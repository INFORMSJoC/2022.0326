# UnitCommitment.jl: Optimization Package for Security-Constrained Unit Commitment
# Copyright (C) 2020-2024, UChicago Argonne, LLC. All rights reserved.
# Released under the modified BSD license. See COPYING.md for more details.

using ADMM, CPLEX, JuMP, UnitCommitment, LinearAlgebra, TimerOutputs, MPI

mutable struct IsfSubproblem
    sp::AdmmSubproblem
    mip
    inj_vars
    w_vars
    e_max_vars
    e_min_vars
    added_constraints
end

function build_subproblem(this_zone::Zone,
                          all_zones::Array{Zone,1},
                          T::Int;
                          gap::Float64,
                          relax::Bool,
                          solver_threads::Int,
                          verbose::Bool)

    optimizer = CPLEX.Optimizer(CPX_PARAM_SCRIND = (verbose ? 1 : 0),
                                CPX_PARAM_THREADS = solver_threads,
                                CPX_PARAM_EPGAP = gap)
    time = 1:T
    buses = [this_zone.buses_int; this_zone.buses_bnd]
    other_zones = setdiff(all_zones, [this_zone])
    model = uc_model(this_zone.generators, buses, T=T, optimizer)
    inj = model.inj_vars

    # Remove centralized power balance equations
    delete.(model.mip, model.power_balance_eq)

    # Add virtual injection variables
    @variable(model.mip, w[z in [zone.index for zone in all_zones],
                           [b.index for b in all_zones[z].buses_bnd],
                           time])

    @variable(model.mip, transfer[time])

    # Safety zone (used for external contingencies)
    @variable(model.mip, 0 <= e_max[[l.index for l in this_zone.lines_int], time] <= 0)
    @variable(model.mip, 0 <= e_min[[l.index for l in this_zone.lines_int], time] <= 0)

    # Add zonal power balance equations
    @constraint(model.mip,
                [t in time],
                sum(inj[b.index,t] for b in this_zone.buses_int) + transfer[t] == 0)

    @constraint(model.mip,
                [t in time],
                sum(w[this_zone.index, b.index, t] for b in this_zone.buses_bnd) == transfer[t])

    # Add definition of virtual injection
    for k in other_zones
        for b in 1:length(k.buses_bnd), t in time
            if this_zone.is_neighbor_zone[k.index]
                @constraint(model.mip,
                            w[k.index, k.buses_bnd[b].index, t] ==
                                    - sum(k.link_base[b, findfirst(isequal(c), k.buses_ext)] * inj[c.index, t]
                                          for c in this_zone.buses_int)
                                    - sum(k.link_base[b, findfirst(isequal(c), k.buses_ext)] * w[this_zone.index, c.index, t]
                                          for c in this_zone.buses_bnd if !(c in k.buses_bnd)))
             else
                 @constraint(model.mip, w[k.index, k.buses_bnd[b].index, t] == 0)
             end
        end
    end

#     if this_zone.index == 2
#         println(model.mip)
#     end

    if relax
        relax_bin_vars!(model.mip)
    end

    vars = [w[k.index, b.index, t] for k in all_zones for b in k.buses_bnd, t in time]
    vars = [vars; [transfer[t] for t in time]]

    ext_weight = length(this_zone.lines_int) < 100 ? 1.0 : 0.0
    weights = [k == this_zone ? 1.0 : ext_weight for k in all_zones for b in k.buses_bnd, t in time]
    weights = [weights; [1.0 for t in time]]

    values = [0.0 for k in all_zones for b in k.buses_bnd, t in time]
    values = [values; [0.0 for t in time]]

    n_vars = length(all_variables(model.mip))

    return IsfSubproblem(AdmmSubproblem(model.mip,
                                        model.obj_var,
                                        vars,
                                        weights=weights,
                                        values=values),
                         model.mip, inj, w, e_max, e_min, Set())
end

function solve_uc_isf(name::String;
                      demand_scale::Float64 = 1.0,
                      limit_scale::Float64 = 1.0,
                      T::Int = 24,
                      ρ::Float64 = 2e-2,
                      transmission::Bool = true,
                      security::Bool = true,
                      relax::Bool = false,
                      max_time::Int = 300,
                      min_feasibility::Float64 = 1e-2,
                      max_iterations = 5,
                      gap::Float64 = 1e-3,
                      solver_threads::Int = 8,
                      verbose_solvers::Array{Int,1} = Array{Int}(undef,0),
                      careful::Bool = false) :: AdmmResult

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm) + 1
    root = (rank == 1)
    nprocs = MPI.Comm_size(comm)

    if root @info "Loading instance..." end
    instance = load_standard_uc_instance(name)
    scale_demands!(instance, demand_scale)
    scale_limits!(instance, limit_scale)

    if root @info "Computing sensitivity factors..." end
    isf = injection_shift_factors(instance.lines)
    lodf = line_outage_factors(instance.lines, isf)
    Z = maximum([l.zone for l in instance.lines])
    if security && Z != 2
        error("security constraints only available for 2-zone systems")
    end

    if root @info "Extracting zones..." end
    zones = [extract_zone(instance, z, isf, lodf, security=security) for z in 1:Z]

    if root
        for z in 1:Z
            @info "Zone $z:"
            @info @sprintf("%8d generators", length(zones[z].generators))
            @info @sprintf("%8d internal lines", length(zones[z].lines_int))
            @info @sprintf("%8d internal buses", length(zones[z].buses_int))
            @info @sprintf("%8d boundary buses", length(zones[z].buses_bnd))
        end
    end

    model = build_subproblem(zones[rank],
                             zones,
                             T,
                             gap=gap,
                             relax = relax,
                             solver_threads = solver_threads,
                             verbose = (rank in verbose_solvers))

    timer = TimerOutput()

    prev_w_base = Dict(t => Array{Float64}(undef,0) for t in 1:T)

    function post_solve()
        if !transmission return end
        this_zone = zones[rank]

        lodf_int = lodf[[l.index for l in this_zone.lines_int],
                        [l.index for l in this_zone.lines_int]]

        L = length(this_zone.lines_int)
        BI = length(this_zone.buses_int)
        BB = length(this_zone.buses_bnd)
        BE = length(this_zone.buses_ext)
        inj_int = value.(model.inj_vars)
        w = value.(model.w_vars)

        isf_int, isf_bnd = this_zone.isf_int, this_zone.isf_bnd
        buses_int, buses_bnd = this_zone.buses_int, this_zone.buses_bnd
        e_max, e_min = zeros(L,T), zeros(L,T)

        isf_int[abs.(isf_int) .< 0.005] .= 0
        isf_bnd[abs.(isf_bnd) .< 0.005] .= 0

        if security
            other_zone = zones[3 - this_zone.index]
            @timeit timer "compute safety zones" begin
                @timeit timer "allreduce inj" begin
                    inj = [0.0 for b in instance.buses, t in 1:T]
                    send = [0.0 for b in instance.buses, t in 1:T]
                    for b in this_zone.buses_int, t in 1:T
                        send[b.index, t] = value(model.inj_vars[b.index, t])
                    end
                    MPI.Allreduce!(send, inj, MPI.SUM, comm)
                end
                for t in 1:T
                    @timeit timer "compute w_base" begin
                        inj_ext_t = [inj[this_zone.buses_ext[c].index, t] for c in 1:BE]
                        w_base = this_zone.link_base * inj_ext_t
                    end
                    if length(prev_w_base[t]) <= 0 || norm(prev_w_base[t] - w_base) > 10.0
                        prev_w_base[t] = w_base
                        @timeit timer "compute diffs" begin
                            diffs = Array{Float64}(undef, BB, 0)
                            for outage in other_zone.lines_int
                                outage.is_vulnerable || continue
                                w_outage = this_zone.link_outage[outage.index] * inj_ext_t
                                diff = w_base - w_outage
                                norm(diff) > 10.0 || continue
                                diffs = [diffs diff]
                            end
                        end
                        @timeit timer "compute e_max, e_min" begin
                            _, D = size(diffs)
                            if D > 0
                                for l in 1:L e_max[l,t], e_min[l,t] = -Inf, Inf end
                                post_flow_diff = isf_bnd * diffs
                                max_diff::Array{Float64,1} = maximum(post_flow_diff, dims=2)[:,1]
                                min_diff::Array{Float64,1} = minimum(post_flow_diff, dims=2)[:,1]
                                e_max[:,t] = max(e_max[:,t], max_diff)
                                e_min[:,t] = min(e_min[:,t], min_diff)
                            end
                        end
                        @timeit timer "update model" begin
                            for l in 1:L
                                line = this_zone.lines_int[l]
                                fix(model.e_max_vars[line.index,t], e_max[l,t], force=true)
                                fix(model.e_min_vars[line.index,t], e_min[l,t], force=true)
                            end
                        end
                    end
                end
            end
        end
        
        @timeit timer "find violations" begin
            for t in 1:T
                max_violation = 0.
                max_violated_line = nothing
                max_violated_scenario = nothing
                inj_int_t = [inj_int[b.index, t] for b in buses_int]
                inj_bnd_t = [w[this_zone.index, b.index, t] for b in buses_bnd]
                pre_flow = isf_int * inj_int_t + isf_bnd * inj_bnd_t
                @timeit timer "pre-contingency" begin
                    for l in 1:L
                        line = this_zone.lines_int[l]
                        violation = max(0.,
                                        pre_flow[l] - line.normal_limit + e_max[l,t],
                                        -pre_flow[l] - line.normal_limit - e_min[l,t])
                        if violation > max_violation
                            max_violation = violation
                            max_violated_line = l
                            max_violated_scenario = nothing
                        end
                    end
                end
                if security
                    @timeit timer "post-contingency" begin
                        limits = [line.normal_limit for line in this_zone.lines_int]
                        for outage in 1:L
                            this_zone.lines_int[outage].is_vulnerable || continue
                            @timeit timer "compute post-flow" begin
                                post_flow = pre_flow + lodf_int[:,outage] * pre_flow[outage]
                            end
                            @timeit timer "evaluate post-flow" begin
                                violations = max.(zeros(L), post_flow-limits, -post_flow-limits)
                                value, monitored = findmax(violations)
                                if value > max_violation
                                    max_violation = value
                                    max_violated_line = monitored
                                    max_violated_scenario = outage
                                end
                            end
                        end
                    end
                end
                if max_violation > 1e-3
                    @timeit timer "add constraints" begin
                        monitored = max_violated_line
                        outage = max_violated_scenario
                        line = this_zone.lines_int[monitored]
                        if (t, monitored, outage) in model.added_constraints continue end
                        push!(model.added_constraints, (t, monitored, outage))
                        flow_monitored = @variable(model.mip)
                        eq = @constraint(model.mip,
                                    flow_monitored == sum(isf_int[monitored,b] * model.inj_vars[buses_int[b].index, t] for b in 1:BI) +
                                                      sum(isf_bnd[monitored,b] * model.w_vars[this_zone.index, buses_bnd[b].index, t] for b in 1:BB))
                        #@show eq
                        if max_violated_scenario === nothing
                            @info(@sprintf("Constraint: line %-8d time %-8d violation %-12.2f",
                                           line.index, t, max_violation))
                            @constraint(model.mip, flow_monitored + model.e_max_vars[line.index, t] <= line.normal_limit)
                            @constraint(model.mip, flow_monitored + model.e_min_vars[line.index, t] >= -line.normal_limit)
                        else
                            @info(@sprintf("Constraint: line %-8d time %-8d violation %-12.2f outage %-8d",
                                           line.index, t, max_violation, max_violated_scenario))
                            flow_outage = @variable(model.mip)
                            @constraint(model.mip,
                                        flow_outage == sum(isf_int[outage,b] * model.inj_vars[buses_int[b].index, t] for b in 1:BI) +
                                                       sum(isf_bnd[outage,b] * model.w_vars[this_zone.index, buses_bnd[b].index, t] for b in 1:BB))
                            @constraint(model.mip, flow_monitored + lodf_int[monitored,outage] * flow_outage <= line.normal_limit)
                            @constraint(model.mip, flow_monitored + lodf_int[monitored,outage] * flow_outage >= -line.normal_limit)
                        end
                    end
                end
            end
        end
    end

    result = sharing_admm(model.sp,
                          comm,
                          ρ = ρ,
                          λ_default = 10.0,
                          post_solve_callback = post_solve,
                          min_improvement = 100.0,
                          min_feasibility = min_feasibility,
                          max_time = max_time,
                          max_iterations = max_iterations,
                          timer = timer,
                          obj_change_tolerance = careful ? 1e-3 : 1e-2,
                          infeas_improv_tolerance = careful ? 1e-4 : 1e-3)


    inj = [0.0 for b in instance.buses, t in 1:T]
    send = [0.0 for b in instance.buses, t in 1:T]
    for b in zones[rank].buses_int, t in 1:T
        send[b.index, t] = value(model.inj_vars[b.index, t])
    end
    MPI.Allreduce!(send, inj, MPI.SUM, comm)

#     if rank == 1
#         @show result.λ
#         @show result.values
#         @show inj
#         @show isf * inj
#         @show isf[2,:]
#     end

    if rank == 1
        @info("Verifying solution...")
        for t in 1:T
            UnitCommitment.verify_power_balance(inj[:,t], t)
            if transmission
                UnitCommitment.verify_flows(inj[:,t], isf, lodf, instance.lines, security=security)
            end
        end

        @info @sprintf("%s,isf,%.2f,%.2f,%d,%d,%.6e,%.6e,%d,%.2f,%.2f",
                       instance.name, demand_scale, limit_scale, transmission,
                       security, result.obj, result.infeasibility, result.iterations,
                       result.wallclock_time, result.wallclock_time / result.iterations)
    end

    return result
end

export solve_uc_isf
