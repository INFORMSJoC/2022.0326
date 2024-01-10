# UnitCommitment.jl: Optimization Package for Security-Constrained Unit Commitment
# Copyright (C) 2020-2024, UChicago Argonne, LLC. All rights reserved.
# Released under the modified BSD license. See COPYING.md for more details.

using JuMP, Requires, Printf, MathOptFormat, TimerOutputs, LinearAlgebra
import Base.Threads: @threads, nthreads


struct InfeasibleProblemException <: Exception end


struct SolveStats
    wallclock_time::Float64
    iterations::Int
    termination_status
end


function relax_bin_vars!(model)
    bin_vars = [v for v in all_variables(model) if is_binary(v)]
    unset_binary.(bin_vars)
    set_upper_bound.(bin_vars, 1.0)
    set_lower_bound.(bin_vars, 0.0)
end


function add_constraint(instance::UnitCommitmentInstance,
                        model::UnitCommitmentModel,
                        v::Violation,
                        isf::Array{Float64,2},
                        lodf::Array{Float64,2};
                        output_prefix::String = "")

    model.violations = [model.violations; v]
    if v.time == 1
        if v.outage_line == v.monitored_line
            @info @sprintf("%s    %8.3f MW overflow in line %-5d (pre-contingency)",
                          output_prefix, v.amount, v.monitored_line)
        else
            @info @sprintf("%s    %8.3f MW overflow in line %-5d (outage: line %d)",
                          output_prefix, v.amount, v.monitored_line, v.outage_line)
        end
    end
    flow = @variable(model.mip)
    set_lower_bound(flow, -v.limit)
    set_upper_bound(flow, v.limit)
    inj = model.inj_vars
    buses = instance.buses
    B = length(buses)
    if v.outage_line == v.monitored_line
        @constraint(model.mip, flow ==
            sum(isf[v.monitored_line, b] * inj[buses[b].index, v.time] for b in 1:B))
    else
        @constraint(model.mip, flow ==
            sum((isf[v.monitored_line, b] + lodf[v.monitored_line, v.outage_line] * isf[v.outage_line, b]) * inj[buses[b].index, v.time]
                for b in 1:B))
    end
end

function find_violations(instance::UnitCommitmentInstance,
                         inj;
                         hours::Int,
                         isf::Array{Float64,2},
                         lodf::Array{Float64,2},
                         max_per_line::Int = 1,
                         max_per_period::Int = 5,
                         security::Bool = true,
                        ) :: Array{Violation}
    L = length(instance.lines)
    normal_limits = [l.normal_limit for l in instance.lines]
    emergency_limits = [l.emergency_limit for l in instance.lines]
    local_filters = Dict(t => ViolationFilter(L,
                                              max_per_period=max_per_period,
                                              max_per_line=max_per_line)
                         for t in 1:hours)
    @threads for t in 1:hours
        injt::Array{Float64} = [inj[b.index, t] for b in instance.buses]
        pre_flow::Array{Float64} = isf * injt
        for m in 1:L
            vt::Float64 = max(0, pre_flow[m] - normal_limits[m], - pre_flow[m] - normal_limits[m])
            if vt > 1e-5
                offer(local_filters[t], Violation(t, m, m, vt, normal_limits[m]))
            end
        end
        if security
            for c in instance.lines
                c.is_vulnerable || continue
                for m in 1:L
                    post_flow_m::Float64 = pre_flow[m] + lodf[m, c.index] * pre_flow[c.index]
                    vt::Float64 = max(0, post_flow_m - emergency_limits[m], - post_flow_m - emergency_limits[m])
                    if vt > 1e-5
                        offer(local_filters[t], Violation(t, m, c.index, vt, emergency_limits[m]))
                    end
                end
            end
        end
    end
    global_filter = ViolationFilter(length(instance.lines),
                                    max_per_period=max_per_period,
                                    max_per_line=max_per_line)
    for t in 1:hours
        for v in query(local_filters[t])
            offer(global_filter, v)
        end
    end
    violations = query(global_filter)
    return violations
end



function update_constraints(instance::UnitCommitmentInstance,
                            model::UnitCommitmentModel;
                            isf::Array{Float64,2},
                            lodf::Array{Float64,2},
                            max_per_line::Int = 1,
                            max_per_period::Int = 5,
                            output_prefix::String = "",
                            security::Bool = true,
                            should_add_contraints::Bool = true,
                            T::Int = 24,
                            timer = TimerOutput(),
                           )

    @info "$(output_prefix)Verifying flow constraints..."
    inj = value.(model.inj_vars)
    violations = find_violations(instance,
                                 inj,
                                 hours = model.hours,
                                 isf = isf,
                                 lodf = lodf,
                                 max_per_line = max_per_line,
                                 max_per_period = max_per_period,
                                 security = security)
    @info "$(output_prefix)    $(length(violations)) violations"
    if should_add_contraints
        for t in 1:T, v in violations
            vt = Violation(t, v.monitored_line, v.outage_line, v.amount, v.limit)
            add_constraint(instance, model, vt, isf, lodf, output_prefix=output_prefix)
        end
    end
    return length(violations)
end


"""
    solve(instance_name; kw_args)

Solves a benchmark instance of the Unit Commitment Problem using the baseline solver.

## Keyword Arguments
- `demand_scale::Float64`   Scale all bus demands linearly by the given amount [default: 0.6]
- `gap::Float64`            Relative MIP integrality gap [default: 1e-3]
- `hours::Int`              Number of hours in the planning horizon [default: 24]
- `isf_cutoff::Float64`     Truncate ISF entries with magnitude smaller than the given amount [default: 1e-3]
- `limit_scale::Float64`    Scale all thermal limits linearly by the given amount [default: 1.0]
- `lodf_cutoff::Float64`    Truncate LODF entries with magnitude smaller than the given amount [default: 1e-4]
- `relax::Bool`             Solve linear relaxation of the problem [default: false]
- `security::Bool`          Enforce N-1 security constraints [default: false]
- `threads::Int`            Number of CPU threads to use [default: 4]
- `transmission::Bool`      Enforce transmission constraints [default: false]
- `reserve::Bool`           Amount of system-wide spinning reserves, proportional to the demand [default: 0.01]
- `verbose::Bool`           Show solver log [default: false]
"""
function solve(instance_name::String;
               demand_scale::Float64 = 0.6,
               gap::Float64 = 1e-3,
               hours::Int = 24,
               isf_cutoff::Float64 = 0.005,
               limit_scale::Float64 = 1.0,
               lodf_cutoff::Float64 = 0.001,
               relax::Bool = false,
               reserve::Float64 = 0.01,
               security::Bool = false,
               threads::Int = 4,
               transmission::Bool = false,
               verbose::Bool = false,
               export_lp_filename::Union{Nothing, String} = nothing,
               export_mps_filename::Union{Nothing, String} = nothing,
              ) :: Tuple{UnitCommitmentInstance,
                         UnitCommitmentModel,
                         UnitCommitmentSolution,
                         SolveStats}

    @info "Loading instance: $instance_name"
    instance = load_standard_uc_instance(instance_name)
    print_summary(instance)

    @info "Scaling problem ($demand_scale demands, $limit_scale limits)..."
    scale_demands!(instance, demand_scale)
    scale_limits!(instance, limit_scale)

    @info "Computing sensitivity factors ($isf_cutoff ISF cutoff, $lodf_cutoff LODF cutoff)..."
    isf = injection_shift_factors(instance.lines)
    lodf = line_outage_factors(instance.lines, isf)
    isf[abs.(isf) .< isf_cutoff] .= 0
    lodf[abs.(lodf) .< lodf_cutoff] .= 0

    model = build_model(instance,
	                 	threads = threads,
	                    reserve = reserve,
	                    verbose = verbose,
    	                gap = gap,
    	                hours = hours)

    solution, stats =  solve(instance,
    	                     model,
    	                     hours = hours,
    	                     isf = isf,
    	                     lodf = lodf,
    	                     relax = relax,
    	                     security = security,
    	                     transmission = transmission,
                             export_lp_filename = export_lp_filename,
                             export_mps_filename = export_mps_filename,
                             gap = gap,
                            )

    return instance, model, solution, stats
end


function solve(instance::UnitCommitmentInstance,
               model::UnitCommitmentModel;
               hours::Int = 24,
               relax::Bool = false,
               isf::Array{Float64,2},
               lodf::Array{Float64,2},
               security::Bool = false,
               transmission::Bool = false,
               output_prefix::String = "",
               export_lp_filename::Union{Nothing, String} = nothing,
               export_mps_filename::Union{Nothing, String} = nothing,
               enable_two_phase_gap::Bool = true,
               max_iterations::Int = 1000,
               timer = TimerOutput(),
               before_optimize_callback = nothing,
               after_optimize_callback = nothing,
               gap::Float64,
	          )::Tuple{UnitCommitmentSolution, SolveStats}

    if relax
        @info "$(output_prefix)Relaxing integrality..."
        relax_bin_vars!(model.mip)
    end

    iterations = 0
    solution_time = 0
    large_gap = false
    relaxed_gap  = 1e-2
    original_gap = gap

    if enable_two_phase_gap && two_phase_gap_available && (transmission || security)
        large_gap = true
        @info "$(output_prefix)Relaxing gap to $(relaxed_gap)..."
        set_gap!(model.mip, relaxed_gap)
    end

    solution_time = @elapsed begin
        @timeit timer "Solve Unit Commitment" begin
            while true
                iterations += 1
                @timeit timer "Optimize MILP" begin
                    @info "$(output_prefix)Optimizing..."
                    before_optimize_callback == nothing || before_optimize_callback()
                    mip_time = @elapsed optimize!(model.mip)
                    termination_status(model.mip) == MOI.OPTIMAL || throw(InfeasibleProblemException())
                    obj = objective_value(model.mip)
                    @info @sprintf("%s    obj: %.8e", output_prefix, obj)
                    @info @sprintf("%s   time: %.2f s", output_prefix, mip_time)
                    after_optimize_callback == nothing || after_optimize_callback()
                end

                @timeit timer "Contingency screening" begin
                    n_violations = 0
                    if transmission || security
                        n_violations = update_constraints(instance,
                                                          model,
                                                          isf = isf,
                                                          lodf = lodf,
                                                          T = hours,
                                                          security = security,
                                                          output_prefix = output_prefix,
                                                          should_add_contraints = (iterations < max_iterations),
                                                          timer = timer)
                    end
                end

                if iterations >= max_iterations
                    @info "$(output_prefix)Maximum number of iterations reached. Stopping."
                    break
                end

                if n_violations == 0
                    if large_gap == false
                        break
                    else
                        large_gap = false
                        @info "$(output_prefix)Restoring gap to $(original_gap)..."
                        set_gap!(model.mip, original_gap)
                    end
                end
            end
        end
    end

    if export_lp_filename != nothing
        @info "Exporting: $export_lp_filename"
        lp_model = MathOptFormat.LP.Model()
        MOI.copy_to(lp_model, backend(model.mip))
        MOI.write_to_file(lp_model, export_lp_filename)
    end

    if export_mps_filename != nothing
        @info "Exporting: $export_mps_filename"
        mps_model = MathOptFormat.MPS.Model()
        MOI.copy_to(mps_model, backend(model.mip))
        MOI.write_to_file(mps_model, export_mps_filename)
    end

    @info @sprintf("%sSolved in %.2f seconds", output_prefix, solution_time)
    stats = SolveStats(solution_time, iterations, termination_status(model.mip))

    return get_solution(instance, model), stats
end

function build_model(instance::UnitCommitmentInstance;
                     gap::Float64 = 1e-3,
                     hours::Int = 24,
                     output_prefix::String = "",
                     reserve::Float64 = 0.01,
                     seed::Int = 42,
                     threads::Int = 4,
                     verbose::Bool = false,
                    )::UnitCommitmentModel

    @info "$(output_prefix)Using $solver_name as MILP solver ($gap gap, $threads threads)"
    solver = get_solver(verbose=verbose, threads=threads, gap=gap, seed=seed)

    @info "$(output_prefix)Building MILP model ($hours hours, $reserve reserve)..."
    model = uc_model(instance.generators, instance.buses, T=hours, solver, reserve_pct=reserve)
end

export update_constraints, UnitCommitmentSolution
