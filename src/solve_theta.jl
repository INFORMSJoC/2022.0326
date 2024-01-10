# UnitCommitment.jl: Optimization Package for Security-Constrained Unit Commitment
# Copyright (C) 2020-2024, UChicago Argonne, LLC. All rights reserved.
# Released under the modified BSD license. See COPYING.md for more details.

using ADMM, CPLEX, JuMP, UnitCommitment, LinearAlgebra, TimerOutputs, MPI

mutable struct ThetaSubproblem
    sp::AdmmSubproblem
    inj_vars
end

function build_subproblem(zone::Zone,
                          T::Int;
                          gap::Float64,
                          relax::Bool,
                          solver_threads::Int,
                          verbose::Bool) :: ThetaSubproblem

    optimizer = CPLEX.Optimizer(CPXPARAM_Advance = 0,
                                CPX_PARAM_BAREPCOMP = 1e-3,
                                CPXPARAM_Emphasis_Numerical = 1,
                                CPX_PARAM_SCRIND = (verbose ? 1 : 0),
                                CPX_PARAM_THREADS = solver_threads,
                                CPX_PARAM_EPGAP = gap)

    time = 1:T
    buses = [zone.buses_int; zone.buses_bnd]
    model = uc_model(zone.generators, buses, T=T, optimizer)
    inj = model.inj_vars

    # Remove centralized power balance equations
    delete.(model.mip, model.power_balance_eq)

    @variable(model.mip, theta[buses, time])
    @variable(model.mip, z[zone.buses_bnd, time])
    @variable(model.mip, w[zone.buses_bnd, time])
    @variable(model.mip, -l.normal_limit <= flow[l in zone.lines_int, time] <= l.normal_limit)

    # Definition of flow
    @constraint(model.mip,
                [l in zone.lines_int, t in time],
                flow[l, t] == l.susceptance * (theta[l.source, t] - theta[l.target, t]))

    # Flow conservation (interior buses)
    @constraint(model.mip,
                [b in zone.buses_int, t in time],
                - sum(flow[l, t] for l in zone.lines_int if l.source == b)
                + sum(flow[l, t] for l in zone.lines_int if l.target == b)
                + inj[b.index,t] == 0)

    # Flow conservation (boundary buses)
    @constraint(model.mip,
                [b in zone.buses_bnd, t in time],
                - sum(flow[l, t] for l in zone.lines_int if l.source == b)
                + sum(flow[l, t] for l in zone.lines_int if l.target == b)
                + w[b,t] + inj[b.index,t] == 0)

    @constraint(model.mip,
                [b in zone.buses_bnd, t in time],
                z[b,t] == theta[b,t] * (zone.index == 1 ? 1 : -1))

    if relax
        relax_bin_vars!(model.mip)
    end

    vars = [w[b,t] for b in zone.buses_bnd for t in time]
    vars = [vars; [z[b,t] for b in zone.buses_bnd for t in time]]

    return ThetaSubproblem(AdmmSubproblem(model.mip, model.obj_var, vars),
                           inj)
end

function solve_uc_theta(name::String;
                        demand_scale::Float64 = 1.0,
                        limit_scale::Float64 = 1.0,
                        T::Int = 24,
                        ρ::Float64 = 1e-2,
                        relax::Bool = false,
                        max_time::Int = 300,
                        min_feasibility::Float64 = 10.0,
                        gap::Float64 = 1e-2,
                        solver_threads::Int = 8,
                        verbose_solvers=[]) :: AdmmResult

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm) + 1
    root = (rank == 1)
    nprocs = MPI.Comm_size(comm)

    instance = load_standard_uc_instance(name)
    scale_demands!(instance, demand_scale)
    scale_limits!(instance, limit_scale)
    isf = injection_shift_factors(instance.lines)
    lodf = line_outage_factors(instance.lines, isf)
    Z = maximum([l.zone for l in instance.lines])
    if Z != 2
        @error("Instance must have exactly 2 zones. Found $(Z) zones instead.")
        return
    end

    zones = [extract_zone(instance, z, isf, lodf, security=false) for z in 1:Z]
    model = build_subproblem(zones[rank],
                              T,
                              gap=gap,
                              relax = relax,
                              solver_threads = solver_threads,
                              verbose = (rank in verbose_solvers))

    result = sharing_admm(model.sp,
                          comm,
                          ρ = ρ,
                          ρ_max = 1.0,
                          ρ_update_interval = 25,
                          λ_default = 10.0,
                          min_feasibility = min_feasibility,
                          max_time = max_time,
                          obj_change_tolerance = 1e-3,
                          infeas_improv_tolerance = 1e-4)

    inj = [0.0 for b in instance.buses, t in 1:T]
    send = [0.0 for b in instance.buses, t in 1:T]
    for b in zones[rank].buses_int, t in 1:T
        send[b.index, t] = value(model.inj_vars[b.index, t])
    end
    MPI.Allreduce!(send, inj, MPI.SUM, comm)

    if rank == 1
        @info("Verifying solution...")
        for t in 1:T
            UnitCommitment.verify_power_balance(inj[:,t], t)
            UnitCommitment.verify_flows(inj[:,t], isf, lodf, instance.lines, security=false)
        end

        @info @sprintf("%s,theta,%.2f,%.2f,%d,%d,%.6e,%.6e,%d,%.2f,%.2f",
                       instance.name, demand_scale, limit_scale, 1, 0, result.obj,
                       result.infeasibility, result.iterations, result.wallclock_time,
                       result.wallclock_time / result.iterations)
   end

    return result
end

export solve_uc_theta
