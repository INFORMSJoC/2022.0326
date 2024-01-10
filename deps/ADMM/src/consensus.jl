# UnitCommitment.jl: Optimization Package for Security-Constrained Unit Commitment
# Copyright (C) 2020-2024, UChicago Argonne, LLC. All rights reserved.
# Released under the modified BSD license. See COPYING.md for more details.

using Dates, JuMP, LinearAlgebra, Printf, Statistics, TimerOutputs

function consensus_admm(subproblems::Array{AdmmSubproblem,1},
                        comm;
                        callbacks::AdmmCallbacks = AdmmCallbacks(),
                        print_interval::Int = 1,
                        stop_conditions::AdmmStopConditions = AdmmStopConditions(),
                        λ_default::Float64 = 0.0,
                        ρ_default::Float64 = 0.1,
                        ρ_max::Float64 = 1e4,
                        ρ_multiplier::Float64 = 2.0,
                        ρ_update_interval::Int = 25,
                        timer = TimerOutput(),
                       )::AdmmResult

    mpi = MpiInfo(comm)
    nprobs = length(subproblems)
    nvars = length(subproblems[1].vars)
    iterations = Array{IterationInfo,1}(undef, 0)
    params = [AdmmSubparams(ρ = ρ_default,
                            λ = [λ_default for v in 1:nvars],
                            target = [0.0 for v in 1:nvars])
              for n in 1:nprobs]
    print_header(mpi)
    while true
        elapsed_time = @elapsed begin
            solutions = Array{AdmmSubsolution,1}(undef, nprobs)
            for n in 1:nprobs
                @timeit timer "before_solve callback" begin
                  callbacks.before_solve_subproblem(length(iterations) + 1, n)
                end
                @timeit timer "solve" begin
                  solutions[n] = solve_subproblem(subproblems[n], params[n])
                end
                @timeit timer "after_solve callback" begin
                  callbacks.after_solve_subproblem(length(iterations) + 1, n)
                end
                @timeit timer "after_solve barrier" begin
                  MPI.Barrier(mpi.comm)
                end
            end
            global_obj = compute_global_objective(mpi, solutions)
            global_target = compute_global_target(solutions, mpi)
            for n in 1:nprobs 
                update_λ_and_residuals!(solutions[n], params[n], global_target)
            end
            global_infeas = compute_global_infeasibility(global_target, solutions, mpi)
            global_consensus = compute_global_consensus(solutions, mpi)
            if has_numerical_issues(global_target) break end
        end
        total_time = sum([it.iteration_time for it in iterations]) + elapsed_time
        it = IterationInfo(number = length(iterations) + 1,
                           global_obj = global_obj,
                           global_infeas = global_infeas,
                           global_consensus = global_consensus,
                           iteration_time = elapsed_time,
                           total_time = total_time)
        iterations = [iterations; it]
        print_progress(mpi, it, print_interval)
        if should_stop(mpi, iterations, stop_conditions) break end
        for n in 1:nprobs 
            update_ρ!(it, params[n], ρ_update_interval, ρ_max, ρ_multiplier)
        end
        callbacks.after_iteration(it)
    end
    last_iteration = iterations[length(iterations)]
    return AdmmResult(last_iteration.global_obj,
                      nothing,
                      nothing,
                      last_iteration.global_infeas,
                      last_iteration.number,
                      last_iteration.total_time)
end

function compute_global_objective(mpi::MpiInfo,
                                  solutions::Array{AdmmSubsolution,1},
                                 )::Float64

    local_obj = sum([s.obj for s in solutions])
    global_obj = MPI.Allreduce(local_obj, MPI.SUM, mpi.comm)
    global_obj /= (length(solutions) * mpi.nprocs)
    return global_obj
end


function compute_global_target(solutions::Array{AdmmSubsolution,1},
                               mpi::MpiInfo
                              )::Array{Float64,1}

    n_subproblems = length(solutions) * mpi.nprocs
    send = sum([s.values for s in solutions])
    target = MPI.Allreduce(send, MPI.SUM, mpi.comm)
    target /= n_subproblems
    return target
end


function compute_global_consensus(solutions,
                                  mpi::MpiInfo,
                                 )::Float64

    n_vars = length(solutions[1].values)
    local_residual_sum = sum([abs.(s.residuals) for s in solutions])
    global_residual_sum = MPI.Allreduce(local_residual_sum, MPI.SUM, mpi.comm)
    return sum([r <= 1e-3 for r in global_residual_sum]) / n_vars
end


function compute_global_infeasibility(target::Array{Float64,1},
                                      solutions::Array{AdmmSubsolution,1},
                                      mpi::MpiInfo)::Float64

    local_infeasibility = sum([norm(s.residuals) for s in solutions])
    global_infeas = MPI.Allreduce(local_infeasibility, MPI.SUM, mpi.comm)
    return global_infeas
end


function solve_subproblem(sp::AdmmSubproblem,
                          params::AdmmSubparams,
                         )::AdmmSubsolution

    G = length(sp.vars)
    if norm(params.λ) < 1e-3
        @objective(sp.mip, Min, sp.obj)
    else
        @objective(sp.mip, Min,
                   sp.obj
                   + sum(sp.weights[g] * params.λ[g] * (sp.vars[g] - params.target[g]) for g in 1:G)
                   + (params.ρ / 2) * sum(sp.weights[g] * (sp.vars[g] - params.target[g])^2 for g in 1:G))
    end
    optimize!(sp.mip)
    obj = objective_value(sp.mip)
    values = value.(sp.vars)
    return AdmmSubsolution(obj = obj,
                           values = values,
                           residuals = zeros(G))
end


function update_λ_and_residuals!(solution::AdmmSubsolution, 
                                 params::AdmmSubparams,
                                 global_target::Array{Float64,1},
                                )::Nothing

    n_vars = length(solution.values)
    params.target = global_target
    for n in 1:n_vars
        solution.residuals[n] = solution.values[n] - params.target[n]
        params.λ[n] += params.ρ * solution.residuals[n]
    end
end


function update_ρ!(it::IterationInfo,
                   params::AdmmSubparams,
                   ρ_update_interval::Int,
                   ρ_max::Float64,
                   ρ_multiplier::Float64,
                  )::Nothing

    if it.number % ρ_update_interval == 0
        params.ρ = min(ρ_max, params.ρ * ρ_multiplier)
    end
    return 
end


function print_header(mpi::MpiInfo)::Nothing
    if !mpi.root return end
    @info "Solving via Consensus ADMM:"
    @info @sprintf("%8s %20s %20s %14s %8s %8s", "iter", "obj", "infeas", "consensus", "time-it", "time")
end


function print_progress(mpi::MpiInfo,
                        iteration::IterationInfo,
                        print_interval,
                       )::Nothing

    if !mpi.root return end
    if iteration.number % print_interval != 0 return end
    @info @sprintf("%8d %20.6e %20.6e %12.2f %% %8.2f %8.2f",
                  iteration.number,
                  iteration.global_obj,
                  iteration.global_infeas,
                  iteration.global_consensus * 100,
                  iteration.iteration_time,
                  iteration.total_time)
end


function has_numerical_issues(target::Array{Float64,1})::Bool
    if maximum(target) == NaN
        @warn "Numerical issues detected. Stopping."
        return true
    end
    return false
end


function should_stop(mpi::MpiInfo,
                     iterations::Array{IterationInfo,1},
                     conditions::AdmmStopConditions,
                    )::Bool

    if length(iterations) >= conditions.max_iterations
        if mpi.root @info "Iteration limit reached. Stopping." end
        return true
    end

    if length(iterations) < conditions.min_iterations
        return false
    end

    curr_it = iterations[length(iterations)]
    prev_it = iterations[length(iterations) - 1]

    if curr_it.global_infeas < conditions.min_feasibility
        obj_change = abs(prev_it.global_obj - curr_it.global_obj)
        if obj_change < conditions.min_improvement
            if mpi.root @info "Feasibility limit reached. Stopping." end
            return true
        end
    end

    return false
end


export consensus_admm
