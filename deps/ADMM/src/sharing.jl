# UnitCommitment.jl: Optimization Package for Security-Constrained Unit Commitment
# Copyright (C) 2020-2024, UChicago Argonne, LLC. All rights reserved.
# Released under the modified BSD license. See COPYING.md for more details.

using MPI, Dates, JuMP, LinearAlgebra, Printf, Statistics, TimerOutputs, CPLEX
const MOI = JuMP.MathOptInterface

struct AdmmResult
    obj::Float64
    λ
    values
    infeasibility::Float64
    iterations::Int
    wallclock_time::Float64
end

struct AdmmSubproblemResult
    obj::Float64
    values::Array{Float64,1}
end

function switch_to_qp(sp::AdmmSubproblem, bin_vars::Array{VariableRef})
    cpx = backend(sp.mip).inner
    values = value.(bin_vars)
    for b in 1:length(bin_vars)
        unset_binary(bin_vars[b])
        set_upper_bound(bin_vars[b], round(values[b]))
        set_lower_bound(bin_vars[b], round(values[b]))
    end
    CPLEX.set_prob_type!(cpx, :QP)
    cpx.has_int = false
    cpx.has_qc = true
end

function switch_to_miqp(sp::AdmmSubproblem, bin_vars::Array{VariableRef})
    cpx = backend(sp.mip).inner
    CPLEX.set_prob_type!(cpx, :MIQP)
    cpx.has_int = true
    cpx.has_qc = true
    set_binary.(bin_vars)
end

function optimize(mip)::Bool
    try
        optimize!(mip)
        return true
    catch
        @warn "An exception occurred while optimizing"
        return false
    end
end

function sharing_admm(sp::AdmmSubproblem,
                      comm;
                      ρ::Float64 = 1.0,
                      ρ_update_interval::Int = 5,
                      ρ_multiplier::Float64 = 2.0,
                      ρ_max::Float64 = 1e3,
                      λ::Array{Float64,1} = [0.0],
                      λ_default::Float64 = 0.0,
                      target::Array{Float64,1} = [0.0],
                      max_iterations::Int = 1000000,
                      min_iterations::Int = 10,
                      min_feasibility::Float64 = 1e-3,
                      min_improvement::Float64 = 1e-3,
                      print_interval::Int = 1,
                      max_time::Int = 900,
                      post_solve_callback = (()->()),
                      timer=TimerOutput(),
                      obj_change_tolerance::Float64 = 1e-2,
                      infeas_improv_tolerance::Float64 = 1e-3,
                      verbose=true,
                     ) :: AdmmResult

    rank = MPI.Comm_rank(comm)+1
    root = (rank == 1)
    N = MPI.Comm_size(comm)
    remaining_time = max_time

    if root
        @info "Solving via Sharing ADMM:"
        @info @sprintf("%8s %20s %20s %8s %8s %8s", "iter", "obj", "infeas", "time-it", "time", "type")
    end

     @timeit timer "sharing ADMM" begin
        @timeit timer "initialization" begin
            G::Int = length(sp.vars)
            if length(λ) == 1 λ = [λ_default for g in 1:G] end
            if length(target) == 1 target = [0. for g in 1:G] end
            values::Array{Float64,1} = [target[g] for g in 1:G]
            total_obj::Float64, prev_obj::Float64 = 0.0, 0.0
            infeas::Float64, prev_infeas::Float64 = 0.0, 0.0
            elapsed_time::Float64 = 0.0
            obj::Float64 = 0.0
            iteration::Int = 0
        end

        bin_vars = [v for v in all_variables(sp.mip) if is_binary(v)]
        switch_to_miqp(sp, bin_vars)
        ptype = "MIQP"
        
        while iteration < max_iterations

            iteration += 1
            initial_solve_time::DateTime = now()

            @timeit timer "solve subproblem" begin
                if ptype == "MIQP"
                    @objective(sp.mip, Min,
                               sp.obj
                               + sum(sp.weights[g] * λ[g] * sp.vars[g] for g in 1:G)
                               + (ρ / 2) * sum(sp.weights[g] * (sp.vars[g] - values[g] + target[g])^2 for g in 1:G))
                else
                    @objective(sp.mip, Min,
                               sp.obj
                               + sum(sp.weights[g] * λ[g] * sp.vars[g] for g in 1:G)
                               + (ρ / 2) * sum(1.0 * (sp.vars[g] - values[g] + target[g])^2 for g in 1:G))
                end

                cpx = backend(sp.mip).inner
                CPLEX.set_param!(cpx.env, "CPX_PARAM_TILIM", remaining_time)

                solve_time = @elapsed begin
                    is_successful = optimize(sp.mip) 
                end
                if is_successful
                    status = termination_status(sp.mip)
                    if status == MOI.OPTIMAL || status == MOI.OTHER_ERROR
                        obj = objective_value(sp.mip)
                        values = value.(sp.vars)
                    else
                        @warn "Subproblem $rank not solved to optimality ($status). Using previous solution."
                    end
                end
            end

            @timeit timer "post-optimize barrier" begin
                MPI.Barrier(comm)
            end

            @timeit timer "post-solve callback" begin
                post_solve_callback()
            end

            @timeit timer "reduce values, obj & time" begin
                solve_time = MPI.Allreduce(solve_time, MPI.MAX, comm)
                remaining_time -= solve_time
                prev_obj, total_obj = total_obj, MPI.Allreduce(obj, MPI.SUM, comm)
                MPI.Allreduce!(values, target, MPI.SUM, comm)
                target /= N
            end

            if maximum(target) == NaN
                if root @warn "Numerical issues detected. Stopping." end
                break
            end

            @timeit timer "update multipliers" begin
                prev_infeas, infeas = infeas, norm(target)
                for g in 1:G λ[g] += ρ * target[g] end
            end

            if root
                if iteration % print_interval == 0
                    @info @sprintf("%8d %20.6e %20.6e %8.2f %8.2f %8s",
                                   iteration, total_obj, infeas, solve_time,
                                   max_time - remaining_time, ptype)
                end
            end

            if remaining_time <= 0
                if root @info "Maximum time limit reached. Stopping." end
                break
            end

            if iteration > min_iterations && infeas < min_feasibility
                if root @info "Feasibility tolerance reached. Stopping." end
                break
            end

            if iteration % ρ_update_interval == 0
                ρ = min(ρ_max, ρ * ρ_multiplier)
            end

            if ptype == "MIQP"
                obj_change = abs((prev_obj - total_obj) / total_obj)
                if obj_change < obj_change_tolerance
                    if root @info "Objective value stagnated. Switching to QP." end
                    switch_to_qp(sp, bin_vars)
                    ptype = "QP"
                end
            else
                infeas_improv = abs((prev_infeas - infeas) / infeas)
                if infeas_improv < infeas_improv_tolerance
                    if root @info "Infeasibility stagnated. Switching to MIQP." end
                    switch_to_miqp(sp, bin_vars)
                    ptype = "MIQP"
                end
            end
        end
    end

    if iteration == max_iterations
        if root @info "Iteration limit reached. Stopping." end
    end

    # for n in 1:N
    #     if n == rank println(timer) end
    #     MPI.Barrier(comm)
    # end

    return AdmmResult(total_obj, λ, values, infeas, iteration, max_time - remaining_time)
end

export sharing_admm, AdmmResult
