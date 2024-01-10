# UnitCommitment.jl: Optimization Package for Security-Constrained Unit Commitment
# Copyright (C) 2020-2024, UChicago Argonne, LLC. All rights reserved.
# Released under the modified BSD license. See COPYING.md for more details.

using JuMP, TimerOutputs


mutable struct AdmmSubproblem
    mip::JuMP.Model
    obj::VariableRef
    vars::Array{VariableRef,1}
    weights::Array{Float64,1}
    values::Array{Float64,1}
    timer::TimerOutput
end


function AdmmSubproblem(mip::JuMP.Model,
                        obj::VariableRef,
                        vars::Array{VariableRef, 1};
                        weights::Array{Float64, 1} = [1.0 for v in vars],
                        values::Array{Float64, 1} = [0.0 for v in vars],
                        timer::TimerOutput = TimerOutput(),
                       )::AdmmSubproblem

    return AdmmSubproblem(mip,
                          obj,
                          vars,
                          weights,
                          values,
                          timer)
end


struct AdmmSubsolution
    obj::Float64
    values::Array{Float64,1}
    residuals::Array{Float64,1}
end


function AdmmSubsolution(;
                         obj::Float64,
                         values::Array{Float64,1},
                         residuals::Array{Float64,1},
                        )::AdmmSubsolution

    return AdmmSubsolution(obj, values, residuals)
end


mutable struct AdmmSubparams
    ρ::Float64
    λ::Array{Float64,1}
    target::Array{Float64,1}
end


function AdmmSubparams(;
                       ρ::Float64,
                       λ::Array{Float64,1},
                       target::Array{Float64,1},
                      )::AdmmSubparams

    return AdmmSubparams(ρ, λ, target)
end


struct MpiInfo
    comm
    rank::Int
    root::Bool
    nprocs::Int
end


struct AdmmStopConditions
    max_iterations::Int
    max_time::Int
    min_feasibility::Float64
    min_improvement::Float64
    min_iterations::Int
end


function AdmmStopConditions(;
                            max_iterations::Int = 1000,
                            max_time::Int = 900,
                            min_feasibility::Float64 = 1e-3,
                            min_improvement::Float64 = 1e-3,
                            min_iterations::Int = 2,
                           )::AdmmStopConditions

    return AdmmStopConditions(max_iterations,
                              max_time,
                              min_feasibility,
                              min_improvement,
                              min_iterations)
end


mutable struct IterationInfo
    number::Int
    global_obj::Float64
    global_infeas::Float64
    global_consensus::Float64
    iteration_time::Float64
    total_time::Float64
end


function IterationInfo(;
                       number::Int,
                       global_obj::Float64,
                       global_infeas::Float64,
                       global_consensus::Float64,
                       iteration_time::Float64,
                       total_time::Float64,
                      )::IterationInfo

    return IterationInfo(number,
                         global_obj,
                         global_infeas,
                         global_consensus,
                         iteration_time,
                         total_time)
end


struct AdmmCallbacks
    before_solve_subproblem
    after_solve_subproblem
    after_iteration
end


function AdmmCallbacks(;
                       before_solve_subproblem = (iteration::Int, problem::Int) -> nothing,
                       after_solve_subproblem  = (iteration::Int, problem::Int) -> nothing,
                       after_iteration         = (it::IterationInfo)-> nothing,
                      )::AdmmCallbacks

    return AdmmCallbacks(before_solve_subproblem,
                         after_solve_subproblem,
                         after_iteration)
end 


function MpiInfo(comm)
    rank = MPI.Comm_rank(comm)+1
    is_root = (rank == 1)
    nprocs = MPI.Comm_size(comm)
    return MpiInfo(comm, rank, is_root, nprocs)
end


export AdmmSubproblem, AdmmSubsolution, AdmmSubparams, MpiInfo, IterationInfo,
       AdmmCallbacks, AdmmStopConditions
