# UnitCommitment.jl: Optimization Package for Security-Constrained Unit Commitment
# Copyright (C) 2020-2024, UChicago Argonne, LLC. All rights reserved.
# Released under the modified BSD license. See COPYING.md for more details.

module UnitCommitment
    include("csv.jl")
    include("instance.jl")
    include("sensitivity.jl")
    include("model.jl")
    include("dc_flow.jl")
    include("checker.jl")
    include("solution.jl")
    include("solve.jl")

    using Cbc, Requires
    solver_name = "Cbc"
    two_phase_gap_available = false
    function get_solver(;verbose::Bool, gap::Float64, threads::Int, seed::Int)
        return with_optimizer(Cbc.Optimizer, logLevel=(verbose ? 1 : 0), ratioGap=gap)
    end

    function __init__()
        @require SCIP="82193955-e24f-5292-bf16-6f2c5261a85f" begin
            solver_name = "SCIP"
            function get_solver(;verbose::Bool, gap::Float64, threads::Int, seed::Int)
                return SCIP.Optimizer(display_verblevel=(verbose ? 5 : 1),
                                      numerics_feastol = 1e-5,
                                      limits_gap=gap)
            end
        end
        @require CPLEX="a076750e-1247-5638-91d2-ce28b192dca0" begin
            solver_name = "CPLEX"
            two_phase_gap_available = true
            function get_solver(;verbose::Bool, gap::Float64, threads::Int, seed::Int)
                return CPLEX.Optimizer(CPX_PARAM_SCRIND=(verbose ? 1 : 0),
                                       CPX_PARAM_MIPDISPLAY=4,
                                       CPX_PARAM_THREADS=threads,
                                       CPX_PARAM_EPGAP=gap,
                                       CPX_PARAM_RANDOMSEED=seed)
            end
            function set_gap!(model, gap::Float64)
                CPLEX.set_param!(backend(model).inner.env, "CPX_PARAM_EPGAP", gap)
            end
        end
    end
end
