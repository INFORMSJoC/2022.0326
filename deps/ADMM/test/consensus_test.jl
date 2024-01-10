# UnitCommitment.jl: Optimization Package for Security-Constrained Unit Commitment
# Copyright (C) 2020-2024, UChicago Argonne, LLC. All rights reserved.
# Released under the modified BSD license. See COPYING.md for more details.

using ADMM, CPLEX, JuMP, MPI, Test

function build_subproblem(rank::Int)
    opt_factory = CPLEX.Optimizer(CPX_PARAM_SCRIND=0,
                                  CPX_PARAM_THREADS=1)
    model = direct_model(opt_factory)
    if rank == 1
        @variable(model, 0 <= x1 <= 2)
        @variable(model, 0 <= y1 <= 2)
        @variable(model, obj)
        @constraint(model, obj == x1 - y1)
        return AdmmSubproblem(model, obj, [x1, y1])
    elseif rank == 2
        @variable(model, 1 <= x2 <= 3)
        @variable(model, 1 <= y2 <= 3)
        @variable(model, obj)
        @constraint(model, obj == 0)
        return AdmmSubproblem(model, obj, [x2, y2])
    elseif rank == 3
        @variable(model, 0 <= x3 <= 3)
        @variable(model, 0 <= y3 <= 3)
        @variable(model, obj)
        @constraint(model, obj == 0)
        return AdmmSubproblem(model, obj, [x3, y3])
    end
end

function main()
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    nprocs = MPI.Comm_size(comm)

    sp1 = build_subproblem(1)
    sp2 = build_subproblem(2)
    sp3 = build_subproblem(3)
    result = consensus_admm([sp1, sp2, sp3], comm)

    if rank == 0
        @testset "Consensus ADMM" begin
            @test round(result.obj, digits=3) == -0.333
            #@test round(result.values[1], digits=3) == 1.0
            #@test round(result.values[2], digits=3) == 2.0
        end
    end
end

MPI.Init()
main()
MPI.Finalize()
