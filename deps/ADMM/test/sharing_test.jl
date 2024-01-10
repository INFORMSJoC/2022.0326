# UnitCommitment.jl: Optimization Package for Security-Constrained Unit Commitment
# Copyright (C) 2020-2024, UChicago Argonne, LLC. All rights reserved.
# Released under the modified BSD license. See COPYING.md for more details.

using ADMM, CPLEX, JuMP, MPI, Test

function build_subproblem(rank::Int)
    opt_factory = CPLEX.Optimizer(CPX_PARAM_SCRIND=0,
                                  CPX_PARAM_THREADS=1)
    if rank == 0
        mip = direct_model(opt_factory)
        @variable(mip, 0 <= x1 <= 2)
        @variable(mip, 0 <= y1 <= 2)
        @variable(mip, z1)
        @variable(mip, obj1)
        @constraint(mip, obj1 == x1 - y1)
        return AdmmSubproblem(mip, obj1, [x1, y1, z1])
    elseif rank == 1
        mip = direct_model(opt_factory)
        @variable(mip, 1 <= x2 <= 3)
        @variable(mip, 1 <= y2 <= 3)
        @variable(mip, z2)
        @variable(mip, nx2)
        @variable(mip, ny2)
        @variable(mip, obj2)
        @constraint(mip, obj2 == 0)
        @constraint(mip, nx2 == - x2)
        @constraint(mip, ny2 == - y2)
        return AdmmSubproblem(mip, obj2, [nx2, ny2, z2])
    end
end

function main()
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    nprocs = MPI.Comm_size(comm)

    sp = build_subproblem(rank)
    result = sharing_admm(sp, comm)

    if rank == 0
        @testset "Sharing ADMM" begin
            @test round(result.obj, digits=3) ≈ -1.0
        end
    end
end

MPI.Init()
main()
MPI.Finalize()
