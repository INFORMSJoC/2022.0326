# UnitCommitment.jl: Optimization Package for Security-Constrained Unit Commitment
# Copyright (C) 2020-2024, UChicago Argonne, LLC. All rights reserved.
# Released under the modified BSD license. See COPYING.md for more details.

using ADMM, CPLEX, JuMP, MPI, Test # Precompile all packages
 
function run_mpi(filename)
    run(`mpiexec -n 2 julia --color=yes $filename`)
end

@testset "ADMM" begin
    @testset "MPI" begin
        run_mpi("sharing_test.jl")
        run_mpi("consensus_test.jl")
    end
end