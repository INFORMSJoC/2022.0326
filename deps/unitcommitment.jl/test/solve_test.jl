# UnitCommitment.jl: Optimization Package for Security-Constrained Unit Commitment
# Copyright (C) 2020-2024, UChicago Argonne, LLC. All rights reserved.
# Released under the modified BSD license. See COPYING.md for more details.

using UnitCommitment, LinearAlgebra

@testset "Solve" begin
    instance, model, solution = UnitCommitment.solve("ieee_rts/case14",
    	                                             hours = 2,
					    							 gap = 0.0,
					    	                         security = true,
					    	                         transmission = true)
    
    push!(solution.violations, Violation(1, 30, 50, 0.0, 0.0))
    push!(solution.violations, Violation(2, 30, 50, 0.0, 0.0))
    push!(solution.violations, Violation(1, 20, 25, 0.0, 0.0))
    push!(solution.violations, Violation(2, 20, 25, 0.0, 0.0))
    push!(solution.violations, Violation(1, 60, 25, 0.0, 0.0))
    push!(solution.violations, Violation(2, 60, 25, 0.0, 0.0))

    tmpdir = mktempdir()
    csv_filename = "$tmpdir/case14_solution.csv"
    write_csv(csv_filename, instance, model, solution)
    actual_csv_contents = read(csv_filename, String)
    expected_csv_contents = read("fixtures/case14_solution.csv", String)
    @test actual_csv_contents == expected_csv_contents
    recovered_solution = read_uc_solution_csv(csv_filename, instance, T=2)[1]

    ϵ = 0.1
    @test recovered_solution.cost ≈ solution.cost
    @test recovered_solution.is_on == solution.is_on
    @test norm(recovered_solution.production - solution.production) < ϵ
    @test norm(recovered_solution.injection - solution.injection) < ϵ
    @test recovered_solution.violations == solution.violations

    rm(csv_filename)
    rm(tmpdir)
end
