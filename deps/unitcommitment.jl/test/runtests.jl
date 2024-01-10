# UnitCommitment.jl: Optimization Package for Security-Constrained Unit Commitment
# Copyright (C) 2020-2024, UChicago Argonne, LLC. All rights reserved.
# Released under the modified BSD license. See COPYING.md for more details.

using Test

@testset "UnitCommitment" begin
    include("csv_test.jl")
    include("instance_test.jl")
    include("sensitivity_test.jl")
    include("solve_test.jl")
end