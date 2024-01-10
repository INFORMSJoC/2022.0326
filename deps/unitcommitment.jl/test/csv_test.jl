# UnitCommitment.jl: Optimization Package for Security-Constrained Unit Commitment
# Copyright (C) 2020-2024, UChicago Argonne, LLC. All rights reserved.
# Released under the modified BSD license. See COPYING.md for more details.

using UnitCommitment

@testset "Read CSV" begin
    csv = read_csv("../instances/ieee_rts/case14/lines.csv")
    @test length(csv.column_names) == 8
    @test csv.column_names[1] == "Line"
    @test csv.column_names[3] == "Target"
    @test csv.column_names[8] == "Zone"

    @test csv.data[1, 1] == "1"
    @test csv.data[1, 2] == "1"
    @test csv.data[1, 3] == "2"
    @test csv.data[1, 4] == "0.05917000"
    @test csv.data[1, 5] == "99999"
    @test csv.data[1, 6] == "99999"
    @test csv.data[1, 7] == "1"
    @test csv.data[1, 8] == "1"

    @test csv.data[20, 1] == "20"
    @test csv.data[20, 2] == "13"
    @test csv.data[20, 3] == "14"
    @test csv.data[20, 4] == "0.34802000"
    @test csv.data[20, 5] == "99999"
    @test csv.data[20, 6] == "99999"
    @test csv.data[20, 7] == "1"
    @test csv.data[20, 8] == "1"
end
