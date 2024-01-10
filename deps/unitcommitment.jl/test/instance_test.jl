# UnitCommitment.jl: Optimization Package for Security-Constrained Unit Commitment
# Copyright (C) 2020-2024, UChicago Argonne, LLC. All rights reserved.
# Released under the modified BSD license. See COPYING.md for more details.

using UnitCommitment, LinearAlgebra

function test_file_equals(f1::String, f2::String)
    open(f1) do original
        open(f2) do generated
            original_lines = readlines(original)
            generated_lines = readlines(generated)
            L = length(original_lines)
            @test length(generated_lines) == L
            for l in 1:L @test original_lines[l] == generated_lines[l] end
        end
    end
end

@testset "Read instance" begin
    instance = load_standard_uc_instance("ieee_rts/case14")
    @test instance.name == "case14"
    @test length(instance.lines) == 20
    @test length(instance.buses) == 14
    @test length(instance.generators) == 5

    @test instance.lines[5].index == 5
    @test instance.lines[5].source == instance.buses[2]
    @test instance.lines[5].target == instance.buses[5]
    @test instance.lines[5].reactance ≈ 0.17388
    @test instance.lines[5].susceptance ≈ 10.037550333
    @test instance.lines[5].normal_limit == 99999
    @test instance.lines[5].emergency_limit == 99999
    @test instance.lines[5].is_vulnerable == true
    @test instance.lines[5].zone == 1

    @test instance.buses[9].index == 9
    @test instance.buses[9].zone == 1
    @test norm(instance.buses[9].demand - [58.94, 55.42, 52.79, 51.91, 51.91, 52.79,
                                           65.1, 75.66, 83.58, 84.46, 84.46, 83.58,
                                           83.58, 83.58, 81.82, 82.7, 87.1, 87.98,
                                           87.98, 84.46, 80.06, 73.02, 64.22, 55.42]) < 0.1

    @test instance.generators[3].index == 3
    @test instance.generators[3].min_power ≈ 0.0
    @test instance.generators[3].max_power ≈ 99.99
    @test instance.generators[3].ramp_up_limit ≈ 70.0
    @test instance.generators[3].ramp_down_limit ≈ 70.0
    @test instance.generators[3].shutdown_ramp_limit ≈ 70.0
    @test instance.generators[3].startup_ramp_limit ≈ 70.0
    @test instance.generators[3].initial_state == 0
    @test instance.generators[3].bus == instance.buses[3]
    @test instance.generators[3].is_always_on == false
    @test instance.generators[3].minimum_uptime == 4
    @test instance.generators[3].minimum_downtime == 4
    @test instance.generators[3].cost_min_power ≈ 0.0
    @test instance.generators[3].cost_segment ≈ [36.28, 40.16, 46.74]
    @test instance.generators[3].offer_segment ≈ [33.33, 33.33, 33.33]
    @test instance.generators[3].startup_cost ≈ 0.0

    @test installed_capacity(instance) ≈ 772.38
end

# @testset "Write instance" begin
#     temp = tempdir()
#     base = "$(dirname(@__FILE__))/../instances/pegase/case89pegase/"
#     instance = load_standard_uc_instance("pegase/case89pegase")
#     write_uc_instance(instance, tempdir())
#     test_file_equals("$base/generators.csv", "$temp/generators.csv")
#     test_file_equals("$base/lines.csv", "$temp/lines.csv")
#     test_file_equals("$base/buses.csv", "$temp/buses.csv")
# end
