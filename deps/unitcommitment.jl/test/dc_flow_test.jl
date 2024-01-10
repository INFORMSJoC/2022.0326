# UnitCommitment.jl: Optimization Package for Security-Constrained Unit Commitment
# Copyright (C) 2020-2024, UChicago Argonne, LLC. All rights reserved.
# Released under the modified BSD license. See COPYING.md for more details.

using UnitCommitment, Test

@testset "Pre-contingency flow" begin
    instance = load_standard_uc_instance("split/13b-2z")
    isf = injection_shift_factors(instance.lines)
    B = length(instance.buses)
    inj = [1.0 for b in 1:B]
    inj[1] = -sum(inj[2:B])

    expected = [-7.03861, -4.96139, -4.15444, -1.88417, -3.96139,
                -0.0019305, -1.27606, -1.87645, -2.2722, -2.57336,
                -1.27413, -1.57336, -1.55019, -1.2249, -1.2249,
                -0.275097, -0.275097, -0.5, 0.5]

    actual = pre_contingency_flow(inj, isf)
    @test norm(expected - actual) < 1e-4
end

@testset "Post-contingency flow" begin
    instance = load_standard_uc_instance("split/13b-2z")
    isf = injection_shift_factors(instance.lines)
    lodf = line_outage_factors(instance.lines, isf)
    B = length(instance.buses)
    inj = [1.0 for b in 1:B]
    inj[1] = -sum(inj[2:B])
    outage = 6

    instance.lines[outage].susceptance = 0.0
    isf_outage = injection_shift_factors(instance.lines)
    expected = isf_outage * inj

    pre_flow = pre_contingency_flow(inj, isf)
    actual = post_contingency_flow(inj, isf, lodf, pre_flow, outage)
    @test norm(expected - actual) < 1e-4
end

@testset "Compute violations" begin
    instance = load_standard_uc_instance("split/13b-2z")
    L = length(instance.lines)
    flow = [500.0 for l in 1:L]
    limits = [l.limit for l in instance.lines]
    filter = ViolationFilter(24, L, max_per_line=99, max_per_period=99)
    @test norm(flow - limits) < 1e-8

    flow[1] = 600.0
    flow[2] = -750.0
    flow[5] = 1000.0

    expected = [Violation(1, 1, -1, 100.0, 500.0),
                Violation(1, 2, -1, 250.0, 500.0),
                Violation(1, 5, -1, 500.0, 500.0)]
    find_violations(filter, flow, limits, time=1, scenario=-1)
    actual = query(filter)
    @test actual == expected
end

@testset "Violation filter" begin
    filter = ViolationFilter(2, 2, max_per_line=1, max_per_period=2)

    # Time 1
    offer(filter, Violation(1, 1, 1, 300., 100.))
    offer(filter, Violation(1, 1, 5, 500., 100.))
    offer(filter, Violation(1, 1, 4, 400., 100.))
    offer(filter, Violation(1, 2, 1, 200., 100.))
    offer(filter, Violation(1, 2, 8, 100., 100.))

    # Time 2
    offer(filter, Violation(2, 1, 1, 500., 100.))
    offer(filter, Violation(2, 1, 2, 400., 100.))
    offer(filter, Violation(2, 1, 3, 300., 100.))
    offer(filter, Violation(2, 2, 1, 900., 100.))

    actual = query(filter)
    expected = [Violation(1, 2, 1, 200., 100.),
                Violation(1, 1, 5, 500., 100.),
                Violation(2, 1, 1, 500., 100.),
                Violation(2, 2, 1, 900., 100.)]
    @test actual == expected
end
