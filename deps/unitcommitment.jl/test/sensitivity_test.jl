# UnitCommitment.jl: Optimization Package for Security-Constrained Unit Commitment
# Copyright (C) 2020-2024, UChicago Argonne, LLC. All rights reserved.
# Released under the modified BSD license. See COPYING.md for more details.

using UnitCommitment, Test, LinearAlgebra

@testset "Injection Shift Factors" begin
    instance = load_standard_uc_instance("split/13b-2z")
    isf = injection_shift_factors(instance.lines)
    expected = [0. -0.73 -0.27 -0.64 -0.54 -0.6  -0.58 -0.62 -0.61 -0.61 -0.61 -0.61 -0.61;
                0. -0.27 -0.73 -0.36 -0.46 -0.4  -0.42 -0.38 -0.39 -0.39 -0.39 -0.39 -0.39;
                0.  0.09 -0.09 -0.56 -0.18 -0.4  -0.31 -0.47 -0.44 -0.46 -0.45 -0.45 -0.45;
                0.  0.18 -0.18 -0.08 -0.37 -0.2  -0.27 -0.15 -0.17 -0.16 -0.16 -0.16 -0.16;
                0. -0.27  0.27 -0.36 -0.46 -0.4  -0.42 -0.38 -0.39 -0.39 -0.39 -0.39 -0.39;
                0.  0.04 -0.04  0.19 -0.08 -0.31 -0.   -0.03  0.07  0.02  0.04  0.04  0.04;
                0.  0.02 -0.02  0.11 -0.04 -0.12 -0.05 -0.42 -0.06 -0.24 -0.15 -0.15 -0.15;
                0.  0.03 -0.03  0.14 -0.06  0.02 -0.25 -0.02 -0.45 -0.24 -0.34 -0.34 -0.34;
                0. -0.06  0.06 -0.28  0.11 -0.51 -0.04 -0.36 -0.2  -0.28 -0.24 -0.24 -0.24;
                0. -0.03  0.03 -0.16  0.07 -0.09 -0.65 -0.17 -0.36 -0.27 -0.31 -0.31 -0.31;
                0. -0.02  0.02 -0.09  0.03  0.19 -0.05 -0.39 -0.13 -0.26 -0.19 -0.19 -0.19;
                0. -0.03  0.03 -0.16  0.07 -0.09  0.35 -0.17 -0.36 -0.27 -0.31 -0.31 -0.31;
                0.  0.   -0.    0.02 -0.01  0.07 -0.1   0.2  -0.19 -0.5  -0.35 -0.35 -0.35;
                0. -0.    0.   -0.01  0.   -0.03  0.05 -0.1   0.1  -0.25 -0.49 -0.33 -0.16;
                0. -0.    0.   -0.01  0.   -0.03  0.05 -0.1   0.1  -0.25 -0.16 -0.33 -0.49;
                0.  0.   -0.    0.01 -0.    0.03 -0.05  0.1  -0.1   0.25 -0.34 -0.17 -0.01;
                0.  0.   -0.    0.01 -0.    0.03 -0.05  0.1  -0.1   0.25 -0.01 -0.17 -0.34;
                0. -0.   -0.   -0.   -0.    0.   -0.    0.    0.   -0.    0.17 -0.5  -0.17;
                0.  0.    0.    0.    0.    0.    0.    0.    0.    0.    0.17  0.5  -0.17]
    @test norm(round.(isf, digits=2) - expected) ≈ 0.

    change_slack!(isf, 1, 6)
    expected = [  0.6  -0.13  0.33 -0.04  0.06 0.  0.02 -0.02 -0.01 -0.01 -0.01 -0.01 -0.01;
                  0.4   0.13 -0.33  0.04 -0.06 0. -0.02  0.02  0.01  0.01  0.01  0.01  0.01;
                  0.4   0.49  0.31 -0.15  0.22 0.  0.09 -0.07 -0.04 -0.05 -0.05 -0.05 -0.05;
                  0.2   0.38  0.02  0.12 -0.17 0. -0.07  0.05  0.03  0.04  0.03  0.03  0.03;
                  0.4   0.13  0.67  0.04 -0.06 0. -0.02  0.02  0.01  0.01  0.01  0.01  0.01;
                  0.31  0.34  0.27  0.5   0.23 0.  0.3   0.27  0.37  0.32  0.35  0.35  0.35;
                  0.12  0.14  0.1   0.22  0.08 0.  0.07 -0.3   0.06 -0.12 -0.03 -0.03 -0.03;
                 -0.02  0.01 -0.05  0.12 -0.08 0. -0.27 -0.04 -0.47 -0.26 -0.36 -0.36 -0.36;
                  0.51  0.45  0.56  0.23  0.62 0.  0.46  0.15  0.31  0.23  0.27  0.27  0.27;
                  0.09  0.06  0.12 -0.07  0.16 0. -0.56 -0.08 -0.27 -0.18 -0.22 -0.22 -0.22;
                 -0.19 -0.2  -0.17 -0.27 -0.15 0. -0.24 -0.57 -0.32 -0.45 -0.38 -0.38 -0.38;
                  0.09  0.06  0.12 -0.07  0.16 0.  0.44 -0.08 -0.27 -0.18 -0.22 -0.22 -0.22;
                 -0.07 -0.07 -0.07 -0.05 -0.08 0. -0.17  0.13 -0.26 -0.57 -0.42 -0.42 -0.42;
                  0.03  0.03  0.04  0.03  0.04 0.  0.08 -0.06  0.13 -0.22 -0.46 -0.29 -0.13;
                  0.03  0.03  0.04  0.03  0.04 0.  0.08 -0.06  0.13 -0.22 -0.13 -0.29 -0.46;
                 -0.03 -0.03 -0.04 -0.03 -0.04 0. -0.08  0.06 -0.13  0.22 -0.37 -0.21 -0.04;
                 -0.03 -0.03 -0.04 -0.03 -0.04 0. -0.08  0.06 -0.13  0.22 -0.04 -0.21 -0.37;
                 -0.   -0.   -0.   -0.   -0.   0. -0.    0.    0.   -0.    0.17 -0.5  -0.17;
                 -0.    0.    0.    0.    0.   0.  0.    0.    0.    0.    0.17  0.5  -0.17]
    @test norm(round.(isf, digits=2) - expected) ≈ 0.
end

@testset "Line Outage Distribution Factors" begin
    instance = load_standard_uc_instance("split/13b-2z")
    isf_before = injection_shift_factors(instance.lines)
    lodf = line_outage_factors(instance.lines, isf_before)

    for outage in 1:length(instance.lines)
        prev_susceptance = instance.lines[outage].susceptance
        instance.lines[outage].susceptance = 0.0
        isf_after = injection_shift_factors(instance.lines)
        instance.lines[outage].susceptance = prev_susceptance

        for l in 1:length(instance.lines)
            expected = isf_after[l, :]
            actual = isf_before[l, :] + lodf[l, outage] * isf_before[outage, :]
            @test norm(expected - actual) < 1e-6
        end
    end
end

@testset "Post-Contingency ISF" begin
    instance = load_standard_uc_instance("split/13b-2z")
    isf_before = injection_shift_factors(instance.lines)
    lodf = line_outage_factors(instance.lines, isf_before)

    for outage in 1:length(instance.lines)
        prev_susceptance = instance.lines[outage].susceptance
        instance.lines[outage].susceptance = 0.0
        expected = injection_shift_factors(instance.lines)
        instance.lines[outage].susceptance = prev_susceptance

        actual = post_contingency_isf(isf_before, lodf, outage)
        @test norm(expected - actual) < 1e-6
    end
end
