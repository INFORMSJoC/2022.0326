# UnitCommitment.jl: Optimization Package for Security-Constrained Unit Commitment
# Copyright (C) 2020-2024, UChicago Argonne, LLC. All rights reserved.
# Released under the modified BSD license. See COPYING.md for more details.

using JuMP

function relax_bin_vars!(model::JuMP.Model)
    bin_vars = [v for v in all_variables(model) if is_binary(v)]
    unset_binary.(bin_vars)
    set_upper_bound.(bin_vars, 1.0)
    set_lower_bound.(bin_vars, 0.0)
end
