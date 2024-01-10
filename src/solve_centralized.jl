# UnitCommitment.jl: Optimization Package for Security-Constrained Unit Commitment
# Copyright (C) 2020-2024, UChicago Argonne, LLC. All rights reserved.
# Released under the modified BSD license. See COPYING.md for more details.

using CPLEX, JuMP, UnitCommitment, LinearAlgebra, Printf

function solve_uc_centralized(name::String;
                             demand_scale=1.0,
                             limit_scale=1.0,
                             T=24,
                             transmission=true,
                             security=true,
                             relax=false
                            )
    
    instance, model, solution, stats = UnitCommitment.solve(name,
                                                            demand_scale = demand_scale,
                                                            limit_scale = limit_scale,
                                                            hours = T,
                                                            transmission = transmission,
                                                            security = security,
                                                            relax = relax)

    @info @sprintf("%s,centralized,%.2f,%.2f,%d,%d,%.6e,%.6e,%d,%.2f,%.2f",
                   instance.name, demand_scale, limit_scale, transmission, security, solution.cost,
                   0, 1, stats.wallclock_time, stats.wallclock_time)

    return solution.cost
end

export solve_uc_centralized
