# UnitCommitment.jl: Optimization Package for Security-Constrained Unit Commitment
# Copyright (C) 2020-2024, UChicago Argonne, LLC. All rights reserved.
# Released under the modified BSD license. See COPYING.md for more details.

using JuMP, MathOptInterface

struct Violation
    time::Int
    monitored_line::Int
    outage_line::Int
    amount::Float64
    limit::Float64
end

mutable struct UnitCommitmentModel
    mip::JuMP.Model
    prod_vars::JuMP.Containers.DenseAxisArray
    inj_vars::JuMP.Containers.DenseAxisArray
    is_on_vars::JuMP.Containers.DenseAxisArray
    switch_on_vars::JuMP.Containers.DenseAxisArray
    switch_off_vars::JuMP.Containers.DenseAxisArray
    reserve_vars::JuMP.Containers.DenseAxisArray
    obj_var::JuMP.VariableRef
    power_balance_eq::JuMP.Containers.DenseAxisArray
    inj_eq::JuMP.Containers.DenseAxisArray
    violations::Array{Violation,1}
    hours::Int
end

function uc_model(generators::Array{Generator,1},
                  buses::Array{Bus,1},
                  optimizer;
                  T::Int = 24,
                  reserve_pct::Float64 = 0.0,
                  named::Bool = true) :: UnitCommitmentModel
    time = 1:T
    segments = 1:3

    if isa(optimizer, JuMP.OptimizerFactory)
        mip = Model(optimizer)
    else
        mip = direct_model(optimizer)
    end

    @variable(mip, segprod[generators, segments, time] >= 0)
    @variable(mip, prod[[g.index for g in generators], time] >= 0)
    @variable(mip, reserve[[g.index for g in generators], time] >= 0)
    @variable(mip, inj[[b.index for b in buses], time])
    @variable(mip, is_on[generators, time], Bin)
    @variable(mip, switch_on[generators, time], Bin)
    @variable(mip, switch_off[generators, time], Bin)
    @variable(mip, original_obj)

    # Original objective function
    @objective(mip, Min, original_obj)
    @constraint(mip, original_obj ==
        # Startup costs
        + sum(switch_on[g, t] * g.startup_cost for g in generators, t in time)

        # Minimum production costs
        + sum(is_on[g, t] * g.cost_min_power for g in generators, t in time)

        # Piecewise-linear production costs
        + sum(segprod[g, k, t] * g.cost_segment[k] for g in generators, k in segments, t in time))

    # Definition of net injection
    @constraint(mip,
                inj_eq[b in buses, t in time],
                inj[b.index, t] == - b.demand[t] +
                        sum(prod[g.index, t] for g in generators if g.bus == b))

    # Definition of total production
    @constraint(mip,
                [g in generators, t in time],
                prod[g.index, t] == is_on[g, t] * g.min_power + sum(segprod[g, k, t] for k in segments))

    # Production limits
    @constraint(mip,
                [g in generators, t in time, k in segments],
                segprod[g, k, t] <= g.offer_segment[k])

    for g in generators
        if g.minimum_uptime > 1
            @constraint(mip,
                        [t in 1:T-1],
                        prod[g.index, t] + reserve[g.index, t] <=
                                g.max_power * is_on[g, t]
                                - (g.max_power - g.ramp_up_limit) * switch_on[g, t]
                                - (g.max_power - g.ramp_down_limit) * switch_off[g, t + 1])
        else
            @constraint(mip,
                        [t in 1:T-1],
                        prod[g.index, t] + reserve[g.index, t] <=
                                g.max_power * is_on[g, t]
                                - (g.max_power - g.ramp_up_limit) * switch_on[g, t]
                                - (max(g.ramp_up_limit - g.ramp_down_limit, 0)) * switch_off[g, t + 1])
            @constraint(mip,
                        [t in 1:T-1],
                        prod[g.index, t] + reserve[g.index, t] <=
                                g.max_power * is_on[g, t]
                                - max(g.ramp_down_limit - g.ramp_up_limit, 0) * switch_on[g, t]
                                - (g.max_power - g.ramp_down_limit) * switch_off[g, t + 1])
        end
        @constraint(mip,
                    prod[g.index, T] + reserve[g.index, T] <=
                            g.max_power * is_on[g, T]
                            - (g.max_power - g.ramp_up_limit) * switch_on[g, T])
    end

    # Link binary variables
    @constraint(mip,
                [g in generators, t in 2:T],
                is_on[g, t] - is_on[g, t-1] == switch_on[g, t] - switch_off[g, t])

    # Avoid switching on and off at the same time
    @constraint(mip,
                [g in generators, t in time],
                switch_on[g, t] + switch_off[g, t] <= 1)

    # Power balance equations
    @constraint(mip,
                power_balance_eq[t in time],
                sum(inj[b.index, t] for b in buses) == 0)

    # Ramp up and down limits
    @constraint(mip,
                [g in generators, t in 2:T],
                prod[g.index, t] <= prod[g.index, t - 1] + g.ramp_up_limit)

    @constraint(mip,
                [g in generators, t in 2:T],
                prod[g.index, t] >= prod[g.index, t - 1] - g.ramp_down_limit)

    # Minimum up-time
    @constraint(mip,
                [g in generators, t in time],
                sum(switch_on[g, i] for i in (t - g.minimum_uptime + 1):t if i >= 1) <=
                        is_on[g, t])

    # Minimum down-time
    @constraint(mip,
                [g in generators, t in time],
                sum(switch_on[g, i] for i in (t - g.minimum_downtime + 1):t if i >= 1) <=
                        1 - is_on[g, max(1, t - g.minimum_downtime)])

    # Minimum spinning reserve
    @constraint(mip,
                [t in time],
                sum(reserve[g.index, t] for g in generators) >= reserve_pct * sum(b.demand[t] for b in buses))

    return UnitCommitmentModel(mip,
                               prod,
                               inj,
                               is_on,
                               switch_on,
                               switch_off,
                               reserve,
                               original_obj,
                               power_balance_eq,
                               inj_eq,
                               Array{Violation,1}(undef,0),
                               T)
end

export UnitCommitmentModel, uc_model, Violation
