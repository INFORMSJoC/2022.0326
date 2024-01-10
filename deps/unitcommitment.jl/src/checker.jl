# UnitCommitment.jl: Optimization Package for Security-Constrained Unit Commitment
# Copyright (C) 2020-2024, UChicago Argonne, LLC. All rights reserved.
# Released under the modified BSD license. See COPYING.md for more details.

function verify_flows(inj::Array{Float64,1},
                      isf::Array{Float64,2},
                      lodf::Array{Float64,2},
                      lines::Array{TransmissionLine,1};
                      security::Bool = true)
                      
    L = length(lines)
    pre_flow = isf * inj
    for monitored in 1:L
        abs(pre_flow[monitored]) < lines[monitored].normal_limit + 0.1 ||
                @warn("Transmission violation: monitored=$monitored " *
                      "flow=$(pre_flow[monitored]) " *
                      "limit=$(lines[monitored].normal_limit)")
    end
    if security
        for outage in 1:L
            lines[outage].is_vulnerable || continue
            post_flow = pre_flow + lodf[:,outage] * pre_flow[outage]
            for monitored in 1:L
                abs(post_flow[monitored]) < lines[monitored].emergency_limit + 0.1 ||
                        @warn("Security violation: monitored=$monitored " *
                              "outage=$outage flow=$(post_flow[monitored]) " *
                              "limit=$(lines[monitored].emergency_limit)")
            end
        end
    end
end

function verify_power_balance(inj::Array{Float64,1}, time::Int)
    balance = sum(inj)
    if abs(balance) > 0.1
        @warn("Power balance violation: $balance MW at time $time")
    end
end
