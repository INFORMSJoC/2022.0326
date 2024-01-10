# UnitCommitment.jl: Optimization Package for Security-Constrained Unit Commitment
# Copyright (C) 2020-2024, UChicago Argonne, LLC. All rights reserved.
# Released under the modified BSD license. See COPYING.md for more details.

using DataStructures

mutable struct ViolationFilter
    L::Int
    max_per_line::Int
    max_per_period::Int
    queues::Dict{Int,PriorityQueue{Violation,Float64}}
end

# DEPRECATED
function pre_contingency_flow(inj::Array{Float64,1},
                              isf::Array{Float64,2},
                             )::Array{Float64,1}
    return isf * inj
end

# DEPRECATED
function post_contingency_flow(inj::Array{Float64,1},
                               isf::Array{Float64,2},
                               lodf::Array{Float64,2},
                               pre_flow::Array{Float64,1},
                               outage::Int,
                              )::Array{Float64,1}
    return pre_flow + lodf[:, outage] * pre_flow[outage]
end

# DEPRECATED
function find_violations(filter::ViolationFilter,
                         flow::Array{Float64,1},
                         limits::Array{Float64,1};
                         time=1,
                         outage_line=-1)
    L = length(limits)
    vt = abs.(max.(zeros(L), flow - limits, - flow - limits))
    for line in 1:L
        if vt[line] < 1e-5 continue end
        offer(filter, Violation(time, line, outage_line < 0 ? line : outage_line, vt[line], limits[line]))
    end
end


function ViolationFilter(L::Int;
                         max_per_line::Int=1,
                         max_per_period::Int=5,
                        )::ViolationFilter
    return ViolationFilter(L, max_per_line, max_per_period,
                           Dict(l => PriorityQueue{Violation,Float64}()
                                for l in 1:L))
end

function offer(filter::ViolationFilter, v::Violation)
    q = filter.queues[v.monitored_line]
    if length(q) < filter.max_per_line
        enqueue!(q, v, v.amount)
    else
        if v.amount > DataStructures.peek(q)[1].amount
            dequeue!(q)
            enqueue!(q, v, v.amount)
        end
    end
end

function query(filter::ViolationFilter)
    violations = Array{Violation,1}()
    time_queue = PriorityQueue{Violation,Float64}()
    for l in 1:filter.L
        line_queue = filter.queues[l]
        while length(line_queue) > 0
            v = dequeue!(line_queue)
            if length(time_queue) < filter.max_per_period
                enqueue!(time_queue, v, v.amount)
            else
                if v.amount > DataStructures.peek(time_queue)[1].amount
                    dequeue!(time_queue)
                    enqueue!(time_queue, v, v.amount)
                end
            end
        end
    end
    while length(time_queue) > 0
        violations = [violations; dequeue!(time_queue)]
    end
    return violations
end

export ViolationFilter, pre_contingency_flow, post_contingency_flow, find_violations, offer, query
