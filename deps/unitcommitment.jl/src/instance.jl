# UnitCommitment.jl: Optimization Package for Security-Constrained Unit Commitment
# Copyright (C) 2020-2024, UChicago Argonne, LLC. All rights reserved.
# Released under the modified BSD license. See COPYING.md for more details.

using Printf

mutable struct Bus
    index::Int
    demand::Array{Float64,1}
    zone::Int
end

mutable struct TransmissionLine
    index::Int
    source::Bus
    target::Bus
    reactance::Float64
    susceptance::Float64
    normal_limit::Float64
    emergency_limit::Float64
    is_vulnerable::Bool
    zone::Int
end

mutable struct Generator
    index::Int
    min_power::Float64
    max_power::Float64
    ramp_down_limit::Float64
    ramp_up_limit::Float64
    shutdown_ramp_limit::Float64
    startup_ramp_limit::Float64
    initial_state::Int
    bus::Bus
    is_always_on::Bool
    minimum_uptime::Int
    minimum_downtime::Int
    cost_min_power::Float64
    cost_segment::Array{Float64,1}
    offer_segment::Array{Float64,1}
    startup_cost::Float64
end

mutable struct UnitCommitmentInstance
    name::String
    variation::String
    lines::Array{TransmissionLine,1}
    buses::Array{Bus,1}
    generators::Array{Generator, 1}
end

Base.show(io::IO, x::TransmissionLine) = print(io, "line($(x.source),$(x.target))")
Base.show(io::IO, x::Generator) = print(io, "gen$(x.index)")
Base.show(io::IO, x::Bus) = print(io, "bus$(x.index)")

function load_standard_uc_instance(name::String) :: UnitCommitmentInstance
    basedir = dirname(@__FILE__)
    return read_uc_instance("$basedir/../instances/$name")
end

function read_uc_instance(path::String)::UnitCommitmentInstance
    path = replace(path, r"/$" => "")
    name = replace(path, r"[^/]*/" => "")

    csv = read_csv("$(path)/buses.csv").data
    n_buses = size(csv)[1]
    buses = [Bus(parse(Int, csv[b, 1]),
                 [parse(Float64, csv[b, t+1]) for t in 1:24],
                 parse(Int, csv[b, 26]))
             for b in 1:n_buses]

    csv = read_csv("$(path)/lines.csv").data
    n_lines = size(csv)[1]
    n_lines, n_cols = size(csv)
    if n_cols == 7
        lines = [TransmissionLine(parse(Int, csv[l, 1]), # index
                                  buses[parse(Int, csv[l, 2])], # source
                                  buses[parse(Int, csv[l, 3])], # target
                                  parse(Float64, csv[l, 4]), # reactance
                                  100 * (pi / 180) / parse(Float64, csv[l, 4]), # susceptance
                                  parse(Float64, csv[l, 5]), # normal_limit
                                  parse(Float64, csv[l, 5]), # emergency_limits
                                  parse(Bool, csv[l, 6]), # is_vulnerable
                                  parse(Int, csv[l, 7])) # zone
                 for l in 1:n_lines]
    elseif n_cols == 8
        lines = [TransmissionLine(parse(Int, csv[l, 1]), # index 
                                  buses[parse(Int, csv[l, 2])], # source
                                  buses[parse(Int, csv[l, 3])], # target
                                  parse(Float64, csv[l, 4]), # reactance
                                  100 * (pi / 180) / parse(Float64, csv[l, 4]), # susceptance
                                  parse(Float64, csv[l, 5]), # normal_limit
                                  parse(Float64, csv[l, 6]), # emergency_limit
                                  parse(Bool, csv[l, 7]), # is_vulnerable
                                  parse(Int, csv[l, 8])) # zone
                 for l in 1:n_lines]
    else
        @error "wrong number of columns in lines.csv: $nlines"
    end

    for l in lines
        if l.zone == 0
            l.zone = l.source.zone
        end
    end

    csv = read_csv("$(path)/generators.csv").data
    n_generators = size(csv)[1]
    generators = [Generator(parse(Int, csv[g, 1]),
                            parse(Float64, csv[g, 2]),
                            parse(Float64, csv[g, 3]),
                            parse(Float64, csv[g, 4]),
                            parse(Float64, csv[g, 5]),
                            parse(Float64, csv[g, 6]),
                            parse(Float64, csv[g, 7]),
                            parse(Int, csv[g, 8]),
                            buses[parse(Int, csv[g, 9])],
                            parse(Bool, csv[g, 10]),
                            parse(Int, csv[g, 11]),
                            parse(Int, csv[g, 12]),
                            parse(Float64, csv[g, 13]),
                            [parse(Float64, csv[g, 14+k]) for k in 0:2],
                            [parse(Float64, csv[g, 17+k]) for k in 0:2],
                            parse(Float64, csv[g, 20]))
                  for g in 1:n_generators]

    for g in generators
        g.max_power = g.min_power + sum(g.offer_segment[k] for k in 1:3)
    end
    variation = ""
    return UnitCommitmentInstance(name, variation, lines, buses, generators)
end

function write_uc_instance(instance::UnitCommitmentInstance, path::String)
    open("$path/generators.csv", "w") do file
        write(file, "Unit,Min Power (MW),Max Power (MW),Ramp-Down Limit (MW)," *
                    "Ramp-Up Limit (MW),Shutdown Ramp Limit (MW)," *
                    "Startup Ramp Limit (MW),Initial State,Bus,Always On," *
                    "Minimum Uptime (h),Minimum Downtime (h),Price Min Power (\$)," *
                    "Price Segment 1 (\$/MW),Price Segment 2 (\$/MW)," *
                    "Price Segment 3 (\$/MW),Offer Segment 1 (MW)," *
                    "Offer Segment 2 (MW),Offer Segment 3 (MW),Startup Cost (\$)\n")
        for g in instance.generators
            write(file, @sprintf("%d,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%d,%d,%d,%d,%d,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f\n",
                                 g.index, g.min_power, g.max_power, g.ramp_down_limit, g.ramp_up_limit, g.shutdown_ramp_limit,
                                 g.startup_ramp_limit, g.initial_state, g.bus.index, g.is_always_on, g.minimum_uptime, g.minimum_downtime,
                                 g.cost_min_power, g.cost_segment[1], g.cost_segment[2], g.cost_segment[3],
                                 g.offer_segment[1], g.offer_segment[2], g.offer_segment[3], g.startup_cost))
        end
    end

    open("$path/lines.csv", "w") do file
        write(file, "Line,Source,Target,Reactance,Normal Flow Limit (MW),Emergency Flow Limit (MW),Vulnerable?,Zone\n")
        for l in instance.lines
            write(file, @sprintf("%d,%d,%d,%.8f,%.0f,%.0f,%d,%d\n",
                                 l.index, l.source.index, l.target.index,
                                 l.reactance, l.normal_limit, l.emergency_limit,
                                 l.is_vulnerable, l.zone))
        end
    end

    open("$path/buses.csv", "w") do file
        write(file, "Bus," * join(["Demand $t (MW)," for t in 1:24]) * "Zone\n")
        for bus in instance.buses
            d = join([@sprintf("%.8f,", dt) for dt in bus.demand])
            write(file, @sprintf("%d,%s%d\n", bus.index, d, bus.zone))
        end
    end
end

function print_summary(instance::UnitCommitmentInstance)
    n_buses = size(instance.buses)[1]
    n_generators = size(instance.generators)[1]
    n_lines_total = size(instance.lines)[1]
    @info @sprintf("%12d units", n_generators)
    @info @sprintf("%12d buses", n_buses)
    @info @sprintf("%12d lines", n_lines_total)
end

function scale_demands!(instance::UnitCommitmentInstance, factor::Float64)
    for b in instance.buses
        b.demand *= factor
    end
end

function scale_limits!(instance::UnitCommitmentInstance, factor::Float64)
    for l in instance.lines
        l.normal_limit *= factor
        l.emergency_limit *= factor
    end
end

function scale_startup_costs!(instance::UnitCommitmentInstance, factor::Float64)
    for g in instance.generators
        g.startup_cost *= factor
    end
end

function installed_capacity(instance::UnitCommitmentInstance)::Float64
    return sum([g.max_power for g in instance.generators])
end

function system_demand(instance::UnitCommitmentInstance; time::Int)::Float64
    return sum([b.demand[time] for b in instance.buses])
end

function peak_demand(instance::UnitCommitmentInstance)::Float64
    T = length(instance.buses[1].demand)
    return maximum([system_demand(instance, time=t) for t in 1:T])
end

export load_standard_uc_instance,
       installed_capacity,
       peak_demand,
       read_uc_instance,
       scale_demands!,
       scale_limits!,
       scale_startup_costs!,
       system_demand,
       write_uc_instance,
       Bus,
       Generator,
       TransmissionLine,
       UnitCommitmentInstance
