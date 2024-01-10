# UnitCommitment.jl: Optimization Package for Security-Constrained Unit Commitment
# Copyright (C) 2020-2024, UChicago Argonne, LLC. All rights reserved.
# Released under the modified BSD license. See COPYING.md for more details.

using Dates, UnitCommitment, JuMP, LinearAlgebra, Base.Threads, Printf,
      TimerOutputs, JLD, MPI

struct Zone
    index::Int
    buses_int::Array{Bus, 1}
    buses_bnd::Array{Bus, 1}
    buses_ext::Array{Bus, 1}

    BI::Array{Int,1}
    BIN::Array{Int,1}
    BN::Array{Int,1}
    BNE::Array{Int,1}
    BE::Array{Int,1}
    is_neighbor_zone::Array{Bool,1}

    lines_int::Array{TransmissionLine, 1}
    lines_ext::Array{TransmissionLine, 1}
    generators::Array{Generator, 1}
    isf_int::Array{Float64, 2}
    isf_bnd::Array{Float64, 2}
    isf_ext::Array{Float64, 2}
    link_base::Array{Float64, 2}
    link_outage
end

Base.show(io::IO, x::Zone) = print(io, "zone$(x.index)")

function find_zones!(instance::UnitCommitmentInstance; ϵ::Float64=0.1,
                     gap::Float64=0.5, max_size::Int=10000000)
    find_zones!(instance.buses, instance.lines, instance.generators, ϵ=ϵ,
                gap=gap, max_size=max_size)
end

function find_split(buses::Array{Bus,1},
                    lines::Array{TransmissionLine,1},
                    generators::Array{Generator,1},
                    force_internal::Dict{Bus,Bool};
                    ϵ::Float64, gap::Float64)

    lines_at = Dict(b => [] for b in buses)
    for line in lines
        append!(lines_at[line.source], [line])
        append!(lines_at[line.target], [line])
    end
    generators_int = [g for g in generators if g.bus in buses]
    G = length(generators_int)
    has_generator = Dict(b => false for b in buses)
    for g in generators_int
        has_generator[g.bus] = true
    end

    optimizer = CPLEX.Optimizer(CPX_PARAM_SCRIND=1,
                                CPX_PARAM_THREADS=8,
                                CPX_PARAM_EPGAP=gap)

    model = direct_model(optimizer)
    cpx = backend(model).inner
    cpx.has_int = true

    @variable(model, is_bnd_bus[buses], Bin)
    @variable(model, is_int_bus[buses], Bin)
    @variable(model, is_internal_line[lines], Bin)

    # Minimize number of boundary buses
    @objective(model, Min, sum(is_bnd_bus))

    # Internal buses cannot be adjacent to both zones at same time
    for b in buses
        for l1 in lines_at[b], l2 in lines_at[b]
            if l1 == l2 continue end
            @constraint(model, is_internal_line[l1] + (1 - is_internal_line[l2]) <= 1 + is_bnd_bus[b])
            @constraint(model, (1 - is_internal_line[l1]) + is_internal_line[l2] <= 1 + is_bnd_bus[b])
        end
    end

    # If a bus is incident to an internal (resp. external) line, then it is
    # either internal (resp. external) or boundary.
    for b in buses
        for l in lines_at[b]
            @constraint(model, is_int_bus[b] + is_bnd_bus[b] >= is_internal_line[l])
            @constraint(model, (1 - is_int_bus[b]) + is_bnd_bus[b] >= (1 - is_internal_line[l]))
        end
    end

    # Partition needs to be somewhat balanced
    @constraint(model, sum(is_internal_line) <= length(lines) * (0.5 + ϵ))
    @constraint(model, sum(is_internal_line) >= length(lines) * (0.5 - ϵ))

    for b in buses
        if force_internal[b]
            @constraint(model, is_int_bus[b] == 1)
            @constraint(model, is_bnd_bus[b] == 0)
        end
    end

    # # Number of generators on each side must be somewhat balanced
    # @constraint(model, sum(is_int_bus[b] for b in buses if has_generator[b]) <= G * (0.5 + ϵ))
    # @constraint(model, sum(is_int_bus[b] for b in buses if has_generator[b]) >= G * (0.5 - ϵ))

    # Boundary buses cannot have generators
    for g in generators
        if g.bus in buses
            @constraint(model, is_bnd_bus[g.bus] == 0)
        end
    end

    optimize!(model)
    return round.(value.(is_internal_line)), round.(value.(is_bnd_bus))
end

function find_zones!(buses::Array{Bus,1},
                     all_lines::Array{TransmissionLine,1},
                     generators::Array{Generator,1};
                     ϵ::Float64, gap::Float64, max_size::Int)

    stack = [all_lines]
    next_zone_number = 1
    force_internal = Dict(b => false for b in buses)

    while !isempty(stack)
        lines = pop!(stack)
        next_zone_number += 1

        is_internal, is_boundary = find_split(buses, lines, generators, force_internal, ϵ=ϵ, gap=gap)
        lines_int = [l for l in lines if is_internal[l] == 1]
        lines_ext = [l for l in lines if is_internal[l] == 0]
        for line in lines_ext line.zone = next_zone_number end

        for bus in buses
            if is_boundary[bus] == 1
                bus.demand .= 0
                force_internal[bus] = true
            end
        end

        if length(lines_int) > max_size push!(stack, lines_int) end
        if length(lines_ext) > max_size push!(stack, lines_ext) end
    end
end


function compute_link_outage(isf, lodf, outage, L, BB, BEE)
    n_lines, n_buses = size(isf)
    isf_outage = [isf[l,b] + lodf[l, outage] * isf[outage, b]
                  for l in 1:n_lines, b in 1:n_buses]
    ignored, lo, ignored2 = LAPACK.gels!('N', isf_outage[L,BB], isf_outage[L,BEE])
    return lo
end

function extract_zone(instance::UnitCommitmentInstance,
                      zone::Int,
                      isf::Array{Float64, 2},
                      lodf::Array{Float64, 2};
                      security=true) :: Zone

    Z = maximum([line.zone for line in instance.lines])

    # Set of transmission lines incident to bus b
    lines_at = Dict(b => [] for b in instance.buses)
    for line in instance.lines
      append!(lines_at[line.source], [line])
      append!(lines_at[line.target], [line])
    end

    # Is bus b located in zone z?
    is_bus_in_zone = Dict((b,z) => false for b in instance.buses, z in 1:Z)
    for b in instance.buses
        for line in lines_at[b]
            is_bus_in_zone[b, line.zone] = true
        end
    end

    # Number of zones each bus belongs to
    num_bus_zones = Dict(b => sum(is_bus_in_zone[b,z] for z in 1:Z)
                         for b in instance.buses)

    BI  = [b.index for b in instance.buses if is_bus_in_zone[b, zone] && num_bus_zones[b] == 1]
    BIN = [b.index for b in instance.buses if is_bus_in_zone[b, zone] && num_bus_zones[b] > 1]

    # Is a given zone neighbor to the current one? That is, do they share a bus?
    is_neighbor_zone = [false for z in 1:Z]
    for b in BIN
        for line in lines_at[instance.buses[b]]
            if line.zone == zone continue end
            is_neighbor_zone[line.zone] = true
        end
    end

    is_bus_in_neighbor_zone = Dict(b => false for b in instance.buses)
    for z in 1:Z
        if !is_neighbor_zone[z] continue end
        for b in instance.buses
            if is_bus_in_zone[b,z]
                is_bus_in_neighbor_zone[b] = true
            end
        end
    end

    BN  = [b.index for b in instance.buses if is_bus_in_neighbor_zone[b] && num_bus_zones[b] == 1]
    BNE = [b.index for b in instance.buses if !is_bus_in_zone[b, zone] && is_bus_in_neighbor_zone[b] && num_bus_zones[b] > 1]
    BE  = [b.index for b in instance.buses if !is_bus_in_zone[b, zone] && !is_bus_in_neighbor_zone[b]]

    buses_int = [instance.buses[b] for b in BI]
    buses_bnd = [instance.buses[b] for b in BIN]
    buses_ext = [b for b in instance.buses if !is_bus_in_zone[b, zone]]

    lines_int = [l for l in instance.lines if l.zone == zone]
    lines_ext = [l for l in instance.lines if l.zone != zone]

    generators = [g for g in instance.generators if g.bus in buses_int]

    L = [l.index for l in lines_int]
    isf_base = copy(isf)
    change_slack!(isf_base, 1, buses_int[1].index)
    BI = [b.index for b in buses_int]
    BB = [b.index for b in buses_bnd]
    BEE = [b.index for b in buses_ext]
    ignored, link_base, ignored2 = LAPACK.gels!('N', isf_base[L, BB], isf_base[L, BEE])

    link_outage::Dict{Int64, Array{Float64,2}} = Dict(line.index => Array{Float64}(undef,0,0) for line in lines_ext)
    rank = MPI.Comm_rank(MPI.COMM_WORLD) + 1
    if security
        for outage in lines_ext
            cache_filename = @sprintf("cache/%s/%d/%06d.jld", instance.name, zone, outage.index)
            if isfile(cache_filename)
                link_outage[outage.index] = load(cache_filename, "link_outage")
            else
                link_outage[outage.index] = compute_link_outage(isf_base, lodf, outage.index, L, BB, BEE)
                if length(lines_ext) > 100 && rank == 1
                    mkpath("cache/$(instance.name)/$zone/")
                    save(cache_filename, "link_outage", link_outage[outage.index])
                end
            end
        end
    end

    return Zone(zone,
                buses_int, buses_bnd, buses_ext,
                BI, BIN, BN, BNE, BE, is_neighbor_zone,
                lines_int, lines_ext,
                generators, isf_base[L, BI], isf_base[L, BB], isf_base[L, BEE],
                link_base, link_outage)
end

function zone_topology(instance::UnitCommitmentInstance)
    Z = maximum([line.zone for line in instance.lines])
    adj = zeros(Z,Z)

    lines_at = Dict(b => [] for b in instance.buses)
    for line in instance.lines
        append!(lines_at[line.source], [line])
        append!(lines_at[line.target], [line])
    end

    for b in instance.buses
        for l1 in lines_at[b], l2 in lines_at[b]
            if l1 == l2 continue end
            adj[l1.zone, l2.zone] = 1.0
        end
    end

    for z1 in 1:(Z-1)
        for z2 in (z1+1):Z
            if adj[z1,z2] > 0
                @info z1, z2
            end
        end
    end
end

export Zone, find_zones!, extract_zone, zone_topology
