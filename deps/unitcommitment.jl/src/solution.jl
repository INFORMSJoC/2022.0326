# UnitCommitment.jl: Optimization Package for Security-Constrained Unit Commitment
# Copyright (C) 2020-2024, UChicago Argonne, LLC. All rights reserved.
# Released under the modified BSD license. See COPYING.md for more details.

struct UnitCommitmentSolution
    cost::Float64
    is_on::Array{Bool,2}             # is_on[generator, time]
    production::Array{Float64,2}     # production[generator, time]
    reserve::Array{Float64,2}        # reserve[generator, time]
    injection::Array{Float64,2}      # injection[bus, time]
    violations::Array{Violation, 1}
end


function get_solution(instance::UnitCommitmentInstance,
                      model::UnitCommitmentModel,
                     )::UnitCommitmentSolution
    is_on_values = value.(model.is_on_vars)
    prod_values = value.(model.prod_vars)
    reserve_values = value.(model.reserve_vars)
    inj_values = value.(model.inj_vars)
    cost = value(model.obj_var)
    n_buses = length(instance.buses)

    is_on = [round(is_on_values[g,t]) for g in instance.generators, t in 1:model.hours]
    production = [prod_values[g.index,t] for g in instance.generators, t in 1:model.hours]
    reserve = [reserve_values[g.index,t] for g in instance.generators, t in 1:model.hours]
    injection = [inj_values[b,t] for b in 1:n_buses, t in 1:model.hours]

    return UnitCommitmentSolution(cost,
                                  is_on,
                                  production,
                                  reserve,
                                  injection,
                                  model.violations)
end


function write_csv(path::String,
	               instance::UnitCommitmentInstance,
	               model::UnitCommitmentModel,
	               sol::UnitCommitmentSolution
	              ) :: Nothing

    should_print_header = !isfile(path)
    io = open(path, "a")

    G = length(instance.generators)
    B = length(instance.buses)
    T = model.hours

    if should_print_header
    	print(io, "instance")
    	print(io, ",variation")
    	print(io, ",cost")
    	for g in 1:G, t in 1:T
    		print(io, ",is_on[$g:$t]")
    	end
        for g in 1:G, t in 1:T
            print(io, ",prod[$g:$t]")
        end
        for g in 1:G, t in 1:T
            print(io, ",reserve[$g:$t]")
        end
    	for b in 1:B, t in 1:T
    		print(io, ",inj[$b:$t]")
    	end
    	print(io, ",violations")
    	println(io)
    end

    print(io, instance.name)
    print(io, ",$(instance.variation)")
    print(io, @sprintf(",%.8f", sol.cost))
	for g in 1:G, t in 1:T
		print(io, @sprintf(",%d", sol.is_on[g,t] ? 1 : 0))
	end
    for g in 1:G, t in 1:T
        print(io, @sprintf(",%.8f", sol.production[g,t]))
    end
    for g in 1:G, t in 1:T
        print(io, @sprintf(",%.8f", sol.reserve[g,t]))
    end
	for b in 1:B, t in 1:T
		print(io, @sprintf(",%.8f", sol.injection[b,t]))
	end

	print(io, ",")
	for v in sol.violations
        v.time == 1 || continue
        print(io, @sprintf("%d:%d ", v.monitored_line, v.outage_line))
    end

    println(io)

    close(io)
    return
end


function read_uc_solution_csv(path::String,
                              instance::UnitCommitmentInstance;
                              T::Int = 24,
                             )::Array{UnitCommitmentSolution,1}
    
    csv_data = read_csv(path).data
    n_rows, n_cols = size(csv_data)
    G = length(instance.generators)
    B = length(instance.buses)

    results = Array{UnitCommitmentSolution}(undef,0)
    for row in 1:n_rows
        is_on = zeros(Bool, G, T)
        production = zeros(G, T)
        reserve = zeros(G, T)
        injection = zeros(B, T)
        violations = Array{Violation,1}(undef,0)

        col = 3
        cost = parse(Float64, csv_data[row, col])
        col += 1
        for g in 1:G, t in 1:T
            is_on[g,t] = parse(Bool, csv_data[row, col])
            col += 1
        end
        for g in 1:G, t in 1:T
            production[g,t] = parse(Float64, csv_data[row, col])
            col += 1
        end
        for g in 1:G, t in 1:T
            reserve[g,t] = parse(Float64, csv_data[row, col])
            col += 1
        end
        for b in 1:B, t in 1:T
            injection[b,t] = parse(Float64, csv_data[row, col])
            col += 1
        end
        violation_txt = strip(csv_data[row, col])
        if length(violation_txt) > 0
            for v in split(violation_txt, " ")
                tokens = parse.(Int, split(v, ":"))
                for t in 1:T
                    push!(violations, Violation(t, tokens[1], tokens[2], 0.0, 0.0))
                end
            end
        end
        push!(results, UnitCommitmentSolution(cost,
                                              is_on,
                                              production,
                                              reserve,
                                              injection,
                                              violations))
    end
    return results
end


export write_csv, read_uc_solution_csv