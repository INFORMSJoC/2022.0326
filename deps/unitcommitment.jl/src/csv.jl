# UnitCommitment.jl: Optimization Package for Security-Constrained Unit Commitment
# Copyright (C) 2020-2024, UChicago Argonne, LLC. All rights reserved.
# Released under the modified BSD license. See COPYING.md for more details.

struct CSV
    column_names::Array{String,1}
    data::Array{String,2}
end

function read_csv(path::String;
                  delimiter::String = ",") :: CSV

    ncols::Int = 0
    nlines::Int = countlines(path)
    data::Array{String,2} = Array{String}(undef, 0, 0)
    column_names::Array{String,1} = Array{String}(undef, 0)

    current_line = 0
    for line in eachline(path)
        current_line += 1
        tokens = split(line, delimiter)

        if current_line == 1
            ncols = length(tokens)
            column_names = tokens
            data = Array{String}(undef, nlines - 1, ncols)
        else
            length(tokens) == ncols || @error "Line $current_line has wrong number of columns; $ncols expected, $(length(tokens)) found"
            data[current_line - 1, :] = tokens
        end
    end

    return CSV(column_names, data)
end

export read_csv
