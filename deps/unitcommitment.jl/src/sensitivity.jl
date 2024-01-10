# UnitCommitment.jl: Optimization Package for Security-Constrained Unit Commitment
# Copyright (C) 2020-2024, UChicago Argonne, LLC. All rights reserved.
# Released under the modified BSD license. See COPYING.md for more details.

using SparseArrays, Base.Threads

function reduced_incidence_matrix(lines::Array{TransmissionLine, 1}) :: SparseMatrixCSC{Float64}
    n_buses = maximum([max(l.source.index, l.target.index) for l in lines])
    n_lines = size(lines)[1]
    matrix = spzeros(Float64, n_lines, n_buses)
    for line in lines
        matrix[line.index, line.source.index] = 1
        matrix[line.index, line.target.index] = -1
    end
    return matrix[1:n_lines, 2:n_buses]
end

function susceptance_matrix(lines::Array{TransmissionLine, 1}) :: SparseMatrixCSC{Float64}
    n_lines = size(lines)[1]
    matrix = spzeros(Float64, n_lines, n_lines)
    for line in lines
        matrix[line.index, line.index] = line.susceptance
    end
    return matrix
end

function injection_shift_factors(lines::Array{TransmissionLine, 1}) :: Array{Float64, 2}
    susceptance = susceptance_matrix(lines)
    incidence = reduced_incidence_matrix(lines)
    laplacian = transpose(incidence) * susceptance * incidence
    isf = susceptance * incidence * inv(Array(laplacian))
    return [zeros(length(lines)) isf]
end

function change_slack!(isf::Array{Float64, 2},
                       prev_slack::Int,
                       new_slack::Int)
    m, n = size(isf)
    slack_col = isf[:, new_slack]
    for i in 1:n
        isf[:, i] -= slack_col
    end
end

function line_outage_factors(lines::Array{TransmissionLine, 1},
                             isf::Array{Float64,2}) :: Array{Float64,2}

    n_lines, n_buses = size(isf)
    incidence = Array(reduced_incidence_matrix(lines))
    lodf::Array{Float64,2} = isf[:, 2:n_buses] * transpose(incidence)
    m, n = size(lodf)
    for i in 1:n
        lodf[:, i] *= 1.0 / (1.0 - lodf[i, i])
        lodf[i, i] = -1
    end
    return lodf
end

function post_contingency_isf(isf::Array{Float64, 2},
                              lodf::Array{Float64, 2},
                              outage::Int) :: Array{Float64, 2}
    L, B = size(isf)
    return [isf[l,b] + lodf[l, outage] * isf[outage, b] for l in 1:L, b in 1:B]
end

export reduced_incidence_matrix, susceptance_matrix, injection_shift_factors,
    change_slack!, line_outage_factors, post_contingency_isf
