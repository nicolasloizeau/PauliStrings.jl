
using SparseArrays
"""
    frustration_graph(o::Operator)

Construct the frustration graph, or anticommutation graph of an operator o. Returns a sparse adjacency matrix.
Each vertex represents a pauli string of the operator, and two vertices are connected by an edge if the corresponding pauli strings anticommute.
The vertices are labeled by the index of the pauli string in the operator. Use [`get_pauli`](@ref) to access individual strings.
"""
function frustration_graph(o::Operator)
    o = Operator(o)
    n = length(o)
    G = spzeros(Int, n, n)
    for i in 1:n
        for j in 1:n
            v = o.v[i] ⊻ o.v[j]
            w = o.w[i] ⊻ o.w[j]
            k = (-1)^count_ones(o.v[i] & o.w[j]) - (-1)^count_ones(o.w[i] & o.v[j])
            if k != 0
                G[i, j] = 1
            end
        end
    end
    return G
end
