
using Graphs

function frustration_graph(o)
    G = SimpleGraph(length(o))
    for i in 1:length(o)
        for j in 1:length(o)
            v = o.v[i] ⊻ o.v[j]
            w = o.w[i] ⊻ o.w[j]
            k = (-1)^count_ones(o.v[i] & o.w[j]) - (-1)^count_ones(o.w[i] & o.v[j])
            if k != 0
                add_edge!(G, i, j)
            end
        end
    end
    return G
end
