
"""random 2-local operator"""
function rand_local2(N::Int)
    o = Operator(N)
    for i in 1:N
        for j in 1:N
            if i != j
                for k in ['X', 'Y', 'Z']
                    for l in ['X', 'Y', 'Z']
                        o += (randn(rng, Float64), k, i, l, j)
                    end
                end
            end
        end
    end
    return compress(o)
end

"""random 1-local operator"""
function rand_local1(N::Int)
    o = Operator(N)
    for i in 1:N
        for k in ['X', 'Y', 'Z']
            o += (randn(rng, Float64), k, i)
        end
    end
    return compress(o)
end
