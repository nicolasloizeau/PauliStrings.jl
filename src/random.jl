
"""
    rand_local2(N::Int)

Random 2-local operator
"""
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

"""
    rand_local2(N::Int)

Random 1-local operator
"""
function rand_local1(N::Int)
    o = Operator(N)
    for i in 1:N
        for k in ['X', 'Y', 'Z']
            o += (randn(rng, Float64), k, i)
        end
    end
    return compress(o)
end


"""
    rand_local1_TS1D(N::Int)

Random 1-local OperatorTS1D
"""
function rand_local1_TS1D(N::Int)
    o = Operator(N)
    for k in ['X', 'Y', 'Z']
        o += (randn(rng, Float64), k, 1)
    end
    return OperatorTS1D(o; full=false)
end

"""
    rand_local2_TS1D(N::Int)

Random 2-local OperatorTS1D
"""
function rand_local2_TS1D(N::Int)
    o = Operator(N)
    for i in 2:N
        for k in ['X', 'Y', 'Z']
            for l in ['X', 'Y', 'Z']
                o += (randn(rng, Float64), k, 1, l, i)
            end
        end
    end
    o = compress(o)
    return OperatorTS1D(o; full=false)
end
