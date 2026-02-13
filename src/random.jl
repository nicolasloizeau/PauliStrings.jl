
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
    rand_local2(N::Int, M::Int)

Random 2-local operator with N sites and M terms
"""
function rand_local2(N::Int, M::Int)
    o = Operator(N)
    for i in 1:M
        k = rand(['X', 'Y', 'Z'])
        l = rand(['X', 'Y', 'Z'])
        i = rand(1:N)
        j = i
        while j == i
            j = rand(1:N)
        end
        o += (randn(rng, Float64), k, i, l, j)
    end
    return compress(o)
end




"""
    rand_local1(N::Int)

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

Random 1-local OperatorTS in one dimension
"""
function rand_local1_TS1D(N::Int)
    o = Operator(N)
    for k in ['X', 'Y', 'Z']
        o += (randn(rng, Float64), k, 1)
    end
    return OperatorTS{(N,)}(o)
end

"""
    rand_local2_TS1D(N::Int)

Random 2-local OperatorTS in one dimension
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
    return OperatorTS{(N,)}(o)
end



function randstring(N)
    return join(rand(["1", "X", "Y", "Z"], N))
end


"""
    random_string(N::Int)

Random PauliString of length N
"""
function rand_pauli_string(N::Int)
    symbols = ['1', 'X', 'Y', 'Z']
    s = join([symbols[rand(1:4)] for _ in 1:N])
    return PauliString(s)
end
