
using ProgressBars
using Random

rng = MersenneTwister(0)

function construction_example1()
    A = Operator(4)
    A += "X", 1
    A += "Y", 2
    A += 2, "Z", 3
    return A
end

function construction_example2()
    A = Operator(8)
    A += 2, "X", 1, "X", 2
    A += 1.5, "Y", 2, "Z", 3
    A += "Z", 3, "Z", 8
    return A
end

function construction_example3()
    A = Operator(6)
    A += "XYZ111"
    A += 5, "11XX11"
    return A
end

"""XX hamiltonian https://arxiv.org/abs/1812.08657"""
function XX(N)
    H = Operator(N)
    for j in 1:(N-1)
        H += ('X', j, 'X', j + 1)
        H += ('Z', j, 'Z', j + 1)
    end
    H += ('X', N, 'X', 1)
    H += ('Z', N, 'Z', 1)
    return H
end


"""X operator on all sites"""
function X(N)
    H = Operator(N)
    for j in 1:N
        H += ('X', j)
    end
    return H
end

"""Y operator on all sites"""
function Y(N)
    H = Operator(N)
    for j in 1:N
        H += ('Y', j)
    end
    return H
end

"""Z operator on all sites"""
function Z(N)
    H = Operator(N)
    for j in 1:N
        H += ('Z', j)
    end
    return H
end

"""Z on site i"""
function Zi(N, i)
    H = Operator(N)
    H += ('Z', i)
    return H
end


# build a chaotic spin chain with periodic bc
function chaotic_chain(N::Int)
    H = ps.Operator(N)
    # XX interractions
    for j in 1:(N-1)
        H += "X", j, "X", j + 1
    end
    H += "X", 1, "X", N # close the chain
    # fields
    for j in 1:N
        H += -1.05, "Z", j
        H += 0.5, "X", j
    end
    return H
end



function ising1D(N, g)
    H = Operator(N)
    for j in 1:(N-1)
        H += "Z", j, "Z", j + 1
    end
    H += "Z", 1, "Z", N #periodic boundary condition
    for j in 1:N
        H += g, "X", j
    end
    return -H
end

# 2d ising model with periodic boundary conditions
function ising2D(L1, L2, g)
    H = Operator(L1 * L2)
    for x in 1:L1
        for y in 1:L2
            H += string_2d(("Z", x, y, "Z", x + 1, y), L1, L2, pbc=true) # horizontal
            H += string_2d(("Z", x, y, "Z", x, y + 1), L1, L2, pbc=true) # vertical
            H += g * string_2d(("X", x, y), L1, L2, pbc=true) # transverse field
        end
    end
    return H
end


"""
random 2-local operator with N sites and M terms for use in tests
"""
function rand_local2_M(N::Int, M::Int)
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
random 1-local operator with N sites and M terms for use in tests
"""
function rand_local1_M(N::Int, M::Int)
    o = Operator(N)
    for i in 1:M
        k = rand(['X', 'Y', 'Z'])
        i = rand(1:N)
        o += (randn(rng, Float64), k, i)
    end
    return compress(o)
end


function randstring(N)
    return join(rand(["1", "X", "Y", "Z"], N))
end

function random_string(N::Int)
    symbols = ['1', 'X', 'Y', 'Z']
    s = join([symbols[rand(1:4)] for _ in 1:N])
    return PauliString(s)
end

"""
XXZ model with periodic boundary conditions
"""
function XXZ(N, J, Δ; ts::Bool=false)
    H = Operator(N)
    i = 1
    H += J, "X", i, "X", mod1(i + 1, N)
    H += J, "Y", i, "Y", mod1(i + 1, N)
    H += J * Δ, "Z", i, "Z", mod1(i + 1, N)

    if ts
        return OperatorTS1D(H, full=false)
    else
        return resum(OperatorTS1D(H, full=false))
    end
end
