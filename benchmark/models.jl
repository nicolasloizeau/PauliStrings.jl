

module models

using PauliStrings
include("lattices.jl")

"""
1D XXZ + X field from https://arxiv.org/pdf/1905.08266 eq 2
"""
function XXZh(N::Int, n, m, h)
    Δ = cos(pi * n / m)
    H = Operator(N)
    H += 'X', 1, 'X', 2
    H += 'Y', 1, 'Y', 2
    H += Δ, 'Z', 1, 'Z', 2
    H += h, 'Z', 1
    return OperatorTS1D(H, full=false)
end
XXZh(N::Int) = XXZh(N, 2, 3, 2)

"""
first order of the dynamical symmetry in the 1D XXZ + field
"""
function XXZo3(N::Int)
    H = Operator(N)
    H += "S+", 1, "S+", 2, "S+", 3
    H += "S-", 1, "S-", 2, "S-", 3
    return OperatorTS1D(H, full=false)
end

"""
hamiltonian from fig 4b of https://arxiv.org/pdf/1812.08657
"""
function chaotic_ising(N::Int, h::Real)
    H = Operator(N)
    H += 'X', 1, 'X', 2
    H += -1.05, 'Z', 1
    H += h, 'X', 1
    return OperatorTS1D(H, full=false)
end
chaotic_ising(N::Int) = chaotic_ising(N, 0.5)


"""
operator from fig 4b of https://arxiv.org/pdf/1812.08657
"""
function chaotic_ising_op(N::Int)
    H = Operator(N)
    H += 1.05, 'X', 1, 'X', 2
    H += 'Z', 1
    return OperatorTS1D(H, full=false)
end


"""
1D XXZ next nearest neigbor from https://arxiv.org/pdf/1808.08977 eq 3
"""
function XXZnnn(N::Int)
    Δ = 2
    γ = 1 / 2
    H = Operator(N)
    H += "X", 1, "X", 2
    H += "Y", 1, "Y", 2
    H += Δ, "Z", 1, "Z", 2
    H += γ, "X", 1, "X", 3
    H += γ, "Y", 1, "Y", 3
    H += γ * Δ, "Z", 1, "Z", 3
    return OperatorTS1D(H, full=false)
end

"""
Z on site 1
"""
function Z1(N::Int)
    H = Operator(N)
    H += "Z", 1
    return H
end


function XXZ(N::Int, delta::Real)
    H = Operator(N)
    H += "X", 1, "X", 2
    H += "Y", 1, "Y", 2
    H += delta, "Z", 1, "Z", 2
    return OperatorTS1D(H, full=false)
end
XXZ(N) = XXZ(N, 0.5)

"""
Many-body localized hamiltonian: XXZ + random Z field
"""
function mbl(N::Int, W::Real)
    H = Operator(XXZ(N))
    for i in 1:N
        hi = W * (rand() - 0.5) * 2
        H += hi, "Z", i
    end
    return H
end
mbl(N::Int) = mbl(N, 10.0)


function XXZ2D(n::Int, m::Int, delta::Real; bc="closed")
    H = Operator(n * m)
    for (i, j) in square_lattice(n, m; bc=bc)
        H += "X", i, "X", j
        H += "Y", i, "Y", j
        H += delta, "Z", i, "Z", j
    end
    return H
end
XXZ2D(n::Int, m::Int) = XXZ2D(n, m, 0.5; bc="closed")
function XXZ2D(N::Int)
    n = floor(Int, sqrt(N))
    m = div(N, n)
    @assert n * m == N "N must be a perfect square"
    return XXZ2D(n, m, 0.5; bc="closed")
end

end
