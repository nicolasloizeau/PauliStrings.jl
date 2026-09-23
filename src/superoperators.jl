# Matrix representations of superoperators (linear maps on operators) in a Pauli string basis.

using SparseArrays

"""
    liouvillian(H::Operator, basis::AbstractVector{<:PauliString}; show_progress=false)
    liouvillian(H::Operator; show_progress=false)

Compute the matrix of the Liouvillian superoperator ``\\mathcal{L}(O) = i[H, O]``
in a basis of Pauli strings:
``L_{ij} = i\\,\\mathrm{tr}(P_i [H, P_j]) / 2^N``.

`H` must be Hermitian, in which case `L` is real and antisymmetric.
If no basis is given, the complete basis of all ``4^N`` strings is used.
Returns a sparse matrix. The coefficient vector `c` of an operator ``O = \\sum_j c_j P_j``
evolves as `dc/dt = L * c`.
"""
liouvillian(H::Operator; kwargs...) = liouvillian(H, complete_basis(qubitlength(H)); kwargs...)
function liouvillian(H::Operator, basis::AbstractVector{<:PauliString}; show_progress=false)
    n = length(basis)
    idx = Dict(p => i for (i, p) in enumerate(basis))
    I, J, V = Int[], Int[], Float64[]
    progress = show_progress ? ProgressBar : identity
    for (j, Pj) in progress(enumerate(basis))
        C = commutator(H, Operator(Pj))
        for (Pi, c) in pairs(C)
            i = get(idx, Pi, 0)
            i == 0 && continue
            v = real(im * c / (1im)^ycount(Pi))
            iszero(v) && continue
            push!(I, i); push!(J, j); push!(V, v)
        end
    end
    return sparse(I, J, V, n, n)
end


"""
    depolarizing_noise(N::Int, γ::Real, basis::AbstractVector{<:PauliString})
    depolarizing_noise(N::Int, γ::Real)

Compute the matrix of the depolarizing noise generator in a basis of Pauli strings.
It is diagonal with entries ``D_{ii} = -γ\\, w_i``, where ``w_i`` is the number of
non-identity Pauli operators in ``P_i``.
Together with [`liouvillian`](@ref), the coefficients of an operator evolve as
`dc/dt = (L + D) * c`, and `exp(D * t)` is the channel implemented by [`add_noise`](@ref)
with amplitude `γ * t`.
If no basis is given, the complete basis of all ``4^N`` strings is used.
Returns a sparse matrix.
"""
depolarizing_noise(N::Int, γ::Real) = depolarizing_noise(N, γ, complete_basis(N))
function depolarizing_noise(N::Int, γ::Real, basis::AbstractVector{<:PauliString})
    all(p -> qubitlength(p) == N, basis) || throw(ArgumentError("basis strings must have $N qubits"))
    return spdiagm([-γ * pauli_weight(p) for p in basis])
end


"""
    depolarizing_noise_channel(N::Int, g::Real, basis::AbstractVector{<:PauliString})
    depolarizing_noise_channel(N::Int, g::Real)

Compute the matrix of the depolarizing noise channel in a basis of Pauli strings.
It is diagonal with entries ``e^{-g\\, w_i}``, where ``w_i`` is the number of
non-identity Pauli operators in ``P_i``. Applying it to the coefficient vector of an
operator is equivalent to [`add_noise`](@ref) with amplitude `g`.
It is `exp(depolarizing_noise(N, g, basis))`, see [`depolarizing_noise`](@ref).
If no basis is given, the complete basis of all ``4^N`` strings is used.
Returns a sparse matrix.
"""
depolarizing_noise_channel(N::Int, g::Real) = depolarizing_noise_channel(N, g, complete_basis(N))
function depolarizing_noise_channel(N::Int, g::Real, basis::AbstractVector{<:PauliString})
    all(p -> qubitlength(p) == N, basis) || throw(ArgumentError("basis strings must have $N qubits"))
    return spdiagm([exp(-g * pauli_weight(p)) for p in basis])
end
