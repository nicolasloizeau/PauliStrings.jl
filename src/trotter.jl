"""
    TrotterGate{P,T}

One factor in a product-form Trotter step: on each Pauli string, the Heisenberg update
`exp(im * theta * G/2) * P * exp(-im * theta * G/2)` with generator `G`, implemented via [`pauli_rotation`](@ref).
The `theta` field follows that convention.
"""
struct TrotterGate{P <: AbstractPauliString, T <: Real}
    generator::P
    theta::T
end

function _trotter_theta(coeff::Number, dt::Real, hbar::Real, heisenberg::Bool)
    alpha = real(coeff) * dt / hbar
    theta = 2 * alpha
    heisenberg || (theta = -theta)
    return theta
end

function _lie_gates(H::Operator, dt::Real, hbar::Real, heisenberg::Bool)
    gates = TrotterGate{paulistringtype(H), Float64}[]
    for (c, p) in zip(get_coeffs(H), keys(H))
        push!(gates, TrotterGate(p, _trotter_theta(c, dt, hbar, heisenberg)))
    end
    return gates
end


function _strang_gates(H::Operator, dt::Real, hbar::Real, heisenberg::Bool)
    L = length(H)
    gates = TrotterGate{paulistringtype(H), Float64}[]
    for j in 1:(L - 1)
        c, p = H[j]
        push!(gates, TrotterGate(p, _trotter_theta(c, dt, hbar, heisenberg) / 2))
    end
    c, p = H[L]
    push!(gates, TrotterGate(p, _trotter_theta(c, dt, hbar, heisenberg)))
    for j in (L - 1):-1:1
        c, p = H[j]
        push!(gates, TrotterGate(p, _trotter_theta(c, dt, hbar, heisenberg) / 2))
    end
    return gates
end

"""
    trotterize(H::Operator, dt::Real; order=2, heisenberg=true, hbar=1)

Build a first-order (`order=1`, Lie) or second-order (`order=2`, Strang) Trotter list that approximates
`exp(im * H * dt / hbar)` (Heisenberg) or the conjugate sequence (Schrödinger / density matrix).
Each gate uses [`pauli_rotation`](@ref) with the returned `theta` field.

For `H::Operator{<:PauliStringTS}`, see the specialized [`trotterize`](@ref) that calls [`resum`](@ref) first.
"""
function trotterize(H::Operator, dt::Real; order::Integer = 2, heisenberg::Bool = true, hbar::Real = 1)
    order ∈ (1, 2) || throw(ArgumentError("order must be 1 or 2, got $order"))
    n = qubitlength(H)
    norm(H - H') > 1.0e-10 && throw(ArgumentError("Hamiltonian must be Hermitian for Trotter splitting"))
    if length(H) == 0
        return TrotterGate{paulistringtype(n), Float64}[]
    end
    if order == 1 || length(H) == 1
        return _lie_gates(H, dt, hbar, heisenberg)
    else
        return _strang_gates(H, dt, hbar, heisenberg)
    end
end

"""
trotter_step!(O::Operator, gates::AbstractVector{<:TrotterGate}; truncation::Function=identity, truncate_every::Int=1)

Apply one Trotter step in place. Gates must be listed in matrix-multiply order `U = V1 * V2 * ... * Vn`;
conjugation `O -> U * O * U'` applies factors `Vn, ..., V1` successively (reverse of the list).
Each Pauli string uses the same coefficient convention as [`Matrix`](@ref)(`O`) (weights include division by `im` to the number of `Y` factors).
"""
function trotter_step!(O::Operator, gates::AbstractVector{<:TrotterGate}; truncation::Function = identity, truncate_every::Int = 1)
    isempty(gates) && return O
    N = qubitlength(O)
    d = emptydict(O)
    # U = V1*V2*...*VL => U O U' = V1*...*VL * O * VL'*...*V1' ; apply VL,...,V1 (reverse list order).
    for (i, g) in enumerate(Iterators.reverse(gates))
        n = length(O.strings)
        empty!(d)
        G = g.generator
        stheta, ctheta = sincos(g.theta)
        phase = (1.0im)^ycount(G)
        for (P, c) in pairs(O)
            C, k = commutator(G, P)
            setwith!(+, d, P, c) # [G, P] as an Operator                     # commuting case
            if k != 0
                # these are much smaller
                setwith!(+, d, P, c * ctheta - c)
                setwith!(+, d, C, c * (1im * stheta / 2) * phase * k)
            end
        end
        ks = Vector{keytype(d)}(undef, length(d))
        vs = Vector{valtype(d)}(undef, length(d))
        for (j, (k, v)) in enumerate(pairs(d))
            @inbounds ks[j] = k
            @inbounds vs[j] = v
        end
        O2 = Operator{keytype(d), valtype(d)}(ks, vs)
        (i % truncate_every == 0) && (O2 = truncation(O2))
        empty!(O.strings)
        empty!(O.coeffs)
        append!(O.strings, O2.strings)
        append!(O.coeffs, O2.coeffs)
    end
    return O
end
