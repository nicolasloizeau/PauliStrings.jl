"""
    TrotterGate{P,T}

One factor in a product-form Trotter step: on each Pauli string, the Heisenberg update
`exp(im * theta * G/2) * P * exp(-im * theta * G/2)` with generator `G`, implemented via [`pauli_rotation`](@ref).
The `theta` field follows that convention.
"""
struct TrotterGate{P<:AbstractPauliString,T<:Real}
    generator::P
    theta::T
end

function _trotter_theta(coeff::Number, dt::Real, hbar::Real, heisenberg::Bool)
    alpha = real(coeff) * dt / hbar
    theta = 2 * alpha
    heisenberg || (theta = -theta)
    return theta
end

function _trotter_pairs(Hc::Operator)
    pairs = Tuple{eltype(Hc.strings),Float64}[]
    for (p, c) in zip(Hc.strings, Hc.coeffs)
        p == one(p) && continue
        lambda = c / (1im)^ycount(p)
        abs(imag(lambda)) > 1e-12 && throw(ArgumentError("Hamiltonian must be Hermitian (real Pauli coefficients); got im coeff $c"))
        push!(pairs, (p, Float64(real(lambda))))
    end
    return pairs
end

function _lie_gates(pairs::Vector{<:Tuple}, dt::Real, hbar::Real, heisenberg::Bool)
    isempty(pairs) && return TrotterGate{paulistringtype(0),Float64}[]
    P0 = first(pairs)[1]
    TrotterGate{typeof(P0),Float64}[TrotterGate(p, _trotter_theta(c, dt, hbar, heisenberg)) for (p, c) in pairs]
end

function _strang_gates(pairs::Vector{<:Tuple}, dt::Real, hbar::Real, heisenberg::Bool)
    isempty(pairs) && return TrotterGate{paulistringtype(0),Float64}[]
    P0 = first(pairs)[1]
    L = length(pairs)
    if L == 1
        p, c = pairs[1]
        return TrotterGate{typeof(P0),Float64}[TrotterGate(p, _trotter_theta(c, dt, hbar, heisenberg))]
    end
    gates = TrotterGate{typeof(P0),Float64}[]
    for j in 1:(L - 1)
        p, c = pairs[j]
        push!(gates, TrotterGate(p, _trotter_theta(c, dt, hbar, heisenberg) / 2))
    end
    pL, cL = pairs[L]
    push!(gates, TrotterGate(pL, _trotter_theta(cL, dt, hbar, heisenberg)))
    for j in (L - 1):-1:1
        p, c = pairs[j]
        push!(gates, TrotterGate(p, _trotter_theta(c, dt, hbar, heisenberg) / 2))
    end
    return gates
end

"""
    trotterize(H::Operator, dt::Real; order=2, heisenberg=true, hbar=1)

Build a first-order (`order=1`, Lie) or second-order (`order=2`, Strang) Trotter list that approximates
`exp(im * H * dt / hbar)` (Heisenberg) or the conjugate sequence (Schrödinger / density matrix).
Each gate uses [`pauli_rotation`](@ref) with the returned `theta` field.
"""
function trotterize(H::Operator, dt::Real; order::Integer=2, heisenberg::Bool=true, hbar::Real=1)
    order ∈ (1, 2) || throw(ArgumentError("order must be 1 or 2, got $order"))
    Hc = compress(H)
    n = qubitlength(Hc)
    norm(Hc - dagger(Hc)) > 1e-10 && throw(ArgumentError("Hamiltonian must be Hermitian for Trotter splitting"))
    pairs = _trotter_pairs(Hc)
    if isempty(pairs)
        return TrotterGate{paulistringtype(n),Float64}[]
    end
    if order == 1
        return _lie_gates(pairs, dt, hbar, heisenberg)
    else
        return _strang_gates(pairs, dt, hbar, heisenberg)
    end
end

"""
    trotter_step!(O::Operator, gates; M=2^20, keep=Operator(0))

Apply one Trotter step in place. Gates must be listed in matrix-multiply order `U = V1 * V2 * ... * Vn`;
conjugation `O -> U * O * U'` applies factors `Vn, ..., V1` successively (reverse of the list).
Each Pauli string uses the same coefficient convention as [`Matrix`](@ref)(`O`) (weights include division by `im` to the number of `Y` factors).
"""
function trotter_step!(O::Operator, gates::AbstractVector{<:TrotterGate}; M::Int=2^20, keep::Operator=Operator(0))
    qubitlength(keep) == 0 && (keep = Operator(qubitlength(O)))
    isempty(gates) && return O
    N = qubitlength(O)
    # U = V1*V2*...*VL => U O U' = V1*...*VL * O * VL'*...*V1' ; apply VL,...,V1 (reverse list order).
    for g in Iterators.reverse(gates)
        acc = Operator(N)
        n = length(O.strings)
        acc_strings = sizehint!(similar(O.strings, 0), 2n)
        acc_coeffs  = sizehint!(similar(O.coeffs, 0), 2n)
        for (P, c) in zip(O.strings, O.coeffs)
            G = g.generator
            C, k = commutator(G, P)              # [G, P] as an Operator
            if k == 0                           # commuting case
                push!(acc_strings, P)
                push!(acc_coeffs, c)
            else
                stheta, ctheta = sincos(g.theta)
                push!(acc_strings, P)
                push!(acc_strings, C)
                push!(acc_coeffs, c*ctheta)
                push!(acc_coeffs, c*(1im * stheta / 2) * k*(1im)^ycount(C)/(1im)^ycount(P))
            end
        end
        acc = Operator(acc_strings, acc_coeffs)
        acc = compress(acc)
        acc = trim(acc, M; keep=keep)
        empty!(O.strings)
        empty!(O.coeffs)
        append!(O.strings, acc.strings)
        append!(O.coeffs, acc.coeffs)
    end
    return O
end

"""
    trotter_evolve(H::Operator, O::Operator, dt::Real, nsteps; order=2, heisenberg=true, hbar=1, gates=nothing, M=2^20, keep=Operator(0))

Apply `trotter_step!` `nsteps` times. If `gates` is passed, it is reused (from a prior [`trotterize`](@ref));
otherwise gates are built from `H` each call would be wasteful — pass `gates=trotterize(H, dt; ...)` or let this function build once:

When `gates === nothing`, `trotterize(H, dt; order, heisenberg, hbar)` is called once and reused for all steps.
"""
function trotter_evolve(
    H::Operator, O::Operator, dt::Real, nsteps::Integer;
    gates=nothing,
    order::Integer=2,
    heisenberg::Bool=true,
    hbar::Real=1,
    M::Int=2^20,
    keep::Operator=Operator(0),
)
    nsteps < 0 && throw(ArgumentError("nsteps must be non-negative"))
    qubitlength(H) == qubitlength(O) || throw(DimensionMismatch("H and O must act on the same number of qubits"))
    g = gates === nothing ? trotterize(H, dt; order, heisenberg, hbar) : gates
    for _ in ProgressBar(1:nsteps)
        trotter_step!(O, g; M, keep)
    end
    return O
end
