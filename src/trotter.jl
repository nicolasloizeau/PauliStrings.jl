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

function _lie_gates(H::Operator, dt::Real, hbar::Real, heisenberg::Bool)
    gates = TrotterGate{paulistringtype(H),Float64}[]
    for (c, p) in zip(get_coeffs(H), H.strings)
        push!(gates, TrotterGate(p, _trotter_theta(c, dt, hbar, heisenberg)))
    end
    return gates
end


function _strang_gates(H::Operator, dt::Real, hbar::Real, heisenberg::Bool)
    L = length(H)
    gates = TrotterGate{paulistringtype(H),Float64}[]
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

For `H::Operator{<:PauliStringTS}`, see the specialized [`trotterize`](@ref) that builds gates from representative terms only.
"""
function trotterize(H::Operator, dt::Real; order::Integer=2, heisenberg::Bool=true, hbar::Real=1)
    order ∈ (1, 2) || throw(ArgumentError("order must be 1 or 2, got $order"))
    n = qubitlength(H)
    norm(H - H') > 1e-10 && throw(ArgumentError("Hamiltonian must be Hermitian for Trotter splitting"))
    if length(H) == 0
        return TrotterGate{paulistringtype(n),Float64}[]
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
function trotter_step!(O::Operator, gates::AbstractVector{<:TrotterGate}; truncation::Function=identity, truncate_every::Int=1)
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
        for (P, c) in zip(O.strings, O.coeffs)
            C, k = commutator(G, P)
            setwith!(+, d, P, c)# [G, P] as an Operator                     # commuting case
            if k != 0
                # these are much smaller
                setwith!(+, d, P, c*ctheta-c)
                setwith!(+, d, C, c*(1im * stheta / 2) * phase*k )
            end
        end
        ks = Vector{keytype(d)}(undef, length(d))
        vs = Vector{valtype(d)}(undef, length(d))
        for (j, (k, v)) in enumerate(pairs(d))
            @inbounds ks[j] = k
            @inbounds vs[j] = v
        end
        O2 = Operator{keytype(d),valtype(d)}(ks, vs)
        (i%truncate_every == 0) && (O2 = truncation(O2))
        empty!(O.strings)
        empty!(O.coeffs)
        append!(O.strings, O2.strings)
        append!(O.coeffs, O2.coeffs)
    end
    return O
end

# ──────────────────────────────────────────────────
# Translation-symmetric Trotter
# ──────────────────────────────────────────────────

"""
    trotterize(H::Operator{<:PauliStringTS}, dt::Real; order=2, heisenberg=true, hbar=1)

Build a Trotter gate list from the **representative** (unit-cell) terms of the
translation-symmetric Hamiltonian `H`. This produces `M` gates (order 1) or
`2M-1` gates (order 2) instead of `N×M`, where `M = length(H)` is the number
of unique representative terms and `N` is the number of translation shifts.

The returned gates have plain `PauliString` generators (the representatives).
Pass them to the [`trotter_step!`](@ref) method for `OperatorTS`, which
internally applies every translation of each gate.
"""
function trotterize(H::Operator{<:PauliStringTS}, dt::Real; order::Integer=2, heisenberg::Bool=true, hbar::Real=1)
    order ∈ (1, 2) || throw(ArgumentError("order must be 1 or 2, got $order"))
    Hr = representative(H)
    n = qubitlength(H)
    norm(Hr - Hr') > 1e-10 && throw(ArgumentError("Hamiltonian must be Hermitian for Trotter splitting"))
    if length(Hr) == 0
        return TrotterGate{paulistringtype(n),Float64}[]
    end
    if order == 1 || length(Hr) == 1
        return _lie_gates(Hr, dt, hbar, heisenberg)
    else
        return _strang_gates(Hr, dt, hbar, heisenberg)
    end
end


"""
    trotter_step!(O::Operator{<:PauliStringTS}, gates::AbstractVector{<:TrotterGate}; truncation, truncate_every)

Apply one Trotter step to a translation-symmetric operator **in place**.
`gates` should come from [`trotterize`](@ref) on an `OperatorTS` — they contain
only the representative generators. This method internally applies every
translation of each gate and re-symmetrizes the result.
"""
function trotter_step!(O::Operator{<:PauliStringTS}, gates::AbstractVector{<:TrotterGate};
                       truncation::Function=identity, truncate_every::Int=1)
    isempty(gates) && return O
    Ls = qubitsize(O)
    Ps = periodicflags(O)

    # Unwrap to representative (plain Operator with PauliString entries)
    Or = representative(O)
    d = emptydict(Or)

    gate_count = 0
    # U = V1*V2*...*VL => conjugation applies VL,...,V1 (reverse)
    for g in Iterators.reverse(gates)
        G_rep = g.generator
        stheta, ctheta = sincos(g.theta)

        # Apply every translation of this representative gate sequentially
        for s in all_shifts(Ls, Ps)
            gate_count += 1
            G = shift(G_rep, Ls, Ps, s)
            phase = (1.0im)^ycount(G)

            empty!(d)
            for (P, c) in zip(Or.strings, Or.coeffs)
                C, k = commutator(G, P)
                setwith!(+, d, P, c)
                if k != 0
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
            Or = Operator{keytype(d),valtype(d)}(ks, vs)
            (gate_count % truncate_every == 0) && (Or = truncation(Or))
        end
    end

    # Re-symmetrize back into OperatorTS and update O in place
    Ots = OperatorTS{Ls,Ps}(Or)
    empty!(O.strings)
    empty!(O.coeffs)
    append!(O.strings, Ots.strings)
    append!(O.coeffs, Ots.coeffs)
    return O
end