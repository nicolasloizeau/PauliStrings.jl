

using Base.Iterators


"""
    trace_product(o1::Operator, o2::Operator; scale=0)
    trace_product(o1::OperatorTS, o2::OperatorTS; scale=0)

Efficiently compute `trace(o1*o2)`. This is much faster than doing `trace(o1*o2)`.
If `scale` is not 0, then the result is normalized such that trace(identity)=scale.
"""
function trace_product(o1::Operator, o2::Operator; scale=0)
    # operation is symmetric but more efficient if o1 is the largest collection
    (length(o1) < length(o2)) && return trace_product(o2, o1; scale)

    checklength(o1, o2)
    N = qubitlength(o1)
    tr = zero(scalartype(o1))

    # ensure `@inbounds` is safe
    length(o1.strings) == length(o1.coeffs) || throw(DimensionMismatch("strings and coefficients must have the same length"))
    length(o2.strings) == length(o2.coeffs) || throw(DimensionMismatch("strings and coefficients must have the same length"))

    # trace of product contributes only if product is 1, which only happens when strings are equal
    # this amounts to `indexin`, which we hijack/reimplement here for efficiency
    d = emptydict(o2)
    @inbounds for (p, c) in pairs(o2)
        insert!(d, p, c)
    end

    @inbounds for (p1, c1) in pairs(o1)
        c2 = get(d, p1, nothing)
        # TODO: verify if c2 = zero(c1) without branch is faster implementation
        isnothing(c2) && continue
        p, k = prod(p1, p1)
        tr += c1 * c2 * k
    end

    (scale == 0) && (scale = 2.0^N)
    return tr * scale
end

function trace_product(o1::Operator{<:PauliStringTS}, o2::Operator{<:PauliStringTS}; scale=0)
    checklength(o1, o2)
    Ls = qubitsize(o1)
    Ps = periodicflags(o1)
    tr = zero(scalartype(o1))

    # see above
    d = emptydict(o2)
    for (p2, c2) in pairs(o2)
        insert!(d, p2, c2)
    end

    for (p1, c1) in pairs(o1)
        c2 = get(d, p1, nothing)
        isnothing(c2) && continue
        rep1 = representative(p1)
        p, k = prod(rep1, rep1)
        f = c1 * c2 * k
        for s in all_shifts(Ls, Ps)
            shifted = shift(rep1, Ls, Ps, s)
            if shifted == rep1
                tr += f
            end
        end
    end
    (iszero(scale)) && (scale = 2.0^Base.prod(Ls))
    # Calculate the number of translations: product of lengths for periodic dimensions only
    num_translations = Base.prod(L for (L, p) in zip(Ls, Ps) if p)
    return tr * scale * num_translations
end

Base.@deprecate oppow(o::AbstractOperator, k::Int) o^k

"""
    Base.:^(o::Operator, k::Int)

kth power of o.
"""
Base.:^(o::AbstractOperator, k::Int) = Base.power_by_squaring(o, k)

"""
    trace_product(A::Operator, k::Int, B::Operator, l::Int; scale=0)

Efficiently compute `trace(A^k*B^l)`. This is much faster than doing `trace(A^k*B^l)`.

If `scale` is not 0, then the result is normalized such that trace(identity)=scale.
"""
function trace_product(A::AbstractOperator, k::Int, B::AbstractOperator, l::Int; scale=0)
    @assert typeof(A) == typeof(B)
    m = div(k + l, 2)
    n = k + l - m
    if k < m
        C = A^k * B^(m - k)
        D = B^n
    elseif k > m
        C = A^m
        D = A^(k - m) * B^l
    else
        C = A^k
        D = B^l
    end
    return trace_product(C, D; scale=scale)
end


"""
    trace_product(A::AbstractOperator; scale=0)

Compute `trace(A*A)`. This is much faster than doing `trace(A*A)`.

If `scale` is not 0, then the result is normalized such that trace(identity)=scale.
"""
function trace_product(A::Operator; scale=0)
    c = get_coeffs(A)
    N = qubitlength(A)
    return sum(c.^2) * (iszero(scale) ? 2.0^N : scale)
end



"""
    trace_product(A::Operator{<:PauliStringTS}; scale=0)

Compute `trace(A*A)`. This is much faster than doing `trace(A*A)`.

If `scale` is not 0, then the result is normalized such that trace(identity)=scale.
"""
trace_product(A::Operator{<:PauliStringTS}; scale=0) = trace_product(A, A; scale=scale)



"""
    trace_product(A::AbstractOperator, k::Int; scale=0)

Efficiently compute `trace(A^k)`. This is much faster than doing `trace(A^k)`.

If `scale` is not 0, then the result is normalized such that trace(identity)=scale.
"""
function trace_product(A::AbstractOperator, k::Int; scale=0)
    m = div(k, 2)
    n = k - m
    C = A^m
    (k%2 == 0) && (return trace_product(C; scale=scale))
    D = A^n
    return trace_product(C, D; scale=scale)
end

"""
    trace_product_z(o1::AbstractOperator, o2::AbstractOperator; scale=0)

Efficiently compute `<0|o1*o2|0>`.
If `scale` is not 0, then the result is normalized such that `trace(identity) = scale`.
"""
function trace_product_z(o1::AbstractOperator, o2::AbstractOperator; scale=0)
    scale = iszero(scale) ? 2.0^qubitlength(o1) : scale
    tr = zero(scalartype(o1))

    for i in eachindex(o1.strings)
        p1, c1 = o1.strings[i], o1.coeffs[i]
        for j in eachindex(o2.strings)
            p2, c2 = o2.strings[j], o2.coeffs[j]

            p, k = prod(p1, p2)
            if xcount(p) == ycount(p) == 0
                tr += c1 * c2 * k
            end
        end
    end

    return tr * scale
end

"""
    moments(H::AbstractOperator, kmax::Int; start=1, scale=0)

Compute the first kmax moments of H.
start is the first moment to start from.

If scale is not 0, then the result is normalized such that trace(identity)=scale.
"""
function moments(H::AbstractOperator, kmax::Int; start=1, scale=0)
    return [trace_product(H, k; scale=scale) for k in start:kmax]
end


# ============================================================================
# trace_moment for translation-symmetric operators  (issue #80)
# ----------------------------------------------------------------------------
# tr(H^k) = Σ_{l1..lk} c_{l1}…c_{lk} tr(P_{l1}…P_{lk}). A term contributes only if the
# XOR of the chosen Pauli strings is the identity (v=0,w=0); whether it vanishes depends
# only on the *multiset* of Paulis, not their order. We therefore enumerate the
# identity-yielding multisets and sum the ordering-phase analytically, never building
# H^{k/2}. For a translation-symmetric operator we first `resum` it to its degree-1 dense
# form (cheap — it is the Hamiltonian, not a power), then run the multiset enumeration on
# that small alphabet. This is the OperatorTS counterpart of the plain-Operator moment
# method and is much faster than `trace_product(o,k)` (which builds H^{k/2}) for high/odd k.
#
# Phase factorization: reordering two Paulis flips the product phase by their (anti)commutation
# sign, so   Σ_{orderings} σ(ordering) = K_ref(M) · D̃(M),
# where K_ref is the sign of one canonical ordering (O(k)) and D̃(M) depends only on the
# anticommutation graph of the distinct terms: every term commuting with all others peels
# off as a binomial coefficient, leaving a small "anticommuting core" handled by a
# mixed-radix multiplicity DP. The i^ycount Clifford phases are carried inside the stored
# coefficients, so all phases here are real ±1 (same convention as `trace_product`).
# ============================================================================

# true iff PauliStrings a and b anticommute (odd number of non-commuting sites)
@inline _anticommutes(a::P, b::P) where {P<:PauliString} =
    isodd(count_ones(a.v & b.w) + count_ones(a.w & b.v))

mutable struct _TSMomentWS{P,Cc,Tv,C}
    terms::Vector{P}      # M' alphabet (resummed Pauli strings), sorted by lowest active site
    coeffs::Vector{Cc}    # matching (physical) coefficients
    n::Int                # number of alphabet terms
    reach::Vector{Tv}     # reach[i] = OR of (v|w) over terms[i:n]  (suffix support union)
    maxsupp::Int          # max pauli_weight over the alphabet
    tv::Tv                # target XOR (v) the chosen multiset must reach (= anchor rep for folding)
    tw::Tv                # target XOR (w)
    anchorv::Tv           # v-mask of the anchored first factor (for the canonical sign prefix)
    anchorc::Cc           # coefficient prefactor of the anchored first factor
    idx::Vector{Int}      # distinct term indices on the DFS stack (length k)
    mult::Vector{Int}     # their multiplicities (length k)
    acc::C                # accumulator
    keep::Vector{Bool}    # scratch: terms not yet peeled (length k)
    core::Vector{Int}     # scratch: anticommuting-core positions (length k)
    strides::Vector{Int}  # scratch: mixed-radix strides (length k)
    digits::Vector{Int}   # scratch: mixed-radix digits (length k)
    cbuf::Vector{Int128}  # scratch: DP buffer
end

# allocate a fresh workspace (own scratch) sharing the read-only alphabet arrays
function _new_ws(sterms::Vector{P}, scoeffs::Vector{Cc}, n::Int, reach::Vector{Tv},
        maxsupp::Int, cap::Int, ::Type{C}) where {P,Cc,Tv,C}
    return _TSMomentWS{P,Cc,Tv,C}(
        sterms, scoeffs, n, reach, maxsupp,
        zero(Tv), zero(Tv), zero(Tv), one(Cc),
        Vector{Int}(undef, cap), Vector{Int}(undef, cap), zero(C),
        Vector{Bool}(undef, cap), Vector{Int}(undef, cap),
        Vector{Int}(undef, cap), Vector{Int}(undef, cap), Int128[],
    )
end

# build the moment-enumeration alphabet: sort terms by lowest active site (so the reachability
# prune resolves low sites early), the suffix-union `reach`, and the max single-term support.
function _moment_ws_data(strings::Vector{P}, coeffs::Vector{Cc}) where {P,Cc}
    n = length(strings)
    perm = sortperm(strings; by = p -> _minsite(p))
    sterms = strings[perm]
    scoeffs = coeffs[perm]
    Tv = typeof(sterms[1].v)
    reach = Vector{Tv}(undef, n)
    acc_reach = zero(Tv)
    @inbounds for i in n:-1:1
        acc_reach |= (sterms[i].v | sterms[i].w)
        reach[i] = acc_reach
    end
    maxsupp = 0
    @inbounds for p in sterms
        maxsupp = max(maxsupp, pauli_weight(p))
    end
    maxsupp == 0 && (maxsupp = 1)     # alphabet is all-identity; avoid 0 in the prune bound
    return sterms, scoeffs, n, reach, maxsupp
end

# sign (±1) of the ordered product of the canonical ordering of the current multiset,
# with the anchored first factor prepended (its v-mask seeds the running product).
# (canonical = anchor, then idx[1] repeated mult[1] times, then idx[2] …, the DFS discovery order)
@inline function _canonical_sign(ws::_TSMomentWS{P,Cc,Tv,C}, r::Int) where {P,Cc,Tv,C}
    Av = ws.anchorv          # the anchor contributes phase 1 (Av starts at 0), then leaves Av=anchor.v
    sgn = 1
    @inbounds for j in 1:r
        q = ws.terms[ws.idx[j]]
        qv = q.v
        qw = q.w
        for _ in 1:ws.mult[j]
            sgn *= 1 - ((count_ones(Av & qw) & 1) << 1)
            Av ⊻= qv
        end
    end
    return sgn
end

# D̃(M): signed sum over distinct orderings of (-1)^(anticommuting inversions vs canonical).
# Commuting terms peel off as binomials; the anticommuting core is summed by a mixed-radix DP.
function _dtilde!(ws::_TSMomentWS{P,Cc,Tv,C}, r::Int) where {P,Cc,Tv,C}
    @inbounds for j in 1:r
        ws.keep[j] = true
    end
    factor = Int128(1)
    ktot = 0
    @inbounds for j in 1:r
        ktot += ws.mult[j]
    end
    # peel terms that commute with all other kept terms
    changed = true
    while changed
        changed = false
        @inbounds for j in 1:r
            ws.keep[j] || continue
            anti = false
            qj = ws.terms[ws.idx[j]]
            for c in 1:r
                (c == j || !ws.keep[c]) && continue
                if _anticommutes(qj, ws.terms[ws.idx[c]])
                    anti = true
                    break
                end
            end
            if !anti
                factor *= Int128(binomial(ktot, ws.mult[j]))
                ktot -= ws.mult[j]
                ws.keep[j] = false
                changed = true
            end
        end
    end
    # collect anticommuting core
    s = 0
    @inbounds for j in 1:r
        if ws.keep[j]
            s += 1
            ws.core[s] = j
        end
    end
    s == 0 && return factor
    return factor * _core_phase_dp!(ws, s)
end

# mixed-radix DP over the anticommuting core of size s (core positions in ws.core[1:s])
function _core_phase_dp!(ws::_TSMomentWS{P,Cc,Tv,C}, s::Int) where {P,Cc,Tv,C}
    total = 1
    @inbounds for a in 1:s
        ws.strides[a] = total
        total *= (ws.mult[ws.core[a]] + 1)
    end
    length(ws.cbuf) < total && resize!(ws.cbuf, total)
    c = ws.cbuf
    @inbounds c[1] = Int128(1)      # code 0 = empty
    @inbounds for code in 1:total-1
        tmp = code
        for a in 1:s
            mj1 = ws.mult[ws.core[a]] + 1
            ws.digits[a] = tmp % mj1
            tmp ÷= mj1
        end
        val = Int128(0)
        for a in 1:s
            ws.digits[a] > 0 || continue
            qa = ws.terms[ws.idx[ws.core[a]]]
            e = 0
            for b in (a+1):s
                if _anticommutes(qa, ws.terms[ws.idx[ws.core[b]]])
                    e += ws.digits[b]
                end
            end
            sgn = iseven(e) ? Int128(1) : Int128(-1)
            val += sgn * c[code + 1 - ws.strides[a]]
        end
        c[code+1] = val
    end
    return c[total]
end

# contribution of the (anchor ∪ multiset) configuration currently on the stack
@inline function _leaf_contribution(ws::_TSMomentWS{P,Cc,Tv,C}, r::Int) where {P,Cc,Tv,C}
    w = ws.anchorc           # coefficient of the anchored first factor (1 when unfolded)
    @inbounds for j in 1:r
        w *= ws.coeffs[ws.idx[j]]^ws.mult[j]
    end
    phase = _canonical_sign(ws, r) * _dtilde!(ws, r)
    return C(w) * Float64(phase)
end

# accumulate the contribution of the configuration currently on the stack
@inline function _ts_leaf!(ws::_TSMomentWS{P,Cc,Tv,C}, r::Int) where {P,Cc,Tv,C}
    ws.acc += _leaf_contribution(ws, r)
    return nothing
end

# pruned DFS over the alphabet: choose a multiplicity for term i, then advance.
# The chosen multiset M' must reach the target XOR (ws.tv, ws.tw).
function _ts_moment_dfs!(ws::_TSMomentWS{P,Cc,Tv,C}, i::Int, Rv::Tv, Rw::Tv, rem::Int, depth::Int) where {P,Cc,Tv,C}
    if rem == 0
        (Rv == ws.tv && Rw == ws.tw) && _ts_leaf!(ws, depth)
        return nothing
    end
    i > ws.n && return nothing
    # deficit = sites where the running XOR still differs from the target
    deficit = (Rv ⊻ ws.tv) | (Rw ⊻ ws.tw)
    # prune: the deficit sites must be reachable by terms i:n
    @inbounds if (deficit & ~ws.reach[i]) != zero(Tv)
        return nothing
    end
    # prune: need at least ceil(count/maxsupp) more terms to clear the deficit
    cnt = count_ones(deficit)
    cnt > rem * ws.maxsupp && return nothing
    # branch 1: skip term i
    _ts_moment_dfs!(ws, i + 1, Rv, Rw, rem, depth)
    # branch 2: use term i with multiplicity t = 1..rem
    @inbounds ws.idx[depth+1] = i
    vi = ws.terms[i].v
    wi = ws.terms[i].w
    @inbounds for t in 1:rem
        ws.mult[depth+1] = t
        if isodd(t)
            _ts_moment_dfs!(ws, i + 1, Rv ⊻ vi, Rw ⊻ wi, rem - t, depth + 1)
        else
            _ts_moment_dfs!(ws, i + 1, Rv, Rw, rem - t, depth + 1)
        end
    end
    return nothing
end

# enumerate only the M' whose smallest used term index is exactly `i0` (forces term i0 used).
# Partitions the search space into independent subtrees for multithreading.
function _ts_seed_dfs!(ws::_TSMomentWS{P,Cc,Tv,C}, i0::Int, rem::Int) where {P,Cc,Tv,C}
    i0 > ws.n && return nothing
    @inbounds ws.idx[1] = i0
    vi = ws.terms[i0].v
    wi = ws.terms[i0].w
    z = zero(Tv)
    @inbounds for t in 1:rem
        ws.mult[1] = t
        if isodd(t)
            _ts_moment_dfs!(ws, i0 + 1, vi, wi, rem - t, 1)
        else
            _ts_moment_dfs!(ws, i0 + 1, z, z, rem - t, 1)
        end
    end
    return nothing
end

# worker (function barrier) for one independent seed (target/anchor, smallest-used-term i0).
# Builds its own workspace so it is safe to call concurrently from different threads.
# `tv,tw` = target XOR the multiset must reach; `anchorv` = v-mask seeding the canonical sign
# (equal to the target for the translation-folded case, zero otherwise); `anchorc` = prefactor.
function _ts_seed_partial(sterms::Vector{P}, scoeffs::Vector{Cc}, n::Int, reach::Vector{Tv},
        maxsupp::Int, cap::Int, tv::Tv, tw::Tv, anchorv::Tv, anchorc::Cc, i0::Int, rem::Int, ::Type{C}) where {P,Cc,Tv,C}
    ws = _new_ws(sterms, scoeffs, n, reach, maxsupp, cap, C)
    ws.tv = tv
    ws.tw = tw
    ws.anchorv = anchorv
    ws.anchorc = anchorc
    _ts_seed_dfs!(ws, i0, rem)
    return ws.acc
end

# unpruned DFS enumerating ALL size-`rem` multisets of the alphabet; records each multiset's
# XOR target R=(Rv,Rw) and accumulated contribution into `table` (for the 4-arg trace(A^k B^l)).
function _xor_table_dfs!(ws::_TSMomentWS{P,Cc,Tv,C}, table::Dict{Tuple{Tv,Tv},C},
        i::Int, Rv::Tv, Rw::Tv, rem::Int, depth::Int) where {P,Cc,Tv,C}
    if rem == 0
        if depth > 0
            f = _leaf_contribution(ws, depth)
            key = (Rv, Rw)
            table[key] = get(table, key, zero(C)) + f
        end
        return nothing
    end
    i > ws.n && return nothing
    _xor_table_dfs!(ws, table, i + 1, Rv, Rw, rem, depth)
    @inbounds ws.idx[depth+1] = i
    vi = ws.terms[i].v
    wi = ws.terms[i].w
    @inbounds for t in 1:rem
        ws.mult[depth+1] = t
        if isodd(t)
            _xor_table_dfs!(ws, table, i + 1, Rv ⊻ vi, Rw ⊻ wi, rem - t, depth + 1)
        else
            _xor_table_dfs!(ws, table, i + 1, Rv, Rw, rem - t, depth + 1)
        end
    end
    return nothing
end

# tabulate one operator side: Dict mapping XOR target R => Σ (∏coeff · K_ref · D̃) over its size-`p` multisets
function _xor_table(sterms::Vector{P}, scoeffs::Vector{Cc}, n::Int, reach::Vector{Tv},
        maxsupp::Int, p::Int, ::Type{C}) where {P,Cc,Tv,C}
    ws = _new_ws(sterms, scoeffs, n, reach, maxsupp, max(p, 1), C)
    table = Dict{Tuple{Tv,Tv},C}()
    _xor_table_dfs!(ws, table, 1, zero(Tv), zero(Tv), p, 0)
    return table
end

"""
    trace_moment(o::Operator{<:PauliStringTS}, k::Int; scale=0)

Efficiently compute `trace(o^k)` for a translation-symmetric operator `o` by enumerating
only the Pauli-string multisets whose product is the identity and summing each ordering
phase analytically — without ever constructing `o^(k/2)`.

This is the translation-symmetric counterpart of the moment method. It exploits translation
invariance by anchoring the first factor to a stored representative and multiplying by the
number of translations, so it is typically much faster than [`trace_product`](@ref)`(o, k)`
for higher and odd moments. The result matches `trace_product(o, k)` to machine precision.

If `scale` is not 0, then the result is normalized such that `trace(identity)=scale`.

The keyword `fold` (default `true`) toggles the translation-folding optimization; `fold=false`
runs the equivalent un-folded enumeration over the resummed operator and is used for testing.
`multithreaded` distributes the (folded) search over `Threads.nthreads()` threads; it defaults
to `true` whenever Julia is started with more than one thread. The result is independent of the
number of threads (the partial sums are reduced in a fixed order).
"""
function trace_moment(o::Operator{<:PauliStringTS}, k::Int; scale=0, fold::Bool=true, multithreaded::Bool=Threads.nthreads() > 1)
    k >= 0 || throw(ArgumentError("k must be >= 0, got $k"))
    N = qubitlength(o)
    scl = iszero(scale) ? 2.0^N : float(scale)
    T = scalartype(o)
    k == 0 && return T(scl)

    Ls = qubitsize(o)
    Ps = periodicflags(o)
    num_translations = Base.prod(L for (L, p) in zip(Ls, Ps) if p)

    C = complex(float(T))
    H = resum(o)                      # degree-1 dense operator (cheap: the Hamiltonian)
    isempty(H.strings) && return T(zero(C))   # empty operator: tr(0^k)=0 for k>=1

    sterms, scoeffs, n, reach, maxsupp = _moment_ws_data(H.strings, H.coeffs)
    Tv = typeof(sterms[1].v)
    Cc = eltype(scoeffs)
    cap = max(k, 1)

    if !fold
        # un-folded reference: enumerate identity multisets of size k over the resummed operator
        ws = _new_ws(sterms, scoeffs, n, reach, maxsupp, cap, C)
        _ts_moment_dfs!(ws, 1, zero(Tv), zero(Tv), k, 0)
        return ws.acc * scl
    end

    # folded: anchor the first factor to each stored representative (shift = identity); the
    # remaining k-1 factors form a multiset M' over the full alphabet with XOR = anchor.
    reps = keys(o)
    cf = values(o)
    nanchor = length(reps)

    if multithreaded && k >= 2
        # independent seeds: (anchor a, smallest used term i0) — each runs in its own workspace.
        # Precompute the per-seed anchor data so the threaded loop only calls a function barrier.
        nseed = nanchor * n
        seed_r0v = Vector{Tv}(undef, nseed)
        seed_r0w = Vector{Tv}(undef, nseed)
        seed_c = Vector{Cc}(undef, nseed)
        seed_i0 = Vector{Int}(undef, nseed)
        s = 0
        for a in 1:nanchor
            r0 = representative(reps[a])
            ca = Cc(cf[a])
            for i0 in 1:n
                s += 1
                seed_r0v[s] = r0.v
                seed_r0w[s] = r0.w
                seed_c[s] = ca
                seed_i0[s] = i0
            end
        end
        partials = zeros(C, nseed)
        Threads.@threads for q in 1:nseed
            partials[q] = _ts_seed_partial(sterms, scoeffs, n, reach, maxsupp, cap,
                seed_r0v[q], seed_r0w[q], seed_r0v[q], seed_c[q], seed_i0[q], k - 1, C)
        end
        return sum(partials) * scl * num_translations
    end

    ws = _new_ws(sterms, scoeffs, n, reach, maxsupp, cap, C)
    @inbounds for a in 1:nanchor
        r0 = representative(reps[a])
        ws.tv = r0.v
        ws.tw = r0.w
        ws.anchorv = r0.v
        ws.anchorc = Cc(cf[a])
        _ts_moment_dfs!(ws, 1, zero(Tv), zero(Tv), k - 1, 0)
    end
    return ws.acc * scl * num_translations
end

# number of size-p multisets over m items (BigInt to avoid overflow in the cost heuristic)
_multiset_count(m::Int, p::Int) = binomial(big(m + p - 1), big(p))

"""
    trace_moment(o::Operator, k::Int; scale=0)

Efficiently compute `trace(o^k)` by enumerating only the Pauli-string multisets whose product
is the identity and summing each ordering phase analytically — without ever constructing
`o^(k/2)`. This is the moment method of issue #80 and is typically much faster than
[`trace_product`](@ref)`(o, k)` for higher and odd moments; the result matches it to machine
precision.

If `scale` is not 0, the result is normalized such that `trace(identity)=scale`.
`multithreaded` distributes the search over `Threads.nthreads()` threads (default: on when
Julia is started with more than one thread). The result is independent of the thread count.
"""
function trace_moment(o::Operator{<:PauliString}, k::Int; scale=0, multithreaded::Bool=Threads.nthreads() > 1)
    k >= 0 || throw(ArgumentError("k must be >= 0, got $k"))
    N = qubitlength(o)
    scl = iszero(scale) ? 2.0^N : float(scale)
    T = scalartype(o)
    k == 0 && return T(scl)
    C = complex(float(T))
    isempty(o.strings) && return T(zero(C))   # empty operator: tr(0^k)=0 for k>=1

    sterms, scoeffs, n, reach, maxsupp = _moment_ws_data(o.strings, o.coeffs)
    Tv = typeof(sterms[1].v)
    Cc = eltype(scoeffs)
    cap = max(k, 1)

    if multithreaded && k >= 2
        # partition by the smallest used term i0; each seed runs in its own workspace
        partials = zeros(C, n)
        Threads.@threads for i0 in 1:n
            partials[i0] = _ts_seed_partial(sterms, scoeffs, n, reach, maxsupp, cap,
                zero(Tv), zero(Tv), zero(Tv), one(Cc), i0, k, C)
        end
        return sum(partials) * scl
    end

    ws = _new_ws(sterms, scoeffs, n, reach, maxsupp, cap, C)
    _ts_moment_dfs!(ws, 1, zero(Tv), zero(Tv), k, 0)
    return ws.acc * scl
end

"""
    trace_moment(A::Operator, k::Int, B::Operator, l::Int; scale=0)

Efficiently compute `trace(A^k * B^l)` by the moment method, without constructing the powers.
The cheaper side is tabulated by XOR target `R` (summing its multiset contributions), and the
other side is enumerated for each `R`; the two blocks combine through a `(-1)^{ycount(R)}`
cross phase. Matches [`trace_product`](@ref)`(A, k, B, l)` to machine precision.

If `scale` is not 0, the result is normalized such that `trace(identity)=scale`.
"""
function trace_moment(A::Operator{<:PauliString}, k::Int, B::Operator{<:PauliString}, l::Int; scale=0, multithreaded::Bool=Threads.nthreads() > 1)
    checklength(A, B)
    (k >= 0 && l >= 0) || throw(ArgumentError("k,l must be >= 0, got k=$k, l=$l"))
    N = qubitlength(A)
    scl = iszero(scale) ? 2.0^N : float(scale)
    T = promote_type(scalartype(A), scalartype(B))
    C = complex(float(T))
    k == 0 && return trace_moment(B, l; scale=scale, multithreaded=multithreaded)
    l == 0 && return trace_moment(A, k; scale=scale, multithreaded=multithreaded)
    (isempty(A.strings) || isempty(B.strings)) && return T(zero(C))

    As, Acf, An, Areach, Amax = _moment_ws_data(A.strings, Vector{C}(A.coeffs))
    Bs, Bcf, Bn, Breach, Bmax = _moment_ws_data(B.strings, Vector{C}(B.coeffs))

    # tabulate the cheaper side (fewer size-p multisets); target-DFS the other.
    if _multiset_count(Bn, l) <= _multiset_count(An, k)
        table = _xor_table(Bs, Bcf, Bn, Breach, Bmax, l, C)
        Ls, Lcf, Ln, Lreach, Lmax, lp = As, Acf, An, Areach, Amax, k
    else
        table = _xor_table(As, Acf, An, Areach, Amax, k, C)
        Ls, Lcf, Ln, Lreach, Lmax, lp = Bs, Bcf, Bn, Breach, Bmax, l
    end
    Tv = typeof(Ls[1].v)
    cap = max(lp, 1)

    if multithreaded && lp >= 2 && !isempty(table)
        keys_R = collect(keys(table))
        seeds = Tuple{Int,Int}[]
        for ri in eachindex(keys_R), i0 in 1:Ln
            push!(seeds, (ri, i0))
        end
        partials = zeros(C, length(seeds))
        Threads.@threads for s in eachindex(seeds)
            ri, i0 = seeds[s]
            Rv, Rw = keys_R[ri]
            yR = count_ones(Rv & Rw)
            ac = table[(Rv, Rw)] * (iseven(yR) ? one(C) : -one(C))
            partials[s] = _ts_seed_partial(Ls, Lcf, Ln, Lreach, Lmax, cap,
                Rv, Rw, zero(Tv), ac, i0, lp, C)
        end
        return sum(partials) * scl
    end

    ws = _new_ws(Ls, Lcf, Ln, Lreach, Lmax, cap, C)
    for (R, fR) in table
        Rv, Rw = R
        yR = count_ones(Rv & Rw)
        ws.tv = Rv
        ws.tw = Rw
        ws.anchorv = zero(Tv)
        ws.anchorc = fR * (iseven(yR) ? one(C) : -one(C))
        _ts_moment_dfs!(ws, 1, zero(Tv), zero(Tv), lp, 0)
    end
    return ws.acc * scl
end

# lowest active site of a Pauli string (1-based); identity sorts last
@inline _minsite(p::PauliString) = (u = p.v | p.w; iszero(u) ? typemax(Int) : trailing_zeros(u) + 1)


# Oerations between Operator and PauliString
# ----------------------------------------------------


function trace_product(o::Operator, p::PauliString; scale=0)
    checklength(o, p)
    c = get_coeff(o, p)
    N = qubitlength(o)
    (scale == 0) && (scale = 2.0^N)
    return c * scale
end

trace_product(p::PauliString, o::Operator; scale=0) = trace_product(o, p; scale=scale)



# Operations between OperatorTS and PauliStringTS
# ----------------------------------------------------

function trace_product(o1::Operator{<:PauliStringTS}, o2::PauliStringTS; scale=0)
    checklength(o1, o2)
    Ls = qubitsize(o1)
    Ps = periodicflags(o1)
    tr = zero(scalartype(o1))
    i = findfirst(==(o2), o1.strings)
    isnothing(i) && return tr
    rep1 = representative(o2)
    p, k = prod(rep1, rep1)
    c1 = o1.coeffs[i]
    c2 = (1im)^ycount(o2)
    f = c1 * c2 * k
    for s in all_shifts(Ls, Ps)
        shifted = shift(rep1, Ls, Ps, s)
        if shifted == rep1
            tr += f
        end
    end
    (iszero(scale)) && (scale = 2.0^Base.prod(Ls))
    num_translations = Base.prod(L for (L, p) in zip(Ls, Ps) if p)
    return tr * scale * num_translations
end

trace_product(p::PauliStringTS, o::Operator{<:PauliStringTS}; scale=0) = trace_product(o, p; scale=scale)



# Operations between PauliString and PauliString
# ----------------------------------------------

function trace_product(s1::P, s2::P; scale=0) where {P<:PauliString}
    N = qubitlength(s1)
    if s1 == s2
        (iszero(scale)) && (scale = 2.0^N)
        return scale
    else
        return 0
    end
end

function trace_product(s1::P, s2::P; scale=0) where {P<:PauliStringTS}
    N = qubitlength(s1)
    if s1 == s2
        (scale == 0) && (scale = 2.0^N)
        Ls = qubitsize(s1)
        Ps = periodicflags(s1)
        num_translations = Base.prod(L for (L, p) in zip(Ls, Ps) if p)
        return scale * num_translations
    else
        return 0
    end
end
