

using Base.Iterators


"""
    trace_product(o1::Operator, o2::Operator; scale=0)
    trace_product(o1::OperatorTS, o2::OperatorTS; scale=0)

Efficiently compute `trace(o1*o2)`. This is much faster than doing `trace(o1*o2)`.
If `scale` is not 0, then the result is normalized such that trace(identity)=scale.
"""
function trace_product(o1::Operator, o2::Operator; scale=0)
    # operation is symmetric but more efficient if o1 is the largest collection
    (length(o1.strings) < length(o2.strings)) && return trace_product(o2, o1; scale)

    checklength(o1, o2)
    N = qubitlength(o1)
    tr = zero(scalartype(o1))

    # ensure `@inbounds` is safe
    length(o1.strings) == length(o1.coeffs) || throw(DimensionMismatch("strings and coefficients must have the same length"))
    length(o2.strings) == length(o2.coeffs) || throw(DimensionMismatch("strings and coefficients must have the same length"))

    # trace of product contributes only if product is 1, which only happens when strings are equal
    # this amounts to `indexin`, which we hijack/reimplement here for efficiency
    d = emptydict(o2)
    @inbounds for i in eachindex(o2.strings)
        insert!(d, o2.strings[i], o2.coeffs[i])
    end

    @inbounds for i in eachindex(o1.strings)
        p1, c1 = o1.strings[i], o1.coeffs[i]
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
    for (p2, c2) in zip(o2.strings, o2.coeffs)
        insert!(d, p2, c2)
    end

    for (p1, c1) in zip(o1.strings, o1.coeffs)
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

# Moments via identity-yielding multiset enumeration
# ===================================================
#
# For `H = sum_l c_l P_l`, expanding `tr(H^k) = sum_{l_1...l_k} c_{l_1}...c_{l_k} tr(P_{l_1}...P_{l_k})`,
# a term contributes only when the (order-independent) product of Pauli strings is proportional to
# the identity, i.e. when `P_{l_1} XOR ... XOR P_{l_k} = 1`. The trace then depends on the *multiset*
# of Pauli strings, not their order. Grouping ordered tuples by their multiset `M` (with
# multiplicities `m_i` and distinct stored coefficients `a_i`) gives
#
#   tr(H^k) = scale * sum_{M : XOR(M)=1} prod_i a_i^{m_i} * K_ref(M) * Dtilde(M),
#
# where `K_ref(M)` is the +-1 phase of the product taken in one fixed (canonical) order and
# `Dtilde(M) = sum over distinct orderings of (-1)^(# anticommuting out-of-order pairs)` only
# depends on the anticommutation graph of the distinct terms. Reordering two Pauli strings flips
# the phase by the sign of their (anti)commutation, which is why the whole order dependence factors
# into `K_ref * Dtilde`. This removes the `k!` factor of the naive sum.
#
# `Dtilde(M)` is evaluated cheaply:
#  * every term that commutes with all other terms of `M` factors out as a binomial coefficient,
#    leaving only an "anticommuting core";
#  * the core is summed with the multiplicity dynamic program
#       Dtilde(n) = sum_{j: n_j>0} (-1)^{sum_{c>j, anticommute(j,c)} n_c} * Dtilde(n - e_j),
#    which costs `prod_{core} (m_j + 1)` (usually tiny for local operators).

# +-1 sign of prod(a, b): the phase of P_a * P_b = ksign * P_{a XOR b}
@inline _ksign(av::Unsigned, bw::Unsigned) = 1 - ((count_ones(av & bw) & 1) << 1)

# whether two Pauli strings anticommute
@inline _anticommute(a::PauliString, b::PauliString) =
    isodd(count_ones(a.v & b.w) + count_ones(a.w & b.v))

# C(n, r) as an exact integer (the intermediate phase quantities are integers; this keeps the
# result exact for any `k` up to ~20, beyond which the moment itself overflows Float64 anyway)
function _binomial_int(n::Int, r::Int)
    (r < 0 || r > n) && return zero(Int64)
    r = min(r, n - r)
    b = one(Int64)
    for i in 1:r
        b = (b * (n - r + i)) ÷ i
    end
    return b
end

# lowest active site of a Pauli string (1-based), or a large sentinel for the identity
@inline function _minsite(p::PauliString)
    u = p.v | p.w
    return iszero(u) ? typemax(Int) : trailing_zeros(u) + 1
end

# Reusable workspace for the moment enumeration, holding the (site-sorted) operator data,
# pruning helpers, the current multiset stack, and preallocated buffers for the phase program.
mutable struct _MomentWS{P,Cc,T,C}
    strings::Vector{P}
    coeffs::Vector{Cc}
    reach::Vector{T}     # reach[i] = union of supports of strings i:n
    maxsup::Int          # max support (number of non-identity sites) of any term
    n::Int
    idx::Vector{Int}     # distinct term indices in the current multiset
    mult::Vector{Int}    # their multiplicities
    depth::Int           # number of distinct terms currently on the stack
    keep::Vector{Bool}   # phase-program scratch: still in the core after peeling
    core::Vector{Int}    # core color positions
    strides::Vector{Int}
    digits::Vector{Int}
    cbuf::Vector{Int64}
    acc::Base.RefValue{C}
end

# +-1 sign of the product of the multiset taken in its canonical (stacked) order
function _canonical_sign(ws::_MomentWS{P,Cc,T,C}, r::Int) where {P,Cc,T,C}
    Sv = zero(T)
    K = 1
    @inbounds for j in 1:r
        q = ws.strings[ws.idx[j]]
        for _ in 1:ws.mult[j]
            K *= 1 - ((count_ones(Sv & q.w) & 1) << 1)
            Sv ⊻= q.v
        end
    end
    return K
end

# Signed sum over distinct orderings, Dtilde(M), via the anticommuting-core dynamic program.
function _dtilde!(ws::_MomentWS{P,Cc,T,C}, r::Int) where {P,Cc,T,C}
    # Peel off every term that commutes with all remaining terms; each contributes a binomial.
    @inbounds for j in 1:r
        ws.keep[j] = true
    end
    factor = one(Int64)
    ktot = 0
    @inbounds for j in 1:r
        ktot += ws.mult[j]
    end
    changed = true
    @inbounds while changed
        changed = false
        for j in 1:r
            ws.keep[j] || continue
            commutes_all = true
            qj = ws.strings[ws.idx[j]]
            for c in 1:r
                if c != j && ws.keep[c] && _anticommute(qj, ws.strings[ws.idx[c]])
                    commutes_all = false
                    break
                end
            end
            if commutes_all
                factor *= _binomial_int(ktot, ws.mult[j])
                ktot -= ws.mult[j]
                ws.keep[j] = false
                changed = true
            end
        end
    end

    # collect the anticommuting core
    s = 0
    @inbounds for j in 1:r
        if ws.keep[j]
            s += 1
            ws.core[s] = j
        end
    end
    s == 0 && return factor

    # multiplicity dynamic program over the core (canonical order = stacked order)
    total = 1
    @inbounds for a in 1:s
        ws.strides[a] = total
        total *= (ws.mult[ws.core[a]] + 1)
    end
    length(ws.cbuf) < total && resize!(ws.cbuf, total)
    c = ws.cbuf
    c[1] = one(Int64)
    @inbounds for code in 1:total-1
        tmp = code
        for a in 1:s
            mj1 = ws.mult[ws.core[a]] + 1
            ws.digits[a] = tmp % mj1
            tmp = tmp ÷ mj1
        end
        val = zero(Int64)
        for a in 1:s
            da = ws.digits[a]
            if da > 0
                qa = ws.strings[ws.idx[ws.core[a]]]
                e = 0
                for b in a+1:s
                    if _anticommute(qa, ws.strings[ws.idx[ws.core[b]]])
                        e += ws.digits[b]
                    end
                end
                sgn = iseven(e) ? one(Int64) : -one(Int64)
                val += sgn * c[code+1-ws.strides[a]]
            end
        end
        c[code+1] = val
    end
    return factor * c[total]
end

# accumulate the contribution of one identity-yielding multiset
function _moment_leaf!(ws::_MomentWS{P,Cc,T,C}) where {P,Cc,T,C}
    r = ws.depth
    w = one(Cc)
    @inbounds for j in 1:r
        w *= ws.coeffs[ws.idx[j]]^ws.mult[j]
    end
    kref = _canonical_sign(ws, r)
    dt = _dtilde!(ws, r)
    ws.acc[] += w * (kref * dt)
    return nothing
end

# depth-first enumeration of size-`rem` multisets of terms `i:n` whose running XOR can still
# return to the identity. Two prunes keep the search close to the number of valid multisets:
#  (1) the remaining terms i:n can only touch sites in reach[i];
#  (2) clearing the s active sites of R needs at least ceil(s / maxsup) more terms.
function _moment_dfs!(ws::_MomentWS{P,Cc,T,C}, i::Int, Rv::T, Rw::T, rem::Int) where {P,Cc,T,C}
    if rem == 0
        (iszero(Rv) && iszero(Rw)) && _moment_leaf!(ws)
        return nothing
    end
    i > ws.n && return nothing
    @inbounds reach_i = ws.reach[i]
    !iszero((Rv | Rw) & ~reach_i) && return nothing
    rem * ws.maxsup < count_ones(Rv | Rw) && return nothing
    @inbounds s = ws.strings[i]
    # use term i zero times
    _moment_dfs!(ws, i + 1, Rv, Rw, rem)
    # use term i one or more times
    ws.depth += 1
    d = ws.depth
    @inbounds ws.idx[d] = i
    Rvt = Rv
    Rwt = Rw
    @inbounds for t in 1:rem
        Rvt ⊻= s.v
        Rwt ⊻= s.w
        ws.mult[d] = t
        _moment_dfs!(ws, i + 1, Rvt, Rwt, rem - t)
    end
    ws.depth -= 1
    return nothing
end

"""
    trace_moment(o::Operator, k::Int; scale=0)

Efficiently compute `trace(o^k)` by enumerating only the multisets of Pauli strings whose product
is proportional to the identity and summing their phase contributions analytically.

Unlike `trace_product(o, k)`, this never constructs `o^(k/2)`, so for sparse (e.g. local) operators
it is typically much faster at high moments `k` and uses very little memory.

If `scale` is not 0, then the result is normalized such that `trace(identity)=scale`.

For very large `k` (beyond ~20) the integer phase bookkeeping can overflow; such high moments are
rarely useful since the value itself exceeds `Float64` range.

# Example
```julia
H = Operator(4)
H += "Z", 1, "Z", 2
H += 0.5, "X", 1
trace_moment(H, 4) ≈ trace_product(H, 4)   # true
```
"""
function trace_moment(o::Operator, k::Int; scale=0)
    k < 0 && throw(ArgumentError("k must be non-negative"))
    N = qubitlength(o)
    scaleval = iszero(scale) ? 2.0^N : scale
    C = scalartype(o)
    k == 0 && return one(C) * scaleval
    n = length(o.strings)
    n == 0 && return zero(C) * scaleval
    length(o.strings) == length(o.coeffs) || throw(DimensionMismatch("strings and coefficients must have the same length"))
    P = eltype(o.strings)
    T = typeof(o.strings[1].v)
    Cc = eltype(o.coeffs)

    # Sort terms by lowest active site so the reach prune forces low sites to be resolved early.
    perm = sortperm(o.strings; by=_minsite)
    strings = o.strings[perm]
    coeffs = o.coeffs[perm]

    # suffix support and max single-term support
    reach = Vector{T}(undef, n + 1)
    reach[n+1] = zero(T)
    maxsup = 0
    @inbounds for i in n:-1:1
        reach[i] = reach[i+1] | (strings[i].v | strings[i].w)
        maxsup = max(maxsup, count_ones(strings[i].v | strings[i].w))
    end
    maxsup = max(maxsup, 1)

    ws = _MomentWS{P,Cc,T,C}(strings, coeffs, reach, maxsup, n,
        Vector{Int}(undef, k), Vector{Int}(undef, k), 0,
        Vector{Bool}(undef, k), Vector{Int}(undef, k), Vector{Int}(undef, k), Vector{Int}(undef, k),
        Vector{Int64}(undef, 2), Base.RefValue{C}(zero(C)))

    _moment_dfs!(ws, 1, zero(T), zero(T), k)
    return ws.acc[] * scaleval
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
