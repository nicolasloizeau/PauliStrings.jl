

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
function _trace_product_power_halves(A::AbstractOperator, k::Int; scale=0)
    m = div(k, 2)
    n = k - m
    C = A^m
    (k%2 == 0) && (return trace_product(C; scale=scale))
    D = A^n
    return trace_product(C, D; scale=scale)
end

trace_product(A::AbstractOperator, k::Int; scale=0) = _trace_product_power_halves(A, k; scale=scale)

const _TRACE_PRODUCT_MULTISET_MAX_REACHABLE = 2_000_000

# Enumerate the odd-count support of identity-yielding multisets.  A multiset
# count vector `n` produces the identity iff the XOR of strings with odd `n_i`
# is the identity.  Once those odd supports are fixed, any additional pairs of
# strings preserve the identity, so `_trace_product_multisets` only has to
# distribute `(k - length(odd_support)) ÷ 2` pairs.
function _identity_parity_patterns(strings::Vector{P}, k::Int; maxstates::Int=_TRACE_PRODUCT_MULTISET_MAX_REACHABLE) where {P<:PauliString}
    M = length(strings)
    id = one(P)
    maxodd = min(k, M)
    reachable = [Set{P}() for _ in 1:(M + 1), _ in 1:(maxodd + 1)]
    push!(reachable[M + 1, 1], id)
    totalstates = 1

    for i in M:-1:1
        p = strings[i]
        for r in 0:maxodd
            dst = reachable[i, r + 1]
            union!(dst, reachable[i + 1, r + 1])
            if r > 0
                for q in reachable[i + 1, r]
                    push!(dst, p ⊻ q)
                end
            end
            totalstates += length(dst)
            totalstates > maxstates && return nothing
        end
    end

    patterns = Vector{Vector{Int}}()
    current = Int[]

    function visit(i::Int, r::Int, target::P)
        if r == 0
            target == id && push!(patterns, copy(current))
            return nothing
        end
        i > M && return nothing
        r > M - i + 1 && return nothing

        if target in reachable[i + 1, r + 1]
            visit(i + 1, r, target)
        end

        nexttarget = target ⊻ strings[i]
        if nexttarget in reachable[i + 1, r]
            push!(current, i)
            visit(i + 1, r - 1, nexttarget)
            pop!(current)
        end
        return nothing
    end

    for r in (k & 1):2:maxodd
        visit(1, r, id)
        length(patterns) > maxstates && return nothing
    end

    return patterns
end

function _multiset_phase_sum(strings::Vector{P}, indices::Vector{Int}, counts0::Vector{Int}) where {P<:PauliString}
    isempty(indices) && return 1
    counts = copy(counts0)
    memo = Dict{Tuple{Vararg{Int}},Int}()

    function visit(current_v::Unsigned)
        key = Tuple(counts)
        cached = get(memo, key, nothing)
        isnothing(cached) || return cached

        total = 0
        anyleft = false
        for a in eachindex(indices)
            counts[a] == 0 && continue
            anyleft = true
            p = strings[indices[a]]
            counts[a] -= 1
            total += _prod_phase_sign(current_v, p.w) * visit(current_v ⊻ p.v)
            counts[a] += 1
        end
        anyleft || (total = 1)
        memo[key] = total
        return total
    end

    return visit(zero(typeof(strings[1].v)))
end

function _trace_product_multisets(A::Operator{P,T}, k::Int; scale=0) where {P<:PauliString,T}
    k < 0 && throw(DomainError(k, "power must be non-negative"))
    N = qubitlength(A)
    trscale = iszero(scale) ? 2.0^N : scale
    k == 0 && return one(T) * trscale
    isempty(A.strings) && return zero(T) * trscale

    strings = A.strings
    coeffs = A.coeffs
    length(strings) == length(coeffs) || throw(DimensionMismatch("strings and coefficients must have the same length"))

    parity_patterns = _identity_parity_patterns(strings, k)
    isnothing(parity_patterns) && return nothing

    M = length(strings)
    coeffpowers = [[coeffs[i]^m for m in 0:k] for i in 1:M]
    odd = falses(M)
    indices = Int[]
    counts = Int[]
    moment = Ref(zero(T))
    multiset_count = Ref(0)
    aborted = Ref(false)

    function visit_counts(i::Int, remaining_pairs::Int, coeffprod)
        aborted[] && return nothing
        if i > M
            remaining_pairs == 0 || return nothing
            multiset_count[] += 1
            if multiset_count[] > _TRACE_PRODUCT_MULTISET_MAX_REACHABLE
                aborted[] = true
                return nothing
            end
            phase_sum = _multiset_phase_sum(strings, indices, counts)
            moment[] += coeffprod * phase_sum
            return nothing
        end

        basecount = odd[i] ? 1 : 0
        for pairs in 0:remaining_pairs
            count = basecount + 2pairs
            if count == 0
                visit_counts(i + 1, remaining_pairs - pairs, coeffprod)
            else
                push!(indices, i)
                push!(counts, count)
                visit_counts(i + 1, remaining_pairs - pairs, coeffprod * coeffpowers[i][count + 1])
                pop!(counts)
                pop!(indices)
            end
            aborted[] && return nothing
        end
        return nothing
    end

    for pattern in parity_patterns
        fill!(odd, false)
        for i in pattern
            odd[i] = true
        end
        pairs = (k - length(pattern)) ÷ 2
        visit_counts(1, pairs, one(T))
        aborted[] && return nothing
    end

    return moment[] * trscale
end

function trace_product(A::Operator{P}, k::Int; scale=0) where {P<:PauliString}
    k < 0 && throw(DomainError(k, "power must be non-negative"))
    k <= 2 && return _trace_product_power_halves(A, k; scale=scale)
    moment = _trace_product_multisets(A, k; scale=scale)
    isnothing(moment) && return _trace_product_power_halves(A, k; scale=scale)
    return moment
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
