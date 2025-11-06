

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
        for s in all_shifts(Ls)
            shifted = shift(rep1, Ls, s)
            if shifted == rep1
                tr += f
            end
        end
    end
    (scale == 0) && (scale = 2.0^Base.prod(Ls))
    return tr * scale * Base.prod(Ls)
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
    trace_product(A::AbstractOperator, k::Int; scale=0)

Efficiently compute `trace(A^k)`. This is much faster than doing `trace(A^k)`.

If `scale` is not 0, then the result is normalized such that trace(identity)=scale.
"""
function trace_product(A::AbstractOperator, k::Int; scale=0)
    m = div(k, 2)
    n = k - m
    C = A^m
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
    tr = zero(scalartype(o1))
    i = findfirst(==(o2), o1.strings)
    isnothing(i) && return tr
    rep1 = representative(o2)
    p, k = prod(rep1, rep1)
    c1 = o1.coeffs[i]
    c2 = (1im)^ycount(o2)
    f = c1 * c2 * k
    for s in all_shifts(Ls)
        shifted = shift(rep1, Ls, s)
        if shifted == rep1
            tr += f
        end
    end
    (scale == 0) && (scale = 2.0^Base.prod(Ls))
    return tr * scale * Base.prod(Ls)
end

trace_product(p::PauliStringTS, o::Operator{<:PauliStringTS}; scale=0) = trace_product(o, p; scale=scale)
