

using Base.Iterators


"""
    trace_product(o1::Operator, o2::Operator; scale=0)
    trace_product(o1::OperatorTS1D, o2::OperatorTS1D; scale=0)

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
function trace_product(o1::OperatorTS1D, o2::OperatorTS1D; scale=0)
    checklength(o1, o2)
    N = qubitlength(o1)
    (scale == 0) && (scale = 2.0^N)
    tr = zero(scalartype(o1))

    # ensure `@inbounds` is safe
    length(o1.strings) == length(o1.coeffs) || throw(DimensionMismatch("strings and coefficients must have the same length"))
    length(o2.strings) == length(o2.coeffs) || throw(DimensionMismatch("strings and coefficients must have the same length"))

    d = emptydict(o2)
    @inbounds for i in eachindex(o2.strings)
        insert!(d, o2.strings[i], o2.coeffs[i])
    end

    @inbounds for i in eachindex(o1.strings)
        p1, c1 = o1.strings[i], o1.coeffs[i]
        for n in 0:(N-1)
            p3 = rotate_lower(p1, n)
            c2 = get(d, p3, nothing)
            isnothing(c2) && continue
            p, k = prod(p3, p3)
            tr += c1 * c2 * k
        end
    end
    (scale == 0) && (scale = 2.0^N)
    return tr * scale*N
end


function trace_product(o1::OperatorTS2D, o2::OperatorTS2D; scale=0)
    checklength(o1, o2)
    N = qubitlength(o1)
    L1 = extent(o1)
    L2 = N ÷ L1
    (scale == 0) && (scale = 2.0^N)
    tr = zero(scalartype(o1))

    # ensure `@inbounds` is safe
    length(o1.strings) == length(o1.coeffs) || throw(DimensionMismatch("strings and coefficients must have the same length"))
    length(o2.strings) == length(o2.coeffs) || throw(DimensionMismatch("strings and coefficients must have the same length"))

    @inbounds for i1 in eachindex(o1.strings)
        p1, c1 = o1.strings[i1], o1.coeffs[i1]
        for i2 in eachindex(o2.strings)
            p2, c2 = o2.strings[i2], o2.coeffs[i2]
            for i in 0:L1-1
                for j in 0:L2-1
                    p3 = rotate_lower(p2, i, j, L1)
                    p, k = prod(p1, p3)
                    if isone(p)
                        tr += c1 * c2 * k
                    end
                end
            end
        end
    end
    return tr * scale * N
end

Base.@deprecate oppow(o::AbstractOperator, k::Int) o^k

"""
    Base.:^(o::Operator, k::Int)

kth power of o.
"""
Base.:^(o::AbstractOperator, k::Int) = Base.power_by_squaring(o, k)

"""
    trace_product(A::Operator, k::Int, B::Operator, l::Int; scale=0)
    trace_product(A::OperatorTS1D, k::Int, B::OperatorTS1D, l::Int; scale=0)
    trace_product(A::OperatorTS2D, k::Int, B::OperatorTS2D, l::Int; scale=0)

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
