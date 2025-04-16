

using Base.Iterators


"""
    trace_product(o1::Operator, o2::Operator; scale=0)
    trace_product(o1::OperatorTS1D, o2::OperatorTS1D; scale=0)

Efficiently compute `trace(o1*o2)`. This is much faster than doing `trace(o1*o2)`.
If `scale` is not 0, then the result is normalized such that trace(identity)=scale.
"""
function trace_product(o1::Operator, o2::Operator; scale=0)
    checklength(o1, o2)
    N = qubitlength(o1)
    (scale == 0) && (scale = 2.0^N)
    tr = zero(scalartype(o1))

    # ensure `@inbounds` is safe
    length(o1.strings) == length(o1.coeffs) || throw(DimensionMismatch("strings and coefficients must have the same length"))
    length(o2.strings) == length(o2.coeffs) || throw(DimensionMismatch("strings and coefficients must have the same length"))

    @inbounds for i in eachindex(o1.strings)
        p1, c1 = o1.strings[i], o1.coeffs[i]
        for j in eachindex(o2.strings)
            p2, c2 = o2.strings[j], o2.coeffs[j]

            p, k = prod(p1, p2)
            if isone(p)
                tr += c1 * c2 * k
            end
        end
    end
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

    @inbounds for i in eachindex(o1.strings)
        p1, c1 = o1.strings[i], o1.coeffs[i]
        for j in eachindex(o2.strings)
            p2, c2 = o2.strings[j], o2.coeffs[j]

            for n in 0:(N-1)
                p3 = rotate_lower(p2, n)
                p, k = prod(p1, p3)
                if isone(p)
                    tr += c1 * c2 * k
                end
            end
        end
    end
    return tr * scale * N
end

@doc """
    oppow(o::Operator, k::Int)
    oppow(o::OperatorTS1D, k::Int)

kth power of o. Same as `^`.
""" oppow
oppow(o::AbstractOperator, k::Int) = o^k

"""
    Base.:^(o::Operator, k::Int)

kth power of o. Same as `oppow`.
"""
Base.:^(o::AbstractOperator, k::Int) = Base.power_by_squaring(o, k)

"""
    trace_product(A::Operator, k::Int, B::Operator, l::Int; scale=0)
    trace_product(A::OperatorTS1D, k::Int, B::OperatorTS1D, l::Int; scale=0)

Efficiently compute `trace(A^k*B^l)`. This is much faster than doing `trace(A^k*B^l)`.

If `scale` is not 0, then the result is normalized such that trace(identity)=scale.
"""
function trace_product(A::AbstractOperator, k::Int, B::AbstractOperator, l::Int; scale=0)
    @assert typeof(A) == typeof(B)
    m = div(k + l, 2)
    n = k + l - m
    if k < m
        C = oppow(A, k) * oppow(B, m - k)
        D = oppow(B, n)
    elseif k > m
        C = oppow(A, m)
        D = oppow(A, k - m) * oppow(B, l)
    else
        C = oppow(A, k)
        D = oppow(B, l)
    end
    return trace_product(C, D; scale=scale)
end

"""
    trace_product(A::Operator, k::Int; scale=0)

Efficiently compute `trace(A^k)`. This is much faster than doing `trace(A^k)`.

If `scale` is not 0, then the result is normalized such that trace(identity)=scale.
"""
function trace_product(A::AbstractOperator, k::Int; scale=0)
    m = div(k, 2)
    n = k - m
    C = oppow(A, m)
    D = oppow(A, n)
    return trace_product(C, D; scale=scale)
end


"""
    moments(H::Operator, kmax::Int; start=1, scale=0)
    moments(H::OperatorTS1D, kmax::Int; start=1, scale=0)

Compute the first kmax moments of H.
start is the first moment to start from.

If scale is not 0, then the result is normalized such that trace(identity)=scale.
"""
function moments(H::Operator, kmax::Int; start=1, scale=0)
    return [trace_product(H, k; scale=scale) for k in start:kmax]
end
