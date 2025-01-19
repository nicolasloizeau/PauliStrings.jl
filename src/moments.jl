

using Base.Iterators


"""
    trace_product(o1::Operator, o2::Operator; scale=0)
    trace_product(o1::OperatorTS1D, o2::OperatorTS1D; scale=0)

Efficiently compute `trace(o1*o2)`. This is much faster than doing `trace(o1*o2)`.
If `scale` is not 0, then the result is normalized such that trace(identity)=scale.
"""
function trace_product(o1::Operator, o2::Operator; scale=0)
    (scale == 0) && (scale = 2^o1.N)
    tr = 0
    T = uinttype(o1)
    o1v::Vector{T} = o1.v
    o2v::Vector{T} = o2.v
    o1w::Vector{T} = o1.w
    o2w::Vector{T} = o2.w
    for i in eachindex(o1v)
        o1vi = o1v[i]
        for j in eachindex(o2v)
            if (o1vi == o2v[j] && o1w[i] == o2w[j])
                tr += o1.coef[i] * o2.coef[j] * (-1)^count_ones(o1v[i] & o2w[j])
            end
        end
    end
    return tr * scale
end
function trace_product(o1::OperatorTS1D, o2::OperatorTS1D; scale=0)
    (scale == 0) && (scale = 2^o1.N)
    tr = 0
    N = o1.N
    T = uinttype(o1)
    o1v::Vector{T} = o1.v
    o2v::Vector{T} = o2.v
    o1w::Vector{T} = o1.w
    o2w::Vector{T} = o2.w
    for k in 0:N-1
        for j in eachindex(o2v)
            o2vj = o2v[j]
            o2wj = o2w[j]
            v2 = rotate_lower(o2vj, N, k)
            w2 = rotate_lower(o2wj, N, k)
            for i in eachindex(o1v)
                if (o1v[i] == v2) && (o1w[i] == w2)
                    tr += o1.coef[i] * o2.coef[j] * (-1)^count_ones(o1v[i] & w2)
                end
            end
        end
    end
    return tr * scale * N
end


"""
    oppow(o::Operator, k::Int)
    oppow(o::OperatorTS1D, k::Int)

kth power of o. Same as `^`.
"""
function oppow(o::Operator, k::Int)
    # divide and conqueer is not faster
    r = eye(o.N)
    for i in 1:k
        r *= o
    end
    return r
end
function oppow(o::OperatorTS1D, k::Int)
    # divide and conqueer is not faster
    r = OperatorTS1D(eye(o.N))
    for i in 1:k
        r *= o
    end
    return r
end


"""
    Base.:^(o::Operator, k::Int)

kth power of o. Same as `oppow`.
"""
Base.:^(o::Operator, k::Int) = oppow(o::Operator, k::Int)
Base.:^(o::OperatorTS1D, k::Int) = oppow(o::OperatorTS1D, k::Int)

"""
    trace_product(A::Operator, k::Int, B::Operator, l::Int; scale=0)
    trace_product(A::OperatorTS1D, k::Int, B::OperatorTS1D, l::Int; scale=0)

Efficiently compute `trace(A^k*B^l)`. This is much faster than doing `trace(A^k*B^l)`.

If `scale` is not 0, then the result is normalized such that trace(identity)=scale.
"""
function trace_product(A::Operator, k::Int, B::Operator, l::Int; scale=0)
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
function trace_product(A::Operator, k::Int; scale=0)
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
