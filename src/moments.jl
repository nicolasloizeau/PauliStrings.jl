

using Base.Iterators


"""
    trace_product(o1::Operator, o2::Operator; scale=0)

Efficiently compute `trace(o1*o2)`. This is much faster than doing `trace(o1*o2)`.
If `scale` is not 0, then the result is normalized such that trace(identity)=scale.
"""
function trace_product(o1::Operator, o2::Operator; scale=0)
    (scale==0)&&(scale=2^o1.N)
    tr = 0
    o1v::Vector{Int} = o1.v
    o2v::Vector{Int} = o2.v
    o1w::Vector{Int} = o1.w
    o2w::Vector{Int} = o2.w
    for i in eachindex(o1v)
        o1vi = o1v[i]
        for j in eachindex(o2v)
            v = o1vi ⊻ o2v[j]
            if (v==0)
                w = o1w[i] ⊻ o2w[j]
                if (w==0)
                    tr += o1.coef[i] * o2.coef[j]* (-1)^count_ones(o1v[i] & o2w[j])
                end
            end
        end
    end
    return tr*scale
end


"""
    oppow(o::Operator, k::Int)

kth power of o. Same as `^`.
"""
function oppow(o::Operator, k::Int)
    # divide and conqueer is not faster
    r = eye(o.N)
    for i in 1:k
        r*=o
    end
    return r
end

"""
    Base.:^(o::Operator, k::Int)

kth power of o. Same as `oppow`.
"""
Base.:^(o::Operator, k::Int) = oppow(o::Operator, k::Int)


"""
    trace_product(A::Operator, k::Int, B::Operator, l::Int; scale=0)

Efficiently compute `trace(A^k*B^l)`. This is much faster than doing `trace(A^k*B^l)`.

If `scale` is not 0, then the result is normalized such that trace(identity)=scale.
"""
function trace_product(A::Operator, k::Int, B::Operator, l::Int; scale=0)
    m = div(k+l, 2)
    n = k+l-m
    if k<m
        C = oppow(A, k)*oppow(B, m-k)
        D = oppow(B, n)
    elseif k>m
        C = oppow(A, m)
        D = oppow(A, k-m)*oppow(B, l)
    else
        C = oppow(A, k)
        D = oppow(B, l)
    end
    return trace_product(C,D; scale=scale)
end

"""
    trace_product(A::Operator, k::Int; scale=0)

Efficiently compute `trace(A^k)`. This is much faster than doing `trace(A^k)`.

If `scale` is not 0, then the result is normalized such that trace(identity)=scale.
"""
function trace_product(A::Operator, k::Int; scale=0)
    m = div(k, 2)
    n = k-m
    C = oppow(A, m)
    D = oppow(A, n)
    return trace_product(C,D; scale=scale)
end


"""
    moments(H::Operator, kmax::Int; start=1, scale=0)

Compute the first kmax moments of H.
start is the first moment to start from.

If scale is not 0, then the result is normalized such that trace(identity)=scale.
"""
function moments(H::Operator, kmax::Int; start=1, scale=0)
    mus = []
    for k in start:kmax
        push!(mus, trace_product(H, k;scale=scale))
    end
    return mus
end