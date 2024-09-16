

using Base.Iterators


function trace_product(o1::Operator, o2::Operator)
    tr = 0
    for i in 1:length(o1.v)
        for j in 1:length(o2.v)
            v = o1.v[i] ⊻ o2.v[j]
            w = o1.w[i] ⊻ o2.w[j]
            if (v==0 && w==0)
                tr += o1.coef[i] * o2.coef[j]* (-1)^count_ones(o1.v[i] & o2.w[j])
            end
        end
    end
    return tr*2^o1.N
end
"""
Power of an operator
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
compute trace(A^k*B^l)
"""
function trace_product(A::Operator, k::Int, B::Operator, l::Int)
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
    return trace_product(C,D)
end

function trace_product(A::Operator, k::Int)
    m = div(k, 2)
    n = k-m
    C = oppow(A, m)
    D = oppow(A, n)
    return trace_product(C,D)
end



function multrings(vs, ws, s)
    v = 0
    w = 0
    c = 0
    for i in 1:s
        c += count_ones(v & ws[i])
        v = v ⊻ vs[i]
        w = w ⊻ ws[i]
    end
    return (-1)^c, v, w
end

"""
Trace(exp(-beta H))
"""
function trace_exp(H::Operator, beta::Real, n::Int)
    r = 0
    for k in 1:n
        r += (-1)^k*beta^k*trace_product(H, k)/factorial(k)
    end
    return r
end

"""
Trace(O*exp(-beta H))
"""
function trace_exp(H::Operator, O::Operator, beta::Real, n::Int)
    r = 0
    for k in 1:n
        r += (-1)^k*beta^k*trace_product(O, 1, H, k)/factorial(k)
    end
    return r
end


# function trace_product(os::NTuple{3, Operator})
#     tr = 0
#     l = 3
#     for i in 1:length(os[1])
#     for j in 1:length(os[2])
#     for k in 1:length(os[3])
#     vs = (os[1].v[i], os[2].v[j], os[3].v[k])
#     ws = (os[1].w[i], os[2].w[j], os[3].w[k])
#     c, v, w = multrings(vs, ws, l)
#     if (v==0 && w==0)
#     tr += c*os[1].coef[i]*os[2].coef[j]*os[3].coef[k]
#     end
#     end
#     end
#     end
#     return tr*2^os[1].N
# end

function trace_product(os::NTuple{3, Operator})
    tr = 0
    l = 3
    for i in 1:length(os[1])
    for j in 1:length(os[2])
    for k in 1:length(os[3])
    vs = (os[1].v[i], os[2].v[j], os[3].v[k])
    ws = (os[1].w[i], os[2].w[j], os[3].w[k])
    c, v, w = multrings(vs, ws, l)
    if (v==0 && w==0)
    tr += c*os[1].coef[i]*os[2].coef[j]*os[3].coef[k]
    end
    end
    end
    end
    return tr*2^os[1].N
end


function trace_product(os::NTuple{4, Operator})
    (o1, o2, o3, o4) = os
    tr = 0
    h = 4
    for i in eachindex(o1.v, o1.w)
    for j in eachindex(o2.v, o2.w)
    for k in eachindex(o3.v, o3.w)
    for l in eachindex(o4.v, o4.w)
    vs = (o1.v[i], o2.v[j], o3.v[k], o4.v[l])
    ws = (o1.w[i], o2.w[j], o3.w[k], o4.w[l])
    c, v, w = multrings(vs, ws, h)
    if (v==0 && w==0)
    tr += c*o1.coef[i]*o2.coef[j]*o3.coef[k]*o4.coef[l]
    end
    end
    end
    end
    end
    return tr*2^o1.N
end


function trace(args)
    body
end


# function trace_product(os::{4, Operator})
# function trace_product(os::Vararg{Operator,L}) where {L}
#     tr = 0
#     p = [length(o) for o in os]
#     ranges = [1:i for i in p]
#     vs = [os[i].v for i in 1:L]
#     ws = [os[i].w for i in 1:L]
#     productv = Iterators.product(eachindex.(vs)...)
#     productw = Iterators.product(eachindex.(ws)...)
#     for (indsv,indsw) in zip(productv, productw)
#
#         v = 0
#         w = 0
#         c = 0
#         for i in eachindex(os, indsv)
#             c += count_ones(v & os[i].w[indsw[i]])
#             v = v ⊻ vs[i][indsv[i]]
#             w = w ⊻ ws[i][indsw[i]]
#         end
#         c = (-1)^c
#
#         if (v==0 && w==0)
#             coefs = [os[i].coef[indsv[i]] for i in 1:L]
#             c = reduce(*, coefs) * c
#             tr += c
#         end
#     end
#     return tr*2^os[1].N
# end

# # function trace_product(os::{4, Operator})
# function trace_product(os::Vararg{Operator,L}) where {L}
#     tr = 0
#     p = [length(o) for o in os]
#     ranges = [1:i for i in p]
#     vs = [os[i].v for i in 1:L]
#     ws = [os[i].w for i in 1:L]
#     # productv = Iterators.product(eachindex.(vs)...)
#     # productw = Iterators.product(eachindex.(ws)...)
#     for inds in Iterators.product(eachindex.(vs, ws)...)
#
#         vs = [os[i].v[inds[i]] for i in eachindex(os, inds)]
#         ws = [os[i].w[inds[i]] for i in eachindex(os, inds)]
#
#
#
#
#
#
#         c, v, w = multrings(vs, ws)
#         if (v==0 && w==0)
#             coefs = [os[i].coef[inds[i]] for i in 1:L]
#             c = reduce(*, coefs) * c
#             tr += c
#         end
#     end
#     return tr*2^os[1].N
# end
