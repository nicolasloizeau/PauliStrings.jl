
using Dictionaries


"""
    OperatorTS1D(o::Operator; full=true)

Initialize a 1D translation invariant operator from an Operator
\$O=\\sum_i o_i O_i\$ where \$O_i=T_i(O_0)\$ and \$T_i\$ is the i-sites translation operator.
Set full=true if passing \$O\$, an Operator that is supported on the whole chain (i.e converting from a translation symmetric [`Operator`](@ref))
Set full=false if passing \$O_0\$,a local term o such that the full operator is \$O=\\sum_i o_i T_i(O_0)\$
"""
function OperatorTS1D(o::Operator; full=true)
    if full && !is_ts(o)
        error("o is not translation symmetric. If you want to initialize an OperatorTS1D only with its local part H_0, then set full=false")
    end
    o2 = shift_left(o)
    full && (o2 /= qubitlength(o))
    o3 = OperatorTS1D(o2.N, o2.v, o2.w, o2.coef)
end



"""
    Operator(o::OperatorTS1D)

Convert an OperatorTS1D to an Operator
"""
function Operator(o::OperatorTS1D; rs=true)
    rs && (o = resum(o))
    return Operator(o.N, o.v, o.w, o.coef)
end




"""
    is_ts(o::Operator)

return true if o is translation symmetric
"""
function is_ts(o::Operator)
    for i in 1:qubitlength(o)
        if opnorm(o - shift(o, i)) / opnorm(o) > 1e-10
            return false
        end
    end
    return true
end



"""rotate left the first n bits of x by r"""
function rotate_lower(x::Unsigned, n::Int, r::Int)
    @assert r <= n
    mask = (one(x) << n) - one(x)
    lower_bits = x & mask
    rotated_bits = (lower_bits >> r) | (lower_bits << (n - r))
    rotated_bits &= mask
    return (x & ~mask) | rotated_bits
end
function rotate_lower(p::PauliString{N,T}, r::Int) where {N,T}
    return PauliString{N,T}(rotate_lower(p.v, N, r), rotate_lower(p.w, N, r))
end



"""
    rotate(o::Operator, r::Int)

Rotate (translate/shift) left the qubits of `O` by `r`
"""
rotate(o::Operator, r::Int) = shift(o, r)

shift(o::Operator, r::Int) = Operator(shift.(o.strings, r), copy(o.coeffs))

"""shift the string v,w so it starts on site 1"""
function shift_left(p::PauliString)
    v = p.v
    w = p.w
    N = qubitlength(p)
    l = (2 * one(v))^(N + 1)
    v2 = v
    w2 = w
    for i in 0:N
        v3 = rotate_lower(v, N, i)
        w3 = rotate_lower(w, N, i)
        if v3 | w3 + v3 & w3 / 100000 < l
            v2 = v3
            w2 = w3
            l = v3 | w3 + v3 & w3 / 100000
        end
    end
    return typeof(p)(v2, w2)
end

"""
    shift_left(O::Operator)

Shift evey string left so they start on site 1. This usefull for using translation symmetry in 1D systems
# Example
```julia
A = Operator(4)
A += "XYZ1"
A += "11ZZ"
A += "1XZY"
A += "ZZ11"
```
```julia
julia> shift_left(A)
(1.0 - 0.0im) XZY1
(1.0 - 0.0im) XYZ1
(2.0 + 0.0im) ZZ11
```
"""
shift_left(O::Operator) = Operator(shift_left.(O.strings), copy(O.coeffs))

shift1(O::Operator) = shift_left(O)


function resum(o::OperatorTS1D)
    o2 = Operator(o.N)
    for i in 1:o.N
        oi = Operator(o.N, o.v, o.w, o.coef)
        o2 += shift(oi, i)
    end
    return o2
end


Base.:+(a::Number, o::OperatorTS1D) = OperatorTS1D(Operator(o, rs=false) + a / o.N; full=false)
Base.:+(o::OperatorTS1D, a::Number) = a + o


function Base.:*(o1::OperatorTS1D, o2::OperatorTS1D)
    if o1.N != o2.N
        error("Multiplying OperatorTS1D of different dimention")
    end
    N = o1.N
    d = emptydict(o1)
    for i in 1:length(o1.v)
        for j in 1:length(o2.v)
            for k in 0:N-1
                v2 = rotate_lower(o2.v[j], N, k)
                w2 = rotate_lower(o2.w[j], N, k)
                v = o1.v[i] ⊻ v2
                w = o1.w[i] ⊻ w2
                c = o1.coef[i] * o2.coef[j] * (-1)^count_ones(o1.v[i] & w2)
                if isassigned(d, (v, w))
                    d[(v, w)] += c
                else
                    insert!(d, (v, w), c)
                end
            end
        end
    end
    o = op_from_dict(d, N, OperatorTS1D)
    return OperatorTS1D(compress(shift_left(o)); full=false)
end



trace(o::OperatorTS1D) = trace(Operator(o; rs=false)) * o.N
opnorm(o::OperatorTS1D) = opnorm(Operator(o; rs=false)) * sqrt(o.N)


# diag(o::OperatorTS1D) = OperatorTS1D(diag(o.o); full=false)
# compress(o::OperatorTS1D) = OperatorTS1D(compress(o.o); full=false)


function com(o1::OperatorTS1D, o2::OperatorTS1D; epsilon::Real=0, maxlength::Int=1000, anti=false)
    s = 1
    anti && (s = -1)
    if o1.N != o2.N
        error("Multiplying OperatorTS1D of different dimention")
    end
    N = o1.N
    d = emptydict(o1)
    for i in 1:length(o1.v)
        for j in 1:length(o2.v)
            for k in 0:N-1
                v2 = rotate_lower(o2.v[j], N, k)
                w2 = rotate_lower(o2.w[j], N, k)
                v = o1.v[i] ⊻ v2
                w = o1.w[i] ⊻ w2
                k = (-1)^count_ones(o1.v[i] & w2) - s * (-1)^count_ones(o1.w[i] & v2)
                c = o1.coef[i] * o2.coef[j] * k
                if (k != 0) && (abs(c) > epsilon) && pauli_weight(v, w) < maxlength
                    if isassigned(d, (v, w))
                        d[(v, w)] += c
                    else
                        insert!(d, (v, w), c)
                    end
                end
            end
        end
    end
    o = op_from_dict(d, N, Operator)
    o = shift_left(o)
    return OperatorTS1D(o.N, o.v, o.w, o.coef)

    # return OperatorTS1D(compress(shift_left(o)); full=false)
end
