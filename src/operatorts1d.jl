
using Dictionaries


mutable struct OperatorTS1D
    o::Operator
    N::Int
    @doc raw"""
        OperatorTS1D(o::Operator; full=true)

    Initialize a 1D translation invariant operator from an Operator
    $O=\sum_i o_i O_i$ where $O_i=T_i(O_0)$ and $T_i$ is the translation by 1 site
    set full=true if passing $O$, an Operator that is supported on the whole chain
        full=false if passing $O_0$,a local term o such that the full operator is $O=\sum_i o_i T_i(O_0)$
    """
    function OperatorTS1D(o::Operator; full=true)
        full && !is_ts(o) && error("o needs to be translation symmetric. If you want to initialize a OperatorTS1D only with its local part H_0, then set full=false")
        o2 = shift_left(o)
        full && (o2/=o.N)
        new(o2, o.N)
    end
end

"""
    Operator(o::OperatorTS1D)

Convert an OperatorTS1D to an Operator
"""
Operator(o::OperatorTS1D) = resum(o)


"""
return true if o is translation symmetric
"""
function is_ts(o::Operator)
    for i in 1:o.N
        if opnorm(o-shift(o,i))/opnorm(o) > 1e-10
            return false
        end
    end
    return true
end



"""rotate left the first n bits of x by r"""
function rotate_lower(x::Int, n::Int, r::Int)
    @assert r<=n
    mask = (1 << n) - 1
    lower_bits = x & mask
    rotated_bits = (lower_bits >> r) | (lower_bits << (n - r))
    rotated_bits &= mask
    return (x & ~mask) | rotated_bits
end


"""rotate left the qubits of O by r"""
function shift(o::Operator, r::Int)
    o2 = deepcopy(o)
    for i in 1:length(o)
        o2.v[i]  = rotate_lower(o2.v[i], o2.N, r)
        o2.w[i]  = rotate_lower(o2.w[i], o2.N, r)
    end
    return compress(o2)
end

"""shift the string v,w so it starts on site 1"""
function shift_left(v,w,N)
    l = 2^(N+1)
    v2 = v
    w2 = w
    for i in 0:N
        v3 = rotate_lower(v, N, i)
        w3 = rotate_lower(w, N, i)
        if v3|w3+v3&w3/100000 < l
            v2 = v3
            w2 = w3
            l = v3|w3+v3&w3/100000
        end
    end
    return (v2,w2)
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
function shift_left(O::Operator)
    d = UnorderedDictionary{Tuple{Int, Int}, Complex{Float64}}()
    for i in 1:length(O)
        v,w  = shift_left(O.v[i],O.w[i],O.N)
        c = O.coef[i]
        if isassigned(d, (v,w))
            d[(v,w)] += c
        else
            insert!(d, (v,w), c)
        end
    end
    return op_from_dict(d, O.N)
end

shift1(O::Operator) = shift_left(O)





function Base.length(o::OperatorTS1D)
    return length(o.o.v)
end

Base.show(io::IO, o::OperatorTS1D) = println(io, o.o)


function resum(o::OperatorTS1D)
    o2 = Operator(o.N)
    for i in 1:o.N
        o2 += shift(o.o, i)
    end
    return o2
end


function Base.:*(o1::OperatorTS1D, o2::OperatorTS1D)
    o1 = o1.o
    o2 = o2.o
    if o1.N != o2.N
        error("Multiplying OperatorTS1D of different dimention")
    end
    N = o1.N
    d = UnorderedDictionary{Tuple{Int, Int}, Complex{Float64}}()
    for i in 1:length(o1.v)
        for j in 1:length(o2.v)
            for k in 0:N-1
                v2 = rotate_lower(o2.v[j], N, k)
                w2 = rotate_lower(o2.w[j], N, k)
                v = o1.v[i] ⊻ v2
                w = o1.w[i] ⊻ w2
                c = o1.coef[i] * o2.coef[j] * (-1)^count_ones(o1.v[i] & w2)
                if isassigned(d, (v,w))
                    d[(v,w)] += c
                else
                    insert!(d, (v,w), c)
                end
            end
        end
    end
    o = op_from_dict(d, N)
    return OperatorTS1D(compress(shift_left(o)); full=false)
end


Base.:+(o1::OperatorTS1D, o2::OperatorTS1D) = OperatorTS1D(o1.o+o2.o; full=false)
Base.:*(o::OperatorTS1D, a::Number) = OperatorTS1D(o.o*a; full=false)
Base.:/(o::OperatorTS1D, a::Number) = OperatorTS1D(o.o/a; full=false)
Base.:*(a::Number, o::OperatorTS1D) = o*a
Base.:+(a::Number, o::OperatorTS1D) = OperatorTS1D(o.o+a/o.N; full=false)
Base.:+(o::OperatorTS1D, a::Number) = a+o
Base.:-(o::OperatorTS1D) = -1*o
Base.:-(o1::OperatorTS1D, o2::OperatorTS1D) = o1+(-o2)
Base.:-(o::OperatorTS1D, a::Number) = o+(-a)
Base.:-(a::Number, o::OperatorTS1D) = a+(-o)

trace(o::OperatorTS1D) = trace(o.o)*o.N
opnorm(o::OperatorTS1D) = sqrt(o.N)*opnorm(o.o)
dagger(o::OperatorTS1D) = OperatorTS1D(dagger(o.o))


diag(o::OperatorTS1D) = OperatorTS1D(diag(o.o))
compress(o::OperatorTS1D) = OperatorTS1D(compress(o.o))


function com(o1::OperatorTS1D, o2::OperatorTS1D; anti=false)
    s = 1
    anti && (s=-1)
    o1 = o1.o
    o2 = o2.o
    if o1.N != o2.N
        error("Multiplying OperatorTS1D of different dimention")
    end
    N = o1.N
    d = UnorderedDictionary{Tuple{Int, Int}, Complex{Float64}}()
    for i in 1:length(o1.v)
        for j in 1:length(o2.v)
            for k in 0:N-1
                v2 = rotate_lower(o2.v[j], N, k)
                w2 = rotate_lower(o2.w[j], N, k)
                v = o1.v[i] ⊻ v2
                w = o1.w[i] ⊻ w2
                k = (-1)^count_ones(o1.v[i] & w2) - s*(-1)^count_ones(o1.w[i] & v2)
                c = o1.coef[i] * o2.coef[j] * k
                if isassigned(d, (v,w))
                    d[(v,w)] += c
                else
                    insert!(d, (v,w), c)
                end
            end
        end
    end
    o = op_from_dict(d, N)
    return OperatorTS1D(compress(shift_left(o)); full=false)
end


# TRUNCATION
# OperatorTS1D versions of functions from truncation.jl
PauliStrings.truncate(o::OperatorTS1D, N::Int; keepnorm::Bool = false) = OperatorTS1D(truncate(o.o, N; keepnorm=keepnorm); full=false)
k_local_part(o::OperatorTS1D, k::Int) = OperatorTS1D(k_local_part(o.o, k); full=false)
trim(o::OperatorTS1D, N::Int; keepnorm::Bool = false, keep::Operator=Operator(0)) = OperatorTS1D(trim(o.o, N; keepnorm=keepnorm, keep=keep); full=false)
cutoff(o::OperatorTS1D, epsilon::Real; keepnorm::Bool = false) = OperatorTS1D(cutoff(o.o, epsilon; keepnorm=keepnorm); full=false)
add_noise(o::OperatorTS1D, g::Real) = OperatorTS1D(add_noise(o.o, g); full=false)
