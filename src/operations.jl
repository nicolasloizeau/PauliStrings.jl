


function uinttype(o::Operator)
    if typeof(o) == Operator64 || typeof(o) == OperatorTS1D64
        return UInt64
    elseif typeof(o) == Operator128 || typeof(o) == OperatorTS1D128
        return UInt128
    else
        error("Type not recognized")
    end
end

function emptydict(o::Operator)
    T = uinttype(o)
    return UnorderedDictionary{Tuple{T,T},Complex{Float64}}()
end



"""
    add(o1::Operator, o2::Operator)
    Base.:+(o1::Operator, o2::Operator)
    Base.:+(o::Operator, a::Number)
    Base.:+(a::Number, o::Operator)

Add two operators together or add a number to an operator

# Example
```julia
A = Operator(4)
A += "XYZ1"
A += 1, "Y", 4
B = Operator(4)
B += 2, "Y", 2, "Y", 4
B += 1, "Z", 3
```

```
julia> A
(1.0 - 0.0im) 111Y
(1.0 - 0.0im) XYZ1

julia> B
(1.0 + 0.0im) 11Z1
(2.0 - 0.0im) 1Y1Y

julia> A+B
(1.0 + 0.0im) 11Z1
(2.0 - 0.0im) 1Y1Y
(1.0 - 0.0im) 111Y
(1.0 - 0.0im) XYZ1

julia> A+5
(1.0 - 0.0im) 111Y
(1.0 - 0.0im) XYZ1
(5.0 + 0.0im) 1111
```
"""
function Base.:+(o1::Operator, o2::Operator)
    @assert o1.N == o2.N "Adding operators of different dimention"
    @assert typeof(o1) == typeof(o2) "Adding operators of different types"
    o3 = typeof(o1)(o1.N)
    o3.v = vcat(o1.v, o2.v)
    o3.w = vcat(o1.w, o2.w)
    o3.coef = vcat(o1.coef, o2.coef)
    return compress(o3)
end


function Base.:+(o::Operator, a::Number)
    o1 = deepcopy(o)
    i = ione(o)
    if i >= 0
        o1.coef[ione(o)] += a
    else
        push!(o1.coef, a)
        push!(o1.v, 0)
        push!(o1.w, 0)
    end
    return o1
end

Base.:+(a::Number, o::Operator) = o + a

"""
    Base.:*(o1::Operator, o2::Operator)
    Base.:*(o::Operator, a::Number)
    Base.:*(a::Number, o::Operator)

Multiply two operators together or an operator with a number

# Example
```julia
A = Operator(4)
A += "XYZ1"
A += 1, "Y", 4
B = Operator(4)
B += 2, "Y", 2, "Y", 4
B += 1, "Z", 3
```

```
julia> A
(1.0 - 0.0im) 111Y
(1.0 - 0.0im) XYZ1


julia> B
(1.0 + 0.0im) 11Z1
(2.0 - 0.0im) 1Y1Y

julia> A*B
(2.0 - 0.0im) X1ZY
(1.0 - 0.0im) 11ZY
(2.0 - 0.0im) 1Y11
(1.0 - 0.0im) XY11

julia> A*5
(5.0 - 0.0im) 111Y
(5.0 - 0.0im) XYZ1
```
"""
function Base.:*(o1::Operator, o2::Operator)
    @assert o1.N == o2.N "Multiplying operators of different dimention"
    @assert typeof(o1) == typeof(o2) "Multiplying operators of different types"
    d = emptydict(o1)
    for i in 1:length(o1.v)
        for j in 1:length(o2.v)
            v = o1.v[i] ⊻ o2.v[j]
            w = o1.w[i] ⊻ o2.w[j]
            c = o1.coef[i] * o2.coef[j] * (-1)^count_ones(o1.v[i] & o2.w[j])
            if isassigned(d, (v, w))
                d[(v, w)] += c
            else
                insert!(d, (v, w), c)
            end
        end
    end
    return op_from_dict(d, o1.N)
end


function op_from_dict(d::UnorderedDictionary{Tuple{T,T},Complex{Float64}}, N::Int) where {T<:Unsigned}
    o = Operator(N)
    for (v, w) in keys(d)
        push!(o.v, v)
        push!(o.w, w)
        push!(o.coef, d[(v, w)])
    end
    return o
end



function Base.:*(o::Operator, a::Number)
    o1 = deepcopy(o)
    o1.coef .*= a
    return o1
end

Base.:*(a::Number, o::Operator) = o * a


"""
    Base.:/(o::Operator, a::Number)

Divide an operator by a number
"""
function Base.:/(o::Operator, a::Number)
    o1 = deepcopy(o)
    o1.coef ./= a
    return o1
end

"""
    Base.:-(o::Operator)
    Base.:-(o1::Operator, o2::Operator)
    Base.:-(o::Operator, a::Real)
    Base.:-(a::Real, o::Operator)

Subtraction between operators and numbers
"""
Base.:-(o::Operator) = -1 * o
Base.:-(o1::Operator, o2::Operator) = o1 + (-o2)
Base.:-(o::Operator, a::Number) = o + (-a)
Base.:-(a::Number, o::Operator) = a + (-o)


"""
    com(o1::Operator, o2::Operator; epsilon::Real=0, maxlength::Int=1000)
    com(o1::OperatorTS1D, o2::OperatorTS1D; anti=false)

Commutator of two operators. Set anti=true to compute the anti-commutator.

# Example
```
julia> A = Operator(4)
julia> A += "X111"
julia> B = Operator(4)
julia> B += "Z111"
julia> B += "XYZ1"
julia> com(A,B)
(0.0 - 2.0im) Y111
```
"""
function com(o1::Operator, o2::Operator; epsilon::Real=0, maxlength::Int=1000, anti=false)
    s = 1
    anti && (s = -1)
    @assert o1.N == o2.N "Commuting operators of different dimention"
    @assert typeof(o1) == typeof(o2) "Commuting operators of different types"
    o3 = Operator(o1.N)
    d = emptydict(o1)
    for i in 1:length(o1.v)
        for j in 1:length(o2.v)
            v = o1.v[i] ⊻ o2.v[j]
            w = o1.w[i] ⊻ o2.w[j]
            k = (-1)^count_ones(o1.v[i] & o2.w[j]) - s * (-1)^count_ones(o1.w[i] & o2.v[j])
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
    for (v, w) in keys(d)
        push!(o3.v, v)
        push!(o3.w, w)
        push!(o3.coef, d[(v, w)])
    end
    return o3
end



"""
    com(v1::Unsigned, w1::Unsigned, v2::Unsigned, w2::Unsigned)

Commutator of two pauli strings in integer representation
Return k,v,w
"""
function com(v1::Unsigned, w1::Unsigned, v2::Unsigned, w2::Unsigned)
    v = v1 ⊻ v2
    w = w1 ⊻ w2
    k = (-1)^count_ones(v1 & w2) - (-1)^count_ones(w1 & v2)
    return k, v, w
end


"""
    compress(o::Operator)

Accumulate repeated terms and remove terms with a coeficient smaller than 1e-16
"""
function compress(o::Operator)
    o1 = typeof(o)(o.N)
    T = uinttype(o)
    vw = Set{Tuple{T,T}}(zip(o.v, o.w))
    d = UnorderedDictionary{Tuple{T,T},Complex{Float64}}(vw, zeros(length(vw)))
    for i in 1:length(o)
        v = o.v[i]
        w = o.w[i]
        d[(v, w)] += o.coef[i]
    end
    for (v, w) in keys(d)
        if abs(d[(v, w)]) > 1e-16
            push!(o1.v, v)
            push!(o1.w, w)
            push!(o1.coef, d[(v, w)])
        end
    end
    return o1
end


"""return the index of the 1 string"""
function ione(o::Operator)
    for i in 1:length(o)
        if o.v[i] == 0 && o.w[i] == 0
            return i
        end
    end
    return -1
end

"""
    trace(o::Operator)
    trace(o::OperatorTS1D)

Trace of an operator

# Example
```
julia> A = Operator(4)
julia> A += 2,"1111"
julia> A += 3,"XYZ1"
julia> trace(A)
32.0 + 0.0im
```
"""
function trace(o::Operator)
    t = 0
    for i in 1:length(o.v)
        if o.v[i] == 0 && o.w[i] == 0
            t += o.coef[i]
        end
    end
    return t * 2^o.N
end


"""
    diag(o::Operator)
    diag(o::OperatorTS1D)

Diagonal of an operator. Keep the strings that only contain 1's or Z's.
Return another operator.

# Example
```
julia> A = Operator(4)
julia> A += 2,"1111"
julia> A += 3,"XYZ1"
julia> A += 3,"Z11Z"
julia> diag(A)
(2.0 + 0.0im) 1111
(3.0 + 0.0im) Z11Z
```
"""
function diag(o::Operator)
    o2 = Operator(o.N)
    for i in 1:length(o)
        v = o.v[i]
        w = o.w[i]
        if xcount(v, w) == 0 && ycount(v, w) == 0
            push!(o2.coef, o.coef[i])
            push!(o2.v, v)
            push!(o2.w, w)
        end
    end
    return o2
end

"""
    ycount(v::Unsigned, w::Unsigned)

Count the number of Y in a string
"""
ycount(v::Unsigned, w::Unsigned) = count_ones(v & w)

"""
    zcount(v::Unsigned, w::Unsigned)

Count the number of Z in a string
"""
zcount(v::Unsigned, w::Unsigned) = count_ones(v & ~w)

"""
    xcount(v::Unsigned, w::Unsigned)

Count the number of X in a string
"""
xcount(v::Unsigned, w::Unsigned) = count_ones(~v & w)



"""
    opnorm(o::Operator)
    opnorm(o::OperatorTS1D)

Frobenius norm

# Example
```
julia> A = Operator(4)
julia> A += 2,"X",2
julia> A += 1,"Z",1,"Z",3
julia> opnorm(A)
8.94427190999916
```
"""
function opnorm(o::Operator)
    return norm(o.coef) * sqrt(2^o.N)
end


"""
    dagger(o::Operator)
    dagger(o::OperatorTS1D)

Conjugate transpose

# Example
```julia
A = Operator(3)
A += 1im,"X",2
A += 1,"Z",1,"Z",3
```
```
julia> A

(1.0 + 0.0im) Z1Z
(0.0 + 1.0im) 1X1


julia> dagger(A)
(1.0 - 0.0im) Z1Z
(0.0 - 1.0im) 1X1
```
"""
function dagger(o::Operator)
    o1 = deepcopy(o)
    for i in 1:length(o1)
        s = (-1)^count_ones(o1.v[i] & o1.w[i])
        o1.coef[i] = s * conj(o1.coef[i])
    end
    return o1
end




"""
v,w encode a string.
return true if at least one index of keep is non unit in vw
"""
function tokeep(v::Unsigned, w::Unsigned, keep::Vector{Int})
    for i in keep
        if bit(v | w, i)
            return true
        end
    end
    return false
end

"""
    ptrace(o::Operator, keep::Vector{Int})

Partial trace.

`keep` is list of qubits indices to keep starting at 1
note that this still returns an operator of size N and doesnt permute the qubits
this only gets rid of Pauli strings that have no support on keep
and add their coeficient*2^N to the identity string

# Example
```julia
A = Operator(5)
A += "XY1XZ"
A += "XY11Z"
```
```
julia> ptrace(A, [3,4])
(1.0 - 0.0im) XY1XZ
(8.0 - 0.0im) 11111

julia> ptrace(A, [1,5])
(1.0 - 0.0im) XY1XZ
(1.0 - 0.0im) XY11Z
```
"""
function ptrace(o::Operator, keep::Vector{Int})
    o2 = Operator(o.N)
    NA = length(keep)
    NB = o.N - NA
    for i in 1:length(o)
        if tokeep(o.v[i], o.w[i], keep)
            push!(o2.v, o.v[i])
            push!(o2.w, o.w[i])
            push!(o2.coef, o.coef[i])
        else
            o2 += o.coef[i] * 2^NB / (1im)^count_ones(o.v[i] & o.w[i])
        end
    end
    return o2
end


"""
    vw_in_o(v::Unsigned, w::Unsigned, o::Operator)

Return true is string (v,w) is in o
"""
function vw_in_o(v::Unsigned, w::Unsigned, o::Operator)
    for i in 1:length(o)
        if v == o.v[i] && w == o.w[i]
            return true
        end
    end
    return false
end
