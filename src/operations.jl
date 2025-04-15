# PauliString operations
# ======================

@inline Base.:(==)(p1::P, p2::P) where {P<:PauliString} = (p1.v == p2.v) & (p1.w == p2.w)

# magic number for the Fibonacci hash function in UInt
const fib_magic_64 = 0x9e3779b97f4a7c15
Base.hash(p::PauliString{N,UInt64}, h::UInt) where {N} = hash(p.v * fib_magic_64 + p.w, h)

# unary operations
# ----------------
"""
    xcount(p::PauliString)

Count the number of X operators in a string.
"""
xcount(p::PauliString) = count_ones(~p.v & p.w)

"""
    ycount(p::PauliString)

Count the number of Y operators in a string.
"""
ycount(p::PauliString) = count_ones(p.v & p.w)

"""
    zcount(p::PauliString)

Count the number of Z operators in a string.
"""
zcount(p::PauliString) = count_ones(p.v & ~p.w)

"""
    pauli_weight(p::PauliString)

Count the number of non unit operators in a string.
"""
pauli_weight(p::PauliString) = count_ones(p.v | p.w)

# TODO: do we want to name this Base.circshift?
"""
    shift(p::PauliString, i::Int)

Rotate the Pauli string `p` by `i` qubits to the left.
"""
function shift(p::PauliString, i::Int)
    N = qubitlength(p)
    return typeof(p)(rotate_lower(p.v, N, i), rotate_lower(p.w, N, i))
end


# binary operations
# -----------------
Base.xor(p1::P, p2::P) where {P<:PauliString} = P(p1.v ⊻ p2.v, p1.w ⊻ p2.w)

function commutator(p1::P, p2::P) where {P<:PauliString}
    p = p1 ⊻ p2
    k = ((count_ones(p2.v & p1.w) & 1) << 1) - ((count_ones(p1.v & p2.w) & 1) << 1)
    return p, k
end

function anticommutator(p1::P, p2::P) where {P<:PauliString}
    p = p1 ⊻ p2
    k = 2 - (((count_ones(p1.v & p2.w) & 1) << 1) + ((count_ones(p1.w & p2.v) & 1) << 1))
    return p, k
end

function prod(p1::P, p2::P) where {P<:PauliString}
    p = p1 ⊻ p2
    k = 1 - ((count_ones(p1.v & p2.w) & 1) << 1)
    return p, k
end



emptydict(o::AbstractOperator) = UnorderedDictionary{eltype(o.strings),eltype(o.coeffs)}()



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
function Base.:+(o1::O, o2::O) where {O<:AbstractOperator}
    checklength(o1, o2)

    d = emptydict(o1)

    # add the first operator
    ps, cs = o1.strings, o1.coeffs
    length(ps) == length(cs) || throw(DimensionMismatch("strings and coefficients must have the same length"))
    @inbounds for i in eachindex(ps)
        setwith!(+, d, ps[i], cs[i])
    end
    # add the second operator
    ps, cs = o2.strings, o2.coeffs
    length(ps) == length(cs) || throw(DimensionMismatch("strings and coefficients must have the same length"))
    @inbounds for i in eachindex(ps)
        setwith!(+, d, ps[i], cs[i])
    end

    # assemble output
    return typeof(o1)(collect(keys(d)), collect(values(d)))
end

function Base.:-(o1::O, o2::O) where {O<:AbstractOperator}
    checklength(o1, o2)

    d = emptydict(o1)

    # add the first operator
    ps, cs = o1.strings, o1.coeffs
    length(ps) == length(cs) || throw(DimensionMismatch("strings and coefficients must have the same length"))
    @inbounds for i in eachindex(ps)
        setwith!(+, d, ps[i], cs[i])
    end
    # subtract the second operator
    ps, cs = o2.strings, o2.coeffs
    length(ps) == length(cs) || throw(DimensionMismatch("strings and coefficients must have the same length"))
    @inbounds for i in eachindex(ps)
        setwith!(-, d, ps[i], cs[i])
    end

    # assemble output
    return typeof(o1)(collect(keys(d)), collect(values(d)))
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
    binary_kernel(f, A::Operator, B::Operator; epsilon::Real=0, maxlength::Int=1000)

Compute-kernel of applying a function `f` to all pairs of strings in two operators `A` and `B`,
reducing the result to a new operator.
"""
function binary_kernel(f, A::Operator, B::Operator; epsilon::Real=0, maxlength::Int=1000)
    checklength(A, B)

    d = emptydict(A) # reducer
    p1s, c1s = A.strings, A.coeffs
    p2s, c2s = B.strings, B.coeffs

    # check boundaries to safely use `@inbounds`
    length(p1s) == length(c1s) || throw(DimensionMismatch("strings and coefficients must have the same length"))
    length(p2s) == length(c2s) || throw(DimensionMismatch("strings and coefficients must have the same length"))

    # core kernel logic
    @inbounds for i1 in eachindex(p1s)
        p1, c1 = p1s[i1], c1s[i1]
        for i2 in eachindex(p2s)
            p2, c2 = p2s[i2], c2s[i2]
            p, k = f(p1, p2)
            c = c1 * c2 * k
            if (k != 0) && abs(c) > epsilon && pauli_weight(p) < maxlength
                setwith!(+, d, p, c)
            end
        end
    end

    # assemble output
    return Operator{keytype(d),valtype(d)}(collect(keys(d)), collect(values(d)))
end

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
    return binary_kernel(prod, o1, o2)
end

function commutator(o1::Operator, o2::Operator)
    return binary_kernel(commutator, o1, o2)
end

function anticommutator(o1::Operator, o2::Operator)
    return binary_kernel(anticommutator, o1, o2)
end



Base.:*(o::Operator, a::Number) = Operator(copy(o.strings), o.coeffs * a)
Base.:*(a::Number, o::AbstractOperator) = o * a

"""
    Base.:/(o::AbstractOperator, a::Number)

Divide an operator by a number
"""
Base.:/(o::AbstractOperator, a::Number) = o * inv(a)
Base.:\(a::Number, o::AbstractOperator) = o * inv(a)

"""
    Base.:-(o::Operator)
    Base.:-(o1::Operator, o2::Operator)
    Base.:-(o::Operator, a::Real)
    Base.:-(a::Real, o::Operator)

Subtraction between operators and numbers
"""
Base.:-(o::AbstractOperator) = -1 * o
Base.:-(o1::AbstractOperator, o2::AbstractOperator) = o1 + (-o2)
Base.:-(o::AbstractOperator, a::Number) = o + (-a)
Base.:-(a::Number, o::AbstractOperator) = a + (-o)


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
    compress(o::AbstractOperator)

Accumulate repeated terms and remove terms with a coeficient smaller than 1e-16
"""
function compress(o::AbstractOperator)
    d = emptydict(o)
    ps, cs = o.strings, o.coeffs
    length(ps) == length(cs) || throw(DimensionMismatch("strings and coefficients must have the same length"))
    @inbounds for i in eachindex(ps)
        setwith!(+, d, ps[i], cs[i])
    end
    return typeof(o)(collect(keys(d)), collect(values(d)))
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
    trace(o::Operator; normalize=false)
    trace(o::OperatorTS1D)

Trace of an operator. If normalize is true, return the trace divided by `2^N`.

# Example
```
julia> A = Operator(4)
julia> A += 2,"1111"
julia> A += 3,"XYZ1"
julia> trace(A)
32.0 + 0.0im
```
"""
function trace(o::Operator; normalize=false)
    t = 0
    for i in 1:length(o.v)
        if o.v[i] == 0 && o.w[i] == 0
            t += o.coef[i]
        end
    end
    if normalize
        return t
    end
    return t * 2.0^o.N
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
    return norm(o.coeffs) * (2.0^(qubitlength(o) / 2))
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
