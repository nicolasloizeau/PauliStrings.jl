"""
    AbstractOperator

Abstract supertype for operators that can be represented in terms of Pauli strings.
"""
abstract type AbstractOperator end

@doc """
    paulistringtype(x::AbstractOperator)
    paulistringtype(::Type{<:AbstractOperator})
    paulistringtype(N::Integer)

Returns the type of the Pauli strings used in the operator, ie the type of a single string.
Alternatively, returns the default type for Pauli strings on `N` qubits.
""" paulistringtype
paulistringtype(x::AbstractOperator) = paulistringtype(typeof(x))
# avoid infinite recursion:
paulistringtype(T::Type{<:AbstractOperator}) = throw(MethodError(paulistringtype, T))

@doc """
    qubitlength(x::AbstractOperator)
    qubitlength(::Type{<:AbstractOperator})

Returns the number of qubits the operator acts on.
""" qubitlength
qubitlength(x::AbstractOperator) = qubitlength(typeof(x))
qubitlength(x::Type{<:AbstractOperator}) = qubitlength(paulistringtype(x))

"""
    scalartype(x::AbstractOperator)
    scalartype(::Type{<:AbstractOperator})

Returns the type of the coefficients used in the operator.
"""
scalartype
scalartype(o::AbstractOperator) = scalartype(typeof(o))
# avoid infinite recursion:
scalartype(::Type{<:AbstractOperator}) = throw(MethodError(scalartype, T))

Base.one(x::AbstractOperator) = one(typeof(x))
Base.zero(x::AbstractOperator) = zero(typeof(x))

"""
    AbstractPauliString

Abstract supertype for Pauli strings, ie strings of Pauli operators (I, X, Y, Z) acting on qubits.
"""
abstract type AbstractPauliString <: AbstractOperator end

paulistringtype(::Type{P}) where {P<:AbstractPauliString} = P

"""
    PauliString{N,T<:Unsigned} <: AbstractPauliString

A concrete type representing a Pauli string on `N` qubits as a pair of unsigned integers, as
described in https://journals.aps.org/pra/abstract/10.1103/PhysRevA.68.042318.
"""
struct PauliString{N,T<:Unsigned} <: AbstractPauliString
    v::T
    w::T
end

qubitlength(::Type{<:PauliString{N}}) where {N} = N

function uinttype(N::Integer)
    N < 0 && throw(DomainError(N, "N must be non-negative"))
    return N ≤ 8 ? UInt8 : N ≤ 16 ? UInt16 : N ≤ 32 ? UInt32 : N ≤ 64 ? UInt64 : N ≤ 128 ? UInt128 : throw(DomainError(N, "N must be <= 128"))
end
paulistringtype(N::Integer) = PauliString{N,uinttype(N)}

PauliString(pauli::AbstractString) = PauliString{length(pauli)}(pauli)
PauliString{N}(pauli::AbstractString) where {N} = PauliString{N,uinttype(N)}(pauli)

function PauliString{N,T}(pauli::AbstractString) where {N,T}
    length(pauli) == N || throw(ArgumentError("pauli string length must be $N"))
    v = zero(T)
    w = zero(T)
    two = T(2)

    for (k, p) in enumerate(pauli)
        if p == 'X'
            w += two^(k - 1)
        elseif p == 'Z'
            v += two^(k - 1)
        elseif p == 'Y'
            w += two^(k - 1)
            v += two^(k - 1)
        elseif (p != 'I') && (p != '1')
            @show p typeof(p)
            throw(ArgumentError("Invalid character in pauli string: $p"))
        end
    end

    return PauliString{N,T}(v, w)
end

Base.one(::Type{PauliString{N,T}}) where {N,T} = PauliString{N,T}(zero(T), zero(T))
Base.copy(p::PauliString{N,T}) where {N,T} = PauliString{N,T}(copy(p.v), copy(p.w))

# TODO: should we use `Char` instead of `Symbol`?
# TODO: should we use `:I` or `Symbol(1)` for identity?

@inline function Base.getindex(p::PauliString, i::Int)
    @boundscheck checkbounds(p, i)
    vi = p.v >> (i - 1) & 1
    wi = p.w >> (i - 1) & 1

    if vi & wi == 1
        return :Y
    elseif vi == 1
        return :Z
    elseif wi == 1
        return :X
    else
        return Symbol(1)
    end
end

@inline function Base.setindex(p::PauliString, val::Symbol, i::Int)
    @boundscheck checkbounds(p, i)
    val in (:X, :Y, :Z, :I, Symbol(1)) || throw(ArgumentError("Invalid value for Pauli string: `:$val`"))
    bitmask = 1 << (i - 1)
    v = if (val === :Y || val === :Z)
        p.v |= bitmask # set bit
    else
        p.v &= ~bitmask # clear bit
    end
    w = if (val === :X || val === :Y)
        p.w |= bitmask # set bit
    else
        p.w &= ~bitmask # clear bit
    end
    return typeof(p)(v, w)
end

function Base.string(x::PauliString)
    N = qubitlength(x)
    iob = IOBuffer(; sizehint=qubitlength(x))
    for i in 1:N
        print(iob, x[i])
    end
    return String(take!(iob))
end

function Base.unsigned(p::PauliString{N,T}) where {N,T}
    return (widen(p.v) << N + p.w)
end

"""
    Operator{P<:PauliString,T<:Number} <: AbstractOperator

A concrete type representing an operator as a sum of Pauli strings.
The operator is represented as a vector of Pauli strings and their corresponding coefficients.
The type parameters `P` and `T` specify the type of the Pauli strings and the type of the coefficients, respectively.
"""
struct Operator{P<:PauliString,T<:Number} <: AbstractOperator
    strings::Vector{P}
    coeffs::Vector{T}
end

"""
    Operator(N::Integer)

Initialize a zero operator on `N` qubits.
"""
Operator(N::Integer) = Operator{paulistringtype(N),ComplexF64}()
Operator{P,T}() where {P,T} = Operator{P,T}(P[], T[])

function Operator(N::Int, v::Vector{T}, w::Vector{T}, coef::Vector{Complex{Float64}}) where {T<:Unsigned}
    length(v) == length(w) == length(coef) || error("v, w, and coef must have the same length")
    P = paulistringtype(N)
    strings = P.(v, w)
    return Operator{P,ComplexF64}(strings, coef)
end

Operator(pauli::AbstractString) = Operator{paulistringtype(length(pauli))}(pauli)
Operator{P}(pauli::AbstractString) where {P} = Operator{P,ComplexF64}(pauli)
function Operator{P,T}(pauli::AbstractString) where {P,T}
    s = P(pauli)
    c = T((1.0im)^ycount(s))
    return Operator{P,T}([s], [c])
end

Operator(o::Operator) = Operator(copy(o.strings), copy(o.coeffs))

paulistringtype(::Type{<:Operator{P}}) where {P} = P
scalartype(::Type{Operator{P,T}}) where {P,T} = T

Base.one(::Type{O}) where {O<:Operator} = O([one(paulistringtype(O))], [one(scalartype(O))])
Base.zero(::Type{O}) where {O<:Operator} = O()
Base.copy(o::Operator) = typeof(o)(copy(o.strings), copy(o.coeffs))

"""
    OperatorTS1D{P<:PauliString,T<:Number} <: AbstractOperator

A concrete type representing a 1D translationally invariant operator as a sum of Pauli strings.
The operator is represented as a representative vector of Pauli strings and their corresponding coefficients,
which are implicitly repeated to form the full operator.
The type parameters `P` and `T` specify the type of the Pauli strings and the type of the coefficients, respectively.
"""
struct OperatorTS1D{P<:PauliString,T<:Number} <: AbstractOperator
    strings::Vector{P}
    coeffs::Vector{T}
end

"""
    OperatorTS1D(N::Integer)

Initialize a zero 1D translation-invariant operator on `N` qubits.
"""
OperatorTS1D(N::Integer) = OperatorTS1D{paulistringtype(N),ComplexF64}()
OperatorTS1D{P,T}() where {P,T} = OperatorTS1D{P,T}(P[], T[])

function OperatorTS1D(N::Int, v::Vector{T}, w::Vector{T}, coef::Vector{Complex{Float64}}) where {T<:Unsigned}
    length(v) == length(w) == length(coef) || error("v, w, and coef must have the same length")
    P = paulistringtype(N)
    strings = P.(v, w)
    return OperatorTS1D{P,ComplexF64}(strings, coef)
end

OperatorTS1D(pauli::AbstractString) = OperatorTS1D{paulistringtype(length(pauli))}(pauli)
OperatorTS1D{P}(pauli::AbstractString) where {P} = OperatorTS1D{P,ComplexF64}(pauli)
function OperatorTS1D{P,T}(pauli::AbstractString) where {P,T}
    s = P(pauli)
    c = T((1.0im)^ycount(s))
    return OperatorTS1D{P,T}([s], [c])
end

OperatorTS1D(o::OperatorTS1D) = OperatorTS1D(copy(o.strings), copy(o.coeffs))

paulistringtype(::Type{<:OperatorTS1D{P}}) where {P} = P
scalartype(::Type{OperatorTS1D{P,T}}) where {P,T} = T

Base.one(::Type{O}) where {O<:OperatorTS1D} = O([one(paulistringtype(O))], [one(scalartype(O)) / qubitlength(O)])
Base.zero(::Type{O}) where {O<:OperatorTS1D} = O()


"""
    OperatorTS2D{P<:PauliString,T<:Number,L1} <: AbstractOperator

A concrete type representing a 2D translationally invariant operator as a sum of Pauli strings.
The operator is represented as a representative vector of Pauli strings and their corresponding coefficients,
which are implicitly repeated to form the full operator.
The type parameters `P` and `T` specify the type of the Pauli strings and the type of the coefficients, respectively.
The parameter `L1` will keep track of the extent of the lattice in the a1 direction following the column-major representation,
that is, the flattened index is `i + j*L1`, for a lattice of size L1 * L2.
"""
struct OperatorTS2D{P<:PauliString,T<:Number,L1} <: AbstractOperator
    strings::Vector{P}
    coeffs::Vector{T}
end

"""
    OperatorTS2D(N::Integer, L1::Integer)

Initialize a zero 2D translation-invariant operator on `N` qubits, with extent of `L1` in the ``a_1`` direction.
"""
OperatorTS2D(N::Integer, L1::Integer) = (N % L1 == 0) ? OperatorTS2D{paulistringtype(N),ComplexF64,L1}() : error("N must be divisible by L1")
OperatorTS2D{P,T,L1}() where {P,T,L1} = OperatorTS2D{P,T,L1}(P[], T[])

function OperatorTS2D(N::Int, L1::Int, v::Vector{T}, w::Vector{T}, coef::Vector{Complex{Float64}}) where {T<:Unsigned}
    length(v) == length(w) == length(coef) || error("v, w, and coef must have the same length")
    P = paulistringtype(N)
    strings = P.(v, w)
    return OperatorTS2D{P,ComplexF64,L1}(strings, coef)
end

OperatorTS2D(pauli::AbstractString, L1::Integer) = OperatorTS2D{paulistringtype(length(pauli)),L1}(pauli)
OperatorTS2D{P,L1}(pauli::AbstractString) where {P,L1} = OperatorTS2D{P,ComplexF64,L1}(pauli)
function OperatorTS2D{P,T,L1}(pauli::AbstractString) where {P,T,L1}
    s = P(pauli)
    c = T((1.0im)^ycount(s))
    return OperatorTS2D{P,T,L1}([s], [c])
end

OperatorTS2D(o::OperatorTS2D) = typeof(o)(copy(o.strings), copy(o.coeffs))
extent(o::OperatorTS2D{P,T,L1}) where {P,T,L1} = L1

paulistringtype(::Type{<:OperatorTS2D{P,T,L1}}) where {P,T,L1} = P
scalartype(::Type{OperatorTS2D{P,T,L1}}) where {P,T,L1} = T

Base.one(::Type{O}) where {O<:OperatorTS2D} = O([one(paulistringtype(O))], [one(scalartype(O)) / qubitlength(O)])
Base.zero(::Type{O}) where {O<:OperatorTS2D} = O()

Base.copy(o::Union{OperatorTS1D,OperatorTS2D}) = typeof(o)(copy(o.strings), copy(o.coeffs))

"""
    Base.length(o::Operator)
    Base.length(o::OperatorTS1D)
    Base.length(o::OperatorTS2D)

Number of pauli strings in an operator

# Example
```
julia> A = Operator(4)
julia> A += "X111"
julia> A += "XYZ1"
julia> A += 2, "Y", 4
julia> length(A)
3
```
"""
Base.length(o::Union{Operator,OperatorTS1D,OperatorTS2D}) = length(o.strings)

"""
    eye(N::Int)

Identity operator on N qubits
"""
function eye(N::Int)
    O = Operator(N)
    return O + 1
end


# Utility functions
# -----------------
Base.checkbounds(::Type{Bool}, p::PauliString, i::Int) = (0 < i <= qubitlength(p))
Base.checkbounds(p::PauliString, i::Int) = checkbounds(Bool, p, i) || throw(BoundsError(p, i))

checklength(::Type{Bool}, o1::AbstractOperator) = true
checklength(::Type{Bool}, o1::AbstractOperator, o2::AbstractOperator) = qubitlength(o1) == qubitlength(o2)
checklength(::Type{Bool}, o1::AbstractOperator, o2::AbstractOperator...) = checklength(Bool, o1, first(o2)) & checklength(Bool, o1, tail(o2)...)

checklength(o1::AbstractOperator, o2::AbstractOperator...) = checklength(Bool, o1, o2...) || throw(DimensionMismatch("Operators must have the same number of qubits"))
