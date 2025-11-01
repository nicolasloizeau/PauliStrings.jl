
"""
    AbstractOperator

Abstract supertype for operators that can be represented in terms of Pauli strings.
"""
abstract type AbstractOperator end


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
PauliString{N}(v::Integer, w::Integer) where {N} = PauliString{N,uinttype(N)}(v, w)

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

Base.string(x::PauliString) = join([x[i] for i = 1:qubitlength(x)])

Base.unsigned(p::PauliString{N,T}) where {N,T} = (widen(p.v) << N + p.w)
