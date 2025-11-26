
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
    Operator{P<:PauliString,T<:Number} <: AbstractOperator

A concrete type representing an operator as a sum of Pauli strings.
The operator is represented as a vector of Pauli strings and their corresponding coefficients.
The type parameters `P` and `T` specify the type of the Pauli strings and the type of the coefficients, respectively.
"""
struct Operator{P<:AbstractPauliString,T<:Number} <: AbstractOperator
    strings::Vector{P}
    coeffs::Vector{T}
end

"""
    Operator(N::Integer)

Initialize a zero operator on `N` qubits.
"""
Operator(N::Integer) = Operator{paulistringtype(N),ComplexF64}()
Operator{P,T}() where {P,T} = Operator{P,T}(P[], T[])
# Operator(strings::Vector{<:AbstractPauliString}, coeffs::Vector{<:Number}) = Operator{eltype(strings),eltype(coeffs)}(strings, coeffs)

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


Operator(pauli::PauliString) = Operator{typeof(pauli),ComplexF64}([pauli], [(1.0im)^ycount(pauli)])

paulistringtype(::Type{<:Operator{P}}) where {P} = P
scalartype(::Type{Operator{P,T}}) where {P,T} = T

Base.one(::Type{O}) where {O<:Operator} = O([one(paulistringtype(O))], [one(scalartype(O))])
Base.zero(::Type{O}) where {O<:Operator} = O()
Base.copy(o::Operator) = typeof(o)(copy(o.strings), copy(o.coeffs))

"""
    Base.length(o::Operator)

Number of Pauli strings in an operator

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
Base.length(o::Union{Operator}) = length(o.strings)


"""
    eye(N::Int)

Identity operator on N qubits
"""
eye(N::Int) = Operator(N) + 1


# Utility functions
# -----------------
Base.checkbounds(::Type{Bool}, p::PauliString, i::Int) = (0 < i <= qubitlength(p))
Base.checkbounds(p::PauliString, i::Int) = checkbounds(Bool, p, i) || throw(BoundsError(p, i))

checklength(::Type{Bool}, o1::AbstractOperator) = true
checklength(::Type{Bool}, o1::AbstractOperator, o2::AbstractOperator) = qubitlength(o1) == qubitlength(o2)
checklength(::Type{Bool}, o1::AbstractOperator, o2::AbstractOperator...) = checklength(Bool, o1, first(o2)) & checklength(Bool, o1, tail(o2)...)

checklength(o1::AbstractOperator, o2::AbstractOperator...) = checklength(Bool, o1, o2...) || throw(DimensionMismatch("Operators must have the same number of qubits"))
checklength(::Type{Bool}, o1::AbstractOperator, o2::PauliString) = qubitlength(o1) == qubitlength(o2)
checklength(::Type{Bool}, o1::PauliString, o2::AbstractOperator) = qubitlength(o1) == qubitlength(o2)
