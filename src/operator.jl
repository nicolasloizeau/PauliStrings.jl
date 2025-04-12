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

"""
    PauliString{N,T<:Unsigned} <: AbstractPauliString

A concrete type representing a Pauli string on `N` qubits as a pair of unsigned integers, as
described in https://journals.aps.org/pra/abstract/10.1103/PhysRevA.68.042318.
"""
struct PauliString{N,T<:Unsigned} <: AbstractPauliString
    v::T
    w::T
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


function uinttype(N::Int)
    if N ≤ 64
        return UInt64
    elseif N ≤ 128
        return UInt128
    else
        error("N needs to be <= 128 qubits")
    end
end


"""
operator as a sum of pauli string encoded like in
https://journals.aps.org/pra/abstract/10.1103/PhysRevA.68.042318
intialized as : O=Operator(N)
where N is the number of qubits
"""



function uinttype(o::Operator)
    if typeof(o) == Operator64 || typeof(o) == OperatorTS1D64
        return UInt64
    elseif typeof(o) == Operator128 || typeof(o) == OperatorTS1D128
        return UInt128
    else
        error("Type not recognized")
    end
end

mutable struct Operator64 <: Operator
    N::Int
    v::Vector{UInt64}
    w::Vector{UInt64}
    coef::Vector{Complex{Float64}}
end

mutable struct Operator128 <: Operator
    N::Int
    v::Vector{UInt128}
    w::Vector{UInt128}
    coef::Vector{Complex{Float64}}
end

mutable struct OperatorTS1D64 <: OperatorTS1D
    N::Int
    v::Vector{UInt64}
    w::Vector{UInt64}
    coef::Vector{Complex{Float64}}
end

mutable struct OperatorTS1D128 <: OperatorTS1D
    N::Int
    v::Vector{UInt128}
    w::Vector{UInt128}
    coef::Vector{Complex{Float64}}
end


Operator64(N::Int) = Operator64(N, UInt64[], UInt64[], Complex{Float64}[])
Operator128(N::Int) = Operator128(N, UInt128[], UInt128[], Complex{Float64}[])
OperatorTS1D64(N::Int) = OperatorTS1D64(N, UInt64[], UInt64[], Complex{Float64}[])
OperatorTS1D128(N::Int) = OperatorTS1D128(N, UInt128[], UInt128[], Complex{Float64}[])


function Operator(N::Int, v::Vector{T}, w::Vector{T}, coef::Vector{Complex{Float64}}) where {T<:Unsigned}
    if N <= 64
        return Operator64(N, v, w, coef)
    elseif N <= 128
        return Operator128(N, v, w, coef)
    else
        error("N needs to be <= 128 qubits")
    end
end

"""
    Operator(N::Int)

Initialize an empty operator on N qubits
"""
function Operator(N::Int)
    if N <= 64
        return Operator64(N)
    elseif N <= 128
        return Operator128(N)
    else
        error("N needs to be <= 128 qubits")
    end
end

function OperatorTS1D(N::Int, v::Vector{T}, w::Vector{T}, coef::Vector{Complex{Float64}}) where {T<:Unsigned}
    if N <= 64
        return OperatorTS1D64(N, v, w, coef)
    elseif N <= 128
        return OperatorTS1D128(N, v, w, coef)
    else
        error("N needs to be <= 128 qubits")
    end
end


function OperatorTS1D(N::Int)
    if N <= 64
        return OperatorTS1D64(N)
    elseif N <= 128
        return OperatorTS1D128(N)
    else
        error("N needs to be <= 128 qubits")
    end
end


function Operator(pauli::String)
    N = length(pauli)
    v, w = string_to_vw(pauli)
    return Operator(N, [v], [w], [(1.0im)^ycount(v, w)])
end


Operator(o::Operator) = deepcopy(o)

"""
    Base.length(o::Operator)
    Base.length(o::OperatorTS1D)

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
function Base.length(o::Operator)
    return length(o.v)
end

"""
    eye(N::Int)

Identity operator on N qubits
"""
function eye(N::Int)
    O = Operator(N)
    return O + 1
end

"""number of non unit paulis in a string encoded by v,w"""
function pauli_weight(v::Unsigned, w::Unsigned)
    return count_ones(v | w)
end


"""returns the position of (v,w) in O. return 0 if (v,w) not in O"""
function posvw(v, w, O)
    for i in 1:length(O)
        if O.v[i] == v && O.w[i] == w
            return i
        end
    end
    return 0
end
