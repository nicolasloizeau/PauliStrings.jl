
abstract type Operator end
abstract type OperatorTS1D <: Operator end

"""
operator as a sum of pauli string encoded like in
https://journals.aps.org/pra/abstract/10.1103/PhysRevA.68.042318
intialized as : O=Operator(N)
where N is the number of qubits
"""

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
