
abstract type Operator end
abstract type OperatorTS1D <: Operator end


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
