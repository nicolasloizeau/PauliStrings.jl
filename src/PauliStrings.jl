module PauliStrings

export Operator
export trace, opnorm, eye, dagger, com, add, compress, ptrace, shift_left
export diag, xcount, ycount, zcount
export truncate, trim, cutoff, prune, add_noise, k_local_part
export rand_local1, rand_local2
export lanczos, rk4, norm_lanczos, rotate_lower
export op_to_strings, vw_to_string

using Random
using LinearAlgebra
using ProgressBars
using Dictionaries

rng = MersenneTwister(0)

"""
operator as a sum of pauli string encoded like in
https://journals.aps.org/pra/abstract/10.1103/PhysRevA.68.042318
intialized as : O=Operator(N)
where N is the number of qubits
"""
mutable struct Operator
    N::Int
    v::Vector{Int}
    w::Vector{Int}
    coef::Vector{Complex{Float64}}
    @doc """
        Operator(N::Int)

    Initialize an empty operator on N qubits
    """
    function Operator(N::Int)
        new(N, Int[], Int[], Complex{Float64}[])
    end
    function Operator(N::Int, v::Vector{Int}, w::Vector{Int}, coef::Vector{Complex{Float64}})
        new(N, v, w, coef)
    end
end

include("io.jl")
include("operations.jl")



"""
    Base.length(o::Operator)

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
    return O+1
end

"""number of non unit paulis in a string encoded by v,w"""
function pauli_weight(v::Int,w::Int)
    return count_ones(v | w)
end


"""returns the position of (v,w) in O. return 0 if (v,w) not in O"""
function posvw(v,w,O)
    for i in 1:length(O)
        if O.v[i] == v && O.w[i] == w
            return i
        end
    end
    return 0
end


include("lanczos.jl")
include("truncation.jl")
include("random.jl")
include("time_evolution.jl")
end
