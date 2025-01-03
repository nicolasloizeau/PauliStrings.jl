module PauliStrings

export Operator
export trace, opnorm, eye, dagger, com, add, compress, ptrace, shift_left, shift, com
export diag, xcount, ycount, zcount
export truncate, trim, cutoff, prune, add_noise, k_local_part, participation
export rand_local1, rand_local2
export lanczos, rk4, norm_lanczos, rotate_lower
export op_to_strings, vw_to_string, string_to_vw, tring_to_dense, op_to_dense, get_pauli, push!, vw_in_o
export get_coefs, get_coef
export trace_product, oppow, trace_product_pow, trace_exp, moments
export OperatorTS1D, resum, rand_local1_TS1D, rand_local2_TS1D, is_ts
export all_strings, set_coefs, all_z, all_x, all_y, all_k_local
export equivalence_class
export frustration_graph

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
        N > 64 && error("N needs to be <= 64 qubits")
        new(N, Int[], Int[], Complex{Float64}[])
    end
    function Operator(N::Int, v::Vector{Int}, w::Vector{Int}, coef::Vector{Complex{Float64}})
        N > 64 && error("N needs to be <= 64 qubits")
        new(N, v, w, coef)
    end
    function Operator(pauli::String)
        N = length(pauli)
        v, w = string_to_vw(pauli)
        new(N, [v], [w], [(1im)^ycount(v, w)])
    end
end

include("operatorts1d.jl")
include("io.jl")
include("operations.jl")



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
function pauli_weight(v::Int, w::Int)
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

include("lanczos.jl")
include("truncation.jl")
include("random.jl")
include("time_evolution.jl")
include("moments.jl")
include("construction.jl")
include("equivalence.jl")
include("graph.jl")
end
