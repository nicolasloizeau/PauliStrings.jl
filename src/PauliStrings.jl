module PauliStrings

export Operator, OperatorTS1D, Operator64, Operator128, OperatorTS1D64, OperatorTS1D128
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
export renyi_entropy
export Circuits


using Random
using LinearAlgebra
using ProgressBars
using Dictionaries

rng = MersenneTwister(0)

include("operator.jl")
include("operatorts1d.jl")
include("io.jl")
include("operations.jl")
include("lanczos.jl")
include("truncation.jl")
include("random.jl")
include("time_evolution.jl")
include("moments.jl")
include("construction.jl")
include("equivalence.jl")
include("graph.jl")
include("entropy.jl")
include("circuits.jl")

end
