using PauliStrings
import PauliStrings as ps
using Test
using LinearAlgebra: norm



ishermitian(H::Operator) = opnorm(H - dagger(H)) < 1e-10
isunitary(U::Operator) = opnorm(U * dagger(U) - one(U)) < 1e-10
isidentity(U::Operator) = opnorm(U - one(U)) < 1e-10


include("examples.jl")

include("basic_tests.jl")
include("io.jl")
include("operator.jl")
include("operations_strings.jl")
include("algorithms.jl")
include("operatorts1d.jl")
include("operatorts2d.jl")
include("construction.jl")
include("circuits.jl")
include("states.jl")
include("evolution.jl")
include("truncation.jl")
include("symbolics.jl")
