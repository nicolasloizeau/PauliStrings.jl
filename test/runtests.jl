using PauliStrings
import PauliStrings as ps
using Test
using LinearAlgebra: norm



ishermitian(H::Operator) = opnorm(H-dagger(H)) < 1e-10
isunitary(U::Operator) = opnorm(U * dagger(U) - eye(U.N)) < 1e-10
isidentity(U::Operator) = opnorm(U - eye(U.N)) < 1e-10


include("examples.jl")

include("operator.jl")
include("algorithms.jl")
include("operatorts1d.jl")
include("circuits.jl")
