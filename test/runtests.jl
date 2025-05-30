using PauliStrings
using PauliStrings: paulistringtype
import PauliStrings as ps
using Test
using LinearAlgebra: norm
using Symbolics
using Symbolics: Num, @variables, simplify
import Symbolics


ishermitian(H::Operator) = opnorm(H - dagger(H)) < 1e-10
isunitary(U::Operator) = opnorm(U * dagger(U) - one(U)) < 1e-10
isidentity(U::Operator) = opnorm(U - one(U)) < 1e-10


include("examples.jl")

include("io.jl")
include("operator.jl")
include("algorithms.jl")
include("operatorts1d.jl")
include("construction.jl")
include("circuits.jl")
include("states.jl")
include("symbolic_paulis.jl")
include("evolution.jl")
include("truncation.jl")
