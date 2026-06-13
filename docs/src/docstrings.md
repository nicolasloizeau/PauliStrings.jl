
# Documentation



## Basics
```@docs
Operator(N::Int)
OperatorTS{Ls}(o)
qubitlength
```

## Truncation
```@docs
truncate(o::Operator, max_lenght::Int; keepnorm::Bool = false)
k_local_part
trim
prune
cutoff
xpart
ypart
zpart
```

## Noise
```@docs
add_noise
add_dephasing_noise
```



## Time evolution
```@docs
evolve
EvolutionResult
AbstractEvolutionMethod
RK4
DOPRI5
Trotter
TrotterTS
Exact
TrotterGate
trotterize
trotter_step!
pauli_rotation
rk4
rk4_lindblad
```


## Other algorithms
```@docs
lanczos
lioms
```



## Operations
```@docs
Base.:+(o1::O, o2::O) where {O<:AbstractOperator}
Base.:-(o1::O, o2::O) where {O<:AbstractOperator}
Base.:*(o1::Operator, o2::Operator; kwargs...)
Base.:^(o::Operator, k::Int)
commutator
anticommutator
Base.:/(o::AbstractOperator, a::Number)
compress(o::AbstractOperator)
trace(o::Operator; normalize=false)
LinearAlgebra.diag(o::AbstractOperator)
LinearAlgebra.norm(o::AbstractOperator; normalize=false)
Base.adjoint(o::AbstractOperator)
ptrace(o::AbstractOperator, keep::Vector{Int})
```



## Power and moments
```@docs
Base.:^(o::AbstractOperator, k::Int)
trace_product
trace_product_z
moments
```


## Random operators
```@docs
rand_local1
rand_local2
rand_local1_TS1D
rand_local2_TS1D
```

## Construction
```@docs
Base.:+(o::Operator, args::Tuple{Number,Vararg{Any}})
complete_basis
k_local_basis
k_local_basis_1d
x_basis
y_basis
z_basis
majorana
string_2d
all_strings
```


## Tranlation symmetry
```@docs
PauliStringTS
OperatorTS
OperatorTS1D
periodicflags
representative
resum
is_ts
is_ts2d
```

## States
```@docs
trace_zpart(o::Operator)
expect(o::Operator, state::String)
expect(o::Operator, in_state::String, out_state::String)
expect_product(o1::Operator, o2::Operator, state::String)
```

## Circuits

```@docs
Circuit
push!(c::Circuit, gate::String, sites::Real...)
pushfirst!(c::Circuit, gate::String, sites::Real...)
XGate
UGate
RXGate
PhaseGate
CPhaseGate
CXGate
SwapGate
CSXGate
MCZGate
XXPlusYYGate
grover_diffusion
compile
expect(c::Circuit, state::String)
expect(c::Circuit, in_state::String, out_state::String)
```

## I/O and conversion
```@docs
op_to_strings(o::Operator)
get_coeffs
Matrix(o::Operator)
SparseArrays.sparse(pauli::PauliString)
SparseArrays.sparse(o::Operator)
get_coeff
get_pauli
Base.string(x::PauliString)
```

## Symbolics 
```@docs
OperatorSymbolics(N::Int)
simplify_operator(o::Operator{P,Complex{Num}}) where {P}
substitute_operator(o::Operator{P,Complex{Num}}, dict::Dict) where {P}
```


## Other tools
```@docs
compress
xcount
ycount
zcount
pauli_weight
support
```

## Index

```@index
```
