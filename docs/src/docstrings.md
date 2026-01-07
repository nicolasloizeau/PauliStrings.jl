
# Documentation



## Basics
```@docs
Operator(N::Int)
OperatorTS{Ls}(o)
qubitlength
```

## Truncation
```@docs
truncate(o::AbstractOperator, max_lenght::Int; keepnorm::Bool=false)
k_local_part(o::AbstractOperator, k::Int; atmost=false)
trim(o::AbstractOperator, max_strings::Int; keepnorm::Bool=false, keep::Operator=Operator(0))
prune(o::AbstractOperator, alpha::Real; keepnorm::Bool=false)
cutoff(o::AbstractOperator, epsilon::Real; keepnorm::Bool=false)
xpart(o::AbstractOperator)
ypart(o::AbstractOperator)
zpart(o::AbstractOperator)
```

## Noise
```@docs
add_noise(o::AbstractOperator, g::Real)
add_noise(o::AbstractOperator, g::AbstractVector{<:Real})
add_dephasing_noise(o::AbstractOperator, g::Real; basis::Symbol=:Z)
add_dephasing_noise(o::AbstractOperator, g::AbstractVector{<:Real})
```



## Algorithms
```@docs
rk4(H::AbstractOperator, O::AbstractOperator, dt::Real; hbar::Real=1, heisenberg=true, M=2^20, keep::Operator=Operator(0))
rk4(H::Function, O::AbstractOperator, dt::Real, t::Real; hbar::Real=1, heisenberg=true, M=2^20, keep::Operator=Operator(0))
rk4_lindblad(H::AbstractOperator, O::AbstractOperator, dt::Real, L; hbar::Real=1, heisenberg=true, M=2^20, keep::Operator=Operator(0), gamma=[])
lanczos(H::AbstractOperator, O::AbstractOperator, steps::Int, nterms::Int; keepnorm=true, maxlength=1000, returnOn=false, observer=false)
lioms(H::AbstractOperator, support::Vector{T}; threshold::Real=1e-14, f::Function=(H,O)->im * commutator(H,O)) where {T<:AbstractOperator}
```

## Operations
```@docs
Base.:+(o1::O, o2::O) where {O<:AbstractOperator}
Base.:-(o1::O, o2::O) where {O<:AbstractOperator}
Base.:*(o1::Operator, o2::Operator; kwargs...)
Base.:^(o::Operator, k::Int)
commutator(o1::Operator, o2::Operator; kwargs...)
anticommutator(o1::Operator, o2::Operator; kwargs...)
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
trace_product(o1::Operator, o2::Operator; scale=0)
trace_product(A::AbstractOperator, k::Int, B::AbstractOperator, l::Int; scale=0)
trace_product(A::AbstractOperator, k::Int; scale=0)
trace_product_z(o1::AbstractOperator, o2::AbstractOperator; scale=0)
moments(H::AbstractOperator, kmax::Int; start=1, scale=0)
```


## Random operators
```@docs
rand_local1(N::Int)
rand_local2(N::Int)
rand_local1_TS1D(N::Int)
rand_local2_TS1D(N::Int)
```

## Construction
```@docs
Base.:+(o::Operator, args::Tuple{Number,Vararg{Any}})
complete_basis(N::Int)
k_local_basis(N::Int, k::Int; atmost=false)
k_local_basis_1d(N::Int, k::Int; translational_symmetry::Bool=false)
x_basis(N::Int)
y_basis(N::Int)
z_basis(N::Int)
majorana(N::Int, k::Int)
string_2d(args::Tuple{Vararg{Any}}, L1::Int, L2::Int; pbc=false)
```


## Tranlation symmetry
```@docs
PauliStringTS{Ls}(p::PauliString)
OperatorTS{Ls}(o::Operator)
OperatorTS1D(o::Operator; full=true)
representative(o::OperatorTS)
representative(p::PauliStringTS)
resum(o::OperatorTS)
is_ts(o::Operator)
is_ts(o::Operator, Ls::Tuple)
is_ts2d(o::Operator, L1)
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
Circuit(N::Int; max_strings=2^30, noise_amplitude=0)
push!(c::Circuit, gate::String, sites::Real...)
pushfirst!(c::Circuit, gate::String, sites::Real...)
XGate(N::Int, i::Int)
UGate(N::Int, i::Int, theta::Real, phi::Real, lam::Real)
RXGate(N::Int, i::Int, theta::Real)
PhaseGate(N::Int, i::Int, theta::Real)
CPhaseGate(N::Int, i::Int, j::Int, theta::Real)
CXGate(N::Int, i::Int, j::Int)
SwapGate(N::Int, i::Int, j::Int)
CSXGate(N::Int, i::Int, j::Int)
CCXGate(N::Int, i::Int, j::Int, k::Int)
MCZGate(N::Int, sites::Int...)
XXPlusYYGate(N::Int, i::Int, j::Int, theta::Real, beta::Real)
grover_diffusion(N::Int, sites::Int...)
compile(c::Circuit)
expect(c::Circuit, state::String)
expect(c::Circuit, in_state::String, out_state::String)
```

## I/O and conversion
```@docs
op_to_strings(o::Operator)
get_coeffs(o::AbstractOperator)
set_coeffs(o::AbstractOperator, coeffs::Vector{T}) where {T<:Number}
Matrix(o::Operator)
SparseArrays.sparse(pauli::PauliString)
SparseArrays.sparse(o::Operator)
get_coeff(o::Operator{P}, p::P) where {P}
get_pauli(o::Operator, i::Int)
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
compress(o::Operator)
xcount(p::PauliString)
ycount(p::PauliString)
zcount(p::PauliString)
pauli_weight(p::PauliString)
```

## Index

```@index
```
