
# Documentation



## Basics
```@docs
Operator(N::Int)
OperatorTS1D(o::Operator; full=true)
OperatorTS2D(o::Operator, L1::Int; full=true)
Operator(o::OperatorTS1D)
Operator(o::OperatorTS2D)
```


## Truncation and noise
```@docs
add_noise(o::AbstractOperator, g::Real)
add_noise(o::AbstractOperator, g::AbstractVector{<:Real})
add_dephasing_noise(o::AbstractOperator, g::Real)
add_dephasing_noise(o::AbstractOperator, g::AbstractVector{<:Real})
truncate(o::AbstractOperator, max_lenght::Int; keepnorm::Bool=false)
k_local_part(o::AbstractOperator, k::Int; atmost=false)
trim(o::AbstractOperator, max_strings::Int; keepnorm::Bool=false, keep::Operator=Operator(0))
prune(o::AbstractOperator, alpha::Real; keepnorm::Bool=false)
cutoff(o::AbstractOperator, epsilon::Real; keepnorm::Bool=false)
```

## Algorithms
```@docs
rk4(H::AbstractOperator, O::AbstractOperator, dt::Real; hbar::Real=1, heisenberg=true, M=2^20, keep::Operator=Operator(0))
rk4(H::Function, O::AbstractOperator, dt::Real, t::Real; hbar::Real=1, heisenberg=true, M=2^20, keep::Operator=Operator(0))
rk4_lindblad(H::AbstractOperator, O::AbstractOperator, dt::Real, L; hbar::Real=1, heisenberg=true, M=2^20, keep::Operator=Operator(0), gamma=[])
lanczos(H::AbstractOperator, O::AbstractOperator, steps::Int, nterms::Int; keepnorm=true, maxlength=1000, returnOn=false, observer=false)
```

## Operations
```@docs
Base.:+(o1::O, o2::O) where {O<:AbstractOperator}
Base.:-(o1::O, o2::O) where {O<:AbstractOperator}
Base.:*(o1::Operator, o2::Operator; kwargs...)
commutator(o1::Operator, o2::Operator; kwargs...)
anticommutator(o1::Operator, o2::Operator; kwargs...)
Base.:/(o::AbstractOperator, a::Number)
compress(o::AbstractOperator)
trace(o::Operator; normalize=false)
diag(o::AbstractOperator)
opnorm(o::AbstractOperator; normalize=false)
dagger(o::AbstractOperator)
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
all_strings(N::Int)
all_k_local(N::Int, k::Int)
all_x(N::Int)
all_y(N::Int)
all_z(N::Int)
set_coeffs(o::AbstractOperator, coeffs::Vector{T}) where {T<:Number}
majorana(N::Int, k::Int)
string_2d(args::Tuple{Vararg{Any}}, L1::Int, L2::Int; pbc=false)
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

## Other tools
```@docs
compress(o::Operator)
op_to_strings(o::Operator)
get_coeffs(o::Operator)
op_to_dense(o::Operator)
Matrix(o::Operator)
shift_left(O::Operator)
rotate(o::Operator, r::Int)
xcount(p::PauliString)
ycount(p::PauliString)
zcount(p::PauliString)
is_ts(o::Operator)
is_ts2d(o::Operator,L1::Int)
get_coeff(o::Operator{P}, p::P) where {P}
get_pauli(o::Operator, i::Int)
```

## Index

```@index
```
