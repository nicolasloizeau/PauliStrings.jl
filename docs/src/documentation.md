
# Documentation


## Basics
```@docs
Operator(N::Int)
```

```@docs
Operator(o::OperatorTS1D)
```

```@docs
OperatorTS1D(o::Operator; full=true)
```

```@docs
Base.length(o::Operator)
```

```@docs
eye(N::Int)
```

```@docs
Base.:+(o::Operator, args::Tuple{Number, Vararg{Any}})
Base.:+(o::Operator, args::Tuple{Vararg{Any}})
Base.:+(o::Operator, term::Tuple{Number, String})
Base.:+(o::Operator, term::String)
```


## Operations


```@docs
Base.:+(o1::Operator, o2::Operator)
Base.:+(o::Operator, a::Number)
Base.:+(a::Number, o::Operator)
```


```@docs
Base.:*(o1::Operator, o2::Operator)
Base.:*(o::Operator, a::Number)
Base.:*(a::Number, o::Operator)
```

```@docs
Base.:-(o::Operator)
Base.:-(o1::Operator, o2::Operator)
Base.:-(o::Operator, a::Real)
Base.:-(a::Real, o::Operator)
```

```@docs
com(o1::Operator, o2::Operator; epsilon::Real=0, maxlength::Int=1000)
```

```@docs
diag(o::Operator)
```

```@docs
trace(o::Operator)
```

```@docs
opnorm(o::Operator)
```

```@docs
dagger(o::Operator)
```

```@docs
ptrace(o::Operator, keep::Vector{Int})
```


## Power and moments
```@docs
Base.:^(o::Operator, k::Int)
```
```@docs
oppow(o::Operator, k::Int)
```
```@docs
trace_product(o1::Operator, o2::Operator; scale=0)
```
```@docs
trace_product(A::Operator, k::Int, B::Operator, l::Int; scale=0)
```
```@docs
trace_product(A::Operator, k::Int; scale=0)
```
```@docs
moments(H::Operator, kmax::Int; start=1, scale=0)
```




## Random operators
```@docs
rand_local1(N::Int)
```
```@docs
rand_local2(N::Int)
```
```@docs
rand_local1_TS1D(N::Int)
```
```@docs
rand_local2_TS1D(N::Int)
```



## Truncation and noise
```@docs
truncate(o::Operator, N::Int; keepnorm::Bool = false)
```
```@docs
trim(o::Operator, N::Int; keepnorm::Bool = false, keep::Operator=Operator(N))
```
```@docs
prune(o::Operator, alpha::Real; keepnorm::Bool = false)
```
```@docs
cutoff(o::Operator, epsilon::Real; keepnorm::Bool = false)
```
```@docs
add_noise(o::Operator, g::Real)
```
```@docs
k_local_part(o::Operator, k::Int; atmost=false)
```

## Algorithms
```@docs
lanczos(H::Operator, O::Operator, steps::Int, nterms::Int; keepnorm=true, maxlength=1000, returnOn=false)
lanczos(H::OperatorTS1D, O::OperatorTS1D, steps::Int, nterms::Int; keepnorm=true, maxlength=1000, returnOn=false)
```

```@docs
rk4(H::Operator, O::Operator, dt::Real; hbar::Real=1, heisenberg=false, M=2^20, keep::Operator=Operator(N))
```

```@docs
rk4(H::Function, O::Operator, dt::Real, t::Real; hbar::Real=1, heisenberg=false)
```

```@docs
equivalence_class(A1::Union{Operator64,Operator128}, H::Union{Operator64,Operator128})
```

```@docs
frustration_graph(o::Operator)
```
## Construction
```@docs
all_strings(N::Int)
```
```@docs
all_k_local(N::Int, k::Int)
```
```@docs
all_x(N::Int)
```
```@docs
all_y(N::Int)
```
```@docs
all_z(N::Int)
```
```@docs
set_coefs(o::Operator, coefs::Vector{T}) where T <: Number
```

## Circuits

```@docs
Circuit(N::Int; max_strings=2^30, noise_amplitude=0)
```
```@docs
push!(c::Circuit, gate::String, sites::Real...)
```
```@docs
pushfirst!(c::Circuit, gate::String, sites::Real...)
```

```@docs
XGate(N::Int, i::Int)
```
```@docs
PhaseGate(N::Int, i::Int, theta::Real)
```

```@docs
CXGate(N::Int, i::Int, j::Int)
```

```@docs
SwapGate(N::Int, i::Int, j::Int)
```
```@docs
CSXGate(N::Int, i::Int, j::Int)
```

```@docs
CCXGate(N::Int, i::Int, j::Int, k::Int)
```
```@docs
MCZGate(N::Int, sites::Int...)
```
```@docs
grover_diffusion(N::Int, sites::Int...)
```
```@docs
compile(c::Circuit)
```
```@docs
expect(c::Circuit, state::String)
```
```@docs
expect(c::Circuit, in_state::String, out_state::String)
```




## Tools
```@docs
compress(o::Operator)
```
```@docs
op_to_strings(o::Operator)
```
```@docs
get_coefs(o::Operator)
```
```@docs
op_to_dense(o::Operator)
```
```@docs
shift_left(O::Operator)
```
```@docs
xcount(v::Unsigned, w::Unsigned)
```
```@docs
ycount(v::Unsigned, w::Unsigned)
```
```@docs
zcount(v::Unsigned, w::Unsigned)
```

```@docs
is_ts(o::Operator)
```

```@docs
get_coef(o::Operator, v::Unsigned, w::Unsigned)
```

```@docs
get_pauli(o::Operator, i::Int)
```


## Low level

```@docs
com(v1::Unsigned, w1::Unsigned, v2::Unsigned, w2::Unsigned)
```

```@docs
string_to_vw(pauli::String)
```

```@docs
vw_to_string(v::Unsigned, w::Unsigned, N::Int)
```

```@docs
vw_in_o(v::Unsigned, w::Unsigned, o::Operator)
```

```@docs
Base.push!(o::Operator, c::Number, v::Unsigned, w::Unsigned)
```

## Index

```@index
```
