
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
add(o1::Operator, o2::Operator)
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
equivalence_class(A1::Operator, H::Operator)
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
xcount(v::Int, w::Int)
```
```@docs
ycount(v::Int, w::Int)
```
```@docs
zcount(v::Int, w::Int)
```

```@docs
is_ts(o::Operator)
```

```@docs
get_coef(o::Operator, v::Int, w::Int)
```

```@docs
get_pauli(o::Operator, i::Int)
```


## Low level

```@docs
com(v1::Int, w1::Int, v2::Int, w2::Int)
```

```@docs
string_to_vw(pauli::String)
```

```@docs
vw_to_string(v::Int, w::Int, N::Int)
```

```@docs
vw_in_o(v::Int, w::Int, o::Operator)
```

```@docs
Base.push!(o::Operator, c::Number, v::Int, w::Int)
```

## Index

```@index
```
