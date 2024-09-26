
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
com(o1::OperatorTS1D, o2::OperatorTS1D; anti=false)
```

```@docs
diag(o::Operator)
diag(o::OperatorTS1D)
```

```@docs
trace(o::Operator)
trace(o::OperatorTS1D)
```

```@docs
opnorm(o::Operator)
opnorm(o::OperatorTS1D)
```

```@docs
dagger(o::Operator)
dagger(o::OperatorTS1D)
```

```@docs
ptrace(o::Operator, keep::Vector{Int})
```


## Power and moments
```@docs
Base.:^(o::Operator, k::Int)
Base.:^(o::OperatorTS1D, k::Int)
```
```@docs
oppow(o::Operator, k::Int)
oppow(o::OperatorTS1D, k::Int)
```
```@docs
trace_product(o1::Operator, o2::Operator; scale=0)
trace_product(o1::OperatorTS1D, o2::OperatorTS1D; scale=0)
```
```@docs
trace_product(A::Operator, k::Int, B::Operator, l::Int; scale=0)
trace_product(A::OperatorTS1D, k::Int, B::OperatorTS1D, l::Int; scale=0)
```
```@docs
trace_product(A::Operator, k::Int; scale=0)
trace_product(A::OperatorTS1D, k::Int; scale=0)
```
```@docs
moments(H::Operator, kmax::Int; start=1, scale=0)
moments(H::OperatorTS1D, kmax::Int; start=1, scale=0)
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
k_local_part(o::Operator, k::Int)
```

## Algorithms
```@docs
lanczos(H::Operator, O::Operator, steps::Int, nterms::Int; keepnorm=true, maxlength=1000, localop=false)
```

```@docs
rk4(H::Operator, O::Operator, dt::Real; hbar::Real=1, heisenberg=false, M=2^20, keep::Operator=Operator(N))
```

```@docs
rk4(H::Function, O::Operator, dt::Real, t::Real; hbar::Real=1, heisenberg=false)
```

## Tools
```@docs
compress(o::Operator)
compress(o::OperatorTS1D)
```
```@docs
op_to_strings(o::Operator)
```
```@docs
op_to_dense(o::Operator)
```
```@docs
shift_left(O::Operator)
```

## Index

```@index
```
