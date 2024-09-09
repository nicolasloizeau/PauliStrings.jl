
# Documentation


## Basics
```@docs
Operator(N::Int)
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

## Random operators
```@docs
rand_local1(N::Int)
```
```@docs
rand_local2(N::Int)
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
```
```@docs
op_to_strings(o::Operator)
```

```@docs
shift_left(O::Operator)
```

## Index

```@index
```
