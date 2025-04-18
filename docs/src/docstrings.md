
# Documentation



## Basics
```@docs
Operator(N::Int)
OperatorTS1D(o::Operator; full=true)
Operator(o::OperatorTS1D)
```


## Truncation
```@docs
add_noise(o::AbstractOperator, g::Real)
truncate(o::AbstractOperator, max_lenght::Int; keepnorm::Bool=false)
k_local_part(o::AbstractOperator, k::Int; atmost=false)
trim(o::AbstractOperator, max_strings::Int; keepnorm::Bool=false, keep::Operator=Operator(0))
prune(o::AbstractOperator, alpha::Real; keepnorm::Bool=false)
cutoff(o::AbstractOperator, epsilon::Real; keepnorm::Bool=false)
```

## Algorithms
```@docs
rk4(H::Operator, O::Operator, dt::Real; hbar::Real=1, heisenberg=true, M=2^20, keep::Operator=Operator(0))
rk4(H::Function, O::Operator, dt::Real, t::Real; hbar::Real=1, heisenberg=true)
rk4_lindblad(H::Operator, O::Operator, dt::Real, L; hbar::Real=1, heisenberg=true, M=2^20, keep::Operator=Operator(0), gamma=[])
lanczos(H::AbstractOperator, O::AbstractOperator, steps::Int, nterms::Int; keepnorm=true, maxlength=1000, returnOn=false, observer=false)
```

## Operations
```@docs
Base.:+(o1::O, o2::O) where {O<:AbstractOperator}
Base.:-(o1::O, o2::O) where {O<:AbstractOperator}
Base.:*(o1::Operator, o2::Operator; kwargs...)
Base.:/(o::AbstractOperator, a::Number)
trace_product(o1::Operator, o2::Operator; scale=0)
compress(o::AbstractOperator)
trace(o::Operator; normalize=false)
diag(o::AbstractOperator)
opnorm(o::AbstractOperator; normalize=false)
dagger(o::AbstractOperator)
ptrace(o::AbstractOperator, keep::Vector{Int})
```
## Index

```@index
```
