# Symbolics

In this tutorial, we show how to work with symbolic operators using the [Symbolics.jl](https://docs.sciml.ai/Symbolics/stable/) package. Make sure you have `Symbolics.jl` installed.


## Construction and basic operations

First import the necessary packages:
```julia
using Symbolics
using PauliStrings
```


Initialize an empty 2-qubit symbolic operator:
```julia
N = 2
H = OperatorSymbolic(N)
println(typeof(H))
```
```
Operator{PauliString{2, UInt8}, Complex{Num}}
```
Now let's create a two-site Ising model, where the value of the transverse field is a symbolic variable `h`
```julia
@variables h
H += "Z", 1, "Z", 2
H += h, "X", 1
H += h, "X", 2
```
We can perform the usual operations supported for [`Operator`](@ref)
```julia
println(H)
println(H + H + 0.5)
println(trace_product(H, 4))
```
```
(h) 1X
(h) X1
(1.0) ZZ

(0.5) 11
(2h) 1X
(2h) X1
(2.0) ZZ

4.0(4(h^4) + (1 + 2(h^2))^2)
```

## Simplification of symbolic operators
`simplify_operator` can be used to simplify the coefficients of a symbolic operator.

Let's apply this to `H^2`:
```julia
H2 = H^2
println(H2)
println(simplify_operator(H2))
```
```
(1 + 2(h^2)) 11
(0.0) ZY
(0.0) YZ
(2(h^2)) XX

(1 + 2(h^2)) 11
(2(h^2)) XX
```
However, keep in mind that `simplify` doesn't reduce all the expressions:
```julia
H3 = H^3
println(H3)
println(simplify_operator(H3))
```
```
(2(h^3) + h*(1 + 2(h^2))) 1X
(-2.0(h^2)) YY
(2(h^3) + h*(1 + 2(h^2))) X1
(1 + 2(h^2)) ZZ

(2(h^3) + h*(1 + 2(h^2))) 1X
(-2.0(h^2)) YY
(2(h^3) + h*(1 + 2(h^2))) X1
(1 + 2(h^2)) ZZ
```
As a final example, let's calculate some commutators with another operator `O`
```julia
O1 = OperatorSymbolic(N)
O1 += "X", 1
O2 = commutator(H, O1)
O3 = commutator(H, O2)
println(O2)
println(O3)
```
```
(2.0im) YZ

(4.0h) YY
(4.0) X1
(-4h) ZZ
```

## Substituting variables with numericald values

To substitute the variables for concrete numerical values we can use `substitute_operator`
```julia
O = substitute_operator(O3, Dict(h => 0.5))
println(typeof(O))
println(O)
```
```
Operator{PauliString{2, UInt8}, ComplexF64}
(2.0 - 0.0im) YY
(4.0 + 0.0im) X1
(-2.0 + 0.0im) ZZ
```
