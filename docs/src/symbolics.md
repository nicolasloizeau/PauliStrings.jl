# Symbolics

In this tutorial, we show how to work with symbolic operators using the [Symbolics.jl](https://docs.sciml.ai/Symbolics/stable/) package.

Before running any calculation, we first must install `Symbolics.jl`, and import it:
```julia
using Symbolics
using PauliStrings
```
To define a symbolic operator, we use the [`OperatorSymbolic`](@ref) constructor:
```julia
N = 2
H = OperatorSymbolic(N)
```
Now let's create a two-site Ising model, where the value of the transverse field is a symbolic variable `h`
```julia
@variables h
H += "Z",1,"Z",2
H += h, "X", 1
H += h, "X", 2
```
We can perform the usual operatons supported for [`Operator`](@ref)
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
To simplify the expressions for each coefficient, we use `Symbolics.simplify` via [`simplify_op`](@ref)
```julia
H2 = H^2
println(H2)
println(simplify_op(H2))
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
println(simplify_op(H3))
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
As a final example, let's perform one step of time evolution using [`rk4`](@ref) 
```julia
O = OperatorSymbolic(N)
O += "X", 1 # Operator to evolve
dt = 0.01
O = rk4(H, O, dt; heisenberg=true, M=16)
println(O)
```
```
(2.666666666666667e-6(h^2)) ZY
(0.0016666666666666666(-7.9996 + 0.0008(h^2) + 2(-1.9998 + 0.0004(h^2)))) YZ
(1.3333333333333335e-8(h^2)) 1X
(1 + 0.0016666666666666666(-0.08 + 0.02(-1.9998 + 0.0004(h^2)))) X1
(-0.0016666666666666666(0.08h - 8.000000000000001e-6(h^3) - 0.02h*(-1.9998 + 0.0004(h^2)))) YY
(0.0016666666666666666(0.08h - 8.000000000000001e-6(h^3) - 0.02h*(-1.9998 + 0.0004(h^2)))) ZZ
```
To substitute the variables for concrete numerical values we use [`substitute_op`](@ref)
```julia
op = substitute_op(O, Dict(h=>0.5))
println(typeof(op))
println(op)
```
```
(6.667e-7 - 0.0im) ZY
(-0.019998 - 0.0im) YZ
(3.3e-9 + 0.0im) 1X
(0.99980001 + 0.0im) X1
(-9.99933e-5 - 0.0im) YY
(9.99933e-5 + 0.0im) ZZ
```
