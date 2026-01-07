# Symbolics

In this tutorial, we show how to work with symbolic operators using  [MathLink.jl](https://github.com/JuliaInterop/MathLink.jl) and [Symbolics.jl](https://docs.sciml.ai/Symbolics/stable/). The main idea is that the coefficients of the `Operator` can be arbitrary `Number` types, we can thus use symbolic numbers from `Symbolics.jl` and expressions from `MathLink.jl` to create symbolic operators.
Two `PauliStrings.jl` extensions are provided to facilitate the construction and manipulation of symbolic operators using `MathLink.jl` and `Symbolics.jl`.


## Mathlink.jl (Wolfram Mathematica)
If you have access to Wolfram Mathematica, we recommend using `MathLink.jl` to create symbolic operators. The `MathLink` `PauliStrings` extension provides a `MathLinkNumber` type that wraps Mathematica symbolic expressions into `Numbers`.

Let's create a simple symbolic operator:

```julia
using MathLink
using PauliStrings
a = W`A`
b = W`B`
O = OperatorMathLink(2)
O += a, "X", 1
O += b, "Z", 1, "Z", 2
println(typeof(O))
```
```
Operator{PauliString{2, UInt8}, MathLinkPauliStringsExt.MathLinkNumber}
```

We can perform the usual operations supported for [`Operator`](@ref)
```julia
println(O)
println(O + O + 1)
println(trace_product(O, 4))
```
```julia
(A) X1
(B) ZZ

(Complex[1.0, 0.0]) 11
(Times[2, A]) X1
(Times[2, B]) ZZ

Times[4.0, Power[Plus[Power[A, 2], Power[B, 2]], 2]]
```

It is possible to simplify a `MathLinkNumber` using `simplify` with assumptions:
```julia
println(simplify(norm(O^2)))
println(simplify(norm(O^2); assumptions = W`Assumptions -> {A>0, B>0}`))
```
Operators can be similarly simplified using `simplify_operator`.

Lets now compute a few lanczos coefficients for the initial operator $Z_1$ in the transverse field Ising model with symbolic transverse field `h`:
Define the Ising Hamiltonian function:
```julia
function ising(h)
    H = OperatorMathLink(N)
    for i in 1:N
        H += h, "X", i
    end
    for i in 1:N
        H += "Z", i, "Z", mod1(i + 1, N)
    end
    return H
end
```
Initialize two operators :
```julia
N = 10 #system size

# Initial operator O = X_1
O = OperatorMathLink(N)
O += 1, "X", 1

# Initialize a symbolic hamiltonian with a symbolic parameter h
h = W`h`
H = ising(h)
```

Now we can compute the lanczos coefficients with simplification assumptions:
```julia
# Define assumptions for simplification
assumptions = W`Assumptions -> h > 0`

# Compute Lanczos coefficients
bn = lanczos(H, O, 5; assumptions=assumptions)
for b in bn
    println(b)
end
```
```julia
Times[2, Power[2, Rational[1, 2]]]
Power[Plus[8, Times[8, Power[h, 2]]], Rational[1, 2]]
Times[2, h, Power[Times[Power[Plus[1, Power[h, 2]], -1], Plus[5, Times[2, Power[h, 2]]]], Rational[1, 2]]]
Power[Times[Power[Plus[5, Times[7, Power[h, 2]], Times[2, Power[h, 4]]], -1], Plus[64, Times[96, Power[h, 2]], Times[68, Power[h, 4]]]], Rational[1, 2]]
Times[2, Power[Times[Power[Plus[80, Times[152, Power[h, 2]], Times[133, Power[h, 4]], Times[34, Power[h, 6]]], -1], Plus[64, Times[516, Power[h, 2]], Times[587, Power[h, 4]], Times[231, Power[h, 6]], Times[96, Power[h, 8]]]], Rational[1, 2]]]
```
These expressions can be easilly imported back into Mathematica for further manipulation, or conversion into LaTeX with [`TeXForm`](https://reference.wolfram.com/language/ref/TeXForm.html):

```math
\left\{2 \sqrt{2},\sqrt{8 h^2+8},2 h \sqrt{\frac{2 h^2+5}{h^2+1}},\sqrt{\frac{68 h^4+96 h^2+64}{2 h^4+7 h^2+5}},2
   \sqrt{\frac{96 h^8+231 h^6+587 h^4+516 h^2+64}{34 h^6+133 h^4+152 h^2+80}}\right\}
```

The `MathLink` expression can be extracted from a `MathLinkNumber` doing `x.expression` where `x` is a `MathLinkNumber`.


## Symbolics.jl
The second option is to use `Symbolics.jl`.

### Construction and basic operations

First import the necessary packages:
```julia
using Symbolics
using PauliStrings
```


Initialize an empty 2-qubit symbolic operator:
```julia
N = 2
H = OperatorSymbolics(N)
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

### Simplification of symbolic operators
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
O1 = OperatorSymbolics(N)
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

### Substituting variables with numericald values

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
