# Symbolic Calculations with PauliStrings (@id paulistrings_symbolics.md)

This tutorial demonstrates how to integrate [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl) into
[PauliStrings.jl](https://github.com/nicolasloizeau/PauliStrings.jl) for fully symbolic Pauli-operator
algebra—without introducing hard dependencies on [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl).

## Setup

Install `Symbolics.jl` only if you need symbolic functionality:

```julia
julia> using PauliStrings
julia> using Symbolics
```

## Defining Symbolic Operators

```julia
using PauliStrings # hide
using PauliStrings: paulistringtype, commutator # hide
using Symbolics; import Symbolics: simplify, show
a, b, γ, Δ = @variables a b γ Δ
N = 2
o = Operator{paulistringtype(N), Complex{Num}}()
o += a,  "ZX"
(a) ZX

o += b,  "XY"
(a) ZX
(b) XY
```


## Symbolic Arithmetic

All symbolic terms are preserved—no numeric cutoffs or boolean branches:

```julia
o2 = o + o
(2a) ZX
(2b) XY

o3 = o * o
(a^2 + b^2) 11
(-2a*b) YZ
```

## Simplification of Coefficients

Use `simplify` to apply `Symbolics.simplify` per coefficient and combine like Pauli strings:

```julia
o4 = a*b*o*(a*o + b*(o*o) - a^2*(o*o) - a*o)
(0.0) 11
((a^2)*((a^2 + b^2)*b - (a^2)*(a^2 + b^2))*b - a*(-2a*(b^2) + 2(a^3)*b)*(b^2)) ZX
(0.0) YZ
(-(a^2)*(-2a*(b^2) + 2(a^3)*b)*b + a*((a^2 + b^2)*b - (a^2)*(a^2 + b^2))*(b^2)) XY

os = simplify(o4)
((a^4)*(b^2) + 3(a^2)*(b^4) - (a^6)*b - 3(a^4)*(b^3)) ZX
(3(a^3)*(b^3) + a*(b^5) - 3(a^5)*(b^2) - (a^3)*(b^4)) XY
```

## Parameterized Hamiltonian

Construct a 3-qubit Hamiltonian, compute its commutator and simplify:

```julia
a, b, c = @variables a b c
H = Operator{paulistringtype(3), Complex{Num}}()
H += a, "X11"
H += b, "Y1Z"
H += c, "Z11"

C = commutator(H, H)
(0.0) 111
(0.0) Z1Z
(0.0) Y11
(0.0) X1Z

isempty(simplify(C).coeffs)
true
```
