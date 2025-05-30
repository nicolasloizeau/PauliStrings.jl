# Symbolic Calculations with PauliStrings (@id paulistrings_symbolics.md)

This tutorial demonstrates how to integrate [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl) into
[PauliStrings.jl](https://github.com/nicolasloizeau/PauliStrings.jl) for fully symbolic Pauli-operator
algebra—without introducing hard or weak dependencies on [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl).

## Setup

Install `Symbolics.jl` only if you need symbolic functionality:

```julia
julia> using PauliStrings
julia> using Symbolics
```

## Defining Symbolic Operators

```jldoctest symbolicexamples
julia> using PauliStrings # hide

julia> using PauliStrings: paulistringtype, commutator # hide

julia> using Symbolics; import Symbolics: simplify, show

julia> a, b, γ, Δ = @variables a b γ Δ
4-element Vector{Num}:
 a
 b
 γ
 Δ

julia> N = 2
2

julia> o = Operator{paulistringtype(N), Complex{Num}}();

julia> o += a,  "ZX"
(a) ZX

julia> o += b,  "XY"
(a) ZX
(b) XY
```


## Symbolic Arithmetic

All symbolic terms are preserved—no numeric cutoffs or boolean branches:

```jldoctest symbolicexamples
julia> o2 = o + o
(2a) ZX
(2b) XY

julia> o3 = o * o
(a^2 + b^2) 11
(-2a*b) YZ
```

## Simplification of Coefficients

Use `simplify` to apply `Symbolics.simplify` per coefficient and combine like Pauli strings:

```jldoctest symbolicexamples
julia> o4 = a*b*o*(a*o + b*(o*o) - a^2*(o*o) - a*o)
(0.0) 11
((a^2)*((a^2 + b^2)*b - (a^2)*(a^2 + b^2))*b - a*(-2a*(b^2) + 2(a^3)*b)*(b^2)) ZX
(0.0) YZ
(-(a^2)*(-2a*(b^2) + 2(a^3)*b)*b + a*((a^2 + b^2)*b - (a^2)*(a^2 + b^2))*(b^2)) XY


julia> os = simplify(o4)
(0.0) 11
((a^2)*((a^2 + b^2)*b - (a^2)*(a^2 + b^2))*b - a*(-2a*(b^2) + 2(a^3)*b)*(b^2)) ZX
(0.0) YZ
(-(a^2)*(-2a*(b^2) + 2(a^3)*b)*b + a*((a^2 + b^2)*b - (a^2)*(a^2 + b^2))*(b^2)) XY
```

## Parameterized Hamiltonian

Construct a 3-qubit Hamiltonian, compute its commutator and simplify:

```jldoctest symbolicexamples
julia> a, b, c = @variables a b c
3-element Vector{Num}:
 a
 b
 c

julia> H = Operator{paulistringtype(3), Complex{Num}}();

julia> H += a, "X11"
(a) X11

julia> H += b, "Y1Z"
(a) X11
(b) Y1Z

julia> H += c, "Z11"
(a) X11
(b) Y1Z
(c) Z11

julia> C = commutator(H, H)
(0.0) Z1Z
(0.0) Y11
(0.0) X1Z


julia> iszero(simplify(C).coeffs)
true
```
