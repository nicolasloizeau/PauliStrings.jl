# Symbolic Calculations with PauliStrings (@id paulistrings_symbolics.md)

This tutorial demonstrates how to integrate [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl) into
[PauliStrings.jl](https://github.com/nicolasloizeau/PauliStrings.jl) for fully symbolic Pauli-operator
algebra—*without* introducing hard dependencies on Symbolics.

## 1. Setup

Install Symbolics.jl** only if you need symbolic functionality:

```julia
julia> using PauliStrings
julia> using Symbolics
```

## 2. Defining Symbolic Operators

```jldoctest symbolicspauli
julia> using PauliStrings # hide

julia> using PauliStrings: paulistringtype, commutator # hide

julia> using Symbolics; import Symbolics: simplify

julia> a, b, γ, Δ = @variables a b γ Δ
4-element Vector{Num}:
 a
 b
 γ
 Δ

julia> N = 2
2

julia> o = Operator{paulistringtype(N), Complex{Num}}()

julia> o += a,  "ZX" # ZX has 0 Y’s → (1i)^0 = 1
a ZX

julia> o += b,  "XY" # XY has 1 Y   → (1i)^1 = im
a ZX
b*im XY
```


## 3. Symbolic Arithmetic

All symbolic terms are preserved—no numeric cutoffs or boolean branches:

```jldoctest symbolicspauli
julia> o2 = o + o
2a ZX
2b*im XY


julia> o3 = o * o
a^2 + b^2 11
-2a*b*im YZ
```

## 4. Simplification of Coefficients

Use `simplify` to apply `Symbolics.simplify` per coefficient and combine like Pauli strings:

```jldoctest symbolicspauli
julia> o4 = a*b*o*(a*o + b*(o*o) - a^2*(o*o) - a*o)
0 11
(a^2)*((a^2 + b^2)*b - (a^2)*(a^2 + b^2))*b - a*(-2a*(b^2) + 2(a^3)*b)*(b^2) ZX
0 YZ
(-(a^2)*(-2a*(b^2) + 2(a^3)*b)*b + a*((a^2 + b^2)*b - (a^2)*(a^2 + b^2))*(b^2))*im XY


julia> os = simplify(o4)
(a^4)*(b^2) + 3(a^2)*(b^4) - (a^6)*b - 3(a^4)*(b^3) ZX
(3(a^3)*(b^3) + a*(b^5) - 3(a^5)*(b^2) - (a^3)*(b^4))*im XY
```

## 5. Parameterized Hamiltonian

Construct a 3-qubit Hamiltonian, compute its commutator and simplify:

```jldoctest symbolicspauli
julia> a, b, c = @variables a b c
3-element Vector{Num}:
 a
 b
 c

julia> H = Operator{paulistringtype(3), Complex{Num}}()


julia> H += a, "X11"
a X11


julia> H += b, "Y1Z"
a X11
b*im Y1Z


julia> H += c, "Z11"
a X11
b*im Y1Z
c Z11


julia> C = commutator(H, H)
0 111
0 Z1Z
0 Y11
0 X1Z


julia> isempty(simplify(C).coeffs)
true
```

