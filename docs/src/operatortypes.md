# Tutorial on Operator and Pauli-String Types

[PauliStrings.jl](https://github.com/nicolasloizeau/PauliStrings.jl) is a Julia
library designed for efficient representation and manipulation of quantum operators
built from [`PauliString`](@ref). It supports generic operators ([`Operator`](@ref)),
as well as translationally invariant ones for 1D systems ([`OperatorTS1D`](@ref)).
This tutorial highlights the different types of operators and Pauli strings supported
by the library.

## Type Hierarchy

```@raw html
<div class="mermaid">
flowchart TD
    subgraph Abstract Types
        AO[AbstractOperator]
        APS[AbstractPauliString]
    end

    subgraph Concrete Types
        PS[PauliString]
        OP[Operator]
        OP1D[OperatorTS1D]
    end

    AO --> APS
    APS --> PS
    AO --> OP
    AO --> OP1D
</div>
```

- **[`AbstractOperator`](@ref)**: Interface for all quantum operators, defines: [`paulistringtype`](@ref), [`qubitlength`](@ref), and [`scalartype`](@ref)

- **[`AbstractPauliString`](@ref)**: Specialization of [`AbstractOperator`](@ref) representing a single sequence of Pauli matrices.

- **[`PauliString`](@ref)**: Concrete N-qubit [`PauliString`](@ref), parameterized by an unsigned integer type `T`.

- **[`Operator`](@ref)**: Arbitrary linear combination of [`PauliString`](@ref)s with coefficients.

- **[`OperatorTS1D`](@ref)**: Translationally invariant 1D operator: repeats a generating set across a periodic lattice.

## PauliString{N,T}

A [`PauliString`](@ref) 

``P = \\bigotimes_{i=1}^{N} P_i,\\quad P_i \\in \\{I,\\,X,\\,Y,\\,Z\\}``

encodes an `N`-qubit tensor product of Pauli matrices (`I`, `X`, `Y`, `Z`).

Internally it is encoded by two bitmasks `v,w ∈T`, where

``v_i =\\begin{cases}
1, & \\text{if } P_i \\in \\{Z,\\,Y\\},\\\\
0, & \\text{otherwise},
\\end{cases}
\\quad
w_i =
\\begin{cases}
1, & \\text{if } P_i \\in \\{X,\\,Y\\},\\\\
0, & \\text{otherwise}.
\\end{cases}``

Here `i` indexes the `(i−1)`-th bit of the integer.

### Construction

```jldoctest example
julia> using PauliStrings # hide

julia> ps1 = PauliString("XZYI")
XZY1

julia> ps2 = PauliString{4,UInt16}("YXXZ")
YXXZ

julia> ps3 = PauliString{4,UInt8}(0b1010, 0b0101)
XZXZ

julia> qubitlength(ps1)
4

julia> ps1[2]
:Z

julia> p = PauliString("YIXZ")
Y1XZ

julia> p[1]
:Y
```

##  Operator

[`Operator`](@ref) encapsulates a sum of Pauli strings:

``O = \\sum_{k=1}^{M} \\alpha_k\\,P^{(k)},\\quad
P^{(k)} \\in \\text{PauliString},\\;
\\alpha_k \\in T``


This type supports arbitrary operator algebra on `N`-qubit Hilbert spaces.

### Constructors

```jldoctest example
julia> opA = Operator{PauliString{3,UInt8},ComplexF64}();

julia> op = Operator(3);

julia> op += 1.0, "XIZ"
(1.0 + 0.0im) X1Z


julia> op += -0.5im, "YXI"
(1.0 + 0.0im) X1Z
(0.0 - 0.5im) YX1


julia> v1, w1, v2, w2 = 0b010, 0b001, 0b011, 0b111;

julia> α1, α2 = 0.5, -1.2im;

julia> op2 = Operator(3, [v1,v2], [w1,w2], [α1,α2])
(0.5 + 0.0im) XZ1
(-0.0 + 1.2im) YYX
```

### Example
```jldoctest example
julia> A = Operator(3);

julia> A += 1.0, "XZI";

julia> A += -0.5im, "YIX"
(1.0 + 0.0im) XZ1
(0.0 - 0.5im) Y1X


julia> length(A)
2
```

## OperatorTS1D

A 1D translationally invariant operator on N qubits stores a generating
set of terms {P(k)} and repeats each term with periodic boundary conditions:

``O_{\\text{TI}}
= \\sum_{k=1}^{m} \\sum_{j=1}^{N}
  \\alpha_k \\;T^j\\bigl(P^{(k)}\\bigr)``

where `T` is the translation operator on the lattice.

```jldoctest example
julia> ti0 = OperatorTS1D(5);

julia> ti1 = OperatorTS1D("XY")
(1.0 - 0.0im) XY

julia> N, n, m, h  = 40, 2, 3, 2;

julia> Δ = cos(pi * n / m);

julia> H = Operator(N);

julia> H += 'X', 1, 'X', 2
(1.0 + 0.0im) XX11111111111111111111111111111111111111


julia> H += 'Y', 1, 'Y', 2
(1.0 - 0.0im) YY11111111111111111111111111111111111111
(1.0 + 0.0im) XX11111111111111111111111111111111111111


julia> H += Δ, 'Z', 1, 'Z', 2
(-0.5 + 0.0im) ZZ11111111111111111111111111111111111111
(1.0 - 0.0im) YY11111111111111111111111111111111111111
(1.0 + 0.0im) XX11111111111111111111111111111111111111


julia> H += h, 'Z', 1
(-0.5 + 0.0im) ZZ11111111111111111111111111111111111111
(2.0 + 0.0im) Z111111111111111111111111111111111111111
(1.0 - 0.0im) YY11111111111111111111111111111111111111
(1.0 + 0.0im) XX11111111111111111111111111111111111111


julia> op = OperatorTS1D(H, full=false)
(-0.5 + 0.0im) 11111111111111111111111111111111111111ZZ
(2.0 + 0.0im) 111111111111111111111111111111111111111Z
(1.0 + 0.0im) 11111111111111111111111111111111111111XX
(1.0 - 0.0im) 11111111111111111111111111111111111111YY
```
