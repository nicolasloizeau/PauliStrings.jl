# Advanced usage and data structures

This page documents the internal representation used by PauliStrings.jl.
Most users can work with [`Operator`](@ref), `PauliString`, and the
high-level constructors, but understanding the storage model is useful when
writing performance-sensitive code or contributing new operations.

## Binary encoding

Each `PauliString` stores a string of Pauli operators as two unsigned
integers, `v` and `w`. Site `i` is encoded in bit `i - 1` of both integers,
so the first qubit is stored in the least significant bit.

The current code uses the following convention:

| Pauli | `v_i` | `w_i` |
|:------|:------|:------|
| `1` or `I` | `0` | `0` |
| `X` | `0` | `1` |
| `Z` | `1` | `0` |
| `Y` | `1` | `1` |

For example, `PauliString("Y1Z1")` has `v = 0b0101` and `w = 0b0001`.
The leftmost printed operator is site 1, which corresponds to the least
significant bit:

```@example advanced
using PauliStrings

p = PauliString("Y1Z1")
(p.v, p.w)
```

The same representation makes common predicates cheap bit operations:

| Quantity | Bit expression |
|:---------|:---------------|
| X count | `count_ones(~p.v & p.w)` |
| Y count | `count_ones(p.v & p.w)` |
| Z count | `count_ones(p.v & ~p.w)` |
| Pauli weight | `count_ones(p.v | p.w)` |

## Integer types and size limits

The integer type is part of the concrete Pauli string type:

```@example advanced
typeof(PauliString("XYZ1"))
```

The default type for `N` qubits is selected by `paulistringtype(N)`. Small
systems use Julia's native unsigned integers, while larger systems use
`BitIntegers`-provided fixed-width integers:

```@example advanced
paulistringtype(8)
paulistringtype(64)
paulistringtype(128)
paulistringtype(256)
```

The package documentation advertises support up to 1024 qubits. Within that
range, the backing width is the next power of two that can hold all sites
(`UInt8`, `UInt16`, `UInt32`, `UInt64`, `UInt128`, `UInt256`, and so on).

You can request a concrete string type explicitly when needed:

```@example advanced
P = paulistringtype(64)
P("X" * "1"^63)
```

## Stored phases and coefficients

The matrix identity

```math
Y = i X Z
```

is reflected in the storage convention. A bare `PauliString` has no
user coefficient, but when it is viewed through the `AbstractOperator`
interface its value is `(1im)^ycount(p)`.

[`Operator`](@ref) therefore stores coefficients in `O.coeffs` with the
internal Y phase included. The user-facing coefficient can be recovered with
[`get_coeffs`](@ref), which divides that phase back out.

```@example advanced
O = Operator("Y")

O.coeffs
get_coeffs(O)
```

This is why code that needs physical coefficients should prefer
[`get_coeffs`](@ref), [`get_coeff`](@ref), or [`op_to_strings`](@ref) instead
of reading `O.coeffs` directly. Direct access is useful for low-level kernels,
where keeping the phase in the coefficient avoids recomputing it at every
operation.

## Pauli algebra as boolean algebra

Multiplication of the Pauli labels is mostly XOR on the two bitstrings. If
`p1` and `p2` share the same concrete Pauli string type, the product support is

```julia
v = p1.v ⊻ p2.v
w = p1.w ⊻ p2.w
```

The remaining work is a sign or commutation factor. For instance, `X` and `Z`
anticommute on the same site:

```@example advanced
x = PauliString("X")
z = PauliString("Z")

x ⊻ z
commutator(x, z)
```

The first value returned by `commutator(x, z)` is the XOR-produced Pauli
string, and the second value is the integer prefactor used by the operator
kernel. For full operators, use the high-level [`commutator`](@ref) function:

```@example advanced
commutator(Operator("X"), Operator("Z"))
```

This bit representation is also why the specialized commutator is faster than
forming `A * B - B * A`: it computes only the terms that survive the
commutator and accumulates them directly.

## Type hierarchy

PauliStrings.jl uses a small hierarchy so that ordinary strings, sums of
strings, and translation-symmetric strings can share algebra kernels.

| Type | Role | Main stored fields |
|:-----|:-----|:-------------------|
| `AbstractOperator` | Supertype for objects that can be iterated as Pauli strings with coefficients. | Interface only. |
| `AbstractPauliString` | Supertype for a single Pauli string-like object. It is also an `AbstractOperator` with length 1. | Interface only. |
| `PauliString{N,T}` | A concrete Pauli string on `N` sites. | `v::T`, `w::T`. |
| `Operator{P,T}` | A sum of Pauli strings of type `P` with scalar type `T`. | `strings::Vector{P}`, `coeffs::Vector{T}`. |
| `PauliStringTS{Ls,Ps,T}` | A representative of a translation-symmetric sum on a lattice with size tuple `Ls` and periodic flags `Ps`. | `v::T`, `w::T`. |
| `OperatorTS{Ls,Ps,U,T}` | Alias for `Operator{PauliStringTS{Ls,Ps,U},T}`. | Same storage as `Operator`. |

Use `PauliString` when you need a single string and want direct access
to bit operations. Use [`Operator`](@ref) for generic sums, Hamiltonians, and
observables. Use [`PauliStringTS`](@ref) or [`OperatorTS`](@ref) when the
problem has translation symmetry and you want to store only one representative
per translation orbit.

The helper functions `paulistringtype`, [`qubitlength`](@ref), and
`scalartype` let generic code discover the concrete string type,
number of qubits, and scalar type without reaching into fields directly.
