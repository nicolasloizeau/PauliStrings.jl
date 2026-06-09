# Advanced usage, datatypes, and data structures

This page describes the internal representation used by PauliStrings.jl. It is
intended for users who want to inspect stored strings directly, write extension
methods, or reason about performance-critical operations.

## Binary encoding

A `PauliString` stores a length-`N` Pauli string as two unsigned integers, `v`
and `w`:

```julia
struct PauliString{N,T<:Unsigned} <: AbstractPauliString
    v::T
    w::T
end
```

The bit at site `i` is read from `v` and `w`. In the current implementation,
`v` stores the Z component and `w` stores the X component:

| Pauli | `v_i` | `w_i` |
|---|---:|---:|
| `1` or `I` | 0 | 0 |
| `X` | 0 | 1 |
| `Z` | 1 | 0 |
| `Y` | 1 | 1 |

For example, `"Y1Z1"` has a Y on site 1 and a Z on site 3, so
`v = 0b0101` and `w = 0b0001`:

```julia
using PauliStrings

p = PauliString{4}(0b0101, 0b0001)
println(p)      # Y1Z1
println(p.v)    # 5
println(p.w)    # 1
```

The site index is one-based, and site 1 is stored in the least significant bit.
That convention makes local bit operations cheap, but it also means that binary
literals are visually read from right to left with respect to site order.

The backing integer type is selected by `uinttype`, through `paulistringtype`.
For small systems the package uses Julia's built-in unsigned integer types; for
larger systems it uses `BitIntegers` integer types. Common size thresholds are:

| Qubits `N` | Default integer type |
|---:|---|
| `N <= 8` | `UInt8` |
| `9 <= N <= 16` | `UInt16` |
| `17 <= N <= 32` | `UInt32` |
| `33 <= N <= 64` | `UInt64` |
| `65 <= N <= 128` | `UInt128` |
| `129 <= N <= 256` | `BitIntegers.UInt256` |

The actual rule is to use the next power-of-two width that can hold `N` bits.
This keeps `xor`, `&`, `|`, shifts, `count_ones`, and hashing fast while still
supporting systems larger than 128 qubits.

## Phase and coefficient storage

An `Operator` stores a sum of Pauli strings as parallel vectors:

```julia
struct Operator{P<:AbstractPauliString,T<:Number} <: AbstractOperator
    strings::Vector{P}
    coeffs::Vector{T}
end
```

The public coefficient of a term is the number printed in front of its Pauli
string. Internally, `O.coeffs` also includes the phase associated with the number
of Y operators in the corresponding string:

```julia
internal_coeff = public_coeff * (1im)^ycount(string)
```

That is why `O.coeffs` can differ from [`get_coeffs`](@ref):

```julia
using PauliStrings

O = Operator(4)
O += 2, "1XXY"
O += 3, "11Z1"

O.coeffs       # internal storage, including the Y phase
get_coeffs(O) # user-facing coefficients
```

For `"1XXY"`, `ycount` is 1, so the stored coefficient is multiplied by `1im`.
When an operator is printed, exported with [`op_to_strings`](@ref), or queried
with [`get_coeffs`](@ref), this phase is divided out again.

This convention lets algebra kernels work with the compact `(v, w)` encoding
while preserving the usual Pauli matrix phases in user-facing coefficients.

## Pauli algebra as Boolean algebra

Most single-string algebra reduces to bit operations on `v` and `w`.

The Pauli support of a product is the bitwise xor of the two input strings:

```julia
p = xor(p1, p2)
```

The remaining scalar factor is a sign or zero computed from parity checks. The
core multiplication rule used by `prod(p1, p2)` is:

```julia
p = xor(p1, p2)
k = 1 - ((count_ones(p1.v & p2.w) & 1) << 1)
```

So `k` is `1` or `-1` depending on the parity of Z bits in the left operand that
cross X bits in the right operand. Operator multiplication then combines this
with the stored coefficients.

Commutators and anticommutators use the same xor result and only change the
scalar factor:

```julia
commutator(p1, p2)     # returns p, k where k is -2, 0, or 2
anticommutator(p1, p2) # returns p, k where k is 0 or 2
```

For example, on one qubit, `X` and `Z` anticommute:

```julia
X = PauliString("X")
Z = PauliString("Z")

p, k = commutator(X, Z)
println(p) # Y
println(k) # 2
```

As operators, this becomes the usual matrix-algebra result after the stored
Y phase is applied:

```julia
commutator(Operator("X"), Operator("Z")) # -2im * Y
```

For a quick mental model:

- `xor` decides the output Pauli string.
- `count_ones` parity decides whether terms commute or anticommute.
- `ycount` and `O.coeffs` carry the phase convention needed for printed and
  exported coefficients.

## Type hierarchy and storage

The main data model has two levels: single strings and sums of strings.
Translation-symmetric types reuse the same storage pattern, but canonicalize
strings under translations.

| Type | Supertype | Storage | Use when |
|---|---|---|---|
| `AbstractOperator` | none | Abstract interface for Pauli-string-like operators | Writing methods that should accept both single strings and sums |
| `AbstractPauliString` | `AbstractOperator` | Abstract interface for one encoded Pauli string | Writing methods specialized to individual Pauli strings |
| `PauliString{N,T}` | `AbstractPauliString` | Two unsigned integers, `v::T` and `w::T` | Representing one concrete N-qubit Pauli string |
| `Operator{P,T}` | `AbstractOperator` | `strings::Vector{P}` and `coeffs::Vector{T}` | Representing a linear combination of Pauli strings |
| `PauliStringTS{Ls,Ps,T}` | `AbstractPauliString` | Two unsigned integers for a canonical translation representative | Representing one translation-symmetric Pauli-string orbit |
| `OperatorTS{Ls,Ps,U,T}` | alias of `Operator{PauliStringTS{Ls,Ps,U},T}` | Vectors of translation-symmetric string representatives and coefficients | Representing a translation-symmetric linear combination |

In tree form:

```text
AbstractOperator
|-- AbstractPauliString
|   |-- PauliString{N,T}
|   `-- PauliStringTS{Ls,Ps,T}
`-- Operator{P,T}
    `-- OperatorTS{Ls,Ps,U,T} = Operator{PauliStringTS{Ls,Ps,U},T}
```

### `PauliString`

Use `PauliString` when you need one concrete string and want cheap bit-level
operations:

```julia
p = PauliString("X1YZ")
support(p)
pauli_weight(p)
```

The type parameter `N` stores the number of qubits at the type level, and the
integer type `T` stores the chosen bit width.

### `Operator`

Use `Operator` for sums, Hamiltonians, observables, and any expression with
more than one term:

```julia
H = Operator(4)
H += "X", 1, "X", 2
H += -1.05, "Z", 1
```

The strings in `H.strings` all have the same qubit length. Repeated strings can
be accumulated with [`compress`](@ref), and user-facing coefficients can be read
with [`get_coeffs`](@ref).

### `PauliStringTS`

Use `PauliStringTS` for a translation-symmetric orbit of a string. The type
parameter `Ls` gives the lattice shape, while `Ps` marks periodic dimensions:

```julia
p = PauliStringTS{(4,)}("X1Y1")
qubitsize(p)      # (4,)
periodicflags(p)  # (true,)
representative(p)
```

Only a canonical representative is stored. The full translated sum is handled
lazily by translation-symmetric operations.

### `OperatorTS`

`OperatorTS` is a type alias for an `Operator` whose string type is
`PauliStringTS`. Use it for translation-invariant Hamiltonians or observables:

```julia
H0 = Operator(4)
H0 += "X", 1, "X", 2

H = OperatorTS1D(H0; full=false)
representative(H)
resum(H)
```

[`representative`](@ref) returns the stored canonical terms. [`resum`](@ref)
expands the translation-symmetric representation into an ordinary `Operator` by
summing over translations.

## Practical guidance

- Use public constructors and `get_coeffs` unless you specifically need the
  internal representation.
- If you inspect `O.coeffs` directly, remember that Y phases are already folded
  into those values.
- Prefer `PauliString` for individual strings and `Operator` for sums.
- Prefer `OperatorTS1D`, `OperatorTS2D`, or `OperatorTS{Ls}` when translation
  symmetry is part of the mathematical object.
- When extending algebra operations, follow the existing kernels: xor for the
  output support, parity checks for signs/commutation, and `compress` or
  `setwith!` to accumulate repeated output strings.
