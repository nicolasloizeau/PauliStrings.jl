# Advanced usage and internals

This page is for users who already know how to build and manipulate
`Operator`s and want to understand the representation choices behind
PauliStrings.jl. The public constructors should be enough for most code, but
knowing the internals is useful when writing new algorithms, inspecting stored
terms, or contributing new operator types.

## Binary encoding of Pauli strings

A [`PauliString`](@ref) stores a string of Pauli matrices as two unsigned
integers, named `v` and `w`. Each qubit contributes one bit to each integer.
For qubit `i`, the pair `(v_i, w_i)` encodes the local Pauli matrix:

| Pauli | `v_i` | `w_i` |
|:------|:-----:|:-----:|
| `1` or `I` | 0 | 0 |
| `X` | 0 | 1 |
| `Y` | 1 | 1 |
| `Z` | 1 | 0 |

The first qubit is stored in the least significant bit. For example:

```@example advanced
using PauliStrings

p = PauliString("Y1Z1")
println(p)
println("v = ", p.v)
println("w = ", p.w)
```

The concrete integer type is chosen from the number of qubits. Small systems
use Julia's built-in unsigned integer types, while larger systems use wider
integer types provided through BitIntegers. For example, the default string
types use `UInt64` up to 64 qubits, `UInt128` up to 128 qubits, and wider
integer types such as `UInt256` for larger systems.

This compact representation makes equality, hashing, multiplication, and
commutation fast because they operate on machine words instead of on strings or
arrays of symbols. PauliStrings.jl currently supports strings up to 1024 qubits
through the default constructors.

## Coefficients and stored phases

An [`Operator`](@ref) stores terms in two vectors:

- `O.strings`: the Pauli strings.
- `O.coeffs`: the stored coefficients.

The stored coefficients include the phase convention used by the binary Pauli
encoding. This means `O.coeffs` is an internal representation, while
[`get_coeffs`](@ref) returns the coefficients as users normally write them in
front of printed Pauli strings.

```@example advanced
O = Operator(3)
O += 2, "Y", 1
O += -3im, "X", 2, "Z", 3

println("printed operator:")
println(O)
println("stored coeffs: ", O.coeffs)
println("user coeffs:   ", get_coeffs(O))
```

When you need to inspect or export coefficients, prefer [`get_coeffs`](@ref)
or [`get_coeff`](@ref). When you are extending low-level routines, keep the
stored phase convention intact and use the existing constructors and helpers
instead of modifying `coeffs` by hand.

## Pauli algebra as boolean algebra

Ignoring the phase, multiplication of two Pauli strings is bitwise XOR of the
`v` and `w` bitstrings:

```@example advanced
p = PauliString("X1Z")
q = PauliString("ZY1")
r = xor(p, q)

println(p)
println(q)
println(r)
```

The missing piece is the phase. Pauli multiplication is not commutative: for
example, `X * Z = -im * Y` while `Z * X = im * Y`. PauliStrings.jl's
commutator and anticommutator kernels compute both the XOR result and the
phase/sign factor:

```@example advanced
x = PauliString("X")
z = PauliString("Z")

println(commutator(Operator(x), z))
println(anticommutator(Operator(x), z))
```

The same rule extends term-by-term to operators:

```@example advanced
A = Operator("X1") + Operator("1Z")
B = Operator("ZX")

println(commutator(A, B))
```

For performance-sensitive code, prefer [`commutator`](@ref) and
[`anticommutator`](@ref) over expanding expressions such as `A * B - B * A`.
The dedicated kernels combine the boolean representation and coefficient phase
updates directly.

## Type hierarchy and when to use each type

PauliStrings.jl separates single strings from sums of strings, and also has
translation-symmetric variants.

| Type | Stores | Use when |
|:-----|:-------|:---------|
| [`AbstractOperator`](@ref) | Supertype for all Pauli-string based operators | Writing generic methods that work on full operators and single strings |
| [`AbstractPauliString`](@ref) | Supertype for single Pauli-string-like objects | Writing methods specialized to one string or one translation class |
| [`PauliString`](@ref) | One Pauli string on `N` qubits as `(v, w)` | You need a single concrete string, for example a basis term or local generator |
| [`Operator`](@ref) | A vector of strings and a vector of coefficients | You need a general finite sum of Pauli strings |
| [`PauliStringTS`](@ref) | One representative of a translation-symmetric sum | Your model is translation invariant and you want to store an orbit compactly |
| [`OperatorTS`](@ref) | A sum of translation-symmetric strings | You need a translation-symmetric operator and want operations to stay in that representation |

The most common workflow is to use [`Operator`](@ref) for arbitrary Hamiltonians
or observables:

```@example advanced
H = Operator(4)
H += "Z", 1, "Z", 2
H += "Z", 2, "Z", 3
H += 0.5, "X", 4

println(typeof(H))
println(H)
```

Use [`PauliString`](@ref) when a single term is the natural unit of work:

```@example advanced
p = PauliString("XX1Z")
println(typeof(p))
println(support(p))
```

Use [`PauliStringTS`](@ref) and [`OperatorTS`](@ref) when the calculation should
identify terms that differ only by lattice translations:

```@example advanced
pts = PauliStringTS{(4,)}("XX11")
ots = pts + PauliStringTS{(4,)}("Z111")

println(typeof(pts))
println(typeof(ots))
println(ots)
```

Translation-symmetric types are especially useful for periodic spin chains and
lattices, where many translated terms share the same coefficient and algebraic
role.

## Practical guidelines for contributors

- Use constructors such as `Operator(N)`, `Operator("X1Z")`, and
  `PauliString{N}(...)` to preserve internal invariants.
- Use [`get_coeffs`](@ref) for user-facing coefficients, and read `coeffs`
  directly only in low-level code that understands the stored phase.
- Keep methods generic over [`AbstractOperator`](@ref) when possible, and
  specialize on [`Operator`](@ref), [`PauliString`](@ref), or
  [`OperatorTS`](@ref) only when the representation matters.
- Prefer the existing algebra kernels (`commutator`, `anticommutator`,
  `compress`, `trace_product`) to manual loops unless the new loop is part of a
  measured optimization.
