# Manipulating Pauli strings

Here is a short tutorial on the `PauliString` type that encodes a single Pauli string.
It's more efficient to use this type than an `Operator` that stores a single string..

## Initializing single Pauli strings

The type `PauliStrings{N, T}` represents a single Pauli string acting on `N` spins. `T` is the integer type used to store the binary representation of the string.
A single Pauli string can be initialized from a string:
```@example singlestrings
using PauliStrings
s = PauliString("1Z1Y")
println(typeof(s))
println(s)
```

or from a list of operators and their site index. In this case, we need to specify the number of sites `N`, here `N=20`:
```@example singlestrings
s = PauliString{20}("Z", 2, "Y", 10)
println(typeof(s))
println(s)
```

THe binary representation consists of two integers  (bitstrings) $v$, $w$. The couple $v_i$, $w_i$  encodes the operator acting on site `i`. The mapping is as follows: X $\to$ (1,0), Y $\to$ (1,1), Z $\to$ (0,1), 1 $\to$ (0,0).

For example, the string "Y1Z1" is represented by $v = 0101_2 = 5$ and $w = 0001_2 = 1$:
```@example singlestrings
s = PauliString{4}(5, 1)
println(s)
println(s.v, " ", s.w)
```
(note that the strings are printed with the least significant site on the left while the least significant bit is on the right in the binary representation).

## Iterate over all strings of an operator
The strings are stored in `o.strings` and the coefficients in `o.coeffs`. Note that `o.coeffs` also stores a phase that counts the number of Y in the string, to get a list of coefficients, use [`get_coeffs`](@ref). We can iterate over all strings of an operator `o` as follows:
```@example singlestrings
o = rand_local2(2)
for string in o.strings
    println(string)
end
```


## Operations between operators and single strings

A sum of strings is an operator:
```@example singlestrings
o = PauliString("1Z1Y") + PauliString("1Z1Z")
println(typeof(o))
println(o)
```

The main operations, like `+`, `-` `*`, `commutator`, `trace_product` are defined between operators and single strings. For example:
```@example singlestrings
N = 4
o = Operator(4)
o += "X1Y1"
o += "Y1Z1"
s = PauliString("ZZZZ")
println(commutator(s,o))
```

# Translation symmetric Pauli strings
The type `PauliStringTS{Ls}` represents a translation symmetric sum of Pauli strings where the tuple `Ls` specifies the period in each dimension. For example `PauliStringTS{(4,)}` represents a translation symmetric sum of Pauli strings on a 1D lattice of length 4, while `PauliStringTS{(2,2)}` represents a translation symmetric sum of Pauli strings on a 2D lattice of size 2x2.
When printed, only a representative of the equivalence class under translation is shown.
As above, most operations are defined between translation symmetric strings and other translation symmetric strings or translation symmetric operators.
We can construct a `PauliStringTS` out of a `String` or a `PauliString`:
```@example singlestrings
s1 = PauliString("X1Y1")
s1ts = PauliStringTS{(4,)}(s1)
println(s1ts)
println(typeof(s1ts))
```
```@example singlestrings
s2ts = PauliStringTS{(4,)}("XX1X")
println(s2ts)
```
Note a `PauliStringTS` represent a sum of strings, so when computing the commutation of two `PauliStringTS`, the result is not necessarilly a single string and is returned as an `OperatorTS`:
```@example singlestrings
sts1 = PauliStringTS{(4,)}("XX11")
sts2 = PauliStringTS{(4,)}("ZZ11")
o = commutator(sts1, sts2)
println(o)
println(typeof(o))
```
