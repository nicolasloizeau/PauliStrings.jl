

# Getting started

PauliStrings.jl is a Julia package for many-body quantum mechanics with Pauli string represented as binary integers (as in [https://journals.aps.org/pra/abstract/10.1103/PhysRevA.68.042318](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.68.042318)).

## Initializing an operator

Import the library and initialize a operator of 4 qubits
```julia
using PauliStrings
import PauliStrings as ps
H = ps.Operator(4)
```

Add a Pauli strings to the operator

```julia
H += "XYZ1"
H += "1YZY"
```

```
julia> H
(1.0 - 0.0im) XYZ1
(1.0 - 0.0im) 1YZY
```

Add a Pauli string with a coeficient
```julia
H += -1.2,"XXXZ" #coeficient can be complex
```

Add a 2-qubit string coupling qubits i and j with X and Y:
```julia
H += 2, "X", i, "Y", j # with a coeficient=2
H += "X", i, "Y", j # with a coeficient=1
```

Add a 1-qubit string on site 1
```julia
H += 2, "Z", i # with a coeficient=2
H += "Z", i # with a coeficient=1
H += "S+", i
```

Supported sites operators are `X`, `Y`, `Z`, `Sx`$=X/2$, `Sy`$=Y/2$, `Sz`$=Z/2$, `S+`$=(X+iY)/2$, `S-`$=(X-iY)/2$.

## Basic Algebra
The Operator type supports the +,-,* operators with other Operators and Numbers:
```julia
H3 = H1*H2
H3 = H1+H2
H3 = H1-H2
H3 = H1+2 # adding a scalar is equivalent to adding the unit times the scalar
H = 5*H # multiply operator by a scalar
```
Trace : `ps.trace(H)`

Frobenius norm : `ps.opnorm(H)`

Conjugate transpose : `ps.dagger(H)`

Number of terms: `length(H)`

Commutator: `ps.com(H1, H2)`. This is much faster than `H1*H2-H2*H1`
