
[![Build Status](https://github.com/nicolasloizeau/PauliStrings.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/nicolasloizeau/PauliStrings.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://nicolasloizeau.github.io/PauliStrings.jl/dev)

# PauliStrings.jl
PauliStrings.jl is a Julia package for many-body quantum mechanics with Pauli string represented as binary integers (as in [https://journals.aps.org/pra/abstract/10.1103/PhysRevA.68.042318](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.68.042318)).

## [Documentation](https://nicolasloizeau.github.io/PauliStrings.jl/dev/)

## Installation
`using Pkg; Pkg.add(url="https://github.com/nicolasloizeau/PauliStrings.jl")` or `] add https://github.com/nicolasloizeau/PauliStrings.jl`

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

Add a 1-qubit string:
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


## Print and export
`print` shows a list of terms with coeficients e.g :
```julia
julia> println(H)
(10.0 - 0.0im) 1ZZ
(5.0 - 0.0im) 1Z1
(15.0 + 0.0im) XYZ
(5.0 + 0.0im) 1YY
```

Export a list of strings with coeficients:
```julia
coefs, strings = ps.op_to_strings(H)
```

## Truncate, Cutoff, Trim, Noise
`ps.truncate(H,M)` removes Pauli strings longer than M (returns a new Operator)
`ps.cutoff(H,c)` removes Pauli strings with coeficient smaller than c in absolute value (returns a new Operator)
`ps.trim(H,N)` keeps the first N trings with higest weight (returns a new Operator)
`ps.prune(H,alpha)` keeps terms with probability 1-exp(-alpha*abs(c)) (returns a new Operator)

`ps.add_noise(H,g)` adds depolarizing noise that make each strings decay like $e^{gw}$ where $w$ is the lenght of the string. This is usefull when used with `trim` to keep the number of strings manageable during time evolution.


## Time evolution

`ps.rk4(H, O, dt; hbar=1, heisenberg=false)` performs a step of Runge Kutta and returns the new updated O(t+dt)

H can be an Operator, or a function that takes a time and return an Operator. In case H is a function, a time also needs to be passed to `rk4(H, O, dt, t)`. O is an Observable or a density matrix to time evolve.
If evolving an observable in the heisenberg picture, set `heisenberg=true`.

An example is in `time_evolve_example.jl`.
The following will time evolve O in the Heisenberg picture. At each step, we add depolarizing noise and trim the operator to keep the number of strings manageable
```julia
function evolve(H, O, M, times, noise)
    dt = times[2]-times[1]
    for t in times
        O = ps.rk4(H, O, dt; heisenberg=true, M=M) #preform one step of rk4, keep only M strings
        O = ps.add_noise(O, noise*dt) #add depolarizingn noise
        O = ps.trim(O, M) # keep the M strings with the largest weight
    end
    return O
end
```

Time evolution of the spin correlation function $\textup{Tr}(Z_1(0)Z_1(t))$ in the chaotic spin chain.
Check time_evolve_example.jl to reproduce the plot.
<!-- ![plot](./time_evolve_example.png) -->

## Lanczos
Compute lanczos coeficients
```julia
bs = ps.lanczos(H, O, steps, nterms)
```
`H` : Hamiltonian

`O` : starting operator

`nterms` : maximum number of terms in the operator. Used by trim at every step

Results for X in XX from https://journals.aps.org/prx/pdf/10.1103/PhysRevX.9.041017 :

<!-- ![plot](./lanczos_example.png) -->
