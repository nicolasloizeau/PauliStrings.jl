# Constructing operators

Start by importing PauliStrings:
```@example constructing
using PauliStrings
import PauliStrings as ps
```

To construct an operator we first need to declare an empty operator of $N$ qubits:
```julia
H = Operator(N)
```
For the moment, `PauliStrings.jl` supports a maximum of 64 qubits.

We can add a term of the form $J X_i$ by doing
```julia
H += J, "X", i
```

and a term of the form $J X_iX_j$ by doing
```julia
H += J, "X", i, "X", j
```

Similarly, we add a term of the form $J X_iX_jX_k$ by doing
```julia
H += J, "X", i, "X", j, "X", k
```
etc.

## 1D transverse Ising model
Let's construct the Hamiltonian of a [1D transverse Ising model](https://en.wikipedia.org/wiki/Transverse-field_Ising_model)
$$H=-J(\sum_{<i,j>}Z_i Z_j +g \sum_i X_i)$$

```@example constructing
function ising1D(N, J, g)
    H = Operator(N)
    for j in 1:(N - 1)
        H += "Z",j,"Z",j+1
    end
    H += "Z",1,"Z",N # periodic boundary condition
    for j in 1:N
        H += g,"X",j
    end
    return -J*H
end
```
Note that the first qubit starts at index 1, following Julia's 1-based index.

Operators can be printed in strings format with the `println` function:
```@example constructing
println(ising1D(3, 1, 0.5))
```


## 2D transverse Ising model


Here we construct a 2D Ising model on a square lattice of L*L sites, with no periodic boundary conditions.


```julia
function ising2D(L, J, g)
    H = ps.Operator(L*L)
    for x in 1:L-1
        for y in 1:L
            # convert x,y to qubit index
            i = L*(y-1)+x
            j = L*(y-1)+(x+1)
            # horizontal interaction terms
            H += ('Z',i,'Z',j)
            # convert x,y to qubit index
            i = L*(x-1)+y
            j = L*x+y
            # vertical interaction terms
            H += ('Z',i,'Z',j)
        end
    end
    for j in 1:L*L
        H += g,"X",j
    end
    return -J*H
end
```
