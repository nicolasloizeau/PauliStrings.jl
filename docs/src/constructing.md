# Constructing operators

To construct an operator we first need to declare an empty operator of $N$ qubits :
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

## 1D transverse ising model
Lets construct the Hamiltonian of a [1D transverse ising model](https://en.wikipedia.org/wiki/Transverse-field_Ising_model)
$$H=-J(\sum_{<i,j>}Z_i Z_j +g \sum_i X_i)$$

```julia
function ising1D(N, J, g)
    H = Operator(N)
    for j in 1:(N - 1)
        H += "Z",j,"Z",j+1
    end
    H += "Z",1,"Z",N #periodic boundary condition
    for j in 1:N
        H += g,"X",j
    end
    return -J*H
end
```
Note that the first qubit start at index 1, following Julia's 1-based index.

## 2D transverse ising model
