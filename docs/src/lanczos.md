# Lanczos


Here we show how to use PauliStrings.jl to run the recursion method and compute [`lanczos`](@ref) coefficients.
We will focus on reproducing the X in XX results from figure 2 of [Parker 2019](https://journals.aps.org/prx/abstract/10.1103/PhysRevX.9.041017).

Start by importing PauliStrings:
```julia
using PauliStrings
import PauliStrings as ps
```

Define the XX Hamiltonian $$H = \sum_i X_iX_{i+1}+Y_iY_{i+1}$$ on a 1D chain with periodic boundary conditions
```julia
function XX(N)
    H = ps.Operator(N)
    for j in 1:(N - 1)
        H += "X",j,"X",j+1
        H += "Z",j,"Z",j+1
    end
    H += "X",1,"X",N
    H += "Z",1,"Z",N
    return H
end
```
and the X operator $$\sum_{i=1}^N X_i$$
```julia
function X(N)
    H = ps.Operator(N)
    for j in 1:N
        H += "X",j
    end
    return H
end
```

Initialize a Hamiltonian and an operator:
```julia
N = 50 # system size
H = XX(N) # Hamiltonian
O = X(N) # operator
```

Compute and plot the [`lanczos`](@ref) coefficients for different truncations.
For each level of truncation, we keep only $$2^p$$ terms with the highest weight at each step of the [`lanczos`](@ref) algorithm.
```julia
ioff() # pyplot
# nterms is the max pauli string length
for p in (14,16,18,20)
    @time bs = ps.lanczos(H, O, 20, 2^p; keepnorm=true)
    plot(bs, label="trim: 2^$p")
end
legend()
ylabel(L"$b_n$")
xlabel(L"$n$")
title("X in XX, N=$N spins")
savefig("lanczos_example.png")
show()
```
![plot](./lanczos_example.png)


# Tacking advantage of translation symmetry
If the problem is 1D and has translation symmetry, we can take advantage of the symmetry to save time and memory. Just pass $$O$$ and $$H$$ as [`OperatorTS1D`](@ref) to  [`lanczos`](@ref). For example :
```julia
Hts = OperatorTS1D(H)
Ots = OperatorTS1D(O)
bs = ps.lanczos(Hts, Ots, 20, 2^p; keepnorm=true)
```
Check the [translation symmetry tutorial](@ref translation).
