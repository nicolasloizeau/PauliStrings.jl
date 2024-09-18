# Lanczos


Here we show how to use PauliStrings.jl to run the recursion method and compute [`lanczos`](@ref) coefficients.
We will focus on reproducing the X in XX results from figure 2 of [Parker 2019](https://journals.aps.org/prx/abstract/10.1103/PhysRevX.9.041017).

Start by importing PauliStrings :
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
and the X operator $\sum_{i=1}^N X_i$
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
H = XX(N) #hamiltonian
O = X(N) #operator
```

Compute and plot the [`lanczos`](@ref) coefficients for different truncation.
For each level of truncation, we keep only 2^p with the highest weight at each step of the [`lanczos`](@ref) algorithm.
```julia
ioff()#pyplot
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
If the problem is 1D and has translation symmetry, we can take advantage of the symmetry to save memory. Instead of storing the full operator at each step of the [`lanczos`](@ref) algorithm, we only store Pauli strings starting on the first site. In order for this to be possible, at each step, we need to shift all the strings so they start on the first site. This is done using the [`shift_left`](@ref) function :
```julia
A = Operator(4)
A += "XYZ1"
A += "11ZZ"
A += "1XZY"
A += "ZZ11"
```
```julia
julia> shift_left(A)
(1.0 - 0.0im) XZY1
(1.0 - 0.0im) XYZ1
(2.0 + 0.0im) ZZ11
```
Note that the [`shift_left`](@ref) function accumulates identical terms. Here it accumulated "11ZZ" and "ZZ11".
This idea is already implemented in the [`lanczos`](@ref) function. The user just needs to set `localop=true` when calling [`lanczos`](@ref). Now we can define the initial X operator just on a single site:
```julia
function X(N)
    H = ps.Operator(N)
    H += "X",1
    return H
end
```
Note that we can also still use the full operator defined on every site. If so, then the `lanczos` function will shift if left for us before starting iterating.  
