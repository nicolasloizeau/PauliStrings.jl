# [Translation symmetry in 1D with [`OperatorTS1D`](@ref)] (@id translation)

Here we will show how to take advantage of translation symmetry to save time and memory in PauliStrings.jl.
Consider the 1D Ising Hamiltonian with periodic boundary conditions
$$H=-J(\sum_{i}Z_i Z_{i+1} +g \sum_i X_i)$$.
There is no need to actually store all the Pauli strings in this case. H is fully specified by just $$-JZ_1Z_2$$ and $$-Jg X_1$$ and the fact that it's translation symmetric.

In general, a 1D translation symmetric operator can be written as $$\sum_i T_i H_0$$ where $$T_i$$ is the i-sites translation operator and $$H_0$$ is the local operator that generates $$H$$. $$H_0$$ can be chosen so that it's only composed of Pauli strings that start on the first site.

In PauliStrings.jl, the structure [`OperatorTS1D`](@ref) lets you manipulate operators in this format. If your problem is 1D translation symmetric, [`OperatorTS1D`](@ref) will be much faster than [`Operator`](@ref).

## Construction

There are two ways of constructing an [`OperatorTS1D`](@ref) :
### From the full translation symmetric $$H$$
First we construct the full [`Operator`](@ref):
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
H = ising1D(N, J, g)
```
then convert it to an [`OperatorTS1D`](@ref)
```julia
Hts = OperatorTS1D(H)
```
If `H` is not translation symmetric, then an error will be returned
### From the local generator $$H_0$$
Construct $$H_0$$ using [`Operator`](@ref):
```julia
H = Operator(N)
H += -J, "Z",1,"Z",2
H += -J*g,"X",1
```
then convert it to an [`OperatorTS1D`](@ref)
```julia
Hts = OperatorTS1D(H, full=false)
```
Note that here, we need to set `full=false` in order to specify that we are not passing the full Hamiltonian but just its local generator.

## Manipulation
All the operations defined on [`Operator`](@ref) are also defined on [`OperatorTS1D`](@ref).

Construct a simple translation invariant operator on 4 sites:
```julia
H = Operator(N)
H += "X", 1
H += "Z", 1,"Z", 2

Hts = OperatorTS1D(H, full=false)

println(H)
println(Operator(Hts))
```
```
(1.0 + 0.0im) X111
(1.0 + 0.0im) ZZ11

(1.0 + 0.0im) 1ZZ1
(1.0 + 0.0im) 1X11
(1.0 + 0.0im) X111
(1.0 + 0.0im) Z11Z
(1.0 + 0.0im) 11ZZ
(1.0 + 0.0im) 11X1
(1.0 + 0.0im) ZZ11
(1.0 + 0.0im) 111X
```
Note that only the local generator is printed when printing [`OperatorTS1D`](@ref), not the full operator.

Multiplication:
```
julia> Hts*Hts
(1.0 + 0.0im) X1X1
(2.0 + 0.0im) Z1Z1
(2.0 + 0.0im) 1111
(1.0 + 0.0im) ZZZZ
(2.0 + 0.0im) XZZ1
(2.0 + 0.0im) XX11
(2.0 + 0.0im) ZZX1
```

Addition:
```
julia> Hts+Hts
(2.0 + 0.0im) X111
(2.0 + 0.0im) ZZ11
```

etc.

## Example: computing $$Tr(H^k)$$
As an example of performance gains of using [`OperatorTS1D`](@ref) instead of [`Operator`](@ref) we compute the 8th moment ($$Tr(H^8)$$) of a 30 spin system.

Using the function defined above we construct a Ising Hamiltonian and convert it to an [`OperatorTS1D`](@ref):
```
N = 30
k = 8
H = ising1D(N, 1)
Hts = OperatorTS1D(H)
```
then we compute the kth moment ([`trace_product`](@ref)) using both [`Operator`](@ref) and [`OperatorTS1D`](@ref) :
```
julia> @time println(trace_product(H, k))
1.1904927790006272e18 + 0.0im
 80.697013 seconds (28.91 k allocations: 111.213 MiB, 0.07% gc time, 0.04% compilation time)

julia> @time println(trace_product(Hts, k))
1.190492779000627e18 + 0.0im
  1.951678 seconds (37.09 k allocations: 36.165 MiB, 2.00% gc time, 2.01% compilation time)
```
[`OperatorTS1D`](@ref) is 40 times faster in this case.

## Translation symmetry in 2D

Similar to the 1D case, a 2D translation symmetric operator can be constructed by either specifying an operator on each site, or by using a local generator. Here we construct a 2D transverse field ising model on a rectangular lattice with periodic boundary conditions with the help of the [`string_2d`](@ref) function:
```julia
function ising2D(L1, L2, g)
    H = Operator(L1 * L2)
    for x in 1:L1
        for y in 1:L2
            H += string_2d(("Z", x, y, "Z", x + 1, y), L1, L2, pbc=true) # horizontal
            H += string_2d(("Z", x, y, "Z", x, y + 1), L1, L2, pbc=true) # vertical
            H += g * string_2d(("X", x, y), L1, L2, pbc=true) # transverse field
        end
    end
    return H
end
L1 = 3
L2 = 2
```
```
julia> H = ising2D(L1, L2, 0.5)
(1.0 + 0.0im) Z1Z111
(1.0 + 0.0im) 111Z1Z
(1.0 + 0.0im) 1ZZ111
(1.0 + 0.0im) 1111ZZ
(0.5 + 0.0im) 1X1111
(2.0 + 0.0im) 11Z11Z
(0.5 + 0.0im) X11111
(2.0 + 0.0im) Z11Z11
(2.0 + 0.0im) 1Z11Z1
(1.0 + 0.0im) 111ZZ1
(0.5 + 0.0im) 1111X1
(0.5 + 0.0im) 11111X
(0.5 + 0.0im) 11X111
(1.0 + 0.0im) ZZ1111
(0.5 + 0.0im) 111X11
```
In general, if you have a lattice with extents $L_1, L_2$ in the $a_1$ and $a_2$ directions respectively, a `PauliString` is written in column-major order (similar to how a matrix is flattened in Julia), that is, $L_2$ concatenated chunks of length $L_1$ each.

To convert to an [`OperatorTS2D`](@ref), we have to specify the extent of the lattice in the $a_1$ direction.
```
julia> Hts = OperatorTS2D(H, L1)
(0.5 + 0.0im)   X11
                111
(1.0 + 0.0im)   Z11
                Z11
(1.0 + 0.0im)   ZZ1
                111
```
Alternatively, you can also specify the local generator only:
```julia
H0 = Operator(L1 * L2)
H0 += 1.0 * string_2d(("Z", 1, 1, "Z", 2, 1), L1, L2)
H0 += 1.0 * string_2d(("Z", 1, 1, "Z", 1, 2), L1, L2)
H0 += 0.5 * string_2d(("X", 1, 1), L1, L2)
```
```
julia> Hts = OperatorTS2D(H0, L1, full=false)
(0.5 + 0.0im)   X11
                111
(1.0 + 0.0im)   Z11
                Z11
(1.0 + 0.0im)   ZZ1
                111
```
