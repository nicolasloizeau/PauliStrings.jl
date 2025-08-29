# [Translation symmetry with [`OperatorTS`](@ref)] (@id translation)

Here we will show how to take advantage of translation symmetry to save time and memory in PauliStrings.jl.
Consider the 1D Ising Hamiltonian with periodic boundary conditions
$$H=-J(\sum_{i}Z_i Z_{i+1} +g \sum_i X_i)$$.
There is no need to actually store all the Pauli strings in this case. H is fully specified by just $$-JZ_1Z_2$$ and $$-Jg X_1$$ and the fact that it's translation symmetric.

In general, a 1D translation symmetric operator can be written as $$\sum_i T_i H_0$$ where $$T_i$$ is the i-sites translation operator and $$H_0$$ is the local operator that generates $$H$$. $$H_0$$ can be chosen so that it's only composed of Pauli strings that start on the first site.

In PauliStrings.jl, the structure [`OperatorTS`](@ref) lets you manipulate operators in this format. If your problem is translation symmetric, [`OperatorTS`](@ref) will be much faster than [`Operator`](@ref).

## Construction

The constructor of [`OperatorTS`](@ref) takes a normal operator and turns it into a lazy sum over all possible translations. There are two ways of using this constructor.
### From the full translation symmetric $$H$$
First we construct the full [`Operator`](@ref):
```jldoctest 1d; output=false, setup=:(using PauliStrings)
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
N = 4
J = -1.0
g = 1.0

H = ising1D(N, J, g)

# output

(1.0 + 0.0im) X111
(1.0 + 0.0im) 1X11
(1.0 + 0.0im) 11X1
(1.0 + 0.0im) 111X
(1.0 + 0.0im) ZZ11
(1.0 + 0.0im) 1ZZ1
(1.0 + 0.0im) Z11Z
(1.0 + 0.0im) 11ZZ
```
then convert it to an [`OperatorTS`](@ref)
```jldoctest 1d
Hts = OperatorTS{(N,)}(H)/N

# output

(1.0 + 0.0im) 111X
(1.0 + 0.0im) 11ZZ
```
Note the normalization factor `1/N` that normalizes the symmetric sum. We can use [`is_ts`](@ref) to check if `H` is actually symmetric.
### From the local generator $$H_0$$
Construct $$H_0$$ using [`Operator`](@ref):
```jldoctest 1d; output=false
H = Operator(N)
H += -J, "Z",1,"Z",2
H += -J*g,"X",1

# output

(1.0 + 0.0im) X111
(1.0 + 0.0im) ZZ11
```
then convert it to an [`OperatorTS`](@ref)
```jldoctest 1d
Hts = OperatorTS{(N,)}(H)

# output

(1.0 + 0.0im) 111X
(1.0 + 0.0im) 11ZZ
```
In this case, no normalization factor is needed.
## Manipulation
All the operations defined on [`Operator`](@ref) are also defined on [`OperatorTS`](@ref).

Construct a simple translation invariant operator on 4 sites:
```jldoctest 1d
H = Operator(N)
H += "X", 1
H += "Z", 1,"Z", 2

Hts = OperatorTS{(N,)}(H)

println(H)
println(resum(Hts))

# output

(1.0 + 0.0im) X111
(1.0 + 0.0im) ZZ11

(1.0 + 0.0im) X111
(1.0 + 0.0im) 1X11
(1.0 + 0.0im) 11X1
(1.0 + 0.0im) 111X
(1.0 + 0.0im) ZZ11
(1.0 + 0.0im) 1ZZ1
(1.0 + 0.0im) Z11Z
(1.0 + 0.0im) 11ZZ
```
Note that only the local generator is printed when printing [`OperatorTS`](@ref), not the full operator.

Multiplication:
```jldoctest 1d
julia> Hts*Hts
(1.0 + 0.0im) 1X1X
(1.0 + 0.0im) ZZZZ
(2.0 + 0.0im) 1111
(2.0 + 0.0im) 11XX
(2.0 + 0.0im) 1Z1Z
(2.0 + 0.0im) X1ZZ
(2.0 + 0.0im) 1XZZ
```

Addition:
```jldoctest 1d
julia> Hts+Hts
(2.0 + 0.0im) 111X
(2.0 + 0.0im) 11ZZ
```

etc.

## Example: computing $$Tr(H^k)$$
As an example of performance gains of using [`OperatorTS`](@ref) instead of [`Operator`](@ref) we compute the 8th moment ($$Tr(H^8)$$) of a 30 spin system.

Using the function defined above we construct a Ising Hamiltonian and convert it to an [`OperatorTS`](@ref):
```jldoctest 1d;output=false
N = 30
k = 8
H = ising1D(N, J, g)
Hts = OperatorTS{(N,)}(H)/N

# output

(1.0 + 0.0im) 11111111111111111111111111111X
(1.0 + 0.0im) 1111111111111111111111111111ZZ
```
then we compute the kth moment ([`trace_product`](@ref)) using both [`Operator`](@ref) and [`OperatorTS`](@ref) :
```julia
julia> @time trace_product(H, k)
  0.688497 seconds (451 allocations: 135.925 MiB, 54.67% gc time)
1.1904927790006272e18 + 0.0im
julia> @time trace_product(Hts, k)
  0.038475 seconds (235 allocations: 7.023 MiB, 27.41% gc time)
1.1904927790006272e18 + 0.0im
```
[`OperatorTS`](@ref) is 18 times faster in this case.

## Translation symmetry in 2D

Similar to the 1D case, a 2D translation symmetric operator can be constructed by either specifying an operator on each site, or by using a local generator. Here we construct a 2D transverse field ising model on a rectangular lattice with periodic boundary conditions with the help of the [`string_2d`](@ref) function:
```jldoctest 2d; setup=:(using PauliStrings), output = false
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

# output

2
```
```jldoctest 2d
julia> H = ising2D(L1, L2, 0.5)
(0.5 + 0.0im) X11111
(0.5 + 0.0im) 1X1111
(0.5 + 0.0im) 11X111
(0.5 + 0.0im) 111X11
(0.5 + 0.0im) 1111X1
(0.5 + 0.0im) 11111X
(1.0 + 0.0im) ZZ1111
(1.0 + 0.0im) Z1Z111
(1.0 + 0.0im) 1ZZ111
(1.0 + 0.0im) 111ZZ1
(1.0 + 0.0im) 111Z1Z
(1.0 + 0.0im) 1111ZZ
(2.0 + 0.0im) Z11Z11
(2.0 + 0.0im) 1Z11Z1
(2.0 + 0.0im) 11Z11Z
```
In general, if you have a lattice with extents $L_1, L_2$ in the $a_1$ and $a_2$ directions respectively, a `PauliString` is written in column-major order (similar to how a matrix is flattened in Julia), that is, $L_2$ concatenated chunks of length $L_1$ each.

To convert to an [`OperatorTS`](@ref), in higher dimensions we have to specify the size of the lattice.
```jldoctest 2d
julia> Hts = OperatorTS{(L1, L2)}(H)/(L1*L2)
(0.5 + 0.0im) 111
              11X
(1.0 + 0.0im) 11Z
              11Z
(1.0 + 0.0im) 111
              1ZZ
```
Alternatively, you can also specify the local generator only:
```jldoctest 2d; output=false
H0 = Operator(L1 * L2)
H0 += 1.0 * string_2d(("Z", 1, 1, "Z", 2, 1), L1, L2)
H0 += 1.0 * string_2d(("Z", 1, 1, "Z", 1, 2), L1, L2)
H0 += 0.5 * string_2d(("X", 1, 1), L1, L2)

# output

(0.5 + 0.0im) X11111
(1.0 + 0.0im) ZZ1111
(1.0 + 0.0im) Z11Z11
```
```jldoctest 2d
julia> Hts = OperatorTS{(L1,L2)}(H0)
(0.5 + 0.0im) 111
              11X
(1.0 + 0.0im) 11Z
              11Z
(1.0 + 0.0im) 111
              1ZZ
```
