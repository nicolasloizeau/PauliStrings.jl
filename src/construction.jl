

"""
    all_strings(N::Int)

Return the sum of all the strings supported on N spins, with coeficients 1
"""
function all_strings(N::Int)
    @warn "`all_strings` is deprecated — use `complete_basis`"
    O = Operator(N)
    for i in 0:2^N-1
        for j in 0:2^N-1
            p = paulistringtype(O)(i, j)
            push!(O.strings, p)
            push!(O.coeffs, (1im)^ycount(p))
        end
    end
    return O
end


"""
    complete_basis(N::Int)

Return a vector of all the strings supported on N spins.
"""
function complete_basis(N::Int)
    strings = PauliString{N}[]
    for i in 0:2^N-1
        for j in 0:2^N-1
            push!(strings, PauliString{N}(i, j))
        end
    end
    return strings
end



"""
    all_k_local(N::Int, k::Int; atmost=false)

Return the sum of all the k-local strings supported on N spins, with coeficients 1.
These are k-local only strings, not including strings shorter than k.

# Example
```
julia> all_k_local(2, 1)
(1.0 + 0.0im) X1
(1.0 + 0.0im) 1X
(1.0 + 0.0im) Z1
(1.0 - 0.0im) Y1
(1.0 + 0.0im) 1Z
(1.0 - 0.0im) 1Y
```
"""
function all_k_local(N::Int, k::Int; atmost=false)
    @warn "`all_k_local` is deprecated — use `k_local_basis`"
    if atmost
        types = [0, 1, 2, 3]
    else
        types = [1, 2, 3]
    end
    O = Operator(N)
    for sites in combinations(1:N, k)
        ranges = [types for _ in 1:k]
        for ptypes in Iterators.product(ranges...)
            p = zeros(Int, N)
            p[sites] .= ptypes
            O += string_from_inds(p)
        end
    end
    coeffs = ones(length(O))
    set_coeffs(O, coeffs)
    return O
end


"""
    k_local_basis(N::Int, k::Int; atmost=false)

Return a vector of all the k-local strings supported on N spins.
# Example
```
julia> k_local_basis(3,2)
27-element Vector{PauliString}:
 XX1
 YX1
 ZX1
 XY1
 ⋮
 1ZY
 1XZ
 1YZ
 1ZZ
```
"""
function k_local_basis(N::Int, k::Int; atmost=false)
    if atmost
        types = [0, 1, 2, 3]
    else
        types = [1, 2, 3]
    end
    strings = PauliString[]
    for sites in combinations(1:N, k)
        ranges = [types for _ in 1:k]
        for ptypes in Iterators.product(ranges...)
            p = zeros(Int, N)
            p[sites] .= ptypes
            push!(strings, PauliString(string_from_inds(p)))
        end
    end
    return strings
end


"""
    z_basis(N::Int)

Return a vector of all the strings supported on N spins with only z and identity.
# Example
```
julia> z_basis(2)
4-element Vector{PauliString}:
 11
 Z1
 1Z
 ZZ
```
"""
function z_basis(N::Int)
    strings = PauliString[]
    for i in 0:2^N-1
        p = PauliString{N}(i, 0)
        push!(strings, p)
    end
    return strings
end

"""
    x_basis(N::Int)

Return a vector of all the strings supported on N spins with only x and identity.
# Example
```
julia> x_basis(2)
4-element Vector{PauliString}:
 11
 X1
 1X
 XX
```
"""
function x_basis(N::Int)
    strings = PauliString[]
    for i in 0:2^N-1
        p = PauliString{N}(0, i)
        push!(strings, p)
    end
    return strings
end

"""
    y_basis(N::Int)

Return a vector of all the strings supported on N spins with only y and identity.
# Example
```
julia> y_basis(2)
4-element Vector{PauliString}:
 11
 Y1
 1Y
 YY
```
"""
function y_basis(N::Int)
    strings = PauliString[]
    for i in 0:2^N-1
        p = PauliString{N}(i, i)
        push!(strings, p)
    end
    return strings
end



"""
    all_z(N::Int)

Return the sum of all the strings supported on N spins with only z and with coeficients 1

# Example
```
julia> all_z(2)
(1.0 + 0.0im) 11
(1.0 + 0.0im) Z1
(1.0 + 0.0im) 1Z
(1.0 + 0.0im) ZZ
```
"""
function all_z(N::Int)
    @warn "`all_z(N)` is deprecated — use `z_basis(N)`"
    O = Operator(N)
    for i in 0:2^N-1
        p = paulistringtype(O)(i, 0)
        push!(O.strings, p)
        push!(O.coeffs, 1)
    end
    return O
end

function all_z(N::Int, bits::Vector{Int})
    O = Operator(N)
    mask = sum(1 << (b - 1) for b in bits)
    for i in 0:2^N-1
        if (i & ~mask) == 0
            p = paulistringtype(O)(i, 0)
            push!(O.strings, p)
            push!(O.coeffs, 1)
        end
    end
    return O
end


"""
    all_x(N::Int)

Return the sum of all the strings supported on N spins with only x and with coeficients 1

# Example
```
julia> all_x(2)
(1.0 + 0.0im) 11
(1.0 + 0.0im) X1
(1.0 + 0.0im) 1X
(1.0 + 0.0im) XX
```
"""
function all_x(N::Int)
    @warn "`all_x(N)` is deprecated — use `x_basis(N)`"
    O = Operator(N)
    for i in 0:2^N-1
        p = paulistringtype(O)(0, i)
        push!(O.strings, p)
        push!(O.coeffs, 1)
    end
    return O
end


"""
    all_y(N::Int)

Return the sum of all the strings supported on N spins with only y and with coeficients 1

# Example
```
julia> all_y(2)
(1.0 + 0.0im) 11
(1.0 - 0.0im) Y1
(1.0 - 0.0im) 1Y
(1.0 - 0.0im) YY
```
"""
function all_y(N::Int)
    @warn "`all_y(N)` is deprecated — use `y_basis(N)`"
    O = Operator(N)
    for i in 0:2^N-1
        p = paulistringtype(O)(i, i)
        push!(O.strings, p)
        push!(O.coeffs, (1im)^ycount(p))
    end
    return O
end

"""
    Base.push!(o::Operator, c::Number, v::Unsigned, w::Unsigned)

Add string c,(v,w) to operator o. Note that c is not the coeficient in front of the pauli string
but the coeficient in front of the real string. A factor (1im)^ycount(v,w) relate the two.
"""
function Base.push!(o::Operator, c::Number, v::Unsigned, w::Unsigned)
    p = paulistringtype(o)(v, w)
    push!(o.strings, p)
    push!(o.coeffs, c)
end




"""
    majorana(N::Int, k::Int)

Return the k-th Majorana operator on N spins.
There are 2N Majoranas supported on N spins.
They all anticomute :
```math
\\{\\gamma_i, \\gamma_j\\} = 2\\delta_{ij}
```

# Example
```
julia> majorana(4,1)
(1.0 + 0.0im) X111

julia> majorana(4,2)
(1.0 - 0.0im) Y111

julia> majorana(4,3)
(1.0 + 0.0im) ZX11

julia> majorana(4,4)
(1.0 - 0.0im) ZY11
```
"""
function majorana(N::Int, k::Int)
    @assert k <= 2 * N && k >= 1 "k must be between 1 and 2N"
    k -= 1
    O = Operator(N)
    nz = k ÷ 2
    n1 = N - nz - 1
    if k % 2 == 0
        s = "Z"^nz * "X" * "1"^n1
    else
        s = "Z"^nz * "Y" * "1"^n1
    end
    add_string(O, s, 1)
    return O
end

"""
    string_2d(args::Tuple{Vararg{Any}}, L1::Int, L2::Int; pbc=false)

Helper functions to construct 2d pauli strings. The tuple is composed of triplets the form `(P,x,y,...)` where P is one of
"X", "Y", "Z", "Sx", "Sy", "Sz", "S-", "S+", and `x, y` label the position in the lattice. If `pbc = true`, the x and y
coordinates will always be brough back to the range \$[1, L_1]\$ and \$[1, L_2]\$ respectively.


Example:
```julia
L1 = L2 = 2
A = Operator(L1 * L2)
A += 0.5 * string_2d(("Z", 1, 1, "Z", 2, 1), L1, L2) # Horizontal interaction
A += 1.0 * string_2d(("Z", 1, 1, "Z", 1, 2), L1, L2) # Vertical interaction
```
```
julia> A
(1.0 + 0.0im) Z1Z1
(0.5 + 0.0im) ZZ11
```
"""
function string_2d(args::Tuple{Vararg{Any}}, L1::Int, L2::Int; pbc=false)
    o = one(Operator(L1 * L2))
    for t in 1:3:length(args)
        o2 = zero(o)
        P = args[t]::String
        x = args[t+1]::Int
        y = args[t+2]::Int
        if pbc
            idx = mod1(x, L1) + L1 * (mod1(y, L2) - 1)
        else
            @assert (1 <= x && x <= L1) && (1 <= y && y <= L2)
            idx = x + L1 * (y - 1)
        end
        o2 += P, idx
        o *= o2
    end
    return o
end
