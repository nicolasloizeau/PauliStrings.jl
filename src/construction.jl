
"""
    all_strings(N::Int)

Return the sum of all the strings supported on N spins, with coeficients 1
"""
function all_strings(N::Int)
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
    all_k_local(N::Int, k::Int)

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
function all_k_local(N::Int, k::Int)
    O = Operator(N)
    for i in 0:2^N-1
        for j in 0:2^N-1
            p = paulistringtype(O)(i, j)
            if pauli_weight(p) == k
                push!(O.strings, p)
                push!(O.coeffs, (1im)^ycount(p))
            end
        end
    end
    return O
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
