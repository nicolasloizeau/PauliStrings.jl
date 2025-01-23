
"""
    all_strings(N::Int)

Return the sum of all the strings supported on N spins, with coeficients 1
"""
function all_strings(N::Int)
    O = Operator(N)
    for i in 0:2^N-1
        for j in 0:2^N-1
            push!(O.v, i)
            push!(O.w, j)
            push!(O.coef, (1im)^ycount(i, j))
        end
    end
    return O
end


"""
    all_k_local(N::Int, k::Int)

Return the sum of all the k-local strings supported on N spins, with coeficients 1.
These are k-local only strings, not including strings shorter than k.

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
            if pauli_weight(i, j) == k
                push!(O.v, i)
                push!(O.w, j)
                push!(O.coef, (1im)^ycount(i, j))
            end
        end
    end
    return O
end


"""
    all_z(N::Int)

Return the sum of all the strings supported on N spins with only z and with coeficients 1

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
        push!(O.v, i)
        push!(O.w, 0)
        push!(O.coef, 1)
    end
    return O
end

function all_z(N::Int, bits::Vector{Int})
    O = Operator(N)
    mask = sum(1 << (b-1) for b in bits)
    for i in 0:2^N-1
        if (i & ~mask) == 0
            push!(O.v, i)
            push!(O.w, 0)
            push!(O.coef, 1)
        end
    end
    return O
end



"""
    all_x(N::Int)

Return the sum of all the strings supported on N spins with only x and with coeficients 1

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
        push!(O.v, 0)
        push!(O.w, i)
        push!(O.coef, 1)
    end
    return O
end


"""
    all_y(N::Int)

Return the sum of all the strings supported on N spins with only y and with coeficients 1

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
        push!(O.v, i)
        push!(O.w, i)
        push!(O.coef, (1im)^ycount(i, i))
    end
    return O
end

"""
    Base.push!(o::Operator, c::Number, v::Unsigned, w::Unsigned)

Add string c,(v,w) to operator o. Note that c is not the coeficient in front of the pauli string
but the coeficient in front of the real string. A factor (1im)^ycount(v,w) relate the two.
"""
function Base.push!(o::Operator, c::Number, v::Unsigned, w::Unsigned)
    push!(o.v, v)
    push!(o.w, w)
    push!(o.coef, c)
end
