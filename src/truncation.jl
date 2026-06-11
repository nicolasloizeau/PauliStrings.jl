
"""
    truncate(o::Operator, max_lenght::Int; keepnorm::Bool = false)

Remove all terms of length > max_lenght.
Keep all terms of length <= max_lenght.
i.e remove all M-local terms with M>max_lenght

# Example
```julia
A = Operator(4)
A += "X",1,"X",2
A += "Z",1,"Z",2,"Z",4
```
```
julia> A
(1.0 + 0.0im) ZZ1Z
(1.0 + 0.0im) XX11

julia> ps.truncate(A,2)
(1.0 + 0.0im) XX11
```
"""
function PauliStrings.truncate(o::AbstractOperator, max_lenght::Int; keepnorm::Bool=false)
    ks = paulistringtype(o)[]
    vs = scalartype(o)[]
    for (p, c) in pairs(o)
        if pauli_weight(p) <= max_lenght
            push!(ks, p)
            push!(vs, c)
        end
    end
    o2 = typeof(o)(ks, vs)
    if keepnorm
        return o2 * norm(o) / norm(o2)
    end
    return o2
end



"""
    k_local_part(o::Operator, k::Int; atmost=false)

Return the k-local part of o. I.e all the strings of lenght k.
Set `atmost=true` to return the 'at most' k-local part, i.e all the strings of length <= k.

# Example
```julia
A = Operator(4)
A += "X",1,"X",2
A += "Z",1,"Z",2,"Z",4
A += "1X11"
```
```
julia> A
(1.0 + 0.0im) ZZ1Z
(1.0 + 0.0im) 1X11
(1.0 + 0.0im) XX11

julia> k_local_part(A,2)
(1.0 + 0.0im) XX11
```
"""
function k_local_part(o::AbstractOperator, k::Int; atmost=false)
    ks = paulistringtype(o)[]
    vs = scalartype(o)[]
    for (p, c) in pairs(o)
        l = pauli_weight(p)
        if l == k || (atmost && l <= k)
            push!(ks, p)
            push!(vs, c)
        end
    end
    return typeof(o)(ks, vs)
end




"""
    trim(o::Operator, max_strings::Int; keepnorm::Bool = false, keep::Operator=Operator(N))

Keep the first `max_strings` terms with largest coeficients.

`keepnorm` is set to true to keep the norm of o.

`keep` is an operator that specify a set of strings that cannot be removed

# Example
```julia
A = Operator(4)
A += 1,"XXXX"
A += 2,"XX11"
A += 3,"XX1X"
A += 4,"ZZXX"
B = Operator(4)
B += 1,"XX11"
B += 1,"XX1X"
```
```
julia> trim(A,2)
(4.0 + 0.0im) ZZXX
(3.0 + 0.0im) XX1X

julia> trim(A,2;keep=B)
(4.0 + 0.0im) ZZXX
(3.0 + 0.0im) XX1X
(2.0 + 0.0im) XX11
```
"""
function trim(o::AbstractOperator, max_strings::Int; keepnorm::Bool=false, keep::AbstractOperator=Operator(0))
    if length(o) <= max_strings
        return deepcopy(o)
    end
    # position-aligned strings/coeffs
    ks = collect(keys(o))
    vs = collect(values(o))
    # keep the N first indices
    i = collect(partialsortperm(abs.(vs), 1:max_strings; rev=true))
    # add the string to keep in case there was a specified string to keep
    if length(keep) > 0
        for p_keep in keys(keep) # for each string in the keep operator
            # we check if it is in o and has been removed
            j = findfirst(==(p_keep), ks)
            if !isnothing(j) && !(j in i)
                push!(i, j)
            end
        end
    end

    o1 = typeof(o)(ks[i], vs[i])

    if keepnorm
        return o1 * norm(o) / norm(o1)
    end
    return o1
end

"""
     prune(o::Operator, alpha::Real; keepnorm::Bool = false)

Keep terms with probability 1-exp(-alpha*abs(c)) where c is the weight of the term
"""
function prune(o::AbstractOperator, alpha::Real; keepnorm::Bool=false)
    ks = paulistringtype(o)[]
    vs = scalartype(o)[]
    for (s, c) in pairs(o)
        p = 1 - exp(-alpha * abs(c))
        if rand() < p
            push!(ks, s)
            push!(vs, c)
        end
    end
    o1 = typeof(o)(ks, vs)
    if keepnorm
        return o1 * norm(o) / norm(o1)
    end
    return o1
end

"""
    cutoff(o::Operator, epsilon::Real; keepnorm::Bool = false)

Remove all terms with weight < epsilon

# Example
```julia
A = Operator(4)
A += 1,"XXXX"
A += 2,"XX11"
A += 3,"XX1X"
A += 4,"ZZXX"
```
```
julia> cutoff(A, 2.5)
(3.0 + 0.0im) XX1X
(4.0 + 0.0im) ZZXX
```
"""
function cutoff(o::AbstractOperator, epsilon::Real; keepnorm::Bool=false)
    ks = paulistringtype(o)[]
    vs = scalartype(o)[]
    for (p, c) in pairs(o)
        if abs(c) > epsilon
            push!(ks, p)
            push!(vs, c)
        end
    end
    o2 = typeof(o)(ks, vs)
    if keepnorm
        return o2 * norm(o) / norm(o2)
    end
    return o2
end



function participation(o::Operator)
    return sum(o.coeffs .^ 4) / sum(o.coeffs .^ 2)^2
end


"""
    zpart(o::AbstractOperator)

Keep strings with only Z and identity
"""
zpart(o::AbstractOperator) = diag(o)


"""
    xpart(o::AbstractOperator)

Keep strings with only X and identity
"""
function xpart(o::AbstractOperator)
    ks = paulistringtype(o)[]
    vs = scalartype(o)[]
    for (p, c) in pairs(o)
        if zcount(p) == 0 && ycount(p) == 0
            push!(ks, p)
            push!(vs, c)
        end
    end
    return typeof(o)(ks, vs)
end


"""
    ypart(o::AbstractOperator)

Keep strings with only Y and identity
"""
function ypart(o::AbstractOperator)
    ks = paulistringtype(o)[]
    vs = scalartype(o)[]
    for (p, c) in pairs(o)
        if zcount(p) == 0 && xcount(p) == 0
            push!(ks, p)
            push!(vs, c)
        end
    end
    return typeof(o)(ks, vs)
end
