
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
    o2 = typeof(o)()
    for i in 1:length(o)
        p = o.strings[i]
        if pauli_weight(p) <= max_lenght
            push!(o2.coeffs, o.coeffs[i])
            push!(o2.strings, p)
        end
    end
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
    o2 = typeof(o)()
    for i in 1:length(o)
        p = o.strings[i]
        l = pauli_weight(p)
        if l == k || (atmost && l <= k)
            push!(o2.coeffs, o.coeffs[i])
            push!(o2.strings, p)
        end
    end
    return o2
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
function trim(o::AbstractOperator, max_strings::Int; keepnorm::Bool=false, keep::Operator=Operator(0))
    if length(o) <= max_strings
        return deepcopy(o)
    end

    # keep the N first indices
    i = collect(partialsortperm(o.coeffs, 1:max_strings; rev=true, by=abs))

    # add the string to keep in case there was a specified string to keep
    if length(keep) > 0
        for tau in 1:length(keep) #for each string tau in the keep operator
            # we check if tau is in o and has been removed
            p_keep = keep.strings[tau]
            j = findfirst(==(p_keep), o.strings)
            if !isnothing(j) && !(j in i)
                push!(i, j)
            end
        end
    end

    o1 = typeof(o)(o.strings[i], o.coeffs[i])

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
    i = Int[]
    for k in 1:length(o)
        p = 1 - exp(-alpha * abs(o.coeffs[k]))
        if rand() < p
            push!(i, k)
        end
    end
    o1 = typeof(o)(o.strings[i], o.coeffs[i])
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
    o2 = zero(o)
    for i in 1:length(o)
        if abs(o.coeffs[i]) > epsilon
            push!(o2.coeffs, o.coeffs[i])
            push!(o2.strings, o.strings[i])
        end
    end
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
    I = findall(p -> zcount(p) == 0 && ycount(p) == 0, o.strings)
    return typeof(o)(o.strings[I], o.coeffs[I])
end


"""
    ypart(o::AbstractOperator)

Keep strings with only Y and identity
"""
function ypart(o::AbstractOperator)
    I = findall(p -> zcount(p) == 0 && xcount(p) == 0, o.strings)
    return typeof(o)(o.strings[I], o.coeffs[I])
end
