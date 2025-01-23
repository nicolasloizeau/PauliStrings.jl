
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
function PauliStrings.truncate(o::Operator, max_lenght::Int; keepnorm::Bool=false)
    o2 = typeof(o)(o.N)
    for i in 1:length(o)
        v = o.v[i]
        w = o.w[i]
        if pauli_weight(v, w) <= max_lenght
            push!(o2.coef, o.coef[i])
            push!(o2.v, v)
            push!(o2.w, w)
        end
    end
    if keepnorm
        return o2 * opnorm(o) / opnorm(o2)
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
function k_local_part(o::Operator, k::Int; atmost=false)
    o2 = typeof(o)(o.N)
    for i in 1:length(o)
        v = o.v[i]
        w = o.w[i]
        l = pauli_weight(v, w)
        if l == k || (atmost && l <= k)
            push!(o2.coef, o.coef[i])
            push!(o2.v, v)
            push!(o2.w, w)
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
function trim(o::Operator, max_strings::Int; keepnorm::Bool=false, keep::Operator=Operator(0))
    if length(o) <= max_strings
        return deepcopy(o)
    end
    # keep the N first indices
    i = sortperm(abs.(o.coef), rev=true)[1:max_strings]
    # add the string to keep in case there was a specified string to keep
    if length(keep) > 0
        for tau in 1:length(keep) #for each string tau in the keep operator
            # we check if tau is in o and has been removed
            j = posvw(keep.v[tau], keep.w[tau], o)
            if !(j in i) && j != 0
                push!(i, j)
            end
        end
    end
    o1 = typeof(o)(o.N, o.v[i], o.w[i], o.coef[i])
    if keepnorm
        return o1 * opnorm(o) / opnorm(o1)
    end
    return o1
end

"""
     prune(o::Operator, alpha::Real; keepnorm::Bool = false)

Keep terms with probability 1-exp(-alpha*abs(c)) where c is the weight of the term
"""
function prune(o::Operator, alpha::Real; keepnorm::Bool=false)
    i = Int[]
    for k in 1:length(o)
        p = 1 - exp(-alpha * abs(o.coef[k]))
        if rand() < p
            push!(i, k)
        end
    end
    o1 = typeof(o)(o.N, o.v[i], o.w[i], o.coef[i])
    if keepnorm
        return o1 * opnorm(o) / opnorm(o1)
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
function cutoff(o::Operator, epsilon::Real; keepnorm::Bool=false)
    o2 = typeof(o)(o.N)
    for i in 1:length(o)
        if abs(o.coef[i]) > epsilon
            push!(o2.coef, o.coef[i])
            push!(o2.v, o.v[i])
            push!(o2.w, o.w[i])
        end
    end
    if keepnorm
        return o2 * opnorm(o) / opnorm(o2)
    end
    return o2
end



"""
    add_noise(o::Operator, g::Real)

Add depolarizing noise that make the long string decays. `g` is the noise amplitude.
# Example
```julia
A = add_noise(A, 0.1)
```
# Reference
[https://arxiv.org/pdf/2407.12768](https://arxiv.org/pdf/2407.12768)
"""
function add_noise(o::Operator, g::Real)
    o2 = deepcopy(o)
    for i in 1:length(o)
        o2.coef[i] *= exp(-pauli_weight(o.v[i], o.w[i]) * g)
    end
    return o2
end


function participation(o::Operator)
    return sum(o.coef .^ 4) / sum(o.coef .^ 2)^2
end
