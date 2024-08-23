
"""
remove all terms of length > N
keep all terms of length <= N
i.e remove all M-local terms with M>N
"""
function truncate(o::Operator, N::Int; keepnorm::Bool = false)
    o2 = Operator(o.N)
    for i in 1:length(o)
        v = o.v[i]
        w = o.w[i]
        if pauli_weight(v,w)<=N
            push!(o2.coef, o.coef[i])
            push!(o2.v, v)
            push!(o2.w, w)
        end
    end
    if keepnorm
        return o2*opnorm(o)/opnorm(o2)
    end
    return o2
end


"""
keep the first N terms with larger coeficients
keepnorm : set to true to keep the norm of o
keep : an operator that specify a set of strings that cannot be removed
"""
function trim(o::Operator, N::Int; keepnorm::Bool = false, keep::Operator=Operator(N))
    if length(o)<=N
        return deepcopy(o)
    end
    # keep the N first indices
    i = sortperm(abs.(o.coef), rev=true)[1:N]
    # add the string to keep in case there was a specified string to keep
    if length(keep) > 0
        for tau in 1:length(keep) #for each string tau in the keep operator
            # we check if tau is in o and has been removed
            j = posvw(keep.v[tau], keep.w[tau], o)
            if !(j in i) && j!=0
                push!(i, j)
            end
        end
    end
    o1 = Operator(o.N, o.v[i], o.w[i], o.coef[i] )
    if keepnorm
        return o1*opnorm(o)/opnorm(o1)
    end
    return o1
end

"""
keep terms with probability 1-exp(-alpha*abs(c)) where c is the weight of the term
"""
function prune(o::Operator, alpha::Real; keepnorm::Bool = false)
    i = Int[]
    for k in 1:length(o)
        p = 1-exp(-alpha*abs(o.coef[k]))
        if rand() < p
            push!(i,k)
        end
    end
    o1 = Operator(o.N, o.v[i], o.w[i], o.coef[i] )
    if keepnorm
        return o1*opnorm(o)/opnorm(o1)
    end
    return o1
end

"""remove all terms with coef < epsilon"""
function cutoff(o::Operator, epsilon::Real; keepnorm::Bool = false)
    o2 = Operator(o.N)
    for i in 1:length(o)
        if abs(o.coef[i])>epsilon
            push!(o2.coef, o.coef[i])
            push!(o2.v, o.v[i])
            push!(o2.w, o.w[i])
        end
    end
    if keepnorm
        return o2*opnorm(o)/opnorm(o2)
    end
    return o2
end



"""
add noise that make the long string decays
https://arxiv.org/pdf/2407.12768
"""
function add_noise(o::Operator, g::Real)
    o2 = deepcopy(o)
    for i in 1:length(o)
        o2.coef[i] *= exp(-pauli_weight(o.v[i],o.w[i])*g)
    end
    return o2
end
