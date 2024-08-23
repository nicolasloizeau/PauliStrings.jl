
function add(o1::Operator, o2::Operator)
    if o1.N != o2.N
        error("Adding operators of different dimention")
    end
    d1 = Dict{Tuple{Int, Int}, Complex{Float64}}(zip(zip(o1.v, o1.w), o1.coef))
    d2 = Dict{Tuple{Int, Int}, Complex{Float64}}(zip(zip(o2.v, o2.w), o2.coef))
    d = mergewith(+, d1, d2)
    o3 = Operator(o1.N)
    vw =  collect(keys(d))
    o3.v = [i[1] for i in vw]
    o3.w = [i[2] for i in vw]
    o3.coef =  collect(values(d))
    return o3
end



function Base.:+(o::Operator, a::Number)
    o1 = deepcopy(o)
    i = ione(o)
    if i >=0
        o1.coef[ione(o)] += a
    else
        push!(o1.coef, a)
        push!(o1.v, 0)
        push!(o1.w, 0)
    end
    return o1
end

Base.:+(a::Number, o::Operator) = o+a


function Base.:*(o1::Operator, o2::Operator)
    if o1.N != o2.N
        error("Multiplying operators of different dimention")
    end
    o3 = Operator(o1.N)
    d = UnorderedDictionary{Tuple{Int, Int}, Complex{Float64}}()
    for i in 1:length(o1.v)
        for j in 1:length(o2.v)
            v = o1.v[i] ⊻ o2.v[j]
            w = o1.w[i] ⊻ o2.w[j]
            c = o1.coef[i] * o2.coef[j] * (-1)^count_ones(o1.v[i] & o2.w[j])
            if isassigned(d, (v,w))
                d[(v,w)] += c
            else
                insert!(d, (v,w), c)
            end
        end
    end
    for (v,w) in keys(d)
        push!(o3.v, v)
        push!(o3.w, w)
        push!(o3.coef, d[(v,w)])
    end
    return o3
end


function Base.:*(o::Operator, a::Number)
    o1 = deepcopy(o)
    o1.coef .*= a
    return o1
end

Base.:*(a::Number, o::Operator) = o*a


function Base.:/(o::Operator, a::Number)
    o1 = deepcopy(o)
    o1.coef ./= a
    return o1
end

Base.:-(o::Operator) = -1*o
Base.:-(o1::Operator, o2::Operator) = o1+(-o2)
Base.:-(o::Operator, a::Real) = o+(-a)
Base.:-(a::Real, o::Operator) = a+(-o)


"""
Commutator. This is much faster than doing o1*o2-o2*o1
"""
function com(o1::Operator, o2::Operator; epsilon::Real=0, maxlength::Int=1000)
    if o1.N != o2.N
        error("Commuting operators of different dimention")
    end
    o3 = Operator(o1.N)
    d = UnorderedDictionary{Tuple{Int, Int}, Complex{Float64}}()
    for i in 1:length(o1.v)
        for j in 1:length(o2.v)
            v = o1.v[i] ⊻ o2.v[j]
            w = o1.w[i] ⊻ o2.w[j]
            k = (-1)^count_ones(o1.v[i] & o2.w[j]) - (-1)^count_ones(o1.w[i] & o2.v[j])
            c = o1.coef[i] * o2.coef[j] * k
            if (k != 0) && (abs(c)>epsilon) && pauli_weight(v,w)<maxlength
                if isassigned(d, (v,w))
                    d[(v,w)] += c
                else
                    insert!(d, (v,w), c)
                end
            end
        end
    end
    for (v,w) in keys(d)
        push!(o3.v, v)
        push!(o3.w, w)
        push!(o3.coef, d[(v,w)])
    end
    return o3
end



"""
accumulate terms with the same pauli string
"""
function compress(o::Operator)
    o1 = Operator(o.N)
    vw = Set{Tuple{Int, Int}}(zip(o.v, o.w))
    d = UnorderedDictionary{Tuple{Int, Int}, Complex{Float64}}(vw, zeros(length(vw)))
    for i in 1:length(o)
        v = o.v[i]
        w = o.w[i]
        d[(v,w)] += o.coef[i]
    end
    for (v,w) in keys(d)
        if abs(d[(v,w)])>1e-20
            push!(o1.v, v)
            push!(o1.w, w)
            push!(o1.coef, d[(v,w)])
        end
    end
    return o1
end


"""return the index of the 1 string"""
function ione(o::Operator)
    for i in 1:length(o)
        if o.v[i]==0 && o.w[i]==0
            return i
        end
    end
    return -1
end


function trace(o::Operator)
    t = 0
    for i in 1:length(o.v)
        if o.v[i]==0 && o.w[i]==0
            t += o.coef[i]
        end
    end
    return t*2^o.N
end

"""frobenius norm"""
function opnorm(o::Operator)
    return norm(o.coef)*sqrt(2^o.N)
end


"""conjugate transpose"""
function dagger(o::Operator)
    o1 = deepcopy(o)
    for i in 1:length(o1)
        s = (-1)^count_ones(o1.v[i] & o1.w[i])
        o1.coef[i] = s*conj(o1.coef[i])
    end
    return o1
end



"""
v,w encode a string.
return true if at least one index of keep is non unit in vw
"""
function tokeep(v::Int, w::Int, keep::Vector{Int})
    for i in keep
        if bit(v|w, i)
            return true
        end
    end
    return false
end

"""
partial trace
keep : list of qubits indices to keep starting at 1
note that this still returns an operator of size N and doesnt permute the qubits
this only gets rid of Pauli strings that have no support on keep
and add their coeficient*2^NB to the identity string
"""
function ptrace(o::Operator, keep::Vector{Int})
    o2 = Operator(o.N)
    NA = length(keep)
    NB = o.N-NA
    for i in 1:length(o)
        if tokeep(o.v[i], o.w[i], keep)
            push!(o2.v, o.v[i])
            push!(o2.w, o.w[i])
            push!(o2.coef, o.coef[i])
        else
            o2 += o.coef[i]*2^NB/(1im)^count_ones(o.v[i] & o.w[i])
        end
    end
    return o2
end
