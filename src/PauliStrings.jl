module PauliStrings


using Random
using LinearAlgebra
using ProgressBars
using Dictionaries

rng = MersenneTwister(0)

"""
operator as a sum of pauli string encoded like in
https://journals.aps.org/pra/abstract/10.1103/PhysRevA.68.042318
intialized as : o=Operator(N)
where N is the number of qubits
"""
mutable struct Operator
    N::Int
    v::Vector{Int}
    w::Vector{Int}
    coef::Vector{Complex{Float64}}
    function Operator(N::Int)
        new(N, Int[], Int[], Complex{Float64}[])
    end
    function Operator(N::Int, v::Vector{Int}, w::Vector{Int}, coef::Vector{Complex{Float64}})
        new(N, v, w, coef)
    end
end

"""number of Y in a pauli string"""
function getdelta(pauli::String)
    return count(c -> c == 'Y', pauli)
end

"""
get v and w from ref
https://journals.aps.org/pra/abstract/10.1103/PhysRevA.68.042318
"""
function getvw(pauli::String)
    v::Int = 0
    w::Int = 0
    for k in 1:length(pauli)
        if pauli[k] == 'X'
            w += 2^(k-1)
        end
        if pauli[k] == 'Z'
            v += 2^(k-1)
        end
        if pauli[k] == 'Y'
            w += 2^(k-1)
            v += 2^(k-1)
        end
    end
    return v,w
end

"""
add a pauli string term to an operator
"""
function add_string(o::Operator, pauli::String, J::Number)
    v,w = getvw(pauli)
    c = (1im)^getdelta(pauli)*J
    push!(o.v, v)
    push!(o.w, w)
    push!(o.coef, c)
end

function string_from_inds(ind::Vector{Int})
    l::Vector{Char} = []
    paulis = ['1', 'X', 'Y', 'Z']
    for i in ind
        push!(l, paulis[i+1])
    end
    return join(l)
end


function Base.:+(o::Operator, term::Tuple{Number, Char, Int, Char, Int})
    o1 = deepcopy(o)
    J, Pi, i, Pj, j = term
    pauli = fill('1', o1.N)
    pauli[i] = Pi
    pauli[j] = Pj
    pauli = join(pauli)
    add_string(o1, pauli, J)
    return compress(o1)
end


function Base.:+(o::Operator, term::Tuple{Number, Char, Int})
    o1 = deepcopy(o)
    J, Pi, i = term
    pauli = fill('1', o1.N)
    pauli[i] = Pi
    pauli = join(pauli)
    add_string(o1, pauli, J)
    return compress(o1)
end

Base.:+(o::Operator, term::Tuple{Char, Int, Char, Int}) = o+(1, term...)
Base.:+(o::Operator, term::Tuple{Char, Int}) = o+(1, term...)


"""
Main function to contruct spin operators in a user friendly way
# Examples
```
O = Operator(4)
O += 1, "X",1,"X",2
O += 2, "S-",1,"S+",4
```
"""
function Base.:+(o::Operator, args::Tuple{Number, Vararg{Any}})
    term = Operator(o.N)+1
    c = args[1]
    for i in 2:2:length(args)
        o2 = Operator(o.N)
        symbol = args[i]::String
        site = args[i+1]::Int
        if occursin(symbol, "XYZ")
            o2 += only(symbol), site
        elseif occursin(symbol, "SxSySz")
            o2 += 0.5, only(uppercase(symbol[2])), site
        elseif symbol == "S+"
            o2 += 0.5, 'X', site
            o2 += 0.5im, 'Y', site
        elseif symbol == "S-"
            o2 += 0.5, 'X', site
            o2 += -0.5im, 'Y', site
        else
            error("Allowed operators: X,Y,Z,Sx,Sy,Sz,S-,S+")
        end
        term *= o2
    end
    return compress(o+c*term)
end
Base.:+(o::Operator, args::Tuple{Vararg{Any}}) = o+(1, args...)
Base.:-(o::Operator, args::Tuple{Number, Vararg{Any}}) = o+(-args[1], args[2:end]...)
Base.:-(o::Operator, args::Tuple{Vararg{Any}}) = o+(-1, args...)



function Base.:+(o::Operator, term::Tuple{Number, String})
    o1 = deepcopy(o)
    c, pauli = term
    if o1.N != length(pauli)
        error("The string needs to be of the same size as the operator")
    end
    add_string(o1, pauli, c)
    return compress(o1)
end

Base.:+(o::Operator, term::String) = o+(1,term)
Base.:-(o::Operator, term::String) = o+(-1,term)
Base.:-(o::Operator, term::Tuple{Number, String}) = o+(-term[1], term[2])

Base.:+(o::Operator, term::Tuple{Number, Vector{Int}}) = o+(term[1], string_from_inds(term[2]))
Base.:-(o::Operator, term::Tuple{Number, Vector{Int}}) = o-(term[1], string_from_inds(term[2]))

Base.:+(o::Operator, term::Vector{Int}) = o+(1, string_from_inds(term))
Base.:-(o::Operator, term::Vector{Int}) = o-(1, string_from_inds(term))


"""true if bit i of n is set"""
function bit(n::Integer, i::Integer)
    return (n & (1 << (i - 1))) != 0
end

"""convert v,w to a string and a phase"""
function vw_to_string(v::Int, w::Int, N::Int)
    string::String = ""
    phase::Complex{Float64} = 1
    for i in 1:N
        if !bit(v,i) && !bit(w,i)
            string *= '1'
        end
        if !bit(v,i) && bit(w,i)
            string *= 'X'
        end
        if bit(v,i) && !bit(w,i)
            string *= 'Z'
        end
        if bit(v,i) && bit(w,i)
            string *= 'Y'
            phase *= 1im
        end
    end
    return string, phase
end

"""print an operator"""
function Base.show(io::IO, o::Operator)
    for i in 1:length(o.v)
        pauli, phase = vw_to_string(o.v[i],o.w[i],o.N)
        c = o.coef[i]/phase
        println(io, "($(round(c, digits=10))) ", pauli)
    end
end


function Base.:+(o1::Operator, o2::Operator)
    if o1.N != o2.N
        error("Adding operators of different dimention")
    end
    o3 = Operator(o1.N)
    o3.v = vcat(o1.v, o2.v)
    o3.w = vcat(o1.w, o2.w)
    o3.coef = vcat(o1.coef, o2.coef)
    return compress(o3)
end


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

"""random 2-local operator"""
function rand_local2(N::Int)
    o = Operator(N)
    for i in 1:N
        for j in 1:N
            for k in ['X', 'Y', 'Z']
                for l in ['X', 'Y', 'Z']
                    o += (randn(rng, Float64), k, i, l, j)
                end
            end
        end
    end
    return compress(o)
end

"""random 1-local operator"""
function rand_local1(N::Int)
    o = Operator(N)
    for i in 1:N
        for k in ['X', 'Y', 'Z']
            o += (randn(rng, Float64), k, i)
        end
    end
    return compress(o)
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

"""number of pauli strings in an operator"""
function Base.length(container::Operator)
    return length(container.v)
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
takes an operator,
return (coefs, strings)
where coefs is a list of numbers and strings is a list of pauli string
coefs[i] multiply strings[i]
"""
function op_to_strings(o::Operator)
    strings::Vector{String} = []
    coefs::Vector{Complex{Float64}} = []
    for i in 1:length(o)
        pauli,phase = vw_to_string(o.v[i], o.w[i], o.N)
        push!(coefs, o.coef[i]/phase)
        push!(strings, pauli)
    end
    return coefs, strings
end

"""number of non unit paulis in a string encoded by v,w"""
function pauli_weight(v::Int,w::Int)
    return count_ones(v | w)
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

"""returns the position of (v,w) in O. return 0 if (v,w) not in O"""
function posvw(v,w,O)
    for i in 1:length(O)
        if O.v[i] == v && O.w[i] == w
            return i
        end
    end
    return 0
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
remove all terms of length >= N
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

"""
Runge–Kutta-4 with time independant Hamiltonian
heisenberg : set to true if evolving an observable
return : rho(t+dt)
"""
function rk4(H::Operator, O::Operator, dt::Real; hbar::Real=1, heisenberg=false, M=2^20, keep::Operator=Operator(N))
    s = 1
    if heisenberg
        s = -1
    end
    k1 = -s*1im/hbar*com(H, O)
    k1 = trim(k1, M; keep=keep)
    k2 = -s*1im/hbar*com(H, O+dt*k1/2)
    k2 = trim(k2, M; keep=keep)
    k3 = -s*1im/hbar*com(H, O+dt*k2/2)
    k3 = trim(k3, M; keep=keep)
    k4 = -s*1im/hbar*com(H, O+dt*k3)
    k4 = trim(k4, M; keep=keep)
    return O+(k1+2*k2+2*k3+k4)*dt/6
end

"""
Runge–Kutta-4 with time dependant Hamiltonian
H : function that takes a time and returns an operator
heisenberg : set to true if evolving an observable
return : rho(t+dt)
"""
function rk4(H::Function, O::Operator, dt::Real, t::Real; hbar::Real=1, heisenberg=false)
    s = 1
    if heisenberg
        s = -1
    end
    k1 = -s*1im/hbar*com(H(t), O)
    k2 = -s*1im/hbar*com(H(t+dt/2), O+dt*k1/2)
    k3 = -s*1im/hbar*com(H(t+dt/2), O+dt*k2/2)
    k4 = -s*1im/hbar*com(H(t+dt), O+dt*k3)
    return O+(k1+2*k2+2*k3+k4)*dt/6
end

function norm_lanczos(O::Operator)
    return opnorm(O)/sqrt(2^O.N)
end

"
https://journals.aps.org/prx/pdf/10.1103/PhysRevX.9.041017 equation 4
H : hamiltonian MPO
O : operator MPO
steps : numer of lanczos steps
nterms : maximum number of terms in the operator. Used by trim at every step
"
function lanczos(H::Operator, O::Operator, steps::Int, nterms::Int; keepnorm=true, maxlength=1000)
    N = H.N
    O0 = deepcopy(O)
    b = norm_lanczos(com(H, O0))
    O1 = com(H, O0)/b
    bs = [b]
    for n in ProgressBar(0:steps-2)
        LHO = com(H, O1; maxlength=maxlength)
        A = LHO-b*O0
        b = norm_lanczos(A)
        O = A/b
        O = trim(O, nterms; keepnorm=keepnorm)
        O0 = deepcopy(O1)
        O1 = deepcopy(O)
        push!(bs, b)
    end
    return bs
end



"
https://journals.aps.org/prx/pdf/10.1103/PhysRevX.9.041017 equation 4
H : hamiltonian MPO
O : operator MPO
steps : numer of lanczos steps
nterms : maximum number of terms in the operator. Used by trim at every step
"
function lanczoscutoff(H::Operator, O::Operator, steps::Int, epsilon::Real; keepnorm=true, nterms=2^28)
    N = H.N
    O0 = deepcopy(O)
    b = norm_lanczos(com(H, O0; epsilon=epsilon))
    O1 = com(H, O0; epsilon=epsilon)/b
    bs = [b]
    for n in ProgressBar(0:steps-2)
        LHO = com(H, O1; epsilon=epsilon)
        A = LHO-b*O0
        b = norm_lanczos(A)
        O = A/b
        O = trim(O, nterms; keepnorm=keepnorm)
        O0 = deepcopy(O1)
        O1 = deepcopy(O)
        push!(bs, b)
    end
    return bs
end

end
