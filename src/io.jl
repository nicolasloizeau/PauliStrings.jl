

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


"""print an operator"""
function Base.show(io::IO, o::Operator)
    for i in 1:length(o.v)
        pauli, phase = vw_to_string(o.v[i],o.w[i],o.N)
        c = o.coef[i]/phase
        println(io, "($(round(c, digits=10))) ", pauli)
    end
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
