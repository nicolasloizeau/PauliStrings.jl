
import LinearAlgebra as la


"""number of Y in a pauli string"""
function getdelta(pauli::String)
    return count(c -> c == 'Y', pauli)
end

"""
get v and w from ref
https://journals.aps.org/pra/abstract/10.1103/PhysRevA.68.042318
"""
function getvw(pauli::String)
    T = uinttype(length(pauli))
    v::T = 0
    w::T = 0
    two::T = 2
    for k in 1:length(pauli)
        if pauli[k] == 'X'
            w += two^(k - 1)
        end
        if pauli[k] == 'Z'
            v += two^(k - 1)
        end
        if pauli[k] == 'Y'
            w += two^(k - 1)
            v += two^(k - 1)
        end
    end
    return v, w
end

"""
    string_to_vw(pauli::String)

Convert a string to two integers (v,w)
"""
string_to_vw(pauli::String) = getvw(pauli)


"""
    set_coeffs(o::Operator, coefs::Vector{T}) where T <: Number

Sets the coefficient of `o` to `coefs`. Inplace.

```julia
A = Operator(4)
A += 2, "1XXY"
A += 3, "11Z1"
```
```
julia> A
(3.0 + 0.0im) 11Z1
(2.0 - 0.0im) 1XXY
julia> set_coeffs(A, [5,6])
julia> A
(5.0 + 0.0im) 11Z1
(6.0 - 0.0im) 1XXY
```
"""
function set_coeffs(o::AbstractOperator, coeffs::Vector{T}) where {T<:Number}
    length(o) != length(coeffs) && error("length(o) != length(coefs)")
    for i in 1:length(o)
        o.coeffs[i] = (1im)^ycount(o.strings[i]) * coeffs[i]
    end
end

set_coefs(o::AbstractOperator, coeffs::Vector{T}) where {T<:Number} = set_coeffs(o, coefs)

"""
add a pauli string term to an operator
"""
function add_string(o::Operator, pauli::String, J::Number)
    p = paulistringtype(o)(pauli)
    c = (1im)^ycount(p) * J
    push!(o.strings, p)
    push!(o.coeffs, c)
    return o
end

function string_from_inds(ind::Vector{Int})
    l::Vector{Char} = []
    paulis = ['1', 'X', 'Y', 'Z']
    for i in ind
        push!(l, paulis[i+1])
    end
    return join(l)
end


function Base.:+(o::Operator, term::Tuple{Number,Char,Int,Char,Int})
    o1 = deepcopy(o)
    J, Pi, i, Pj, j = term
    pauli = fill('1', qubitlength(o1))
    pauli[i] = Pi
    pauli[j] = Pj
    pauli = join(pauli)
    add_string(o1, pauli, J)
    return compress(o1)
end


function Base.:+(o::Operator, term::Tuple{Number,Char,Int})
    o1 = deepcopy(o)
    J, Pi, i = term
    pauli = fill('1', qubitlength(o1))
    pauli[i] = Pi
    pauli = join(pauli)
    add_string(o1, pauli, J)
    return compress(o1)
end

Base.:+(o::Operator, term::Tuple{Char,Int,Char,Int}) = o + (1, term...)
Base.:+(o::Operator, term::Tuple{Char,Int}) = o + (1, term...)


"""
    Base.:+(o::Operator, args::Tuple{Number, Vararg{Any}})
    Base.:+(o::Operator, args::Tuple{Vararg{Any}})
    Base.:+(o::Operator, term::Tuple{Number, String})
    Base.:+(o::Operator, term::String)

Main functions to contruct spin operators.
Identical signatures are available for `-`.

# Examples
k-local terms can be added by adding a tuple to the operator.
The first element of the tuple is an optional coefficient.
The other element are couples (symbol,site) where symbol can be "X", "Y", "Z", "Sx", "Sy", "Sz", "S+", "S-" and site is an integer specifying the site on wich the symbol is acting.

```julia
A = Operator(4)
A += 2, "X",1,"X",2
A += 3, "Y",1,"X",2
A += "X",3,"X",4
A += 4,"Z",3
A += 5.2,"X",1,"Y",2,"Z",3
```
```
julia> A
(4.0 + 0.0im) 11Z1
(3.0 - 0.0im) YX11
(1.0 + 0.0im) 11XX
(2.0 + 0.0im) XX11
(5.2 - 0.0im) XYZ1
```

Full strings can also be added:
```julia
A = Operator(4)
A += 2, "1XXY"
A += 2im, "11Z1"
```
```
julia> A
(0.0 + 2.0im) 11Z1
(2.0 - 0.0im) 1XXY
```
"""
function Base.:+(o::Operator, args::Tuple{Number,Vararg{Any}})
    term = one(o)
    c = args[1]
    for i in 2:2:length(args)
        o2 = zero(o)
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
    return compress(o + c * term)
end
Base.:+(o::Operator, args::Tuple{Vararg{Any}}) = o + (1, args...)
Base.:-(o::Operator, args::Tuple{Number,Vararg{Any}}) = o + (-args[1], args[2:end]...)
Base.:-(o::Operator, args::Tuple{Vararg{Any}}) = o + (-1, args...)



function Base.:+(o::Operator, term::Tuple{Number,String})
    o1 = deepcopy(o)
    c, pauli = term
    if qubitlength(o1) != length(pauli)
        error("The string needs to be of the same size as the operator")
    end
    add_string(o1, pauli, c)
    return compress(o1)
end

Base.:+(o::Operator, term::String) = o + (1, term)
Base.:-(o::Operator, term::String) = o + (-1, term)
Base.:-(o::Operator, term::Tuple{Number,String}) = o + (-term[1], term[2])

Base.:+(o::Operator, term::Tuple{Number,Vector{Int}}) = o + (term[1], string_from_inds(term[2]))
Base.:-(o::Operator, term::Tuple{Number,Vector{Int}}) = o - (term[1], string_from_inds(term[2]))

Base.:+(o::Operator, term::Vector{Int}) = o + (1, string_from_inds(term))
Base.:-(o::Operator, term::Vector{Int}) = o - (1, string_from_inds(term))


"""true if bit i of n is set"""
function bit(n::Unsigned, i::Int)
    return (n & (one(n) << (i - one(n)))) != 0
end

"""
    vw_to_string(v::Int, w::Unsigned, N::Unsigned)

convert v,w to a string and a phase
"""
function vw_to_string(v::Unsigned, w::Unsigned, N::Int)
    string::String = ""
    phase::Complex{Float64} = 1
    T = uinttype(N)
    for i in 1:N
        if !bit(v, i) && !bit(w, i)
            string *= '1'
        end
        if !bit(v, i) && bit(w, i)
            string *= 'X'
        end
        if bit(v, i) && !bit(w, i)
            string *= 'Z'
        end
        if bit(v, i) && bit(w, i)
            string *= 'Y'
            phase *= 1im
        end
    end
    return string, phase
end



function Base.show(io::IO, o::AbstractOperator)
    N = qubitlength(o)
    t = eltype(o.coeffs)
    o = sort(o)
    for (p,c) in zip(o.strings, o.coeffs)
        phase = 1im^ycount(p)
        c /= phase

        pauli = string(p)

        if t == ComplexF64
            prefix = "($(round(c, digits=10))) "
        else
            prefix = "($c) "
        end

        space = "\n"*repeat(" ", length(prefix))
        print(io, prefix)
        join(io, split(pauli, '\n'), space)
        println(io)
    end
end

Base.show(io::IO, o::AbstractPauliString) = print(io, string(o))

"""
    get_coeffs(o::Operator)

Return the list of coefficient in front of each strings.
"""
function get_coeffs(o::Operator)
    return [o.coeffs[i] / (1im)^ycount(o.strings[i]) for i in 1:length(o)]
end

get_coefs(o) = get_coeffs(o)
get_coef(o) = get_coeff(o)

"""
    op_to_strings(o::Operator)

takes an operator,
return (coefs, strings)
where coefs is a list of numbers and strings is a list of pauli string
coefs[i] multiply strings[i]
"""
function op_to_strings(o::Operator)
    strings::Vector{String} = []
    coefs::Vector{Complex{Float64}} = []
    for i in 1:length(o)
        pauli, phase = vw_to_string(o.strings[i].v, o.strings[i].w, qubitlength(o))
        push!(coefs, o.coeffs[i] / phase)
        push!(strings, pauli)
    end
    return coefs, strings
end

s0 = [1 0; 0 1]
sx = [0 1; 1 0]
sy = [0 -im; im 0]
sz = [1 0; 0 -1]
pdict = Dict('1' => s0, 'X' => sx, 'Y' => sy, 'Z' => sz)

function string_to_dense(v, w, N)
    pauli, phase = vw_to_string(v, w, N)
    tau = 1
    for s in pauli
        tau = la.kron(tau, pdict[s])
    end
    return tau, phase
end
string_to_dense(p::PauliString) = string_to_dense(p.v, p.w, qubitlength(p))

"""
    op_to_dense(o::Operator)

Convert an operator to a dense matrix.
"""
function op_to_dense(o::Operator)
    N = qubitlength(o)
    dense = zeros(Complex, 2^N, 2^N)
    for i in 1:length(o)
        p = o.strings[i]
        tau, phase = string_to_dense(p)
        dense .+= tau * o.coeffs[i] / phase
    end
    return dense
end


"""
    Matrix(o::Operator)

Convert an operator to a dense matrix.
"""
Base.Matrix(o::Operator) = op_to_dense(o)


"""
    get_coeff(o::Operator{P}, p::P) where {P}

Return the coefficient of the string p in o.
"""
function get_coeff(o::Operator{P}, p::P) where {P}
    id = findfirst(==(p), o.strings)
    return isnothing(id) ? zero(scalartype(o)) : (o.coeffs[id] / (1im)^ycount(o.strings[id]))
end


"""
    get_pauli(o::Operator, i::Int)

Return an operator that represent the i-th pauli string of `o'.
Does not return the string multiplied by the coefficient. Only the string.
"""
function get_pauli(o::Operator, i::Int)
    p = o.strings[i]
    return typeof(o)([p], [(1.0im)^ycount(p)])
end


op_to_dense(o::OperatorTS) = op_to_dense(resum(o))

"""
    p"paulistring"

String macro to create a pauli string.
"""
macro p_str(pauli)
    return PauliString(pauli)
end


"""
    vw_in_o(v::Unsigned, w::Unsigned, o::Operator)

Return true is string (v,w) is in o
"""
function vw_in_o(v::Unsigned, w::Unsigned, o::Operator)
    for string in o.strings
        if v == string.v && w == string.w
            return true
        end
    end
    return false
end



function Base.sort(o::Operator)
    i = sortperm(eachindex(o.coeffs), by=i->(abs(o.coeffs[i]), o.strings[i]))
    o2 = typeof(o)(o.strings[i], o.coeffs[i])
    return o2
end
