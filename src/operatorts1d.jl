
using Dictionaries

"""
    OperatorTS1D{P<:PauliString,T<:Number} <: AbstractOperator

A concrete type representing a 1D translationally invariant operator as a sum of Pauli strings.
The operator is represented as a representative vector of Pauli strings and their corresponding coefficients,
which are implicitly repeated to form the full operator.
The type parameters `P` and `T` specify the type of the Pauli strings and the type of the coefficients, respectively.
"""
struct OperatorTS1D{P<:PauliString,T<:Number} <: AbstractOperator
    strings::Vector{P}
    coeffs::Vector{T}
end

"""
    OperatorTS1D(N::Integer)

Initialize a zero 1D translation-invariant operator on `N` qubits.
"""
OperatorTS1D(N::Integer) = OperatorTS1D{paulistringtype(N),ComplexF64}()
OperatorTS1D{P,T}() where {P,T} = OperatorTS1D{P,T}(P[], T[])

function OperatorTS1D(N::Int, v::Vector{T}, w::Vector{T}, coef::Vector{Complex{Float64}}) where {T<:Unsigned}
    length(v) == length(w) == length(coef) || error("v, w, and coef must have the same length")
    P = paulistringtype(N)
    strings = P.(v, w)
    return OperatorTS1D{P,ComplexF64}(strings, coef)
end

OperatorTS1D(pauli::AbstractString) = OperatorTS1D{paulistringtype(length(pauli))}(pauli)
OperatorTS1D{P}(pauli::AbstractString) where {P} = OperatorTS1D{P,ComplexF64}(pauli)
function OperatorTS1D{P,T}(pauli::AbstractString) where {P,T}
    s = P(pauli)
    c = T((1.0im)^ycount(s))
    return OperatorTS1D{P,T}([s], [c])
end

OperatorTS1D(o::OperatorTS1D) = OperatorTS1D(copy(o.strings), copy(o.coeffs))

paulistringtype(::Type{<:OperatorTS1D{P}}) where {P} = P
scalartype(::Type{OperatorTS1D{P,T}}) where {P,T} = T

Base.one(::Type{O}) where {O<:OperatorTS1D} = O([one(paulistringtype(O))], [one(scalartype(O)) / qubitlength(O)])
Base.zero(::Type{O}) where {O<:OperatorTS1D} = O()


"""
    OperatorTS1D(o::Operator; full=true)

Initialize a 1D translation invariant operator from an Operator
\$O=\\sum_i o_i O_i\$ where \$O_i=T_i(O_0)\$ and \$T_i\$ is the i-sites translation operator.
Set full=true if passing \$O\$, an Operator that is supported on the whole chain (i.e converting from a translation symmetric [`Operator`](@ref))
Set full=false if passing \$O_0\$,a local term o such that the full operator is \$O=\\sum_i o_i T_i(O_0)\$
"""
function OperatorTS1D(o::Operator; full=true)
    if full && !is_ts(o)
        error("o is not translation symmetric. If you want to initialize an OperatorTS1D only with its local part H_0, then set full=false")
    end
    o2 = shift_left(o)
    full && (o2 /= qubitlength(o))
    o3 = OperatorTS1D(o2.strings, o2.coeffs)
end



"""
    Operator(o::OperatorTS1D)

Convert an OperatorTS1D to an Operator
"""
function Operator(o::OperatorTS1D; rs=true)
    rs && (o = resum(o))
    return Operator(o.strings, o.coeffs)
end




"""
    is_ts(o::Operator)

return true if o is translation symmetric
"""
function is_ts(o::Operator)
    for i in 1:qubitlength(o)
        if opnorm(o - shift(o, i)) / opnorm(o) > 1e-10
            return false
        end
    end
    return true
end



"""rotate left the first n bits of x by r"""
function rotate_lower(x::Unsigned, n::Int, r::Int)
    @assert r <= n
    mask = (one(x) << n) - one(x)
    lower_bits = x & mask
    rotated_bits = (lower_bits >> r) | (lower_bits << (n - r))
    rotated_bits &= mask
    return (x & ~mask) | rotated_bits
end
function rotate_lower(p::PauliString{N,T}, r::Int) where {N,T}
    return PauliString{N,T}(rotate_lower(p.v, N, r), rotate_lower(p.w, N, r))
end



"""
    rotate(o::Operator, r::Int)

Rotate (translate/shift) left the qubits of `O` by `r`
"""
rotate(o::Operator, r::Int) = shift(o, r)

shift(o::Operator, r::Int) = Operator(shift.(o.strings, r), copy(o.coeffs))

"""shift the string v,w so it starts on site 1"""
shift_left(p::PauliString) = maximum(Base.Fix1(rotate_lower, p), 0:qubitlength(p)-1)

"""
    shift_left(O::Operator)

Shift evey string left so they start on site 1. This usefull for using translation symmetry in 1D systems
# Example
```julia
A = Operator(4)
A += "XYZ1"
A += "11ZZ"
A += "1XZY"
A += "ZZ11"
```
```julia
julia> shift_left(A)
(1.0 - 0.0im) XZY1
(1.0 - 0.0im) XYZ1
(2.0 + 0.0im) ZZ11
```
"""
shift_left(O::Operator) = compress(Operator(shift_left.(O.strings), copy(O.coeffs)))
shift_left(O::OperatorTS1D) = compress(OperatorTS1D(shift_left.(O.strings), copy(O.coeffs)))

shift1(O::Operator) = shift_left(O)


function resum(o::OperatorTS1D)
    o2 = Operator(similar(o.strings, 0), similar(o.coeffs, 0))
    for i in 1:qubitlength(o)
        oi = Operator(copy(o.strings), copy(o.coeffs))
        o2 += shift(oi, i)
    end
    return o2
end


Base.:+(a::Number, o::OperatorTS1D) = OperatorTS1D(Operator(o, rs=false) + a / qubitlength(o); full=false)
Base.:+(o::OperatorTS1D, a::Number) = a + o

function binary_kernel(op, A::OperatorTS1D, B::OperatorTS1D; epsilon::Real=0, maxlength::Int=1000)
    checklength(A, B)
    N = qubitlength(A)

    d = emptydict(A)
    p1s, c1s = A.strings, A.coeffs
    p2s, c2s = B.strings, B.coeffs

    # check lengths to safely use `@inbounds`
    length(p1s) == length(c1s) || throw(DimensionMismatch("strings and coefficients must have the same length"))
    length(p2s) == length(c2s) || throw(DimensionMismatch("strings and coefficients must have the same length"))

    # core kernel logic
    @inbounds for i1 in eachindex(p1s)
        p1, c1 = p1s[i1], c1s[i1]
        for i2 in eachindex(p2s)
            p2, c2 = p2s[i2], c2s[i2]
            for k in 0:(N-1)
                p, k = op(p1, rotate_lower(p2, k))
                c = c1 * c2 * k
                if (k != 0) && (abs(c) > epsilon) && pauli_weight(p) < maxlength
                    setwith!(+, d, shift_left(p), c)
                end
            end
        end
    end

    return OperatorTS1D(collect(keys(d)), collect(values(d)))
end

Base.:*(o1::OperatorTS1D, o2::OperatorTS1D; kwargs...) = binary_kernel(prod, o1, o2; kwargs...)
commutator(o1::OperatorTS1D, o2::OperatorTS1D; kwargs...) = binary_kernel(commutator, o1, o2; kwargs...)
anticommutator(o1::OperatorTS1D, o2::OperatorTS1D; kwargs...) = binary_kernel(anticommutator, o1, o2; kwargs...)


trace(o::OperatorTS1D) = trace(Operator(o; rs=false)) * qubitlength(o)
opnorm(o::OperatorTS1D) = opnorm(Operator(o; rs=false)) * sqrt(qubitlength(o))
