using Dictionaries

"""
    OperatorTS2D{P<:PauliStringTS,T<:Number,L1} <: AbstractOperator

A concrete type representing a 2D translationally invariant operator as a sum of Pauli strings.
The operator is represented as a representative vector of Pauli strings and their corresponding coefficients,
which are implicitly repeated to form the full operator.
The type parameters `P` and `T` specify the type of the Pauli strings and the type of the coefficients, respectively.
The parameter `L1` will keep track of the extent of the lattice in the a1 direction following the column-major representation,
that is, the flattened index is `i + j*L1`, for a lattice of size L1 * L2.
"""
struct OperatorTS2D{P,T} <: AbstractOperator
    strings::Vector{P}
    coeffs::Vector{T}
end

"""
    OperatorTS2D(N::Integer, L1::Integer)

Initialize a zero 2D translation-invariant operator on `N` qubits, with extent of `L1` in the ``a_1`` direction.
"""
OperatorTS2D(N::Integer, L1::Integer) = (N % L1 == 0) ? OperatorTS2D{paulistringtype((L1, N÷L1)),ComplexF64}() : error("N must be divisible by L1")
OperatorTS2D{P,T}() where {P,T} = OperatorTS2D{P,T}(P[], T[])

function OperatorTS2D(N::Int, L1::Int, v::Vector{T}, w::Vector{T}, coef::Vector{Complex{Float64}}) where {T<:Unsigned}
    length(v) == length(w) == length(coef) || error("v, w, and coef must have the same length")
    P = paulistringtype(N)
    strings = P.(v, w)
    return OperatorTS2D{P,ComplexF64,L1}(strings, coef)
end

OperatorTS2D(pauli::AbstractString, L1::Integer) = OperatorTS2D{paulistringtype((L1,length(pauli)÷L1)), ComplexF64}(pauli)
OperatorTS2D{P}(pauli::AbstractString) where {P} = OperatorTS2D{P,ComplexF64}(pauli)
function OperatorTS2D{P,T}(pauli::AbstractString) where {P,T}
    s = P(pauli)
    c = T((1.0im)^ycount(s))
    return OperatorTS2D{P,T}([s], [c])
end

OperatorTS2D(o::OperatorTS2D) = typeof(o)(copy(o.strings), copy(o.coeffs))

@doc raw"""
    OperatorTS2D(o::Operator, L1::Int; full=true)

Initialize a 2D translation invariant operator from an Operator
$O=%%sum_{i,j} o_{i,j} O_{i,j}$ where $O_{i,j}=T_{a_1}^i T_{a_2}^j(O_0)$ is the application of translation operators in the direction of the basis vectors $a_1$ and $a_2$,
$i$ and $j$ times respectively.
`L1` is the number of sites in the $a_1$ direction. For example, if the lattice is 2x3, `L1 = 2` and `L2 = 3`, and the lattice is flattened in column-major order:
(x,y) -> (1,1),(2,1),(1,2),(2,2),(1,3),(2,3).

Set `full=true` if passing $O$, an Operator that is supported on the whole lattice (i.e converting from a translation symmetric [`Operator`](@ref))
Set `full=false` if passing $O_0$, a local term o such that the full operator is $O=%%sum_{i,j} o_{i,j} T_a^i T_b^j(O_0)$.
"""
function OperatorTS2D(o::Operator, L1::Int; full=true)
    if full && !is_ts2d(o, L1)
        error("o is not 2d translation symmetric. If you want to initialize an OperatorTS2D only with its local part H_0, then set full=false")
    end

    L2 = qubitlength(o)÷L1

    periodic_strings = PauliStringTS{(L1,L2)}.(o.strings)
    coeffs = full ? o.coeffs/qubitlength(o) : copy(o.coeffs)

    return compress(OperatorTS2D{eltype(periodic_strings),eltype(coeffs)}(periodic_strings, coeffs))
end

# deprecate this?
extent(::OperatorTS2D{<:PauliStringTS{Ls},T}) where {Ls,T} = Ls[1]

paulistringtype(::Type{<:OperatorTS2D{P,T}}) where {P,T} = P
scalartype(::Type{OperatorTS2D{P,T}}) where {P,T} = T

qubitsize(::Type{<:OperatorTS2D{P}}) where {P} = qubitsize(P)
qubitsize(op::OperatorTS2D) = qubitsize(typeof(op))

Base.one(::Type{O}) where {O<:OperatorTS2D} = O([one(paulistringtype(O))], [one(scalartype(O)) / qubitlength(O)])
Base.zero(::Type{O}) where {O<:OperatorTS2D} = O()

Base.copy(o::Union{OperatorTS1D,OperatorTS2D}) = typeof(o)(copy(o.strings), copy(o.coeffs))

"""
    Base.length(o::Operator)
    Base.length(o::OperatorTS1D)
    Base.length(o::OperatorTS2D)

Number of Pauli strings in an operator

# Example
```
julia> A = Operator(4)
julia> A += "X111"
julia> A += "XYZ1"
julia> A += 2, "Y", 4
julia> length(A)
3
```
"""
Base.length(o::Union{Operator,OperatorTS1D,OperatorTS2D}) = length(o.strings)


"""
    Operator(o::OperatorTS2D)

Convert an OperatorTS2D to an Operator
"""
Operator(o::OperatorTS2D; rs=true) = rs ? resum(o) : Operator(representative.(o.strings), o.coeffs)

"""
    is_ts2d(o::Operator, L1::Int)

Return true if o is 2d translation symmetric.
"""
function is_ts2d(o::Operator, L1::Int)
    @assert qubitlength(o) % L1 == 0
    L2 = qubitlength(o) ÷ L1
    for i in 0:L1-1
        for j in 0:L2-1
            if opnorm(o - shift(o, (L1, L2), (i, j))) / opnorm(o) > 1e-10
                return false
            end
        end
    end
    return true
end

shift(o::Operator, Ls::Tuple, shifts::Tuple) = shift(o, Val(Ls), shifts)
shift(o::Operator, ::Val{Ls}, shifts::Tuple) where {Ls} =
    compress(Operator(shift.(o.strings, Val(Ls), Ref(shifts)), copy(o.coeffs)))

function shift(o::Operator, r1::Integer, r2::Integer, L1::Integer)
    Base.depwarn("use shift(o, (L1, L2), (r1, r2)) instead", :shift)
    shift(o, (L1, qubitlength(o)÷L1), (r1, r2))
end

function resum(o::OperatorTS2D)
    L1, L2 = qubitsize(o)
    rep_op = Operator(representative.(o.strings), copy(o.coeffs))

    op2 = Operator(similar(rep_op.strings, 0), similar(rep_op.coeffs, 0))
    for i in 0:L1-1
        for j in 0:L2-1
            op2 += shift(rep_op, Val((L1,L2)), (i, j))
        end
    end
    return op2
end


Base.:+(a::Number, o::OperatorTS2D) = OperatorTS2D(Operator(o, rs=false) + a / qubitlength(o), extent(o); full=false)
Base.:+(o::OperatorTS2D, a::Number) = a + o

function binary_kernel(op, A::OperatorTS2D, B::OperatorTS2D; epsilon::Real=0, maxlength::Int=1000)
    checklength(A, B)
    L1, L2 = qubitsize(A)

    d = emptydict(A)
    p1s, c1s = A.strings, A.coeffs
    p2s, c2s = B.strings, B.coeffs

    # check lengths to safely use `@inbounds`
    length(p1s) == length(c1s) || throw(DimensionMismatch("strings and coefficients must have the same length"))
    length(p2s) == length(c2s) || throw(DimensionMismatch("strings and coefficients must have the same length"))

    # core kernel logic
    @inbounds for (p1, c1) in zip(p1s, c1s)
        rep1 = representative(p1)
        for (p2, c2) in zip(p2s, c2s)
            rep2 = representative(p2)
            for i in 0:L1-1
                for j in 0:L2-1
                    p, k = op(rep1, shift(rep2, Val((L1,L2)), (i,j)))
                    c = c1 * c2 * k
                    if (k != 0) && (abs(c) > epsilon) && pauli_weight(p) < maxlength
                        setwith!(+, d, PauliStringTS{(L1,L2)}(p), c)
                    end
                end
            end
        end
    end

    return typeof(A)(collect(keys(d)), collect(values(d)))
end

Base.:*(o1::OperatorTS2D, o2::OperatorTS2D; kwargs...) = binary_kernel(prod, o1, o2; kwargs...)
commutator(o1::OperatorTS2D, o2::OperatorTS2D; kwargs...) = binary_kernel(commutator, o1, o2; kwargs...)
anticommutator(o1::OperatorTS2D, o2::OperatorTS2D; kwargs...) = binary_kernel(anticommutator, o1, o2; kwargs...)


function trace(o::OperatorTS2D)
    r = zero(scalartype(o))

    for (p, c) in zip(o.strings, o.coeffs)
        if isone(representative(p))
            r += c
        end
    end
    return r * qubitlength(o) * 2^qubitlength(o)
end

opnorm(o::OperatorTS2D) = sqrt(trace_product(dagger(o),o))
