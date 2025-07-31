"""
    PauliStringTS{Ls, T<:Unsigned} <: AbstractPauliString

Type representing a translation symmetric sum of a Pauli string. The tuple `Ls` specifies the period in each dimension.
"""
struct PauliStringTS{Ls, T <: Unsigned} <: AbstractPauliString
    v::T
    w::T
end

periodicpaulistringtype(Ls::NTuple{<:Any, Integer}) = PauliStringTS{Ls, uinttype(Base.prod(Ls))}
qubitlength(::Type{<:PauliStringTS{Ls}}) where {Ls} = Base.prod(Ls)

"""
    qubitsize(::Type{<:PauliStringTS})
    qubitsize(p::PauliStringTS)

Get the tuple of periods of a given `PauliStringTS`.
"""
qubitsize(::Type{<:PauliStringTS{Ls}}) where {Ls} = Ls
qubitsize(p::PauliStringTS) = qubitsize(p)

for count in [:xcount, :ycount, :zcount, :pauli_weight]
    eval(:($count(p::PauliStringTS) = $count(representative(p))))
end

function Base.string(p::PauliStringTS)
    Ls = qubitsize(typeof(p))

    function str(slice)
        if ndims(slice) == 1
            return join(slice)
        elseif ndims(slice) == 2
            return join(join.(eachcol(slice)), "\n")
        end

        segments = str.(eachslice(slice, dims = ((3:ndims(slice))...,)))

        return join(["[:, :, $(join(Tuple(I), ", "))] =\n" * segments[I] for I in CartesianIndices(segments)], "\n\n")
    end
    σij = reshape(collect(vw_to_string(p.v, p.w, Base.prod(Ls))[1]), Ls)
    return str(σij)
end

Base.one(::Type{<:PauliStringTS{Ls, T}}) where {Ls, T} = PauliStringTS{Ls}(one(PauliString{Base.prod(Ls), T}))

"""
    PauliStringTS{Ls}(p::PauliString)

Construct a translation symmetric sum of a Pauli string `p`. `Ls` is a tuple that specifies the periodicity in each dimension.
"""
function PauliStringTS{Ls}(p::PauliString) where {Ls}
    if !(Ls isa Tuple)
        error("Cannot construct PauliStringTS{$Ls}: $Ls is not a Tuple")
    end

    N = qubitlength(p)
    if Base.prod(Ls) != N
        error("Cannot construct PauliStringTS{$Ls} from PauliString{$N}: $(join(Ls, "×")) != $N.")
    end

    rep = find_representative(p, Ls)
    return PauliStringTS{Ls, typeof(rep.v)}(rep.v, rep.w)
end

"""
    representative(p::PauliStringTS) -> PauliString

Returns a unique representative string of the translation symmetric sum of the Pauli string `p`.
"""
representative(p::PauliStringTS{Ls, T}) where {Ls, T} = PauliString{Base.prod(Ls), T}(p.v, p.w)

@inline function find_representative(p::PauliString, Ls)
    # It is crucial for performance that this function, along with shift
    # gets inlined so that Ls get constant-propagated eliminating some integer divisions.
    #
    # This broke when using maximum instead of the loop.
    pmax = p
    for shifts in all_shifts(Ls)
        pshift = shift(p, Ls, shifts)
        if pshift > pmax
            pmax = pshift
        end
    end
    return pmax
end

@inline all_shifts(Ls) = Iterators.product(map(L -> 1:L, Ls)...)
@inline all_shifts(::Type{<:PauliStringTS{Ls}}) where {Ls} = all_shifts(Ls)

"""
    shift(p::PauliString, Ls::Tuple, shifts::Tuple)
    shift(p::Operator, Ls::Tuple, shifts::Tuple)

Interpret a PauliString as a multidimensional array whose size is given by the tuple `Ls` and apply a tuple of periodic `shifts` in the different dimensions to it.
"""
@inline shift(p::PauliString{N, T}, Ls::Tuple, shifts::Tuple) where {N, T} =
    PauliString{N, T}(_shift(p.v, Ls, shifts), _shift(p.w, Ls, shifts))

shift(p::Operator, Ls::Tuple, shifts::Tuple) = Operator(shift.(p.strings, Ref(Ls), Ref(shifts)), copy(p.coeffs))

"""
    shift(p::PauliString, r::Integer)
    shift(p::Operator, r::Integer)

Shift `p` by `r` bits, assuming it lives in 1D.
"""
shift(p::Union{PauliString, Operator}, r::Integer) = shift(p, (qubitlength(p),), (r,))

"""
    _shift(x::Unsigned, Ls::Tuple, shifts::Tuple)

Interpret `x` as a multidimensional array of bits whose size is given by the tuple `Ls` and apply a tuple of periodic `shifts` to it.
"""
@inline function _shift(x::Unsigned, Ls::Tuple, distances::Tuple)
    stride = 1
    for (distance, L) in zip(distances, Ls)
        x = _shift(x, stride * distance, stride * L)
        stride *= L
    end
    return x & ~(~zero(x) << Base.prod(Ls))
end

"""
    _shift(x, shift, stride)

Periodically shift `x` by `s` bits within periodic windows of length `stride`.
"""
@inline function _shift(x::Unsigned, s, stride)
    @assert 0 <= s <= stride
    _, s, stride = unsigned.(promote(x, s, stride))
    rightzeros = 8 * sizeof(x) % stride
    rshift = stride - s

    # masks for the region that crosses the periodic boundary after the shift
    magicmask = (~(~zero(x) << s)) * ((~zero(x) ÷ ((one(x) << stride) - one(x))) >> rightzeros << rshift)

    # mostly we can just shift x by s, but for bits crossing the boundary,
    # we shift them all the way to the other side.
    return (x << s) & (~(magicmask >> rshift)) | ((x & magicmask) >> rshift)
end

const OperatorTS{Ls, U, T} = Operator{PauliStringTS{Ls, U}, T}

function OperatorTS{Ls}(o::Operator) where {Ls}
    periodic_strings = PauliStringTS{Ls}.(o.strings)
    coeffs = copy(o.coeffs)

    return compress(Operator{eltype(periodic_strings), eltype(coeffs)}(periodic_strings, coeffs))
end

qubitsize(::Type{<:OperatorTS{Ls}}) where {Ls} = Ls
qubitsize(op::Operator{<:PauliStringTS}) = qubitsize(typeof(op))

function representative(o::OperatorTS)
    return Operator(representative.(o.strings), copy(o.coeffs))
end

function resum(o::OperatorTS)
    Ls = qubitsize(o)
    rep_op = representative(o)

    op = Operator(similar(rep_op.strings, 0), similar(rep_op.coeffs, 0))
    for s in all_shifts(paulistringtype(o))
        op += shift(rep_op, Ls, s)
    end
    return op
end

function binary_kernel(op, A::Operator{<:PauliStringTS}, B::Operator{<:PauliStringTS}; epsilon::Real = 0, maxlength::Int = 1000)
    checklength(A, B)
    Ls = qubitsize(A)

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
            for s in all_shifts(paulistringtype(A))
                p, k = op(rep1, shift(rep2, Ls, s))
                c = c1 * c2 * k
                if (k != 0) && (abs(c) > epsilon) && pauli_weight(p) < maxlength
                    setwith!(+, d, PauliStringTS{Ls}(p), c)
                end
            end
        end
    end

    return typeof(A)(collect(keys(d)), collect(values(d)))
end

function trace(o::Operator{<:PauliStringTS})
    r = zero(scalartype(o))

    for (p, c) in zip(o.strings, o.coeffs)
        if isone(representative(p))
            r += c
        end
    end
    return r * qubitlength(o) * 2^qubitlength(o)
end

opnorm(o::Operator{<:PauliStringTS}) = sqrt(real(trace_product(dagger(o), o)))


"""
    is_ts(o::Operator)
    
return true if `o` is translation symmetric in one dimension

    is_ts(o::Operator, Ls::Tuple)

return true if `o` is translation symmetric on a hypercube with side lengths `Ls`.
"""
is_ts(o::Operator) = is_ts(o, (qubitlength(o),))
function is_ts(o::Operator, Ls::Tuple)
    for s in all_shifts(Ls)
        if opnorm(o - shift(o, Ls, s)) / opnorm(o) > 1.0e-10
            return false
        end
    end
    return true
end

"""
    is_ts2d(o::Operator, L1)

return true if `o` is translation symmetric on a rectangle with sidelengths `L1` × `qubitlength(o)÷L1`.
"""
is_ts2d(o, L1) = is_ts(o, (L1, qubitlength(o) ÷ L1))


Base.:+(a::Number, o::Operator{<:PauliStringTS}) = OperatorTS{qubitsize(o)}(representative(o) + a / qubitlength(o))
Base.:+(o::Operator{<:PauliStringTS}, a::Number) = a + o

# deprecated

"""
    OperatorTS2D(N::Integer, L1::Integer)

Initialize a zero 2D translation-invariant operator on `N` qubits, with extent of `L1` in the ``a_1`` direction.
"""
OperatorTS2D(N::Integer, L1::Integer) = (N % L1 == 0) ? Operator{periodicpaulistringtype((L1, N ÷ L1)), ComplexF64}() : error("N must be divisible by L1")

function OperatorTS2D(N::Int, L1::Int, v::Vector{T}, w::Vector{T}, coef::AbstractVector) where {T <: Unsigned}
    length(v) == length(w) == length(coef) || error("v, w, and coef must have the same length")
    P = periodicpaulistringtype((L1, N ÷ L1))
    strings = P.(v, w)
    return Operator(strings, coef)
end

OperatorTS2D(pauli::AbstractString, L1::Integer) = Operator{periodicpaulistringtype((L1, length(pauli) ÷ L1)), ComplexF64}(pauli)
function OperatorTS2D(op::Operator, L1::Integer; full = true)
    L2 = qubitlength(op) ÷ L1
    if full && !is_ts(op, (L1, L2))
        error("o is not translation symmetric. If you want to initialize an OperatorTS1D only with its local part H_0, then set full=false")
    end

    full && (op /= qubitlength(op))
    return OperatorTS{(L1, L2)}(op)
end

OperatorTS2D(op::OperatorTS) = typeof(op)(copy(op.strings), copy(op.coeffs))

"""
    OperatorTS1D(N::Integer)

Initialize a zero 1D translation-invariant operator on `N` qubits.
"""
OperatorTS1D(N::Integer) = OperatorTS1D{paulistringtype(N), ComplexF64}()

function OperatorTS1D(N::Int, v::Vector{T}, w::Vector{T}, coef::AbstractVector) where {T <: Unsigned}
    length(v) == length(w) == length(coef) || error("v, w, and coef must have the same length")
    P = periodicpaulistringtype((N,))
    strings = P.(v, w)
    return OperatorTS1D{P, ComplexF64}(strings, coef)
end

"""
    OperatorTS1D(o::Operator; full=true)

Initialize a 1D translation invariant operator from an Operator
\$O=\\sum_i o_i O_i\$ where \$O_i=T_i(O_0)\$ and \$T_i\$ is the i-sites translation operator.
Set full=true if passing \$O\$, an Operator that is supported on the whole chain (i.e converting from a translation symmetric [`Operator`](@ref))
Set full=false if passing \$O_0\$,a local term o such that the full operator is \$O=\\sum_i o_i T_i(O_0)\$
"""
function OperatorTS1D(o::Operator; full = true)
    if full && !is_ts(o)
        error("o is not translation symmetric. If you want to initialize an OperatorTS1D only with its local part H_0, then set full=false")
    end
    full && (o /= qubitlength(o))
    return o3 = OperatorTS{(qubitlength(o),)}(o)
end

function Operator(o::Operator{<:PauliStringTS}; rs = true)
    Base.depwarn("Operator(o::OperatorTS; rs) is deprecated. If rs=false, use representative(o), if rs=true, use resum(o)", :Operator)
    rs && (o = resum(o))
    return Operator(copy(o.strings), copy(o.coeffs))
end


@deprecate rotate_lower(p, r) shift(p, r)

@deprecate rotate(p, r) shift(p, r)

@deprecate shift_left(p::PauliString) find_representative(p, (qubitlength(p),))
@deprecate shift_left(O::Operator) representative(OperatorTS{(qubitlength(O),)}(O))
@deprecate shift1(O::Operator) representative(OperatorTS{(qubitlength(O),)}(O))
@deprecate shift_origin(O::Operator, L1) representative(OperatorTS{(L1, qubitlength(O) ÷ L1)}(O))
