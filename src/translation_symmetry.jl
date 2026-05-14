"""
    PauliStringTS{Ls, Ps, T<:Unsigned} <: AbstractPauliString

Type representing a translation symmetric sum of a Pauli string. The tuple `Ls` specifies the period in each dimension.
The tuple `Ps` specifies which dimensions are periodic (true) or open (false). By default all dimensions are periodic.
"""
struct PauliStringTS{Ls,Ps,T<:Unsigned} <: AbstractPauliString
    v::T
    w::T
end

periodicpaulistringtype(Ls::NTuple{<:Any,Integer}) = PauliStringTS{Ls,ntuple(i -> true, length(Ls)),uinttype(Base.prod(Ls))}
periodicpaulistringtype(Ls::NTuple{<:Any,Integer}, Ps::NTuple{<:Any,Bool}) = PauliStringTS{Ls,Ps,uinttype(Base.prod(Ls))}
qubitlength(::Type{<:PauliStringTS{Ls}}) where {Ls} = Base.prod(Ls)

"""
    qubitsize(::Type{<:PauliStringTS})
    qubitsize(p::PauliStringTS)

Get the tuple of periods of a given `PauliStringTS`.
"""
qubitsize(::Type{<:PauliStringTS{Ls}}) where {Ls} = Ls
qubitsize(p::PauliStringTS) = qubitsize(typeof(p))

"""
    periodicflags(::Type{<:PauliStringTS})
    periodicflags(p::PauliStringTS)

Get the tuple of periodic flags of a given `PauliStringTS`.
"""
periodicflags(::Type{<:PauliStringTS{Ls,Ps}}) where {Ls,Ps} = Ps
periodicflags(p::PauliStringTS) = periodicflags(typeof(p))

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

        segments = str.(eachslice(slice, dims=((3:ndims(slice))...,)))

        return join(["[:, :, $(join(Tuple(I), ", "))] =\n" * segments[I] for I in CartesianIndices(segments)], "\n\n")
    end
    σij = reshape(collect(vw_to_string(p.v, p.w, Base.prod(Ls))[1]), Ls)
    return str(σij)
end

Base.one(::Type{<:PauliStringTS{Ls,Ps,T}}) where {Ls,Ps,T} = PauliStringTS{Ls,Ps}(one(PauliString{Base.prod(Ls),T}))
Base.isless(a::PauliStringTS, b::PauliStringTS) = isless(representative(a), representative(b))

"""
    PauliStringTS{Ls}(p::PauliString)
    PauliStringTS{Ls,Ps}(p::PauliString)

Construct a translation symmetric sum of a Pauli string `p`. `Ls` is a tuple that specifies the periodicity in each dimension.
`Ps` is an optional tuple of Bool values indicating which dimensions are periodic (default: all true).
"""
function PauliStringTS{Ls}(p::PauliString) where {Ls}
    return PauliStringTS{Ls,ntuple(i -> true, length(Ls))}(p)
end

function PauliStringTS{Ls,Ps}(p::PauliString) where {Ls,Ps}
    if !(Ls isa Tuple)
        error("Cannot construct PauliStringTS{$Ls}: $Ls is not a Tuple")
    end
    if !(Ps isa Tuple)
        error("Cannot construct PauliStringTS{$Ls,$Ps}: $Ps is not a Tuple")
    end
    if length(Ls) != length(Ps)
        error("Cannot construct PauliStringTS{$Ls,$Ps}: Ls and Ps must have the same length")
    end

    N = qubitlength(p)
    if Base.prod(Ls) != N
        error("Cannot construct PauliStringTS{$Ls,$Ps} from PauliString{$N}: $(join(Ls, "×")) != $N.")
    end

    rep = find_representative(p, Ls, Ps)
    return PauliStringTS{Ls,Ps,typeof(rep.v)}(rep.v, rep.w)
end

PauliStringTS{Ls}(pauli::AbstractString) where {Ls} = PauliStringTS{Ls}(PauliString(pauli))


"""
    representative(p::PauliStringTS)

Returns a unique representative string of the translation symmetric sum of the Pauli string `p`.
"""
representative(p::PauliStringTS{Ls,Ps,T}) where {Ls,Ps,T} = PauliString{Base.prod(Ls),T}(p.v, p.w)

@inline function find_representative(p::PauliString, Ls)
    return find_representative(p, Ls, ntuple(i -> true, length(Ls)))
end

@inline function find_representative(p::PauliString, Ls, Ps)
    # It is crucial for performance that this function, along with shift
    # gets inlined so that Ls get constant-propagated eliminating some integer divisions.
    #
    # This broke when using maximum instead of the loop.
    pmax = p
    for shifts in all_shifts(Ls, Ps)
        pshift = shift(p, Ls, Ps, shifts)
        if pshift > pmax
            pmax = pshift
        end
    end
    return pmax
end

@inline all_shifts(Ls) = all_shifts(Ls, ntuple(i -> true, length(Ls)))
@inline function all_shifts(Ls, Ps)
    # For each dimension, if periodic, generate all shifts 1:L, otherwise just 1 (identity)
    return Iterators.product(map((L, p) -> p ? (1:L) : (1:1), Ls, Ps)...)
end
@inline all_shifts(::Type{<:PauliStringTS{Ls,Ps}}) where {Ls,Ps} = all_shifts(Ls, Ps)

"""
    shift(p::PauliString, Ls::Tuple, shifts::Tuple)
    shift(p::PauliString, Ls::Tuple, Ps::Tuple, shifts::Tuple)
    shift(p::Operator, Ls::Tuple, shifts::Tuple)
    shift(p::Operator, Ls::Tuple, Ps::Tuple, shifts::Tuple)

Interpret a PauliString as a multidimensional array whose size is given by the tuple `Ls` and apply a tuple of periodic `shifts` in the different dimensions to it.
`Ps` is an optional tuple of Bool values indicating which dimensions are periodic (default: all true).
"""
@inline shift(p::PauliString{N,T}, Ls::Tuple, shifts::Tuple) where {N,T} =
    shift(p, Ls, ntuple(i -> true, length(Ls)), shifts)

@inline shift(p::PauliString{N,T}, Ls::Tuple, Ps::Tuple, shifts::Tuple) where {N,T} =
    PauliString{N,T}(_shift(p.v, Ls, Ps, shifts), _shift(p.w, Ls, Ps, shifts))

shift(p::Operator, Ls::Tuple, shifts::Tuple) = shift(p, Ls, ntuple(i -> true, length(Ls)), shifts)
shift(p::Operator, Ls::Tuple, Ps::Tuple, shifts::Tuple) = Operator(shift.(p.strings, Ref(Ls), Ref(Ps), Ref(shifts)), copy(p.coeffs))

"""
    shift(p::PauliString, r::Integer)
    shift(p::Operator, r::Integer)

Shift `p` by `r` bits, assuming it lives in 1D.
"""
shift(p::Union{PauliString,Operator}, r::Integer) = shift(p, (qubitlength(p),), (r,))

"""
    _shift(x::Unsigned, Ls::Tuple, shifts::Tuple)
    _shift(x::Unsigned, Ls::Tuple, Ps::Tuple, shifts::Tuple)

Interpret `x` as a multidimensional array of bits whose size is given by the tuple `Ls` and apply a tuple of periodic `shifts` to it.
`Ps` is an optional tuple of Bool values indicating which dimensions are periodic (default: all true).
"""
@inline _shift(x::Unsigned, Ls::Tuple, distances::Tuple) = _shift(x, Ls, ntuple(i -> true, length(Ls)), distances)

@inline function _shift(x::Unsigned, Ls::Tuple, Ps::Tuple, distances::Tuple)
    stride = 1
    for (distance, L, p) in zip(distances, Ls, Ps)
        if p  # Only shift if this dimension is periodic
            x = _shift(x, stride * distance, stride * L)
        end
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
    stride_int = stride
    _, s, stride = unsigned.(promote(x, s, stride))
    rightzeros = 8 * sizeof(x) % stride_int
    rshift = stride - s

    # masks for the region that crosses the periodic boundary after the shift
    magicmask = (~(~zero(x) << s)) * ((~zero(x) ÷ ((one(x) << stride) - one(x))) >> rightzeros << rshift)

    # mostly we can just shift x by s, but for bits crossing the boundary,
    # we shift them all the way to the other side.
    return (x << s) & (~(magicmask >> rshift)) | ((x & magicmask) >> rshift)
end

const OperatorTS{Ls,Ps,U,T} = Operator{PauliStringTS{Ls,Ps,U},T}

@doc raw"""
    OperatorTS{Ls}(o::Operator)
    OperatorTS{Ls,Ps}(o::Operator)

Construct an ``n``-dimensional translation symmetric operator from `o` where `Ls` is a tuple of integers `(L1, L2, ...)`
and `Ps` is an optional tuple of Bool values indicating which dimensions are periodic (default: all true).
The resulting operator is equivalent to

```math
O_\mathrm{TS} = \sum_T T^\dag O T
```
where ``T`` are all translations on the `L1`×`L2`×… hypercube. So if you feed it an operator that is already a sum, you should afterwards normalize it by the number of sites.

To get a dense operator from this lazy sum representation, see [`resum`](@ref). To get a single term, see [`representative`](@ref).
"""
function OperatorTS{Ls}(o::Operator) where {Ls}
    return OperatorTS{Ls,ntuple(i -> true, length(Ls))}(o)
end

function OperatorTS{Ls,Ps}(o::Operator) where {Ls,Ps}
    periodic_strings = PauliStringTS{Ls,Ps}.(o.strings)
    coeffs = copy(o.coeffs)

    return compress(Operator{eltype(periodic_strings),eltype(coeffs)}(periodic_strings, coeffs))
end


OperatorTS{Ls}(pauli::PauliString) where {Ls} = OperatorTS{Ls}(Operator(pauli))
function OperatorTS(pauli::PauliStringTS)
    periodic_strings = [pauli]
    coeffs = [1.0im^ycount(pauli)]
    return compress(Operator{eltype(periodic_strings),eltype(coeffs)}(periodic_strings, coeffs))
end
OperatorTS{Ls}(pauli::AbstractString) where {Ls} = OperatorTS{Ls}(PauliString(pauli))


qubitsize(::Type{<:OperatorTS{Ls}}) where {Ls} = Ls
qubitsize(op::Operator{<:PauliStringTS}) = qubitsize(typeof(op))

periodicflags(::Type{<:OperatorTS{Ls,Ps}}) where {Ls,Ps} = Ps
periodicflags(op::Operator{<:PauliStringTS}) = periodicflags(typeof(op))

"""
    representative(o::OperatorTS)

Returns a unique term of the symmetric sum represented by `o`.
"""
representative(o::OperatorTS) = Operator(representative.(o.strings), copy(o.coeffs))


"""
    resum(o::OperatorTS)

Perform the symmetric sum represented by `o` to yield a dense `Operator` containing unsymmetrized PauliStrings.
"""
function resum(o::OperatorTS)
    Ls = qubitsize(o)
    Ps = periodicflags(o)
    rep_op = representative(o)

    op = Operator(similar(rep_op.strings, 0), similar(rep_op.coeffs, 0))
    for s in all_shifts(paulistringtype(o))
        op += shift(rep_op, Ls, Ps, s)
    end
    return op
end

function binary_kernel(op, A::Operator{<:PauliStringTS}, B::Operator{<:PauliStringTS}; epsilon::Real=0, maxlength::Int=1000)
    checklength(A, B)
    Ls = qubitsize(A)
    Ps = periodicflags(A)

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
                p, k = op(rep1, shift(rep2, Ls, Ps, s))
                c = c1 * c2 * k
                if (k != 0) && pauli_weight(p) < maxlength
                    setwith!(+, d, PauliStringTS{Ls,Ps}(p), c)
                end
            end
        end
    end

    o = typeof(A)(collect(keys(d)), collect(values(d)))
    return (eltype(o.coeffs) == ComplexF64) ? cutoff(o, epsilon) : o
end

function trace(o::Operator{<:PauliStringTS})
    r = zero(scalartype(o))

    for (p, c) in zip(o.strings, o.coeffs)
        if isone(representative(p))
            r += c
        end
    end
    # Calculate the number of translations: product of lengths for periodic dimensions only
    Ls = qubitsize(o)
    Ps = periodicflags(o)
    num_translations = Base.prod(L for (L, p) in zip(Ls, Ps) if p)
    return r * num_translations * 2^qubitlength(o)
end

Base.@deprecate opnorm(o::Operator{<:PauliStringTS}) LinearAlgebra.norm(o::Operator{<:PauliStringTS})

function LinearAlgebra.norm(o::Operator{<:PauliStringTS}; normalize=false)
    n = sqrt(real(trace_product(o', o)))
    if normalize
        return n / (2.0^(qubitlength(o) / 2))
    else
        return n
    end
end

function LinearAlgebra.norm(p::PauliStringTS; normalize=false)
    normalize && return sqrt(qubitlength(p))
    return sqrt(2.0^qubitlength(p) * qubitlength(p))
end


"""
    is_ts(o::Operator)

return true if `o` is translation symmetric in one dimension
"""
is_ts(o::Operator) = is_ts(o, (qubitlength(o),), (true,))
"""
    is_ts(o::Operator, Ls::Tuple)
    is_ts(o::Operator, Ls::Tuple, Ps::Tuple)

return true if `o` is translation symmetric on a hypercube with side lengths `Ls`.
`Ps` is an optional tuple of Bool values indicating which dimensions are periodic (default: all true).
"""
is_ts(o::Operator, Ls::Tuple) = is_ts(o, Ls, ntuple(i -> true, length(Ls)))
function is_ts(o::Operator, Ls::Tuple, Ps::Tuple)
    for s in all_shifts(Ls, Ps)
        if norm(o - shift(o, Ls, Ps, s)) / norm(o) > 1.0e-10
            return false
        end
    end
    return true
end

"""
    is_ts2d(o::Operator, L1)
    is_ts2d(o::Operator, L1, Ps::NTuple{2,Bool})

return true if `o` is translation symmetric on a rectangle with sidelengths `L1` × `qubitlength(o)÷L1`.
`Ps` is an optional tuple of Bool values indicating which dimensions are periodic (default: both true).
"""
is_ts2d(o, L1) = is_ts(o, (L1, qubitlength(o) ÷ L1), (true, true))
is_ts2d(o, L1, Ps::NTuple{2,Bool}) = is_ts(o, (L1, qubitlength(o) ÷ L1), Ps)


Base.:+(a::Number, o::Operator{<:PauliStringTS}) = begin
    # Calculate the number of translations: product of lengths for periodic dimensions only
    Ls = qubitsize(o)
    Ps = periodicflags(o)
    num_translations = Base.prod(L for (L, p) in zip(Ls, Ps) if p)
    OperatorTS{qubitsize(o),periodicflags(o)}(representative(o) + a / num_translations)
end
Base.:+(o::Operator{<:PauliStringTS}, a::Number) = a + o

# deprecated

"""
    OperatorTS2D(N::Integer, L1::Integer)

Initialize a zero 2D translation-invariant operator on `N` qubits, with extent of `L1` in the ``a_1`` direction.
"""
OperatorTS2D(N::Integer, L1::Integer) = (N % L1 == 0) ? Operator{periodicpaulistringtype((L1, N ÷ L1)),ComplexF64}() : error("N must be divisible by L1")

function OperatorTS2D(N::Int, L1::Int, v::Vector{T}, w::Vector{T}, coef::AbstractVector) where {T<:Unsigned}
    length(v) == length(w) == length(coef) || error("v, w, and coef must have the same length")
    P = periodicpaulistringtype((L1, N ÷ L1))
    strings = P.(v, w)
    return Operator(strings, coef)
end

OperatorTS2D(pauli::AbstractString, L1::Integer) = Operator{periodicpaulistringtype((L1, length(pauli) ÷ L1)),ComplexF64}(pauli)
function OperatorTS2D(op::Operator, L1::Integer; full=true, periodic::NTuple{2,Bool}=(true, true))
    L2 = qubitlength(op) ÷ L1
    if full && !is_ts(op, (L1, L2), periodic)
        error("o is not translation symmetric. If you want to initialize an OperatorTS1D only with its local part H_0, then set full=false")
    end

    if full
        # Divide by the number of translations: product of lengths for periodic dimensions only
        num_translations = Base.prod(L for (L, p) in zip((L1, L2), periodic) if p)
        op /= num_translations
    end
    return OperatorTS{(L1, L2),periodic}(op)
end

OperatorTS2D(op::OperatorTS) = typeof(op)(copy(op.strings), copy(op.coeffs))


"""
Initialize a zero 1D translation-invariant operator on `N` qubits.
"""
OperatorTS1D(N::Integer) = OperatorTS1D{paulistringtype(N),ComplexF64}()

function OperatorTS1D(N::Int, v::Vector{T}, w::Vector{T}, coef::AbstractVector) where {T<:Unsigned}
    length(v) == length(w) == length(coef) || error("v, w, and coef must have the same length")
    P = periodicpaulistringtype((N,))
    strings = P.(v, w)
    return OperatorTS1D{P,ComplexF64}(strings, coef)
end

"""
    OperatorTS1D(o::Operator; full=true)

Initialize a 1D translation invariant operator from an Operator
\$O=\\sum_i o_i O_i\$ where \$O_i=T_i(O_0)\$ and \$T_i\$ is the i-sites translation operator.
Set full=true if passing \$O\$, an Operator that is supported on the whole chain (i.e converting from a translation symmetric [`Operator`](@ref))
Set full=false if passing \$O_0\$,a local term o such that the full operator is \$O=\\sum_i o_i T_i(O_0)\$
"""
function OperatorTS1D(o::Operator; full=true, periodic::Bool=true)
    if full && !is_ts(o, (qubitlength(o),), (periodic,))
        error("o is not translation symmetric. If you want to initialize an OperatorTS1D only with its local part H_0, then set full=false")
    end
    if full
        # Divide by the number of translations: product of lengths for periodic dimensions only
        # For 1D: if periodic, divide by qubitlength(o), else divide by 1 (no division)
        num_translations = periodic ? qubitlength(o) : 1
        o /= num_translations
    end
    return OperatorTS{(qubitlength(o),),(periodic,)}(o)
end

function Operator(o::Operator{<:PauliStringTS}; rs=true)
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
