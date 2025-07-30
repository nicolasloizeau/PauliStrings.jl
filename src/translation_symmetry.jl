"""
    PauliStringTS{Ls, T<:Unsigned} <: AbstractPauliString

Type representing a translation symmetric sum of a Pauli string. The tuple `Ls` specifies the period in each dimension. The sum is normalized by `1/prod(Ls)`.
"""
struct PauliStringTS{Ls, T <: Unsigned} <: AbstractPauliString
    v::T
    w::T
end

paulistringtype(Ls::NTuple{<:Any, Integer}) = PauliStringTS{Ls, uinttype(Base.prod(Ls))}
qubitlength(::Type{<:PauliStringTS{Ls}}) where {Ls} = Base.prod(Ls)

"""
    qubitsize(::Type{<:PauliStringTS})
    qubitsize(p::PauliStringTS)

Get the tuple of periods of a given `PauliStringTS`.
"""
qubitsize(::Type{<:PauliStringTS{Ls}}) where {Ls} = Ls
qubitsize(p::PauliStringTS) = qubitsize(p)

for count in [:xcount, :ycount, :zcount]
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
    N = qubitlength(p)
    if Base.prod(Ls) != N
        error("Cannot construct PauliString{$Ls} from PauliString{$N}: $(join(Ls, "×")) != $N.")
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
    pmax = p
    for shifts in Iterators.product(map(L -> 1:L, Ls)...)
        pshift = shift(p, Ls, shifts)
        if pshift > pmax
            pmax = pshift
        end
    end
    return pmax
end

"""
    shift(p::PauliString, Ls, shifts)

Interpret a PauliString as a multidimensional array whose size is given by the tuple `Ls` and apply a tuple of periodic `shifts` in the different dimensions to it.
"""
@inline shift(p::PauliString{N, T}, Ls, shifts) where {N, T} =
    PauliString{N, T}(shift(p.v, Ls, shifts), shift(p.w, Ls, shifts))

"""
    shift(x::Unsigned, Ls::Tuple, shifts::Tuple)

Interpret `x` as a multidimensional array of bits whose size is given by the tuple `Ls` and apply a tuple of periodic `shifts` to it.
"""
@inline function shift(x::Unsigned, Ls::Tuple, distances::Tuple)
    stride = 1
    for (distance, L) in zip(distances, Ls)
        x = shift(x, stride * distance, stride * L)
        stride *= L
    end
    return x & ~(~zero(x) << Base.prod(Ls))
end

"""
    shift(x, shift, stride)

Periodically shift `x` by `s` bits within periodic windows of length `stride`.
"""
@inline function shift(x::Unsigned, s, stride)
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
