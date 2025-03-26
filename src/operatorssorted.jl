# Single Site
# -----------
"""
    @enum[UInt8] PauliBasis I Z X Y

Elementary single-site pauli operators, ordered for easy conversion to and from
the bit representation: `PauliBasis(v, w) = v + 2w`
"""
@enum PauliBasis::UInt8 I = 0b0000 Z = 0b0001 X = 0b0010 Y = 0b0011

function Base.convert(::Type{PauliBasis}, c::Char)
    return (c == '1' || c == 'I') ? I : c == 'X' ? X : c == 'Y' ? Y : c == 'Z' ? Z :
           throw(DomainError(c, "PauliBasis must be one of `'1'`, `'I'`, `'X'`, `'Y'`, `'Z'`"))
end

# @enum ExtendedPauliBasis I X Y Z Sx Sy Sz S⁻ S⁺

# Single String
# -------------
"""
    PauliString{N,B}

Single string of Pauli operators, of length `N` with two bitstrings `v` and `w`.
"""
struct PauliString{N,B} <: AbstractVector{PauliBasis}
    v::B
    w::B
end
PauliString{N}(v::B, w::B) where {N,B} = PauliString{N,B}(v, w)

storagetype(N::Int) = storagetype(Val(N))
function storagetype(::Val{N}) where {N}
    return N < 1 ? throw(DomainError(N, "String length must be strictly positive")) :
           N ≤ 32 ? UInt32 : N ≤ 64 ? UInt64 : N ≤ 128 ? UInt128 : BitVector
end

Base.length(::Type{<:PauliString{N}}) where {N} = N

PauliString(x::Union{String,Vector}) = PauliString{length(x)}(x)
PauliString{N}(x::Union{String,Vector}) where {N} = PauliString{N,storagetype(Val(N))}(x)

# enable syntax p"IXXI"
macro p_str(x)
    return PauliString(x)
end

function Base.show(io::IO, x::PauliString)
    print(io, "p\"")
    join(io, x)
    print(io, "\"")
end
function Base.show(io::IO, ::MIME"text/plain", x::PauliString)
    print(io, "p\"")
    join(io, x)
    print(io, "\"")
end

function PauliString{N,B}(x::Union{String,Vector}) where {N,B}
    length(x) == N || throw(ArgumentError("Invalid input length"))
    p = one(PauliString{N,B})
    @inbounds for i in eachindex(x)
        p = Base.setindex(p, x[i], i)
    end
    return p
end

Base.sign(p::PauliString) = prod(x -> x == Y ? im : 1, p; init=complex(1))

Base.one(::Type{PauliString{N,B}}) where {N,B} = PauliString{N,B}(zero(B), zero(B))

Base.size(::PauliString{N}) where {N} = (N,)

@inline function Base.getindex(p::PauliString, i::Int)
    @boundscheck checkbounds(p, i)
    vi = ((p.v >> (i - 1)) & 1) % UInt8
    wi = ((p.w >> (i - 2)) & 2) % UInt8
    return PauliBasis(vi | wi)
end

@inline function Base.setindex(p::PauliString, b::PauliBasis, i::Int)
    @boundscheck checkbounds(p, i)
    v = p.v | ((UInt8(b) & 1) << (i - 1))
    w = p.w | ((UInt8(b) & 2) << (i - 2))
    return typeof(p)(v, w)
end
@inline Base.setindex(p::PauliString, b::Char, i::Int) = Base.setindex(p, convert(PauliBasis, b), i)

Base.xor(p1::P, p2::P) where {P<:PauliString} = P(p1.v ⊻ p2.v, p1.w ⊻ p2.w)
Base.:*(p1::P, p2::P) where {P} = p1 ⊻ p2, 1 - ((count_ones(p1.v & p2.w) & 1) << 1)
Base.isequal(x::P, y::P) where {P<:PauliString} = x.v == y.v && x.w == y.w
Base.:(==)(x::P, y::P) where {P<:PauliString} = isequal(x, y)

@inline Base.hash(x::P, h::UInt) where {P<:PauliString} = hash((x.v, x.w), h)
@inline Base.hash(x::PauliString{N,UInt32}, h::UInt) where {N} = hash(UInt(x.v) << sizeof(UInt) + UInt(x.w), h)

function commutator(p1::P, p2::P) where {P<:PauliString}
    p = p1 ⊻ p2
    k = ((count_ones(p2.v & p1.w) & 1) << 1) - ((count_ones(p1.v & p2.w) & 1) << 1)
    # k2 = -((count_ones((p2.v & p1.w) ⊻ (p1.v & p2.w)) & 1) << 1)
    # @show k2
    return p, k
end

function anticommutator(p1::P, p2::P) where {P<:PauliString}
    p = p1 ⊻ p2
    k = 2 - (((count_ones(p1.v & p2.w) & 1) << 1) + ((count_ones(p1.w & p2.v) & 1) << 1))
    return p, k
end

pauli_weight(p::PauliString) = count_ones(p.v | p.w)

Base.:+(x::P, y::P, z::P...) where {P<:PauliString} = OperatorSorted(collect((x, y, z...)))

# Sum of strings
# --------------
"""
    OperatorSorted

New datastructure for linear combinations of Pauli strings, kept in sorted order.
"""
struct OperatorSorted{P<:PauliString,C<:Complex} <: Operator
    paulistrings::Vector{Pair{P,C}}
    OperatorSorted{P,C}(x::Vector{Pair{P,C}}) where {P<:PauliString,C<:Complex} = new{P,C}(x)
end

function OperatorSorted(x::Vector{Pair{P,C}}) where {P<:PauliString,C<:Complex}
    return OperatorSorted{P,C}(x)
end

OperatorSorted(x::Vector{P}) where {P<:PauliString} = OperatorSorted{P,Complex{Bool}}(x)
function OperatorSorted{P,C}(x::Vector{P}) where {P<:PauliString,C<:Number}
    return OperatorSorted{P,C}(map(p -> (p => one(C)), x))
end

function OperatorSorted(N::Int, vs::Vector{T}, ws::Vector{T}, coeff::Vector{C}) where {T<:Unsigned,C<:Complex}
    return OperatorSorted(map(vs, ws, coeff) do v, w, c
        return PauliString{N,T}(v, w) => c
    end)
end
function OperatorSorted(N::Int)
    B = storagetype(N)
    P = PauliString{N,B}
    return OperatorSorted{P,Complex{Bool}}(Vector{Pair{P,Complex{Bool}}}(undef, 0))
end

function Base.show(io::IO, o::OperatorSorted)
    for (p, c) in o
        pstr = string(p)
        cstr = round(c / sign(p), digits=10)
        println(io, "($cstr) $pstr")
    end
end

stringtype(O::OperatorSorted) = stringtype(typeof(O))
stringtype(::Type{OperatorSorted{P,C}}) where {P,C} = P
scalartype(O::OperatorSorted) = scalartype(typeof(O))
scalartype(::Type{OperatorSorted{P,C}}) where {P,C} = C

Base.length(o::OperatorSorted) = length(o.paulistrings)
Base.zero(o::OperatorSorted) = OperatorSorted(similar(o.paulistrings, 0))
Base.zero(::Type{OperatorSorted{P,C}}) where {P,C} = OperatorSorted{P,C}(Pair{P,C}[])
Base.one(o::OperatorSorted) = OperatorSorted([one(stringtype(o)) => one(scalartype(o))])

function Base.similar(::OperatorSorted{P,C}, ::Type{T}) where {P,C,T<:Complex}
    return zero(OperatorSorted{P,T})
end
Base.copy(o::OperatorSorted) = OperatorSorted(copy(o.paulistrings))

function Base.:(==)(o1::OperatorSorted, o2::OperatorSorted)
    return o1.paulistrings == o2.paulistrings
end
function Base.isapprox(o1::OperatorSorted, o2::OperatorSorted; kwargs...)
    return first.(o1.paulistrings) == first.(o2.paulistrings) && isapprox(last.(o1.paulistrings), last.(o2.paulistrings); kwargs...)
end

@inline Base.getindex(o::OperatorSorted, i::Int) = o.paulistrings[i]
@inline function Base.getindex(o::OperatorSorted, p::PauliString)
    length(p) == length(stringtype(o)) || error()
    i = searchsortedfirst(o.paulistrings, p; by=first)
    return i > length(o) ? zero(scalartype(o)) : last(@inbounds o[i])
end

@inline Base.iterate(o::OperatorSorted, args...) = iterate(o.paulistrings, args...)

function Base.merge!(o::OperatorSorted)
    m = length(o)
    m == 0 && return o

    paulistrings = sort!(o.paulistrings; by=first)

    current = processed = 1
    @inbounds current_str, current_val = paulistrings[current]
    @inbounds while current < m
        current += 1
        next_str, next_val = paulistrings[current]
        if current_str == next_str
            current_val = paulistrings[processed].second
            paulistrings[processed] = current_str => (current_val + next_val)
        else
            processed += 1
            current_str = next_str
            current_val = next_val
            if processed < current
                paulistrings[processed] = current_str => current_val
            end
        end
    end
    processed < m && resize!(paulistrings, processed)

    return o
end
function Base.merge(o1::OperatorSorted, o2::OperatorSorted)
    paulistrings = merge_sorted(o1.paulistrings, o2.paulistrings; by=first)
    m = length(paulistrings)
    current = processed = 1
    @inbounds current_str, current_val = paulistrings[current]
    @inbounds while current < m
        current += 1
        next_str, next_val = paulistrings[current]
        if current_str == next_str
            current_val = paulistrings[processed].second
            paulistrings[processed] = current_str => (current_val + next_val)
        else
            processed += 1
            current_str = next_str
            current_val = next_val
            if processed < current
                paulistrings[processed] = current_str => current_val
            end
        end
    end
    processed < m && resize!(paulistrings, processed)

    return OperatorSorted(paulistrings)
end

function Base.merge!(o::OperatorSorted, os::OperatorSorted...)
    for o2 in os
        append!(o.paulistrings, o2.paulistrings)
    end
    return merge!(o)
end

function Base.:+(o1::O, o2::O) where {O<:OperatorSorted}
    paulistrings = O(vcat(o1.paulistrings, o2.paulistrings))
    return merge!(paulistrings)
end

Base.:-(o::OperatorSorted) = -one(scalartype(o)) * o
Base.:-(o1::O, o2::O) where {O<:OperatorSorted} = o1 + (-o2)

function Base.:*(o::OperatorSorted, λ::Number)
    return OperatorSorted(map(o.paulistrings) do (str, val)
        return str => (val * λ)
    end)
end
Base.:*(λ::Number, o::OperatorSorted) = o * λ
Base.:/(o::OperatorSorted, λ::Number) = o * inv(λ)
Base.:\(λ::Number, o::OperatorSorted) = o * inv(λ)

const MAX_BUFFER = 2^18

function Base.:*(o1::O, o2::O) where {O<:OperatorSorted}
    T = Base.promote_op(*, scalartype(o1), scalartype(o2))
    result = similar(o1, T)
    buffer = similar(o1, T)
    for (p1, c1) in o1, (p2, c2) in o2
        p, k = p1 * p2
        if k != 0
            push!(buffer.paulistrings, p => k * c1 * c2)
            if length(buffer) > MAX_BUFFER
                # @debug "temporary merge"
                result = merge(result, merge!(buffer))
                resize!(buffer.paulistrings, 0)
            end
        end
    end
    return merge(result, merge!(buffer))
end

for f in (:commutator, :anticommutator)
    @eval function $f(o1::O, o2::O; maxlength::Int=1000, epsilon::Real=0) where {O<:OperatorSorted}
        result = zero(o1)
        buffer = zero(o1)
        for (p1, c1) in o1, (p2, c2) in o2
            p, k = $f(p1, p2)
            c = k * c1 * c2
            (k != 0 && pauli_weight(p) < maxlength && abs(c) > epsilon) || continue
            push!(buffer.paulistrings, p => c)
            if length(buffer) > MAX_BUFFER
                # @debug "temporary merge"
                result = merge(result, merge!(buffer))
                resize!(buffer.paulistrings, 0)
            end
        end
        return merge(result, merge!(buffer))
    end
end

function opnorm(o::OperatorSorted)
    return norm(last.(o.paulistrings)) / (2.0^(length(stringtype(o)) / 2))
end

function norm_lanczos(o::OperatorSorted)
    return norm(last.(o.paulistrings))
end

function trim(o::OperatorSorted, max_strings::Int; keepnorm::Bool=false, keep::Operator=zero(o))
    length(o) ≤ max_strings && return copy(o)
    I = partialsortperm(o.paulistrings, 1:max_strings; by=(abs ∘ last), rev=true)

    @assert length(keep) == 0 "TBA"

    o′ = OperatorSorted(sort!(o.paulistrings[I]; by=first))
    return keepnorm ? o′ * (opnorm(o) / opnorm(o′)) : o′
end

# Unsorted version
# ----------------
"""
    OperatorUnsorted

New datastructure for linear combinations of Pauli strings, kept in unsorted order.
"""
struct OperatorUnsorted{P<:PauliString,C<:Complex} <: Operator
    paulistrings::Dictionary{P,C}
    function OperatorUnsorted{P,C}(x::Dictionary{P,C}) where {P<:PauliString,C<:Complex}
        return new{P,C}(x)
    end
end

function OperatorUnsorted(x::Dictionary{P,C}) where {P<:PauliString,C<:Complex}
    return OperatorUnsorted{P,C}(x)
end

function OperatorUnsorted(x::Vector{Pair{P,C}}) where {P<:PauliString,C<:Complex}
    dict = Dictionary{P,C}(first.(x), last.(x))
    return OperatorUnsorted(dict)
end

OperatorUnsorted(x::Vector{P}) where {P<:PauliString} = OperatorUnsorted{P,Complex{Bool}}(x)
function OperatorUnsorted{P,C}(x::Vector{P}) where {P<:PauliString,C<:Number}
    return OperatorUnsorted(map(p -> (p => one(C)), x))
end

function OperatorUnsorted(N::Int)
    P = PauliString{N,storagetype(N)}
    return OperatorUnsorted(Dictionary{P,ComplexF64}())
end

# Construction Utility
# --------------------
Base.:+(o::OperatorUnsorted, a::Number) = o + a * one(o)
Base.:+(a::Number, o::OperatorUnsorted) = o + a * one(o)
Base.:-(o::OperatorUnsorted, a::Number) = o + (-a) * one(o)
Base.:-(a::Number, o::OperatorUnsorted) = a * one(o) - o

function Base.:+(o::OperatorUnsorted, args::Tuple{Number,Vararg{Any}})
    term = one(o)
    c = args[1]
    for i in 2:2:length(args)
        symbol = args[i]::String
        site = args[i+1]::Int
        if length(symbol) == 1
            b = convert(PauliBasis, symbol[1])
            p = Base.setindex(one(stringtype(o)), b, site)
            term *= p
        else
            o2 = zero(o)
            if occursin(symbol, "SxSySz")
                o2 += 0.5, only(uppercase(symbol[2])), site
            elseif symbol == "S+" || symbol == "S-"
                o2 += 0.5, 'X', site
                o2 += (symbol == "S+" ? 0.5im : -0.5im), 'Y', site
            else
                error("Allowed operators: X,Y,Z,Sx,Sy,Sz,S-,S+")
            end
            term *= o2
        end
    end
    return o + c * term
end

function Base.show(io::IO, o::OperatorUnsorted)
    for (p, c) in o
        pstr = string(p)
        cstr = round(c / sign(p), digits=10)
        println(io, "($cstr) $pstr")
    end
end

stringtype(O::OperatorUnsorted) = stringtype(typeof(O))
stringtype(::Type{OperatorUnsorted{P,C}}) where {P,C} = P
scalartype(O::OperatorUnsorted) = scalartype(typeof(O))
scalartype(::Type{OperatorUnsorted{P,C}}) where {P,C} = C

Base.length(o::OperatorUnsorted) = length(o.paulistrings)
Base.zero(o::OperatorUnsorted) = zero(typeof(o))
Base.zero(::Type{OperatorUnsorted{P,C}}) where {P,C} = OperatorUnsorted{P,C}(Dictionary{P,C}())
Base.one(o::OperatorUnsorted) = OperatorUnsorted([one(stringtype(o)) => one(scalartype(o))])

function Base.similar(::OperatorUnsorted{P,C}, ::Type{T}) where {P,C,T<:Complex}
    return zero(OperatorUnsorted{P,T})
end
Base.copy(o::OperatorUnsorted) = OperatorUnsorted(copy(o.paulistrings))

function Base.:(==)(o1::OperatorUnsorted, o2::OperatorUnsorted)
    return isdictequal(o1.paulistrings, o2.paulistrings)
end
function Base.isapprox(o1::OperatorUnsorted, o2::OperatorUnsorted; kwargs...)
    return first.(o1.paulistrings) == first.(o2.paulistrings) && isapprox(last.(o1.paulistrings), last.(o2.paulistrings); kwargs...)
end

# @inline Base.getindex(o::OperatorSorted, i::Int) = o.paulistrings[i]
@inline Base.getindex(o::OperatorUnsorted, p::PauliString) = get(o.paulistrings, p, zero(scalartype(o)))

@inline Base.iterate(o::OperatorUnsorted, args...) = iterate(pairs(o.paulistrings), args...)

function Base.:+(o1::OperatorUnsorted, o2::OperatorUnsorted)
    if length(o1) > length(o2)
        paulistrings = mergewith(+, o1.paulistrings, o2.paulistrings)
    else
        paulistrings = mergewith(+, o2.paulistrings, o1.paulistrings)
    end
    return OperatorUnsorted(paulistrings)
end
Base.:+(o1::OperatorUnsorted{P}, o2::P) where {P<:PauliString} = o1 + OperatorUnsorted([o2 => one(scalartype(o1))])
Base.:+(o1::P, o2::OperatorUnsorted{P}) where {P<:PauliString} = OperatorUnsorted([o1 => one(scalartype(o2))]) + o2

Base.:-(o::OperatorUnsorted) = OperatorUnsorted(map(-, o.paulistrings))
Base.:-(o1::OperatorUnsorted, o2::OperatorUnsorted) = o1 + (-o2)

function Base.:*(o::OperatorUnsorted, λ::Number)
    return OperatorUnsorted(map(Base.Fix2(*, λ), o.paulistrings))
end
Base.:*(λ::Number, o::OperatorUnsorted) = o * λ
Base.:/(o::OperatorUnsorted, λ::Number) = o * inv(λ)
Base.:\(λ::Number, o::OperatorUnsorted) = o * inv(λ)

function Base.:*(o1::O, o2::O) where {O<:OperatorUnsorted}
    T = Base.promote_op(*, scalartype(o1), scalartype(o2))
    result = similar(o1, T)

    @inbounds for (p1, c1) in o1, (p2, c2) in o2
        p, k = p1 * p2
        k == 0 && continue
        c = k * c1 * c2
        hadtoken, token = gettoken!(result.paulistrings, p)
        hadtoken && (c += gettokenvalue(result.paulistrings, token))
        settokenvalue!(result.paulistrings, token, c)
    end
    return result
end
function Base.:*(o1::OperatorUnsorted{P}, o2::P) where {P}
    return o1 * OperatorUnsorted([o2 => one(scalartype(o1))])
end

function commutator(o1::O, o2::O; maxlength::Int=1000, epsilon::Real=0) where {O<:OperatorUnsorted}
    T = Base.promote_op(*, scalartype(o1), scalartype(o2))
    Teps = real(T(epsilon))
    sizehint = max(length(o1), length(o2)) # conservative guess
    paulistrings = Dictionary{stringtype(O),T}(; sizehint)
    keys1 = collect(keys(o1.paulistrings))
    keys2 = collect(keys(o2.paulistrings))
    vals1 = collect(values(o1.paulistrings))
    vals2 = collect(values(o2.paulistrings))
    @inbounds for i in eachindex(keys1), j in eachindex(keys2)
        p, k = commutator(keys1[i], keys2[j])
        c = k * vals1[i] * vals2[j]

        (k != 0 && abs(c) > Teps && pauli_weight(p) < maxlength) || continue
        setwith!(+, paulistrings, p, c)
    end
    # o1_iterator = zip(o1.paulistrings.indices.values, o1.paulistrings.values)
    # o2_iterator = zip(o2.paulistrings.indices.values, o2.paulistrings.values)
    # @inbounds for (p1, c1) in o1_iterator, (p2, c2) in o2_iterator
    #     p, k = commutator(p1, p2)
    #     c = k * c1 * c2

    #     (k != 0 && abs(c) > Teps && pauli_weight(p) < maxlength) || continue
    #     setwith!(+, paulistrings, p, c)
    # end
    return OperatorUnsorted(paulistrings)
end

function anticommutator(o1::O, o2::O; maxlength::Int=1000, epsilon::Real=0) where {O<:OperatorUnsorted}
    T = Base.promote_op(*, scalartype(o1), scalartype(o2))
    sizehint = max(length(o1), length(o2)) # conservative guess
    paulistrings = Dictionary{stringtype(O),T}(; sizehint)
    @inbounds for (p1, c1) in o1, (p2, c2) in o2
        p, k = anticommutator(p1, p2)
        c = k * c1 * c2

        (k != 0 && abs(c) > epsilon && pauli_weight(p) < maxlength) || continue
        setwith!(+, paulistrings, p, c)
    end
    return OperatorUnsorted(paulistrings)
end

function opnorm(o::OperatorUnsorted)
    return norm(values(o.paulistrings)) / (2.0^(length(stringtype(o)) / 2))
end

function norm_lanczos(o::OperatorUnsorted)
    return norm(values(o.paulistrings))
end

function trim(o::OperatorUnsorted, max_strings::Int; keepnorm::Bool=false, keep::Operator=zero(o))
    length(o) ≤ max_strings && return copy(o)

    vals = o.paulistrings.values
    I = partialsortperm(vals, 1:max_strings; by=abs, rev=true)
    ks = @view keys(o.paulistrings).values[I]
    vs = @view vals[I]
    paulistrings = Dictionary(ks, vs)

    @assert length(keep) == 0 "TBA"

    o′ = OperatorUnsorted(paulistrings)
    return keepnorm ? o′ * (opnorm(o) / opnorm(o′)) : o′
end

# Glue
# ----
const NewOperators = Union{OperatorSorted,OperatorUnsorted}
function com(o1::NewOperators, o2::NewOperators; anti::Bool=false, kwargs...)
    return anti ? anticommutator(o1, o2; kwargs...) : commutator(o1, o2; kwargs...)
end

OperatorSorted(x::OperatorSorted) = copy(x)
OperatorUnsorted(x::OperatorUnsorted) = copy(x)

function OperatorSorted(x::Operator)
    if length(x) == 0
        B = eltype(x.v)
        return OperatorSorted(PauliString{x.N,B}[])
    else
        paulistrings = map(x.v, x.w, x.coef) do v, w, c
            return PauliString{x.N}(v, w) => c
        end
        return OperatorSorted(sort!(paulistrings; by=first))
    end
end
function OperatorUnsorted(x::Operator)
    if length(x) == 0
        B = eltype(x.v)
        return OperatorUnsorted(PauliString{x.N,B}[])
    else
        strings = map(x.v, x.w) do v, w
            return PauliString{x.N}(v, w)
        end
        return OperatorUnsorted(Dictionary(strings, x.coef))
    end
end

# Utility
# -------
function merge_sorted(v1, v2; by=identity)
    isempty(v1) && return copy(v2)
    isempty(v2) && return copy(v1)
    v = similar(v1, length(v1) + length(v2))
    i1 = i2 = 1
    @inbounds for i in 1:length(v)
        if i2 > length(v2) || (i1 ≤ length(v1) && by(v1[i1]) ≤ by(v2[i2]))
            v[i] = v1[i1]
            i1 += 1
        else
            v[i] = v2[i2]
            i2 += 1
        end
    end
    # @assert i1 == length(v1) && i2 == length(v2)
    return v
end
