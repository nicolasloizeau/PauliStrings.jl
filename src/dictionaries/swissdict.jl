using Dictionaries
using Base: @propagate_inbounds

# Swiss dict
# ----------
struct SwissDictionary{K,V} <: AbstractDictionary{K,V}
    indices::SwissIndices{K}
    values::Vector{V}

    function SwissDictionary{K,V}(inds::SwissIndices{K}, values::Vector{V}) where {K,V}
        length(values) == length(inds.keys) || throw(DimensionMismatch())
        return new{K,V}(inds, values)
    end

    # slots::VecType{_u8x16}
    # keys::VecType{K}
    # vals::VecType{V}
    # nbfull::Int
    # count::Int
    # age::UInt
    # idxfloor::Int  # an index <= the indices of all used slots

    # function SwissDictionary{K,V}() where {K,V}
    #     new(fill(_expand16(0x00), 1), Vector{K}(undef, 16), Vector{V}(undef, 16), 0, 0, 0, 1)
    # end
    # function SwissDictionary{K,V}(d::SwissDictionary{K,V}) where {K,V}
    #     new(copy(d.slots), copy(d.keys), copy(d.vals), d.nbfull, d.count, d.age,
    #         d.idxfloor)
    # end
    # function SwissDictionary{K,V}(slots, keys, vals, nbfull, count, age, idxfloor) where {K,V}
    #     new(slots, keys, vals, nbfull, count, age, idxfloor)
    # end
end

SwissDictionary(; sizehint=8) = SwissDictionary{Any,Any}(; sizehint)
SwissDictionary{K}(; sizehint=8) where {K} = SwissDictionary{K,Any}(; sizehint)
function SwissDictionary{K,V}(; sizehint::Int=_SIZEHINT) where {K,V}
    indices = SwissIndices{K}(; sizehint)
    values = Vector{V}(undef, length(indices.keys))
    return SwissDictionary{K,V}(indices, values)
end

function SwissDictionary{K,V}(kv) where {K,V}
    h = SwissDictionary{K,V}()
    for (k, v) in kv
        h[k] = v
    end
    return h
end
SwissDictionary{K,V}(p::Pair) where {K,V} = setindex!(SwissDictionary{K,V}(), p.second, p.first)
function SwissDictionary{K,V}(ps::Pair...) where {K,V}
    h = SwissDictionary{K,V}()
    sizehint!(h, length(ps))
    for p in ps
        h[p.first] = p.second
    end
    return h
end
Base.copy(d::SwissDictionary) = SwissDictionary(d)
Base.empty(d::SwissDictionary, ::Type{K}, ::Type{V}) where {K,V} = SwissDictionary{K,V}()

function SwissDictionary(kv)
    try
        dict_with_eltype((K, V) -> SwissDictionary{K,V}, kv, eltype(kv))
    catch e
        if !isiterable(typeof(kv)) || !all(x -> isa(x, Union{Tuple,Pair}), kv)
            throw(ArgumentError("SwissDictionary(kv): kv needs to be an iterator of tuples or pairs"))
        else
            rethrow(e)
        end
    end
end

# Token interface
# ---------------
Base.keys(dict::SwissDictionary) = dict.indices
# Base.values(dict::SwissDictionary) = dict.values


Dictionaries.istokenizable(::SwissDictionary) = true
Dictionaries.tokentype(::SwissDictionary) = Int

Dictionaries.tokenized(d::SwissDictionary) = d.values

function Dictionaries.istokenassigned(dict::SwissDictionary, (_slot, index))
    return isassigned(dict.values, index)
end

function Dictionaries.istokenassigned(dict::SwissDictionary, index::Int)
    return isassigned(dict.values, index)
end

@propagate_inbounds function Dictionaries.gettokenvalue(dict::SwissDictionary, (_slot, index))
    return dict.values[index]
end

@propagate_inbounds function Dictionaries.gettokenvalue(dict::SwissDictionary, index::Int)
    return dict.values[index]
end

Dictionaries.issettable(::SwissDictionary) = true

@propagate_inbounds function Dictionaries.settokenvalue!(dict::SwissDictionary{<:Any,T}, (_slot, index), value::T) where {T}
    dict.values[index] = value
    return dict
end

@propagate_inbounds function Dictionaries.settokenvalue!(dict::SwissDictionary{<:Any,T}, index::Int, value::T) where {T}
    dict.values[index] = value
    return dict
end

# insertion

Dictionaries.isinsertable(::SwissDictionary) = true

function Dictionaries.gettoken!(dict::SwissDictionary{I}, key::I) where {I}
    i0, tag = _hashtag(hash(key))
    found, token = _gettoken!(keys(dict), key, i0, tag)
    if !found
        @inbounds _insert!(keys(dict), key, token, tag)
        if maybe_rehash_grow!(dict)
            _, token = _gettoken!(keys(dict), key, i0, tag)
        end
    end
    return found, token
end

function maybe_rehash_grow!(dict::SwissDictionary)
    h = keys(dict)
    sz = length(h.keys)
    if h.count > sz * SWISS_DICT_LOAD_FACTOR || (h.nbfull - 1) * 10 > sz * 6
        rehash!(dict, sz << 2)
        return true
    end
    return false
end

function rehash!(dict::SwissDictionary, newsz=length(keys(dict).keys))
    h = keys(dict)
    newsz = Base._tablesz(newsz)
    (newsz * SWISS_DICT_LOAD_FACTOR) > length(h) || (newsz <<= 1)

    if isempty(h)
        resize!(h.slots, newsz >> 4)
        fill!(h.slots, _expand16(0x00))
        resize!(h.keys, newsz)
        resize!(dict.values, newsz)
        h.nbfull = 0
        h.idxfloor = 1
        return dict
    end

    dict′ = SwissDictionary{keytype(dict),valtype(dict)}(; sizehint=newsz)
    for (k, v) in pairs(dict)
        insert!(dict′, k, v)
    end
    h′ = keys(dict′)
    h.slots = h′.slots
    h.keys = h′.keys
    h.count = h′.count
    h.nbfull = h′.nbfull
    h.ndel = h′.ndel
    h.idxfloor = h′.idxfloor
    resize!(dict.values, newsz)
    copyto!(dict.values, dict′.values)

    return dict
end


function Dictionaries.deletetoken!(dict::SwissDictionary, (slot, index))
    isbitstype(valtype(dict)) || ccall(:jl_arrayunset, Cvoid, (Any, UInt), dict.values, index - 1)
    deletetoken!(keys(dict), (slot, index), (dict.values,))
    return dict
end

function Base.empty!(dict::SwissDictionary)
    empty!(dict.values)
    empty!(keys(dict))
    return dict
end

# Basic operations

# get the index where a key is stored, or -1 if not present
# ht_keyindex(h::SwissDictionary, key) = ht_keyindex(h::SwissDictionary, key, _hashtag(hash(key))...)
# function ht_keyindex(h::SwissDictionary, key, i0, tag)
#     slots = h.slots
#     keys = h.keys
#     sz = length(slots)
#     i = i0 & (sz - 1)
#     _prefetchr(pointer(h.keys, i * 16 + 1))
#     _prefetchr(pointer(h.vals, i * 16 + 1))
#     #Todo/discuss: _prefetchr(pointer(h.keys, i*16+9))?
#     @inbounds while true
#         msk = slots[i+1]
#         cands, done = _find_candidates(msk, tag)
#         while cands != 0
#             off = trailing_zeros(cands)
#             idx = i * 16 + off + 1
#             isequal(keys[idx], key) && return idx
#             cands = _blsr(cands)
#         end
#         done && break
#         i = (i + 1) & (sz - 1)
#     end
#     return -1
# end

# # get the index where a key is stored, or -pos if not present
# # and the key would be inserted at pos
# # This version is for use by setindex! and get!. It never rehashes.
# ht_keyindex2!(h::SwissDictionary, key) = ht_keyindex2!(h, key, _hashtag(hash(key))...)
# @inline function ht_keyindex2!(h::SwissDictionary, key, i0, tag)
#     slots = h.slots
#     keys = h.keys
#     sz = length(slots)
#     i = i0 & (sz - 1)
#     _prefetchw(pointer(h.keys, i * 16 + 1))
#     _prefetchw(pointer(h.vals, i * 16 + 1))
#     #Todo/discuss: _prefetchr(pointer(h.keys, i*16+9))?
#     @inbounds while true
#         msk = slots[i+1]
#         cands, done = _find_candidates(msk, tag)
#         while cands != 0
#             off = trailing_zeros(cands)
#             idx = i * 16 + off + 1
#             isequal(keys[idx], key) && return idx, tag
#             cands = _blsr(cands)
#         end
#         done && break
#         i = (i + 1) & (sz - 1)
#     end
#     i = i0 & (sz - 1)
#     @inbounds while true
#         msk = slots[i+1]
#         cands = _find_free(msk)
#         if cands != 0
#             off = trailing_zeros(cands)
#             idx = i * 16 + off + 1
#             return -idx, tag
#         end
#         i = (i + 1) & (sz - 1)
#     end
# end

# function _setindex!(h::SwissDictionary, v, key, index, tag)
#     @inbounds h.keys[index] = key
#     @inbounds h.vals[index] = v
#     h.count += 1
#     h.age += 1
#     so = _slotget(h.slots, index)
#     h.nbfull += (iszero(index & 0x0f) & (so == 0x00))
#     _slotset!(h.slots, tag, index)
#     if index < h.idxfloor
#         h.idxfloor = index
#     end
#     maybe_rehash_grow!(h)
# end

# function _delete!(h::SwissDictionary{K,V}, index) where {K,V}
#     # Caller is responsible for maybe shrinking the SwissDictionary after the deletion.
#     isbitstype(K) || isbitsunion(K) || ccall(:jl_arrayunset, Cvoid, (Any, UInt), h.keys, index - 1)
#     isbitstype(V) || isbitsunion(V) || ccall(:jl_arrayunset, Cvoid, (Any, UInt), h.vals, index - 1)
#     isboundary = iszero(index & 0x0f) #boundaries: 16, 32, ...
#     @inbounds _slotset!(h.slots, ifelse(isboundary, 0x01, 0x00), index)
#     h.count -= 1
#     h.age += 1
#     maybe_rehash_shrink!(h)
# end


# # fast iteration over active slots.
# function _iterslots(h::SwissDictionary, start::Int)
#     i0 = ((start - 1) & (length(h.keys) - 1)) >> 4 + 1
#     off = (start - 1) & 0x0f
#     @inbounds sl = _find_free(h.slots[i0>>4+1])
#     sl = ((~sl & 0xffff) >> off) << off
#     return _iterslots(h, (i0, sl))
# end

# function _iterslots(h::SwissDictionary, state)
#     i, sl = state
#     while iszero(sl)
#         i += 1
#         i <= length(h.slots) || return nothing
#         @inbounds msk = h.slots[i]
#         sl = _find_free(msk)
#         sl = (~sl & 0xffff)
#     end
#     return ((i - 1) * 16 + trailing_zeros(sl) + 1, (i, _blsr(sl)))
# end

# # Dictionary resize logic:
# # Guarantee 40% of buckets and 15% of entries free, and at least 25% of entries filled
# # growth when > 85% entries full or > 60% buckets full, shrink when <25% entries full.
# # >60% bucket full should be super rare outside of very bad hash collisions or
# # super long-lived Dictionaries (expected 0.85^16 = 7% buckets full at 85% entries full).
# # worst-case hysteresis: shrink at 25% vs grow at 30% if all hashes collide.
# # expected hysteresis is 25% to 42.5%.
# function maybe_rehash_grow!(h::SwissDictionary)
#     sz = length(h.keys)
#     if h.count > sz * SWISS_DICT_LOAD_FACTOR || (h.nbfull - 1) * 10 > sz * 6
#         rehash!(h, sz << 2)
#     end
# end

# function maybe_rehash_shrink!(h::SwissDictionary)
#     sz = length(h.keys)
#     if h.count * 4 < sz && sz > 16
#         rehash!(h, sz >> 1)
#     end
# end

# function Base.sizehint!(d::SwissDictionary, newsz)
#     newsz = _tablesz(newsz * 2)  # *2 for keys and values in same array
#     oldsz = length(d.keys)
#     # grow at least 25%
#     if newsz < (oldsz * 5) >> 2
#         return d
#     end
#     rehash!(d, newsz)
# end

# function rehash!(h::SwissDictionary{K,V}, newsz=length(h.keys)) where {K,V}
#     olds = h.slots
#     oldk = h.keys
#     oldv = h.vals
#     sz = length(oldk)
#     newsz = _tablesz(newsz)
#     (newsz * SWISS_DICT_LOAD_FACTOR) > h.count || (newsz <<= 1)
#     h.age += 1
#     h.idxfloor = 1
#     if h.count == 0
#         resize!(h.slots, newsz >> 4)
#         fill!(h.slots, _expand16(0x00))
#         resize!(h.keys, newsz)
#         resize!(h.vals, newsz)
#         h.nbfull = 0
#         return h
#     end
#     nssz = newsz >> 4
#     slots = fill(_expand16(0x00), nssz)
#     keys = Vector{K}(undef, newsz)
#     vals = Vector{V}(undef, newsz)
#     age0 = h.age
#     nbfull = 0
#     is = _iterslots(h, 1)
#     count = 0
#     @inbounds while is !== nothing
#         i, s = is
#         k = oldk[i]
#         v = oldv[i]
#         i0, t = _hashtag(hash(k))
#         i = i0 & (nssz - 1)
#         idx = 0
#         while true
#             msk = slots[i+1]
#             cands = _find_free(msk)
#             if cands != 0
#                 off = trailing_zeros(cands)
#                 idx = i * 16 + off + 1
#                 break
#             end
#             i = (i + 1) & (nssz - 1)
#         end
#         _slotset!(slots, t, idx)
#         keys[idx] = k
#         vals[idx] = v
#         nbfull += iszero(idx & 0x0f)
#         count += 1
#         if h.age != age0
#             return rehash!(h, newsz)
#         end
#         is = _iterslots(h, s)
#     end
#     h.slots = slots
#     h.keys = keys
#     h.vals = vals
#     h.nbfull = nbfull
#     @assert h.age == age0
#     @assert h.count == count
#     return h
# end

# Base.isempty(t::SwissDictionary) = (t.count == 0)
# Base.length(t::SwissDictionary) = t.count

# """
#     empty!(collection) -> collection

# Remove all elements from a `collection`.

# # Examples
# ```jldoctest
# julia> A = SwissDictionary("a" => 1, "b" => 2)
# SwissDictionary{String, Int64} with 2 entries:
#   "a" => 1
#   "b" => 2

# julia> empty!(A);


# julia> A
# SwissDictionary{String, Int64}()
# ```
# """
# function Base.empty!(h::SwissDictionary{K,V}) where {K,V}
#     fill!(h.slots, _expand16(0x00))
#     sz = length(h.keys)
#     empty!(h.keys)
#     empty!(h.vals)
#     resize!(h.keys, sz)
#     resize!(h.vals, sz)
#     h.nbfull = 0
#     h.count = 0
#     h.age += 1
#     h.idxfloor = 1
#     return h
# end

# function Base.setindex!(h::SwissDictionary{K,V}, v0, key0) where {K,V}
#     key = convert(K, key0)
#     _setindex!(h, v0, key)
# end

# function _setindex!(h::SwissDictionary{K,V}, v0, key::K) where {K,V}
#     v = convert(V, v0)
#     index, tag = ht_keyindex2!(h, key)

#     if index > 0
#         h.age += 1
#         @inbounds h.keys[index] = key
#         @inbounds h.vals[index] = v
#     else
#         _setindex!(h, v, key, -index, tag)
#     end

#     return h
# end

# """
#     get!(collection, key, default)

# Return the value stored for the given key, or if no mapping for the key is present, store
# `key => default`, and return `default`.

# # Examples
# ```jldoctest
# julia> d = SwissDictionary("a"=>1, "b"=>2, "c"=>3);

# julia> get!(d, "a", 5)
# 1

# julia> get!(d, "d", 4)
# 4

# julia> d
# SwissDictionary{String, Int64} with 4 entries:
#   "a" => 1
#   "b" => 2
#   "c" => 3
#   "d" => 4
# ```
# """
# Base.get!(h::SwissDictionary{K,V}, key0, default) where {K,V} = get!(() -> default, h, key0)

# """
#     get!(f::Function, collection, key)

# Return the value stored for the given key, or if no mapping for the key is present, store
# `key => f()`, and return `f()`.

# This is intended to be called using `do` block syntax:
# ```julia
# get!(dict, key) do
#     # default value calculated here
#     time()
# end
# ```
# """
# function Base.get!(default::Base.Callable, h::SwissDictionary{K,V}, key0) where {K,V}
#     key = convert(K, key0)
#     return _get!(default, h, key)
# end

# function _get!(default::Base.Callable, h::SwissDictionary{K,V}, key::K) where {K,V}
#     index, tag = ht_keyindex2!(h, key)

#     index > 0 && return @inbounds h.vals[index]

#     age0 = h.age
#     v = convert(V, default())
#     if h.age != age0
#         index, tag = ht_keyindex2!(h, key)
#     end
#     if index > 0
#         h.age += 1
#         @inbounds h.keys[index] = key
#         @inbounds h.vals[index] = v
#     else
#         _setindex!(h, v, key, -index, tag)
#     end
#     return v
# end

# function Base.getindex(h::SwissDictionary{K,V}, key) where {K,V}
#     index = ht_keyindex(h, key)
#     @inbounds return (index < 0) ? throw(KeyError(key)) : h.vals[index]::V
# end

# """
#     get(collection, key, default)

# Return the value stored for the given key, or the given default value if no mapping for the
# key is present.

# # Examples
# ```jldoctest
# julia> d = SwissDictionary("a"=>1, "b"=>2);


# julia> get(d, "a", 3)
# 1

# julia> get(d, "c", 3)
# 3
# ```
# """
# function Base.get(h::SwissDictionary{K,V}, key, default) where {K,V}
#     index = ht_keyindex(h, key)
#     @inbounds return (index < 0) ? default : h.vals[index]::V
# end

# """
#     get(f::Function, collection, key)

# Return the value stored for the given key, or if no mapping for the key is present, return
# `f()`.  Use [`get!`](@ref) to also store the default value in the dictionary.

# This is intended to be called using `do` block syntax

# ```julia
# get(dict, key) do
#     # default value calculated here
#     time()
# end
# ```
# """
# function Base.get(default::Base.Callable, h::SwissDictionary{K,V}, key) where {K,V}
#     index = ht_keyindex(h, key)
#     @inbounds return (index < 0) ? default() : h.vals[index]::V
# end

# """
#     haskey(collection, key) -> Bool

# Determine whether a collection has a mapping for a given `key`.

# # Examples
# ```jldoctest
# julia> D = SwissDictionary('a'=>2, 'b'=>3)
# SwissDictionary{Char, Int64} with 2 entries:
#   'a' => 2
#   'b' => 3

# julia> haskey(D, 'a')
# true

# julia> haskey(D, 'c')
# false
# ```
# """
# Base.haskey(h::SwissDictionary, key) = (ht_keyindex(h, key) > 0)
# Base.in(key, v::Base.KeySet{<:Any,<:SwissDictionary}) = (ht_keyindex(v.dict, key) > 0)

# """
#     getkey(collection, key, default)

# Return the key matching argument `key` if one exists in `collection`, otherwise return `default`.

# # Examples
# ```jldoctest
# julia> D = SwissDictionary('a'=>2, 'b'=>3)
# SwissDictionary{Char, Int64} with 2 entries:
#   'a' => 2
#   'b' => 3

# julia> getkey(D, 'a', 1)
# 'a': ASCII/Unicode U+0061 (category Ll: Letter, lowercase)

# julia> getkey(D, 'd', 'a')
# 'a': ASCII/Unicode U+0061 (category Ll: Letter, lowercase)
# ```
# """
# function Base.getkey(h::SwissDictionary{K,V}, key, default) where {K,V}
#     index = ht_keyindex(h, key)
#     @inbounds return (index < 0) ? default : h.keys[index]::K
# end

# function _pop!(h::SwissDictionary, index)
#     @inbounds val = h.vals[index]
#     _delete!(h, index)
#     maybe_rehash_shrink!(h)
#     return val
# end

# """
#     pop!(collection, key[, default])

# Delete and return the mapping for `key` if it exists in `collection`, otherwise return
# `default`, or throw an error if `default` is not specified.

# # Examples
# ```jldoctest
# julia> d = SwissDictionary("a"=>1, "b"=>2, "c"=>3);

# julia> pop!(d, "a")
# 1

# julia> pop!(d, "d")
# ERROR: KeyError: key "d" not found
# [...]

# julia> pop!(d, "e", 4)
# 4
# ```
# """
# function Base.pop!(h::SwissDictionary, key)
#     index = ht_keyindex(h, key)
#     return index > 0 ? _pop!(h, index) : throw(KeyError(key))
# end

# function Base.pop!(h::SwissDictionary, key, default)
#     index = ht_keyindex(h, key)
#     return index > 0 ? _pop!(h, index) : default
# end

# function Base.pop!(h::SwissDictionary)
#     isempty(h) && throw(ArgumentError("SwissDictionary must be non-empty"))
#     is = _iterslots(h, h.idxfloor)
#     @assert is !== nothing
#     idx, s = is
#     @inbounds key = h.keys[idx]
#     @inbounds val = h.vals[idx]
#     _delete!(h, idx)
#     h.idxfloor = idx
#     return key => val
# end

# """
#     delete!(collection, key)

# Delete the mapping for the given key in a collection, and return the collection.

# # Examples
# ```jldoctest
# julia> d = SwissDictionary("a"=>1, "b"=>2)
# SwissDictionary{String, Int64} with 2 entries:
#   "a" => 1
#   "b" => 2

# julia> delete!(d, "b")
# SwissDictionary{String, Int64} with 1 entry:
#   "a" => 1
# ```
# """
# function Base.delete!(h::SwissDictionary, key)
#     index = ht_keyindex(h, key)
#     if index > 0
#         _delete!(h, index)
#     end
#     maybe_rehash_shrink!(h)
#     return h
# end

# Base.@propagate_inbounds function Base.iterate(h::SwissDictionary, state=h.idxfloor)
#     is = _iterslots(h, state)
#     is === nothing && return nothing
#     i, s = is
#     @inbounds p = h.keys[i] => h.vals[i]
#     return (p, s)
# end

# Base.@propagate_inbounds function Base.iterate(v::Union{Base.KeySet{<:Any,<:SwissDictionary},Base.ValueIterator{<:SwissDictionary}}, state=v.dict.idxfloor)
#     is = _iterslots(v.dict, state)
#     is === nothing && return nothing
#     i, s = is
#     return (v isa KeySet ? v.dict.keys[i] : v.dict.vals[i], s)
# end
