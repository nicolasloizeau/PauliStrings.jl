const _SIZEHINT = 16

mutable struct SwissIndices{K} <: AbstractIndices{K}
    "swiss table metadata"
    slots::VecType{_u8x16}

    "open-addressed storage of keys"
    keys::VecType{K}

    "number of filled buckets"
    nbfull::Int

    "number of stored keys"
    count::Int

    "number of deletions"
    ndel::Int

    "index smaller than the first filled slot - speed up iteration"
    idxfloor::Int

    "maximal number of probes before we stop searching"
    maxprobe::Int
end

SwissIndices() = SwissIndices{Any}()

function SwissIndices{K}(; sizehint::Int=_SIZEHINT) where {K}
    sz = Base._tablesz(sizehint)
    slots = fill(_expand16(0x00), sz >> 4)
    keys = VecType{K}(undef, sz)
    return SwissIndices{K}(slots, keys, 0, 0, 0, 1, 16)
end

function SwissIndices(iter)
    if Base.IteratorEltype(iter) == Base.EltypeUnknown()
        iter = collect(iter)
    end
    return SwissIndices{eltype(iter)}(iter)
end

function SwissIndices{T}(iter) where {T}
    if Base.IteratorSize(iter) isa Union{Base.HasLength,Base.HasShape}
        h = SwissIndices{T}(; sizehint=(length(iter) * 5) >>> 0x02)
    else
        h = SwissIndices{T}()
    end

    for item in iter
        insert!(h, item)
    end

    return h
end

Base.length(h::SwissIndices) = h.count

# Token Interface
# ---------------
Dictionaries.istokenizable(::SwissIndices) = true
Dictionaries.tokentype(::SwissIndices) = Int

function Dictionaries.iteratetoken(h::SwissIndices{K}) where {K}
    i0 = ((h.idxfloor - 1) & (length(h.keys) - 1)) >> 4 + 1
    off = (h.idxfloor - 1) & 0x0f
    @inbounds slot = _find_free(h.slots[i0>>4+1])
    slot = ((~slot & 0xffff) >> off) << off
    return Dictionaries.iteratetoken(h, (i0, slot))
end

function Dictionaries.iteratetoken(h::SwissIndices{K}, (i, slot)) where {K}
    while iszero(slot)
        i += 1
        i <= length(h.slots) || return nothing
        @inbounds mask = h.slots[i]
        slot = _find_free(mask)
        slot = (~slot & 0xffff)
    end
    return (i - 1) * 16 + trailing_zeros(slot) + 1, (i, _blsr(slot))
end

# function Dictionaries.iteratetoken_reverse(h::SwissIndices{K}) where {K}
#     start = length(h.keys) - rem(length(h.keys), 16) - 1
#     i0 = ((start - 1) & (length(h.keys) - 1)) >> 4 + 1
#     off = (start - 1) & 0x0f
#     @inbounds slot = _find_free(h.slots[i0>>4+1])
#     slot = ((~slot & 0xffff) >> off) << off
#     return Dictionaries.iteratetoken_reverse(h, (i0, slot))
# end
# function Dictionaries.iteratetoken_reverse(h::SwissIndices{K}, (i, slot)) where {K}
#     while iszero(slot)
#         i -= 1
#         i < 1 && return nothing
#         @inbounds mask = h.slots[i]
#         slot = _find_free(mask)
#         slot = (~slot & 0xffff)
#     end
#     return (i - 1) * 16 + trailing_zeros(slot) + 1, (i, _blsr(slot))
# end
# TODO: iteratetoken_reverse

function Dictionaries.gettoken(h::SwissIndices{K}, key::K) where {K}
    slots = h.slots
    keys = h.keys
    sz = length(slots)

    i0, tag = _hashtag(hash(key))
    _prefetchr(pointer(keys, i0 * 16 + 1))

    # check if the key is already present
    i = i0 & (sz - 1) # sz is a power of 2
    while true                # outer loop - blocks of 16x8
        msk = slots[i+1]
        candidates, done = _find_candidates(msk, tag)
        while !iszero(candidates)       # inner loop - bits of 16
            offset = trailing_zeros(candidates)
            idx = i * 16 + offset + 1
            isequal(keys[idx], key) && return (true, idx)
            candidates = _blsr(candidates)
        end
        done && break
        i = (i + 1) & (sz - 1)
    end

    return false, -1
end

@inline Dictionaries.gettokenvalue(h::SwissIndices, token::Int) = h.keys[token]
@inline Dictionaries.settokenvalue!(h::SwissIndices{K}, token, value::K) where {K} = h.keys[token] = value

function Dictionaries.gettoken!(h::SwissIndices{K}, key::K) where {K}
    i0, tag = _hashtag(hash(key))
    found, token = _gettoken!(h, key, i0, tag)
    if !found
        @inbounds _insert!(h, key, token, tag)
        if maybe_rehash_grow!(h)
            _, token = _gettoken!(h, key, io, tag)
        end
    end
    return found, token
end

function _gettoken!(h::SwissIndices{K}, key::K, i0, tag) where {K}
    slots = h.slots
    keys = h.keys
    sz = length(slots)

    _prefetchr(pointer(keys, i0 * 16 + 1))

    # check if the key is already present
    i = i0 & (sz - 1) # sz is a power of 2
    while true                # outer loop - blocks of 16x8
        msk = slots[i+1]
        candidates, done = _find_candidates(msk, tag)
        while !iszero(candidates)       # inner loop - bits of 16
            offset = trailing_zeros(candidates)
            idx = i * 16 + offset + 1
            isequal(keys[idx], key) && return (true, idx)
            candidates = _blsr(candidates)
        end
        done && break
        i = (i + 1) & (sz - 1)
    end

    # if not found, check for the first free slot
    # TODO: if we disallow deletions/detect no deletions, we don't need to check again?
    i = i0 & (sz - 1) # sz is a power of 2
    while true
        msk = slots[i+1]
        candidates = _find_free(msk)
        if !iszero(candidates)
            offset = trailing_zeros(candidates)
            idx = i * 16 + offset + 1
            return (false, idx)
        end
        i = (i + 1) & (sz - 1)
    end
end

function _insert!(h::SwissIndices{K}, key::K, token::Int, tag) where {K}
    so = _slotget(h.slots, token)
    h.nbfull += (iszero(token & 0x0f) & (so == 0x00))
    _slotset!(h.slots, tag, token)
    h.idxfloor = min(h.idxfloor, token)
    h.count += 1
    h.keys[token] = key
    return h
end

# Rehashing
# ---------
function maybe_rehash_grow!(h::SwissIndices)
    sz = length(h.keys)
    hastogrow = h.count > sz * SWISS_DICT_LOAD_FACTOR || (h.nbfull - 1) * 10 > sz * 6
    hastogrow && rehash!(h, sz << 2)
    return hastogrow
end

function maybe_rehash_shrink!(h::SwissIndices)
    sz = length(h.keys)
    if h.count * 4 < sz && sz > 16
        rehash!(h, sz >> 1)
    end
end

function Base.sizehint!(d::SwissIndices, newsz)
    newsz = _tablesz(newsz)
    oldsz = length(d.keys)
    # grow at least 25%
    if newsz < (oldsz * 5) >> 2
        return d
    end
    rehash!(d, newsz)
end

function rehash!(h::SwissIndices{K}, newsz=length(h.keys)) where {K}
    newsz = Base._tablesz(newsz)
    (newsz * SWISS_DICT_LOAD_FACTOR) > length(h) || (newsz <<= 1)

    if isempty(h)
        resize!(h.slots, newsz >> 4)
        fill!(h.slots, _expand16(0x00))
        resize!(h.keys, newsz)
        h.nbfull = 0
        h.idxfloor = 1
        return h
    end

    h′ = SwissIndices{K}(; sizehint=newsz)
    for k in h
        insert!(h′, k)
    end
    h.slots = h′.slots
    h.keys = h′.keys
    h.count = h′.count
    h.nbfull = h′.nbfull
    h.ndel = h′.ndel
    h.idxfloor = h′.idxfloor

    return h
end
