const SWISS_DICT_LOAD_FACTOR = 0.97
const _u8x16 = NTuple{16,VecElement{UInt8}}

# attempt to be backwards compatible
@static if Core.isdefined(Main, :Memory)
    const VecType{T} = Memory{T}
else
    const VecType{T} = Vector{T}
end

# SIMD utilities
# This enables the swiss dict core concepts: we use a cheap probe to compare 16 values at once
# before we start comparing the keys.
@inline _expand16(u::UInt8) = ntuple(i -> VecElement(u), Val(16))
_blsr(i::UInt32) = i & (i - Int32(1))

@inline _vcmp_eq(u::_u8x16, v::_u8x16) = Core.Intrinsics.llvmcall(("""
%cmp = icmp eq <16 x i8> %0, %1
%cmp16 = bitcast <16 x i1> %cmp to i16
%res = zext i16 %cmp16 to i32
ret i32 %res
"""), UInt32, Tuple{_u8x16,_u8x16}, u, v)

@inline _vcmp_le(u::_u8x16, v::_u8x16) = Core.Intrinsics.llvmcall(("""
%cmp = icmp ule <16 x i8> %0, %1
%cmp16 = bitcast <16 x i1> %cmp to i16
%res = zext i16 %cmp16 to i32
ret i32 %res
"""), UInt32, Tuple{_u8x16,_u8x16}, u, v)

# prefetching memory for reading purposes
@inline function _prefetchr(p::Ptr)
    ccall("llvm.prefetch", llvmcall, Cvoid, (Ref{Int8}, Int32, Int32, Int32), Ptr{Int8}(p), 0, 3, 1)
end

# prefetching memory for writing purposes
@inline function _prefetchw(p::Ptr)
    ccall("llvm.prefetch", llvmcall, Cvoid, (Ref{Int8}, Int32, Int32, Int32), Ptr{Int8}(p), 1, 3, 1)
end

@inline function _hashtag(u::Unsigned)
    #extracts tag between 0x02 and 0xff from lower bits, rotates tag bits to front
    u = u % UInt
    tag = u % UInt8
    if UInt === UInt64
        hi = ((u >> 8) | (u << 56)) % Int
    else
        hi = ((u >> 8) | (u << 24)) % Int
    end
    tag = tag > 1 ? tag : tag + 0x02
    return (hi, tag)
end

Base.@propagate_inbounds function _slotget(slots::VecType{_u8x16}, i::Int)
    @boundscheck 0 < i <= length(slots) * 16 || throw(BoundsError(slots, 1 + (i - 1) >> 4))
    GC.@preserve slots begin
        return unsafe_load(convert(Ptr{UInt8}, pointer(slots)), i)
    end
end

Base.@propagate_inbounds function _slotset!(slots::VecType{_u8x16}, v::UInt8, i::Int)
    @boundscheck 0 < i <= length(slots) * 16 || throw(BoundsError(slots, 1 + (i - 1) >> 4))
    GC.@preserve slots begin
        return unsafe_store!(convert(Ptr{UInt8}, pointer(slots)), v, i)
    end
end

@inline function _find_candidates(v::_u8x16, tag::UInt8)
    match = _vcmp_eq(v, _expand16(tag))
    return (match, v[16].value === 0x00)
end

@inline _find_free(v::_u8x16) = _vcmp_le(v, _expand16(UInt8(1)))
