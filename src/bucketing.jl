# Bucketed, multithreaded binary kernel
# ======================================
#
# For plain `PauliString`s, every binary kernel function (`prod`, `commutator`,
# `anticommutator`) produces the same output string `p = p₁ ⊻ p₂` (only the
# scalar phase `k` differs). If strings are bucketed by a hash `h` that is
# *GF(2)-linear*, i.e. `h(p₁ ⊻ p₂) == h(p₁) ⊻ h(p₂)`, then the output bucket of a
# product is `c = h(p₁) ⊻ h(p₂)` — fully determined by the input buckets.
#
# Consequences:
#   * every output string lands in exactly one bucket, so output buckets are
#     *disjoint* and each can own a private accumulator with no contention;
#   * for a fixed output bucket `c`, the only contributing input-bucket pairs
#     are `(a, a ⊻ c)`, so the `|A|·|B|` work is partitioned into independent
#     chunks that parallelize over output buckets.
#
# The translation-invariant kernel (`PauliStringTS`) does *not* satisfy this: a
# single product is summed over `all_shifts`, emitting many output strings at
# different XOR values, so the single-output-bucket invariant breaks. The
# `is_translation_invariant` trait below gates bucketing off for those.
#
# Note: the Fibonacci `hash` in `operations.jl` is *not* XOR-linear (it uses
# `muladd`), so the bucketing hashes below are implemented as fresh linear maps.

"""
    is_translation_invariant(o::AbstractOperator) -> Bool
    is_translation_invariant(::Type{<:AbstractPauliString}) -> Bool

Trait marking whether an operator's Pauli-string type is translation invariant.
Translation-invariant strings (`PauliStringTS`) sum each product over all lattice
shifts, so the bucketed binary kernel — which assumes a product's output string
is the single `p₁ ⊻ p₂` — cannot be used for them; [`default_strategy`](@ref) and
[`binary_kernel!`](@ref) dispatch on this trait. Evaluated in the type domain, so
it folds to a compile-time constant.
"""
is_translation_invariant(::Type{<:PauliString}) = false
is_translation_invariant(::Type{<:PauliStringTS}) = true
is_translation_invariant(o::AbstractOperator) = is_translation_invariant(paulistringtype(o))

"""
    AbstractBucketStrategy

Supertype for strategies selecting how [`binary_kernel!`](@ref) buckets Pauli
strings to enable contention-free multithreading. A bucketing strategy `s`
defines two methods:

  * `nbits(s)::Int` — there are `2^nbits(s)` buckets;
  * `bucketindex(s, p::PauliString)::Int` — the 0-based bucket of `p`, in
    `0:2^nbits(s)-1`.

The bucket map *must* be GF(2)-linear:

    bucketindex(s, p₁ ⊻ p₂) == bucketindex(s, p₁) ⊻ bucketindex(s, p₂)

This is what makes output buckets disjoint. Use [`is_xor_linear`](@ref) to check
a strategy against this contract.
"""
abstract type AbstractBucketStrategy end

"""
    Serial()

Bucketing strategy selecting the single-threaded reference kernel. This is the
fallback used by [`default_strategy`](@ref) for translation-invariant operators.
"""
struct Serial <: AbstractBucketStrategy end

"""
    nbits(s::AbstractBucketStrategy) -> Int

Number of bucket-index bits; there are `2^nbits(s)` buckets.
"""
function nbits end

"""
    bucketindex(s::AbstractBucketStrategy, p::PauliString) -> Int

0-based bucket index of `p`, in `0:2^nbits(s)-1`.
"""
function bucketindex end

# Hash families
# -------------

"""
    XorVW{B}()

Bucket by the low `B` bits of `v ⊻ w`: `bucketindex = (p.v ⊻ p.w) & (2^B - 1)`.
Cheapest family; mixes both bitstrings. `2^B` buckets.
"""
struct XorVW{B} <: AbstractBucketStrategy end

nbits(::XorVW{B}) where {B} = B
@inline function bucketindex(::XorVW{B}, p::PauliString) where {B}
    mask = (one(p.v) << B) - one(p.v)
    return Int((p.v ⊻ p.w) & mask)
end

# XOR-fold an unsigned word down to `B` bits by XORing its successive `B`-bit
# chunks. Linear because each output bit is a parity of input bits.
@inline function foldbits(x::T, ::Val{B}) where {T <: Unsigned, B}
    mask = (one(T) << B) - one(T)
    acc = zero(T)
    nchunks = cld(8 * sizeof(T), B)
    for k in 0:(nchunks - 1)
        acc ⊻= (x >> (k * B)) & mask
    end
    return acc
end

"""
    Folded{B}()

Bucket by XOR-folding the *whole* `v` and `w` words down to `B` bits (XOR of
successive `B`-bit chunks, then `fold(v) ⊻ fold(w)`). Unlike [`XorVW`](@ref),
high qubits also contribute, which improves load balance for spatially-local
operators where the low qubits dominate. `2^B` buckets.
"""
struct Folded{B} <: AbstractBucketStrategy end

nbits(::Folded{B}) where {B} = B
@inline function bucketindex(::Folded{B}, p::PauliString) where {B}
    return Int(foldbits(p.v, Val(B)) ⊻ foldbits(p.w, Val(B)))
end

"""
    LinearMatrix{B,T}(masks::NTuple{B,Tuple{T,T}})

Most general GF(2)-linear bucketing: output bit `i` is the parity of selected
bits of `v` and `w`,

    bitᵢ = (count_ones(p.v & avᵢ) ⊻ count_ones(p.w & awᵢ)) & 1

with `masks[i] = (avᵢ, awᵢ)`. `2^B` buckets, at a cost of `B` popcounts per
string. Build with [`lowbits_matrix`](@ref), [`random_matrix`](@ref) or
[`mixing_matrix`](@ref).
"""
struct LinearMatrix{B, T <: Unsigned} <: AbstractBucketStrategy
    masks::NTuple{B, Tuple{T, T}}
end

nbits(::LinearMatrix{B}) where {B} = B
@inline function bucketindex(s::LinearMatrix{B, T}, p::PauliString{N, T}) where {B, N, T}
    idx = 0
    @inbounds for i in 1:B
        avᵢ, awᵢ = s.masks[i]
        bit = (count_ones(p.v & avᵢ) ⊻ count_ones(p.w & awᵢ)) & 1
        idx |= bit << (i - 1)
    end
    return idx
end

"""
    lowbits_matrix(N, B) -> LinearMatrix

A [`LinearMatrix`](@ref) whose output bit `i` selects bit `i-1` of `v` (and no
bits of `w`) — i.e. the low `B` bits of `v`. Useful as a simple, reproducible
baseline matrix.
"""
function lowbits_matrix(N::Integer, B::Integer)
    T = uinttype(N)
    masks = ntuple(i -> (one(T) << (i - 1), zero(T)), B)
    return LinearMatrix{B, T}(masks)
end

"""
    random_matrix(N, B; rng) -> LinearMatrix

A [`LinearMatrix`](@ref) with random `v`/`w` masks over the `N` qubits. The map
mixes all bits, which tends to balance bucket occupancy at the cost of `B`
popcounts per string. Pass a seeded `rng` for reproducibility; otherwise the
deterministic [`mixing_matrix`](@ref) is preferred for the default path.
"""
function random_matrix(N::Integer, B::Integer; rng::AbstractRNG = Random.default_rng())
    T = uinttype(N)
    qmask = N >= 8 * sizeof(T) ? typemax(T) : (one(T) << N) - one(T)
    masks = ntuple(_ -> (rand(rng, T) & qmask, rand(rng, T) & qmask), B)
    return LinearMatrix{B, T}(masks)
end

"""
    mixing_matrix(N, B) -> LinearMatrix

A deterministic, full-rank, well-mixing GF(2) bucketing matrix over `N` qubits
with `B ≤ N` output bits. Output bit `i` has a distinct pivot qubit `i` in `v`
(so the rows are independent by construction and all `2^B` buckets are
reachable), plus a broad spread of higher `v` qubits and `w` qubits taken from a
per-row rotation of a mixing constant. This is what [`default_strategy`](@ref)
selects for multithreaded products — like [`random_matrix`](@ref) it balances
bucket occupancy well, but is reproducible.
"""
function mixing_matrix(N::Integer, B::Integer)
    B <= N || throw(ArgumentError("mixing_matrix needs B ≤ N (got B=$B, N=$N)"))
    T = uinttype(N)
    qmask = N >= 8 * sizeof(T) ? typemax(T) : (one(T) << N) - one(T)
    magic = fib_magic_64 % T  # low bits of the golden-ratio constant, in type T
    masks = ntuple(B) do i₁
        i = i₁ - 1
        above = qmask & ~((one(T) << i) - one(T)) & ~(one(T) << i)  # qubits > i
        av = (one(T) << i) | (bitrotate(magic, 7i) & above)         # pivot i + high-v spread
        aw = bitrotate(magic, 7i + 3) & qmask                       # w spread
        (av, aw)
    end
    return LinearMatrix{B, T}(masks)
end

# A single bucket: parallel slices of strings and coefficients, exposing the same
# `keys`/`values`/`pairs` interface as `Operator` so it can be iterated uniformly.
struct Bucket{P, T, VS <: AbstractVector{P}, VC <: AbstractVector{T}}
    strings::VS
    coeffs::VC
end

Base.keys(b::Bucket) = b.strings
Base.values(b::Bucket) = b.coeffs
Base.pairs(b::Bucket) = zip(b.strings, b.coeffs)
Base.length(b::Bucket) = length(b.strings)
Base.isempty(b::Bucket) = isempty(b.strings)

# Counting-sort of an operator's strings into contiguous per-bucket ranges.
# Indexable (1-based) as a vector of `Bucket`s: `buckets[a+1]` is the bucket for
# 0-based hash `a`. `keys`/`values` return the full bucket-ordered arrays.
struct Buckets{P, T} <: AbstractVector{Bucket{P, T, SubArray{P, 1, Vector{P}, Tuple{UnitRange{Int}}, true}, SubArray{T, 1, Vector{T}, Tuple{UnitRange{Int}}, true}}}
    offsets::Vector{Int}   # length nb+1
    strings::Vector{P}
    coeffs::Vector{T}
end

Base.size(b::Buckets) = (length(b.offsets) - 1,)
Base.keys(b::Buckets) = b.strings
Base.values(b::Buckets) = b.coeffs

@inline function Base.getindex(b::Buckets, i::Int)
    @boundscheck checkbounds(b, i)
    @inbounds r = (b.offsets[i]):(b.offsets[i + 1] - 1)
    return Bucket(view(b.strings, r), view(b.coeffs, r))
end

function bucketize(s::AbstractBucketStrategy, o::AbstractOperator)
    ks, vs = keys(o), values(o)
    # take element types from the actual data: a bare `PauliString` operand reports
    # `scalartype` Float64 but yields complex coefficients via `values`.
    P = eltype(ks)
    T = eltype(vs)
    n = length(o)
    nb = 1 << nbits(s)

    counts = zeros(Int, nb)
    @inbounds for i in 1:n
        counts[bucketindex(s, ks[i]) + 1] += 1
    end

    offsets = Vector{Int}(undef, nb + 1)
    offsets[1] = 1
    @inbounds for a in 1:nb
        offsets[a + 1] = offsets[a] + counts[a]
    end

    strings = Vector{P}(undef, n)
    coeffs = Vector{T}(undef, n)
    pos = copy(offsets)  # next write index per (1-based) bucket
    @inbounds for i in 1:n
        a = bucketindex(s, ks[i]) + 1
        p = pos[a]
        strings[p] = ks[i]
        coeffs[p] = vs[i]
        pos[a] = p + 1
    end

    return Buckets{P, T}(offsets, strings, coeffs)
end

"""
    is_xor_linear(s, P; nsamples=1000, rng) -> Bool

Check the GF(2)-linearity contract of strategy `s` for Pauli string type `P` on
`nsamples` random pairs: `bucketindex(s, p₁ ⊻ p₂) == bucketindex(s, p₁) ⊻ bucketindex(s, p₂)`.
"""
function is_xor_linear(
        s::AbstractBucketStrategy, ::Type{PauliString{N, T}};
        nsamples::Int = 1000, rng::AbstractRNG = Random.default_rng()
    ) where {N, T}
    P = PauliString{N, T}
    for _ in 1:nsamples
        p₁ = P(rand(rng, T), rand(rng, T))
        p₂ = P(rand(rng, T), rand(rng, T))
        bucketindex(s, p₁ ⊻ p₂) == (bucketindex(s, p₁) ⊻ bucketindex(s, p₂)) || return false
    end
    return true
end

# Default strategy selection
# --------------------------
# Plain `PauliString` operators are always bucketed (the bucketed kernel runs
# sequentially on a single thread); only translation-invariant operators fall
# back to the serial reference kernel.

"""
    default_strategy(A, B) -> AbstractBucketStrategy

Pick the bucketing strategy backing `*`, `commutator`, and `anticommutator` when
none is given. Returns [`Serial`](@ref) for translation-invariant operators (the
only case the bucketed kernel cannot handle); otherwise a deterministic
[`mixing_matrix`](@ref) with `2^b` buckets sitting a few× above `nthreads()`,
which lets the `:greedy` scheduler balance per-bucket work well across models
(see `benchmark/bucketing/loadbalance.jl`).
"""
function default_strategy(A::AbstractOperator, B::AbstractOperator)
    P = paulistringtype(A)
    is_translation_invariant(P) && return Serial()
    N = qubitlength(A)
    b = min(clamp(ceil(Int, log2(max(Threads.nthreads(), 1))) + 3, 4, 10), N)
    return mixing_matrix(N, b)
end

# Bucketed, multithreaded kernel
# ------------------------------

# Accumulate output bucket `c` = Σ_a A[a]·B[a⊻c] into `d`. Lives in its own
# function so the OhMyThreads task body stays small and type-stable. `nonemptyA`
# is the precomputed list of 0-based nonempty A-bucket keys.
function _accumulate_bucket!(d, f::F, c::Int, α, nonemptyA, bucketsA::Buckets, bucketsB::Buckets, maxlength::Int) where {F}
    @inbounds for a in nonemptyA
        bucketB = bucketsB[(a ⊻ c) + 1]  # output bucket c ⟹ b = a ⊻ c
        isempty(bucketB) && continue
        bucketA = bucketsA[a + 1]
        ksA, vsA = keys(bucketA), values(bucketA)
        ksB, vsB = keys(bucketB), values(bucketB)
        for i in 1:length(bucketA)
            p₁ = ksA[i]
            αc₁ = α * vsA[i]
            for j in 1:length(bucketB)
                p, k = f(p₁, ksB[j])
                (iszero(k) || pauli_weight(p) >= maxlength) && continue
                setwith!(+, d, p, αc₁ * vsB[j] * k)
            end
        end
    end
    return d
end

# Seed `d` with the β·C terms of one output bucket.
@inline function _seed_bucket!(d, bucketC::Bucket, β)
    ks, vs = keys(bucketC), values(bucketC)
    @inbounds for t in 1:length(bucketC)
        setwith!(+, d, ks[t], vs[t] * β)
    end
    return d
end

function binary_kernel!(
        f::F, C::AbstractOperator, A::AbstractOperator, B::AbstractOperator,
        α::Number, β::Number, s::AbstractBucketStrategy;
        maxlength::Int = 1000, epsilon::Real = eps(real(scalartype(C)))
    ) where {F}
    checklength(C, A, B)

    T = scalartype(C)
    P = paulistringtype(C)
    is_translation_invariant(P) && throw(
        ArgumentError(
            "bucketed binary_kernel! supports only plain PauliString; use Serial() for translation-invariant operators"
        )
    )

    nb = 1 << nbits(s)
    bucketsA = bucketize(s, A)
    bucketsB = bucketize(s, B)
    seed = !iszero(β)
    bucketsC = seed ? bucketize(s, C) : bucketsA  # bucketsA is unused when !seed

    # 0-based keys of the nonempty A buckets, hoisted out of the per-bucket loop
    # so each task skips empty A buckets without rescanning.
    nonemptyA = [a for a in 0:(nb - 1) if !isempty(bucketsA[a + 1])]

    # per-bucket output-size estimate: the serial heuristic max(|A|,|B|) for the
    # whole output, spread across the `nb` buckets.
    base_hint = cld(max(length(A), length(B)), nb)

    # one (keys, values) pair per output bucket; keys are disjoint across buckets
    outkeys = Vector{Vector{P}}(undef, nb)
    outvals = Vector{Vector{T}}(undef, nb)

    # Greedy scheduler: OhMyThreads spawns ~nthreads tasks that pull output
    # buckets `c` from a shared queue. The reducer dict `d` is task-local
    # (`@local`) — only ~nthreads dicts total, each reused across the buckets its
    # task handles. `empty!(d)` resets it between buckets (capacity retained, no
    # rehash). Each iteration writes only its own output slot, so it is safe
    # despite the non-deterministic task order.
    @tasks for c in 0:(nb - 1)
        @set scheduler = :greedy
        @local d = UnorderedDictionary{P, T}(; sizehint = base_hint)
        empty!(d)
        seed && _seed_bucket!(d, bucketsC[c + 1], β)
        _accumulate_bucket!(d, f, c, α, nonemptyA, bucketsA, bucketsB, maxlength)
        outkeys[c + 1] = collect(keys(d))   # copy out before the next iteration empties d
        outvals[c + 1] = collect(values(d))
    end

    # Assemble output by concatenating the disjoint accumulators in a single
    # preallocated pass, fusing the epsilon cutoff.
    resize!(C, sum(length, outkeys))
    ksC, vsC = keys(C), values(C)
    i = 1
    if epsilon > 0
        ϵ² = epsilon^2
        @inbounds for b in 1:nb
            ks, vs = outkeys[b], outvals[b]
            for j in eachindex(ks, vs)
                ksC[i] = ks[j]
                vsC[i] = vs[j]
                i += abs2(vs[j]) > ϵ²
            end
        end
        resize!(C, i - 1)
    else
        @inbounds for b in 1:nb
            ks, vs = outkeys[b], outvals[b]
            for j in eachindex(ks, vs)
                ksC[i] = ks[j]
                vsC[i] = vs[j]
                i += 1
            end
        end
    end

    return C
end
