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

# GF(2)-linear bucketing matrix
# -----------------------------

"""
    LinearMatrix{B,T}(masks::NTuple{B,Tuple{T,T}})

GF(2)-linear bucketing: output bit `i` is the parity of selected bits of `v` and
`w`,

    bitᵢ = (count_ones(p.v & avᵢ) ⊻ count_ones(p.w & awᵢ)) & 1

with `masks[i] = (avᵢ, awᵢ)`. `2^B` buckets, at a cost of `B` popcounts per
string. Build with [`mixing_matrix`](@ref).
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
    mixing_matrix(N, B) -> LinearMatrix

A deterministic, full-rank, well-mixing GF(2) bucketing matrix over `N` qubits
with `B ≤ N` output bits. Output bit `i` has a distinct pivot qubit `i` in `v`
(so the rows are independent by construction and all `2^B` buckets are
reachable), plus a broad spread of higher `v` qubits and `w` qubits taken from a
per-row rotation of a mixing constant. This is what [`default_strategy`](@ref)
selects for multithreaded products: it balances bucket occupancy well across all
tested models and is reproducible.
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
    L2_TARGET_BYTES

Target byte size of a single per-bucket accumulator dictionary, used by
[`default_strategy`](@ref) to choose the bucket count so that each thread's local
dictionary stays roughly L2-cache-resident (insertions into a cache-resident dict
are what let the multithreaded kernel scale). Conservative per-core default; tune
to your hardware.
"""
const L2_TARGET_BYTES = 512 * 1024

"""
    default_strategy(A, B) -> AbstractBucketStrategy

Pick the bucketing strategy backing `*`, `commutator`, and `anticommutator` when
none is given. Returns [`Serial`](@ref) for translation-invariant operators (the
only case the bucketed kernel cannot handle); otherwise a deterministic
[`mixing_matrix`](@ref) with `2^b` buckets.

`b` is chosen from the *operator sizes* so that each per-bucket accumulator
dictionary stays within [`L2_TARGET_BYTES`](@ref): with the per-bucket occupancy
estimated as `max(|A|,|B|) / 2^b` entries (the same heuristic as the kernel's dict
`sizehint`), `b` grows with the operands so a single bucket's dict stays
L2-resident. A floor of a few× `nthreads()` buckets keeps the scheduler supplied
with enough chunks to balance load.
"""
function default_strategy(A::AbstractOperator, B::AbstractOperator)
    P = paulistringtype(A)
    is_translation_invariant(P) && return Serial()
    N = qubitlength(A)
    nt = max(Threads.nthreads(), 1)
    T = complex(Base.promote_op(*, scalartype(A), scalartype(B)))
    # bytes per accumulator entry: PauliString key (two words) + complex coeff +
    # an allowance for the dictionary's hash-index overhead.
    bytes_per_entry = sizeof(P) + sizeof(T) + 16
    fit = max(L2_TARGET_BYTES ÷ bytes_per_entry, 1)        # entries keeping a bucket in L2
    entries = max(length(A), length(B))
    # enough buckets to (a) keep each per-bucket dict ≈ L2-resident and (b) give the
    # scheduler several× more chunks than threads to balance load. The L2 term only
    # raises `b` above the parallelism floor once the operands exceed ~fit·8·nthreads.
    nbuckets = max(cld(entries, fit), 8 * nt)
    b = clamp(ceil(Int, log2(nbuckets)), min(4, N), min(N, 16))
    return mixing_matrix(N, b)
end

# Bucketed, multithreaded kernel
# ------------------------------

# Buckets are stored in 1-based vectors where slot `i` holds the GF(2)-hash label
# `i - 1`. The output-bucket relation `label_C = label_A ⊻ label_B` is XOR on the
# 0-based labels, so this is the corresponding group operation on the 1-based slots:
# given A-slot `a` and output slot `c`, it returns the complementary B-slot.
@inline _xorbucket(a::Int, c::Int) = ((a - 1) ⊻ (c - 1)) + 1

# Accumulate output bucket `c` = Σ_a A[a]·B[a⊻c] into `d`. Lives in its own
# function so the parallel task body stays small and type-stable. `nonemptyA` is the
# precomputed list of nonempty A-bucket slots (1-based).
function _accumulate_bucket!(d, f::F, c::Int, α, nonemptyA, bucketsA::Buckets, bucketsB::Buckets, maxlength::Int) where {F}
    @inbounds for a in nonemptyA
        bucketB = bucketsB[_xorbucket(a, c)]
        isempty(bucketB) && continue
        bucketA = bucketsA[a]
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

# Body of the parallel loop: compute output bucket `c` into the task-local dict `d`,
# then collect it (with the epsilon cutoff applied here, so the assembly is a plain
# `copyto!`) into this bucket's output slot. `d` is reused across the buckets a task
# handles, so `empty!` resets it (capacity retained, no rehash) and we copy out before
# the dict is reused. Each call writes only slot `c`, so it is data-race free despite
# the non-deterministic task order. The cutoff is branch-free: every entry is written
# at the running index and the index only advances when the entry is kept.
function _fill_bucket!(
        c::Int, d, outkeys::Vector{Vector{P}}, outvals::Vector{Vector{T}},
        f::F, α, β, seed::Bool, nonemptyA, bucketsA::Buckets, bucketsB::Buckets,
        bucketsC, maxlength::Int, docut::Bool, ϵ²::Real
    ) where {P, T, F}
    empty!(d)
    seed && _seed_bucket!(d, bucketsC[c], β)
    _accumulate_bucket!(d, f, c, α, nonemptyA, bucketsA, bucketsB, maxlength)
    ks = Vector{P}(undef, length(d))
    vs = Vector{T}(undef, length(d))
    n = 0
    @inbounds for (p, coeff) in zip(keys(d), values(d))
        ks[n + 1] = p
        vs[n + 1] = coeff
        n += !docut | (abs2(coeff) > ϵ²)
    end
    resize!(ks, n); resize!(vs, n)
    outkeys[c] = ks; outvals[c] = vs
    return nothing
end

function binary_kernel!(
        f::F, C::AbstractOperator, A::AbstractOperator, B::AbstractOperator,
        α::Number, β::Number, s::AbstractBucketStrategy;
        maxlength::Int = 1000, epsilon::Real = eps(real(scalartype(C))),
        scheduler::Scheduler = GreedyScheduler()
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

    # slots (1-based) of the nonempty A buckets, hoisted out of the per-bucket loop
    # so each task skips empty A buckets without rescanning.
    nonemptyA = [a for a in 1:nb if !isempty(bucketsA[a])]

    # per-bucket output-size estimate: the serial heuristic max(|A|,|B|) for the
    # whole output, spread across the `nb` buckets.
    base_hint = cld(max(length(A), length(B)), nb)
    docut = epsilon > 0
    ϵ² = epsilon^2

    # one (keys, values) pair per output bucket; keys are disjoint across buckets
    outkeys = Vector{Vector{P}}(undef, nb)
    outvals = Vector{Vector{T}}(undef, nb)

    # `scheduler` (default `GreedyScheduler()`) spawns ~nthreads tasks that pull
    # output buckets `c`; pass e.g. `DynamicScheduler()` to compare. The reducer dict
    # is held in a `TaskLocalValue` (the same mechanism `@local` uses) — lazily
    # created once per task and reused across the buckets that task handles, so only
    # ~nthreads dicts exist in total.
    dicts = OhMyThreads.TaskLocalValue{UnorderedDictionary{P, T}}() do
        UnorderedDictionary{P, T}(; sizehint = base_hint)
    end
    tforeach(1:nb; scheduler = scheduler) do c
        _fill_bucket!(c, dicts[], outkeys, outvals, f, α, β, seed, nonemptyA, bucketsA, bucketsB, bucketsC, maxlength, docut, ϵ²)
    end

    # Assemble: the per-bucket slices hold disjoint keys and the cutoff is already
    # applied, so this is a single preallocated `copyto!` pass.
    resize!(C, sum(length, outkeys))
    ksC, vsC = keys(C), values(C)
    i = 1
    @inbounds for c in 1:nb
        ks, vs = outkeys[c], outvals[c]
        n = length(ks)
        copyto!(ksC, i, ks, 1, n)
        copyto!(vsC, i, vs, 1, n)
        i += n
    end

    return C
end
