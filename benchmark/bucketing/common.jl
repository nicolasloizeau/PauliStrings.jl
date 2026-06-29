# Shared helpers for the bucketing benchmarks: representative operands built by
# time-evolving a local operator under various models, plus the bucketing
# strategy set and the per-bucket load-balance metrics.

using PauliStrings
import PauliStrings as ps
using Random

include(joinpath(@__DIR__, "..", "models.jl"))

# A local Z on a central site — the seed operator we evolve.
function central_Z(N::Int)
    O = Operator(N)
    O += "Z", cld(N, 2)
    return O
end

# Heisenberg-evolve `O0` under `H`, truncating to `M` strings, for `nsteps` RK4
# steps. Returns a representative O(t) string collection of at most `M` strings.
function evolve_to(H, O0; M::Int, dt::Real = 0.05, nsteps::Int = 30)
    trunc = O -> ps.trim(O, M)
    O = copy(O0)
    for _ in 1:nsteps
        O = ps.trim(rk4(H, O, dt; heisenberg = true, truncation = trunc), M)
    end
    return O
end

# Models, each as a plain `Operator` Hamiltonian on `N` qubits (TS models are
# converted with `Operator(·)`), paired with the central-Z seed operator.
function model_HO(name::AbstractString, N::Int)
    Random.seed!(0)  # mbl draws a random field; keep operands reproducible
    H = if name == "chaotic_ising"
        Operator(models.chaotic_ising(N))
    elseif name == "XXZnnn"
        Operator(models.XXZnnn(N))
    elseif name == "mbl"
        models.mbl(N)
    elseif name == "XXZ2D"
        models.XXZ2D(N)
    else
        error("unknown model $name")
    end
    return H, central_Z(N)
end

# Models to sweep; override with the BENCH_MODELS env var (comma-separated).
const MODELS_1D = Tuple(split(get(ENV, "BENCH_MODELS", "chaotic_ising,XXZnnn,mbl"), ","))

# Bucketing strategies to compare at a given bucket-bit count `b`. `Mixing` is the
# deterministic full-rank matrix that `default_strategy` selects.
function strategies(N::Int, b::Int)
    return (
        ("Folded{$b}", Folded{b}()),
        ("XorVW{$b}", XorVW{b}()),
        ("lowbits{$b}", ps.lowbits_matrix(N, b)),
        ("Mixing{$b}", mixing_matrix(N, b)),
    )
end

# Load-balance metrics
# --------------------

# Per-output-bucket work for C = A·B: work[c+1] = Σ_a |A[a]|·|B[a⊻c]|. This is
# exactly what one thread owning output bucket `c` must compute.
function bucket_work(s, A, B)
    nb = 1 << ps.nbits(s)
    cntA = zeros(Int, nb)
    cntB = zeros(Int, nb)
    for p in keys(A); cntA[ps.bucketindex(s, p) + 1] += 1; end
    for p in keys(B); cntB[ps.bucketindex(s, p) + 1] += 1; end
    work = zeros(Int, nb)
    for c in 0:(nb - 1), a in 0:(nb - 1)
        work[c + 1] += cntA[a + 1] * cntB[(a ⊻ c) + 1]
    end
    return work
end

# Per-bucket occupancy of a single operator's strings.
function bucket_occupancy(s, O)
    nb = 1 << ps.nbits(s)
    cnt = zeros(Int, nb)
    for p in keys(O); cnt[ps.bucketindex(s, p) + 1] += 1; end
    return cnt
end

# Summarise a per-bucket distribution: imbalance = max/mean (1.0 = perfect),
# p99/mean, and fraction of empty buckets.
function balance_stats(w)
    nb = length(w)
    m = sum(w) / nb
    m == 0 && return (imbalance = 1.0, p99 = 1.0, empty = 1.0)
    sorted = sort(w)
    p99 = sorted[clamp(ceil(Int, 0.99 * nb), 1, nb)]
    return (imbalance = maximum(w) / m, p99 = p99 / m, empty = count(==(0), w) / nb)
end
