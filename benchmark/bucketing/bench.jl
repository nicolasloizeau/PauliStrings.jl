# Single-process timing benchmark of the bucketed `binary_kernel!` for one thread
# count, on representative time-evolved operands. Run via `run_sweep.jl`, which
# launches this once per thread count. Results append as CSV rows to ARGS[1].
#
#   julia --project=benchmark -t <N> benchmark/bucketing/bench.jl <results.csv>
#
# Tunables (env): BENCH_SECONDS (per cell, default 3), BENCH_BITS (bucket bits,
# default "8"), BENCH_N (1D qubits), BENCH_M (target #strings).

include(joinpath(@__DIR__, "common.jl"))
using BenchmarkTools
using Printf

const SECONDS = parse(Float64, get(ENV, "BENCH_SECONDS", "3"))
const BITS = parse.(Int, split(get(ENV, "BENCH_BITS", "8"), ","))
const N1D = parse(Int, get(ENV, "BENCH_N", "18"))
const M = parse(Int, get(ENV, "BENCH_M", "16384"))

results_path = get(ARGS, 1, joinpath(@__DIR__, "bucketing_results.csv"))

# representative operands per model: (H, O(t))
function operands()
    out = Tuple{String,Int,Any,Any}[]
    for name in MODELS_1D
        H, O0 = model_HO(name, N1D)
        push!(out, (name, N1D, H, evolve_to(H, O0; M = M)))
    end
    return out
end

# strategies to time: Serial baseline + the best-mixing families across bit counts
function timing_strategies(N)
    out = Tuple{String,Any}[("Serial", Serial())]
    for b in BITS
        push!(out, ("Folded{$b}", Folded{b}()))
        push!(out, ("Mixing{$b}", mixing_matrix(N, b)))  # the default family
    end
    return out
end

# the two kernel shapes that matter: the RK4 commutator [H,O], and a dense O*O.
kernels(H, O) = (("commutator_HO", () -> commutator(H, O)),
                 ("sqr_OO", () -> O * O))

function run!(io)
    nt = Threads.nthreads()
    for (model, N, H, O) in operands()
        for (kname, _) in kernels(H, O)
            for (sname, s) in timing_strategies(N)
                # build the closure for this strategy
                call = kname == "commutator_HO" ?
                    (() -> ps.binary_kernel(ps.commutator, H, O; strategy = s)) :
                    (() -> ps.binary_kernel(ps.prod, O, O; strategy = s))
                len = length(call())
                bm = @benchmark $call() seconds = SECONDS
                t = minimum(bm).time / 1e9
                mem = minimum(bm).memory
                alloc = minimum(bm).allocs
                @printf(io, "%d,%s,%s,%s,%d,%d,%.6e,%d,%d,%d\n",
                    nt, model, kname, sname, N, M, t, mem, alloc, len)
                flush(io)
                @printf("[t=%2d] %-12s %-14s %-16s  %.3f ms\n", nt, model, kname, sname, t * 1e3)
            end
        end
    end
end

open(results_path, "a") do io
    run!(io)
end
