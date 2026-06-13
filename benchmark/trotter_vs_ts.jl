using PauliStrings
using LinearAlgebra
using Printf

function MFIM(N, h)
    H = Operator(N)
    H += -h, "X", 1
    H += -h / 2, "Z", 2
    H += "Z", 1, "Z", 2
    return OperatorTS{(N,)}(H)
end

function Xtot(N)
    O = Operator(N)
    O += "X", 1
    return OperatorTS{(N,)}(O)
end

# Wrapper to suppress progress bars during benchmarking
function quiet_evolve(args...; kwargs...)
    return evolve(args...; kwargs...)
end

function run_benchmark(N, M_terms, steps, dt, tols)
    h = 0.5
    H = MFIM(N, h)
    O0 = Xtot(N)
    times = 0:dt:(steps*dt)

    # Create explicit non-TS versions for the baseline
    H_exp = Operator(H)
    O0_exp = Operator(O0)

    # Mathematically fair truncation functions
    truncation(O) = trim(O, M_terms * N)   # For Expanded Trotter (ordinary strings)
    truncationTS(O) = trim(O, M_terms)   # For TrotterTS and RK4 (TS orbits)

    # Observable: X magnetization, normalized
    fout(O) = real(trace_product(O0, O)) / (2.0^N)

    @printf("\n\n")
    @printf(" Benchmarking N=%d, M=%d terms, steps=%d, dt=%.2f\n", N, M_terms, steps, dt)
    @printf("=========================================================================\n\n")

    # 1. WARMUP
    _ = quiet_evolve(H_exp, O0_exp, 0:dt:dt; method=Trotter(order=2), truncation=truncation, fout=fout)
    _ = quiet_evolve(H, O0, 0:dt:dt; method=RK4(), truncation=truncationTS, fout=fout)
    for tol in tols
        _ = quiet_evolve(H, O0, 0:dt:dt; method=TrotterTS(order=2, componenttol=tol), truncation=truncationTS, fout=fout)
    end

    # 2. RK4 (TS-Native Ground Truth)
    GC.gc()
    stats_rk4 = @timed quiet_evolve(H, O0, times; method=RK4(), truncation=truncationTS, fout=fout)
    res_true = stats_rk4.value.history
    len_rk4 = length(stats_rk4.value.final)

    # 3. BASELINE (Standard Expanded Trotter)
    GC.gc()
    stats_base = @timed quiet_evolve(H_exp, O0_exp, times; method=Trotter(order=2), truncation=truncation, fout=fout)
    res_base = stats_base.value.history
    len_base = length(stats_base.value.final)

    @printf("%-18s | %8s | %9s | %9s | %10s | %8s\n", "Method", "Time (s)", "Alloc(GB)", "Final Len", "Rel. Error", "Speedup")
    @printf("-------------------|----------|-----------|-----------|------------|----------\n")

    err_base = norm(res_base .- res_true) / norm(res_true)
    @printf("%-18s | %8.3f | %9.3f | %9d | %10.2e | %8s\n", "Expanded Trotter", stats_base.time, stats_base.bytes / 1024^3, len_base, err_base, "1.00x")

    speedup_rk4 = stats_base.time / stats_rk4.time
    @printf("%-18s | %8.3f | %9.3f | %9d | %10s | %7.2fx\n", "RK4 (TS)", stats_rk4.time, stats_rk4.bytes / 1024^3, len_rk4, "-", speedup_rk4)

    # 4. TROTTER TS (Component Batching)
    for tol in tols
        method = TrotterTS(order=2, componenttol=tol)

        GC.gc()
        stats_ts = @timed quiet_evolve(H, O0, times; method=method, truncation=truncationTS, fout=fout)

        res_ts = stats_ts.value.history
        len_ts = length(stats_ts.value.final)

        # Calculate relative error over the time history vs the ground truth (RK4)
        err = norm(res_ts .- res_true) / norm(res_true)
        speedup = stats_base.time / stats_ts.time

        @printf("TrotterTS tol=%-4g | %8.3f | %9.3f | %9d | %10.2e | %7.2fx\n", tol, stats_ts.time, stats_ts.bytes / 1024^3, len_ts, err, speedup)
    end
end

# Ensure code is compiled before timing actual benchmarks
println("Compiling and warming up routines...")
run_benchmark(6, 256, 2, 0.1, [0.99])

# Run actual scaling benchmarks
run_benchmark(16, 2^11, 20, 0.1, [0.99, 0.999, 0.9999, 1.0])
run_benchmark(32, 2^11, 20, 0.1, [0.99, 0.999, 0.9999, 1.0])

# Omit 1.0 for N=64 to avoid extremely long runtimes on exact exponential evaluation
run_benchmark(64, 2^11, 20, 0.1, [0.99, 0.999, 0.9999])
