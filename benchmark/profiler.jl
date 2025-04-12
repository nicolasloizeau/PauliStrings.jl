using BenchmarkTools
using PauliStrings

SUITE = BenchmarkGroup()

# Lanczos benchmarks
# ------------------
# XX hamiltonian: 1d chain with XX + YY interaction
function XX(N)
    H = PauliStrings.Operator(N)
    for j in 1:(N-1)
        H += "X", j, "X", j + 1
        H += "Z", j, "Z", j + 1
    end
    return H
end

# X local operator: X operator on each site
function X(N)
    H = PauliStrings.Operator(N)
    for j in 1:N
        H += "X", j
    end
    return H
end


N = 50
g = addgroup!(SUITE, "lanczos (N=$N)")
H = XX(N)
O = X(N)
steps = 30

p = 20
nterms = 2^p

PauliStrings.lanczos(H, O, steps, nterms; keepnorm=true)

PauliStrings.COM_THREADS[] = 0
bs0 = map(18:24) do p
    nterms = 2^p
    @benchmark PauliStrings.lanczos($H, $O, $steps, $nterms; keepnorm=true)
end

PauliStrings.COM_THREADS[] = 1
bs1 = map(18:24) do p
    nterms = 2^p
    @benchmark PauliStrings.lanczos($H, $O, $steps, $nterms; keepnorm=true)
end

PauliStrings.COM_THREADS[] = 2
bs2 = map(18:24) do p
    nterms = 2^p
    @benchmark PauliStrings.lanczos($H, $O, $steps, $nterms; keepnorm=true)
end

bs0
bs1
bs2

nterms = 2^20
steps = 20
PauliStrings.THREAD_THRESHOLD[] = 2^17
b0, b1, b2 = map(0:2) do nthreads
    PauliStrings.COM_THREADS[] = nthreads
    @benchmark PauliStrings.lanczos($H, $O, $steps, $nterms; keepnorm=true)
end

@profview PauliStrings.lanczos(H, O, steps, nterms; keepnorm=true)



PauliStrings.BLOCK_SIZE[] = 2^20
PauliStrings.lanczos(H, O, steps, nterms; keepnorm=true)
@benchmark PauliStrings.lanczos($H, $O, $steps, $nterms; keepnorm=true)

bs = map(14:2:30) do bsz
    PauliStrings.BLOCK_SIZE[] = 2^bsz
    @benchmark PauliStrings.lanczos($H, $O, $steps, $nterms; keepnorm=true)
end

# for p in 14:2:20
#     nterms = 2^p
#     g["nterms=2^$p"] = @benchmarkable PauliStrings.lanczos($H, $O, $steps, $nterms; keepnorm=true)
# end




