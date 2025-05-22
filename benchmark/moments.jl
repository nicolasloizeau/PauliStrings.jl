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
steps = 20

for p in 14:2:20
    nterms = 2^p
    g["nterms=2^$p"] = @benchmarkable PauliStrings.lanczos($H, $O, $steps, $nterms; keepnorm=true)
end

run(SUITE)
