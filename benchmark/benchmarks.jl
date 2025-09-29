using PauliStrings
using BenchmarkTools
include("models.jl")
include("functions_to_benchmark.jl")


# here are some physically relevant benchmarks
# system size and precision are choosen so that they take ~1min each on a laptop



SUITE = BenchmarkGroup()

# Lanczos benchmarks
# ------------------
# Each model is benchmarked for 2 different precisions to check convergence, and the output is saved in `figures/`
g = addgroup!(SUITE, "lanczos")
# a 1D translation invariant integrable but hard XXZ model + field with a dynamical symmetry
g["1D XXZ"] = @benchmarkable run_lanczos("XXZh", "XXZo3", 20; precision=18)
# a 1D translation invariant chaotic but easy model with linear growth of lanczos coefficients
g["chaotic_ising"] = @benchmarkable run_lanczos("chaotic_ising", "chaotic_ising_op", 30; precision=18)
# a 1D disordered model, not translation invariant, with many-body localization
g["mbl"] = @benchmarkable run_lanczos("mbl", "Z1", 30; precision=18)
# TODO 2D translation invariant



# Time evolution benchmarks
# -------------------------
g = addgroup!(SUITE, "evolution")
# a 1D chaotic model with next nearest neighbor interactions and diffusive transport
# same model as in https://arxiv.org/pdf/2410.09654 figure 1
# H is translation invariant but we dont use OperatorTS1D here because the inital operator is not
g["1D XXZnnn"] = @benchmarkable run_evolution("XXZnnn", "Z1", 10; precision=14, tmax=5.0, noise=0.05)
# a 2D chaotic model with diffusive transport
g["2D XXZ"] = @benchmarkable run_evolution("XXZ2D", "Z1", 5 * 5; precision=14, tmax=5.0, noise=0.05)



# Moments benchmarks
# ------------------
g = addgroup!(SUITE, "moments")
N = 25
kmax = 10 # number of moments to compute
# 1D translation invariant XXZ using OperatorTS1D
H = models.XXZ(N)
g["XXZ OperatorTS1D"] = @benchmarkable moments(H, kmax; scale=1)
# 1D translation invariant XXZ without using OperatorTS1D
H = Operator(models.XXZ(N))
g["XXZ Operator"] = @benchmarkable moments(H, kmax; scale=1)
# TODO 2D translation invariant


# results = run(SUITE["moments"])
# println(results)
