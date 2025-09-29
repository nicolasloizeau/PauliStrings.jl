using PauliStrings
import PauliStrings as ps
import PyPlot as plt
using ProgressBars
include("models.jl")

"""
Run lanczos on a given hamiltonian and operator from models.jl, for system size N.
"""
function run_lanczos(hamiltonian::String, operator::String, N::Int; precision=18)
    H = getfield(models, Symbol(hamiltonian))(N)
    O = getfield(models, Symbol(operator))(N)
    steps = N
    plt.cla()
    plt.clf()
    for p in [precision - 2, precision]
        nterms = 2^p
        b = lanczos(H, O, N, nterms)
        plt.plot(b, label="nterms=2^$p")
    end
    plt.legend()
    plt.xlabel("Lanczos step")
    plt.ylabel("Lanczos coefficient")
    plt.savefig("figures/lanczos-$hamiltonian-$operator-$N.png")
end


function evolve(H, O, M, times, noise)
    N = qubitlength(H)
    correlator = []
    dt = times[2] - times[1]
    O0 = copy(O)
    for t in ProgressBar(times)
        push!(correlator, trace_product(O, O0) / 2^N)
        #one step of rk4, keep only M strings, do not discard ztot
        O = rk4(H, O, dt; heisenberg=true, M=M, keep=O0)
        O = add_noise(O, noise * dt)
        O = ps.trim(O, M; keep=O0)
    end
    return real.(correlator)
end


function run_evolution(hamiltonian::String, operator::String, N::Int; precision=12, tmax=10.0, dt=0.05, noise=0.01)
    H = Operator(getfield(models, Symbol(hamiltonian))(N))
    O = getfield(models, Symbol(operator))(N)
    times = 0:dt:tmax
    plt.cla()
    plt.clf()
    for p in [precision - 2, precision]
        correlator = evolve(H, O, 2^p, times, noise)
        plt.loglog(times, correlator, label="nterms=2^$p")
    end
    plt.xlabel("Time")
    plt.ylabel("tr(O(t) O(0)")
    plt.savefig("figures/evolution-$hamiltonian-$operator-$N.png")
end
