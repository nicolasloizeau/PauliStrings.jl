
# Example comparing Trotter and RK4 evolution of the mixed field Ising model, both using
# the translation-symmetric `OperatorTS` representation. Both should agree.


using PauliStrings
import PyPlot as plt

# Translation invariant mixed field Ising model
function MFIM(N, h)
    H = Operator(N)
    H += -h, "X", 1
    H += -h/2, "Z", 2
    H += "Z", 1, "Z", 2
    return OperatorTS{(N,)}(H)
end

# Total X operator
function Xtot(N)
    H = Operator(N)
    H += "X", 1
    return OperatorTS{(N,)}(H)
end

# init the ising model on 32 qubits
N = 32
H = MFIM(N, 0.5)
O0 = Xtot(N)
dt = 0.1
times = 0:dt:5


fout(O) = trace_product(O0, O) / 2^N

for M in (12, 14)
    truncation(o) = trim(o, 2^M)
    # Trotter and RK4 evolution, both on the translation-symmetric OperatorTS
    res_trotter = evolve(H, O0, times; method=Trotter(), fout=fout, truncation=truncation)
    res_rk4 = evolve(H, O0, times; method=RK4(), fout=fout, truncation=truncation)
    plt.plot(times, real.(res_trotter.history), label="Trotter (M=$M)", lw=2)
    plt.plot(times, real.(res_rk4.history), "--", label="RK4 (M=$M)", lw=2)
end


plt.legend()
plt.xlabel("Time")
plt.ylabel("⟨Xtot(t)⟩")
plt.savefig("./evolve_trotter_ts.png", bbox_inches="tight")
plt.show()
