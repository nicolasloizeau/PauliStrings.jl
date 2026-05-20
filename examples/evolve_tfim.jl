
# Example of evolving the translation invariant transverse field Ising model with dissipation.
# For a non-translation invariant example, see `evolve_chaotic.jl`


using PauliStrings
import PyPlot as plt

# Translation invariant transverse field Ising model
function TFIM(N, h)
    H = Operator(N)
    H += -h, "X", 1
    H += "Z", 1, "Z", 2
    return OperatorTS{(N,)}(H)
end

# Total X operator
function Xtot(N)
    H = Operator(N)
    H += "X", 1
    return OperatorTS{(N,)}(H)
end

# intit the ising model on 32 qubits
N = 32
H = TFIM(N, 0.3)
O0 = Xtot(N)
dt = 0.01

# depolarizing noise with strength 0.05
dissipation(O, dt) = add_noise(O, 0.05*dt)

#function we want to observe at each time step
fout(O) = trace_product(O0, O) / 2^N

times = 0:dt:10
# evolve for different truncation values, plot the results
for M in [10, 12, 14]
    truncation(o) = trim(o, 2^M)
    res = evolve(H, O0, times; method = RK4(), fout=fout, dissipation=dissipation, truncation=truncation)
    plt.plot(times, res.history, label="#strings = 2^$(M)")
end
plt.legend()
plt.xlabel("Time")
plt.ylabel("⟨Xtot(t)⟩")
plt.show()
