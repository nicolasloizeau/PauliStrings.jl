# Example of evolving a chaotic spin chain with dissipation.
# For a translation invariant example, see `evolve_tfim.jl`


using PauliStrings
import PyPlot as plt

# Some chaotic 1D spin chain
function chaotic_chain(N::Int)
    H = Operator(N)
    # XX interractions with periodic bc
    for j in 1:N
        H += "X", j, "X", mod1(j+1, N)
    end
    # fields
    for j in 1:N
        H += -1.05, "Z", j
        H += 0.5, "X", j
    end
    return H
end


# intit the ising model on 32 qubits
N = 32
H = chaotic_chain(N)
O0 = Operator(N) + ("Z", 1) # X on site 1
dt = 0.02

# depolarizing noise with strength 0.05
dissipation(O, dt) = add_noise(O, 0.05*dt)

#function we want to observe at each time step
fout(O) = trace_product(O0, O) / 2^N

times = 0:dt:5
# evolve for different truncation values, plot the results
for M in [10, 12, 14]
    truncation(o) = trim(o, 2^M)
    res = evolve(H, O0, times; method = RK4(), fout=fout, dissipation=dissipation, truncation=truncation)
    plt.loglog(times, res.history, label="#strings = 2^$(M)")
end
plt.legend()
plt.xlabel("Time")
plt.ylabel("⟨Xtot(t)⟩")
plt.show()
