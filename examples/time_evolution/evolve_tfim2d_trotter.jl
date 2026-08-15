
# Same as `evolve_tfim2d.jl` (2D transverse field Ising model at the quantum critical point,
# autocorrelation tr(Ztot(0)*Ztot(t))/2^N), but evolved with Trotter gates instead of RK4.


using PauliStrings
import PyPlot as plt

# Translation invariant 2D transverse field Ising model on an L x L periodic grid
function TFIM2D(L, J, h)
    H = Operator(L * L)
    H += -J * string_2d(("Z", 1, 1, "Z", 2, 1), L, L) # horizontal bond
    H += -J * string_2d(("Z", 1, 1, "Z", 1, 2), L, L) # vertical bond
    H += -h * string_2d(("X", 1, 1), L, L)            # transverse field
    return OperatorTS{(L, L)}(H)
end

# Total Z operator
function Ztot(L)
    O = Operator(L * L)
    O += "Z", 1
    return OperatorTS{(L, L)}(O)
end

# 8x8 grid (64 qubits) at the 2D quantum critical point
L = 8
N = L * L
dt = 0.05

# We don't know how to Trotterize the translation-symmetric representation efficiently,
# so we `resum` the OperatorTS into a plain (non-symmetric) Operator before evolving.
H = resum(TFIM2D(L, 1, 3.04428))
O0 = resum(Ztot(L))

#function we want to observe at each time step
fout(O) = real(trace_product(O0, O)) / 2.0^N

times = 0:dt:2
# evolve for different truncation values, plot the results
for M in [14, 16]
    truncation(o) = trim(o, 2^M)
    res = evolve(H, O0, times; method=Trotter(), fout=fout, truncation=truncation)
    plt.plot(times, res.history, label="#strings = 2^$(M)")
end
plt.legend()
plt.xlabel("Time")
plt.ylabel("tr(Ztot(0)Ztot(t))/2^N")
plt.savefig("./evolve_tfim2d_trotter.png", bbox_inches="tight")
plt.show()
