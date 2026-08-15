
# Example of evolving the 2D transverse field Ising model at the quantum critical point.
# H = -J sum_<ij> Zi Zj - h sum_i Xi on an L x L periodic grid, with h/J = 3.04428
# We compute the autocorrelation tr(Ztot(0)*Ztot(t))/2^N.
# For the 1D version, see `evolve_tfim.jl`


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
H = TFIM2D(L, 1, 3.04428)
O0 = Ztot(L)
dt = 0.02

#function we want to observe at each time step
fout(O) = real(trace_product(O0, O)) / 2.0^N

times = 0:dt:2
# evolve for different truncation values, plot the results
for M in [10, 12, 14]
    truncation(o) = trim(o, 2^M)
    res = evolve(H, O0, times; method=RK4(), fout=fout, truncation=truncation)
    plt.plot(times, res.history, label="#strings = 2^$(M)")
end
plt.legend()
plt.xlabel("Time")
plt.ylabel("tr(Ztot(0)Ztot(t))/2^N")
plt.savefig("./evolve_tfim2d.png", bbox_inches="tight")
plt.show()
