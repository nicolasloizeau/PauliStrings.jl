using PauliStrings
import PauliStrings as ps
using Printf

function XXZ(N, J, Δ)
    H = ps.Operator(N)
    i = 1
    H += J / 2, "S+", i, "S-", mod1(i + 1, N)
    H += J / 2, "S-", i, "S+", mod1(i + 1, N)
    H += J * Δ, "Sz", i, "Sz", mod1(i + 1, N)
    return ps.OperatorTS1D(H, full=false)
end

M = 3
L = 2 * M + 1

# should find 1 LIOM in this basis - the energy current
support = ps.symmetry_adapted_k_local_basis(L, M; time_reversal=:imag, spin_flip=:even, conserve_magnetization=:yes, translational_symmetry=true)
println("Size of basis: $(length(support))")
H = XXZ(L, 1.0, 0.5)
evals, ops = ps.lioms(H, support; return_all=true)

for i in eachindex(evals)
    @printf("Eigenvalue %2d: %1.5e\n", i, evals[i])
    @printf("||[H, O]|| = %1.5e\n", ps.opnorm(ps.commutator(H, ops[i]), normalize=true))
    println()
end

println("Number of LIOMs found: ", count(x -> abs(x) < 1e-10, evals))