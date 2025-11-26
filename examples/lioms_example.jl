# Example demonstrating finding LIOMs in the XXZ model. Follows section III.A of https://arxiv.org/abs/2505.05882

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

# should find 1 LIOM in this basis - the energy current Q3
support = ps.symmetry_adapted_k_local_basis(
    L,
    M;
    time_reversal=:imag,
    spin_flip=:even,
    conserve_magnetization=:yes,
    translational_symmetry=true,
)

println("Size of operator basis: $(length(support))") # should be 3
H = XXZ(L, 1.0, 0.75)
evals, ops = ps.lioms(H, support; threshold=Inf) # returns all eigenmodes

for i in eachindex(evals)
    @printf("Eigenvalue λ_%2d: %1.5e\n", i, evals[i])
    @printf("||[H, O]||^2 = %1.5e\n", ps.opnorm(ps.commutator(H, ops[i]), normalize=true)^2)
    @printf("||O|| = %1.5e\n", ps.opnorm(ops[i], normalize=true))
    println("Operator O_$i in Pauli basis:")
    display(ops[i])
    println()
end

println("Number of LIOMs found: ", count(x -> abs(x) < 1e-10, evals))