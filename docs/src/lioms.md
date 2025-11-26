# Local Integrals of Motion (LIOMs)

Here we show how to use PauliStrings.jl to construct the local integrals of motion ([`lioms`](@ref)) in the 1D XXZ model. We will focus on reproducing Fig. 1 from [Pawłowski 2025](https://arxiv.org/abs/2505.05882), where this algorithm was first introduced.

TODO: outline of the algorithm

Start by setting up PauliStrings:
```@example lioms
using PauliStrings
import PauliStrings as ps
```

Then we define the XXZ Hamiltonian $$H =\sum_i \frac{J}{2}\left(S^+_i S^-_{i+1} + S^-_i S^+_{i+1}\right) + \Delta S^z_i S^z_{i+1}$$ on a 1D chain with periodic boundary conditions:
```@example lioms
function XXZ(N, J, Δ)
    H = ps.Operator(N)
    i = 1
    H += J / 2, "S+", i, "S-", mod1(i + 1, N)
    H += J / 2, "S-", i, "S+", mod1(i + 1, N)
    H += J * Δ, "Sz", i, "Sz", mod1(i + 1, N)
    return ps.OperatorTS1D(H, full=false)
end
```

## General Pauli basis
First, we construct the LIOMs using the full Pauli basis supported on up to $k$ sites (Fig. 1a of [Pawłowski 2025](https://arxiv.org/abs/2505.05882)). It is constructed from all possible Pauli strings acting on up to $k$ contiguous sites, translated across the entire lattice. Each Pauli string is composed of the operators `X`, `Y`, `Z`, and the identity `1`.

```@example lioms
k = 4
L = 2 * k + 1
J = 1.0
Δ = 0.5

H = XXZ(L, J, Δ)
support = ps.k_local_basis(L, k; translational_symmetry=true)
evals, lioms = ps.lioms(H, support; threshold=1e-10)

# 4 zero eigenvalues, corresponding to 4 LIOMs on support up to k=4 sites
@show evals
```
Using the general Pauli basis on up to $k$ sites, we find exactly k LIOMs as expected.


## Symmetry adapted Pauli basis
We can improve the efficiency of the LIOMs construction by using a symmetry-adapted operator basis (Fig. 1b of [Pawłowski 2025](https://arxiv.org/abs/2505.05882)). To this end, we notice that the XXZ Hamiltonian conserves the total magnetization $S^z_\mathrm{tot} = \sum_i S^z_i$. Therefore, we can expect the LIOMs to also commute with $S^z_\mathrm{tot}$, and so we include only such commuting Pauli strings in our basis. This is most easily done by using the raising and lowering operators `S+` and `S-`, along with `Sz` and the identity `1` as the building blocks of the Pauli strings.

Furthermore, we can split the operator basis into sectors with well-defined behavior under spin-flip symmetry, defined by the operator $F = \prod_i X_i$, $F^\dagger O F = \pm O$. 


TODO: finish this and add complete example with a plot


