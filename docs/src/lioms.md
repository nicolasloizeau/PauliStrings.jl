# Local Integrals of Motion (LIOMs)

Here we show how to use PauliStrings.jl to construct the local integrals of motion ([`lioms`](@ref)) in the 1D XXZ model. We will focus on reproducing Fig. 1b and Fig. 2b from [Pawłowski 2025](https://arxiv.org/abs/2505.05882), where this algorithm was first introduced.


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
As a first example, we construct the LIOMs using the full Pauli basis supported on up to $k$ sites. It is constructed from all possible Pauli strings acting on up to $M$ contiguous sites, translated across the entire lattice. Each Pauli string is composed of the operators `X`, `Y`, `Z`, and the identity `1`.

```@example lioms
M = 4
L = 2 * M + 1
J = 1.0
Δ = 0.5

H = XXZ(L, J, Δ)
support = ps.k_local_basis(L, M; translational_symmetry=true)
evals, evecs, lioms = ps.lioms(H, support; threshold=1e-10)

# 4 zero eigenvalues, corresponding to 4 LIOMs on support up to k=4 sites
@show evals
```
Using the general Pauli basis on up to $k$ sites, we find exactly k LIOMs as expected.


## Symmetry adapted Pauli basis
We can improve the efficiency of the LIOMs construction by using a symmetry-adapted operator basis.

- First, we notice that the XXZ Hamiltonian conserves the total magnetization $S^z_\mathrm{tot} = \sum_i S^z_i$. Therefore, we can expect the LIOMs to also commute with $S^z_\mathrm{tot}$, and so we include only such commuting Pauli strings in our basis. This is most easily done by using the raising and lowering operators `S+` and `S-`, along with `Sz` and the identity `1` as the building blocks of the Pauli strings.

- Then, we can split the operator basis into sectors with well-defined behavior under spin-flip symmetry, defined by the operator $F = \prod_i X_i$, $F^\dagger O F = \pm O$.

- And finally, since the XXZ hamiltonian is time-reversal symmetric, we can also split the operator basis into sectors with real matrix elements, $O^s_R = O^s + {O^s}^{\dagger}$, and imaginary matrix elements, $O^s_I = i(O^s - {O^s}^\dagger)$.


```@example symmetry_adapted_lioms
M = 4
support = ps.symmetry_adapted_k_local_basis(
    L,
    M;
    time_reversal=:real,
    spin_flip=:even,
    conserve_magnetization=:yes,
    translational_symmetry=true,
)

evals, evecs, ops = ps.lioms(
    H,
    support;
    threshold=1e-10,
)

# 2 zero eigenvalues, corresponding to 2 LIOMs on support up to M=4 sites
@show evals
```
Using the symmetry-adapted Pauli basis on up to $M$ sites, in the real-even sector, we find exactly $\lfloor M/2 \rfloor$ LIOMs. Remaining $\lceil M/2 \rceil$ LIOMs can be found in the imaginary-even sector.

Repeating the above for $M \in [2, 4, 6]$, we reproduce Fig. 1b and 2b:


![plot](./assets/symmetry_resolved_lioms_xxz.png)

Second panel shows the composition of the $\lambda=0$ LIOMs subspace, in terms of projections on the basis operators. 