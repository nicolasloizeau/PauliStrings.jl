# TS-native Trotter proposal

## Motivation

The current translation-symmetric (TS) Trotter path expands

$$H = \sum_a c_a \sum_r T_r P_a T_r^\dagger$$

with `resum(H)`, then applies a second-order product formula over all translated Pauli terms. For a TS Hamiltonian with `m` representative terms and `N` translations, this creates `O(Nm)` ordinary single-Pauli gates per step. In the mixed-field Ising example with `m = 3`, second-order Strang gives `6N - 1` gates.

The intended TS-native method should instead split over translation orbits, not individual translated terms.

## Orbit-level splitting

Write the Hamiltonian as a sum of orbit Hamiltonians

$$
H = \sum_{a=1}^m H_a^{\mathrm{orb}},
\qquad
H_a^{\mathrm{orb}} = c_a \sum_{r \in G} T_r P_a T_r^\dagger,
$$

where `G` is the translation group. Define the orbit Liouvillian

$$
L_a(O) = \frac{i}{\hbar}[H_a^{\mathrm{orb}}, O].
$$

A second-order TS-native Trotter step is then Strang splitting over orbit Liouvillians:

$$
S_2(\Delta t)
=
e^{\frac{\Delta t}{2} L_1}
e^{\frac{\Delta t}{2} L_2}
\cdots
e^{\Delta t L_m}
\cdots
e^{\frac{\Delta t}{2} L_2}
e^{\frac{\Delta t}{2} L_1}.
$$

This reduces the number of product-formula factors from `O(Nm)` ordinary gates to `O(m)` TS orbit flows.

## Pauli-orbit propagation rule

Let `[Q]` denote the TS Pauli orbit represented by a Pauli string `Q`:

$$
[Q] = \sum_x T_x Q T_x^\dagger.
$$

For one Hamiltonian orbit `[P]`, expand

$$
[[P], [Q]]
= \sum_{r,x} [T_r P T_r^\dagger, T_x Q T_x^\dagger].
$$

Using translation covariance and the relative shift `\delta = r - x`,

$$
[[P], [Q]]
= \sum_x T_x
\left(
  \sum_\delta [T_\delta P T_\delta^\dagger, Q]
\right)
T_x^\dagger.
$$

Thus one orbit Liouvillian acts as a sparse graph on TS Pauli orbits:

```text
input node:  [Q]
edge label:  relative shift δ
output node: canon(T_δ P T_δ† · Q)
weight:      Pauli commutator phase × Hamiltonian coefficient
```

All relative shifts that fold to the same canonical TS representative must be accumulated. This multiplicity is essential because `OperatorTS` represents unnormalized orbit sums.

In code, this graph rule is exposed by:

```julia
orbit_edges(Ha, q)
orbit_liouvillian(Ha, O)
```

where `orbit_edges` gives the outgoing graph edges from a single TS Pauli orbit `q`, and `orbit_liouvillian` sums those edges over an operator `O`.

## Relation to RK4

RK4 approximates the full flow with a polynomial in

$$L_H = \sum_a L_a.$$

The TS-native Trotter proposal instead applies a product of orbit-flow exponentials:

$$e^{\tau L_a}$$

for individual orbit Hamiltonians. This is not RK4 if the orbit flows are treated as exponential actions, e.g. by Krylov or local component exponentials. A Taylor approximation to each orbit flow is useful as a diagnostic but should not be considered the final Trotter primitive.

## Implementation strategy

1. Keep `TrotterTS` as Strang splitting over orbit Liouvillians.
2. Use `orbit_edges` as the generic Pauli-orbit graph generator.
3. Apply orbit exponentials by matrix-free Krylov initially.
4. Investigate local connected components of the orbit graph. If components are small, build local matrices and apply dense exponentials per component.
5. Avoid target-specific assumptions such as 1D nearest-neighbor structure or commuting MFIM orbits in the core method.

## Attempts tried during exploration

We tried a cheap orbit-Strang variant where each orbit flow was replaced by a second-order Taylor approximation, `O + τ L_a O + τ² L_a² O / 2`. This was useful diagnostically and showed that orbit-level splitting can reduce the expanded-gate work, but it is too close in spirit to a commutator-based integrator and should not be considered the final TS Trotter primitive. We also prototyped local connected-component exponentials on the Pauli-orbit graph; small sampled components suggested theoretical promise, but the naive implementation rebuilt components and dense matrices every step, causing excessive allocation and runtime. A 1D active-shift specializer briefly improved hot-loop timings, but it was removed because the core method should remain generic and not be tuned to the target MFIM example.

## Validation goals

For each backend:

1. **Liouvillian correctness**

   ```julia
   orbit_liouvillian(Ha, O) ≈ 1im * commutator(Ha, O) / hbar
   resum(orbit_liouvillian(Ha, O)) ≈ 1im * commutator(resum(Ha), resum(O)) / hbar
   ```

2. **Orbit-flow correctness** on small systems against dense exact evolution under `resum(Ha)`.

3. **Second-order global scaling** of the full Strang method against dense exact evolution.

4. **Work accounting**, not just wall time:
   - expanded ordinary gates versus orbit factors,
   - number of graph edges,
   - component sizes if using component exponentials,
   - allocations,
   - final and intermediate term counts.

## References

- H. F. Trotter, “On the product of semi-groups of operators,” *Proc. AMS* 10, 545–551 (1959).
- G. Strang, “On the Construction and Comparison of Difference Schemes,” *SIAM J. Numer. Anal.* 5, 506–517 (1968).
- N. Hatano and M. Suzuki, “Finding Exponential Product Formulas of Higher Orders,” arXiv:math-ph/0506007.
- A. M. Childs, Y. Su, M. C. Tran, N. Wiebe, and S. Zhu, “A Theory of Trotter Error,” arXiv:1912.08854.
- M. S. Rudolph et al., “Pauli Propagation: A Computational Framework for Simulating Quantum Systems,” arXiv:2505.21606. This is useful context for Pauli-string propagation, though PauliPropagation.jl depends on PauliStrings.jl and is not expected to contain this TS-native implementation.
