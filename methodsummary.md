---
title: "TrotterTS Implementation"
---

# Overview

The TrotterTS method is a translation-symmetric orbit-level product formula for time evolution under translation-invariant Hamiltonians. It replaces the O(Nm) ordinary single-Pauli gates in the baseline Trotter approach with O(m) TS orbit-level operations, where m is the number of Hamiltonian representative terms and N is the system size.

***

# Mathematical Method

## Introduction

### The Evolution Equation (Heisenberg Picture)

The real-time evolution of a translation-symmetric operator $O(t)$ under a translation-invariant Hamiltonian $H$ is governed by the **Heisenberg equation of motion**:

$$\frac{d}{dt} O(t) = \frac{i}{\hbar} [H, O(t)]$$

where $[\cdot, \cdot]$ is the quantum commutator, and $H$ and $O(t)$ are elements of the translation-symmetric operator algebra represented as:

$$H, O(t) \in \mathcal{A}_{\text{TS}}$$

### The Formal Exact Solution

The formal exact solution to this differential equation is given by unitary conjugation under the time-evolution operator:

$$O(t) = e^{\frac{i}{\hbar} H t} O(0) e^{-\frac{i}{\hbar} H t}$$

Equivalently, this can be written in superoperator form by defining the **Liouvillian superoperator** $L_H \equiv \frac{i}{\hbar} [H, \cdot]$:

$$O(t) = e^{t L_H} O(0)$$

### Orbits

Let $G$ be the discrete translation group of the lattice and $Q$ be a representative Pauli string.

The **orbit set** of $Q$ under $G$ is:

$$\text{orbit}(Q) = \{ T_g Q T_g^\dagger \mid g \in G \}$$

In a translation-symmetric operator representation, the **orbit operator** $[Q]$ is represented as the sum of its unique translates:

$$[Q] = \sum_{g \in G / \text{Stab}(Q)} T_g Q T_g^\dagger$$

where $\text{Stab}(Q) = \{ g \in G \mid T_g Q T_g^\dagger = Q \}$ is the stabilizer subgroup of translations that leave $Q$ invariant.

## Orbit-Level Splitting

The Hamiltonian is decomposed into translation orbits rather than individual translated terms:

$$
H = \sum_{a=1}^m H_a^{\text{orb}}, \quad \text{where} \quad H_a^{\text{orb}} = c_a \sum_{r \in G} T_r P_a T_r^{\dagger}
$$

Here $G$ is the translation group, $P_a$ are representative Pauli strings, and $T_r$ denotes translation by $r$.

Operator $O$ undergoes a similar treatment, with its orbit-paulis denoted by $Q$ with necessary subscripts.

## Orbit Liouvillian

For each orbit, define the Liouvillian superoperator:

$$L_a(O) = \frac{i}{\hbar}[H_a^{\text{orb}}, O]$$

## Second-Order Strang Splitting

> Remark: For now, only order 2 trotter is implemented. Higher order update complicate the Pauli graph propogation below

The time evolution is approximated by a Strang product over orbit Liouvillians:

$$S_2(\Delta t) = e^{\frac{\Delta t}{2} L_1} e^{\frac{\Delta t}{2} L_2} \cdots e^{\Delta t L_m} \cdots e^{\frac{\Delta t}{2} L_2} e^{\frac{\Delta t}{2} L_1}$$

This reduces the number of product-formula factors from $O(Nm)$ to $O(m)$.

Now, the complexity of the problem is shifted to applying the Liouvillian of the Pauli orbits, for which there exists no closed form solution, unlike the Pauli string case, where the splitting lemma can be used for pauli rotation updates. In the special case that translated Pauli strings commute, $[P_a, P_b] = 0$, the exponential factorizes, allowing us to apply the Pauli splitting lemma to the active sites. An initial attempt did try to apply a Krylov-Arnoldi update to implement the exponential update, but the cost of building the commutators, inner products, high level OperatorTS allocations for every orbit flow lead to a much worse performance profile. This lead to the next attempt.

TODO: add data: krylov

## Pauli-Orbit Graph Propagation

The action of $L_a$ on a TS Pauli orbit $[Q] = \sum_x T_x Q T_x^\dagger$ generates transitions to other TS Pauli orbits via the commutator:

$$[[P], [Q]] = \sum_{r,x} [T_r P T_r^\dagger, T_x Q T_x^\dagger] = \sum_x T_x \left( \sum_\delta [T_\delta P T_\delta^\dagger, Q] \right) T_x^\dagger$$

where $\delta = r - x$ is the relative shift. This defines a sparse graph on TS Pauli orbits:

- Nodes: TS Pauli orbits $[Q]$
- Edges: Relative shifts $\delta$ with weights from Pauli commutator phases

The reason this is sparse is because both, the Hamiltonian and Operator PauliStrings tend to be local, allowing us to only check active sites and skip the rest. This is the second major optimization that we apply over the full update method currently used in the codebase. In fact, this also can speed up TS-RK4, since the current commutator calc is applied on by expanding all sites, which is not necessary. [TODO: verify]

## Connected Component Exponentiation

The Pauli-orbit graph naturally decomposes into disjoint connected components. These components are precisely the efficient update primitives that exists in a TS native update. They are typically small and update by a single dense matrix exponentiation.

To apply the update for each component:

- Build the Liouvillian matrix $A$ where $A_{i,j} = \text{coefficient of } [Q_i] \text{ in } \frac{i}{\hbar}[H_a^{\text{orb}}, [Q_j]]$
- Apply the matrix exponential: $\text{coeffs\_out} = \exp(dt A) \times \text{coeffs\_in}$

The specific construction is used so $A$ is real. Also, since the number of pauli strings in a component is small, $A$ is also small.

TODO: add data: component size

## Component Weight Freezing

Now the main bottleneck is the number of components that we have to update. This tends to be quite large, and so the only way is to truncate the components to only operate only on high weight components. This mirrors the expansion of PauliString count in the non-TS level, to which the current code applies a dissipation to reduce.

TODO: add data: component count

Components are sorted by weight $w = \sum |\text{coeffs}|^2$. Only components contributing to the top `componenttol` fraction of total weight are evolved; the rest are "frozen" (carried forward unchanged). This exploits the heavy-tailed distribution of component weights. We can also apply a plain truncation, or a "drop" method, but this increased error and norm more significantly than unevolved components.

## Signature Batching (Optimization)

Whereas component truncation does increase the time sufficiently, signature batching allows us to update more components with barely any additional computation.

Signature based component batching requires that multiple disjoint components share the exact same matrix representation under a consistent basis ordering. This allows the joint coefficient matrix $C(t) = [\vec{c}^{(1)}(t), \vec{c}^{(2)}(t), \dots, \vec{c}^{(M)}(t)] \in \mathbb{C}^{n \times M}$ to evolve collectively using a single propagator:

$$C(t) = e^{t A} C(0)$$

There are two primary mathematical formulations to establish the identical matrix representation $A$ needed for this batched propagation:

### 1. Exact Matrix Signature Match
This method groups components whose algebraic representations are identical under their default ordering.

* **Mathematical Condition:** Two components $C_1$ and $C_2$, with ordered bases $B_1 = \{u_1, \dots, u_n\}$ and $B_2 = \{v_1, \dots, v_n\}$, are batched if their matrix representations $A_1$ and $A_2$ satisfy:
  $$A_1 = A_2 \quad \iff \quad (A_1)_{i,j} = (A_2)_{i,j} \quad \forall i, j \in \{1, \dots, n\}$$
* **Symmetry:** This groups components that are translationally equivalent and generated via identical deterministic graph traversals.


### 2. Isomorphic Permuted Match (Canonical Labeling)
This method generalizes batching to isomorphic components that are ordered differently in memory.

* **Mathematical Condition:** Two components $C_1$ and $C_2$ are isomorphic if there exists a permutation matrix $P$ such that their arbitrary matrix representations $A_1$ and $A_2$ satisfy:
  $$A_2 = P A_1 P^T$$
* **Canonicalization Map:** To enable Level 2 batching, we define a canonical labeling function $\pi$ that maps any isomorphic matrix to a unique canonical representative:
  $$\pi(A_1) = \pi(A_2) \equiv A_{\text{canon}}$$
* **The Batching Equation:** Let $P_k$ be the permutation matrix that maps component $k$ to its canonical basis. We permute the initial coefficients before stacking them:
  $$C_{\text{canon}}(0) = \left[ P_1 \vec{c}^{(1)}(0), P_2 \vec{c}^{(2)}(0), \dots, P_M \vec{c}^{(M)}(0) \right]$$
  We then propagate the entire canonicalized batch in one step:
  $$C_{\text{canon}}(t) = e^{t A_{\text{canon}}} C_{\text{canon}}(0)$$
  The physical coefficients are retrieved by applying the inverse permutation:
  $$\vec{c}^{(k)}(t) = P_k^T \vec{c}^{(k)}_{\text{canon}}(t)$$


The code currently implements only Exact signature matching, since the cost of checking for permutations may not be worth it, and the current code is already quite efficient. We might have to study which of these occurs more, 
and decide accordingly, or allow users to pick.

***

Programmatic Method

New Types (src/evolution.jl)

struct TrotterTS <: AbstractEvolutionMethod
    order::Int           # 1 (Lie) or 2 (Strang)
    componenttol::Float64  # fraction of component weight to evolve
end

New Module (src/orbit.jl)

┌───────────────────────────────────────────────┬─────────────────────────────────────────────────────────────┐
│                   Function                    │                           Purpose                           │
├───────────────────────────────────────────────┼─────────────────────────────────────────────────────────────┤
│ _coord(site, Ls, dim)                         │ Convert linear site index to coordinate in dimension dim    │
├───────────────────────────────────────────────┼─────────────────────────────────────────────────────────────┤
│_linear_site(coords, Ls)                      │ Convert coordinates back to linear index                    │
├───────────────────────────────────────────────┼─────────────────────────────────────────────────────────────┤
│_active_shift_tuple(site_p, site_q, Ls, Ps)   │ Compute valid relative shifts between two active sites      │
├───────────────────────────────────────────────┼─────────────────────────────────────────────────────────────┤
│ orbit_edges(A, q)                             │ Generate outgoing edges from TS Pauli orbit q under orbit A │
├───────────────────────────────────────────────┼─────────────────────────────────────────────────────────────┤
│ orbit_liouvillian(A, B)                       │ Full orbit commutator: i[A,B]/ħ                             │
├───────────────────────────────────────────────┼─────────────────────────────────────────────────────────────┤
│ _orbit_component(Ha, seed)                    │ BFS to find connected component containing seed             │
├───────────────────────────────────────────────┼─────────────────────────────────────────────────────────────┤
│ _orbit_component_matrix(Ha, component, index) │ Build Liouvillian matrix for a component                    │
├───────────────────────────────────────────────┼─────────────────────────────────────────────────────────────┤
│_component_plan(cache, Ha, seed)              │ Get or build cached plan for component                      │
├───────────────────────────────────────────────┼─────────────────────────────────────────────────────────────┤
│_component_exp(plan, dt)                      │ Get or compute exp(dt × matrix) with caching                │
├───────────────────────────────────────────────┼─────────────────────────────────────────────────────────────┤
│_component_signature(A)                       │ Hashable signature for grouping identical matrices          │
├───────────────────────────────────────────────┼─────────────────────────────────────────────────────────────┤
│_orbit_terms(H)                               │ Split Hamiltonian into single-term orbit operators          │
└───────────────────────────────────────────────┴─────────────────────────────────────────────────────────────┘

Evolution Pipeline (src/evolution.jl)

_orbit_flow(Ha, O, dt, cache, method, hbar, truncation)
Main orbit evolution function:

1. Build coefficient lookup from O
2. For each TS Pauli string, find its connected component and collect weights
3. Sort components by weight (descending)
4. Group components by matrix signature
5. For each signature group:

- Compute E = exp(dt × A) once
- Apply coeffs_out = E × coeffs_in via batched matrix-matrix multiply

1. Freeze low-weight components (carry forward unchanged)
2. Apply truncation

_trotterts_step(Hterms, caches, O, dt, method, hbar, truncation)
Applies Strang splitting over orbit terms:

- Order 1: sequential exp(dt × L_j) for each orbit
- Order 2: palindromic exp(dt/2 × L_1) ... exp(dt × L_m) ... exp(dt/2 × L_1)

_evolve(::TrotterTS, H::OperatorTS, O::OperatorTS, tspan; ...)
Main evolution loop:

1. Split H into orbit terms via _orbit_terms
2. Initialize per-orbit caches
3. For each timestep:

- Apply _trotterts_step
- Apply dissipation and truncation
- Save to history

Key Data Structures

struct _OrbitComponentPlan{P}
    component::Vector{P}           # Pauli strings in component
    index::Dict{P,Int}             # string → position mapping
    matrix::Matrix{ComplexF64}     # Liouvillian matrix
    exp_cache::Dict{Float64,Matrix}  # dt → exp(dt×matrix) cache
end

mutable struct _OrbitFlowCache{P,O}
    plans::Dict{P,_OrbitComponentPlan}  # component plans by seed
    fallback::Dict{P,O}                 # fallback storage
end

Comparison to Baseline

┌──────────────────────┬──────────────────────────────────┬───────────────────────────────┐
│        Aspect        │         Baseline Trotter         │           TrotterTS           │
├──────────────────────┼──────────────────────────────────┼───────────────────────────────┤
│ Hamiltonian          │ resum(H) → O(Nm) terms           │ m orbit representatives       │
├──────────────────────┼──────────────────────────────────┼───────────────────────────────┤
│ Gates per step       │ 2 × length(resum(H)) - 1 = O(Nm) │ 2 × length(Hterms) - 1 = O(m) │
├──────────────────────┼──────────────────────────────────┼───────────────────────────────┤
│ State representation │ Ordinary Operator                │ OperatorTS (compressed)       │
├──────────────────────┼──────────────────────────────────┼───────────────────────────────┤
│ Core operation       │ Single-Pauli rotation            │ Orbit Liouvillian exponential │
├──────────────────────┼──────────────────────────────────┼───────────────────────────────┤
│ Memory               │ O(N × terms)                     │ O(terms) representatives      │
└──────────────────────┴──────────────────────────────────┴───────────────────────────────┘

***
Files Modified/Added

┌─────────────────────┬───────────────────────────────────────────────────────────────────────────────────────┐
│        File         │                                        Changes                                        │
├─────────────────────┼───────────────────────────────────────────────────────────────────────────────────────┤
│ src/orbit.jl        │ New file (166 lines): orbit graph operations, component finding, matrix exponentials  │
├─────────────────────┼───────────────────────────────────────────────────────────────────────────────────────┤
│ src/evolution.jl    │ Added (133 lines): TrotterTS type, _orbit_flow, _trotterts_step, _evolve(::TrotterTS) │
├─────────────────────┼───────────────────────────────────────────────────────────────────────────────────────┤
│ src/PauliStrings.jl │ Export orbit_edges, orbit_liouvillian                                                 │
├─────────────────────┼───────────────────────────────────────────────────────────────────────────────────────┤
│ test/evolution.jl   │ Added TrotterTS test case                                                             │
├─────────────────────┼───────────────────────────────────────────────────────────────────────────────────────┤
│ benchmark/*.md      │ Theory notes, slowdown analysis, design rationale                                     │
└─────────────────────┴───────────────────────────────────────────────────────────────────────────────────────┘

***
Design Principles

1. Generic method: No assumptions about 1D geometry, nearest-neighbor structure, or commuting orbits
2. Orbit-level, not gate-level: Operates on TS representatives, not expanded translated terms
3. Exact within components: Uses matrix exponentials on connected components, not Taylor approximations
4. Adaptive precision: componenttol parameter trades accuracy for speed by freezing low-weight components
5. Batched computation: Signature grouping reduces redundant matrix exponentials
