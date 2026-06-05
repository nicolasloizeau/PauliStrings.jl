# TS Trotter attempts flow

Internal tracking notes to avoid going in circles.

## Root problem

Current TS Trotter does:

```text
H::OperatorTS -> resum(H) -> ordinary translated gates
O::OperatorTS -> representative(O) -> ordinary workspace
```

For `m` TS Hamiltonian representatives and `N` translations, second-order Strang creates roughly `2mN - 1` ordinary gates. In MFIM, `m=3`, so this is `6N - 1` gates. Profiling showed the main slowdown is the `O(N)` translated-gate sweep and per-gate truncation, not one-time `resum(H)` allocation.

## Attempt tree

```text
1. Lazy resum / fold each full step
   ├─ Idea
   │  Avoid materializing resum(H); generate translated gates lazily and fold O back to OperatorTS after each full step.
   ├─ Result
   │  Correct but not faster. Sometimes slower.
   └─ Why it failed
      It removed a small allocation/setup cost but kept the same O(N) translated-gate algorithmic work.

2. Tune truncate_every
   ├─ Idea
   │  Current expanded Trotter truncates after every gate. Try less frequent truncation.
   ├─ Result
   │  Moderate speedups for values like 4–16, but unstable if too infrequent.
   └─ Why it failed
      This is a band-aid. It does not fix the wrong granularity of O(N) translated gates.

3. Orbit-level Strang splitting
   ├─ Idea
   │  Split H into translation-orbit Hamiltonians H_a^orb and use Strang over orbit Liouvillians L_a(O)=i[H_a^orb,O]/hbar.
   ├─ Result
   │  The right theoretical direction. Reduces product-formula factors from O(Nm) to O(m).
   └─ Remaining problem
      Need an efficient orbit-flow primitive exp(t L_a) acting directly in OperatorTS.

4. Taylor-2 orbit flow
   ├─ Idea
   │  Approximate each orbit flow by O + t L_a O + t^2 L_a^2 O / 2.
   ├─ Result
   │  Fast-ish and showed orbit splitting can reduce work.
   └─ Why it failed
      It is too close to a commutator-based RK-style integrator, not a true Trotter/orbit exponential primitive. Removed from code.

5. Krylov orbit flow
   ├─ Idea
   │  Apply exp(t L_a) with matrix-free Arnoldi/Krylov using orbit_liouvillian.
   ├─ Result
   │  Conceptually correct and now the clean TrotterTS prototype.
   └─ Why it is slow
      Krylov does many OperatorTS operations: ~20 orbit_liouvillian calls, ~100 trace_product inner products, ~50 residual updates per step for m=3,k=4. The backend still uses generic TS commutator scanning all relative shifts.

6. orbit_edges graph interface
   ├─ Idea
   │  Expose Pauli-orbit graph edges q -> r for one orbit generator, to reason about sparse orbit graph structure.
   ├─ Result
   │  Useful for theory and measurement. Shows folded degree around ~4 in sampled MFIM states.
   └─ Why not sufficient
      Generic orbit_edges currently still scans all translations and allocates temporary OperatorTS outputs. It is not yet a fast backend.

7. 1D active-shift specializer
   ├─ Idea
   │  Avoid scanning all shifts by enumerating only support-overlap relative shifts in 1D periodic systems.
   ├─ Result
   │  Could improve hot-loop timings after making it non-allocating.
   └─ Why it was removed
      It was target/example-specific. The core PR should not rely on 1D/MFIM-specific structure.

8. Naive dynamic support-overlap shifts
   ├─ Idea
   │  Generic-ish support-overlap shift filtering per Pauli pair.
   ├─ Result
   │  Correct but much slower.
   └─ Why it failed
      Per-pair allocation of support lists/sets/vectors and bookkeeping dominated any reduction in shift count.

9. Local component exponentials
   ├─ Idea
   │  Build connected components of the orbit graph and apply dense exp on each component.
   ├─ Result
   │  Components looked small in samples, so theory is promising.
   └─ Why naive version failed
      Rebuilt components and dense matrices every orbit flow/time step. Too many allocations and dense exp calls. Removed from production code.

10. Edge caching for Krylov
    ├─ Idea
    │  Cache orbit_edges(Ha,q) by orbit and input Pauli orbit.
    ├─ Result
    │  High hit rate but slower overall.
    └─ Why it failed
       Edge generation was not the main bottleneck in that implementation; cache storage/lookups and temporary edge Operators added overhead.

11. Vector-space Krylov with Dict coefficients
    ├─ Idea
    │  Avoid OperatorTS temporaries by using Dict{PauliStringTS,ComplexF64} vectors.
    ├─ Result
    │  Lower or similar allocations, but much slower.
    └─ Why it failed
       Plain Dict coefficient vectors are too slow. A real vector backend would need indexed bases and sparse edge lists, not ad hoc dictionaries.

12. Linear-combination helper for OperatorTS
    ├─ Idea
    │  Reduce allocations in Krylov residual updates and final axpy by accumulating into one dictionary.
    ├─ Result
    │  Little/no improvement, sometimes slower.
    └─ Why it failed
       Existing OperatorTS primitives are already reasonable. The helper added complexity without clear hot-loop gain.

13. Adjoint-free trace_product helper
    ├─ Idea
    │  Avoid constructing V[i]' in Krylov inner products.
    ├─ Result
    │  Little/no allocation or time improvement.
    └─ Why it failed
       Adjoint construction was not a dominant allocation source. Reverted.

14. Cache Krylov adjoints and norms
    ├─ Idea
    │  Store Vdag[i] and trace_product(Vdag[i], V[i]) once per Krylov vector.
    ├─ Result
    │  Small allocation/time improvement in public step.
    └─ Status
       Kept because it is simple and avoids repeated intermediate adjoints/norms in a hot loop.
```

## Current understanding

The current code is clean and generic:

```text
TrotterTS = Strang over TS orbit Liouvillian exponentials
orbit_liouvillian = i * commutator(A, B) / hbar
orbit_edges = generic graph-inspection primitive
```

The remaining obstacle is that `commutator(OperatorTS, OperatorTS)` is TS-native in representation but not fully TS-native in cost: it scans every relative translation. Thus the orbit flow still pays an `O(N)` relative-shift factor per input orbit.

## Real next target

Develop a generic, dimension-independent active relative-shift generator for Pauli-orbit edges:

```text
current: orbit_liouvillian -> commutator -> all relative shifts
wanted:  orbit_liouvillian -> active orbit-edge generator -> only shifts that can contribute
```

Requirements:

- no MFIM-specific assumptions,
- no 1D-only hot path as the core method,
- no per-pair dynamic allocation like support vectors/sets,
- validate exactly against current commutator/resum,
- measure shift counts, allocations, and timings before trusting final runtime.
