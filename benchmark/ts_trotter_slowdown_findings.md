# TS Trotter slowdown notes

Date: 2026-06-04

## Problem

The current translation-symmetric (TS) Trotter path does not beat RK4 on the TS mixed-field Ising example, unlike the non-TS case.

Initial suspicion was that the slowdown came mainly from expanding `H` with `resum(H)` and evolving outside the TS representation.

## Current code path

Baseline TS Trotter in `src/evolution.jl` does:

```julia
Or = representative(O)
Hr = resum(H)
gates = trotterize(Hr, dt)
trotter_step!(Or, gates; truncation=truncation)
```

For MFIM:

```text
length(Hts) = 3
length(resum(Hts)) = 3N
second-order Strang gates = 2 * length(resum(Hts)) - 1 = 6N - 1
```

So at `N=64`, a 3-term TS Hamiltonian becomes 383 ordinary single-Pauli rotation gates per Trotter step.

## Prototype tried

Implemented a lazy TS gate-list prototype that:

- avoided materializing `resum(H)`,
- generated translated gates on demand,
- folded the evolved operator back into `OperatorTS` after every complete Trotter step,
- applied truncation in TS space.

This prototype was reverted from the working tree and saved as:

```text
stash@{0}: lazy-ts-trotter-prototype-before-ts-native
```

Recover with:

```bash
git stash apply stash@{0}
```

## Lazy prototype benchmark

Same MFIM benchmark, `dt=0.1`, best of 3 after warmup.

| Case | Version | RK4 best | Trotter best | Trotter/RK4 |
|---|---:|---:|---:|---:|
| `N=32, M=12, steps=50` | baseline | `1.173s` | `1.342s` | `1.14x` |
| `N=32, M=12, steps=50` | lazy/fold prototype | `1.169s` | `1.690s` | `1.44x` |
| `N=64, M=14, steps=20` | baseline | `2.285s` | `4.078s` | `1.78x` |
| `N=64, M=14, steps=20` | lazy/fold prototype | `2.231s` | `4.254s` | `1.91x` |

Conclusion: avoiding materialized `resum(H)` is not the main fix. The lazy/fold prototype did not improve the benchmark.

## Component profiling

Short-run component profile:

### `N=64, M=14, steps=8`

| Method | Total | Dominant measured costs |
|---|---:|---|
| RK4 | `0.539s` | `0.453s` in 32 TS commutators, `0.011s` in 40 trims |
| baseline Trotter | `1.503s` | `0.368s` in 3072 trims, 383 gates/step |
| lazy/fold Trotter | `1.909s` | `1.892s` step body, `0.587s` in 3072 trims |

One-time costs were small:

```text
resum(H)    ≈ 0.05s in the longer component run
build gates ≈ 0.11s in the longer component run
```

and much smaller in one-step warmed checks.

## Definitive checks

### 1. Size scaling at fixed TS complexity

Setup:

- MFIM, so `length(Hts) = 3` for every `N`.
- Grow `Ots` to about 2048 TS representative strings.
- Compare one RK4 step to one cached Trotter step using `trotterize(resum(Hts), dt)`.

| N | `length(Hts)` | `length(resum(Hts))` | Strang gates | `length(Ots)` | RK4 step | Trotter step | Ratio |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 16 | 3 | 48 | 95 | 2048 | `0.0105s` | `0.0172s` | 1.64x |
| 24 | 3 | 72 | 143 | 2048 | `0.0128s` | `0.0192s` | 1.51x |
| 32 | 3 | 96 | 191 | 2048 | `0.0147s` | `0.0319s` | 2.17x |
| 48 | 3 | 144 | 287 | 2048 | `0.0193s` | `0.0324s` | 1.68x |
| 64 | 3 | 192 | 383 | 2048 | `0.0226s` | `0.0470s` | 2.08x |

At `N=64`, one-time costs were tiny compared to the step:

```text
resum(H)      ≈ 0.00027s
gate building ≈ 0.000026s
Trotter step  ≈ 0.047s
```

### 2. Prefix gate-list scaling

Setup:

- `N=64`
- fixed `length(Ots)=4096`
- full Strang gate list length = 383
- apply only first `K` ordinary Trotter gates

| Prefix K | Time | Time/gate | Final length |
|---:|---:|---:|---:|
| 8 | `0.00103s` | `0.000129s` | 4704 |
| 16 | `0.00189s` | `0.000118s` | 5039 |
| 32 | `0.00367s` | `0.000115s` | 5724 |
| 64 | `0.00887s` | `0.000139s` | 7548 |
| 128 | `0.0430s` | `0.000336s` | 8192 |
| 256 | `0.0766s` | `0.000299s` | 8192 |
| 383 | `0.107s` | `0.000281s` | 8192 |

This directly shows cost grows with the number of translated ordinary gates. After `K≈128`, truncation and cap effects also contribute.

### 3. `truncate_every` check

Setup:

- `N=64`
- `length(Ots)=4096`
- 383 gates
- trim cap 8192

| `truncate_every` | Time | Final length | Max intermediate length |
|---:|---:|---:|---:|
| 1 | `0.108s` | 8192 | 12169 |
| 2 | `0.102s` | 8192 | 13135 |
| 4 | `0.0958s` | 8192 | 13125 |
| 8 | `0.0926s` | 8790 | 19901 |
| 16 | `0.0918s` | 9636 | 21798 |
| 32 | `0.103s` | 9954 | 32107 |
| 64 | `0.191s` | 14497 | 135942 |

Tuning truncation frequency gives only a modest improvement. It is not the structural fix.

## Finding

The slowdown is not primarily caused by storing `resum(H)` or by final contraction/folding.

The slowdown is caused by the computational primitive:

```text
current TS Trotter = O(N) ordinary translated single-Pauli gate sweeps per step
TS RK4             = a few TS commutator-kernel passes per step
```

For MFIM at `N=64`:

```text
TS Hamiltonian representatives: 3
current second-order gates:     383
```

The prefix-scaling check verifies that runtime tracks the number of translated gates. Therefore, a genuine TS-native Trotter primitive should operate at the level of TS representatives/orbits rather than individual translated ordinary gates.

## Implication

Lazy `resum` generation and `truncate_every` tuning are fine-grained optimizations only.

The likely real fix is a TS-native Trotter/orbit primitive that replaces:

```text
6N - 1 ordinary gate sweeps
```

with a small number of TS-native orbit-level operations, while preserving a justified second-order product formula.
