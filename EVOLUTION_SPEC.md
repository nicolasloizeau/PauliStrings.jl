# PauliStrings.jl — Unified Time Evolution API

## Goal

Create a new file `src/evolution.jl` that defines a unified, extensible
time-evolution API for PauliStrings.jl. This API is the long-term public
interface — the function signature and return type **must not change** in
backward-incompatible ways once shipped, so the design below should be
followed precisely.

The new `evolve` function is a thin orchestrator. It does **not**
reimplement integration logic. The actual per-step work is delegated to
existing single-step routines in `src/runge_kutta.jl` and
`src/trotter.jl`. Discover the relevant single-step functions in those
files and call them from the method-specific `_evolve` implementations.
Do not duplicate their logic.

## Public API contract

Any future change to the items in this section is a breaking change and
must be avoided.

### Function signature

```julia
evolve(H, O, tspan::AbstractVector; kwargs...) -> EvolutionResult
```

Convenience method for fixed-step usage:

```julia
evolve(H, O, dt::Real, nsteps::Integer; kwargs...) -> EvolutionResult
# lowers to evolve(H, O, range(0, step=dt, length=nsteps+1); kwargs...)
```

### Positional arguments

- `H` — the Hamiltonian (or list of terms / gates, depending on method).
  Type is intentionally loose; method-specific `_evolve` implementations
  dispatch on what they need.
- `O` — the initial operator (Heisenberg picture). Type loose for the
  same reason.
- `tspan::AbstractVector` — vector of save times. The integrator may
  take any internal step pattern between save points; values are saved
  at exactly these times. `tspan[1]` is included in the output (the
  initial state is saved before any evolution).

### Keyword arguments

- `method` — algorithm-configuration struct (see below). Default:
  `Trotter()`.
- `truncation::TruncationScheme` — operator-truncation options. Default:
  `TruncationScheme()`.
- `dissipation::Union{Nothing, AbstractDissipation}` — open-system
  dynamics. Default: `nothing` (unitary evolution).
- `fout` — function called as `fout(O, t)` at each save time. Its return
  value is collected into `history`. Default: `copy` (saves the full
  operator trajectory). Set `fout = nothing` to skip history collection
  entirely.
- `hbar::Real` — Planck constant. Default: `1`.

### Return type

```julia
struct EvolutionResult{T, H, O}
    t::Vector{T}      # the save times (== collect(tspan))
    history::H        # Vector of fout returns, or nothing if fout === nothing
    final::O          # the final operator, always returned
end
```

`final` is always populated, regardless of `fout`. This is the key
design choice: `fout` controls *history*, `final` is independent. Users
who want only the final state set `fout = nothing` and read
`result.final`. Users who want a scalar timeseries set
`fout = O -> tr(O*A)` and still get `result.final` for free /
continuation.

The return is a `NamedTuple` constructed via `(; t, history, final)`
**or** a struct `EvolutionResult` — pick the struct version, since it
gives room to add fields (`diagnostics`, `method_state`, ...) without
breaking destructuring later. Define it in `evolution.jl`.

## Required structs

All of the following live in `src/evolution.jl`. Export the user-facing
ones from the package's main module.

### Method types

```julia
abstract type AbstractEvolutionMethod end

struct Trotter <: AbstractEvolutionMethod
    order::Int
    gates::Any        # precomputed gates, or nothing -> built internally
end
Trotter(; order=2, gates=nothing) = Trotter(order, gates)

struct RungeKutta <: AbstractEvolutionMethod
    order::Int        # e.g. 4 for classical RK4; pick what runge_kutta.jl supports
end
RungeKutta(; order=4) = RungeKutta(order)
```

If `runge_kutta.jl` exposes adaptive variants (DOPRI5, etc.), add a
corresponding struct (e.g. `RKAdaptive` or `DP5`) with the relevant
tolerance fields. Do not invent options that the underlying single-step
routines don't support; mirror what's already there.

If there are other existing single-step methods in the codebase
(Krylov, Magnus, etc.), add one struct per method following the same
pattern. **One struct per algorithm.** Each carries only that
algorithm's options.

### Truncation

```julia
struct TruncationScheme
    M::Int            # max number of Pauli strings to keep
    k::Int            # weight cutoff (0 = off)
    trim_every::Int   # apply truncation every N internal steps
end
TruncationScheme(; M=2^20, k=0, trim_every=1) =
    TruncationScheme(M, k, trim_every)
```

Names of fields must match what existing truncation routines in
PauliStrings.jl already accept; if the existing convention uses
different names (e.g. `keep`, `weight_cutoff`), follow the existing
convention rather than inventing new ones. The point is to bundle the
existing scattered kwargs (`M`, `k_truncate`, `trim_every`) into one
struct.

### Dissipation

```julia
abstract type AbstractDissipation end

struct DepolarizingNoise <: AbstractDissipation
    rate::Float64
end

# Add other concrete types as needed by existing noise routines in the
# codebase (dephasing, weight-dependent, generic Lindbladian, ...).
```

`dissipation = nothing` means unitary. `dissipation::AbstractDissipation`
selects an open-system implementation. Mirror whatever noise routines
already exist in the package; do not invent dissipation models that
aren't already implemented somewhere.

## Dispatch structure

```julia
# Public entry point — does input validation and orchestration.
function evolve(H, O, tspan::AbstractVector;
                method::AbstractEvolutionMethod = Trotter(),
                truncation::TruncationScheme = TruncationScheme(),
                dissipation::Union{Nothing, AbstractDissipation} = nothing,
                fout = copy,
                hbar::Real = 1)
    # 1. Validate tspan (non-empty, sorted).
    # 2. Allocate history container based on fout (nothing -> no allocation;
    #    otherwise call fout once on (O, tspan[1]) to determine element type,
    #    then preallocate Vector{T}(undef, length(tspan))).
    # 3. Save initial entry: history[1] = fout(O, tspan[1]) if applicable.
    # 4. Dispatch to _evolve.
    return _evolve(method, H, O, tspan;
                   truncation=truncation,
                   dissipation=dissipation,
                   fout=fout,
                   hbar=hbar)
end

# Convenience.
evolve(H, O, dt::Real, nsteps::Integer; kwargs...) =
    evolve(H, O, range(0, step=dt, length=nsteps+1); kwargs...)

# Method-specific implementations — one per concrete method type.
function _evolve(method::Trotter, H, O, tspan; truncation, dissipation, fout, hbar)
    # Loop over save intervals tspan[i] -> tspan[i+1].
    # Within each interval, repeatedly call the existing single-step Trotter
    # routine from src/trotter.jl with the appropriate dt.
    # Apply truncation every truncation.trim_every internal steps.
    # If dissipation !== nothing, dispatch on its type to apply noise.
    # After reaching tspan[i+1], call fout(O, tspan[i+1]) and store.
    # Return EvolutionResult(t=collect(tspan), history=..., final=O).
end

function _evolve(method::RungeKutta, H, O, tspan; truncation, dissipation, fout, hbar)
    # Same shape, but call the single-step RK routine from src/runge_kutta.jl.
end
```

The orchestration code (history allocation, save-time loop,
truncation cadence, fout handling) should be **shared** between
`_evolve` methods where possible. Factor it into helper functions if
duplication starts to hurt. The only thing that should genuinely
differ between methods is the per-step call.

## Implementation notes

- **Do not reimplement Trotter or RK steps.** Locate the existing
  single-step functions in `src/trotter.jl` and `src/runge_kutta.jl`
  and call them. If their signatures need a small adapter, write the
  adapter; don't rewrite the integrator.
- **`fout` element type.** Call `fout(O, tspan[1])` once before the
  loop to get a sample value, use `typeof(sample)` to allocate
  `Vector{T}(undef, length(tspan))`, and write by index. Avoids
  `push!` overhead and gives concrete element types.
- **`fout = nothing` path.** Skip history allocation entirely; set the
  `history` field of `EvolutionResult` to `nothing`.
- **`fout = copy` default.** If any internal evolution step mutates `O`
  in place, `fout = identity` would store references to the same
  evolving object. `copy` is the safe default.
- **Validate `tspan`.** It must be non-empty and monotonically
  increasing. Throw `ArgumentError` otherwise.
- **`hbar`.** Pass it through to the single-step routines that need it.
  Don't apply it in two places.
- **Truncation cadence.** `trim_every` counts internal steps, not save
  intervals. Track an internal step counter across save intervals.
- **Dissipation integration point.** Decide whether dissipation is
  applied as a separate sub-step (operator splitting) or folded into
  the main step (depending on the method). Mirror whatever pattern
  existing noise routines in the package already use; don't invent a
  new integration scheme.

## What to put where

```
src/
  evolution.jl         <- new file: structs, evolve, _evolve methods
  trotter.jl           <- existing: single-step Trotter, leave alone
  runge_kutta.jl       <- existing: single-step RK, leave alone
  PauliStrings.jl      <- main module: include("evolution.jl") and
                          export evolve, EvolutionResult, Trotter,
                          RungeKutta, TruncationScheme, DepolarizingNoise,
                          (and any other method/dissipation types added)
```

## Backward compatibility

Existing functions like `evolve_trotter` should remain in place
unchanged for now. The new `evolve` lives alongside them. A future PR
can deprecate the old functions in favor of the new `evolve`; that's
out of scope for this change.

## Tests

Add a new test file `test/test_evolution.jl` covering, at minimum:

- Trotter and RungeKutta produce equivalent results (within tolerance)
  on a small, fixed problem.
- `fout = nothing` returns `history = nothing` and a correct `final`.
- `fout = O -> tr(O*A)` for some observable `A` returns a correct
  scalar trajectory and a correct `final`.
- `fout = copy` (default) returns a `history` whose last element equals
  `final` (up to copy).
- The convenience signature `evolve(H, O, dt, nsteps; ...)` produces
  the same result as the explicit `tspan` form.
- `tspan[1]` is included in the output (initial state is saved).
- Validation: empty `tspan`, non-monotonic `tspan`, both raise
  `ArgumentError`.
- A short run with `dissipation = DepolarizingNoise(rate=...)` produces
  results that differ from the unitary run in the expected direction
  (e.g. operator norm decay).

## Style

- Match the existing code style of PauliStrings.jl (indentation,
  naming, docstring conventions). Read a couple of existing files
  first.
- Add docstrings to every exported symbol. The `evolve` docstring
  should include a small worked example of each of: full-trajectory
  default, scalar-observable `fout`, `fout = nothing`, and the
  convenience `(dt, nsteps)` form.
- Keep `evolution.jl` self-contained except for calls into the
  existing single-step routines and existing truncation/noise
  utilities.

## Out of scope

- Adaptive step-size control on top of fixed-step methods. If the
  underlying single-step routine supports adaptivity, expose it via
  the corresponding method struct's fields; otherwise leave it alone.
- A `callback` hook in addition to `fout`. `fout` is the only hook.
  Side effects (printing, logging) go inside the user's `fout`.
- Early termination from `fout`. Not supported in this version.
- Deprecation of `evolve_trotter` and friends. Separate PR.
