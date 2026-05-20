"""
    AbstractEvolutionMethod

Supertype for time-evolution algorithm specifications passed to [`evolve`](@ref).
Each concrete subtype carries only that algorithm's own options. New methods are
new subtypes; the signature of [`evolve`](@ref) does not change.
"""
abstract type AbstractEvolutionMethod end

"""
    Trotter(; order=2, gates=nothing)

Fixed-step product-formula integrator. `order` is 1 (Lie) or 2 (Strang). `gates`
is an optional precomputed gate list (from [`trotterize`](@ref)) matching the
save-interval `dt`; if `nothing`, gates are built internally from `H` and cached
when `tspan` is uniformly spaced.

Only implemented for `H::Operator` and `O::Operator`. For translation-symmetric
operators, use [`evolve_trotter`](@ref) directly.
"""
struct Trotter <: AbstractEvolutionMethod
    order::Int
    gates::Any
end
Trotter(; order::Integer=2, gates=nothing) = Trotter(Int(order), gates)

"""
    RK4()

Classical fixed-step 4th-order Runge–Kutta. Takes one internal step per
save-interval; subdivide by supplying a finer `tspan`.
"""
struct RK4 <: AbstractEvolutionMethod end

"""
    DOPRI5()

Dormand–Prince 5(4) embedded Runge–Kutta. In this version it is used as a
fixed-step method: one internal step per save-interval.
"""
struct DOPRI5 <: AbstractEvolutionMethod end

"""
    Exact()

Reference integrator. Builds dense matrices `Hm = Matrix(H)` and `Om = Matrix(O)`
once, then evolves with `Om <- U * Om * U'` where `U = exp(i Hm dt / hbar)`,
exponentiating per save-interval (cached when `tspan` is uniform).

`history` elements and `final` are dense matrices, not `Operator`s. Intended for
testing/benchmarking small systems against the truncated-Pauli methods; cost and
memory are exponential in the qubit count.
"""
struct Exact <: AbstractEvolutionMethod end

"""
    EvolutionResult{T, H, O}

Result of [`evolve`](@ref).

# Fields
- `t::Vector{T}`: the save times (`collect(tspan)`).
- `history::H`: a `Vector` whose `i`-th element is `fout(O(t[i]), t[i])`, or
  `nothing` if `fout === nothing` was passed.
- `final::O`: the final operator, always populated regardless of `fout`.
"""
struct EvolutionResult{T,H,O}
    t::Vector{T}
    history::H
    final::O
end

# ---------- internal helpers ----------

function _validate_tspan(tspan::AbstractVector)
    isempty(tspan) && throw(ArgumentError("tspan must be non-empty"))
    for i in 1:(length(tspan) - 1)
        tspan[i + 1] > tspan[i] || throw(ArgumentError("tspan must be strictly increasing"))
    end
    return nothing
end

# Allocate history container based on fout. Returns history or nothing.
_alloc_history(::Nothing, O, n) = nothing
function _alloc_history(fout, O, n)
    sample = fout(O)
    history = Vector{typeof(sample)}(undef, n)
    @inbounds history[1] = sample
    return history
end

@inline _save!(::Nothing, ::Any, ::Any, ::Any) = nothing
@inline function _save!(history, fout, O, i)
    @inbounds history[i] = fout(O)
    return nothing
end

# ---------- public entry ----------

"""
    evolve(H::AbstractOperator, O::AbstractOperator, tspan::AbstractVector; kwargs...) -> EvolutionResult
    evolve(H::AbstractOperator, O::AbstractOperator, dt::Real, nsteps::Integer; kwargs...) -> EvolutionResult

Evolve the operator `O` under the Hamiltonian `H` in the Heisenberg picture and
return an [`EvolutionResult`](@ref). Values are saved at the times in `tspan`;
`tspan[1]` is saved before any integration.

The integrator takes one internal step per save-interval. To use a finer
internal step than the spacing at which results are saved, pass a finer `tspan`.

# Keyword arguments
- `method::AbstractEvolutionMethod = RK4()`. One of [`Trotter`](@ref),
  [`RK4`](@ref), [`DOPRI5`](@ref), [`Exact`](@ref).
- `truncation`: function `O -> O` applied after every internal step. Default
  `identity`.
- `dissipation`: function `(O, dt) -> O` applied after every internal step. The
  `dt` argument reflects that dissipation is part of the dynamics (a Lindblad
  semigroup factor) and depends on the step size. Default `(O, dt) -> O`.
- `fout`: function `fout(O)` called at every save time; its return value is
  collected into `history`. Default `(O) -> nothing`. Pass `fout = copy` to
  save the full operator trajectory. Pass `fout = nothing` to skip history
  collection entirely. `final` is populated independently of `fout`.
- `hbar::Real = 1`: Planck constant, forwarded to the underlying step routine.


# Examples

Full-trajectory (`fout = copy`):
```julia
result = evolve(H, O, 0.0:0.05:1.0; fout = copy)
result.history          # Vector of operators saved at each time
result.final            # final operator
```

Scalar observable trajectory:
```julia
result = evolve(H, O, 0.0:0.05:1.0; fout = O -> trace(O * A))
result.history          # Vector{ComplexF64} of ⟨A(t)⟩
result.final            # final operator, free for continuation
```

Skip history (default), keep only the final state:
```julia
result = evolve(H, O, 0.0:0.05:1.0; fout = nothing)
result.history === nothing
result.final            # final operator
```

Convenience `(dt, nsteps)` form:
```julia
evolve(H, O, 0.01, 100)        # 100 steps of dt = 0.01, 101 save points
```

Dissipative run:
```julia
evolve(H, O, 0.0:0.05:1.0; dissipation = (O, dt) -> add_noise(O, 0.1 * dt))
```

Choose an integrator:
```julia
evolve(H, O, 0.0:0.05:1.0; method = RK4())
evolve(H, O, 0.0:0.05:1.0; method = DOPRI5())
```
"""
function evolve(H::AbstractOperator, O::AbstractOperator, tspan::AbstractVector;
                method::AbstractEvolutionMethod = RK4(),
                truncation = identity,
                dissipation = (O, dt) -> O,
                fout = O -> nothing,
                hbar::Real = 1)
    _validate_tspan(tspan)
    return _evolve(method, H, O, tspan;
                   truncation=truncation,
                   dissipation=dissipation,
                   fout=fout,
                   hbar=hbar)
end

evolve(H::AbstractOperator, O::AbstractOperator, dt::Real, nsteps::Integer; kwargs...) =
    evolve(H, O, range(0, step=dt, length=nsteps + 1); kwargs...)




# ---------- method-specific implementations ----------



function _evolve(::RK4, H::AbstractOperator, O::AbstractOperator, tspan;
                 truncation, dissipation, fout, hbar)
    n = length(tspan)
    history = _alloc_history(fout, O, n)
    for i in ProgressBar(1:(n - 1))
        dt = tspan[i + 1] - tspan[i]
        O = rk4(H, O, dt; hbar=hbar, heisenberg=true, truncation=truncation)
        O = dissipation(O, dt)
        O = truncation(O)
        _save!(history, fout, O, i + 1)
    end
    return EvolutionResult(collect(tspan), history, O)
end

function _evolve(::DOPRI5, H::AbstractOperator, O::AbstractOperator, tspan;
                 truncation, dissipation, fout, hbar)
    n = length(tspan)
    history = _alloc_history(fout, O, n)
    for i in ProgressBar(1:(n - 1))
        dt = tspan[i + 1] - tspan[i]
        O = dopri5(H, O, dt; hbar=hbar, heisenberg=true, truncation=truncation)
        O = dissipation(O, dt)
        O = truncation(O)
        _save!(history, fout, O, i + 1)
    end
    return EvolutionResult(collect(tspan), history, O)
end


function _evolve(method::Trotter, H::Operator, O::Operator, tspan;
                 truncation, dissipation, fout, hbar)
    n = length(tspan)
    history = _alloc_history(fout, O, n)
    # `trotter_step!` mutates its operator in place; copy to avoid aliasing the
    # caller's `O`, which would corrupt observables like `fout(O) = trace_product(O0, O)`.
    O = copy(O)

    # Cache gates when the save spacing is uniform; rebuild per step otherwise.
    dt0 = n > 1 ? (tspan[2] - tspan[1]) : zero(eltype(tspan))
    uniform = n > 1 && all(i -> tspan[i + 1] - tspan[i] ≈ dt0, 1:(n - 1))
    gates_cached = method.gates !== nothing ? method.gates :
                   (uniform && n > 1 ?
                    trotterize(H, dt0; order=method.order, heisenberg=true, hbar=hbar) :
                    nothing)

    for i in ProgressBar(1:(n - 1))
        dt = tspan[i + 1] - tspan[i]
        g = gates_cached !== nothing ? gates_cached :
            trotterize(H, dt; order=method.order, heisenberg=true, hbar=hbar)
        trotter_step!(O, g; truncation=truncation)
        O = dissipation(O, dt)
        O = truncation(O)
        _save!(history, fout, O, i + 1)
    end
    return EvolutionResult(collect(tspan), history, O)
end

function _evolve(::Trotter, H::AbstractOperator, O::AbstractOperator, tspan;
                 truncation, dissipation, fout, hbar)
    throw(ArgumentError("Trotter evolution via `evolve` is implemented for `Operator` only, " *
                        "not for `$(typeof(H))`. Use `evolve_trotter` for translation-symmetric operators."))
end


function _evolve(::Exact, H::AbstractOperator, O::AbstractOperator, tspan;
                 truncation, dissipation, fout, hbar)
    n = length(tspan)
    Hm = Matrix(H)
    Om = Matrix(O)
    history = _alloc_history(fout, Om, n)

    # Cache the propagator when the save spacing is uniform.
    dt0 = n > 1 ? (tspan[2] - tspan[1]) : zero(eltype(tspan))
    uniform = n > 1 && all(i -> tspan[i + 1] - tspan[i] ≈ dt0, 1:(n - 1))
    U_cached = uniform && n > 1 ? exp((1im * dt0 / hbar) * Hm) : nothing

    for i in ProgressBar(1:(n - 1))
        dt = tspan[i + 1] - tspan[i]
        U = U_cached !== nothing ? U_cached : exp((1im * dt / hbar) * Hm)
        Om = U * Om * adjoint(U)
        Om = dissipation(Om, dt)
        Om = truncation(Om)
        _save!(history, fout, Om, i + 1)
    end
    return EvolutionResult(collect(tspan), history, Om)
end
