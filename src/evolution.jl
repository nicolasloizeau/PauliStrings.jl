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

Implemented for `H::Operator`/`O::Operator` and for `H::OperatorTS`/`O::OperatorTS`. In
the translation-symmetric case the gates are built from `resum(H)` and applied to the
representative term, so a precomputed `gates` list must come from `trotterize(resum(H), dt)`.
"""
struct Trotter <: AbstractEvolutionMethod
    order::Int
    gates::Any
end
Trotter(; order::Integer=2, gates=nothing) = Trotter(Int(order), gates)

"""
    TrotterTS(; order=2, componenttol=0.9999)

Translation-symmetric orbit-level product formula. For `H::OperatorTS`, splits the
Hamiltonian into its representative translation orbits and applies a first-order
(`order=1`) or second-order (`order=2`, Strang) product over those TS orbit
Liouvillians. Each orbit flow is applied by exponentiating connected components of
the Pauli-orbit graph. Components are processed in descending coefficient weight
until `componenttol` cumulative weight is reached; lower-weight components are
frozen for that orbit flow.
"""
struct TrotterTS <: AbstractEvolutionMethod
    order::Int
    componenttol::Float64
end
TrotterTS(; order::Integer=2, componenttol::Real=0.9999) =
    TrotterTS(Int(order), Float64(componenttol))

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
    for i in 1:(length(tspan)-1)
        tspan[i+1] > tspan[i] || throw(ArgumentError("tspan must be strictly increasing"))
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
  [`TrotterTS`](@ref), [`RK4`](@ref), [`DOPRI5`](@ref), [`Exact`](@ref).
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
    method::AbstractEvolutionMethod=RK4(),
    truncation=identity,
    dissipation=(O, dt) -> O,
    fout=O -> nothing,
    hbar::Real=1)
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
    for i in ProgressBar(1:(n-1))
        dt = tspan[i+1] - tspan[i]
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
    for i in ProgressBar(1:(n-1))
        dt = tspan[i+1] - tspan[i]
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
    uniform = n > 1 && all(i -> tspan[i+1] - tspan[i] ≈ dt0, 1:(n-1))
    gates_cached = method.gates !== nothing ? method.gates :
                   (uniform && n > 1 ?
                    trotterize(H, dt0; order=method.order, heisenberg=true, hbar=hbar) :
                    nothing)

    for i in ProgressBar(1:(n-1))
        dt = tspan[i+1] - tspan[i]
        g = gates_cached !== nothing ? gates_cached :
            trotterize(H, dt; order=method.order, heisenberg=true, hbar=hbar)
        trotter_step!(O, g; truncation=truncation)
        O = dissipation(O, dt)
        O = truncation(O)
        _save!(history, fout, O, i + 1)
    end
    return EvolutionResult(collect(tspan), history, O)
end

function _evolve(method::Trotter, H::Operator{<:PauliStringTS}, O::Operator{<:PauliStringTS}, tspan;
    truncation, dissipation, fout, hbar)
    qubitsize(H) == qubitsize(O) && periodicflags(H) == periodicflags(O) ||
        throw(DimensionMismatch("H and O must share the same translation-symmetry lattice"))
    Ls = qubitsize(O)
    Ps = periodicflags(O)
    n = length(tspan)
    history = _alloc_history(fout, O, n)
    # Evolve the representative as a dense `Operator`: the gates are built from the
    # resummed (translation-unfolded) `H`, so applying them to the representative term
    # reproduces the translation-symmetric dynamics. `representative` returns a fresh
    # `Operator`, so `trotter_step!`'s in-place mutation never aliases the caller's `O`.
    Or = representative(O)
    Hr = resum(H)

    # Cache gates when the save spacing is uniform; rebuild per step otherwise.
    dt0 = n > 1 ? (tspan[2] - tspan[1]) : zero(eltype(tspan))
    uniform = n > 1 && all(i -> tspan[i+1] - tspan[i] ≈ dt0, 1:(n-1))
    gates_cached = method.gates !== nothing ? method.gates :
                   (uniform && n > 1 ?
                    trotterize(Hr, dt0; order=method.order, heisenberg=true, hbar=hbar) :
                    nothing)

    for i in ProgressBar(1:(n-1))
        dt = tspan[i+1] - tspan[i]
        g = gates_cached !== nothing ? gates_cached :
            trotterize(Hr, dt; order=method.order, heisenberg=true, hbar=hbar)
        trotter_step!(Or, g; truncation=truncation)
        Or = dissipation(Or, dt)
        Or = truncation(Or)
        _save!(history, fout, OperatorTS{Ls,Ps}(Or), i + 1)
    end
    return EvolutionResult(collect(tspan), history, OperatorTS{Ls,Ps}(Or))
end

function _gather_components(Ha::Operator{<:PauliStringTS}, O::Operator{<:PauliStringTS}, cache::_OrbitFlowCache, hbar::Real, maxlength::Int)
    coeff_lookup = Dict{eltype(O.strings),eltype(O.coeffs)}()
    sizehint!(coeff_lookup, length(O))
    for (q, c) in zip(O.strings, O.coeffs)
        coeff_lookup[q] = get(coeff_lookup, q, zero(eltype(O.coeffs))) + c
    end

    assigned = Set{eltype(O.strings)}()
    sizehint!(assigned, length(O))

    component_data = Tuple{Float64,Any,Vector{ComplexF64}}[]
    sizehint!(component_data, length(O))
    total_weight = 0.0

    for q0 in O.strings
        q0 in assigned && continue

        # Check cache first
        plan = get(cache.plans, q0, nothing)
        comp_seq, item = if plan !== nothing
            plan.component, plan
        else
            component, index, transitions = _orbit_component_and_transitions(Ha, q0, hbar, maxlength)
            component, (component, index, transitions)
        end

        coeffs = zeros(ComplexF64, length(comp_seq))
        weight = 0.0
        has_active = false
        for (j, q) in enumerate(comp_seq)
            c = get(coeff_lookup, q, nothing)
            if c !== nothing
                coeffs[j] += c
                weight += abs2(c)
                push!(assigned, q)
                has_active = true
            end
        end
        if has_active
            total_weight += weight
            push!(component_data, (weight, item, coeffs))
        end
    end
    return component_data, total_weight
end
function _propagate_batched_groups!(out_d, kept::Dict, dt::Real)
    for group in values(kept)
        plan0 = group[1][1]
        E = _component_exp(plan0, dt)
        coeffmat = Matrix{ComplexF64}(undef, length(plan0.component), length(group))
        for (j, (_, coeffs)) in enumerate(group)
            @inbounds coeffmat[:, j] = coeffs
        end
        coeffmat2 = E * coeffmat  # Real Float64 Matrix * ComplexF64 Matrix
        for (j, (plan, _)) in enumerate(group)
            for (i, q) in enumerate(plan.component)
                val = coeffmat2[i, j]
                if !iszero(val)
                    setwith!(+, out_d, q, val)
                end
            end
        end
    end
end
function _orbit_flow(Ha::Operator{<:PauliStringTS}, O::Operator{<:PauliStringTS}, dt::Real, cache::_OrbitFlowCache, method::TrotterTS, hbar::Real, truncation)
    componenttol = method.componenttol
    0 < componenttol <= 1 || throw(ArgumentError("componenttol must be in (0, 1]"))
    maxlength = 1000

    component_data, total_weight = _gather_components(Ha, O, cache, hbar, maxlength)

    sort!(component_data; by=x -> x[1], rev=true)
    keep_weight = componenttol * total_weight
    accumulated = 0.0
    out_d = emptydict(O)
    #sizehint!(out_d, length(O))

    PType = eltype(O.strings)
    kept = Dict{Any,Vector{Tuple{_OrbitComponentPlan{PType},Vector{ComplexF64}}}}()
    sizehint!(kept, length(component_data) ÷ 4)

    # Tail Freezing
    for (weight, item, coeffs) in component_data
        if accumulated < keep_weight
            plan = if item isa _OrbitComponentPlan
                item
            else
                component, index, transitions = item
                _component_plan!(cache, component, index, transitions)
            end

            sig = plan.signature
            vec = get!(kept, sig) do
                Tuple{_OrbitComponentPlan{PType},Vector{ComplexF64}}[]
            end
            push!(vec, (plan, coeffs))
            accumulated += weight
        else
            # Freeze low-weight components: carry their active coefficients forward unchanged.
            component_nodes = if item isa _OrbitComponentPlan
                item.component
            else
                item[1]
            end
            for (j, q) in enumerate(component_nodes)
                if !iszero(coeffs[j])
                    setwith!(+, out_d, q, coeffs[j])
                end
            end
        end
    end

    _propagate_batched_groups!(out_d, kept, dt)

    out = typeof(O)(collect(keys(out_d)), collect(values(out_d)))
    return truncation(out)
end

function _trotterts_step(Hterms, caches, O::Operator{<:PauliStringTS}, dt::Real, method::TrotterTS, hbar::Real, truncation)
    method.order ∈ (1, 2) || throw(ArgumentError("order must be 1 or 2, got $(method.order)"))
    isempty(Hterms) && return O
    if method.order == 1 || length(Hterms) == 1
        for j in eachindex(Hterms)
            O = _orbit_flow(Hterms[j], O, dt, caches[j], method, hbar, truncation)
        end
    else
        L = length(Hterms)
        for j in 1:(L-1)
            O = _orbit_flow(Hterms[j], O, dt / 2, caches[j], method, hbar, truncation)
        end
        O = _orbit_flow(Hterms[L], O, dt, caches[L], method, hbar, truncation)
        for j in (L-1):-1:1
            O = _orbit_flow(Hterms[j], O, dt / 2, caches[j], method, hbar, truncation)
        end
    end
    return O
end

function _evolve(method::TrotterTS, H::Operator{<:PauliStringTS}, O::Operator{<:PauliStringTS}, tspan;
    truncation, dissipation, fout, hbar)
    qubitsize(H) == qubitsize(O) && periodicflags(H) == periodicflags(O) ||
        throw(DimensionMismatch("H and O must share the same translation-symmetry lattice"))
    n = length(tspan)
    history = _alloc_history(fout, O, n)
    Hterms = _orbit_terms(H)
    caches = [_OrbitFlowCache(paulistringtype(H), typeof(H)) for _ in Hterms]
    O = copy(O)
    for i in ProgressBar(1:(n-1))
        dt = tspan[i+1] - tspan[i]
        O = _trotterts_step(Hterms, caches, O, dt, method, hbar, truncation)
        O = dissipation(O, dt)
        O = truncation(O)
        _save!(history, fout, O, i + 1)
    end
    return EvolutionResult(collect(tspan), history, O)
end

function _evolve(::TrotterTS, H::AbstractOperator, O::AbstractOperator, tspan;
    truncation, dissipation, fout, hbar)
    throw(ArgumentError("TrotterTS evolution via `evolve` is implemented for `OperatorTS` only, not for `$(typeof(H))`."))
end

function _evolve(::Trotter, H::AbstractOperator, O::AbstractOperator, tspan;
    truncation, dissipation, fout, hbar)
    throw(ArgumentError("Trotter evolution via `evolve` is implemented for `Operator` and " *
                        "`OperatorTS` only, not for `$(typeof(H))`."))
end


function _evolve(::Exact, H::AbstractOperator, O::AbstractOperator, tspan;
    truncation, dissipation, fout, hbar)
    n = length(tspan)
    Hm = Matrix(H)
    Om = Matrix(O)
    history = _alloc_history(fout, Om, n)

    # Cache the propagator when the save spacing is uniform.
    dt0 = n > 1 ? (tspan[2] - tspan[1]) : zero(eltype(tspan))
    uniform = n > 1 && all(i -> tspan[i+1] - tspan[i] ≈ dt0, 1:(n-1))
    U_cached = uniform && n > 1 ? exp((1im * dt0 / hbar) * Hm) : nothing

    for i in ProgressBar(1:(n-1))
        dt = tspan[i+1] - tspan[i]
        U = U_cached !== nothing ? U_cached : exp((1im * dt / hbar) * Hm)
        Om = U * Om * adjoint(U)
        Om = dissipation(Om, dt)
        Om = truncation(Om)
        _save!(history, fout, Om, i + 1)
    end
    return EvolutionResult(collect(tspan), history, Om)
end
