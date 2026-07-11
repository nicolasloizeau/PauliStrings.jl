# TODO: introduce to the main package
commutes(p, q) = iseven(count_ones((p.v & q.w) ⊻ (p.w & q.v)))
anticommutes(p, q) = !(commutes(p, q))

#= TODO: 1D optimization using bitrotate is insanely fast

@inline function _fast_1D_periodic_shift(v::T, shift_amount::Int, N::Int) where {T<:Unsigned}
    # Create a bitmask of exactly N ones: e.g., N=4 -> 0000...1111
    # (Note: if N == sizeof(T)*8, shifting by N overflows, so we handle it safely)
    mask = (N == sizeof(T) * 8) ? typemax(T) : (one(T) << N) - one(T)

    # Shift left by d, and wrap the overflow around to the right
    return ((v << shift_amount) | (v >> (N - shift_amount))) & mask
end
=#

@inline function _active_shift_tuple(cp::CartesianIndex{D}, cq::CartesianIndex{D}, Ls::NTuple{D,Int}, Ps::NTuple{D,Bool}) where {D}
    @inbounds for k in 1:D
        if !Ps[k] && cp[k] != cq[k]
            return false, Ls, 0
        end
    end

    shifts = ntuple(Val(D)) do k
        @inbounds begin
            if Ps[k]
                d_shift = mod(cp[k] - cq[k], Ls[k])
                iszero(d_shift) ? Ls[k] : d_shift
            else
                1
            end
        end
    end

    key = 0
    stride = 1
    @inbounds for k in 1:D
        if Ps[k]
            key += (shifts[k] - 1) * stride
        end
        stride *= Ls[k]
    end

    return true, shifts, key
end

function _ts_commutator!(d, prep::PauliString{N,T}, c::Number, qrep::PauliString{N,T}, Ls::NTuple{D,Int}, Ps::NTuple{D,Bool}; maxlength::Int=1000) where {N,T<:Unsigned,D}
    pauli_weight(prep) == 0 && return d
    pauli_weight(qrep) == 0 && return d

    R = CartesianIndices(Ls)

    seen_mask = zero(T)

    m_p = prep.v | prep.w
    m_q_init = qrep.v | qrep.w

    while m_p > 0
        tz_p = trailing_zeros(m_p)
        site_p = tz_p + 1
        cp = R[site_p]

        m_q = m_q_init
        while m_q > 0
            tz_q = trailing_zeros(m_q)
            site_q = tz_q + 1
            cq = R[site_q]

            valid, shifts, key = _active_shift_tuple(cp, cq, Ls, Ps)

            if valid
                bit = one(T) << key
                if iszero(seen_mask & bit)
                    seen_mask |= bit  # Mark as seen
                    q_shifted = shift(qrep, Ls, Ps, shifts)
                    if anticommutes(prep, q_shifted)
                        out, k = commutator(prep, q_shifted)
                        if count_ones(out.v | out.w) < maxlength
                            setwith!(+, d, PauliStringTS{Ls,Ps}(out), c * k)
                        end
                    end
                end
            end
            m_q = m_q ⊻ (one(T) << tz_q)
        end
        m_p = m_p ⊻ (one(T) << tz_p)
    end
    return d
end

"""
    TrotterTSGate{P,T,C}

One factor in a translation-symmetric Trotter step. Evolving a state operator under this
gate corresponds to an orbit-level flow under the generator `generator` (with weight `coeff`) 
for a duration of `dt`. Contains a mutable `cache` to store pre-computed Pauli-orbit plans 
and matrix exponentials.
"""
struct TrotterTSGate{P<:PauliStringTS,T<:Real,C}
    generator::P
    coeff::ComplexF64
    dt::T
    cache::C
end

struct OrbitComponentPlan{P,S}
    component::Vector{P}
    index::Dict{P,Int}
    matrix::Matrix{ComplexF64}
    exp_cache::Dict{Float64,Matrix{ComplexF64}}
    signature::S
end

mutable struct OrbitFlowCache{P,O,C,D,E}
    plans::Dict{P,OrbitComponentPlan{P}}
    fallback::Dict{P,O}
    assigned::Set{P}
    component_data::Vector{Tuple{Float64,Vector{P},Any,Vector{ComplexF64}}}
    out_d::C
    coeff_lookup::D
    edges_d::E
end

function OrbitFlowCache(::Type{P}, ::Type{O}) where {P,O}
    dummy = O(P[], ComplexF64[])
    out_d = emptydict(dummy)
    coeff_lookup = Dict{P,ComplexF64}()
    edges_d = emptydict(dummy)
    return OrbitFlowCache(
        Dict{P,OrbitComponentPlan{P}}(),
        Dict{P,O}(),
        Set{P}(),
        Tuple{Float64,Vector{P},Any,Vector{ComplexF64}}[],
        out_d,
        coeff_lookup,
        edges_d,
    )
end

function orbit_edges!(d, p::PauliStringTS{Ls,Ps,U}, c::Number, q::PauliStringTS{Ls,Ps,U}; maxlength) where {Ls,Ps,U}
    checklength(p, q)
    qrep = representative(q)
    pauli_weight(qrep) == 0 && return d
    prep = representative(p)
    _ts_commutator!(d, prep, c, qrep, Ls, Ps; maxlength)
    return d
end

function orbit_liouvillian(A::Operator{<:PauliStringTS}, B::Operator{<:PauliStringTS}; hbar::Real=1, epsilon::Real=0, maxlength::Int=1000)
    checklength(A, B)
    return 1im * commutator(A, B; epsilon=epsilon, maxlength=maxlength) / hbar
end


function orbit_component_and_transitions(gate::TrotterTSGate, seed::PauliStringTS, hbar::Real, maxlength::Int)
    component = typeof(seed)[seed]
    index = Dict{typeof(seed),Int}(seed => 1)
    transitions = Tuple{Int,Int,ComplexF64}[]
    edges_d = gate.cache.edges_d
    queue_index = 1

    while queue_index <= length(component)
        q = component[queue_index]
        j = queue_index
        queue_index += 1

        empty!(edges_d)
        orbit_edges!(edges_d, gate.generator, gate.coeff, q; maxlength)
        for (r, c) in pairs(edges_d)
            if !haskey(index, r)
                push!(component, r)
                index[r] = length(component)
            end
            i = index[r]
            push!(transitions, (i, j, 1im * c / hbar))
        end
    end
    return component, index, transitions
end

function get_component_plan!(cache::OrbitFlowCache, component, index, transitions)
    # Return already cached plan if another node in this component populated it first
    plan = get(cache.plans, component[1], nothing)
    plan !== nothing && return plan

    n = length(component)
    A = zeros(ComplexF64, n, n)
    for (i, j, val) in transitions # A is small
        A[i, j] += val
    end
    sig = component_signature(transitions, n)
    plan = OrbitComponentPlan(component, index, A, Dict{Float64,Matrix{ComplexF64}}(), sig)
    for q in component
        cache.plans[q] = plan
    end
    return plan
end

function component_exp(plan::OrbitComponentPlan, dt::Real)
    key = Float64(dt)
    E = get(plan.exp_cache, key, nothing)
    if E === nothing
        E = exp(key * plan.matrix)
        plan.exp_cache[key] = E
    end
    return E
end

function component_signature(transitions::Vector{Tuple{Int,Int,ComplexF64}}, n::Int; digits::Int=12)
    entries = Tuple{Int,Int,ComplexF64}[]
    sizehint!(entries, length(transitions))
    for (i, j, val) in transitions
        iszero(val) && continue
        push!(entries, (i, j, round(val, digits=digits)))
    end
    sort!(entries)
    return (n, entries)
end

function build_orbit_components(gate::TrotterTSGate, O::Operator{<:PauliStringTS}, hbar::Real, maxlength::Int)
    cache = gate.cache
    coeff_lookup = cache.coeff_lookup
    empty!(coeff_lookup)
    sizehint!(coeff_lookup, length(O))
    for (q, c) in zip(O.strings, O.coeffs)
        coeff_lookup[q] = c
    end

    component_data = cache.component_data
    empty!(component_data)
    sizehint!(component_data, length(O))

    total_weight = 0.0

    for q0 in O.strings
        haskey(coeff_lookup, q0) || continue

        # Check cache first
        plan = get(cache.plans, q0, nothing)
        if plan === nothing
            component, index, transitions = orbit_component_and_transitions(gate, q0, hbar, maxlength)
            plan = get_component_plan!(cache, component, index, transitions)
        end
        comp_seq = plan.component

        coeffs = zeros(ComplexF64, length(comp_seq))
        weight = 0.0
        has_active = false

        for (j, q) in enumerate(comp_seq)
            c = pop!(coeff_lookup, q, nothing)
            if c !== nothing
                coeffs[j] = c
                weight += abs2(c)
                push!(cache.assigned, q)
                has_active = true
            end
        end
        if has_active
            total_weight += weight
            push!(component_data, (weight, comp_seq, plan, coeffs))
        end
    end
    return component_data, total_weight
end
function propagate_batches!(out_d, kept::Dict, dt::Real)
    for group in values(kept)
        plan0 = group[1][1]
        E = component_exp(plan0, dt)
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
function orbit_update(gate::TrotterTSGate, O::Operator{<:PauliStringTS}, hbar::Real, truncation; componenttol, maxlength)
    0 < componenttol <= 1 || throw(ArgumentError("componenttol must be in (0, 1]"))
    0 < maxlength || throw(ArgumentError("maxlength must be > 0"))

    component_data, total_weight = build_orbit_components(gate, O, hbar, maxlength)

    sort!(component_data; by=x -> x[1], rev=true)
    keep_weight = componenttol * total_weight
    accumulated = 0.0
    out_d = gate.cache.out_d
    empty!(out_d)

    PType = eltype(O.strings)
    kept = Dict{Any,Vector{Tuple{OrbitComponentPlan{PType},Vector{ComplexF64}}}}()
    sizehint!(kept, length(component_data) ÷ 4)

    for (weight, comp_seq, item, coeffs) in component_data
        if accumulated < keep_weight
            plan = if item isa OrbitComponentPlan
                item
            else
                component, index, transitions = item
                get_component_plan!(cache, component, index, transitions)
            end

            sig = plan.signature
            vec = get!(kept, sig) do
                Tuple{OrbitComponentPlan{PType},Vector{ComplexF64}}[]
            end
            push!(vec, (plan, coeffs))
            accumulated += weight
        else
            for (j, q) in enumerate(comp_seq)
                if !iszero(coeffs[j])
                    setwith!(+, out_d, q, coeffs[j])
                end
            end
        end
    end

    propagate_batches!(out_d, kept, gate.dt)

    out = typeof(O)(collect(keys(out_d)), collect(values(out_d)))
    return truncation(out)
end

"""
    ts_trotterize(H::Operator{<:PauliStringTS}, dt::Real; order=2, heisenberg=true, hbar=1, caches=nothing)

Build a first-order (`order=1`) or second-order (`order=2`) translation-symmetric Trotter 
sequence of `TrotterTSGate` factors approximating the evolution. 

To reuse plans and matrix exponentials across multiple time steps, pass a vector of 
pre-allocated caches to the `caches` keyword argument.
"""
function ts_trotterize(H::Operator{<:PauliStringTS}, dt::Real;
    order::Integer=2, heisenberg::Bool=true, hbar::Real=1, caches=nothing)
    order ∈ (1, 2) || throw(ArgumentError("order must be 1 or 2, got $order"))
    L = length(H)
    gates = TrotterTSGate{paulistringtype(H),Float64,typeof(OrbitFlowCache(paulistringtype(H), typeof(H)))}[]
    length(H) == 0 && return gates

    if caches === nothing
        caches = [OrbitFlowCache(paulistringtype(H), typeof(H)) for _ in 1:L]
    end

    dt_eff = heisenberg ? dt : -dt

    if order == 1 || L == 1
        for j in 1:L
            p, c = H.strings[j], H.coeffs[j]
            push!(gates, TrotterTSGate(p, ComplexF64(c), Float64(dt_eff), caches[j]))
        end
    else
        for j in 1:(L-1)
            p, c = H.strings[j], H.coeffs[j]
            push!(gates, TrotterTSGate(p, ComplexF64(c), Float64(dt_eff / 2), caches[j]))
        end
        p, c = H.strings[L], H.coeffs[L]
        push!(gates, TrotterTSGate(p, ComplexF64(c), Float64(dt_eff), caches[L]))
        for j in (L-1):-1:1
            p, c = H.strings[j], H.coeffs[j]
            push!(gates, TrotterTSGate(p, ComplexF64(c), Float64(dt_eff / 2), caches[j]))
        end
    end
    return gates
end

"""
    ts_trotter_step!(O::Operator{<:PauliStringTS}, gates::AbstractVector{<:TrotterTSGate}; 
                  truncation=identity, componenttol=0.9999, maxlength=typemax(Int64), hbar=1)

Apply one translation-symmetric Trotter step in place. Iterates through the gates in reverse 
order, updating `O` by computing orbit flow updates.
"""
function ts_trotter_step!(O::Operator{<:PauliStringTS}, gates::AbstractVector{<:TrotterTSGate};
    truncation::Function=identity, componenttol::Real=0.9999,
    maxlength::Int=typemax(Int64), hbar::Real=1)
    isempty(gates) && return O
    for g in Iterators.reverse(gates)
        O2 = orbit_update(g, O, hbar, truncation; componenttol=componenttol, maxlength=maxlength)
        empty!(O.strings)
        empty!(O.coeffs)
        append!(O.strings, O2.strings)
        append!(O.coeffs, O2.coeffs)
    end
    return O
end

