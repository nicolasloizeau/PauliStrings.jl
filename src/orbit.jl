# TODO: This is actually a better implentation of binary_kernel, so even the RK4 should get a speedup. This is an optimization for a different PR

function _active_shift_tuple(site_p::Integer, site_q::Integer, Ls, Ps)
    R = CartesianIndices(Ls)
    cp = R[site_p]
    cq = R[site_q]

    for k in 1:length(Ls)
        if !Ps[k] && cp[k] != cq[k]
            return false, Ls, 1
        end
    end
    key = 1
    stride = 1
    shifts = ntuple(length(Ls)) do k
        if Ps[k]
            d = mod(cp[k] - cq[k], Ls[k])
            sk = iszero(d) ? Ls[k] : d
            key += (sk - 1) * stride
            stride *= Ls[k]
            sk
        else
            stride *= Ls[k]
            1
        end
    end
    return true, shifts, key
end

function _ts_commutator!(d, prep::PauliString, c::Number, qrep::PauliString, Ls, Ps; seen=Nothing, maxlength::Int=1000)
    pauli_weight(prep) == 0 && return d
    pauli_weight(qrep) == 0 && return d

    if isnothing(seen)
        seen = falses(N)
    else
        fill!(seen, false)
    end

    sup_p = support(prep)
    sup_q = support(qrep)
    for site_p in sup_p, site_q in sup_q
        valid, shifts, key = _active_shift_tuple(site_p, site_q, Ls, Ps)
        valid || continue
        seen[key] && continue
        seen[key] = true

        out, k = commutator(prep, shift(qrep, Ls, Ps, shifts))
        coeff = c * k
        if (k != 0) && pauli_weight(out) < maxlength
            setwith!(+, d, PauliStringTS{Ls,Ps}(out), coeff)
        end
    end
    return d
end

function orbit_edges(A::Operator{<:PauliStringTS}, q::PauliStringTS; epsilon::Real=0, maxlength::Int=1000)
    d = emptydict(A)
    seen = falses(Base.prod(qubitsize(A)))
    _orbit_edges!(d, A, q; seen, maxlength)
    o = typeof(A)(collect(keys(d)), collect(values(d)))
    return (epsilon > 0) ? cutoff(o, epsilon) : o
end

function _orbit_edges!(d, A::Operator{PauliStringTS{Ls, Ps, U}}, q::PauliStringTS{Ls, Ps, U}; seen, maxlength::Int=1000) where {Ls, Ps, U}
    checklength(A, q)
    qrep = representative(q)
    pauli_weight(qrep) == 0 && return typeof(A)(collect(keys(d)), collect(values(d)))

    for (p, c) in zip(A.strings, A.coeffs)
        prep = representative(p)
        _ts_commutator!(d, prep, c, qrep, Ls, Ps; seen, maxlength)
    end
    return d
end

function orbit_liouvillian(A::Operator{<:PauliStringTS}, B::Operator{<:PauliStringTS}; hbar::Real=1, epsilon::Real=0, maxlength::Int=1000)
    checklength(A, B)
    return 1im * commutator(A, B; epsilon=epsilon, maxlength=maxlength) / hbar
end

struct _OrbitComponentPlan{P,S}
    component::Vector{P}
    index::Dict{P,Int}
    matrix::Matrix{ComplexF64}
    exp_cache::Dict{Float64,Matrix{ComplexF64}}
    signature::S
end

mutable struct _OrbitFlowCache{P,O,C,D,E}
    plans::Dict{P,_OrbitComponentPlan{P}}
    fallback::Dict{P,O}
    assigned::Set{P}
    component_data::Vector{Tuple{Float64,Vector{P},Any,Vector{ComplexF64}}}
    out_d::C
    coeff_lookup::D
    edges_d::E
    seen::BitVector
end

function _OrbitFlowCache(::Type{P}, ::Type{O}) where {P,O}
    dummy = O(P[], ComplexF64[])
    out_d = emptydict(dummy)
    coeff_lookup = Dict{P,ComplexF64}()
    edges_d = emptydict(dummy)
    seen = falses(Base.prod(qubitsize(dummy)))
    return _OrbitFlowCache(
        Dict{P,_OrbitComponentPlan{P}}(),
        Dict{P,O}(),
        Set{P}(),
        Tuple{Float64,Vector{P},Any,Vector{ComplexF64}}[],
        out_d,
        coeff_lookup,
        edges_d,
        seen
    )
end

function _orbit_component_and_transitions(Ha::Operator{<:PauliStringTS}, seed::PauliStringTS, hbar::Real, cache::_OrbitFlowCache, maxlength::Int)
    component = typeof(seed)[seed]
    index = Dict{typeof(seed),Int}(seed => 1)
    transitions = Tuple{Int,Int,ComplexF64}[]
    edges_d = cache.edges_d
    seen = cache.seen
    queue_index = 1

    while queue_index <= length(component)
        q = component[queue_index]
        j = queue_index
        queue_index += 1

        empty!(edges_d)
        _orbit_edges!(edges_d, Ha, q; seen, maxlength)
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

function _component_plan!(cache::_OrbitFlowCache, component, index, transitions)
    # Return already cached plan if another node in this component populated it first
    plan = get(cache.plans, component[1], nothing)
    plan !== nothing && return plan

    n = length(component)
    A = zeros(ComplexF64, n, n)
    for (i, j, val) in transitions # A is small
        A[i, j] += val
    end
    sig = _compute_signature(transitions, n)
    plan = _OrbitComponentPlan(component, index, A, Dict{Float64,Matrix{ComplexF64}}(), sig)
    for q in component
        cache.plans[q] = plan
    end
    return plan
end

function _component_exp(plan::_OrbitComponentPlan, dt::Real)
    key = Float64(dt)
    E = get(plan.exp_cache, key, nothing)
    if E === nothing
        E = exp(key * plan.matrix)
        plan.exp_cache[key] = E
    end
    return E
end

function _compute_signature(transitions::Vector{Tuple{Int,Int,ComplexF64}}, n::Int; digits::Int=12)
    entries = Tuple{Int,Int,ComplexF64}[]
    sizehint!(entries, length(transitions))
    for (i, j, val) in transitions
        iszero(val) && continue
        push!(entries, (i, j, round(val, digits=digits)))
    end
    return (n, entries)
end

function _orbit_terms(H::Operator{<:PauliStringTS})
    terms = typeof(H)[]
    for j in 1:length(H)
        push!(terms, typeof(H)([H.strings[j]], [H.coeffs[j]]))
    end
    return terms
end
