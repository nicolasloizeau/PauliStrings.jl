# TODO: This is actually a better implentation of binary_kernel, so even the RK4 should get a speedup. This is an optimization for a different PR

@inline _coord(site::Integer, Ls::Tuple) = Tuple(CartesianIndices(Ls)[site])
@inline _linear_site(coords::Tuple, Ls::Tuple) = LinearIndices(Ls)[coords...]

function _active_shift_tuple(site_p::Integer, site_q::Integer, Ls::Tuple, Ps::Tuple)
    cp = _coord(site_p, Ls)
    cq = _coord(site_q, Ls)
    valid = true
    shifts = ntuple(length(Ls)) do k
        if Ps[k]
            d = mod(cp[k] - cq[k], Ls[k])
            iszero(d) ? Ls[k] : d
        else
            valid &= cp[k] == cq[k]
            1
        end
    end
    return valid, shifts
end

function _ts_commutator!(d, prep::PauliString, c::Number, qrep::PauliString, Ls::Tuple, Ps::Tuple; epsilon::Real=0, maxlength::Int=1000)
    pauli_weight(prep) == 0 && return d
    pauli_weight(qrep) == 0 && return d

    N = Base.prod(Ls)
    seen = falses(N)

    for site_p in support(prep), site_q in support(qrep)
        valid, shifts = _active_shift_tuple(site_p, site_q, Ls, Ps)
        valid || continue
        key = _linear_site(shifts, Ls)
        seen[key] && continue
        seen[key] = true

        out, k = commutator(prep, shift(qrep, Ls, Ps, shifts))
        coeff = c * k
        if (k != 0) && (abs(coeff) > epsilon) && pauli_weight(out) < maxlength
            setwith!(+, d, PauliStringTS{Ls,Ps}(out), coeff)
        end
    end
    return d
end

function orbit_edges!(d, A::Operator{<:PauliStringTS}, q::PauliStringTS; epsilon::Real=0, maxlength::Int=1000)
    checklength(A, q)
    Ls = qubitsize(A)
    Ps = periodicflags(A)
    qrep = representative(q)
    pauli_weight(qrep) == 0 && return typeof(A)(collect(keys(d)), collect(values(d)))

    for (p, c) in zip(A.strings, A.coeffs)
        prep = representative(p) 
        _ts_commutator!(d, prep, c, qrep, Ls, Ps; epsilon=epsilon, maxlength=maxlength)
    end
    return d 
end

function orbit_liouvillian(A::Operator{<:PauliStringTS}, B::Operator{<:PauliStringTS}; hbar::Real=1, epsilon::Real=0, maxlength::Int=1000)
    checklength(A, B)
    return 1im * commutator(A, B; epsilon=epsilon, maxlength=maxlength) / hbar
end

struct _OrbitComponentPlan{P, S}
    component::Vector{P}
    index::Dict{P,Int}
    matrix::Matrix{ComplexF64}
    exp_cache::Dict{Float64,Matrix{ComplexF64}}
    signature::S
end

mutable struct _OrbitFlowCache{P,O}
    plans::Dict{P,_OrbitComponentPlan{P}}
    fallback::Dict{P,O}
end

_OrbitFlowCache(::Type{P}, ::Type{O}) where {P,O} = _OrbitFlowCache{P,O}(Dict{P,_OrbitComponentPlan{P}}(), Dict{P,O}())

function _orbit_component_and_transitions(Ha::Operator{<:PauliStringTS}, seed::PauliStringTS, hbar::Real, maxlength::Int)
    component = typeof(seed)[seed]
    index = Dict{typeof(seed),Int}(seed => 1)
    transitions = Tuple{Int,Int,ComplexF64}[]
    edges_d = emptydict(Ha)
    queue_index = 1

    while queue_index <= length(component)
        q = component[queue_index]
        j = queue_index
        queue_index += 1

        empty!(edges_d)
        orbit_edges!(edges_d, Ha, q; maxlength=maxlength)
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
