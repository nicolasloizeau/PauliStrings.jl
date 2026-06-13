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

function orbit_edges(A::Operator{<:PauliStringTS}, q::PauliStringTS; epsilon::Real=0, maxlength::Int=1000)
    d = emptydict(A)
    seen = falses(Base.prod(qubitsize(A)))
    _orbit_edges!(d, A, q; seen, maxlength)
    o = typeof(A)(collect(keys(d)), collect(values(d)))
    return (epsilon > 0) ? cutoff(o, epsilon) : o
end

function _orbit_edges!(d, A::Operator{PauliStringTS{Ls,Ps,U}}, q::PauliStringTS{Ls,Ps,U}; seen, maxlength::Int=1000) where {Ls,Ps,U}
    checklength(A, q)
    qrep = representative(q)
    pauli_weight(qrep) == 0 && return typeof(A)(collect(keys(d)), collect(values(d)))

    for (p, c) in zip(A.strings, A.coeffs)
        prep = representative(p)
        _ts_commutator!(d, prep, c, qrep, Ls, Ps; maxlength)
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
    sort!(entries)
    return (n, entries)
end

function _orbit_terms(H::Operator{<:PauliStringTS})
    terms = typeof(H)[]
    for j in 1:length(H)
        push!(terms, typeof(H)([H.strings[j]], [H.coeffs[j]]))
    end
    return terms
end
