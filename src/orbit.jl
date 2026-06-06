function _coord(site::Integer, Ls::Tuple, dim::Integer)
    x = site - 1
    for k in 1:(dim - 1)
        x = div(x, Ls[k])
    end
    return mod(x, Ls[dim]) + 1
end

function _linear_site(coords::Tuple, Ls::Tuple)
    site = 1
    stride = 1
    for k in eachindex(Ls)
        site += (coords[k] - 1) * stride
        stride *= Ls[k]
    end
    return site
end

function _active_shift_tuple(site_p::Integer, site_q::Integer, Ls::Tuple, Ps::Tuple)
    valid = true
    shifts = ntuple(length(Ls)) do k
        cp = _coord(site_p, Ls, k)
        cq = _coord(site_q, Ls, k)
        if Ps[k]
            d = mod(cp - cq, Ls[k])
            iszero(d) ? Ls[k] : d
        else
            valid &= cp == cq
            1
        end
    end
    return valid, shifts
end

function orbit_edges(A::Operator{<:PauliStringTS}, q::PauliStringTS; epsilon::Real=0, maxlength::Int=1000)
    checklength(A, q)
    Ls = qubitsize(A)
    Ps = periodicflags(A)
    N = Base.prod(Ls)
    d = emptydict(A)
    qrep = representative(q)
    uq = qrep.v | qrep.w
    iszero(uq) && return typeof(A)(collect(keys(d)), collect(values(d)))

    for (p, c) in zip(A.strings, A.coeffs)
        prep = representative(p)
        up = prep.v | prep.w
        iszero(up) && continue
        seen = falses(N)
        xp = up
        while !iszero(xp)
            site_p = trailing_zeros(xp) + 1
            xp &= xp - one(xp)

            xq = uq
            while !iszero(xq)
                site_q = trailing_zeros(xq) + 1
                xq &= xq - one(xq)

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
        end
    end
    return typeof(A)(collect(keys(d)), collect(values(d)))
end

function orbit_liouvillian(A::Operator{<:PauliStringTS}, B::Operator{<:PauliStringTS}; hbar::Real=1, epsilon::Real=0, maxlength::Int=1000)
    checklength(A, B)
    return 1im * commutator(A, B; epsilon=epsilon, maxlength=maxlength) / hbar
end

struct _OrbitComponentPlan{P}
    component::Vector{P}
    index::Dict{P,Int}
    matrix::Matrix{ComplexF64}
    exp_cache::Dict{Float64,Matrix{ComplexF64}}
end

mutable struct _OrbitFlowCache{P,O}
    plans::Dict{P,_OrbitComponentPlan{P}}
    fallback::Dict{P,O}
end

_OrbitFlowCache(::Type{P}, ::Type{O}) where {P,O} = _OrbitFlowCache{P,O}(Dict{P,_OrbitComponentPlan{P}}(), Dict{P,O}())

function _orbit_component(Ha::Operator{<:PauliStringTS}, seed::PauliStringTS, maxlength::Int)
    component = typeof(seed)[seed]
    index = Dict{typeof(seed),Int}(seed => 1)
    queue_index = 1
    while queue_index <= length(component)
        q = component[queue_index]
        queue_index += 1
        edges = orbit_edges(Ha, q; maxlength=maxlength)
        for r in edges.strings
            haskey(index, r) && continue
            push!(component, r)
            index[r] = length(component)
        end
    end
    return component, index
end

function _orbit_component_matrix(Ha::Operator{<:PauliStringTS}, component, index, hbar::Real, maxlength::Int)
    n = length(component)
    A = zeros(ComplexF64, n, n)
    for (j, q) in enumerate(component)
        edges = orbit_edges(Ha, q; maxlength=maxlength)
        for (r, c) in zip(edges.strings, edges.coeffs)
            i = get(index, r, 0)
            i == 0 && throw(ArgumentError("orbit component is not closed"))
            A[i, j] += 1im * c / hbar
        end
    end
    return A
end

function _component_plan(cache::_OrbitFlowCache, Ha::Operator{<:PauliStringTS}, seed::PauliStringTS, hbar::Real, maxlength::Int)
    plan = get(cache.plans, seed, nothing)
    plan !== nothing && return plan

    component, index = _orbit_component(Ha, seed, maxlength)
    A = _orbit_component_matrix(Ha, component, index, hbar, maxlength)
    plan = _OrbitComponentPlan(component, index, A, Dict{Float64,Matrix{ComplexF64}}())
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

function _component_signature(A::AbstractMatrix{<:Number}; digits::Int=12)
    entries = Tuple{Int,Int,Float64,Float64}[]
    for j in axes(A, 2), i in axes(A, 1)
        v = A[i, j]
        iszero(v) && continue
        push!(entries, (i, j, round(real(v), digits=digits), round(imag(v), digits=digits)))
    end
    return (size(A, 1), entries)
end

function _orbit_terms(H::Operator{<:PauliStringTS})
    terms = typeof(H)[]
    for j in 1:length(H)
        push!(terms, typeof(H)([H.strings[j]], [H.coeffs[j]]))
    end
    return terms
end
