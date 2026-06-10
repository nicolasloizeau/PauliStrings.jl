norm_lanczos(O::AbstractOperator) = norm(O, normalize = true)

"""
    lanczos(H::Operator, O::Operator, steps::Int, nterms::Int; keepnorm=true, maxlength=1000, returnOn=false)

Compute the first `steps` lanczos coeficients for Hamiltonian `H` and initial operator `O`

At every step, the operator is trimed with [`trim`](@ref) and only `nterms` are kept.

Using `maxlength` speeds up the commutator by only keeping terms of length <= `maxlength`

Set `returnOn=true` to save the On's at each step. Then the function returns a pair of lists (bn, On).
The first operators of the list On is O
"""
function lanczos(H::AbstractOperator, O::AbstractOperator, steps::Int, nterms::Int; keepnorm = true, maxlength = 1000, returnOn = false, observer = false, show_progress = true)
    @assert typeof(H) == typeof(O)
    checklength(H, O)
    @assert observer === false || returnOn === false

    # Ôₙ₋₁ starts as Ô₀ = O / ‖O‖.
    Ôₙ₋₁ = scale(O, inv(norm_lanczos(O)))

    Ôₙ = commutator(H, Ôₙ₋₁)
    bₙ = norm_lanczos(Ôₙ)
    Ôₙ = scale!(Ôₙ, inv(bₙ))

    bs = [bₙ]
    returnOn && (Ons = [copy(Ôₙ₋₁), copy(Ôₙ)])
    (observer !== false) && (obs = [observer(Ôₙ₋₁), observer(Ôₙ)])

    progress = show_progress ? ProgressBar : identity
    for n in progress(0:(steps - 2))
        # Lanczos three-term recurrence:
        #   [H, Ôₙ] = bₙ₊₁ Ôₙ₊₁ + bₙ Ôₙ₋₁
        #   ⟹ Ôₙ₊₁ = ([H, Ôₙ] − bₙ Ôₙ₋₁) / bₙ₊₁
        bₙ₊₁Ôₙ₊₁ = commutator!(Ôₙ₋₁, H, Ôₙ, true, -bₙ; maxlength)
        bₙ₊₁ = norm_lanczos(bₙ₊₁Ôₙ₊₁)
        Ôₙ₊₁ = scale!(bₙ₊₁Ôₙ₊₁, inv(bₙ₊₁))
        Ôₙ₊₁ = trim(Ôₙ₊₁, nterms; keepnorm)

        returnOn && push!(Ons, copy(Ôₙ₊₁))
        (observer !== false) && push!(obs, observer(Ôₙ₊₁))
        push!(bs, bₙ₊₁)

        # advance the window for the next step: bₙ ← bₙ₊₁ and (Ôₙ₋₁, Ôₙ) ← (Ôₙ, Ôₙ₊₁);
        bₙ = bₙ₊₁
        Ôₙ₋₁, Ôₙ = Ôₙ, Ôₙ₊₁
    end

    (observer !== false) && return (bs, obs)
    returnOn && (return bs, Ons)
    return bs
end
