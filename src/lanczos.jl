norm_lanczos(O::AbstractOperator) = opnorm(O, normalize=true)
norm_lanczos(a::Number) = abs(a)


inner_lanczos(o1::AbstractOperator, o2::AbstractOperator) = trace_product(o1', o2; scale=1)
"""
    lanczos(H::Operator, O::Operator, steps::Int, nterms::Int; keepnorm=true, maxlength=1000, returnOn=false)
    lanczos(H::OperatorTS1D, O::OperatorTS1D, steps::Int, nterms::Int; keepnorm=true, maxlength=1000, returnOn=false)

Compute the first `steps` lanczos coeficients for Hamiltonian `H` and initial operator `O`

At every step, the operator is trimed with [`trim`](@ref) and only `nterms` are kept.

Using `maxlength` speeds up the commutator by only keeping terms of length <= `maxlength`

Set `returnOn=true` to save the On's at each step. Then the function returns a pair of lists (bn, On).
The first operators of the list On is O
"""
function lanczos(H::AbstractOperator, O::AbstractOperator, steps::Int, nterms::Int; keepnorm=true, maxlength=1000, returnOn=false, observer=false, show_progress=true)
    @assert typeof(H) == typeof(O)
    checklength(H, O)
    @assert observer === false || returnOn === false
    O0 = deepcopy(O)
    O0 /= norm_lanczos(O0)
    LHO = commutator(H, O0)
    b = norm_lanczos(LHO)
    O1 = commutator(H, O0) / b
    bs = [b]
    returnOn && (Ons = [O0, O1])
    (observer !== false) && (obs = [observer(O0), observer(O1)])
    progress = collect
    show_progress && (progress = ProgressBar)
    for n in progress(0:steps-2)
        LHO = commutator(H, O1; maxlength=maxlength)
        O2 = LHO - b * O0
        b = norm_lanczos(O2)
        O2 /= b
        O2 = trim(O2, nterms; keepnorm=keepnorm)
        returnOn && push!(Ons, O2)
        (observer !== false) && push!(obs, observer(O2))
        O0 = deepcopy(O1)
        O1 = deepcopy(O2)
        push!(bs, b)
    end
    (observer !== false) && return (bs, obs)
    returnOn && (return bs, Ons)
    return bs
end


function apply_lindblad(H, noise, O; adjoint=false)
    LH = commutator(H, O)
    LD = noise(O) - O
    if adjoint
        return LH - LD
    else
        return LH + LD
    end
end



"""
    bilanczos(H::Operator, O::Operator, steps::Int, nterms::Int, noise::Function; keepnorm=true, maxlength=1000, returnOn=false, observer=false, show_progress=true)

https://arxiv.org/pdf/1102.3909 fig 2
"""

function bilanczos(H::Operator, O::Operator, steps::Int, nterms::Int, noise::Function; maxlength=1000, returnOn=false, observer=false, show_progress=true)
    @assert typeof(H) == typeof(O)
    progress = collect
    show_progress && (progress = ProgressBar)
    a = Dict()
    b = Dict()
    c = Dict()
    u0 = O / norm_lanczos(O)
    w0 = O / norm_lanczos(O)

    b[0] = norm_lanczos(u0)
    c[0] = trace_product(w0', u0; scale=1)
    s = u0 / b[0]
    st = w0 / c[0]
    t = apply_lindblad(H, noise, s)
    tt = apply_lindblad(H, noise, st; adjoint=true)
    b[1] = 0
    c[1] = 0
    r = 0
    rt = 0
    for i in progress(1:steps)
        a[i] = trace_product(st', t; scale=1)
        t = t - a[i] * s - c[i] * r
        tt = tt - a[i]' * st - b[i]' * rt
        b[i+1] = norm_lanczos(t)
        c[i+1] = trace_product(tt', t; scale=1) / b[i+1]
        r = deepcopy(s)
        rt = deepcopy(st)
        s = t / b[i+1]
        st = tt / c[i+1]
        t = apply_lindblad(H, noise, s)
        tt = apply_lindblad(H, noise, st; adjoint=true)
        t = trim(t, nterms)
        tt = trim(tt, nterms)
    end
    a = [a[k] for k in sort(collect(keys(a)))]
    b = [b[k] for k in sort(collect(keys(b)))][3:end]
    c = [c[k] for k in sort(collect(keys(c)))][3:end]
    return a, b, c
end
