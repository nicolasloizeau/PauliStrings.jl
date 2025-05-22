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
    # O = commutator(H, O)
    if adjoint
        return -im*commutator(H, O)
    else
        return im*commutator(H, O)
    end
    # if adjoint
    #     return O + im * noise(O)
    # end
    # return O - im * noise(O)
    # return O
end



"""
https://arxiv.org/pdf/2405.09628 below eq 250
"""
function bilanczos2(H::Operator, O::Operator, steps::Int, nterms::Int, noise::Function; keepnorm=true, maxlength=1000, returnOn=false, observer=false, show_progress=true)
    @assert typeof(H) == typeof(O)
    checklength(H, O)
    progress = collect
    show_progress && (progress = ProgressBar)
    pm = 0 + 0im
    qm = 0 + 0im
    b = 0 + 0im
    c = 0 + 0im
    p = O / norm_lanczos(O)
    q = O / norm_lanczos(O)
    as = Complex{Float64}[]
    bs = Complex{Float64}[]
    cs = Complex{Float64}[]

    for j in progress(0:steps)
        r = apply_lindblad(H, noise, p; adjoint=true)
        s = apply_lindblad(H, noise, q)
        r = r - b * pm
        s = s - c' * qm
        a = trace_product(dagger(q), r; scale=1)
        r = r - a * p
        s = s - a' * q
        omega = trace_product(dagger(r), s; scale=1)
        cp = sqrt(abs(omega))
        bp = omega' / c
        pp = r / cp
        qp = s / bp'
        # pp = trim(pp, nterms)
        # qp = trim(qp, nterms)

        pm = deepcopy(p)
        p = deepcopy(pp)
        qm = deepcopy(q)
        q = deepcopy(qp)

        c = cp
        b = bp

        push!(as, a)
        push!(bs, b)
        push!(cs, c)

    end
    return as, bs, cs

end


function bilanczos(H::Operator, O::Operator, steps::Int, nterms::Int, noise::Function; keepnorm=true, maxlength=1000, returnOn=false, observer=false, show_progress=true)
    @assert typeof(H) == typeof(O)
    checklength(H, O)
    progress = collect
    show_progress && (progress = ProgressBar)
    # O = O / norm_lanczos(O)
    # P = O
    O1 = O / norm_lanczos(O)
    O2 = 0
    P1 = apply_lindblad(H, noise, O1)
    P1 = P1 - inner_lanczos(P1, O1)*O1
    P1 = P1 / norm_lanczos(P1)
    println(norm_lanczos(P1))
    println(inner_lanczos(P1, O1))
    # P1 = O / norm_lanczos(O)
    P2 = 0
    a1 = 0
    b1 = 0
    c1 = 0
    as = []
    bs = []
    cs = []
    for n in 1:steps
        A = apply_lindblad(H, noise, O1)-a1*O1-c1*O2
        B = apply_lindblad(H, noise, O1; adjoint=true)-a1'*O1-b1*O2
        b = norm_lanczos(A)
        c = inner_lanczos(B, A)/b
        O = A/b
        P = B/c
        a = inner_lanczos(P, apply_lindblad(H, noise, O))
        P1 = deepcopy(P)
        P2 = deepcopy(P1)
        O1 = deepcopy(O)
        O2 = deepcopy(O1)
        a1 = a
        b1 = b
        c1 = c
        push!(as, a)
        push!(bs, b)
        push!(cs, c)
        println(a, " ", b, " ", c)
    end
    return as, bs, cs
end
