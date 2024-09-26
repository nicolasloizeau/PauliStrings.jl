


function norm_lanczos(O::Operator)
    return opnorm(O)/sqrt(2^O.N)
end


"
    lanczos(H::Operator, O::Operator, steps::Int, nterms::Int; keepnorm=true, maxlength=1000, localop=false)

Computer the first `steps` lanczos coeficients for Hamiltonian `H` and initial operator `O`

At every step, the operator is trimed with `PauliStrings.trim` and only `nterms` are kept.

If H and O are 1D-translation-invariant, it is possible to provide a single local term of O and set `localop=true`

Using `maxlength` speeds up the commutator by only keeping terms of length <= `maxlength`
"
function lanczos(H::Operator, O::Operator, steps::Int, nterms::Int; keepnorm=true, maxlength=1000, localop=false)
    c = 1
    N = H.N
    localop && (c = sqrt(N))
    O0 = deepcopy(O)
    localop && (O0 = shift_left(O0))
    O0 /= norm_lanczos(O0)*c
    LHO = com(H, O0)
    localop && (LHO = shift_left(LHO))
    b = norm_lanczos(LHO)*c
    O1 = com(H, O0)/b
    localop && (O1 = shift_left(O1))
    bs = [b]
    for n in ProgressBar(0:steps-2)
        LHO = com(H, O1; maxlength=maxlength)
        localop && (LHO = shift_left(LHO))
        O2 = LHO-b*O0
        localop && (O2 = shift_left(O2))
        b = norm_lanczos(O2)*c
        O2 /= b
        O2 = trim(O2, nterms; keepnorm=keepnorm)
        O0 = deepcopy(O1)
        O1 = deepcopy(O2)
        push!(bs, b)
    end
    return bs
end
