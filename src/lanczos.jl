


function norm_lanczos(O::Operator)
    return opnorm(O) / sqrt(2^O.N)
end


"""
    lanczos(H::Operator, O::Operator, steps::Int, nterms::Int; keepnorm=true, maxlength=1000, returnOn=false)
    lanczos(H::OperatorTS1D, O::OperatorTS1D, steps::Int, nterms::Int; keepnorm=true, maxlength=1000, returnOn=false)

Compute the first `steps` lanczos coeficients for Hamiltonian `H` and initial operator `O`

At every step, the operator is trimed with [`trim`](@ref) and only `nterms` are kept.

Using `maxlength` speeds up the commutator by only keeping terms of length <= `maxlength`

Set `returnOn=false` to save the On's at each step. Then the function returns a pair of lists (bn, On).
The first operators of the list On is O
"""
function lanczos(H::Operator, O::Operator, steps::Int, nterms::Int; keepnorm=true, maxlength=1000, returnOn=false)
    @assert typeof(H) == typeof(O)
    O0 = deepcopy(O)
    O0 /= norm_lanczos(O0)
    LHO = com(H, O0)
    b = norm_lanczos(LHO)
    O1 = com(H, O0) / b
    bs = [b]
    returnOn && (Ons = [O0, O1])
    for n in ProgressBar(0:steps-2)
        LHO = com(H, O1; maxlength=maxlength)
        O2 = LHO - b * O0
        b = norm_lanczos(O2)
        O2 /= b
        O2 = trim(O2, nterms; keepnorm=keepnorm)
        returnOn && push!(Ons, O2)
        O0 = deepcopy(O1)
        O1 = deepcopy(O2)
        push!(bs, b)
    end
    returnOn && (return bs, Ons)
    return bs
end
