


function norm_lanczos(O)
    return opnorm(O)/sqrt(2^O.N)
end


"
    lanczos(H::Operator, O::Operator, steps::Int, nterms::Int; keepnorm=true, maxlength=1000)
    lanczos(H::OperatorTS1D, O::OperatorTS1D, steps::Int, nterms::Int; keepnorm=true, maxlength=1000)

Computer the first `steps` lanczos coeficients for Hamiltonian `H` and initial operator `O`

At every step, the operator is trimed with [`trim`](@ref) and only `nterms` are kept.

Using `maxlength` speeds up the commutator by only keeping terms of length <= `maxlength`
"
function lanczos(H::Union{Operator, OperatorTS1D}, O::Union{Operator, OperatorTS1D}, steps::Int, nterms::Int; keepnorm=true, maxlength=1000)
    @assert typeof(H) == typeof(O)
    O0 = deepcopy(O)
    O0 /= norm_lanczos(O0)
    LHO = com(H, O0)
    b = norm_lanczos(LHO)
    O1 = com(H, O0)/b
    bs = [b]
    for n in ProgressBar(0:steps-2)
        LHO = com(H, O1; maxlength=maxlength)
        O2 = LHO-b*O0
        b = norm_lanczos(O2)
        O2 /= b
        O2 = trim(O2, nterms; keepnorm=keepnorm)
        O0 = deepcopy(O1)
        O1 = deepcopy(O2)
        push!(bs, b)
    end
    return bs
end
