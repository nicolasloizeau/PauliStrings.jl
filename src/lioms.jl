
function f(H::Operator, O::Operator)::Operator
    return im * commutator(H, O)
end


function lioms(H::Operator, support::Vector{Operator}; return_all::Bool=false)::Tuple{Vector{Float64},Vector{Operator}}

    n = length(support)
    support ./= [opnorm(O, normalize=true) for O in support]
    fs = [f(H, O) for O in support]

    Fmat = zeros(Float64, n, n)
    for i in 1:n
        for j in i:n
            Fmat[i, j] = trace_product(fs[i], fs[j]; scale=1)
            i != j && (Fmat[j, i] = Fmat[i, j])
        end
    end

    evals, evecs = eigen(Symmetric(Fmat))

    if return_all
        ops = Vector{Operator}(undef, n)
        for i in 1:n
            ops[i] = cutoff(sum(evecs[:, i] .* support), 1e-10)
        end
        return evals, ops
    else
        n_lioms = count(x -> isapprox(x, 0; atol=1e-10), evals)
        ops = Vector{Operator}(undef, n_lioms)
        for i in 1:n_lioms
            ops[i] = cutoff(sum(evecs[:, i] .* support), 1e-10)
        end
        return evals[1:n_lioms], ops
    end

end