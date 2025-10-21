
@inline function f(H::T, O::T)::T where {T<:AbstractOperator}
    return im * commutator(H, O)
end

function k_local_basis(N::Int, k::Int; translational_symmetry::Bool=false)::Vector{<:AbstractOperator}
    ops = translational_symmetry ? OperatorTS[] : Operator[]

    for i in 1:4^k-1
        op_list = digits(i, base=4, pad=k)

        # no identity on first site
        op_list[1] == 0 && continue

        # all checks passed, building operator
        p = zeros(Int, N)
        sites = 1:k
        p[sites] .= op_list
        op = Operator(N)
        op += string_from_inds(p)

        if translational_symmetry
            push!(ops, OperatorTS1D(op, full=false))
        else
            for s in 0:N-1
                shifted_op = shift(op, s)
                push!(ops, shifted_op)
            end
        end
    end

    return ops
end


function lioms(H::T, support::Vector{T}; return_all::Bool=false, f::Function=f)::Tuple{Vector{Float64},Vector{T}} where {T<:AbstractOperator}

    n = length(support)
    support ./= [opnorm(O, normalize=true) for O in support]
    fs = [f(H, O) for O in support]

    Fmat = zeros(Float64, n, n)
    for i in 1:n
        for j in i:n
            Fmat[i, j] = trace_product(fs[i], fs[j]; scale=1)
        end
    end

    evals, evecs = eigen(Symmetric(Fmat, :U))

    if return_all
        ops = similar(support, n)
        for i in 1:n
            ops[i] = cutoff(sum(evecs[:, i] .* support), 1e-10)
        end
        return evals, ops
    else
        n_lioms = count(x -> isapprox(x, 0; atol=1e-10), evals)
        ops = similar(support, n_lioms)
        for i in 1:n_lioms
            ops[i] = cutoff(sum(evecs[:, i] .* support), 1e-10)
        end
        return evals[1:n_lioms], ops
    end

end

function lioms(H::T, k::Int; return_all::Bool=false, f::Function=f)::Tuple{Vector{Float64},Vector{T}} where {T<:AbstractOperator}
    N = qubitlength(H)
    ts = isa(H, OperatorTS)
    support = k_local_basis(N, k; translational_symmetry=ts)
    return lioms(H, support; return_all=return_all, f=f)
end