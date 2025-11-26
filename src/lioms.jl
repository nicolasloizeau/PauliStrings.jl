
@inline function f(H::T, O::T)::T where {T<:AbstractOperator}
    return im * commutator(H, O)
end

"""
    is_odd_under_spin_flip(op_list::Vector{Int}, time_reversal::Symbol)

Check whether a given operator (represented as a list of integers where 0=1,1=S+,2=Sz,3=S-)
is odd under spin flip.
"""
function is_odd_under_spin_flip(op_list::Vector{Int}, time_reversal::Symbol)
    cnt_z = count(x -> x == 2, op_list)
    sign = (-1)^cnt_z
    if time_reversal == :imag
        return sign == 1
    end
    return sign == -1
end

"""
    k_local_basis(N::Int, k::Int; translational_symmetry::Bool=false)::Vector{<:AbstractOperator}

Generates the basis of all k-site Pauli strings on N qubits, build from X,Y,Z and identity. 
If `translational_symmetry=true`, only the unit cell operators are returned as `OperatorTS1D`.
Otherwise, all translations of each operator are included as `Operator`s.
"""
function k_local_basis(N::Int, k::Int; translational_symmetry::Bool=false)::Vector{<:AbstractOperator}
    ops = translational_symmetry ? OperatorTS[] : Operator[]

    for i in 1:4^k-1
        op_list = digits(i, base=4, pad=k)

        # no identity on first site
        # since this will cause dupklicates after translations
        op_list[1] == 0 && continue

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


"""
    symmetry_adapted_k_local_basis(N::Int, k::Int; translational_symmetry::Bool=false)::Vector{<:AbstractOperator}

Generates the basis of all k-site Pauli strings on N qubits, build from S+,S-,Sz and identity. 
If `translational_symmetry=true`, only the unit cell operators are returned as `OperatorTS1D`.
Otherwise, all translations of each operator are included as `Operator`s. 

If `conserve_magnetization=:yes`, only operators such that [O, S^z_total] = 0 are included. 
If `conserve_magnetization=:no`, then only [O, S^z_total] != 0 are included. 
If `conserve_magnetization=:both`, no restriction is applied.

The spin-flip symmetry is defined via the parity operator P = ∏_i S^x_i. If `spin_flip=:even`, 
only operators such that POP = P are included. If `spin_flip=:odd`, then only POP = -O are included. 

We can form hermitian operators from the non-hermitian S+, S-, Sz basis in two ways, 
O + O† or i(O - O†). If `time_reversal=:real`, only the former are included. 
If `time_reversal=:imaginary`, only the latter are included. 


"""
function symmetry_adapted_k_local_basis(N::Int, k::Int; time_reversal::Symbol=:imag, spin_flip::Symbol=:even, conserve_magnetization::Symbol=:yes, translational_symmetry::Bool=true)
    ops = translational_symmetry ? OperatorTS[] : Operator[]
    base_ops = ["1", "S+", "Sz", "S-"]

    time_reversal ∈ (:real, :imag) || error("time_reversal must be :real or :imag")
    spin_flip ∈ (:even, :odd, :both) || error("spin_flip must be :even, :odd, or :both")
    conserve_magnetization ∈ (:yes, :no, :both) || error("conserve_magnetization must be :yes, :no, or :both")

    for i in 1:4^k-1
        op_list = digits(i, base=4, pad=k)

        # no identity on first site
        op_list[1] == 0 && continue

        # all checks passed, building operator
        op = Operator(N)

        # check if hermitian conjugate already included
        # since we always build ops O+O† or/and i(O-O†)
        op_list_hc = map(x -> x == 1 ? 3 : x == 3 ? 1 : x, op_list)
        op_hc = sum(op_list_hc[i] * 4^(i - 1) for i in 1:length(op_list_hc))
        op_hc < i && continue

        # check for spin flip symmetry
        (spin_flip == :even) && is_odd_under_spin_flip(op_list, time_reversal) && continue
        (spin_flip == :odd) && !is_odd_under_spin_flip(op_list, time_reversal) && continue

        # check for magnetization conservation
        magnetization_check = count(x -> x == 1, op_list) - count(x -> x == 3, op_list)
        (conserve_magnetization == :yes) && magnetization_check != 0 && continue
        (conserve_magnetization == :no) && magnetization_check == 0 && continue

        terms = vcat([[base_ops[o+1], idx] for (idx, o) in enumerate(op_list) if o != 0]...)
        op += (1.0, terms...)

        if time_reversal == :real
            op += dagger(op)
        end

        if time_reversal == :imag
            count(x -> x in (1, 3), op_list) == 0 && continue  # skip if no S+ or S- present
            op -= dagger(op)
            op *= im
        end

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


"""
    lioms(H::T, support::Vector{T}; threshold::Real=1e-14, f::Function=f)::Tuple{Vector{Float64},Vector{T}} where {T<:AbstractOperator}

Algorithm constructing all Local Integrals of Motion (LIOMs) for a Hamiltonian H,
supported on the operator basis from `support`.
Follows definitions in [https://arxiv.org/abs/2505.05882](https://arxiv.org/abs/2505.05882).
'threshold` parameter sets the threshold for eigenvalues above which eigenmodes are discarded. 
By defualt it is 0.0, meaning only exact LIOMs are returned.
Uses a function `f(H,O)` to check for LIOMs, by default the commutator `f(H,O) = im*[H,O]`.
Returns a tuple of eigenvalues and eigenmodes (operators).
"""
function lioms(H::T, support::Vector{T}; threshold::Real=1e-14, f::Function=f)::Tuple{Vector{Float64},Vector{T}} where {T<:AbstractOperator}

    n = length(support)
    L = qubitlength(H)
    scale = isa(H, OperatorTS) ? 1 / (2^L * L) : 1 / (2^L)
    support ./= [opnorm(O; normalize=true) for O in support]
    fs = [f(H, O) for O in support]

    Fmat = zeros(Float64, n, n)
    for i in 1:n
        for j in i:n
            Fmat[i, j] = trace_product(dagger(fs[i]), fs[j]) * scale
        end
    end

    evals, evecs = eigen(Symmetric(Fmat, :U))
    num_to_return = count(x -> x <= threshold, evals)
    ops = similar(support, num_to_return)
    for i in 1:num_to_return
        ops[i] = cutoff(sum(evecs[:, i] .* support), 1e-10)
    end
    return evals[1:num_to_return], ops
end



"""
    lioms(H::T, k::Int; threshold::Real=1e-14, f::Function=f)::Tuple{Vector{Float64},Vector{T}} where {T<:AbstractOperator}

Algorithm constructing all Local Integrals of Motion (LIOMs) for a Hamiltonian H, supported on the 
most general pauli string basis on `k` sites.
Follows definitions in [https://arxiv.org/abs/2505.05882](https://arxiv.org/abs/2505.05882).
`threshold` parameter sets the threshold for eigenvalues above which eigenmodes are discarded. 
By defualt it is 0.0, meaning only exact LIOMs are returned.
Uses a function `f(H,O)` to check for LIOMs, by default the commutator `f(H,O) = im*[H,O]`.
Returns a tuple of eigenvalues and eigenmodes (operators).
"""
function lioms(H::T, k::Int; threshold::Real=1e-14, f::Function=f)::Tuple{Vector{Float64},Vector{T}} where {T<:AbstractOperator}
    N = qubitlength(H)
    ts = isa(H, OperatorTS)
    support = k_local_basis(N, k; translational_symmetry=ts)
    return lioms(H, support; threshold=threshold, f=f)
end