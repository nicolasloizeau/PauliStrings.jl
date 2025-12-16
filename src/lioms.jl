
@inline function f(H::T, O::T)::T where {T<:AbstractOperator}
    return im * commutator(H, O)
end

"""
    is_odd_under_spin_flip(op_list::Vector{Int}, time_reversal::Symbol)

Check whether a given operator is odd under spin flip.

# Arguments
- `op_list::Vector{Int}`: List of integers representing the operator (0=1, 1=S+, 2=Sz, 3=S-)
- `time_reversal::Symbol`: Time reversal symmetry setting (`:real` or `:imag`)

# Returns
- `Bool`: `true` if the operator is odd under spin flip, `false` otherwise
"""
function is_odd_under_spin_flip(op_list::Vector{Int}, time_reversal::Symbol)
    cnt_z = count(==(2), op_list)
    sign = (-1)^cnt_z
    return time_reversal == :imag ? (sign == 1) : (sign == -1)
end

"""
    k_local_basis_1d(N::Int, k::Int; translational_symmetry::Bool=false)::Vector{<:AbstractOperator}

Generates the basis of all k-site Pauli strings on N qubits, build from X,Y,Z and identity.

# Arguments
- `N::Int`: Number of qubits
- `k::Int`: Maximum support size (number of sites)
- `translational_symmetry::Bool=false`: If `true`, return only unit cell operators as `OperatorTS1D`; otherwise return all translations as `Operator`s

# Returns
- `Vector{<:AbstractOperator}`: Basis of k-site Pauli strings
"""
function k_local_basis_1d(N::Int, k::Int; translational_symmetry::Bool=false)::Vector{<:AbstractOperator}
    strings = translational_symmetry ? PauliStringTS[] : PauliString[]

    @inbounds for i in 1:4^k-1
        op_list = digits(i, base=4, pad=k)

        # no identity on first site
        # since this will cause dupklicates after translations
        op_list[1] == 0 && continue

        p = zeros(Int, N)
        sites = 1:k
        p[sites] .= op_list
        string = string_from_inds(p)
        string = PauliString{N}(string)

        if translational_symmetry
            push!(strings, PauliStringTS{(N,)}(string))
        else
            for s in 0:N-1
                push!(strings, shift(string, s))
            end
        end
    end

    return strings
end


"""
    symmetry_adapted_k_local_basis_1d(N::Int, k::Int; time_reversal::Symbol=:imag, spin_flip::Symbol=:even, conserve_magnetization::Symbol=:yes, translational_symmetry::Bool=true)

Generates the basis of all k-site Pauli strings on N qubits, build from S+,S-,Sz and identity.
If `translational_symmetry=true`, only the unit cell operators are returned as `OperatorTS1D`.
Otherwise, all translations of each operator are included as `Operator`s.

The spin-flip symmetry is defined via the parity operator P = ∏_i S^x_i. If `spin_flip=:even`,
only operators such that POP = P are included. If `spin_flip=:odd`, then only POP = -O are included.

We can form hermitian operators from the non-hermitian S+, S-, Sz basis in two ways,
O + O† or i(O - O†). If `time_reversal=:real`, only the former are included.
If `time_reversal=:imaginary`, only the latter are included.

# Arguments
- `N::Int`: Number of qubits
- `k::Int`: Maximum support size (number of sites)
- `time_reversal::Symbol=:imag`: Time reversal symmetry (`:real` for O + O†, `:imag` for i(O - O†))
- `spin_flip::Symbol=:even`: Spin-flip parity (`:even` for POP = O, `:odd` for POP = -O, `:both` for no restriction)
- `conserve_magnetization::Symbol=:yes`: Magnetization conservation (`:yes` for [O, S^z_total] = 0, `:no` for [O, S^z_total] ≠ 0, `:both` for no restriction)
- `translational_symmetry::Bool=true`: If `true`, return only unit cell operators as `OperatorTS1D`; otherwise return all translations

# Returns
- `Vector{<:AbstractOperator}`: Symmetry-adapted basis of k-site Pauli strings
"""
function symmetry_adapted_k_local_basis_1d(N::Int, k::Int; time_reversal::Symbol=:imag, spin_flip::Symbol=:even, conserve_magnetization::Symbol=:yes, translational_symmetry::Bool=true)::Vector{<:AbstractOperator}
    ops = translational_symmetry ? OperatorTS[] : Operator[]
    base_ops = ["1", "S+", "Sz", "S-"]

    time_reversal ∈ (:real, :imag) || error("time_reversal must be :real or :imag")
    spin_flip ∈ (:even, :odd, :both) || error("spin_flip must be :even, :odd, or :both")
    conserve_magnetization ∈ (:yes, :no, :both) || error("conserve_magnetization must be :yes, :no, or :both")

    @inbounds for i in 1:4^k-1
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
        magnetization_check = count(==(1), op_list) - count(==(3), op_list)
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
    lioms(H::AbstractOperator, support::Vector{T}; threshold::Real=1e-14, f::Function=f) where {T<:AbstractOperator}

Algorithm constructing all Local Integrals of Motion (LIOMs) for a Hamiltonian H,
supported on the operator basis from `support`.
Follows definitions in [https://arxiv.org/abs/2505.05882](https://arxiv.org/abs/2505.05882).

# Arguments
- `H::T`: The Hamiltonian operator
- `support::Vector{T}`: Vector of basis operators
- `threshold::Real=1e-14`: Threshold for eigenvalues above which eigenmodes are discarded. By default only exact LIOMs are returned.
- `f::Function=f`: Function to check for LIOMs, by default the commutator `f(H,O) = im*[H,O]`

# Returns
- `evals::Vector{Float64}`: Eigenvalues below threshold
- `evecs::Matrix{Float64}`: Eigenvectors corresponding to eigenvalues
- `ops::Vector{AbstractOperator}`: LIOM operators
"""
function lioms(H::AbstractOperator, support::Vector{T}; threshold::Real=1e-14, f::Function=f)::Tuple{Vector{Float64},Matrix{Float64},Vector{AbstractOperator}} where {T<:AbstractOperator}
    n = length(support)
    support = convert(Vector{Any}, support)
    L = qubitlength(H)
    scale = isa(H, OperatorTS) ? 1 / (2^L * L) : 1 / (2^L)

    norms = map(op -> norm(op; normalize=true), support)
    @inbounds for i in 1:n
        support[i] /= norms[i]
    end

    fs = Vector{AbstractOperator}(undef, n)
    @inbounds for i in 1:n
        fs[i] = f(H, support[i])
    end
    dagger_fs = map(adjoint, fs)

    Fmat = zeros(Float64, n, n)
    @inbounds begin
        if n >= 128 # not sure where is the threshold for threading
            @static if VERSION >= v"1.11"
                Threads.@threads :greedy for i in 1:n
                    di = dagger_fs[i]
                    for j in i:n
                        Fmat[i, j] = trace_product(di, fs[j]) * scale
                    end
                end
            else
                # default schedule, will be :dynamic for 1.8-1.10 and :static for older
                Threads.@threads for i in 1:n
                    di = dagger_fs[i]
                    for j in i:n
                        Fmat[i, j] = trace_product(di, fs[j]) * scale
                    end
                end
            end
        else
            for i in 1:n
                di = dagger_fs[i]
                for j in i:n
                    Fmat[i, j] = trace_product(di, fs[j]) * scale
                end
            end
        end
    end

    evals, evecs = eigen(Symmetric(Fmat, :U))
    num_to_return = count(<=(threshold), evals)
    ops = similar(support, num_to_return)
    @inbounds for i in 1:num_to_return
        ops[i] = cutoff(sum(evecs[:, i] .* support), 1e-10)
    end
    return evals[1:num_to_return], evecs[:, 1:num_to_return], ops
end



"""
    lioms(H::T, k::Int; threshold::Real=1e-14, f::Function=f) where {T<:AbstractOperator}

Algorithm constructing all Local Integrals of Motion (LIOMs) for a Hamiltonian H, supported on the
most general Pauli string basis on `k` sites.
Follows definitions in [https://arxiv.org/abs/2505.05882](https://arxiv.org/abs/2505.05882).

# Arguments
- `H::T`: The Hamiltonian operator
- `k::Int`: Maximum support size for basis operators
- `threshold::Real=1e-14`: Threshold for eigenvalues above which eigenmodes are discarded. By default only exact LIOMs are returned.
- `f::Function=f`: Function to check for LIOMs, by default the commutator `f(H,O) = im*[H,O]`

# Returns
- `evals::Vector{Float64}`: Eigenvalues below threshold
- `evecs::Matrix{Float64}`: Eigenvectors corresponding to eigenvalues
- `ops::Vector{T}`: LIOM operators
"""
function lioms(H::T, k::Int; threshold::Real=1e-14, f::Function=f)::Tuple{Vector{Float64},Matrix{Float64},Vector{T}} where {T<:AbstractOperator}
    N = qubitlength(H)
    ts = isa(H, OperatorTS)
    support = k_local_basis_1d(N, k; translational_symmetry=ts)
    return lioms(H, support; threshold=threshold, f=f)
end
