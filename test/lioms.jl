# tests for lioms functionality


@testset "k_local_basis without translational symmetry" begin
    k = 2
    N = 2 * k + 1

    basis = k_local_basis(N, k; translational_symmetry=false)

    expected_count = 3 * 4^(k - 1) * N
    @test length(basis) == expected_count

    for op in basis
        @test qubitlength(op) == N
    end

end

@testset "k_local_basis with translational symmetry" begin
    k = 2
    N = 2 * k + 1

    basis = k_local_basis(N, k; translational_symmetry=true)

    expected_count = 3 * 4^(k - 1)
    @test length(basis) == expected_count

    # All should be OperatorTS types
    for op in basis
        @test qubitlength(op) == N
    end

end

@testset "k_local_basis for different k values" begin

    k1 = 1
    N = 2 * k1 + 1
    basis_k1 = k_local_basis(N, k1; translational_symmetry=false)
    @test length(basis_k1) == 3 * N

    basis_k1_ts = k_local_basis(N, k1; translational_symmetry=true)
    @test length(basis_k1_ts) == 3

    k2 = 4
    N = 2 * k2 + 1
    basis_k2 = k_local_basis(N, k2; translational_symmetry=false)
    @test length(basis_k2) == 3 * 4^(k2 - 1) * N

    basis_k2_ts = k_local_basis(N, k2; translational_symmetry=true)
    @test length(basis_k2_ts) == 3 * 4^(k2 - 1)
end

@testset "k_local_basis orthogonality" begin
    k = 2
    N = 2 * k + 1
    basis = k_local_basis(N, k; translational_symmetry=false)

    for i in eachindex(basis)
        for j in i+1:length(basis)
            op_i = basis[i]
            op_j = basis[j]
            tp = trace_product(op_i, op_j; scale=1)
            @test isapprox(tp, 0.0; atol=1e-10)
        end
    end

    basis_ts = k_local_basis(N, k; translational_symmetry=true)
    for i in eachindex(basis_ts)
        for j in i+1:length(basis_ts)
            op_i = basis_ts[i]
            op_j = basis_ts[j]
            tp = trace_product(op_i, op_j; scale=1)
            @test isapprox(tp, 0.0; atol=1e-10)
        end
    end
end


@testset "lioms translational symmetry on/off" begin
    k = 3
    N = 2 * k + 1

    H = XXZ(N, 1.0, 0.5; ts=false)
    support = k_local_basis(N, k; translational_symmetry=false)
    evals, ops = lioms(H, support; return_all=false)
    @test length(evals) == k
    @test length(ops) == k

    H_ts = XXZ(N, 1.0, 0.5; ts=true)
    support_ts = k_local_basis(N, k; translational_symmetry=true)
    evals_ts, ops_ts = lioms(H_ts, support_ts; return_all=false)
    @test length(evals_ts) == k
    @test length(ops_ts) == k
end

@testset "lioms return_all option" begin
    k = 3
    N = 2 * k + 1
    H = XXZ(N, 1.0, 0.5)
    support = k_local_basis(N, k; translational_symmetry=false)

    evals_partial, ops_partial = lioms(H, support; return_all=false)
    evals_all, ops_all = lioms(H, support; return_all=true)

    @test length(evals_all) >= length(evals_partial)
    @test length(ops_all) >= length(ops_partial)

    for i in eachindex(evals_partial)
        @test isapprox(evals_all[i], evals_partial[i]; atol=1e-10)
    end
end

@testset "lioms from support size" begin
    k = 3
    N = 2 * k + 1
    H = XXZ(N, 1.0, 0.5; ts=false)

    # Generate support and run lioms with explicit support
    support = k_local_basis(N, k; translational_symmetry=false)
    evals_explicit, ops_explicit = lioms(H, support; return_all=false)

    evals_k, ops_k = lioms(H, k; return_all=false)

    # Should get the same results
    @test length(evals_explicit) == length(evals_k)
    @test length(ops_explicit) == length(ops_k)

    for i in eachindex(evals_explicit)
        @test isapprox(evals_explicit[i], evals_k[i]; atol=1e-10)
    end
end

@testset "lioms properties" begin
    k = 3
    N = 2 * k + 1
    H = XXZ(N, 1.0, 0.5; ts=true)
    support = ps.k_local_basis(N, k; translational_symmetry=true)
    evals, ops = ps.lioms(H, support; return_all=false)

    @test all(abs.(evals) .< 1e-12)
    @test all((commutator(H, O) |> opnorm) .< 1e-12 for O in ops)

    for i in eachindex(ops)
        for j in i+1:length(ops)
            @test abs(trace_product(ops[i], ops[j]; scale=1)) <= 1e-12
        end
    end
end
