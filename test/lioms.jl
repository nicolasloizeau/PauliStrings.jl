# tests for lioms functionality

@testset "symmetry_adapted_k_local_basis basic functionality" begin
    k = 3
    N = 2 * k + 1

    # Test with default parameters
    basis = ps.symmetry_adapted_k_local_basis(N, k)
    @test length(basis) > 0
    @test all(op -> qubitlength(op) == N, basis)

    # All operators should be OperatorTS types by default
    @test all(op -> isa(op, OperatorTS), basis)
end


@testset "symmetry_adapted_k_local_basis time_reversal symmetry" begin
    k = 3
    N = 2 * k + 1

    # Test with real time reversal (O + O†)
    basis_real = ps.symmetry_adapted_k_local_basis(N, k; time_reversal=:real)
    @test length(basis_real) > 0

    # Test with imaginary time reversal (i(O - O†))
    basis_imag = ps.symmetry_adapted_k_local_basis(N, k; time_reversal=:imag)
    @test length(basis_imag) > 0

end


@testset "symmetry_adapted_k_local_basis spin_flip symmetry" begin
    k = 3
    N = 2 * k + 1

    # Test with even spin flip
    basis_even = ps.symmetry_adapted_k_local_basis(N, k; spin_flip=:even)
    @test length(basis_even) > 0

    # Test with odd spin flip
    basis_odd = ps.symmetry_adapted_k_local_basis(N, k; spin_flip=:odd)
    @test length(basis_odd) > 0

    # Test with both
    basis_both = ps.symmetry_adapted_k_local_basis(N, k; spin_flip=:both)
    @test length(basis_both) > 0

    # Both should contain even + odd
    @test length(basis_both) >= length(basis_even)
    @test length(basis_both) >= length(basis_odd)
end


@testset "symmetry_adapted_k_local_basis magnetization conservation" begin
    k = 3
    N = 2 * k + 1

    # Test with magnetization conservation
    basis_yes = ps.symmetry_adapted_k_local_basis(N, k; conserve_magnetization=:yes)
    @test length(basis_yes) > 0

    # Test without magnetization conservation
    basis_no = ps.symmetry_adapted_k_local_basis(N, k; conserve_magnetization=:no)
    @test length(basis_no) > 0

    # Test with both
    basis_both = ps.symmetry_adapted_k_local_basis(N, k; conserve_magnetization=:both)
    @test length(basis_both) > 0

    # Both should be larger
    @test length(basis_both) >= length(basis_yes)
    @test length(basis_both) >= length(basis_no)
end


@testset "symmetry_adapted_k_local_basis translational symmetry" begin
    k = 3
    N = 2 * k + 1

    # With translational symmetry (default)
    basis_ts = ps.symmetry_adapted_k_local_basis(N, k; translational_symmetry=true)
    @test all(op -> isa(op, OperatorTS), basis_ts)

    # Without translational symmetry
    basis_no_ts = ps.symmetry_adapted_k_local_basis(N, k; translational_symmetry=false)
    @test all(op -> isa(op, Operator), basis_no_ts)

    # Without TS should have more operators (all translations)
    @test length(basis_no_ts) > length(basis_ts)
end


@testset "symmetry_adapted_k_local_basis error handling" begin
    k = 3
    N = 2 * k + 1

    # Invalid time_reversal
    @test_throws ErrorException ps.symmetry_adapted_k_local_basis(N, k; time_reversal=:invalid)

    # Invalid spin_flip
    @test_throws ErrorException ps.symmetry_adapted_k_local_basis(N, k; spin_flip=:invalid)

    # Invalid conserve_magnetization
    @test_throws ErrorException ps.symmetry_adapted_k_local_basis(N, k; conserve_magnetization=:invalid)
end


@testset "symmetry_adapted_k_local_basis orthogonality" begin
    k = 4
    N = 2 * k + 1

    basis = ps.symmetry_adapted_k_local_basis(N, k; translational_symmetry=true)

    # Check orthogonality of basis elements
    for i in eachindex(basis)
        for j in i+1:length(basis)
            op_i = basis[i]
            op_j = basis[j]
            tp = trace_product(op_i, op_j; scale=1)
            @test isapprox(tp, 0.0; atol=1e-10)
        end
    end
end


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

@testset "lioms in symmetry adapted basis" begin
    M = 3
    L = 2 * M + 1

    H = XXZ(L, 1.0, 0.5; ts=true)

    # imaginary sector
    support = ps.symmetry_adapted_k_local_basis(L, M; time_reversal=:imag, spin_flip=:even, conserve_magnetization=:yes, translational_symmetry=true)
    evals, ops = ps.lioms(H, support; return_all=false)
    # only energy current expected
    @test length(evals) == 1
    for i in eachindex(ops)
        @test (commutator(H, ops[i]) |> opnorm) < 1e-10
    end

    # real sector
    support_real = ps.symmetry_adapted_k_local_basis(L, M; time_reversal=:real, spin_flip=:even, conserve_magnetization=:yes, translational_symmetry=true)
    evals_real, ops_real = ps.lioms(H, support_real; return_all=false)
    # only hamiltonian expected
    @test length(evals_real) == 1
    for i in eachindex(ops_real)
        @test (commutator(H, ops_real[i]) |> opnorm) < 1e-10
    end
end