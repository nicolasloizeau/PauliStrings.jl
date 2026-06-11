@testset "lanczos" begin
    bs1 = ps.lanczos(XX(10), X(10), 10, 2^14)
    bs2 = ps.lanczos(OperatorTS1D(XX(10)), OperatorTS1D(X(10)), 10, 2^14)
    @test norm(bs1 .- [2.8284271247461903, 4.0, 4.898979485566357, 5.65685424949238, 6.324555320336759, 6.92820323027551, 7.483314773547883, 8.0, 8.487501570645295, 8.957117508915214]) < 1e-13
    @test norm(bs1 .- bs2) < 1e-10
end


@testset "moments" begin
    for N in (10, 70)
        O1 = rand_local2_M(N, 15)
        O2 = rand_local2_M(N, 15)
        @test norm(O1^3 - O1 * O1 * O1) / norm(O1^3) < 1e-9
        @test abs(trace_product(O1, O2) - trace(O1 * O2)) < 1e-9
        @test abs(trace_product(O1, 4) - trace(O1^4)) / abs(trace_product(O1, 4)) < 1e-9
        @test abs(trace_product(O1, 4, O2, 3) - trace(O1^4 * O2^3)) < 1e-8
    end
    N = 6
    O1 = ising1D(N, 0.5)
    O2 = ising1D(N, 1.1)
    O1ts = OperatorTS1D(O1)
    O2ts = OperatorTS1D(O2)
    @test norm(O1ts^3 - O1ts * O1ts * O1ts) < 1e-10
    @test abs(trace_product(O1ts, 4) - trace_product(O1, 4)) < 1e-10
    @test abs(trace_product(O1ts, 4, O2ts, 3) - trace_product(O1, 4, O2, 3)) < 1e-10
end


@testset "trace_moment" begin
    # trace_moment(H, k) must agree with trace_product(H, k) = trace(H^k)
    # k = 0 is the identity trace
    for N in (4, 6, 8)
        @test trace_moment(Operator(N), 0; scale=1) ≈ 1.0
        @test trace_moment(Operator(N), 5; scale=1) ≈ 0.0
    end

    # random real 2-local operators (Hermitian). Both methods accumulate in Float64, so we
    # compare with a relative tolerance.
    for N in (4, 6, 8, 10)
        O = rand_local2_M(N, 12)
        for k in 0:8
            a = trace_moment(O, k; scale=1)
            b = trace_product(O, k; scale=1)
            @test abs(a - b) ≤ 1e-7 * (abs(b) + 1e-9) + 1e-12
        end
    end

    # 1-local operators and the transverse-field Ising model
    for N in (4, 6, 8)
        O = rand_local1_M(N, 8)
        for k in 0:9
            @test abs(trace_moment(O, k; scale=1) - trace_product(O, k; scale=1)) ≤
                  1e-7 * (abs(trace_product(O, k; scale=1)) + 1e-9) + 1e-12
        end
    end
    H = ising1D(8, 0.5)
    for k in 0:10
        a = trace_moment(H, k; scale=1)
        b = trace_product(H, k; scale=1)
        @test abs(a - b) ≤ 1e-7 * (abs(b) + 1e-9) + 1e-12
    end

    # complex coefficients (non-Hermitian) must also match
    O = Operator(5)
    O += (1.0 + 0.3im, "X", 1)
    O += (0.7, "Z", 2, "Z", 3)
    O += (-0.4im, "Y", 4)
    O += (0.9, "Z", 1, "X", 5)
    for k in 0:7
        @test abs(trace_moment(O, k) - trace_product(O, k)) ≤ 1e-8 * (abs(trace_product(O, k)) + 1e-9) + 1e-9
    end

    # the scale keyword behaves like trace_product
    O = rand_local2_M(6, 10)
    @test trace_moment(O, 4) ≈ trace_product(O, 4)
    @test trace_moment(O, 4; scale=2.0^6) ≈ trace_moment(O, 4)
    @test trace_moment(O, 0) ≈ 2.0^6

    @test_throws ArgumentError trace_moment(O, -1)
end


@testset "trace_moment A^k B^l" begin
    # trace_moment(A, k, B, l) must agree with trace_product(A, k, B, l) = trace(A^k B^l)
    for N in (4, 6, 8)
        A = rand_local2_M(N, 10)
        B = rand_local2_M(N, 6)
        for k in 1:4, l in 1:3
            a = trace_moment(A, k, B, l; scale=1)
            b = trace_product(A, k, B, l; scale=1)
            @test abs(a - b) ≤ 1e-7 * (abs(b) + 1e-9) + 1e-12
        end
    end

    # response-style trace(H^k O) with O a single Pauli (the issue #80 follow-up example)
    H = ising1D(8, 0.5)
    O = Operator(8) + ("X", 1)
    for k in 1:8
        a = trace_moment(H, k, O, 1; scale=1)
        b = trace_product(H, k, O, 1; scale=1)
        @test abs(a - b) ≤ 1e-7 * (abs(b) + 1e-9) + 1e-12
    end

    # k=0 / l=0 reduce to plain moments
    A = rand_local2_M(6, 8)
    B = rand_local2_M(6, 5)
    @test trace_moment(A, 0, B, 3; scale=1) ≈ trace_moment(B, 3; scale=1)
    @test trace_moment(A, 3, B, 0; scale=1) ≈ trace_moment(A, 3; scale=1)
    @test_throws ArgumentError trace_moment(A, -1, B, 2)
end


@testset "equivalence" begin
    @test length(equivalence_class(Operator("Y11111"), XX(6))) == 72
    @test length(equivalence_class(Operator("Y111Z1"), XX(6))) == 512
end

@testset "graph" begin
    O1 = ising1D(8, 0.5)
    G = frustration_graph(O1)
    @test sum(G) == 32
end
