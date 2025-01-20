@testset "lanczos" begin
    bs1 = ps.lanczos(XX(10), X(10), 10, 2^14)
    bs2 = ps.lanczos(OperatorTS1D(XX(10)), OperatorTS1D(X(10)), 10, 2^14)
    @test norm(bs1 .- [2.8284271247461903, 4.0, 4.898979485566357, 5.65685424949238, 6.324555320336759, 6.92820323027551, 7.483314773547883, 8.0, 8.487501570645295, 8.957117508915214]) < 1e-14
    @test norm(bs1 .- bs2) < 1e-10
end

@testset "time_evolution" begin
    @test norm(evolve_chaotic() - 0.09208935978092929) < 1e-10
end

@testset "moments" begin
    for N in (10, 70)
        O1 = rand_local2_M(N, 15)
        O2 = rand_local2_M(N, 15)
        @test opnorm(oppow(O1, 3) - O1 * O1 * O1) < 1e-9
        @test abs(trace_product(O1, O2) - trace(O1 * O2)) < 1e-9
        @test abs(trace_product(O1, 4) - trace(oppow(O1, 4))) < 1e-9
        @test abs(trace_product(O1, 4, O2, 3) - trace(oppow(O1, 4) * oppow(O2, 3))) < 1e-8
    end
    N = 6
    O1 = ising1D(N, 0.5)
    O2 = ising1D(N, 1.1)
    O1ts = OperatorTS1D(O1)
    O2ts = OperatorTS1D(O2)
    @test opnorm(oppow(O1ts, 3) - O1ts * O1ts * O1ts) < 1e-10
    @test abs(trace_product(O1ts, 4) - trace_product(O1, 4)) < 1e-10
    @test abs(trace_product(O1ts, 4, O2ts, 3) - trace_product(O1, 4, O2, 3)) < 1e-10
end


@testset "equivalence" begin
    @test length(equivalence_class(Operator("Y11111"), XX(6))) == 72
    @test length(equivalence_class(Operator("Y111Z1"), XX(6))) == 512
end
