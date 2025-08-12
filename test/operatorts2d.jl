
@testset "operatorts2d" begin
    for L1 in 1:4
        L2 = 3
        O1 = ising2D(L1, L2, rand())
        O2 = ising2D(L1, L2, rand())
        O1ts = OperatorTS2D(O1, L1)
        O2ts = OperatorTS2D(O2, L1)
        eps = 1e-8
        @test (opnorm(O1) - opnorm(O1ts)) < eps
        @test (opnorm(O2) - opnorm(O2ts)) < eps
        @test typeof(O1ts + O2ts) == typeof(O1ts)
        @test typeof(O1ts - O2ts) == typeof(O1ts)
        @test typeof(O1ts * O2ts) == typeof(O1ts)
        @test typeof(O1ts - 1) == typeof(O1ts)
        @test typeof(O1ts + 1) == typeof(O1ts)
        @test typeof(O1ts * 2) == typeof(O1ts)
        @test typeof(O1ts / 2) == typeof(O1ts)
        @test typeof(OperatorTS2D(O1ts)) <: OperatorTS2D
        @test opnorm(Operator(O1ts) - O1) < eps
        @test opnorm(Operator(O1ts * O2ts) - O1 * O2) < eps
        @test opnorm(Operator(O1ts * O2ts * O2ts) - O1 * O2 * O2) < eps
        @test opnorm(Operator(O1ts + O2ts) - (O1 + O2)) < eps
        @test opnorm(Operator(O1ts * 7) - O1 * 7) < eps
        @test opnorm(Operator(O1ts + 7) - (O1 + 7)) < eps
        @test is_ts2d(O1, L1)
        @test is_ts2d(O1 * O2, L1)
        @test !is_ts2d(ps.rand_local2(L1 * L2), L1)
        @test trace(O1ts) == trace(dagger(O1))
        @test trace(O1ts) == trace(diag(O1ts))
        k = 1 + L1
        @test isapprox(trace_product(O1ts, k), trace_product(O1, k))
        @test isapprox(trace_product(O1ts, k), trace(O1^k))
        @test isapprox(trace_product(O1ts, k), trace(O1ts^k))
    end
end
