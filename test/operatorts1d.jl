
@testset "operatorts1d" begin
    N = 6
    O1 = ising1D(N, 0.5)
    O2 = ising1D(N, 1.1)
    O1ts = OperatorTS1D(O1)
    O2ts = OperatorTS1D(O2)
    eps = 1e-10
    @test typeof(O1ts + O2ts) == typeof(O1ts)
    @test typeof(O1ts - O2ts) == typeof(O1ts)
    @test typeof(O1ts * O2ts) == typeof(O1ts)
    @test typeof(O1ts - 1) == typeof(O1ts)
    @test typeof(O1ts + 1) == typeof(O1ts)
    @test typeof(O1ts * 2) == typeof(O1ts)
    @test typeof(O1ts / 2) == typeof(O1ts)
    @test typeof(OperatorTS1D(ising1D(70, 1))) == OperatorTS1D128
    @test typeof(OperatorTS1D(ising1D(40, 1))) == OperatorTS1D64
    @test opnorm(Operator(O1ts) - O1) < eps
    @test opnorm(Operator(O2ts) - O2) < eps
    @test opnorm(Operator(O2ts) - O2) < eps
    @test opnorm(Operator(O1ts * O2ts) - O1 * O2) < eps
    @test opnorm(Operator(O1ts * O2ts * O2ts) - O1 * O2 * O2) < eps
    @test opnorm(Operator(O1ts + O2ts) - (O1 + O2)) < eps
    @test opnorm(Operator(O1ts * 7) - O1 * 7) < eps
    @test opnorm(Operator(O1ts + 7) - (O1 + 7)) < eps
    @test is_ts(O1)
    @test !is_ts(ps.rand_local2(N))
    @test trace(O1ts) == trace(dagger(O1))
    @test trace(O1ts) == trace(diag(O1ts))
end
