
@testset "truncation" begin
    N = 8
    O2 = ps.rand_local2(N)
    @test opnorm(ps.truncate(O2, 1)) == 0
    @test opnorm(ps.truncate(O2, 2)) == opnorm(O2)
    @test length(ps.trim(O2, 10)) == 10
    @test length(ps.prune(O2, 2)) <= length(ps.prune(O2, 20))
    @test length(ps.cutoff(XX(N) + 0.1 * X(N), 0.5)) == length(XX(N))
    @test opnorm(ps.cutoff(O2, 0.5)) <= opnorm(O2)
    @test opnorm(ps.add_noise(O2, 0.5)) < opnorm(O2)
    @test opnorm(ps.add_noise(O2, 0.5)) < opnorm(ps.add_noise(O2, 0.1))
    @test opnorm(k_local_part(O2, 1) - ps.truncate(O2, 1)) == 0
    O1 = ising1D(N, 0.5)
    O1ts = OperatorTS1D(O1)
    @test opnorm(ps.truncate(O1ts, 2)) == opnorm(O1)
    @test length(ps.trim(O1ts, 1)) == 1
    @test opnorm(ps.cutoff(O1ts, 0.8)) <= opnorm(O1ts)
    @test opnorm(ps.add_noise(O1ts, 0.5)) < opnorm(O1ts)
    @test opnorm(Operator(ps.add_noise(O1ts, 0.5)) - ps.add_noise(O1, 0.5)) < 1e-10
end
