using PauliStrings
import PauliStrings as ps
using Test
using LinearAlgebra: norm

include("examples.jl")


@testset "io" begin
    @test construction_example1().w == [0, 1, 2]
    @test construction_example2().v == [6, 0, 132]
    @test construction_example2().w == [2, 3, 0]
    @test construction_example3().v == [0, 6]
    @test construction_example3().w == [12, 3]
    @test get_coef(construction_example1(), 4, 0) == 2
    @test get_coef(construction_example1(), 0, 1) == 1
    @test get_coef(construction_example1(), 2, 2) == 1
end

@testset "random" begin
    N = 10
    @test length(ps.rand_local1(N)) == 3 * N
    @test ps.rand_local1(N).N == N
    @test length(ps.rand_local2(N)) == 9 * N * (N - 1) / 2
    @test ps.rand_local2(N).N == N
end

@testset "operators" begin
    @test Operator128 <: Operator
    @test Operator64 <: Operator
    @test OperatorTS1D128 <: OperatorTS1D
    @test OperatorTS1D64 <: OperatorTS1D
    @test OperatorTS1D <: Operator
    @test typeof(Operator(10)) <: Operator
    @test typeof(Operator(70)) <: Operator
    ising10 = ising1D(10, 1)
    ising70 = ising1D(70, 1)
    @test typeof(ising10) <: Operator
    @test typeof(ising70) <: Operator
    @test typeof(ising10) == Operator64
    @test typeof(ising70) == Operator128
    ising10ts = OperatorTS1D(ising10)
    ising70ts = OperatorTS1D(ising70)
    @test typeof(ising10ts) == OperatorTS1D64
    @test typeof(ising70ts) == OperatorTS1D128
    @test typeof(Operator(ising10ts)) == Operator64
    @test typeof(Operator(ising70ts)) == Operator128
end

@testset "operations" begin
    for N in (10, 10)
        O1 = rand_local1_M(N, 20)
        O2 = rand_local2_M(N, 20)
        @test trace(O1 * O2) == 0
        @test trace(O1 * O1) != 0
        @test opnorm(O1) > 0
        @test opnorm(dagger(XX(N)) - XX(N)) == 0
        @test opnorm(dagger(X(N)) - X(N)) == 0
        @test opnorm(dagger(O1) - O1) == 0
        @test opnorm(dagger(O2) - O2) == 0
        @test opnorm(com(O1, O2) - (O1 * O2 - O2 * O1)) <= 1e-10
        @test opnorm(com(O1, O2, anti=true) - (O1 * O2 + O2 * O1)) <= 1e-10
        @test opnorm(com(O1, eye(N))) == 0
        @test opnorm(com(XX(N), Y(N))) == 0
        @test opnorm(com(XX(N), X(N))) == opnorm(com(XX(N), Z(N)))
        @test trace(diag(O1)) == trace(O1)
        @test trace(diag(O2)) == trace(O2)
    end
    O = Operator(6)
    O += "XYZZ1Y"
    @test xcount(O.v[1], O.w[1]) == 1
    @test ycount(O.v[1], O.w[1]) == 2
    @test zcount(O.v[1], O.w[1]) == 2
end

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


@testset "equivalence" begin
    @test length(equivalence_class(Operator("Y11111"), XX(6))) == 72
    @test length(equivalence_class(Operator("Y111Z1"), XX(6))) == 512
end
