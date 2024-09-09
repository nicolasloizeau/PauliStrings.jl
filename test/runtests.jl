using PauliStrings
import PauliStrings as ps
using Test
using LinearAlgebra: norm

include("examples.jl")


@testset "io" begin
    @test construction_example1().w == [0,1,2]
    @test construction_example2().v == [6,0,132]
    @test construction_example2().w == [2,3,0]
    @test construction_example3().v == [0,6]
    @test construction_example3().w == [12,3]
end

@testset "random" begin
    N = 10
    @test length(ps.rand_local1(N)) == 3*N
    @test ps.rand_local1(N).N == N
    @test length(ps.rand_local2(N)) == 9*N*(N-1)/2
    @test ps.rand_local2(N).N == N
end

@testset "operations" begin
    N = 8
    O1 = ps.rand_local1(N)
    O2 = ps.rand_local2(N)
    @test trace(O1*O2)==0
    @test trace(O1*O1)!=0
    @test opnorm(O1)>0
    @test opnorm(dagger(XX(N))-XX(N)) == 0
    @test opnorm(dagger(X(N))-X(N)) == 0
    @test opnorm(dagger(O1)-O1) == 0
    @test opnorm(dagger(O2)-O2) == 0
    @test opnorm(com(O1, eye(N))) == 0
    @test opnorm(com(XX(N), Y(N))) == 0
    @test opnorm(com(XX(N), X(N))) == opnorm(com(XX(N), Z(N)))
    @test trace(diag(O1)) == trace(O1)
    @test trace(diag(O2)) == trace(O2)
    O = Operator(6)
    O += "XYZZ1Y"
    @test xcount(O.v[1], O.w[1]) == 1
    @test ycount(O.v[1], O.w[1]) == 2
    @test zcount(O.v[1], O.w[1]) == 2
end

@testset "truncation" begin
    N = 8
    O2 = ps.rand_local2(N)
    @test opnorm(ps.truncate(O2, 1))==0
    @test opnorm(ps.truncate(O2, 2))==opnorm(ps.truncate(O2, 2))
    @test length(ps.trim(O2, 10)) == 10
    @test length(ps.prune(O2, 2)) <= length(ps.prune(O2, 20))
    @test length(ps.cutoff(XX(N)+0.1*X(N), 0.5)) == length(XX(N))
    @test opnorm(ps.cutoff(O2, 0.5)) <= opnorm(O2)
    @test opnorm(ps.add_noise(O2, 0.5)) < opnorm(O2)
    @test opnorm(ps.add_noise(O2, 0.5)) < opnorm(ps.add_noise(O2, 0.1))
    @test opnorm(k_local_part(O2, 1)-ps.truncate(O2, 1)) == 0
end

@testset "lanczos" begin
    O = Operator(10)
    O += "X", 1
    @test norm(ps.lanczos(XX(10), X(10), 10, 2^14).-[2.8284271247461903, 4.0, 4.898979485566357, 5.65685424949238, 6.324555320336759, 6.92820323027551, 7.483314773547883, 8.0, 8.487501570645295, 8.957117508915214]) <1e-14
    @test norm(ps.lanczos(XX(10), O, 10, 2^14; localop=true).-ps.lanczos(XX(10), X(10), 10, 2^14)) <1e-2
end

@testset "time_evolution" begin
    @test norm(evolve_chaotic()-0.09208935978092929) < 1e-10
end
