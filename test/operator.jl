using Base: normalize_depots_for_relocation
using PauliStrings: paulistringtype

# tests for the Operator type and operations on operators


@testset "random" begin
    N = 10
    @test length(ps.rand_local1(N)) == 3 * N
    @test qubitlength(ps.rand_local1(N)) == N
    @test length(ps.rand_local2(N)) == 9 * N * (N - 1) / 2
    @test qubitlength(ps.rand_local2(N)) == N
end

@testset "operators" begin
    @test Operator <: AbstractOperator
    @test OperatorTS1D <: AbstractOperator
    @test typeof(Operator(10)) <: Operator
    @test typeof(Operator(70)) <: Operator
    ising10 = ising1D(10, 1)
    ising70 = ising1D(70, 1)
    @test typeof(ising10) <: Operator
    @test typeof(ising70) <: Operator
    @test typeof(ising10) == Operator{paulistringtype(10),ComplexF64}
    @test typeof(ising70) == Operator{paulistringtype(70),ComplexF64}
    ising10ts = OperatorTS1D(ising10)
    ising70ts = OperatorTS1D(ising70)
    @test typeof(ising10ts) == OperatorTS1D{paulistringtype(10),ComplexF64}
    @test typeof(ising70ts) == OperatorTS1D{paulistringtype(70),ComplexF64}
    @test typeof(Operator(ising10ts)) == typeof(ising10)
    @test typeof(Operator(ising70ts)) == typeof(ising70)
end

@testset "operations" begin
    for N in (10, 70)
        O1 = rand_local1_M(N, 20)
        O2 = rand_local2_M(N, 20)
        @test trace(O1 * O2) == 0
        @test trace(O1 * O1) != 0
        @test opnorm(O1) > 0
        @test opnorm(dagger(XX(N)) - XX(N)) == 0
        @test opnorm(dagger(X(N)) - X(N)) == 0
        @test opnorm(dagger(O1) - O1) == 0
        @test opnorm(dagger(O2) - O2) == 0
        @test opnorm(commutator(O1, O2) - (O1 * O2 - O2 * O1)) <= 1e-10
        @test opnorm(anticommutator(O1, O2) - (O1 * O2 + O2 * O1)) <= 1e-10
        @test opnorm(commutator(O1, eye(N))) == 0
        @test opnorm(commutator(XX(N), Y(N))) == 0
        @test opnorm(commutator(XX(N), X(N))) == opnorm(commutator(XX(N), Z(N)))
        @test trace(diag(O1)) == trace(O1)
        @test trace(diag(O2)) == trace(O2)
    end
    O = Operator(6)
    O += "XYZZ1Y"
    @test xcount(O.strings[1]) == 1
    @test ycount(O.strings[1]) == 2
    @test zcount(O.strings[1]) == 2


    o = construction_example2()
    o2 = ptrace(o, [3, 4])
    @test trace(o2; normalize=true) == 128
    o3 = typeof(o2)()
    o3 += "11Z1111Z"
    o3 += 128, "11111111"
    o3 += 1.5, "1YZ11111"
    @test opnorm(o2 - o3) == 0
end
