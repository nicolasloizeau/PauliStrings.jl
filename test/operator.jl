using PauliStrings: paulistringtype

# tests for the Operator type and operations on operators
Ns = (10, 70, 500)

@testset "random" begin
    N = 10
    @test length(ps.rand_local1(N)) == 3 * N
    @test qubitlength(ps.rand_local1(N)) == N
    @test length(ps.rand_local2(N)) == 9 * N * (N - 1) / 2
    @test qubitlength(ps.rand_local2(N)) == N
end

@testset "operators" begin
    @test Operator <: AbstractOperator
    for N in Ns
        O = Operator(N)
        @test typeof(O) <: Operator
        @test typeof(O) == Operator{paulistringtype(N),ComplexF64}
        @test qubitlength(O) == N
        ising = ising1D(N, 1)
        @test typeof(ising) <: Operator
        @test typeof(ising) == Operator{paulistringtype(N),ComplexF64}
        isingts = OperatorTS1D(ising)
        @test typeof(Operator(isingts)) == typeof(ising)
    end
end

@testset "operations" begin
    for N in Ns
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
