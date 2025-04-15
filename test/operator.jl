

# tests for the Operator type and operations on operators


function randstring(N)
    return join(rand(["1", "X", "Y", "Z"], N))
end


@testset "io" begin
    o1 = construction_example1()
    getproperty.(o1.strings, :w) == [0, 1, 2]
    o2 = construction_example2()
    getproperty.(o2.strings, :v) == [6, 0, 132]
    getproperty.(o2.strings, :w) == [2, 3, 0]
    o3 = construction_example3()
    getproperty.(o3.strings, :v) == [0, 6]
    getproperty.(o3.strings, :w) == [12, 3]
    @test get_coef(o1, eltype(o1.strings)(4, 0)) == 2
    @test get_coef(o1, eltype(o1.strings)(0, 1)) == 1
    @test get_coef(o1, eltype(o1.strings)(2, 2)) == 1

    for N in (10, 70)
        for i in 1:10
            O = Operator(N)
            st = randstring(N)
            O += st
            @test string(only(O.strings)) == st
            # @test vw_to_string(O.v[1], O.w[1], N)[1] == st
        end
    end
end

@testset "random" begin
    N = 10
    @test length(ps.rand_local1(N)) == 3 * N
    @test qubitlength(ps.rand_local1(N)) == N
    @test length(ps.rand_local2(N)) == 9 * N * (N - 1) / 2
    @test qubitlength(ps.rand_local2(N)) == N
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
