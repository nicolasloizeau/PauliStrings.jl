include("symbolics_ex.jl")

@testset "symbolics" begin
    eps = 1e-8
    N = 2
    @variables g
    dict = Dict(g => 0.55)
    O = OperatorSymbolic(N)
    @test typeof(O + 1) == typeof(O)
    @test typeof(O - 1) == typeof(O)
    @test typeof(O / 0.5) == typeof(O)
    @test typeof(O * 0.5) == typeof(O)
    @test typeof(O * 0.5) == typeof(O)
    @test typeof(O + O) == typeof(O)
    @test typeof(O * O) == typeof(O)
    O += g, "Z", 1
    O += g^2, "X", 1, "X", 2
    Os = substitute_op(O, dict)
    Op = Operator(N)
    Op += dict[g], "Z", 1
    Op += dict[g]^2, "X", 1, "X", 2
    @test typeof(Os) == typeof(Op)
    @test norm(Os - Op) < eps
    @test norm(substitute_op(simplify_op(O), dict) - Op) < eps
    @test norm(substitute_op(0.5 * O^4, dict) - 0.5 * Op^4) < eps
    @test norm(substitute_op(simplify_op(0.5 * O^4), dict) - 0.5 * Os^4) < eps
    @test abs(substitute(trace_product(O, 3), dict) - trace_product(Op, 3)) < eps
    @test abs(substitute(trace_product(O, 4), dict) - trace_product(Os, 4)) < eps
end
