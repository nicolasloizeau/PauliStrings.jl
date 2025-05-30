@testset "Addition" begin
    N = 2
    a, b = @variables a b
    o = Operator{paulistringtype(N), Complex{Num}}()
    o = o + (a, "ZX")
    o = o + (b, "XY")
    @test o.strings[1] == paulistringtype(N)("ZX")
    @test isequal(o.coeffs[1], Complex{Num}(a, 0))
    @test o.strings[2] == paulistringtype(N)("XY")
    @test isequal(o.coeffs[2], Complex{Num}(0, b))
end
@testset "Simplify operator" begin
    N = 1
    a, b = @variables a b
    o = Operator{paulistringtype(N), Complex{Num}}()
    o = o + (a, "X")
    o = o + (b, "X")
    os = simplify(o)
    @test length(os.strings) == 1
    @test string(os.strings[1]) == "X"
    @test isequal(os.coeffs[1], Complex{Num}(a + b, 0))
    H = Operator{paulistringtype(1), Complex{Num}}()
    H = H + (a, "X")
    C = commutator(H, H)
    Cs = simplify(C)
    @test isempty(Cs.strings)
end
@testset "Commutator and anticommutator" begin
    N = 1
    a, b = @variables a b
    X = Operator{paulistringtype(N), Complex{Num}}() + (a, "X")
    Y = Operator{paulistringtype(N), Complex{Num}}() + (b, "Y")
    Z = Operator{paulistringtype(N), Complex{Num}}() + (1, "Z")
    @testset "Commutator [X, Y]" begin
        C = commutator(X, Y)
        Cs = simplify(C)
        @test length(Cs.strings) == 1
        @test Cs.strings[1] == paulistringtype(N)("Z")
        @test isequal(Cs.coeffs[1], Complex{Num}(0, 2a*b))
    end
    @testset "Anticommutator {X, Z}" begin
        A = anticommutator(X, Z)
        As = simplify(A)
        @test isempty(As.strings)
    end
    @testset "Nested operations" begin
        N = 1
        a, b = @variables a b
        X = Operator{paulistringtype(N), Complex{Num}}() + (a, "X")
        Y = Operator{paulistringtype(N), Complex{Num}}() + (b, "Y")
        Z = Operator{paulistringtype(N), Complex{Num}}() + (1, "Z")
        expr = commutator(commutator(X, Y), Z)
        exprs = simplify(expr)
        @test isempty(exprs.strings)
    end
end
