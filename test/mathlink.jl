using Test
using LinearAlgebra
using MathLink
using PauliStrings


function ising(h, N)
    H = OperatorMathLink(N)
    for i in 1:N
        H += h, "X", i
    end
    for i in 1:N
        H += "Z", i, "Z", mod1(i + 1, N)
    end
    return H
end

function ising_ts(h, N)
    H = OperatorMathLinkTS{(N,)}()
    H += h, "X", 1
    H += "Z", 1, "Z", 2
    return H
end

bn_strings = [
"Times[2, Power[2, Rational[1, 2]]]",
"Power[Plus[8, Times[8, Power[h, 2]]], Rational[1, 2]]",
"Times[2, h, Power[Times[Power[Plus[1, Power[h, 2]], -1], Plus[5, Times[2, Power[h, 2]]]], Rational[1, 2]]]",
"Power[Times[Power[Plus[5, Times[7, Power[h, 2]], Times[2, Power[h, 4]]], -1], Plus[64, Times[96, Power[h, 2]], Times[68, Power[h, 4]]]], Rational[1, 2]]",
"Times[2, Power[Times[Power[Plus[80, Times[152, Power[h, 2]], Times[133, Power[h, 4]], Times[34, Power[h, 6]]], -1], Plus[64, Times[516, Power[h, 2]], Times[587, Power[h, 4]], Times[231, Power[h, 6]], Times[96, Power[h, 8]]]], Rational[1, 2]]]"
]


@testset "mathlink operator" begin
    A = OperatorMathLink(1) + (1, "X", 1)
    @test string(norm(A)) == "Power[2, Rational[1, 2]]"
    B = OperatorMathLink(1) + (1, "Y", 1)
    @test string(norm(A*B)) == "Power[2, Rational[1, 2]]"
end


@testset "lanczos with MathLink" begin
    N = 10
    O = OperatorMathLink(N) + (1, "X", 1)
    H = ising(W`h`, N)
    assumptions = W`Assumptions -> h > 0`
    bn = lanczos(H, O, 5; assumptions=assumptions)
    for (b, bn_str) in zip(bn, bn_strings)
        @test string(b) == bn_str
    end

end

@testset "OperatorMathLinkTS construction" begin
    N = 10

    O = Operator(N) + (1, "X", 1)
    Ots = OperatorMathLinkTS{(N,)}(O)
    Oml = OperatorMathLink(N) + (1, "X", 1)
    @test Ots isa Operator{<:PauliStringTS}
    @test length(Ots) == 1
    @test typeof(Ots.coeffs[1]) == typeof(Oml.coeffs[1])

    Hts = ising_ts(W`h`, N)
    @test Hts isa Operator{<:PauliStringTS}
    @test length(Hts) == 2
end

@testset "lanczos with OperatorMathLinkTS" begin
    N = 10
    O = OperatorMathLinkTS{(N,)}()
    O += 1, "X", 1
    H = ising_ts(W`h`, N)
    assumptions = W`Assumptions -> h > 0`
    bn = lanczos(H, O, 5; assumptions=assumptions)
    for (b, bn_str) in zip(bn, bn_strings)
        @test string(b) == bn_str
    end
end
