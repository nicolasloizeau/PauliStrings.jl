using Test
using LinearAlgebra
using MathLink
using PauliStrings


function ising(N, h)
    H = OperatorMathLink(N)
    for i in 1:N
        H += h, "X", i
    end
    for i in 1:N
        H += "Z", i, "Z", mod1(i + 1, N)
    end
    return H
end

function ising_ts(N, h)
    H = OperatorMathLinkTS{(N,)}()
    H += h, "X", 1
    H += "Z", 1, "Z", 2
    return H
end


function xtot(N)
    O = OperatorMathLink(N)
    for i in 1:N
        O += 1, "X", i
    end
    return O
end


function xtot_ts(N)
    O = OperatorMathLinkTS{(N,)}()
    O += 1, "X", 1
    return O
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
    @test string(norm(A * B)) == "Power[2, Rational[1, 2]]"
end


@testset "lanczos with MathLink" begin
    N = 10
    O = OperatorMathLink(N) + (1, "X", 1)
    H = ising(N, W`h`)
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

    Hts = ising_ts(N, W`h`)
    @test Hts isa Operator{<:PauliStringTS}
    @test length(Hts) == 2
end

@testset "lanczos with OperatorMathLinkTS" begin
    N = 12
    O = xtot(N)
    Ots = xtot_ts(N)
    H = ising(N, W`h`)
    H_ts = ising_ts(N, W`h`)
    assumptions = W`Assumptions -> h > 0`
    bn = lanczos(H, O, 12; assumptions=assumptions)
    bn_ts = lanczos(H_ts, Ots, 12; assumptions=assumptions)
    for (b, b_ts) in zip(bn, bn_ts)
        b = weval(W"Simplify"(b.expression, assumptions))
        b = weval(W"Apart"(b))
        b_ts = weval(W"Simplify"(b_ts.expression, assumptions))
        b_ts = weval(W"Apart"(b_ts))
        @test string(weval(W`TrueQ[Simplify[$b - $b_ts == 0, $assumptions]]`)) == "True"
    end
end
