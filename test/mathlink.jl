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
    N = 6

    H = OperatorMathLinkTS{(N,)}(N)
    @test eltype(H.strings) <: PauliStringTS
    @test eltype(H.coeffs) != ComplexF64
    @test length(H) == 0

    @test_throws ErrorException OperatorMathLinkTS{(N,)}(N + 1)

    H += W`h`, "X", 1
    @test length(H) == 1
    @test length(resum(H)) == N

    H += "Z", 1, "Z", 2
    @test length(H) == 2
    @test length(resum(H)) == 2N

    H_dense = ising(W`h`, N)
    H_ts = OperatorMathLinkTS{(N,)}(H_dense)
    @test eltype(H_ts.strings) <: PauliStringTS
    @test eltype(H_ts.coeffs) != ComplexF64
    @test any(occursin("h", string(c)) for c in H_ts.coeffs)
    @test length(H_ts) <= length(H_dense)
end


@testset "OperatorMathLinkTS full=true uses exact rational" begin
    N = 6
    H_dense = OperatorMathLink(N) + (W`h`, "X", 1)
    H_ts = OperatorMathLinkTS{(N,)}(H_dense; full=true)
    @test string(H_ts.coeffs[1]) == "Times[Rational[1, 6], h]"
end


@testset "OperatorMathLinkTS: norm and simplify" begin
    N = 4
    O = OperatorMathLinkTS{(N,)}(N)
    O += 1, "X", 1
    @test string(norm(O)) == "4"

    O_simp = simplify_operator(O)
    @test eltype(O_simp.strings) <: PauliStringTS
    @test length(O_simp) == length(O)
end


@testset "OperatorMathLinkTS: simplify_operator uses assumptions" begin
    N = 4
    O = OperatorMathLinkTS{(N,)}(N)
    O += W`Sqrt[h^2]`, "X", 1

    without = simplify_operator(O)
    @test string(without.coeffs[1]) == "Power[Power[h, 2], Rational[1, 2]]"

    with = simplify_operator(O; assumptions=W`Assumptions -> h > 0`)
    @test string(with.coeffs[1]) == "h"
    @test eltype(with.strings) <: PauliStringTS
end


@testset "OperatorMathLinkTS: lanczos reproduces OperatorMathLink result" begin
    N = 10
    assumptions = W`Assumptions -> h > 0`

    O_dense = OperatorMathLink(N) + (1, "X", 1)
    H_dense = ising(W`h`, N)
    bn_dense = lanczos(H_dense, O_dense, 5; assumptions=assumptions)

    O_ts = OperatorMathLink(N) + (1, "X", 1)
    H_ts = OperatorMathLinkTS{(N,)}(ising(W`h`, N); full=true)
    bn_ts = lanczos(H_ts, O_ts, 5; assumptions=assumptions)

    @test length(bn_ts) == length(bn_dense)
    for (b_dense, b_ts) in zip(bn_dense, bn_ts)
        @test string(b_dense) == string(b_ts)
    end

    for (b, bn_str) in zip(bn_ts, bn_strings)
        @test string(b) == bn_str
    end
end


@testset "OperatorMathLinkTS: issue acceptance workflow" begin
    function xtot(N)
        O = OperatorMathLink(N)
        for i in 1:N
            O += 1, "X", i
        end
        return O
    end

    N = 10
    assumptions = W`Assumptions -> h > 0`

    O = xtot(N)
    O += 1, "X", 1

    bn_ref = lanczos(ising(W`h`, N), O, 5; assumptions=assumptions)

    H_ts = OperatorMathLinkTS{(N,)}(ising(W`h`, N); full=true)
    @test eltype(H_ts.strings) <: PauliStringTS
    bn_ts = lanczos(H_ts, O, 5; assumptions=assumptions)

    @test length(bn_ts) == length(bn_ref)
    for (b_ref, b_ts) in zip(bn_ref, bn_ts)
        @test string(b_ref) == string(b_ts)
    end
end
