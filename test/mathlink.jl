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

# TS version: store only the seed (local terms at site 1).
# resum() will translate these to all N sites, reproducing the full Ising H.
function ising_ts(h, N)
    H = OperatorMathLinkTS{(N,)}(N)
    H += h, "X", 1              # seed X term
    H += "Z", 1, "Z", 2        # seed ZZ bond (wraps via PauliStringTS orbit)
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
    @test length(H) == 0

    # seed: one X term → one orbit representative
    H2 = OperatorMathLinkTS{(N,)}(N)
    H2 += W`h`, "X", 1
    @test length(H2) == 1

    # resum should reproduce the full N-site operator
    H2_full = resum(H2)
    @test length(H2_full) == N

    # Convert from dense OperatorMathLink
    H_dense = ising(W`h`, N)
    H_ts = OperatorMathLinkTS{(N,)}(H_dense; full=true)
    @test eltype(H_ts.strings) <: PauliStringTS
    @test length(H_ts) <= length(H_dense)
end


@testset "OperatorMathLinkTS: norm and simplify" begin
    N = 4
    O = OperatorMathLinkTS{(N,)}(N)
    O += 1, "X", 1
    n = norm(O)
    @test n isa MathLinkNumber

    O_simp = simplify_operator(O)
    @test eltype(O_simp.strings) <: PauliStringTS
end


@testset "OperatorMathLinkTS: lanczos reproduces OperatorMathLink result" begin
    N = 10
    assumptions = W`Assumptions -> h > 0`

    # Dense reference
    O_dense = OperatorMathLink(N) + (1, "X", 1)
    H_dense = ising(W`h`, N)
    bn_dense = lanczos(H_dense, O_dense, 5; assumptions=assumptions)

    # TS: seed-based H, same plain initial operator
    O_ts = OperatorMathLink(N) + (1, "X", 1)
    H_ts = ising_ts(W`h`, N)
    bn_ts = lanczos(H_ts, O_ts, 5; assumptions=assumptions)

    for (b_dense, b_ts) in zip(bn_dense, bn_ts)
        @test string(b_dense) == string(b_ts)
    end
end