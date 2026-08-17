using PauliStrings
import PauliStrings as ps
using Test
using LinearAlgebra

function heisenberg_ring_plain(N)
    H = Operator(N)
    for j in 1:N
        H += "X", j, "X", mod1(j + 1, N)
        H += "Y", j, "Y", mod1(j + 1, N)
        H += "Z", j, "Z", mod1(j + 1, N)
    end
    return H
end

_re(a, b) = abs(a - b) / max(abs(a), abs(b), 1.0)

@testset "trace_moment (plain Operator, issue #80)" begin
    eps = 1e-7

    @testset "trace_moment(o,k) == trace_product == dense" begin
        for N in (4, 6, 8), k in 0:6
            H = ising1D(N, 0.5)
            tm = trace_moment(H, k; scale=1)
            @test _re(tm, trace_moment(H, k; scale=1, multithreaded=true)) < 1e-10
            if k >= 1
                @test _re(tm, trace_product(H, k; scale=1)) < eps
            end
            N <= 6 && @test _re(tm, trace(H^k) / 2.0^N) < eps
        end
        for N in (4, 6), k in 0:6
            H = heisenberg_ring_plain(N)
            tm = trace_moment(H, k; scale=1)
            k >= 1 && @test _re(tm, trace_product(H, k; scale=1)) < eps
            @test _re(tm, trace(H^k) / 2.0^N) < eps
        end
    end

    @testset "complex coeffs, scale, k=0, errors" begin
        Hc = (1.0 + 0.5im) * ising1D(6, 0.5)
        for k in 1:5
            @test _re(trace_moment(Hc, k; scale=1), trace_product(Hc, k; scale=1)) < eps
        end
        H = ising1D(6, 0.5)
        @test _re(trace_moment(H, 4), trace_product(H, 4)) < eps   # default scale = 2^N
        @test trace_moment(H, 0; scale=1) ≈ 1
        @test trace_moment(H, 0) ≈ 2.0^6
        @test_throws ArgumentError trace_moment(H, -1)
    end

    @testset "4-arg trace_moment(A,k,B,l) == trace_product == dense" begin
        for N in (4, 6)
            A = ising1D(N, 0.5)
            B = Operator(N) + ("X", 1)
            for k in 1:5, l in 1:2
                tm = trace_moment(A, k, B, l; scale=1)
                @test _re(tm, trace_product(A, k, B, l; scale=1)) < eps
                @test _re(tm, trace_moment(A, k, B, l; scale=1, multithreaded=true)) < 1e-10
                N <= 6 && @test _re(tm, trace(A^k * B^l) / 2.0^N) < eps
            end
        end
        # two general operators (exercises the alternate tabulation branch / R != identity)
        for N in (4, 6), (k, l) in ((2, 2), (3, 1), (1, 3))
            A = heisenberg_ring_plain(N)
            B = ising1D(N, 0.7)
            @test _re(trace_moment(A, k, B, l; scale=1), trace_product(A, k, B, l; scale=1)) < eps
        end
        # k=0 / l=0 reduce to the single-operator moment
        let H = ising1D(6, 0.5), O = Operator(6) + ("X", 1)
            @test _re(trace_moment(H, 0, O, 3; scale=1), trace_moment(O, 3; scale=1)) < eps
            @test _re(trace_moment(H, 4, O, 0; scale=1), trace_moment(H, 4; scale=1)) < eps
        end
    end
end
