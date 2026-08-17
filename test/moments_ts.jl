using PauliStrings
import PauliStrings as ps
using Test
using LinearAlgebra

# Heisenberg ring (nonzero odd moments, dense anticommuting cores) for the moment tests
function heisenberg_ring(N)
    H = Operator(N)
    for j in 1:N
        H += "X", j, "X", mod1(j + 1, N)
        H += "Y", j, "Y", mod1(j + 1, N)
        H += "Z", j, "Z", mod1(j + 1, N)
    end
    return H
end

_relerr(a, b) = abs(a - b) / max(abs(a), abs(b), 1.0)

@testset "trace_moment (translation symmetric)" begin
    eps = 1e-7

    @testset "1D TFIM: trace_moment == trace_product == dense" begin
        for N in (4, 6, 8), k in 0:6
            Hts = OperatorTS1D(ising1D(N, 0.5))
            tm = trace_moment(Hts, k; scale=1)
            if k >= 1
                @test _relerr(tm, trace_product(Hts, k; scale=1)) < eps
            end
            if N <= 6                       # independent dense (operator-power) oracle
                @test _relerr(tm, trace(resum(Hts)^k) / 2.0^N) < eps
            end
        end
    end

    @testset "1D Heisenberg (nonzero odd moments)" begin
        for N in (4, 6), k in 0:6
            Hts = OperatorTS1D(heisenberg_ring(N))
            tm = trace_moment(Hts, k; scale=1)
            k >= 1 && @test _relerr(tm, trace_product(Hts, k; scale=1)) < eps
            @test _relerr(tm, trace(resum(Hts)^k) / 2.0^N) < eps
        end
    end

    @testset "fold == unfold == multithreaded" begin
        for N in (4, 6), k in 1:6
            Hts = OperatorTS1D(ising1D(N, 0.7))
            f = trace_moment(Hts, k; scale=1, fold=true)
            u = trace_moment(Hts, k; scale=1, fold=false)
            mt = trace_moment(Hts, k; scale=1, multithreaded=true)
            @test _relerr(f, u) < 1e-10
            @test _relerr(f, mt) < 1e-10
        end
        # Heisenberg exercises the anticommuting-core DP in all three paths
        for k in 1:6
            Hts = OperatorTS1D(heisenberg_ring(6))
            f = trace_moment(Hts, k; scale=1, fold=true)
            @test _relerr(f, trace_moment(Hts, k; scale=1, fold=false)) < 1e-10
            @test _relerr(f, trace_moment(Hts, k; scale=1, multithreaded=true)) < 1e-10
        end
    end

    @testset "complex and non-Hermitian coefficients" begin
        Hc = (1.0 + 0.5im) * OperatorTS1D(ising1D(6, 0.5))
        for k in 1:5
            @test _relerr(trace_moment(Hc, k; scale=1), trace_product(Hc, k; scale=1)) < eps
        end
        Hnh = OperatorTS{(6,),(true,)}(Operator("X11111") + 0.3im * Operator("Z11111"))
        for k in 1:5
            @test _relerr(trace_moment(Hnh, k; scale=1), trace(resum(Hnh)^k) / 2.0^6) < eps
        end
    end

    @testset "scale keyword and k=0" begin
        Hts = OperatorTS1D(ising1D(6, 0.5))
        @test _relerr(trace_moment(Hts, 4), trace_product(Hts, 4)) < eps   # default scale = 2^N
        @test trace_moment(Hts, 0; scale=1) ≈ 1
        @test trace_moment(Hts, 0) ≈ 2.0^6
        @test trace_moment(Hts, 0; scale=2.5) ≈ 2.5
    end

    @testset "2D TFIM" begin
        for L1 in 2:3
            L2 = 3
            Hts = OperatorTS2D(ising2D(L1, L2, 0.6), L1)
            for k in 1:4
                tm = trace_moment(Hts, k; scale=1)
                @test _relerr(tm, trace_product(Hts, k; scale=1)) < eps
                @test _relerr(tm, trace(resum(Hts)^k) / 2.0^(L1 * L2)) < eps
            end
        end
    end

    @testset "mixed periodicity (one open dimension)" begin
        Hmix = OperatorTS{(2, 3),(true, false)}(Operator("X11111") + Operator("Z11111"))
        for k in 1:5
            tm = trace_moment(Hmix, k; scale=1)
            @test _relerr(tm, trace_product(Hmix, k; scale=1)) < eps
            @test _relerr(tm, trace(resum(Hmix)^k) / 2.0^6) < eps
        end
    end

    @testset "error handling" begin
        Hts = OperatorTS1D(ising1D(4, 0.5))
        @test_throws ArgumentError trace_moment(Hts, -1)
    end
end
