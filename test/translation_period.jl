using PauliStrings
using Test
using LinearAlgebra

function naive_ts_sum(op::Operator, Ls, Ps, Ks)
    out = Operator(similar(op.strings, 0), similar(op.coeffs, 0))
    for s in PauliStrings.all_shifts(Ls, Ps, Ks)
        out += shift(op, Ls, Ps, s)
    end
    return out
end

@testset "translation period k" begin

    @testset "1D k=2" begin
        N = 6
        k = 2
        Ls = (N,)
        Ps = (true,)
        Ks = (k,)

        local_op = Operator("X11111")
        O_ts = OperatorTS{Ls, Ps, Ks}(local_op)
        O_naive = naive_ts_sum(local_op, Ls, Ps, Ks)

        @test translationperiod(O_ts) == Ks
        @test length(resum(O_ts)) == length(O_naive)
        @test norm(resum(O_ts) - O_naive) < 1e-12

        pts = PauliStringTS{Ls, Ps, Ks}("X11111")
        pts_shift = PauliStringTS{Ls, Ps, Ks}(PauliString("1X1111"))
        @test pts != pts_shift
        @test pts == PauliStringTS{Ls, Ps, Ks}(PauliString("1111X1"))

        @test length(resum(OperatorTS(pts))) == N ÷ k
    end

    @testset "1D k=1 backward compat" begin
        N = 4
        local_op = Operator("XX11")
        O_ts = OperatorTS{(N,)}(local_op)
        O_k1 = OperatorTS{(N,), (true,), (1,)}(local_op)

        @test translationperiod(O_ts) == (1,)
        @test norm(resum(O_ts) - resum(O_k1)) < 1e-12
        @test norm(resum(O_ts) - naive_ts_sum(local_op, (N,), (true,), (1,))) < 1e-12
    end

    @testset "2D k=(2,1)" begin
        Ls = (4, 2)
        Ps = (true, true)
        Ks = (2, 1)

        local_op = Operator("X1111111")
        O_ts = OperatorTS{Ls, Ps, Ks}(local_op)
        O_naive = naive_ts_sum(local_op, Ls, Ps, Ks)

        @test length(resum(O_ts)) == length(O_naive)
        @test norm(resum(O_ts) - O_naive) < 1e-12
        @test PauliStrings.num_translations(Ls, Ps, Ks) == 4
    end

    @testset "trace and products vs naive" begin
        N = 8
        k = 4
        Ls = (N,)
        Ps = (true,)
        Ks = (k,)

        O1 = OperatorTS{Ls, Ps, Ks}(Operator("X1111111"))
        O2 = OperatorTS{Ls, Ps, Ks}(Operator("Z1111111") + 0.5 * Operator("Y1111111"))

        R1 = resum(O1)
        R2 = resum(O2)
        R3 = resum(O1 * O2)

        @test norm(R3 - R1 * R2) < 1e-12
        @test abs(trace_product(O1, O2)) < 1e-12
        @test abs(trace_product(O1, O1) - trace(R1' * R1)) < 1e-12
        @test abs(trace(O1) - trace(R1)) < 1e-12
        @test abs(LinearAlgebra.norm(O1) - LinearAlgebra.norm(R1)) < 1e-12
    end

    @testset "unsigned type constructor compat" begin
        N = 4
        p = PauliString("XX11")
        pts = PauliStringTS{(N,), (true,), UInt64}(p)
        @test translationperiod(pts) == (1,)
        @test string(pts) == string(PauliStringTS{(N,), (true,), (1,)}(p))
    end

    @testset "norm and symmetry guards" begin
        Ls = (6,)
        Ps = (true,)
        Ks = (2,)
        pts = PauliStringTS{Ls, Ps, Ks}("X11111")
        @test norm(pts) ≈ sqrt(2.0^6 * 3)
        @test norm(pts; normalize=true) ≈ sqrt(3)

        O6 = OperatorTS{Ls, Ps, Ks}(Operator("X11111"))
        O8 = OperatorTS{Ls, Ps, (3,)}(Operator("X11111"))
        @test_throws DimensionMismatch O6 * O8
        @test_throws DimensionMismatch trace_product(O6, O8)
        @test is_ts(resum(O6), Ls, Ps, Ks)
        @test !is_ts(Operator("X11111"), Ls, Ps, Ks)

        support = k_local_basis_1d(6, 1; translational_symmetry=true, translation_period=2)
        @test all(translationperiod.(support) .== Ref(Ks))
    end

    @testset "API gaps" begin
        Ls = (6,)
        Ps = (true,)
        Ks = (2,)

        pts = PauliStringTS{Ls, Ps, Ks}("X", 1)
        @test pts == PauliStringTS{Ls, Ps, Ks}("X11111")

        O_ts = OperatorTS1D(Operator("X11111"); full=false, translation_period=2)
        @test translationperiod(O_ts) == Ks

        basis = symmetry_adapted_k_local_basis_1d(6, 1; translational_symmetry=true, translation_period=2)
        @test all(translationperiod.(basis) .== Ref(Ks))

        H = OperatorTS{Ls, Ps, Ks}(Operator("X11111"))
        bad_support = k_local_basis_1d(6, 1; translational_symmetry=true, translation_period=1)
        @test_throws DimensionMismatch lioms(H, bad_support)

        H = Operator(6)
        H += "X", 1
        H += "Z", 1, "Z", 2
        Hrk = OperatorTS{Ls, Ps, Ks}(H)
        Ork = OperatorTS{Ls, Ps, Ks}(Operator("X11111"))
        dt = 0.05
        O1 = evolve(Hrk, Ork, 0:dt:0.2; method=RK4(), fout=identity).final
        @test translationperiod(O1) == Ks
        @test isa(O1, typeof(Ork))
    end

end
