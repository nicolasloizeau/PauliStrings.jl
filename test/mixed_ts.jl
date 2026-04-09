using PauliStrings
using Test
using LinearAlgebra

@testset "mixed ts" begin

    Ls = (2, 2)

    # 1. Test PauliStringTS construction and properties
    p = PauliString("X111")
    pts_per = PauliStringTS{Ls,(true, true)}(p)
    pts_mix = PauliStringTS{Ls,(true, false)}(p)
    pts_open = PauliStringTS{Ls,(false, false)}(p)

    pts_per_shift = PauliStringTS{Ls,(true, true)}(PauliString("1X11"))
    @test pts_per == pts_per_shift
    @test trace_product(pts_per, pts_per_shift) ≈ 2.0^qubitlength(pts_per) * 4

    @test periodicflags(pts_per) == (true, true)
    @test periodicflags(pts_mix) == (true, false)
    @test periodicflags(pts_open) == (false, false)

    # Check number of translations
    # resum should produce:
    # pts_per: 4 terms ( (0,0), (1,0), (0,1), (1,1) )
    # pts_mix: 2 terms ( (0,0), (1,0) )
    # pts_open: 1 term ( (0,0) )

    @test length(resum(OperatorTS(pts_per))) == 4
    @test length(resum(OperatorTS(pts_mix))) == 2
    @test length(resum(OperatorTS(pts_open))) == 1

    # 2. Test operator products with mixed periodicity
    O1 = OperatorTS{Ls,(true, false)}(Operator("X111"))
    O2 = OperatorTS{Ls,(true, false)}(Operator("Z111"))

    # O1 * O2 should be translation symmetric with same periodicity
    O3 = O1 * O2
    @test periodicflags(O3) == (true, false)
    @test qubitsize(O3) == Ls

    # Verify against resummed multiplication
    R1 = resum(O1)
    R2 = resum(O2)
    R3 = resum(O3)
    @test norm(R3 - R1 * R2) < 1e-12

    # 3. Test trace_product
    @test abs(trace_product(O1, O2)) < 1e-12 # X and Z are orthogonal
    @test abs(trace_product(O1, O1) - trace(resum(O1)' * resum(O1))) < 1e-12

    O4 = OperatorTS{Ls,(true, false)}(Operator("X111") + 0.5 * Operator("Y111"))
    O5 = OperatorTS{Ls,(true, false)}(Operator("X111") + 0.3 * Operator("Z111"))
    @test abs(trace_product(O4, O5) - trace(resum(O4)' * resum(O5))) < 1e-12

end
