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

    # 4. Test custom translation periods
    L1 = (6,)
    Ps1 = (true,)
    PsOpen1 = (false,)
    Ks1 = (2,)
    p1 = PauliString("X11111")
    pts_period2 = PauliStringTS{L1,Ps1,Ks1}(p1)
    pts_period2_shift = PauliStringTS{L1,Ps1,Ks1}(shift(p1, L1, Ps1, (2,)))
    pts_period2_other_orbit = PauliStringTS{L1,Ps1,Ks1}(shift(p1, L1, Ps1, (1,)))

    @test pts_period2 == pts_period2_shift
    @test pts_period2 != pts_period2_other_orbit
    @test translationperiods(pts_period2) == Ks1
    @test length(resum(OperatorTS(pts_period2))) == 3
    @test trace_product(pts_period2, pts_period2) ≈ 2.0^qubitlength(pts_period2) * 3
    @test norm(pts_period2) ≈ sqrt(2.0^qubitlength(pts_period2) * 3)
    @test norm(pts_period2; normalize=true) ≈ sqrt(3)
    pts_period2_u16 = PauliStringTS{L1,Ps1,Ks1}(PauliString{6,UInt16}("X11111"))
    @test trace_product(pts_period2, pts_period2_u16) ≈ trace_product(pts_period2, pts_period2)

    @test PauliStringTS{L1,(true,),UInt8}(p1) == PauliStringTS{L1}(p1)
    @test PauliStringTS{L1,PsOpen1,UInt8}(p1) == PauliStringTS{L1,PsOpen1}(p1)
    @test PauliStringTS{L1,UInt8}(p1) == PauliStringTS{L1}(p1)
    @test PauliStringTS{L1,UInt16}(PauliString{6,UInt16}("X11111")) isa PauliStringTS{L1,Ps1,UInt16,UInt16}
    @test PauliStringTS{L1}(p1) isa PauliStringTS{L1,(true,),UInt8}
    @test PauliStringTS{L1,PsOpen1}(p1) isa PauliStringTS{L1,PsOpen1,UInt8}
    @test pts_period2 isa PauliStringTS{L1,Ps1,Ks1,UInt8}
    @test_throws ErrorException PauliStringTS{L1,Ps1,(4,)}(p1)

    O6 = OperatorTS{L1,Ps1,Ks1}(Operator("X11111") + 0.25 * Operator("Y11111"))
    O7 = OperatorTS{L1,Ps1,Ks1}(Operator("Z11111") + 0.5 * Operator("X11111"))
    old_string_dispatch(::PauliStringTS{L1,(true,),UInt8}) = :old_string_default_flags
    old_string_dispatch(::PauliStringTS{L1,PsOpen1,UInt8}) = :old_string_with_flags
    old_operator_dispatch(::OperatorTS{L1,Ps1,UInt8,T}) where {T} = :old_operator_with_flags
    @test old_string_dispatch(PauliStringTS{L1}(p1)) == :old_string_default_flags
    @test old_string_dispatch(PauliStringTS{L1,PsOpen1}(p1)) == :old_string_with_flags
    @test old_string_dispatch(PauliStringTS{L1,UInt8}(p1)) == :old_string_default_flags
    @test OperatorTS{L1,UInt8}(Operator("X11111")) isa OperatorTS{L1,Ps1,UInt8,ComplexF64}
    @test old_operator_dispatch(OperatorTS{L1,Ps1}(Operator("X11111"))) == :old_operator_with_flags
    @test old_operator_dispatch(OperatorTS{L1,UInt8}(Operator("X11111"))) == :old_operator_with_flags
    @test O6 isa OperatorTS{L1,Ps1,Ks1}
    @test translationperiods(O6) == Ks1
    @test length(resum(O6)) == 6
    @test is_ts(resum(O6), L1, Ps1, Ks1)
    @test !is_ts(Operator("X11111"), L1, Ps1, Ks1)
    @test norm(resum(O6 * O7) - resum(O6) * resum(O7)) < 1e-12
    @test abs(trace_product(O6, O7) - trace(resum(O6)' * resum(O7))) < 1e-12

    support_period2 = k_local_basis_1d(6, 1; translational_symmetry=true, translation_period=2)
    @test all(translationperiods.(support_period2) .== Ref(Ks1))

    Ks2 = (3,)
    pts_period3 = PauliStringTS{L1,Ps1,Ks2}(p1)
    O8 = OperatorTS{L1,Ps1,Ks2}(Operator("X11111"))
    @test_throws DimensionMismatch O6 * O8
    @test_throws DimensionMismatch O6 * pts_period3
    @test_throws DimensionMismatch pts_period2 * pts_period3
    @test_throws DimensionMismatch trace_product(O6, O8)
    @test_throws DimensionMismatch trace_product(O6, pts_period3)
    @test_throws DimensionMismatch trace_product(pts_period2, pts_period3)
    @test_throws DimensionMismatch evolve(O6, O8, 0.01, 1; method=Trotter(), fout=identity)

end
