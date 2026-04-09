using PauliStrings
using LinearAlgebra
using Test

"""One Heisenberg step: O -> exp(i H dt / hbar) O exp(-i H dt / hbar)."""
function heisenberg_exact(H::Operator, O::Operator, dt::Real; hbar::Real=1)
    Hm = Matrix(H)
    Om = Matrix(O)
    U = exp((1im * dt / hbar) * Hm)
    return U * Om * adjoint(U)
end

@testset "trotterize / gate list" begin
    H = Operator(2)
    H += 0.5, "X", 1
    H += 0.3, "Z", 2
    g1 = trotterize(H, 0.1; order=1, heisenberg=true)
    g2 = trotterize(H, 0.1; order=2, heisenberg=true)
    @test length(g1) == 2
    @test length(g2) == 3
    @test g2[2].theta == g1[2].theta
    @test g2[1].theta ≈ g1[1].theta / 2
    @test g2[3].theta ≈ g1[1].theta / 2

    H1 = Operator(1)
    H1 += 1.2, "X", 1
    a = trotterize(H1, 0.07; order=1)
    b = trotterize(H1, 0.07; order=2)
    @test length(a) == length(b) == 1
    @test a[1].theta == b[1].theta

    @test isempty(trotterize(Operator(2), 0.05))
    @test_throws ArgumentError trotterize(H, 0.05; order=3)
    @test_throws ArgumentError trotterize(H, 0.05; order=0)

    Hc = Operator(1)
    Hc += 1im, "X", 1
    @test_throws ArgumentError trotterize(Hc, 0.05)

    gh = trotterize(H1, 0.11; heisenberg=true)
    gs = trotterize(H1, 0.11; heisenberg=false)
    @test gh[1].theta == -gs[1].theta
end

@testset "trotter_step! - commuting terms vs exact Heisenberg" begin
    H = Operator(2)
    H += 0.7, "X", 1
    H += 0.3, "Z", 2
    O = Operator(2)
    O += "X", 2
    dt = 0.05
    keep = Operator(2)
    for order in (1, 2)
        g = trotterize(H, dt; heisenberg=true, order=order)
        Ot = copy(O)
        trotter_step!(Ot, g; M=10^6, keep=keep)
        Mex = heisenberg_exact(H, O, dt)
        @test norm(Matrix(Ot) - Mex) < 1e-10
    end
end

@testset "trotter_step! - Strang is closer than Lie for non-commuting H" begin
    LocalH = Operator(1)
    LocalH += 0.4, "X", 1
    LocalH += 0.5, "Z", 1
    O = Operator(1)
    O += "Y", 1
    dt = 0.15
    Mex = heisenberg_exact(LocalH, O, dt)
    keep = Operator(1)
    O1 = copy(O)
    trotter_step!(O1, trotterize(LocalH, dt; heisenberg=true, order=1); M=10^6, keep=keep)
    O2 = copy(O)
    trotter_step!(O2, trotterize(LocalH, dt; heisenberg=true, order=2); M=10^6, keep=keep)
    e1 = norm(Matrix(O1) - Mex)
    e2 = norm(Matrix(O2) - Mex)
    @test e2 < e1
    @test e1 > 1e-4
end

@testset "trotter_step! - empty gate list leaves O unchanged" begin
    O = Operator(2)
    O += "Z", 1
    O0 = copy(O)
    P2 = paulistringtype(2)
    trotter_step!(O, TrotterGate{P2,Float64}[]; M=100, keep=Operator(2))
    @test norm(O - O0) < 1e-14
end

@testset "trotter_evolve" begin
    H = Operator(2)
    H += 0.2, "Z", 1
    H += 0.15, "X", 2
    O = Operator(2)
    O += "X", 1
    O += "Z", 2
    O0 = copy(O)
    trotter_evolve(H, O, 0.02, 0; heisenberg=true, order=2)
    @test norm(O - O0) < 1e-14

    H2 = Operator(1)
    H2 += 0.5, "X", 1
    Oa = Operator(1)
    Oa += "Z", 1
    Ob = copy(Oa)
    g = trotterize(H2, 0.03; heisenberg=true, order=1)
    trotter_evolve(H2, Oa, 0.03, 4; gates=g, heisenberg=true, order=1, M=10^6, keep=Operator(1))
    for _ in 1:4
        trotter_step!(Ob, g; M=10^6, keep=Operator(1))
    end
    @test norm(Matrix(Oa) - Matrix(Ob)) < 1e-12

    H3 = Operator(2)
    O3 = Operator(3)
    O3 += "X", 1
    @test_throws DimensionMismatch trotter_evolve(H3, O3, 0.1, 1)

    @test_throws ArgumentError trotter_evolve(H2, Operator(1), 0.1, -1)
end
