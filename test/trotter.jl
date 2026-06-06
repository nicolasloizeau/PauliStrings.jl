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
    truncation(O) = trim(O, 2^14)
    for order in (1, 2)
        g = trotterize(H, dt; heisenberg=true, order=order)
        Ot = copy(O)
        trotter_step!(Ot, g; truncation=truncation)
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
    truncation(O) = trim(O, 2^14)
    O1 = copy(O)
    trotter_step!(O1, trotterize(LocalH, dt; heisenberg=true, order=1); truncation=truncation)
    O2 = copy(O)
    trotter_step!(O2, trotterize(LocalH, dt; heisenberg=true, order=2); truncation=truncation)
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
    trotter_step!(O, TrotterGate{P2,Float64}[])
    @test norm(O - O0) < 1e-14
end

@testset "trotterize - OperatorTS produces minimal gates" begin
    N = 6
    H = Operator(N)
    H += "Z", 1, "Z", 2
    H += -0.5, "X", 1
    Hts = OperatorTS{(N,)}(H)

    g_ts = trotterize(Hts, 0.1; order=1)
    g_full = trotterize(resum(Hts), 0.1; order=1)

    # TS version has only M representative gates; full has N×M
    @test length(g_ts) == length(Hts)
    @test length(g_full) == length(resum(Hts))
    @test length(g_ts) < length(g_full)

    # order=2
    g_ts2 = trotterize(Hts, 0.1; order=2)
    @test length(g_ts2) == 2 * length(Hts) - 1

    # Empty operator
    Hempty = OperatorTS{(N,)}(Operator(N))
    @test isempty(trotterize(Hempty, 0.1))

    # Invalid order
    @test_throws ArgumentError trotterize(Hts, 0.1; order=3)
end

@testset "trotter_step! - OperatorTS matches resum approach" begin
    N = 6
    H = Operator(N)
    H += "Z", 1, "Z", 2
    H += -0.3, "X", 1
    Hts = OperatorTS{(N,)}(H)
    O = Operator(N)
    O += "X", 1
    Ots = OperatorTS{(N,)}(O)

    dt = 0.05
    truncation(O) = trim(O, 2^14)

    for order in (1, 2)
        # TS path: trotterize on OperatorTS, trotter_step! on OperatorTS
        Ots_copy = copy(Ots)
        g_ts = trotterize(Hts, dt; order=order)
        trotter_step!(Ots_copy, g_ts; truncation=truncation)

        # resum path: trotterize on full Operator, trotter_step! on full Operator
        Or = representative(copy(Ots))
        g_full = trotterize(resum(Hts), dt; order=order)
        trotter_step!(Or, g_full; truncation=truncation)
        Ots_from_full = OperatorTS{(N,)}(Or)

        # Both should give the same result
        @test norm(resum(Ots_copy) - resum(Ots_from_full)) / norm(resum(Ots_from_full)) < 1e-10
    end
end

@testset "trotter_step! - OperatorTS vs exact Heisenberg" begin
    N = 4
    H = Operator(N)
    H += "Z", 1, "Z", 2
    H += -0.5, "X", 1
    Hts = OperatorTS{(N,)}(H)
    O = Operator(N)
    O += "X", 1
    Ots = OperatorTS{(N,)}(O)

    dt = 0.02
    Mex = heisenberg_exact(resum(Hts), resum(Ots), dt)

    Ots_copy = copy(Ots)
    g = trotterize(Hts, dt; order=2)
    trotter_step!(Ots_copy, g)
    @test norm(Matrix(resum(Ots_copy)) - Mex) < 1e-6
end

@testset "trotter_step! - OperatorTS empty gates unchanged" begin
    N = 4
    Ots = OperatorTS{(N,)}(Operator(N) + ("Z", 1))
    Ots0 = copy(Ots)
    P = paulistringtype(N)
    trotter_step!(Ots, TrotterGate{P,Float64}[])
    @test norm(Ots - Ots0) < 1e-14
end