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




function ising(L, hx, hy)
    H = Operator(L )
    H += "Z", 1, "Z", 2
    H += hx, "X", 1
    H += hy, "Y", 1
    return OperatorTS{(L,)}(H)
end

ising(L) = ising(L, 0.5, 0.25)

observer(O::Operator) = sum(get_coeffs(xpart(O)))
observer(O::Operator{<:PauliStringTS}) = sum(get_coeffs(xpart(O)))*qubitlength(O)

function evolve_rk4(H, O, dt, nsteps; M=2^20, noise=0.0, observer=observer)
    result = []
    O0 = deepcopy(O)
    for i in ProgressBar(1:nsteps)
        push!(result, observer(O))
        O = ps.rk4(H, O, dt; heisenberg=true, M=M,  keep=O0)
        O = ps.trim(O, M; keep=O0)
    end
    return real.(result)
end




@testset "evolve_trotter vs rk4" begin
    N = 6
    O = Operator(N) + ("X", 1)
    O = OperatorTS{(N,)}(O)
    H = ising(N)
    M = 14
    dt = 0.02
    tmax = 2
    times = 0:dt:tmax

    res1 = evolve_trotter(H, deepcopy(O), dt, length(times);
        M=2^M, observer=observer)
    res2 = evolve_rk4(H, O, dt, length(times); M=2^M, observer=observer)
    @test norm(res1 - res2) < 0.003
    println("trotter vs rk4 error with M=$M: ", norm(res1 - res2))

    O = resum(O)
    H = resum(H)

    res1 = evolve_trotter(H, deepcopy(O), dt, length(times);
        M=2^M, observer=observer)
    res2 = evolve_rk4(H, O, dt, length(times); M=2^M, observer=observer)
    @test norm(res1 - res2) < 0.003
    println("trotter vs rk4 error with M=$M: ", norm(res1 - res2))

end
