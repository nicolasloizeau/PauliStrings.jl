using PauliStrings
using Test
using LinearAlgebra

function TFIM(N, h)
    H = Operator(N)
    H += -h, "X", 1
    H += "Z", 1, "Z", 2
    return OperatorTS{(N,)}(H)
end

# Total X operator
function Xtot(N)
    H = Operator(N)
    H += "X", 1
    return OperatorTS{(N,)}(H)
end


@testset "evolve" begin
    N = 6
    H = TFIM(N, 0.3)
    O0 = Xtot(N)
    dt = 0.01
    fout(O) = trace_product(O0, O) / 2^N
    fout_trotter(O) = trace_product(resum(O0), O) / 2^N
    fout_exact(O) = tr(Matrix(O0)*O) / 2^N
    times = 0:dt:1
    truncation(O) = trim(O, 2^14)
    res1 = evolve(H, O0, times; method = Exact(), fout=fout_exact).history
    res2 = evolve(H, O0, times; method = RK4(), fout=fout, truncation=truncation).history
    res3 = evolve(H, O0, times; method = DOPRI5(), fout=fout, truncation=truncation).history
    res4 = evolve(resum(H), resum(O0), times; method = Trotter(), fout=fout_trotter, truncation=truncation).history
    res5 = evolve(H, O0, times; method = Trotter(), fout=fout, truncation=truncation).history

    @test norm(res2 .- res1)/norm(res1) < 1e-7
    @test norm(res3 .- res1)/norm(res1) < 1e-7
    @test norm(res4 .- res1)/norm(res1) < 1e-5
    @test norm(res5 .- res1)/norm(res1) < 1e-4

    dissipation(O, dt) = add_noise(O, 0.1*dt)
    res1 = evolve(H, O0, times; method = RK4(), fout=fout, dissipation=dissipation).history
    res2 = evolve(H, O0, times; method = DOPRI5(), fout=fout, dissipation=dissipation).history
    res3 = evolve(H, O0, times; method = RK4(), fout=fout).history
    @test norm(res1 .- res2)/norm(res2) < 1e-7
    @test norm(res1) < norm(res3)
end


# The in-place steppers must reproduce their allocating counterparts (within
# machine precision) while overwriting `O` and returning it.
@testset "rk4! / dopri5! / rk4_lindblad! match allocating versions" begin
    relerr(a, b) = norm(a - b) / norm(b)
    dt = 0.01

    # plain Operator and translation-symmetric OperatorTS
    N = 4
    Hplain = Operator(N)
    Hplain += -0.3, "X", 1
    Hplain += "Z", 1, "Z", 2
    Hplain += 0.7, "Z", 3
    Oplain = Operator(N)
    Oplain += "X", 2
    Oplain += 0.5, "Y", 1, "Y", 3

    Hts = TFIM(6, 0.3)
    Ots = Xtot(6)

    trunc(O) = trim(O, 2^14)

    for (H, O) in ((Hplain, Oplain), (Hts, Ots))
        # heisenberg=true / false, with and without truncation
        for heisenberg in (true, false), truncation in (identity, trunc)
            ref = rk4(H, O, dt; heisenberg=heisenberg, truncation=truncation)
            inplace = rk4!(copy(O), H, dt; heisenberg=heisenberg, truncation=truncation)
            @test relerr(inplace, ref) < 1e-12

            ref = PauliStrings.dopri5(H, O, dt; heisenberg=heisenberg, truncation=truncation)
            inplace = PauliStrings.dopri5!(copy(O), H, dt; heisenberg=heisenberg, truncation=truncation)
            @test relerr(inplace, ref) < 1e-12
        end

        # rk4! returns the same object it was handed (true in-place)
        O2 = copy(O)
        @test rk4!(O2, H, dt) === O2

        # time-dependent Hamiltonian wrapper
        Hfun(t) = H
        @test relerr(rk4!(copy(O), Hfun, dt, 0.0), rk4(Hfun, O, dt, 0.0)) < 1e-12
    end

    # Lindblad path
    L = [Oplain, Operator(N) + ("Z", 1)]
    gamma = [0.4, 0.7]
    ref = rk4_lindblad(Hplain, Oplain, dt, L; gamma=gamma)
    inplace = rk4_lindblad!(copy(Oplain), Hplain, dt, L; gamma=gamma)
    @test relerr(inplace, ref) < 1e-12
end
