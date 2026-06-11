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
