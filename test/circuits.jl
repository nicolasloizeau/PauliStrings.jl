
using PauliStrings.Circuits
using LinearAlgebra


function ccx()
    c = Circuit(3)
    push!(c, "H", 3)
    push!(c, "CNOT", 2, 3)
    push!(c, "Tdg", 3)
    push!(c, "CNOT", 1, 3)
    push!(c, "T", 3)
    push!(c, "CNOT", 2, 3)
    push!(c, "Tdg", 3)
    push!(c, "CNOT", 1, 3)
    push!(c, "T", 2)
    push!(c, "T", 3)
    push!(c, "CNOT", 1, 2)
    push!(c, "H", 3)
    push!(c, "T", 1)
    push!(c, "Tdg", 2)
    push!(c, "CNOT", 1, 2)
    return c
end

function swap()
    c = Circuit(2)
    push!(c, "CNOT", 1, 2)
    push!(c, "CNOT", 2, 1)
    push!(c, "CNOT", 1, 2)
    return compile(c)
end

function swapH()
    c = Circuit(2)
    push!(c, "H", 1)
    push!(c, "CZ", 1, 2)
    push!(c, "H", 1)
    push!(c, "H", 2)
    push!(c, "CZ", 1, 2)
    push!(c, "H", 1)
    push!(c, "H", 2)
    push!(c, "CZ", 1, 2)
    push!(c, "H", 1)
    return compile(c)
end

function cxH()
    c = Circuit(2)
    push!(c, "H", 2)
    push!(c, "CZ", 2, 1)
    push!(c, "H", 2)
    return compile(c)
end

import LinearAlgebra: diag as ldiag
@testset "circuits" begin
    U1 = compress(compile(ccx()))
    U2 = CCXGate(3, 1, 2, 3)
    @test norm(U1 - U2) < 1e-10
    @test isidentity(U1 * dagger(U2))
    @test isunitary(TGate(1, 1))
    @test isunitary(U1)
    @test isunitary(U2)
    @test isidentity(swap() * SwapGate(2, 1, 2))
    @test isidentity(swapH() * SwapGate(2, 1, 2))
    @test isidentity(cxH() * CXGate(2, 1, 2))
    @test norm(CZGate(2, 1, 2) - CZGate(2, 2, 1)) < 1e-10
    @test ldiag(op_to_dense(MCZGate(3))) == [1, 1, 1, 1, 1, 1, 1, -1]
    @test expect(ccx(), "111", "110") ≈ 1
    @test expect(ccx(), "111", "000") ≈ 0
    @test expect(U1, "111", "000") ≈ 0
    @test expect(U1, "111", "110") ≈ 1
    @test isunitary(UGate(2, 1, 0.1, 1.2, 0.3))
    @test isunitary(RYGate(2, 1, 0.1))
    @test isunitary(RXGate(2, 1, 0.1))
    @test norm(ldiag(op_to_dense(RZGate(2, 1, pi))) - [1im, 1im, -1im, -1im]) < 1e-10
    theta, phi, lam = 0.1, 1.2, 0.3
    @test all(op_to_dense(UGate(1, 1, theta, phi, lam)) .≈ [[cos(theta / 2), exp(1im * phi) * sin(theta / 2)] [-exp(1im * lam) * sin(theta / 2), exp(1im * (phi + lam)) * cos(theta / 2)]])
end
