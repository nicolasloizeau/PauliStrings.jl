
using PauliStrings.Circuits

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
    return compile(c)
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
    U1 = compress(ccx())
    U2 = CCXGate(3, 1, 2, 3)
    @test opnorm(U1 - U2) < 1e-10
    @test isidentity(U1 * dagger(U2))
    @test isunitary(TGate(1,1))
    @test isunitary(U1)
    @test isunitary(U2)
    @test isidentity(swap() * SwapGate(2, 1, 2))
    @test isidentity(swapH() * SwapGate(2, 1, 2))
    @test isidentity(cxH() * CXGate(2, 1, 2))
    @test opnorm(CZGate(2, 1, 2)-CZGate(2, 2, 1)) < 1e-10
    @test ldiag(op_to_dense(MCZGate(3))) == [1,1,1,1,1,1,1,-1]
end
