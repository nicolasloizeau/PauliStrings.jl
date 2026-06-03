
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
    @test isidentity(U1 * U2')
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
    qasm = """
    OPENQASM 2.0;
    include "qelib1.inc";
    qreg q[3];
    creg c[3];

    // Bell pair plus a rotation and Toffoli.
    h q[0];
    cx q[0], q[1];
    rz(pi/2) q[1];
    ccx q[0], q[1], q[2];
    sdg q[2];
    u2(pi/4, pi/8) q[0];
    cu1(pi/3) q[0], q[2];
    id q[1];
    barrier q;
    measure q[0] -> c[0];
    """
    c = from_openqasm(qasm)
    @test c.N == 3
    @test first.(c.gates) == ["H", "CX", "RZ", "CCX", "Phase", "U", "CPhase"]
    @test [sites for (_, sites, _) in c.gates] == [[1], [1, 2], [2], [1, 2, 3], [3], [1], [1, 3]]
    @test c.gates[3][3] == Real[pi / 2]
    @test c.gates[5][3] == Real[-pi / 2]
    @test c.gates[6][3] == Real[pi / 2, pi / 4, pi / 8]
    @test c.gates[7][3] == Real[pi / 3]
    @test isunitary(compile(c))

    path = tempname() * ".qasm"
    write(path, qasm)
    @test from_openqasm_file(path).gates == c.gates
    rm(path)

    multi_register = from_openqasm("""
    OPENQASM 2.0;
    qreg a[1];
    qreg b[2];
    x a[0];
    cx a[0], b[1];
    """)
    @test multi_register.N == 3
    @test [sites for (_, sites, _) in multi_register.gates] == [[1], [1, 3]]

    @test_throws ErrorException from_openqasm("""
    OPENQASM 2.0;
    qreg q[1];
    qreg q[1];
    """)

    @test_throws ErrorException from_openqasm("""
    OPENQASM 2.0;
    qreg q[1];
    gate custom a { x a; }
    custom q[0];
    """)
    @test all(op_to_dense(UGate(1, 1, theta, phi, lam)) .≈ [[cos(theta / 2), exp(1im * phi) * sin(theta / 2)] [-exp(1im * lam) * sin(theta / 2), exp(1im * (phi + lam)) * cos(theta / 2)]])
end
