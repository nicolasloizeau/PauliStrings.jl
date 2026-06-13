
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
    @test all(op_to_dense(UGate(1, 1, theta, phi, lam)) .≈ [[cos(theta / 2), exp(1im * phi) * sin(theta / 2)] [-exp(1im * lam) * sin(theta / 2), exp(1im * (phi + lam)) * cos(theta / 2)]])
end

@testset "OpenQASM import" begin
    cat_state = """
    OPENQASM 2.0;
    include "qelib1.inc";
    qreg q[4];
    creg c[4];
    h q[0];
    cx q[0], q[1];
    cx q[1], q[2];
    cx q[2], q[3];
    measure q -> c;
    """

    @test_throws PauliStrings.Circuits.OpenQASMError from_openqasm(cat_state)
    circuit = from_openqasm(cat_state; ignore_measurements = true)
    @test circuit.N == 4
    @test circuit.gates == [
        ("H", [1], Real[]),
        ("CNOT", [1, 2], Real[]),
        ("CNOT", [2, 3], Real[]),
        ("CNOT", [3, 4], Real[]),
    ]
    @test isunitary(compile(circuit))

    register_program = """
    OPENQASM 2.0;
    include "qelib1.inc";
    qreg left[2];
    qreg right[2];
    h left;
    cx left, right;
    rz(pi / 2 + sin(pi / 2)) right[1];
    """
    register_circuit = from_openqasm(register_program)
    @test register_circuit.N == 4
    @test register_circuit.gates[1:4] == [
        ("H", [1], Real[]),
        ("H", [2], Real[]),
        ("CNOT", [1, 3], Real[]),
        ("CNOT", [2, 4], Real[]),
    ]
    @test register_circuit.gates[5][1:2] == ("RZ", [4])
    @test only(register_circuit.gates[5][3]) ≈ pi / 2 + 1

    aliases = from_openqasm("""
    OPENQASM 2.0;
    qreg q[3];
    sdg q[0];
    u2(0, pi) q[1];
    cp(pi / 4) q[1], q[2];
    ccx q[0], q[1], q[2];
    """)
    @test [gate[1] for gate in aliases.gates] == ["Phase", "U", "CPhase", "CCX"]

    temporary_path, io = mktemp()
    try
        write(io, "OPENQASM 2.0; qreg q[1]; x q[0];")
        close(io)
        @test from_openqasm_file(temporary_path).gates == [("X", [1], Real[])]
    finally
        isopen(io) && close(io)
        rm(temporary_path; force = true)
    end

    @test_throws PauliStrings.Circuits.OpenQASMError from_openqasm(
        "OPENQASM 2.0; qreg q[1]; x q[1];",
    )
    @test_throws PauliStrings.Circuits.OpenQASMError from_openqasm(
        "OPENQASM 2.0; qreg q[1]; made_up q[0];",
    )
end
