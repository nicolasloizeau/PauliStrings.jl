using OpenQASM
using Test
using PauliStrings
using PauliStrings.Circuits
using LinearAlgebra

@testset "openqasm" begin
    eps = 1e-12

    @testset "gate set round-trips to a hand-built Circuit" begin
        qasm = """
        OPENQASM 2.0;
        include "qelib1.inc";
        qreg q[3];
        creg c[3];
        x q[0];
        y q[1];
        z q[2];
        h q[0];
        s q[1];
        t q[2];
        sx q[0];
        tdg q[1];
        sdg q[2];
        rx(0.2) q[0];
        ry(0.3) q[1];
        rz(0.4) q[2];
        u1(0.5) q[0];
        p(0.6) q[1];
        u2(0.1,0.2) q[2];
        u3(0.1,0.2,0.3) q[0];
        cx q[0],q[1];
        cy q[1],q[2];
        cz q[0],q[2];
        swap q[0],q[1];
        ccx q[0],q[1],q[2];
        measure q -> c;
        """
        c = parse_qasm(qasm)

        h = Circuit(3)
        push!(h, "X", 1); push!(h, "Y", 2); push!(h, "Z", 3)
        push!(h, "H", 1); push!(h, "S", 2); push!(h, "T", 3)
        push!(h, "SX", 1); push!(h, "Phase", 2, -pi / 4); push!(h, "Phase", 3, -pi / 2)
        push!(h, "RX", 1, 0.2); push!(h, "RY", 2, 0.3); push!(h, "RZ", 3, 0.4)
        push!(h, "Phase", 1, 0.5); push!(h, "Phase", 2, 0.6)
        push!(h, "U", 3, pi / 2, 0.1, 0.2); push!(h, "U", 1, 0.1, 0.2, 0.3)
        push!(h, "CNOT", 1, 2); push!(h, "CY", 2, 3); push!(h, "CZ", 1, 3)
        push!(h, "Swap", 1, 2); push!(h, "CCX", 1, 2, 3)

        @test c.N == 3
        @test norm(Matrix(compile(c)) - Matrix(compile(h))) < eps
    end

    @testset "measure and barrier are ignored" begin
        c = parse_qasm("OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[1];\ncreg c[1];\nh q[0];\nbarrier q;\nmeasure q -> c;\n")
        @test length(c.gates) == 1
        @test c.gates[1][1] == "H"
    end

    @testset "multiple qregs share one index space" begin
        c = parse_qasm("OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg a[2];\nqreg b[1];\nx a[0];\ncx a[1],b[0];\n")
        @test c.N == 3
        @test c.gates[1][2] == [1]
        @test c.gates[2][2] == [2, 3]
    end

    @testset "parameter expressions" begin
        c = parse_qasm("OPENQASM 2.0;\nqreg q[1];\nrz(2*pi) q[0];\nrx(-pi/4) q[0];\nry(pi) q[0];\nrz(sin(0.5)) q[0];\nrx(-0.3) q[0];\nrz(1+2*3) q[0];\n")
        got = [g[3][1] for g in c.gates]
        @test got ≈ [2pi, -pi / 4, pi, sin(0.5), -0.3, 7.0]
    end

    @testset "load_qasm reads from a file" begin
        path = tempname() * ".qasm"
        write(path, "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[2];\nh q[0];\ncx q[0],q[1];\n")
        c = load_qasm(path)
        @test c.N == 2
        @test [g[1] for g in c.gates] == ["H", "CNOT"]
        rm(path)
    end

    @testset "QASMBench-style Deutsch circuit imports and runs" begin
        deutsch = """
        OPENQASM 2.0;
        include "qelib1.inc";
        qreg q[2];
        creg c[2];
        x q[1];
        h q[0];
        h q[1];
        cx q[0],q[1];
        h q[0];
        measure q[0] -> c[0];
        measure q[1] -> c[1];
        """
        c = parse_qasm(deutsch)
        ref = Circuit(2)
        push!(ref, "X", 2); push!(ref, "H", 1); push!(ref, "H", 2)
        push!(ref, "CNOT", 1, 2); push!(ref, "H", 1)
        @test norm(Matrix(compile(c)) - Matrix(compile(ref))) < eps
        @test expect(c, "00") ≈ expect(ref, "00")
    end

    @testset "register broadcasting" begin
        # `h q;` applies H to every qubit of q; `cx q,r;` is element-wise
        c = parse_qasm("OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[3];\nqreg r[3];\nh q;\ncx q,r;\n")
        ref = Circuit(6)
        push!(ref, "H", 1); push!(ref, "H", 2); push!(ref, "H", 3)
        push!(ref, "CNOT", 1, 4); push!(ref, "CNOT", 2, 5); push!(ref, "CNOT", 3, 6)
        @test norm(Matrix(compile(c)) - Matrix(compile(ref))) < eps
        # single qubit broadcasts against a register: `cx q[0], r;`
        c2 = parse_qasm("OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[1];\nqreg r[2];\ncx q[0],r;\n")
        @test [g[2] for g in c2.gates] == [[1, 2], [1, 3]]
    end

    @testset "errors" begin
        @test_throws ArgumentError parse_qasm("OPENQASM 2.0;\n")                                  # no qreg
        @test_throws ArgumentError parse_qasm("OPENQASM 2.0;\nqreg q[1];\nx p[0];\n")             # unknown qreg
        @test_throws ArgumentError parse_qasm("OPENQASM 2.0;\nqreg q[1];\ncrz(0.1) q[0],q[0];\n") # unsupported gate
    end
end
