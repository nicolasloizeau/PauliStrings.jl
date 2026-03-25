

Ns = (10, 70, 200)

@testset "strings operations" begin
    # operations between Operator and PauliString
    for N in Ns
        O = rand_local2_M(N, 40)
        string = random_string(N)
        O2 = Operator(string)
        @test norm(O * string - O * O2) < 1E-10
        @test norm(string * O - O2 * O) < 1E-10
        @test norm(O * string -(string*O')') < 1E-10
        @test norm(commutator(O, string) - commutator(O, O2)) < 1E-10
        @test norm(anticommutator(O, string) - anticommutator(O, O2)) < 1E-10
        @test norm(commutator(O, string) + commutator(string, O)) < 1E-10
        @test norm(O + string - (O + O2)) < 1E-10
        @test norm(O - string - (O - O2)) < 1E-10
        string = rand(O2.strings)
        O2 = Operator(string)
        @test abs(trace_product(O, string; scale=1) - trace_product(O, O2; scale=1)) < 1E-10
    end
    # operations between PauliStringTS and OperatorTS
    for N in Ns
        O1 = compress(OperatorTS{(N,)}(rand_local2_M(N, 40)))
        string2 = PauliStringTS{(N,)}(random_string(N))
        O2 = OperatorTS(string2)
        @test norm(O1 * string2 - O1 * O2) < 1E-10
        @test norm(string2*O1 - O2*O1) < 1E-10
        @test norm(O1 * string2 -(string2*O1')') < 1E-10
        @test norm(commutator(O1, string2) - commutator(O1, O2)) < 1E-10
        @test norm(anticommutator(O1, string2) - anticommutator(O1, O2)) < 1E-10
        @test norm(O1 + string2 - (O1 + O2)) < 1E-10
        @test norm(O1 - string2 - (O1 - O2)) < 1E-10
        string2 = rand(O2.strings)
        O2 = OperatorTS(string2)
        @test abs(trace_product(O1, string2; scale=1) - trace_product(O1, O2; scale=1)) < 1E-10
    end
    # operations between PauliStringTS and PauliStringTS
    for N in Ns
        string1 = PauliStringTS{(N,)}(random_string(N))
        string2 = PauliStringTS{(N,)}(random_string(N))
        O1 = OperatorTS(string1)
        O2 = OperatorTS(string2)
        @test norm(string1 * string2 - (O1 * O2)) < 1E-10
        @test norm(string1 + string2 - (O1 + O2)) < 1E-10
        @test norm(string1 - string2 - (O1 - O2)) < 1E-10
        @test norm(commutator(string1, string2) - commutator(O1, O2)) < 1E-10
        @test norm(anticommutator(string1, string2) - anticommutator(O1, O2)) < 1E-10
        @test norm(trace_product(string1, string2)-trace_product(O1, O2)) < 1E-10
        @test norm(trace_product(string1, string1)-trace_product(O1, O1)) < 1E-10
        @test norm(string1) == norm(O1)
        @test norm(string1; normalize=true) == norm(O1; normalize=true)
    end
    # operations between PauliString and PauliString
    for N in Ns
        string1 = random_string(N)
        string2 = random_string(N)
        O1 = Operator(string1)
        O2 = Operator(string2)
        @test norm(string1 + string2 - (O1 + O2)) < 1E-10
        @test norm(string1 - string2 - (O1 - O2)) < 1E-10
        @test norm(trace_product(string1, string2)-trace_product(O1, O2)) < 1E-10
        @test norm(trace_product(string1, string1)-trace_product(O1, O1)) < 1E-10
        @test norm(string1) == norm(O1)
        @test norm(string1; normalize=true) == norm(O1; normalize=true)
    end
    for N in Ns
        for k in 1:4
            string1 = random_string(N)
            string2 = random_string(N)
            @test abs(trace_product(string1, string2) - trace_product(Operator(string1), Operator(string2))) < 1E-10
        end
    end

    # cis and exp for PauliString should match matrix exponential (cis also unitary)
    paulis = ['1', 'X', 'Y', 'Z']

    # all strings of length 1 and 2
    for N in (1, 2)
        inds = Iterators.product(ntuple(_ -> paulis, N)...)
        for chars in inds
            s = String(collect(chars))
            p = PauliString(s)
            U = cis(p)
            Udense = Matrix(U)
            Pmat = Matrix(p)
            Uref = exp(1im * Pmat)
            @test maximum(abs.(Udense .- Uref)) < 1e-12
            @test maximum(abs.(Udense' * Udense .- I(size(Udense, 1)))) < 1e-12

            E = exp(p)
            Edense = Matrix(E)
            Eref = exp(Pmat)
            @test maximum(abs.(Edense .- Eref)) < 1e-12
        end
    end

    # a few random 3-qubit strings
    for _ in 1:10
        p = random_string(3)
        U = cis(p)
        Udense = Matrix(U)
        Pmat = Matrix(p)
        Uref = exp(1im * Pmat)
        @test maximum(abs.(Udense .- Uref)) < 1e-12
        @test maximum(abs.(Udense' * Udense .- I(size(Udense, 1)))) < 1e-12

        E = exp(p)
        Edense = Matrix(E)
        Eref = exp(Pmat)
        @test maximum(abs.(Edense .- Eref)) < 1e-12
    end

    # pauli_rotation should match dense Heisenberg evolution for small systems
    angles = (0.1, 0.3, 0.7)

    # All nontrivial 1-qubit pairs
    one_qubit = ["X", "Y", "Z"]
    for Gs in one_qubit, Ps in one_qubit, theta in angles
        G = PauliString(Gs)
        P = PauliString(Ps)
        R = pauli_rotation(G, P, theta)
        Rdense = Matrix(R)

        Gmat = Matrix(G)
        Pmat = Matrix(P)
        U = exp(1im * theta * Gmat / 2)
        Rref = U * Pmat * U'
        @test maximum(abs.(Rdense .- Rref)) < 1e-12
    end

    # A few random 2-qubit non-commuting pairs
    for theta in angles
        for _ in 1:10
            G = random_string(2)
            P = random_string(2)
            _, kPG = commutator(P, G)
            kPG == 0 && continue  # skip commuting
            R = pauli_rotation(G, P, theta)
            Rdense = Matrix(R)

            Gmat = Matrix(G)
            Pmat = Matrix(P)
            U = exp(1im * theta * Gmat / 2)
            Rref = U * Pmat * U'
            @test maximum(abs.(Rdense .- Rref)) < 1e-12
        end
    end

end
