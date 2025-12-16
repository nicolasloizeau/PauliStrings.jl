

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

end
