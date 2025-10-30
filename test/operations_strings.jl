

function random_string(N::Int)
    symbols = ['1', 'X', 'Y', 'Z']
    s = join([symbols[rand(1:4)] for _ in 1:N])
    return PauliString(s)
end


@testset "strings operations" begin
    # operations between Operator and PauliString
    for N in (10, 70)
        O = rand_local2_M(N, 40)
        string = random_string(N)
        O2 = Operator(string)
        @test opnorm(O * string - O * O2) < 1E-10
        @test opnorm(commutator(O, string) - commutator(O, O2)) < 1E-10
        @test opnorm(anticommutator(O, string) - anticommutator(O, O2)) < 1E-10
        @test opnorm(commutator(O, string) + commutator(string, O)) < 1E-10
        @test opnorm(O + string - (O + O2)) < 1E-10
        @test opnorm(O - string - (O - O2)) < 1E-10
    end
    # operations between PauliStringTS and OperatorTS
    for N in (10, 70)
        O1 = compress(OperatorTS{(N,)}(rand_local2_M(N, 40)))
        string2 = PauliStringTS{(N,)}(random_string(N))
        O2 = OperatorTS(string2)
        @test opnorm(O1 * string2 - O1 * O2) < 1E-10
        @test opnorm(commutator(O1, string2) - commutator(O1, O2)) < 1E-10
        @test opnorm(anticommutator(O1, string2) - anticommutator(O1, O2)) < 1E-10
        @test opnorm(O1 + string2 - (O1 + O2)) < 1E-10
        @test opnorm(O1 - string2 - (O1 - O2)) < 1E-10
    end
    # operations between PauliStringTS and PauliStringTS
    for N in (10, 70)
        string1 = PauliStringTS{(N,)}(random_string(N))
        string2 = PauliStringTS{(N,)}(random_string(N))
        O1 = OperatorTS(string1)
        O2 = OperatorTS(string2)
        @test opnorm(string1 * string2 - (O1 * O2)) < 1E-10
        @test opnorm(string1 + string2 - (O1 + O2)) < 1E-10
        @test opnorm(string1 - string2 - (O1 - O2)) < 1E-10
        @test opnorm(commutator(string1, string2) - commutator(O1, O2)) < 1E-10
        @test opnorm(anticommutator(string1, string2) - anticommutator(O1, O2)) < 1E-10
    end
    # operations between PauliString and PauliString
    for N in (10, 70)
        string1 = random_string(N)
        string2 = random_string(N)
        O1 = Operator(string1)
        O2 = Operator(string2)
        @test opnorm(string1 + string2 - (O1 + O2)) < 1E-10
        @test opnorm(string1 - string2 - (O1 - O2)) < 1E-10
    end

end
