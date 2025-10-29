


@testset "strings operations" begin
    for N in (10, 70)
        O = rand_local2_M(N, 40)
        string = get_pauli(O, 5)
        O2 = Operator(string)
        @test opnorm(O*string-O*O2) < 1E-10
        @test opnorm(commutator(O, string)-commutator(O, O2)) < 1E-10
        @test opnorm(anticommutator(O, string)-anticommutator(O, O2)) < 1E-10
        @test opnorm(commutator(O, string)+commutator(string, O)) < 1E-10
        @test opnorm(O+string - (O+O2)) < 1E-10
        @test opnorm(O-string - (O-O2)) < 1E-10
    end


end
