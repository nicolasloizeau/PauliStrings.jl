

using LinearAlgebra: I

@testset "fermions" begin
    N = 4
    for i in 1:N*2
        for j in 1:N*2
            m1 = majorana(N, i)
            m2 = majorana(N, j)
            prod = 2 * eye(N) * I[i, j]
            @test opnorm(anticommutator(m1, m2) - prod) == 0
        end
    end
end

@testset "construction" begin
    N = 3
    @test length(all_strings(N)) == 4^N
    N = 6
    k = 3
    @test length(all_k_local(N, k)) == binomial(N, k)*3^k
end
