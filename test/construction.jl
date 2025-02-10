

using LinearAlgebra: I

@testset "fermions" begin
    N = 4
    for i in 1:N*2
        for j in 1:N*2
            m1 = majorana(N, i)
            m2 = majorana(N, j)
            prod = 2 * eye(N) * I[i, j]
            @test opnorm(com(m1, m2; anti=true) - prod) == 0
        end
    end
end
