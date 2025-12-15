

using LinearAlgebra: I

@testset "fermions" begin
    N = 4
    for i in 1:N*2
        for j in 1:N*2
            m1 = majorana(N, i)
            m2 = majorana(N, j)
            prod = 2 * eye(N) * I[i, j]
            @test norm(anticommutator(m1, m2) - prod) == 0
        end
    end
end

@testset "construction" begin
    N = 3
    @test length(all_strings(N)) == 4^N
    N = 6
    k = 3
    @test length(all_k_local(N, k)) == binomial(N, k) * 3^k
    L1 = 3
    L2 = 2
    a = string_2d(("X", 1, 1, "Y", 2, 1, "Z", 3, 1), L1, L2)
    b = string_2d(("Z", 3, 1, "Z", 1, 2, "X", 2, 2, "X", 3, 2), L1, L2)
    c = a * b
    @test string(a.strings[1]) == "XYZ111"
    @test string(b.strings[1]) == "11ZZXX"
    @test string(c.strings[1]) == "XY1ZXX"
end
