"""
 Creates a random 2-local operator with qubit q only supporting 'I' or 'Z'.
"""
function rand_local2_qubit_diag(N::Int, q::Int)
    o = Operator(N)
    for i in 1:N
        for j in 1:N
            if i != j
                for k in ['X', 'Y', 'Z']
                    k = (i == q && k != 'Z') ? rand(['I', 'Z']) : k
                    for l in ['X', 'Y', 'Z']
                        l = (j == q && l != 'Z') ? rand(['I', 'Z']) : l
                        o += (randn(rng, Float64), k, i, l, j)
                    end
                end
            end
        end
    end
    return compress(o)
end

@testset "truncation" begin
    N = 8
    O2 = ps.rand_local2(N)
    @test opnorm(ps.truncate(O2, 1)) == 0
    @test opnorm(ps.truncate(O2, 2)) == opnorm(O2)
    @test length(ps.trim(O2, 10)) == 10
    @test length(ps.prune(O2, 2)) <= length(ps.prune(O2, 20))
    @test length(ps.cutoff(XX(N) + 0.1 * X(N), 0.5)) == length(XX(N))
    @test opnorm(ps.cutoff(O2, 0.5)) <= opnorm(O2)
    @test opnorm(ps.add_noise(O2, 0.5)) < opnorm(O2)
    @test opnorm(ps.add_noise(O2, 0.5)) < opnorm(ps.add_noise(O2, 0.1))
    g = rand(N)
    @test opnorm(ps.add_dephasing_noise(O2, g)) < opnorm(O2)
    @test opnorm(ps.add_dephasing_noise(O2, g)) < opnorm(ps.add_dephasing_noise(O2, 0.2 * g))
    @test opnorm(Operator(ps.add_dephasing_noise(O2, fill(0.5, N)) - ps.add_dephasing_noise(O2, 0.5))) < 1e-10
    @test opnorm(Operator(ps.diag(ps.add_dephasing_noise(O2, g)) - ps.diag(O2))) < 1e-10
    q = rand(1:N) # qubit only supporting 'I' or 'Z'
    g = [i == q ? rand() : 0.0 for i in 1:N] # noise only affecting qubit q
    O2q = rand_local2_qubit_diag(N, q)
    @test opnorm(Operator(ps.add_dephasing_noise(O2q, g) - O2q)) < 1e-10
    @test opnorm(k_local_part(O2, 1) - ps.truncate(O2, 1)) == 0
    O1 = ising1D(N, 0.5)
    O1ts = OperatorTS1D(O1)
    @test opnorm(ps.truncate(O1ts, 2)) == opnorm(O1)
    @test length(ps.trim(O1ts, 1)) == 1
    @test opnorm(ps.cutoff(O1ts, 0.8)) <= opnorm(O1ts)
    @test opnorm(ps.add_noise(O1ts, 0.5)) < opnorm(O1ts)
    @test opnorm(Operator(ps.add_noise(O1ts, 0.5)) - ps.add_noise(O1, 0.5)) < 1e-10
end
