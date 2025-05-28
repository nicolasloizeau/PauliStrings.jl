
@testset "truncation" begin

    function ising1D(N, g)
        H = Operator(N)
        for j in 1:(N-1)
            H += "Z", j, "Z", j + 1
        end
        H += "Z", 1, "Z", N #periodic boundary condition
        for j in 1:N
            H += g, "X", j
        end
        return -H
    end

    function XX(N)
        H = Operator(N)
        for j in 1:(N-1)
            H += ('X', j, 'X', j + 1)
            H += ('Z', j, 'Z', j + 1)
        end
        H += ('X', N, 'X', 1)
        H += ('Z', N, 'Z', 1)
        return H
    end

    """X operator on all sites"""
    function X(N)
        H = Operator(N)
        for j in 1:N
            H += ('X', j)
        end
        return H
    end

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
    @test opnorm(k_local_part(O2, 1) - ps.truncate(O2, 1)) == 0
    O1 = ising1D(N, 0.5)
    O1ts = OperatorTS1D(O1)
    @test opnorm(ps.truncate(O1ts, 2)) == opnorm(O1)
    @test length(ps.trim(O1ts, 1)) == 1
    @test opnorm(ps.cutoff(O1ts, 0.8)) <= opnorm(O1ts)
    @test opnorm(ps.add_noise(O1ts, 0.5)) < opnorm(O1ts)
    @test opnorm(Operator(ps.add_noise(O1ts, 0.5)) - ps.add_noise(O1, 0.5)) < 1e-10
    @test opnorm(ps.add_dephasing_noise(O2, zeros(N)) - O2) < 1e-12
    dephased = ps.add_dephasing_noise(O2, fill(0.5, N))
    @test opnorm(dephased) < opnorm(O2)
    op_low = opnorm(ps.add_dephasing_noise(O2, fill(0.1, N)))
    op_high = opnorm(ps.add_dephasing_noise(O2, fill(1.0, N)))
    @test op_high < op_low
    # Test Z-only strings are unaffected
    O_z = Operator(N)
    O_z += 1.0, "Z"^N
    O_z_dephased = ps.add_dephasing_noise(O_z, fill(0.9, N))
    @test opnorm(O_z_dephased - O_z) < 1e-12
    I_op = Operator(N)
    I_op += 1.0, "I"^N
    I_dephased = ps.add_dephasing_noise(I_op, fill(1.0, N))
    @test opnorm(I_dephased - I_op) < 1e-12
    # X-only string should be exponentially damped
    N = 8
    O_x = Operator(N)
    O_x += 1.0, "X"^N
    O_x_damped = ps.add_dephasing_noise(O_x, fill(1.0, N))
    coeff = O_x_damped.coeffs[1]
    expected = exp(-N)
    @test isapprox(abs(coeff), expected, atol=1e-12)
    # Random noise vector
    gammas = rand(N)
    O_rand = ps.rand_local2(N)
    O_rand_damped = ps.add_dephasing_noise(O_rand, gammas)
    @test opnorm(O_rand_damped) < opnorm(O_rand)
    # Manual check: one coefficient with X and Y should be damped
    O_test = Operator(3)
    O_test += 1.0, "XYZ"
    O_test += 2.0, "ZIZ"
    g = [0.1, 0.2, 0.3]
    O_damped = ps.add_dephasing_noise(O_test, g)
    # locate by PauliString
    ziz_str = ps.PauliString("ZIZ")
    xyz_str = ps.PauliString("XYZ")
    i_ziz = findfirst(==(ziz_str), O_damped.strings)
    i_xyz = findfirst(==(xyz_str), O_damped.strings)
    # Z-only term must be exactly unchanged
    @test isapprox(O_damped.coeffs[i_ziz], 2.0 + 0im, atol=1e-12)
    # XYZ term must be damped by exp(-g1 - g2)
    expected = exp(-g[1]) * exp(-g[2])
    @test isapprox(abs(O_damped.coeffs[i_xyz]), expected, atol=1e-10)
end
