@testset "shift" begin
    @testset for Ls in [(5,), (2,7), (16,8), (4,3,1,2)]
        x = rand(UInt128)

        num_to_tensor(x) = reshape(reverse(digits(x;base=2,pad=prod(Ls))[1:prod(Ls)]),Ls)

        for shifts in Iterators.product(map(L->1:L, Ls)...)
            @test num_to_tensor(ps._shift(x, Ls, shifts)) == circshift(num_to_tensor(x), -1 .* shifts)
        end
    end
end

@testset "operatorts2d" begin
    @testset for L1 in 1:4
        L2 = 3
        O1 = ising2D(L1, L2, rand())
        O2 = 1im * ising2D(L1, L2, rand())
        O1ts = OperatorTS2D(O1, L1)
        O2ts = OperatorTS2D(O2, L1)
        eps = 1e-8
        @test typeof(O1ts + O2ts) == typeof(O1ts)
        @test typeof(O1ts - O2ts) == typeof(O1ts)
        @test typeof(O1ts * O2ts) == typeof(O1ts)
        @test typeof(O1ts - 1) == typeof(O1ts)
        @test typeof(O1ts + 1) == typeof(O1ts)
        @test typeof(O1ts * 2) == typeof(O1ts)
        @test typeof(O1ts / 2) == typeof(O1ts)
        @test opnorm(Operator(O1ts) - O1) < eps
        @test opnorm(Operator(O1ts * O2ts) - O1 * O2) < eps
        @test opnorm(Operator(O1ts * O2ts * O2ts) - O1 * O2 * O2) < eps
        @test opnorm(Operator(O1ts + O2ts) - (O1 + O2)) < eps
        @test opnorm(Operator(O1ts * 7) - O1 * 7) < eps
        @test opnorm(Operator(O1ts + 7) - (O1 + 7)) < eps
        @test is_ts2d(O1,L1)
        @test is_ts2d(O1*O2,L1)
        @test !is_ts2d(ps.rand_local2(L1*L2),L1)
        @test trace(O1ts) == trace(dagger(O1))
        @test isapprox(trace(O1ts^2), trace(O1^2))
        @test trace(O1ts) == trace(diag(O1ts))
        @test isapprox(opnorm(O1ts*O2ts), opnorm(O1*O2))
        k = 1 + L1
        @test isapprox(trace_product(O1ts, k),trace_product(O1, k))
        @test isapprox(trace_product(O1ts, k), trace(O1^k))
        @test isapprox(trace_product(O1ts, k), trace(O1ts^k))
    end
end
