# Tests for the bucketed, multithreaded binary kernel.
#
# These exercise more parallelism when run with multiple threads
# (`julia -t auto`), but they are correct on a single thread too: every
# bucketing strategy must reproduce the `Serial` reference kernel exactly.

@testset "bucketing" begin
    relerr(A, B) = norm(A - B) / max(norm(B), 1)

    # The translation-invariance trait drives all kernel dispatch.
    @testset "is_translation_invariant trait" begin
        @test ps.is_translation_invariant(ps.PauliString{6}) == false
        @test ps.is_translation_invariant(paulistringtype(12)) == false
        @test ps.is_translation_invariant(ps.periodicpaulistringtype((6,))) == true
        @test ps.is_translation_invariant(rand_local2_M(8, 20)) == false
        Hts = OperatorTS1D(rand_local2_M(8, 20); full = false)
        @test ps.is_translation_invariant(Hts) == true
    end

    # (N, #terms) covering UInt8, UInt32, UInt64 backing integers
    cases = ((6, 60), (20, 200), (40, 300))
    ops = (("prod", ps.prod), ("commutator", ps.commutator), ("anticommutator", ps.anticommutator))

    @testset "N=$N" for (N, M) in cases
        P = paulistringtype(N)
        Bs = filter(b -> b < min(N, 8 * sizeof(ps.uinttype(N))), (3, 5, 7))
        strategies = AbstractBucketStrategy[mixing_matrix(N, b) for b in Bs]

        @testset "GF(2)-linearity" begin
            for s in strategies
                @test ps.is_xor_linear(s, P; nsamples = 2000)
            end
        end

        # mixing_matrix is full-rank by construction: every one of its 2^b buckets
        # is reachable, so no bucket is structurally always-empty.
        @testset "mixing_matrix full rank" begin
            for b in Bs
                s = mixing_matrix(N, b)
                seen = Set(ps.bucketindex(s, ps.PauliString{N}(rand(ps.uinttype(N)), rand(ps.uinttype(N)))) for _ in 1:(1 << (b + 6)))
                @test length(seen) == 1 << b
            end
        end

        A = rand_local2_M(N, M)
        B = rand_local2_M(N, M)
        C0 = rand_local1_M(N, M ÷ 4)

        @testset "$name equals Serial" for (name, f) in ops
            ref = ps.binary_kernel(f, A, B; strategy = Serial())
            for s in strategies
                @test relerr(ps.binary_kernel(f, A, B; strategy = s), ref) < 1e-10
            end

            # α scaling + β·C accumulation in place
            α, β = 0.5 + 0.2im, 2.0
            refC = copy(C0)
            ps.binary_kernel!(f, refC, A, B, α, β; strategy = Serial())
            for s in strategies
                Cb = copy(C0)
                ps.binary_kernel!(f, Cb, A, B, α, β; strategy = s)
                @test relerr(Cb, refC) < 1e-10
            end

            # maxlength truncation and epsilon cutoff
            for s in strategies
                @test relerr(
                    ps.binary_kernel(f, A, B; strategy = s, maxlength = 4),
                    ps.binary_kernel(f, A, B; strategy = Serial(), maxlength = 4),
                ) < 1e-10
                @test relerr(
                    ps.binary_kernel(f, A, B; strategy = s, epsilon = 1e-2),
                    ps.binary_kernel(f, A, B; strategy = Serial(), epsilon = 1e-2),
                ) < 1e-10
            end
        end

        # Disjointness/routing: a product lands in the bucket fixed by linearity,
        # bucket(p₁ · p₂) == bucket(p₁) ⊻ bucket(p₂), so output buckets never overlap.
        @testset "output routing" begin
            ksA, ksB = keys(A), keys(B)
            for s in strategies
                ok = true
                for p₁ in ksA, p₂ in ksB
                    p, _ = ps.prod(p₁, p₂)
                    ok &= ps.bucketindex(s, p) == (ps.bucketindex(s, p₁) ⊻ ps.bucketindex(s, p₂))
                end
                @test ok
            end
        end
    end

    # The bucketized operator exposes the keys/values interface, per bucket and
    # as a whole, and partitions the operator exactly.
    @testset "Bucket keys/values interface" begin
        N = 20
        A = rand_local2_M(N, 200)
        for s in (mixing_matrix(N, 6),)
            bk = ps.bucketize(s, A)
            @test length(bk) == 1 << ps.nbits(s)
            @test sum(length, bk) == length(A)
            @test length(keys(bk)) == length(A)        # full bucket-ordered arrays
            @test length(values(bk)) == length(A)
            for (i, bucket) in enumerate(bk)            # iterable as a vector of buckets
                @test keys(bucket) === bucket.strings
                @test values(bucket) === bucket.coeffs
                @test length(bucket) == length(keys(bucket)) == length(values(bucket))
                # every string in bucket i-1 hashes to i-1
                @test all(p -> ps.bucketindex(s, p) == i - 1, keys(bucket))
                @test length(collect(pairs(bucket))) == length(bucket)
            end
        end
    end

    # The public operators dispatch through the adaptive default and must agree
    # with the serial kernel regardless of how `default_strategy` routes them.
    @testset "default dispatch" begin
        N = 24
        A = rand_local2_M(N, 400)
        B = rand_local2_M(N, 400)
        @test ps.default_strategy(A, B) isa LinearMatrix    # plain PauliString: always bucketed
        @test norm(A * B - ps.binary_kernel(ps.prod, A, B; strategy = Serial())) / norm(A * B) < 1e-10
        @test norm(commutator(A, B) - ps.binary_kernel(ps.commutator, A, B; strategy = Serial())) / max(norm(commutator(A, B)), 1) < 1e-10
        @test norm(anticommutator(A, B) - ps.binary_kernel(ps.anticommutator, A, B; strategy = Serial())) / max(norm(anticommutator(A, B)), 1) < 1e-10
    end

    # Translation-invariant operators are not bucketable here: default falls back
    # to Serial, and an explicit bucketing strategy is rejected.
    @testset "translation-invariant fallback" begin
        N = 12
        Hts = OperatorTS1D(rand_local2_M(N, 40); full = false)
        @test ps.default_strategy(Hts, Hts) isa Serial
        @test Hts * Hts isa Operator           # works via the serial path
        @test_throws ArgumentError ps.binary_kernel(ps.prod, Hts, Hts; strategy = mixing_matrix(N, 5))
    end
end
