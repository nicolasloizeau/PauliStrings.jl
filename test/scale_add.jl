using VectorInterface

@testset "scale (allocating) matches *" begin
    for N in (6, 12, 40)
        A = ps.rand_local2(N)
        for a in (2.0, -1.5, 0.5 + 0.25im, im)
            @test norm(scale(A, a) - a * A) < 1e-10
            @test norm(scale(A, a) - A * a) < 1e-10
        end
    end
end

@testset "scale! in place" begin
    for N in (6, 12, 40)
        A = ps.rand_local2(N)
        for a in (2.0, -1.5, 0.5 + 0.25im, im)
            o = copy(A)
            r = scale!(o, a)
            @test r === o                  # returns the same object
            @test norm(o - a * A) < 1e-10
        end
        # scale!(y, x, a) overwrites y with a*x
        y = ps.rand_local2(N)
        r = scale!(y, A, 3.0)
        @test r === y
        @test norm(y - 3.0 * A) < 1e-10
    end
end

@testset "scale!! maybe-in-place" begin
    A = ps.rand_local2(8)
    o = copy(A)
    r = VectorInterface.scale!!(o, 2.0)   # ComplexF64 stays ComplexF64 -> in place
    @test r === o
    @test norm(o - 2.0 * A) < 1e-10
end

@testset "zerovector" begin
    A = ps.rand_local2(8)
    z = VectorInterface.zerovector(A)
    @test length(z) == 0
    @test scalartype(z) == scalartype(A)
    z2 = VectorInterface.zerovector(A, ComplexF32)
    @test scalartype(z2) == ComplexF32
    o = copy(A)
    @test VectorInterface.zerovector!(o) === o
    @test length(o) == 0
end

@testset "add (allocating) matches +/-" begin
    for N in (6, 12, 40)
        A = ps.rand_local2(N)
        B = ps.rand_local2(N)
        @test norm(add(A, B) - (A + B)) < 1e-10            # α = β = 1
        @test norm(add(A, B, -1) - (A - B)) < 1e-10        # α = -1, β = 1
        # general β*A + α*B with real and complex weights
        for (α, β) in ((2.0, 3.0), (-0.5im, 1.0), (1.0, 0.25 + 0.1im))
            @test norm(add(A, B, α, β) - (β * A + α * B)) < 1e-10
        end
    end
end

@testset "add! / add!! in place match add" begin
    for N in (6, 12, 40)
        A = ps.rand_local2(N)
        B = ps.rand_local2(N)
        ref = add(A, B, 2.0, 3.0)
        y = copy(A)
        r = add!(y, B, 2.0, 3.0)
        @test r === y
        @test norm(y - ref) < 1e-10
        # add!! in place when the eltype is preserved (ComplexF64)
        y2 = copy(A)
        r2 = VectorInterface.add!!(y2, B, -1.0, 1.0)
        @test r2 === y2
        @test norm(y2 - (A - B)) < 1e-10
        # self-aliasing: add!(C, C, α, β) matches the non-aliased add(A, A, α, β)
        C = copy(A)
        add!(C, C, 0.7, -1.3)
        @test norm(C - add(A, A, 0.7, -1.3)) < 1e-10
        # exact cancellation collapses to an empty operator
        C = copy(A)
        @test length(add!(C, A, -1)) == 0
    end
end

@testset "scale / add on translation-symmetric operators" begin
    N = 6
    A = OperatorTS1D(ising1D(N, 1.0))
    B = OperatorTS1D(ising1D(N, 0.3))
    @test A isa OperatorTS
    @test norm(scale(A, 2.5) - 2.5 * A) < 1e-10
    @test norm(add(A, B, 0.7, -1.3) - (-1.3 * A + 0.7 * B)) < 1e-10
    y = copy(A)
    @test add!(y, B, 2.0) === y
    @test norm(y - (A + 2.0 * B)) < 1e-10
end

@testset "+/-/* delegate to add/scale" begin
    # exercise the delegated operators end-to-end against a hand-built reference
    A = ps.rand_local2(10)
    B = ps.rand_local2(10)
    @test norm((A + B) - add(A, B)) < 1e-12
    @test norm((A - B) - add(A, B, -1)) < 1e-12
    @test norm((2.5 * A) - scale(A, 2.5)) < 1e-12
    @test norm((A * 2.5) - scale(A, 2.5)) < 1e-12
end
