# reference implementation: L_ij = i tr(P_i [H, P_j]) / 2^N
function liouvillian_reference(H, basis)
    N = qubitlength(H)
    L = zeros(length(basis), length(basis))
    for (j, Pj) in enumerate(basis), (i, Pi) in enumerate(basis)
        L[i, j] = real(im * trace_product(Operator(Pi), commutator(H, Operator(Pj)))) / 2^N
    end
    return L
end

@testset "liouvillian" begin
    N = 3
    H = rand_local2(N)
    H = H + H'
    basis = complete_basis(N)
    L = liouvillian(H, basis)
    @test size(L) == (4^N, 4^N)
    @test L ≈ liouvillian_reference(H, basis)
    @test L ≈ liouvillian(H)
    @test L ≈ -L'
    # dc/dt = L c reproduces i[H, O]
    O = rand_local2(N)
    c = [trace_product(Operator(P), O) / 2^N for P in basis]
    dc = L * c
    dO = im * commutator(H, O)
    @test dc ≈ [trace_product(Operator(P), dO) / 2^N for P in basis]
    # sub-basis
    sub = k_local_basis(N, 1)
    Lsub = liouvillian(H, sub)
    @test size(Lsub) == (length(sub), length(sub))
    @test Lsub ≈ liouvillian_reference(H, sub)
end

@testset "depolarizing_noise" begin
    N = 3
    γ = 0.3
    basis = complete_basis(N)
    D = depolarizing_noise(N, γ)
    @test D ≈ depolarizing_noise(N, γ, basis)
    @test size(D) == (4^N, 4^N)
    @test isdiag(Matrix(D))
    @test all(diag(D) .== [-γ * pauli_weight(p) for p in basis])
    # exp(D t) reproduces add_noise with amplitude γ t
    O = rand_local2(N)
    t = 0.7
    c = [trace_product(Operator(P), O) / 2^N for P in basis]
    c_noisy = exp(Matrix(D) * t) * c
    @test c_noisy ≈ [trace_product(Operator(P), add_noise(O, γ * t)) / 2^N for P in basis]
    sub = k_local_basis(N, 1)
    @test diag(depolarizing_noise(N, γ, sub)) == fill(-γ, length(sub))
    @test_throws ArgumentError depolarizing_noise(N + 1, γ, sub)
end

@testset "depolarizing_noise_channel" begin
    N = 3
    g = 0.3
    basis = complete_basis(N)
    E = depolarizing_noise_channel(N, g)
    @test E ≈ depolarizing_noise_channel(N, g, basis)
    @test isdiag(Matrix(E))
    @test Matrix(E) ≈ exp(Matrix(depolarizing_noise(N, g)))
    O = rand_local2(N)
    c = [trace_product(Operator(P), O) / 2^N for P in basis]
    @test E * c ≈ [trace_product(Operator(P), add_noise(O, g)) / 2^N for P in basis]
    sub = k_local_basis(N, 1)
    @test diag(depolarizing_noise_channel(N, g, sub)) ≈ fill(exp(-g), length(sub))
    @test_throws ArgumentError depolarizing_noise_channel(N + 1, g, sub)
end
