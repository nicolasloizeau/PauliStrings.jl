using PauliStrings

# heisenberg evolution of the operator O using rk4
# return tr(O(0)*O(t))/tr(O(t)^2)
# M is the number of strings to keep at each step
# noise is the amplitude of depolarizing noise
function evolve(H, O, M, times, noise)
    echo = []
    O0 = deepcopy(O)
    dt = times[2] - times[1]
    for t in ProgressBar(times)
        push!(echo, ps.trace(O * ps.dagger(O0)) / ps.trace(O0 * O0))
        #preform one step of rk4, keep only M strings, do not discard O0
        O = ps.rk4(H, O, dt; heisenberg=true, M=M, keep=O0)
        #add depolarizingn oise
        O = ps.add_noise(O, noise * dt)
        # keep the M strings with the largest weight. Do not discard O0
        O = ps.trim(O, M; keep=O0)
    end
    return real.(echo)
end



function test_evolve_chaotic()
    N = 20 # system size
    H = chaotic_chain(N) #hamiltonian
    O = ps.Operator(N) #operator to time evolve
    O += "Z", 1 # Z on site 1
    times = range(0, stop=5, step=0.05)
    noise = 0.01
    S = evolve(H, O, 2^8, times, noise)
    return S[end]
end




function evolve_lindblad(H, L, O, M, times)
    echo = []
    O0 = deepcopy(O)
    dt = times[2] - times[1]
    for t in ProgressBar(times)
        push!(echo, trace(O * O0) / 2^qubitlength(O))
        O = rk4_lindblad(H, O, dt, L; heisenberg=true, M=M)
        O = trim(O, M)
    end
    return real.(echo)
end

function test_evolve_lindblad_xx()
    M = 2^12
    N = 4
    H = XX(N) #hamiltonian
    O = X(N) #operator to time evolve
    L = [Zi(N, i) * 0.5 for i in 1:N] # Lindblad jump operators
    times = 0:0.01:1
    S = evolve_lindblad(H, L, O, M, times)
    return S[end]
end


@testset "time_evolution" begin
    @test_broken norm(test_evolve_chaotic() - 0.09208935978092929) < 1e-10
    @test norm(test_evolve_lindblad_xx() + 0.2493404959734003) < 1e-10
end
