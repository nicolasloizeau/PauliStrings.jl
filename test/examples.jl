
using ProgressBars
using Random

rng = MersenneTwister(0)

function construction_example1()
    A = Operator(4)
    A += "X", 1
    A += "Y", 2
    A += 2, "Z", 3
    return A
end

function construction_example2()
    A = Operator(8)
    A += 2, "X", 1, "X", 2
    A += 1.5, "Y", 2, "Z", 3
    A += "Z", 3, "Z", 8
    return A
end

function construction_example3()
    A = Operator(6)
    A += "XYZ111"
    A += 5, "11XX11"
    return A
end

"""XX hamiltonian https://arxiv.org/abs/1812.08657"""
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

"""Y operator on all sites"""
function Y(N)
    H = Operator(N)
    for j in 1:N
        H += ('Y', j)
    end
    return H
end

"""Z operator on all sites"""
function Z(N)
    H = Operator(N)
    for j in 1:N
        H += ('Z', j)
    end
    return H
end

# build a chaotic spin chain with periodic bc
function chaotic_chain(N::Int)
    H = ps.Operator(N)
    # XX interractions
    for j in 1:(N-1)
        H += "X", j, "X", j + 1
    end
    H += "X", 1, "X", N # close the chain
    # fields
    for j in 1:N
        H += -1.05, "Z", j
        H += 0.5, "X", j
    end
    return H
end



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

function evolve_chaotic()
    N = 20 # system size
    H = chaotic_chain(N) #hamiltonian
    O = ps.Operator(N) #operator to time evolve
    O += "Z", 1 # Z on site 1
    times = range(0, stop=5, step=0.05)
    noise = 0.01
    S = evolve(H, O, 2^8, times, noise)
    return S[end]
end


"""
random 2-local operator with N sites and M terms for use in tests
"""
function rand_local2_M(N::Int, M::Int)
    o = Operator(N)
    for i in 1:M
        k = rand(['X', 'Y', 'Z'])
        l = rand(['X', 'Y', 'Z'])
        i = rand(1:N)
        j = i
        while j == i
            j = rand(1:N)
        end
        o += (randn(rng, Float64), k, i, l, j)
    end
    return compress(o)
end


"""
random 1-local operator with N sites and M terms for use in tests
"""
function rand_local1_M(N::Int, M::Int)
    o = Operator(N)
    for i in 1:M
        k = rand(['X', 'Y', 'Z'])
        i = rand(1:N)
        o += (randn(rng, Float64), k, i)
    end
    return compress(o)
end
