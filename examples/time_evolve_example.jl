
# example of time evolving Z_1 on a chaotic spin chain


using PauliStrings
import PauliStrings as ps
using PyPlot
using ProgressBars


# build a chaotic spin chain with periodic bc
function chaotic_chain(N::Int)
    H = ps.Operator(N)
    # XX interractions
    for j in 1:(N - 1)
        H += "X",j,"X",j+1
    end
    H += "X",1,"X",N # close the chain
    # fields
    for j in 1:N
        H += -1.05,"Z",j
        H += 0.5,"X",j
    end
    return H
end


N = 32 # system size
H = chaotic_chain(N) #hamiltonian
O = ps.Operator(N) #operator to time evolve
O += "Z", 1 # Z on site 1

ioff()


# heisenberg evolution of the operator O using rk4
# return tr(O(0)*O(t))/tr(O(t)^2)
# M is the number of strings to keep at each step
# noise is the amplitude of depolarizing noise
function evolve(H, O, M, times, noise)
    echo = []
    O0 = deepcopy(O)
    dt = times[2]-times[1]
    for t in ProgressBar(times)
        push!(echo, ps.trace(O*ps.dagger(O0))/ps.trace(O0*O0))
        #preform one step of rk4, keep only M strings, do not discard O0
        O = ps.rk4(H, O, dt; heisenberg=true, M=M,  keep=O0)
        #add depolarizingn oise
        O = ps.add_noise(O, noise*dt)
        # keep the M strings with the largest weight. Do not discard O0
        O = ps.trim(O, M; keep=O0)
    end
    return real.(echo)
end

plt.cla()
# time evolve O for different trim values
times = range(0, stop=5, step=0.05)
noise = 0.01
for trim in (10,12,14)
    S = evolve(H, O, 2^trim, times, noise)
    loglog(times, S) #plot S(t)
end

legend()
title("N=$N")
xlabel("t")
ylabel(L"tr$(Z_1(0)*Z_1(t))$")
savefig("time_evolve_example.png")
show()
