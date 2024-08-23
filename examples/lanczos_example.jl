

# lanczos example from :
# https://journals.aps.org/prx/abstract/10.1103/PhysRevX.9.041017
# figure 2, "X in XX"

using PauliStrings
import PauliStrings as ps
using PyPlot

# XX hamiltonian: 1d chain with XX+YY interraction
function XX(N)
    H = ps.Operator(N)
    for j in 1:(N - 1)
        H += "X",j,"X",j+1
        H += "Z",j,"Z",j+1
    end
    return H
end

# X local operator: X operator on each site
function X(N)
    H = ps.Operator(N)
    for j in 1:N
        H += "X",j
    end
    return H
end


N = 50 # system size
H = XX(N) #hamiltonian
O = X(N) #operator
O = O/ps.norm_lanczos(O)

ioff()#pyplot

# nterms is the max pauli string length
for p in (14,16,18,20)
    @time bs = ps.lanczos(H, O, 20, 2^p; keepnorm=true)
    plot(bs, label="trim: 2^$p")
end



legend()
ylabel(L"$b_n$")
xlabel(L"$n$")
title("X in XX, N=$N spins")
savefig("lanczos_example.png")
show()
