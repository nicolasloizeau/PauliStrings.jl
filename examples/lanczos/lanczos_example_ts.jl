

# lanczos example from :
# https://journals.aps.org/prx/abstract/10.1103/PhysRevX.9.041017
# figure 2, "X in XX"

using PauliStrings
using PyPlot

# XX hamiltonian: 1d chain with XX+YY interraction
function XX(N)
    H = Operator(N)
    H += "X",1,"X",2
    H += "Z",1,"Z",2
    return OperatorTS{(N,)}(H)
end

# X local operator: X operator on each site
function X(N)
    H = Operator(N)
    H += "X",1
    return OperatorTS{(N,)}(H)
end


N = 50 # system size
H = XX(N) #hamiltonian
O = X(N) #operator

ioff()#pyplot

# nterms is the max pauli string length
for p in (14,16,18,20)
    @time bs = lanczos(H, O, 20, 2^p; keepnorm=true)
    plot(bs, label="trim: 2^$p")
end



legend()
ylabel(L"$b_n$")
xlabel(L"$n$")
title("X in XX, N=$N spins")
savefig("lanczos_example_ts.png")
show()
