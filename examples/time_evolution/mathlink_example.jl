
# Compute a few Lanczos coefficients symbolocally using MathLink
# ---------------------------------------------------------------

using LinearAlgebra
using MathLink
using PauliStrings


function ising(h)
    H = OperatorMathLink(N)
    for i in 1:N
        H += h, "X", i
    end
    for i in 1:N
        H += "Z", i, "Z", mod1(i + 1, N)
    end
    return H
end

N = 10 #system size

# Initial operator O = X_1
O = OperatorMathLink(N)
O += 1, "X", 1

# Initialize a symbolic hamiltonian with a symbolic parameter h
h = W`h`
H = ising(h)


# Define assumptions for simplification
assumptions = W`Assumptions -> h > 0`

# Compute Lanczos coefficients
bn = lanczos(H, O, 5; assumptions=assumptions)
for b in bn
    println(b)
end
