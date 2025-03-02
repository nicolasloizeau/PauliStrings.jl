using BenchmarkTools
using Test
using PauliStrings
using Profile

# Lanczos benchmarks
# ------------------
# XX hamiltonian: 1d chain with XX + YY interaction
function XX(N)
    H = PauliStrings.Operator(N)
    for j in 1:(N-1)
        H += "X", j, "X", j + 1
        H += "Z", j, "Z", j + 1
    end
    return H
end

# X local operator: X operator on each site
function X(N)
    H = PauliStrings.Operator(N)
    for j in 1:N
        H += "X", j
    end
    return H
end


N = 40
H = XX(N)
O = X(N)
Hs = OperatorSorted(H)
Hu = OperatorUnsorted(H)
Os = OperatorSorted(O)
Ou = OperatorUnsorted(O)

steps = 20
p = 16
nterms = 2^p

@test OperatorSorted(O * O) == Os * Os
@test OperatorUnsorted(O * O) == Ou * Ou
@test OperatorSorted(com(O, O)) == com(Os, Os)
@test OperatorUnsorted(com(O, O)) == com(Ou, Ou)
@test OperatorSorted(H * H) == Hs * Hs
@test OperatorUnsorted(H * H) == Hu * Hu
@test OperatorSorted(com(H, H)) == com(Hs, Hs)
@test OperatorUnsorted(com(H, H)) == com(Hu, Hu)

result = PauliStrings.lanczos(H, O, steps, nterms; keepnorm=true)

# Hs = OperatorSorted(H)
# Os = OperatorSorted(O)
#
# result2 = PauliStrings.lanczos(Hs, Os, steps, nterms; keepnorm=true)

Hu = OperatorUnsorted(H)
Ou = OperatorUnsorted(O)

result3 = PauliStrings.lanczos(Hu, Ou, steps, nterms; keepnorm=true)

@assert !any(isnan, result3)

# b2 = @benchmark PauliStrings.lanczos($Hs, $Os, $steps, $nterms; keepnorm=true)
# display(b2)
# b3 = @benchmark PauliStrings.lanczos($Hu, $Ou, $steps, $nterms; keepnorm=true)
# b1 = @benchmark PauliStrings.lanczos($H, $O, $steps, $nterms; keepnorm=true)
# display(b1)
# display(b3)

Profile.clear()
@bprofile PauliStrings.lanczos($Hu, $Ou, $steps, $nterms; keepnorm=true)
VSCodeServer.view_profile()

