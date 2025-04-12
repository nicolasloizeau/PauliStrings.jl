
"""
    equivalence_class(A1::Operator, H::Operator)

Equivalence class of A1 with respect to Hamiltonian H, as defined in [https://arxiv.org/abs/2406.13026](https://arxiv.org/abs/2406.13026) Algorithm 1.
The equivalence class is returned as an operator with coeficients 1
Use [`get_pauli`](@ref) to access individual strings.
"""
function equivalence_class(A1::Operator, H::Operator)
    length(A1) != 1 && error("A1 needs to be an operator with a single Pauli string")
    A = deepcopy(A1)
    D = 1
    counter = 0
    while counter < D
        p = A.strings[counter+1]
        for i in 1:length(H)
            a, Ap = com(p, H.strings[i])
            if a != 0 && !(Ap in A.strings)
                push!(A.strings, Ap)
                push!(A.coeffs, (1im)^ycount(Ap.v, Ap.w))
                D += 1
            end
        end
        counter += 1
    end
    return A
end
