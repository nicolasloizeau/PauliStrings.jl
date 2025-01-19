
"""
    equivalence_class(A1::Union{Operator64,Operator128}, H::Union{Operator64,Operator128})

Equivalence class of A1 with respect to Hamiltonian H, as defined in [https://arxiv.org/abs/2406.13026](https://arxiv.org/abs/2406.13026) Algorithm 1.
The equivalence class is returned as an operator with coeficients 1
Use [`get_pauli`](@ref) to access individual strings.
"""
function equivalence_class(A1::Union{Operator64,Operator128}, H::Union{Operator64,Operator128})
    length(A1) != 1 && error("A1 needs to be an operator with a single Pauli string")
    A = deepcopy(A1)
    D = 1
    counter = 0
    while counter < D
        Vv = A.v[counter+1]
        Vw = A.w[counter+1]
        for i in 1:length(H)
            a, Av, Aw = com(Vv, Vw, H.v[i], H.w[i])
            if a != 0 && !vw_in_o(Av, Aw, A)
                push!(A, (1im)^ycount(Av, Aw), Av, Aw)
                D += 1
            end
        end
        counter += 1
    end
    return A
end
