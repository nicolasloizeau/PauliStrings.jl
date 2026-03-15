
"""
    add_noise(o::Operator, g::Real)

Add depolarizing noise that make the long string decays. `g` is the noise amplitude.
Each string is multiplied by exp(-g * w), where w is the number of non-identity Pauli operators in the string.
This is equivalent to the compostion of ``N`` single qubit depolarizing channels with Kraus operators
``\\sqrt{1-\\frac{3p}{4}} I_i, \\, \\sqrt{\\frac{p}{4}} X_i, \\, \\sqrt{\\frac{p}{4}} Y_i, \\, \\sqrt{\\frac{p}{4}} Z_i``
and ``p=1-e^{-g}``.

# Example
```julia
A = add_noise(A, 0.1)
```
# Reference
[https://arxiv.org/pdf/2407.12768](https://arxiv.org/pdf/2407.12768)
"""
function add_noise(o::AbstractOperator, g::Real)
    o2 = deepcopy(o)
    for i in 1:length(o)
        o2.coeffs[i] *= exp(-pauli_weight(o.strings[i]) * g)
    end
    return o2
end


"""
    add_noise(o::Operator, g::AbstractVector{<:Real})

Add local depolarizing noise.

If \$g_j\$ is the noise amplitude for site \$j\$, then each string will be multiplied by
\$e^{-\\sum_j g_j}\$, where the sum runs over the sites with non-unit Pauli operators.

"""
function add_noise(o::AbstractOperator, g::AbstractVector{<:Real})
    o2 = deepcopy(o)
    N = qubitlength(o)
    N != length(g) && throw(ArgumentError("length of g ($(length(g))) must be $N"))
    for i in 1:length(o)
        p = o.strings[i]
        noise_bits = p.v | p.w
        o2.coeffs[i] *= exp(-sum(Real[g[j] for j in 1:N if ((noise_bits >> (j - 1)) & 1) == 1]))
    end
    return o2
end



dephasing_weight_z(p::PauliString) = count_ones(p.w)
dephasing_weight_x(p::PauliString) = count_ones(p.v)
dephasing_weight_y(p::PauliString) = count_ones(p.v ⊻ p.w)

dephasing_weights = Dict(:Z => dephasing_weight_z,
                         :X => dephasing_weight_x,
                         :Y => dephasing_weight_y)


"""
    add_dephasing_noise(o::AbstractOperator, g::Real; basis::Symbol=:Z)

Add dephasing noise.

If ``g`` is the noise amplitude, then each string will decay by a factor of
``e^{-gw}``, where ``w`` is the count of Pauli operators in the string that are
either ``X`` or ``Y``.

# Reference
[https://arxiv.org/abs/2306.05804](https://arxiv.org/pdf/2306.05804)
"""
function add_dephasing_noise(o::AbstractOperator, g::Real; basis::Symbol=:Z)
    dephasing_weight = dephasing_weights[basis]
    o2 = deepcopy(o)
    for i in 1:length(o)
        p = o.strings[i]
        o2.coeffs[i] *= exp(-g * dephasing_weight(p)) # p.w gives the bits containing 'X' or 'Y' operators
    end
    return o2
end



"""
    add_dephasing_noise(o::AbstractOperator, g::AbstractVector{<:Real})

Add local dephasing noise.

If ``g_j`` is the noise amplitude of site ``j``, then each string will be multiplied
by ``e^{-\\sum_j g_j}``, where the sum runs over the sites with Pauli operators that
are either ``X`` or ``Y``.

"""
function add_dephasing_noise(o::AbstractOperator, g::AbstractVector{<:Real})
    o2 = deepcopy(o)
    N = qubitlength(o)
    N != length(g) && throw(ArgumentError("length of g ($(length(g))) must be $N"))
    for i in 1:length(o)
        p = o.strings[i]
        noise_bits = p.w # bits containing 'X' or 'Y' operators.
        o2.coeffs[i] *= exp(-sum(Real[g[j] for j in 1:N if ((noise_bits >> (j - 1)) & 1) == 1]))
    end
    return o2
end
