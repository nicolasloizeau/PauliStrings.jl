


f_unitary(H, O, s, hbar) = s * commutator(H, O) / hbar


"""
    rk4(H::AbstractOperator, O::AbstractOperator, dt::Real; hbar::Real=1, heisenberg=true, truncation::Function=copy, f=f_unitary)

Single step of Runge–Kutta-4 with time independant Hamiltonian.
Returns O(t+dt).
`truncation` : function that takes an operator and returns a truncated version of it. By default it is `identity` (no truncation).
"""
function rk4(H::AbstractOperator, O::AbstractOperator, dt::Real;
             hbar::Real=1, heisenberg=true, truncation::Function=identity, f=f_unitary)
    s = heisenberg ? 1im : -1im
    k1 = truncation(f(H, O, s, hbar))
    k2 = truncation(f(H, O + dt * k1 / 2, s, hbar))
    k3 = truncation(f(H, O + dt * k2 / 2, s, hbar))
    k4 = truncation(f(H, O + dt * k3, s, hbar))
    return O + (k1 + 2 * k2 + 2 * k3 + k4) * dt / 6
end



"""
    dopri5(H::AbstractOperator, O::AbstractOperator, dt::Real; hbar::Real=1, heisenberg=true, truncation::Function=identity, f=f_unitary)

Single step of Dormand–Prince-5 with time independant Hamiltonian.
Returns O(t+dt).
`truncation` : function that takes an operator and returns a truncated version of it. By default it is `identity` (no truncation).
"""
function dopri5(H::AbstractOperator, O::AbstractOperator, dt::Real;
                hbar::Real=1, heisenberg=true, truncation::Function=identity, f=f_unitary)
    s = heisenberg ? 1im : -1im
    a21 = 1/5
    a31, a32 = 3/40, 9/40
    a41, a42, a43 = 44/45, -56/15, 32/9
    a51, a52, a53, a54 = 19372/6561, -25360/2187, 64448/6561, -212/729
    a61, a62, a63, a64, a65 = 9017/3168, -355/33, 46732/5247, 49/176, -5103/18656
    b1, b3, b4, b5, b6 = 35/384, 500/1113, 125/192, -2187/6784, 11/84
    k1 = truncation(f(H, O, s, hbar))
    k2 = truncation(f(H, O + dt*a21*k1, s, hbar))
    k3 = truncation(f(H, O + dt*(a31*k1 + a32*k2), s, hbar))
    k4 = truncation(f(H, O + dt*(a41*k1 + a42*k2 + a43*k3), s, hbar))
    k5 = truncation(f(H, O + dt*(a51*k1 + a52*k2 + a53*k3 + a54*k4), s, hbar))
    k6 = truncation(f(H, O + dt*(a61*k1 + a62*k2 + a63*k3 + a64*k4 + a65*k5), s, hbar))
    return O + dt*(b1*k1 + b3*k3 + b4*k4 + b5*k5 + b6*k6)
end



"""
    rk4(H::Function, O::AbstractOperator, dt::Real, t::Real; hbar::Real=1, heisenberg=true, truncation::Function=identity, f=f_unitary)

Single step of Runge–Kutta-4 with time dependant Hamiltonian.
Returns O(t+dt).
`H` is a function that takes time `t` and returns the Hamiltonian at that time.
`truncation` : function that takes an operator and returns a truncated version of it. By default it is `identity` (no truncation).
"""
function rk4(H::Function, O::AbstractOperator, dt::Real, t::Real; hbar::Real=1, heisenberg=true, truncation::Function=identity, f=f_unitary)
    return rk4(H(t), O, dt; hbar=hbar, heisenberg=heisenberg, truncation=truncation, f=f)
end



function f_lindblad(H, O, s, hbar, L, gamma)
    Odot = f_unitary(H, O, s, hbar)
    for i in 1:length(L)
        if s == 1im
            A = L[i]' * O * L[i]
        elseif s == -1im
            A = L[i] * O * L[i]'
        else
            throw(ArgumentError("s should be 1im or -1im"))
        end
        LL = L[i]' * L[i]
        Odot += gamma[i] * (A - anticommutator(LL, O) / 2)
    end
    return Odot
end




"""
    rk4_lindblad(H::AbstractOperator, O::AbstractOperator, dt::Real, L; hbar::Real=1, heisenberg=true, truncation::Function=copy, gamma=[]))

Single step of Runge–Kutta-4 for solving the Lindblad equation

``\\dot{O}=i[H,O]+\\sum_i \\gamma_i \\left(L_i^\\dagger O L_i -\\frac{1}{2} \\{ L_i^\\dagger L_i, O\\} \\right)``

Returns O(t+dt).
`L` is a list of jump operators.
"""
function rk4_lindblad(H::AbstractOperator, O::AbstractOperator, dt::Real, L; hbar::Real=1, heisenberg=true, truncation::Function=copy, gamma=[])
    @assert length(gamma) == length(L) || length(gamma) == 0
    (length(gamma) == 0) && (gamma = ones(length(L)))
    f(H, O, s, hbar) = f_lindblad(H, O, s, hbar, L, gamma)
    return rk4(H, O, dt; hbar=hbar, heisenberg=heisenberg, truncation=truncation, f=f)
end



# ---------------------------------------------------------------------------
# In-place variants
#
# These mirror the functions above but lean on the in-place kernels
# `commutator!`/`anticommutator!` and `scale!` to avoid intermediate
# allocations. They overwrite the operator `O` with `O(t+dt)` and return it.
# ---------------------------------------------------------------------------

# In-place time-derivative interface (mirrors `f`/`f_unitary`, but mutating).
#
#   Ȯ = (s/ħ)·[H, O]                         (unitary, s = ±i)
#
# `f!(C, H, O, s, hbar)` writes the derivative of `O` into `C` and returns it.
# `C` may alias `O`: `commutator!` reads `H` and `O` in full before it
# resizes/overwrites `C`, so a stage buffer `Y = O + c·dt·k` can be turned into
# its own derivative in place via `f!(Y, H, Y, s, hbar)`.
f_unitary!(C::AbstractOperator, H, O, s, hbar) = commutator!(C, H, O, s / hbar)

# Overwrite `O`'s storage with the contents of `R`, returning `O`.
function _store!(O::AbstractOperator, R::AbstractOperator)
    resize!(O, length(R))
    copyto!(keys(O), keys(R))
    copyto!(values(O), values(R))
    return O
end


"""
    rk4!(O::AbstractOperator, H::AbstractOperator, dt::Real; hbar::Real=1, heisenberg=true, truncation::Function=identity, f! = f_unitary!)

In-place single step of Runge–Kutta-4 with a time independant Hamiltonian.
Overwrites `O` with `O(t+dt)` and returns it, reusing the in-place
`commutator!`/`scale!` kernels to avoid intermediate allocations.
`truncation` : function that takes an operator and returns a truncated version of it. By default it is `identity` (no truncation).
See [`rk4`](@ref) for the allocating version.
"""
function rk4!(O::AbstractOperator, H::AbstractOperator, dt::Real;
              hbar::Real=1, heisenberg=true, truncation::Function=identity, f! = f_unitary!)
    s = heisenberg ? 1im : -1im
    # RK4 (Ȯ = f!(·); default Ȯ = (s/ħ)·[H, O]):
    #   k₁ = Ȯ(O)
    #   k₂ = Ȯ(O + dt/2·k₁)
    #   k₃ = Ȯ(O + dt/2·k₂)
    #   k₄ = Ȯ(O + dt·k₃)
    #   O(t+dt) = O + dt/6·(k₁ + 2k₂ + 2k₃ + k₄)

    # k₁: derivative of O into a fresh buffer (must not clobber O, still needed below)
    k₁ = truncation(f!(zero(typeof(O)), H, O, s, hbar))

    # k₂…k₄: build the stage argument Y (fresh `+`), then overwrite it with its
    # own derivative in place — `f!(Y, H, Y, …)` is safe because C may alias O.
    k₂ = O + (dt / 2) * k₁
    k₂ = truncation(f!(k₂, H, k₂, s, hbar))
    k₃ = O + (dt / 2) * k₂
    k₃ = truncation(f!(k₃, H, k₃, s, hbar))
    k₄ = O + dt * k₃
    k₄ = truncation(f!(k₄, H, k₄, s, hbar))

    # Δ = dt/6·(k₁ + 2k₂ + 2k₃ + k₄); scale the now-dead kᵢ in place with scale!
    Δ = k₁ + scale!(k₂, 2) + scale!(k₃, 2) + k₄
    scale!(Δ, dt / 6)

    # O(t+dt) = O + Δ, overwriting O's storage in place
    return _store!(O, O + Δ)
end


"""
    dopri5!(O::AbstractOperator, H::AbstractOperator, dt::Real; hbar::Real=1, heisenberg=true, truncation::Function=identity, f! = f_unitary!)

In-place single step of Dormand–Prince-5 with a time independant Hamiltonian.
Overwrites `O` with `O(t+dt)` and returns it.
`truncation` : function that takes an operator and returns a truncated version of it. By default it is `identity` (no truncation).
See [`dopri5`](@ref) for the allocating version.
"""
function dopri5!(O::AbstractOperator, H::AbstractOperator, dt::Real;
                 hbar::Real=1, heisenberg=true, truncation::Function=identity, f! = f_unitary!)
    s = heisenberg ? 1im : -1im
    a21 = 1/5
    a31, a32 = 3/40, 9/40
    a41, a42, a43 = 44/45, -56/15, 32/9
    a51, a52, a53, a54 = 19372/6561, -25360/2187, 64448/6561, -212/729
    a61, a62, a63, a64, a65 = 9017/3168, -355/33, 46732/5247, 49/176, -5103/18656
    b1, b3, b4, b5, b6 = 35/384, 500/1113, 125/192, -2187/6784, 11/84

    # Each stage argument is a fresh operator; turn it into its own derivative in
    # place with f!(kᵢ, H, kᵢ, s, hbar). k₁ uses a fresh buffer to spare `O`.
    k1 = truncation(f!(zero(typeof(O)), H, O, s, hbar))
    k2 = O + dt*a21*k1
    k2 = truncation(f!(k2, H, k2, s, hbar))
    k3 = O + dt*(a31*k1 + a32*k2)
    k3 = truncation(f!(k3, H, k3, s, hbar))
    k4 = O + dt*(a41*k1 + a42*k2 + a43*k3)
    k4 = truncation(f!(k4, H, k4, s, hbar))
    k5 = O + dt*(a51*k1 + a52*k2 + a53*k3 + a54*k4)
    k5 = truncation(f!(k5, H, k5, s, hbar))
    k6 = O + dt*(a61*k1 + a62*k2 + a63*k3 + a64*k4 + a65*k5)
    k6 = truncation(f!(k6, H, k6, s, hbar))

    # O(t+dt) = O + dt·(b1·k₁ + b3·k₃ + b4·k₄ + b5·k₅ + b6·k₆)
    Δ = scale!(k1, b1) + scale!(k3, b3) + scale!(k4, b4) + scale!(k5, b5) + scale!(k6, b6)
    scale!(Δ, dt)
    return _store!(O, O + Δ)
end


"""
    rk4!(O::AbstractOperator, H::Function, dt::Real, t::Real; hbar::Real=1, heisenberg=true, truncation::Function=identity, f! = f_unitary!)

In-place single step of Runge–Kutta-4 with a time dependant Hamiltonian.
Overwrites `O` with `O(t+dt)` and returns it.
`H` is a function that takes time `t` and returns the Hamiltonian at that time.
`truncation` : function that takes an operator and returns a truncated version of it. By default it is `identity` (no truncation).
"""
function rk4!(O::AbstractOperator, H::Function, dt::Real, t::Real; hbar::Real=1, heisenberg=true, truncation::Function=identity, f! = f_unitary!)
    return rk4!(O, H(t), dt; hbar=hbar, heisenberg=heisenberg, truncation=truncation, f! = f!)
end


# In-place Lindblad derivative:
#   Ȯ = (s/ħ)·[H, O] + Σᵢ γᵢ ( Aᵢ − ½{Lᵢ†Lᵢ, O} ),  Aᵢ = Lᵢ†OLᵢ (s=+i) | LᵢOLᵢ† (s=−i)
# Writes the result into C and returns it; C may alias O.
function f_lindblad!(C::AbstractOperator, H, O, s, hbar, L, gamma)
    # The dissipator is assembled from the *original* O first, because C may alias
    # O and gets overwritten by the unitary part below.
    Odiss = zero(typeof(O))
    for i in eachindex(L)
        if s == 1im
            A = L[i]' * O * L[i]
        elseif s == -1im
            A = L[i] * O * L[i]'
        else
            throw(ArgumentError("s should be 1im or -1im"))
        end
        LL = L[i]' * L[i]
        Odiss += gamma[i] * (A - anticommutator(LL, O) / 2)
    end
    f_unitary!(C, H, O, s, hbar)   # C := (s/ħ)·[H, O]   (safe even if C === O)
    return C + Odiss               # no in-place operator add exists, so this allocates
end


"""
    rk4_lindblad!(O::AbstractOperator, H::AbstractOperator, dt::Real, L; hbar::Real=1, heisenberg=true, truncation::Function=copy, gamma=[])

In-place single step of Runge–Kutta-4 for solving the Lindblad equation

``\\dot{O}=i[H,O]+\\sum_i \\gamma_i \\left(L_i^\\dagger O L_i -\\frac{1}{2} \\{ L_i^\\dagger L_i, O\\} \\right)``

Overwrites `O` with `O(t+dt)` and returns it.
`L` is a list of jump operators.
See [`rk4_lindblad`](@ref) for the allocating version.
"""
function rk4_lindblad!(O::AbstractOperator, H::AbstractOperator, dt::Real, L; hbar::Real=1, heisenberg=true, truncation::Function=copy, gamma=[])
    @assert length(gamma) == length(L) || length(gamma) == 0
    (length(gamma) == 0) && (gamma = ones(length(L)))
    f!(C, H, O, s, hbar) = f_lindblad!(C, H, O, s, hbar, L, gamma)
    return rk4!(O, H, dt; hbar=hbar, heisenberg=heisenberg, truncation=truncation, f! = f!)
end
