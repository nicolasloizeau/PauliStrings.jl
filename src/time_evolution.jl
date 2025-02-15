



function f_unitary(H, O, s, hbar)
    return s * com(H, O) / hbar
end



"""
    rk4(H::Operator, O::Operator, dt::Real; hbar::Real=1, heisenberg=true, M=2^20, keep::Operator=Operator(N))

Single step of Runge–Kutta-4 with time independant Hamiltonian.
Returns O(t+dt).
Set `heisenberg=true` for evolving an observable in the heisenberg picture.
If `heisenberg=false` then it is assumed that O is a density matrix.
`M` is the number of strings to keep.
"""
function rk4(H::Operator, O::Operator, dt::Real; hbar::Real=1, heisenberg=true, M=2^20, keep::Operator=Operator(0))
    (keep.N == 0) && (keep = Operator(O.N))
    s = -1im
    heisenberg && (s = 1im)
    k1 = f_unitary(H, O, s, hbar)
    k1 = trim(k1, M; keep=keep)
    k2 = f_unitary(H, O + dt * k1 / 2, s, hbar)
    k2 = trim(k2, M; keep=keep)
    k3 = f_unitary(H, O + dt * k2 / 2, s, hbar)
    k3 = trim(k3, M; keep=keep)
    k4 = f_unitary(H, O + dt * k3, s, hbar)
    k4 = trim(k4, M; keep=keep)
    return O + (k1 + 2 * k2 + 2 * k3 + k4) * dt / 6
end

"""
    rk4(H::Function, O::Operator, dt::Real, t::Real; hbar::Real=1, heisenberg=false)

Single step of Runge–Kutta-4 with time dependant Hamiltonian.
`H` is a function that takes a number (time) and returns an operator.
"""
function rk4(H::Function, O::Operator, dt::Real, t::Real; hbar::Real=1, heisenberg=true)
    s = -1im
    heisenberg && (s *= -1)
    itime && (s *= 1im)
    k1 = s / hbar * com(H(t), O)
    k2 = s / hbar * com(H(t + dt / 2), O + dt * k1 / 2)
    k3 = s / hbar * com(H(t + dt / 2), O + dt * k2 / 2)
    k4 = s / hbar * com(H(t + dt), O + dt * k3)
    return O + (k1 + 2 * k2 + 2 * k3 + k4) * dt / 6
end



function f_lindblad(H, O, s, hbar, L, gamma)
    Odot = f_unitary(H, O, s, hbar)
    for i in 1:length(L)
        if s == 1im
            A = dagger(L[i]) * O * L[i]
        elseif s == -1im
            A = L[i] * O * dagger(L[i])
        else
            throw(ArgumentError("s should be 1im or -1im"))
        end
        LL = dagger(L[i]) * L[i]
        Odot += gamma[i] * (A - com(LL, O; anti=true) / 2)
    end
    return Odot
end




"""
    rk4_lindblad(H::Operator, O::Operator, dt::Real, L; hbar::Real=1, heisenberg=true, M=2^20, keep::Operator=Operator(0), gamma=[])

Single step of Runge–Kutta-4 for solving the Lindblad equation

``\\dot{O}=i[H,O]+\\sum_i \\gamma_i \\left(L_i^\\dagger O L_i -\\frac{1}{2} \\{ L_i^\\dagger L_i, O\\} \\right)``

Returns O(t+dt).
`L` is a list of jump operators.
Set `heisenberg=true` for evolving an observable in the heisenberg picture.
If `heisenberg=false` then it is assumed that O is a density matrix.
"""
function rk4_lindblad(H::Operator, O::Operator, dt::Real, L; hbar::Real=1, heisenberg=true, M=2^20, keep::Operator=Operator(0), gamma=[])
    @assert length(gamma) == length(L) || length(gamma) == 0
    if length(gamma) == 0
        gamma = ones(length(L))
    end
    (keep.N == 0) && (keep = Operator(O.N))
    s = -1im
    heisenberg && (s = 1im)
    k1 = f_lindblad(H, O, s, hbar, L, gamma)
    k1 = trim(k1, M; keep=keep)
    k2 = f_lindblad(H, O + dt * k1 / 2, s, hbar, L, gamma)
    k2 = trim(k2, M; keep=keep)
    k3 = f_lindblad(H, O + dt * k2 / 2, s, hbar, L, gamma)
    k3 = trim(k3, M; keep=keep)
    k4 = f_lindblad(H, O + dt * k3, s, hbar, L, gamma)
    k4 = trim(k4, M; keep=keep)
    return O + (k1 + 2 * k2 + 2 * k3 + k4) * dt / 6
end
