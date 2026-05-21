


f_unitary(H, O, s, hbar) = s * commutator(H, O) / hbar


"""
    rk4(H::AbstractOperator, O::AbstractOperator, dt::Real; hbar::Real=1, heisenberg=true, M=2^20, trim::Function=trim)

Single step of Runge–Kutta-4 with time independant Hamiltonian.
Returns O(t+dt).
Set `heisenberg=true` for evolving an observable in the heisenberg picture.
If `heisenberg=false` then it is assumed that O is a density matrix.
`M` is the number of strings to keep.
"""
function rk4(H::AbstractOperator, O::AbstractOperator, dt::Real;
             hbar::Real=1, heisenberg=true, truncation::Function=copy, f=f_unitary)
    s = heisenberg ? 1im : -1im
    k1 = truncation(f(H, O, s, hbar))
    k2 = truncation(f(H, O + dt * k1 / 2, s, hbar))
    k3 = truncation(f(H, O + dt * k2 / 2, s, hbar))
    k4 = truncation(f(H, O + dt * k3, s, hbar))
    return O + (k1 + 2 * k2 + 2 * k3 + k4) * dt / 6
end




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
    rk4(H::Function, O::AbstractOperator, dt::Real, t::Real; hbar::Real=1, heisenberg=true, M=2^20, keep::Operator=Operator(0))

Single step of Runge–Kutta-4 with time dependant Hamiltonian.
Returns O(t+dt).
`H` is a function that takes a number (time) and returns an operator.
Set `heisenberg=true` for evolving an observable in the heisenberg picture.
If `heisenberg=false` then it is assumed that O is a density matrix.
`M` is the number of strings to keep.
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
    rk4_lindblad(H::AbstractOperator, O::AbstractOperator, dt::Real, L; hbar::Real=1, heisenberg=true, M=2^20, keep::Operator=Operator(0), gamma=[])

Single step of Runge–Kutta-4 for solving the Lindblad equation

``\\dot{O}=i[H,O]+\\sum_i \\gamma_i \\left(L_i^\\dagger O L_i -\\frac{1}{2} \\{ L_i^\\dagger L_i, O\\} \\right)``

Returns O(t+dt).
`L` is a list of jump operators.
Set `heisenberg=true` for evolving an observable in the heisenberg picture.
If `heisenberg=false` then it is assumed that O is a density matrix.
"""
function rk4_lindblad(H::AbstractOperator, O::AbstractOperator, dt::Real, L; hbar::Real=1, heisenberg=true, truncation::Function=copy, gamma=[])
    @assert length(gamma) == length(L) || length(gamma) == 0
    (length(gamma) == 0) && (gamma = ones(length(L)))
    f(H, O, s, hbar) = f_lindblad(H, O, s, hbar, L, gamma)
    return rk4(H, O, dt; hbar=hbar, heisenberg=heisenberg, truncation=truncation, f=f)
end
