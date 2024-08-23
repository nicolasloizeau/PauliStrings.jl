"""
Runge–Kutta-4 with time independant Hamiltonian
heisenberg : set to true if evolving an observable
return : rho(t+dt)
"""
function rk4(H::Operator, O::Operator, dt::Real; hbar::Real=1, heisenberg=false, M=2^20, keep::Operator=Operator(N))
    s = 1
    if heisenberg
        s = -1
    end
    k1 = -s*1im/hbar*com(H, O)
    k1 = trim(k1, M; keep=keep)
    k2 = -s*1im/hbar*com(H, O+dt*k1/2)
    k2 = trim(k2, M; keep=keep)
    k3 = -s*1im/hbar*com(H, O+dt*k2/2)
    k3 = trim(k3, M; keep=keep)
    k4 = -s*1im/hbar*com(H, O+dt*k3)
    k4 = trim(k4, M; keep=keep)
    return O+(k1+2*k2+2*k3+k4)*dt/6
end

"""
Runge–Kutta-4 with time dependant Hamiltonian
H : function that takes a time and returns an operator
heisenberg : set to true if evolving an observable
return : rho(t+dt)
"""
function rk4(H::Function, O::Operator, dt::Real, t::Real; hbar::Real=1, heisenberg=false)
    s = 1
    if heisenberg
        s = -1
    end
    k1 = -s*1im/hbar*com(H(t), O)
    k2 = -s*1im/hbar*com(H(t+dt/2), O+dt*k1/2)
    k3 = -s*1im/hbar*com(H(t+dt/2), O+dt*k2/2)
    k4 = -s*1im/hbar*com(H(t+dt), O+dt*k3)
    return O+(k1+2*k2+2*k3+k4)*dt/6
end
