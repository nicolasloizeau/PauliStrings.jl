


function norm_lanczos(O::Operator)
    return opnorm(O)/sqrt(2^O.N)
end


"""rotate left the first n bits of x by r"""
function rotate_lower(x::Int, n::Int, r::Int)
    mask = (1 << n) - 1
    lower_bits = x & mask
    rotated_bits = (lower_bits >> r) | (lower_bits << (n - r))
    rotated_bits &= mask
    return (x & ~mask) | rotated_bits
end

"""shift the string v,w so it starts on site 1"""
function shift1(v,w,N)
    l = 2^(N+1)
    v2 = v
    w2 = w
    for i in 1:N
        v3 = rotate_lower(v, N, i)
        w3 = rotate_lower(w, N, i)
        if v3+w3<l
            v2 = v3
            w2 = w3
            l = v3+w3
        end
    end
    return (v2,w2)
end

"""shift evey string so they start on site 1"""
function shift1(O::Operator)
    for i in 1:length(O)
        O.v[i],O.w[i]  = shift1(O.v[i],O.w[i],O.N)
    end
    return compress(O)
end


"
https://journals.aps.org/prx/pdf/10.1103/PhysRevX.9.041017 equation 4
H : hamiltonian MPO
O : operator MPO
steps : numer of lanczos steps
nterms : maximum number of terms in the operator. Used by trim at every step
localop : if H and O are 1d translation invariant, we can provide a single local term of O and set localop=true
"
function lanczos(H::Operator, O::Operator, steps::Int, nterms::Int; keepnorm=true, maxlength=1000, localop=false)
    c = 1
    N = H.N
    localop && (c = sqrt(N))
    O0 = deepcopy(O)
    O0 /= norm_lanczos(O0)*c
    b = norm_lanczos(com(H, O0))*c
    O1 = com(H, O0)/b
    localop && (O1 = shift1(O1))
    bs = [b]
    for n in ProgressBar(0:steps-2)
        LHO = com(H, O1; maxlength=maxlength)
        localop && (LHO = shift1(LHO))
        O2 = LHO-b*O0
        b = norm_lanczos(O2)*c
        O2 /= b
        O2 = trim(O2, nterms; keepnorm=keepnorm)
        O0 = deepcopy(O1)
        O1 = deepcopy(O2)
        push!(bs, b)
    end
    return bs
end
