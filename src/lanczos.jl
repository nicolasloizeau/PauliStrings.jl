


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

"""rotate left the qubits of O by r"""
function shift(O::Operator, r::Int)
    for i in 1:length(O)
        O.v[i]  = rotate_lower(O.v[i], O.N, r)
        O.w[i]  = rotate_lower(O.w[i], O.N, r)
    end
    return compress(O)
end

"""shift the string v,w so it starts on site 1"""
function shift_left(v,w,N)
    l = 2^(N+1)
    v2 = v
    w2 = w
    for i in 0:N
        v3 = rotate_lower(v, N, i)
        w3 = rotate_lower(w, N, i)
        if v3|w3+v3&w3/100000 < l
            v2 = v3
            w2 = w3
            l = v3|w3+v3&w3/100000
        end
    end
    return (v2,w2)
end

"""
    shift_left(O::Operator)

Shift evey string left so they start on site 1. This usefull for using translation symmetry in 1D systems
# Example
```julia
A = Operator(4)
A += "XYZ1"
A += "11ZZ"
A += "1XZY"
A += "ZZ11"
```
```julia
julia> shift_left(A)
(1.0 - 0.0im) XZY1
(1.0 - 0.0im) XYZ1
(2.0 + 0.0im) ZZ11
```
"""
function shift_left(O::Operator)
    for i in 1:length(O)
        O.v[i],O.w[i]  = shift_left(O.v[i],O.w[i],O.N)
    end
    return compress(O)
end

shift1(O::Operator) = shift_left(O)

"
    lanczos(H::Operator, O::Operator, steps::Int, nterms::Int; keepnorm=true, maxlength=1000, localop=false)

Computer the first `steps` lanczos coeficients for Hamiltonian `H` and initial operator `O`

At every step, the operator is trimed with `PauliStrings.trim` and only `nterms` are kept.

If H and O are 1D-translation-invariant, it is possible to provide a single local term of O and set `localop=true`

Using `maxlength` speeds up the commutator by only keeping terms of length <= `maxlength`
"
function lanczos(H::Operator, O::Operator, steps::Int, nterms::Int; keepnorm=true, maxlength=1000, localop=false)
    c = 1
    N = H.N
    localop && (c = sqrt(N))
    O0 = deepcopy(O)
    localop && (O0 = shift_left(O0))
    O0 /= norm_lanczos(O0)*c
    LHO = com(H, O0)
    localop && (LHO = shift_left(LHO))
    b = norm_lanczos(LHO)*c
    O1 = com(H, O0)/b
    localop && (O1 = shift_left(O1))
    bs = [b]
    for n in ProgressBar(0:steps-2)
        LHO = com(H, O1; maxlength=maxlength)
        localop && (LHO = shift_left(LHO))
        O2 = LHO-b*O0
        localop && (O2 = shift_left(O2))
        b = norm_lanczos(O2)*c
        O2 /= b
        O2 = trim(O2, nterms; keepnorm=keepnorm)
        O0 = deepcopy(O1)
        O1 = deepcopy(O2)
        push!(bs, b)
    end
    return bs
end
