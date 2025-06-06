using Dictionaries

"""
    OperatorTS2D(o::Operator, L1::Int; full=true)

Initialize a 2D translation invariant operator from an Operator
\$O=\\sum_{i,j} o_{i,j} O_{i,j}\$ where \$O_{i,j}=T_{a_1}^i T_{a_2}^j(O_0)\$ is the application of translation operators in the direction of the basis vectors \$a_1\$ and \$\a_2\$,
\$i\$ and \$j\$ times respectively.
`L1` is the number of sites in the a1 direction. For example, if the lattice is 2x3, `L1 = 2` and `L2 = 3`, and the lattice is flattened as
[(1,1),(2,1),(1,2),(2,2),(1,3),(2,3)] (column-major order).
Set full=true if passing \$O\$, an Operator that is supported on the whole lattice (i.e converting from a translation symmetric [`Operator`](@ref))
Set full=false if passing \$O_0\$, a local term o such that the full operator is \$O=\\sum_{i,j} o_{i,j} T_a^i T_b^j(O_0)\$.
"""
function OperatorTS2D(o::Operator, L1::Int; full=true)
    if full && !is_ts2d(o, L1)
        error("o is not 2d translation symmetric.")
    end
    o2 = shift_origin(o, L1)
    full && (o2 /= qubitlength(o))

    return OperatorTS2D{eltype(o2.strings),eltype(o2.coeffs),L1}(o2.strings, o2.coeffs)
end

"""
    Operator(o::OperatorTS2D)

Convert an OperatorTS2D to an Operator
"""
function Operator(o::OperatorTS2D; rs=true)
    rs && (o = resum(o))
    return Operator(o.strings, o.coeffs)
end

"""
    is_ts2d(o::Operator)

return true if o is 2d translation symmetric.
"""
function is_ts2d(o::Operator, L1::Int)
    @assert qubitlength(o) % L1 == 0
    L2 = qubitlength(o) ÷ L1
    for i in 0:L1-1
        for j in 0:L2-1
            if opnorm(o - shift(o, i, j, L1)) / opnorm(o) > 1e-10
                return false
            end
        end
    end
    return true
end


"""rotate right the bits in the range [n1, n2) of x by r (this implies left rotation of a PauliString)"""
function rotate_chunk(x::Unsigned, n1::Int, n2::Int, r::Int)
    @assert n2 > n1
    mask = (one(x) << n2) - (one(x) << n1)
    bits_to_rotate = x & mask
    rotated_bits = (bits_to_rotate >> r) | (bits_to_rotate << (n2 - n1 - r))
    rotated_bits &= mask
    return (x & ~mask) | rotated_bits
end

""" 2d generalization of rotate_lower. Rotates the PauliStrng corresponding to a 2d lattice by r1 and r2 sites in the
a1 and a2 directions, respectively.
"""
function rotate_lower(p::PauliString{N,T}, r1::Int, r2::Int, L1::Int) where {N,T}
    L2 = N ÷ L1
    v = rotate_chunk(p.v, 0, N, r2 * L1)
    w = rotate_chunk(p.w, 0, N, r2 * L1)
    for i in 0:L2-1
        v = rotate_chunk(v, i * L1, (i + 1) * L1, r1)
        w = rotate_chunk(w, i * L1, (i + 1) * L1, r1)
    end
    return PauliString{N,T}(v, w)
end

"""Shift the string v, w so it starts on site 1"""
function shift_origin(p::PauliString{N,T}, L1::Int) where {N,T}
    L2 = N ÷ L1
    return maximum(rotate_lower(p, i, j, L1) for i in 0:L1-1, j in 0:L2-1)
end

shift(o::Operator, r1::Int, r2::Int, L1::Int) = Operator(rotate_lower.(o.strings, r1, r2, L1), copy(o.coeffs))

shift_origin(o::Operator, L1::Int) = compress(typeof(o)(shift_origin.(o.strings, L1), copy(o.coeffs)))
function shift_origin(o::OperatorTS2D)
    L1 = extent(o)
    return compress(typeof(o)(shift_origin.(o.strings, L1), copy(o.coeffs)))
end

function resum(o::OperatorTS2D)
    o2 = Operator(similar(o.strings, 0), similar(o.coeffs, 0))
    L1 = extent(o)
    N = qubitlength(o)
    L2 = N ÷ L1
    for i in 0:L1-1
        for j in 0:L2-1
            oij = Operator(copy(o.strings), copy(o.coeffs))
            o2 += shift(oij, i, j, L1)
        end
    end
    return o2
end


Base.:+(a::Number, o::OperatorTS2D) = OperatorTS2D(Operator(o, rs=false) + a / qubitlength(o), extent(o); full=false)
Base.:+(o::OperatorTS2D, a::Number) = a + o

function binary_kernel(op, A::OperatorTS2D, B::OperatorTS2D; epsilon::Real=0, maxlength::Int=1000)
    checklength(A, B)
    N = qubitlength(A)
    L1 = extent(A)
    L2 = N ÷ L1

    d = emptydict(A)
    p1s, c1s = A.strings, A.coeffs
    p2s, c2s = B.strings, B.coeffs

    # check lengths to safely use `@inbounds`
    length(p1s) == length(c1s) || throw(DimensionMismatch("strings and coefficients must have the same length"))
    length(p2s) == length(c2s) || throw(DimensionMismatch("strings and coefficients must have the same length"))

    # core kernel logic
    @inbounds for i1 in eachindex(p1s)
        p1, c1 = p1s[i1], c1s[i1]
        for i2 in eachindex(p2s)
            p2, c2 = p2s[i2], c2s[i2]
            for i in 0:L1-1
                for j in 0:L2-1
                    p, k = op(p1, rotate_lower(p2, i, j, L1))
                    c = c1 * c2 * k
                    if (k != 0) && (abs(c) > epsilon) && pauli_weight(p) < maxlength
                        setwith!(+, d, shift_origin(p, L1), c)
                    end
                end
            end
        end
    end

    return typeof(A)(collect(keys(d)), collect(values(d)))
end

Base.:*(o1::OperatorTS2D, o2::OperatorTS2D; kwargs...) = binary_kernel(prod, o1, o2; kwargs...)
commutator(o1::OperatorTS2D, o2::OperatorTS2D; kwargs...) = binary_kernel(commutator, o1, o2; kwargs...)
anticommutator(o1::OperatorTS2D, o2::OperatorTS2D; kwargs...) = binary_kernel(anticommutator, o1, o2; kwargs...)


trace(o::OperatorTS2D) = trace(Operator(o; rs=false)) * qubitlength(o)
opnorm(o::OperatorTS2D) = opnorm(Operator(o))
