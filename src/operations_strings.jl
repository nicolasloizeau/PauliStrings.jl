


function add(o::Operator, s::PauliString; sign=1)
    o2 = Operator(o)
    c = (1im)^ycount(s)
    if s in o2.strings
        i = findfirst(==(s), o2.strings)
        o2.coeffs[i] += c
    else
        push!(o2.strings, s)
        push!(o2.coeffs, c)
    end
    return o2
end

Base.:+(o::Operator, s::PauliString) = add(o, s; sign=1)
Base.:-(o::Operator, s::PauliString) = add(o, s; sign=-1)
Base.:+(s::PauliString, o::Operator) = o + s
Base.:-(s::PauliString, o::Operator) = o - s

function Base.:+(s1::T, s2::T) where {T<:PauliString}
    O = Operator(qubitlength(s1))
    O += s1
    O += s2
    return O
end


function binary_kernel(f, A::Operator, B::PauliString)
    # checklength(A, B)

    d = emptydict(A) # reducer
    p1s, c1s = A.strings, A.coeffs

    # check boundaries to safely use `@inbounds`
    length(p1s) == length(c1s) || throw(DimensionMismatch("strings and coefficients must have the same length"))

    # core kernel logic
    p2 = B
    c2 = (1im)^ycount(p2)
    @inbounds for i1 in eachindex(p1s)
        p1, c1 = p1s[i1], c1s[i1]
        p, k = f(p1, p2)
        c = c1 * c2 * k
        if (k != 0)
            setwith!(+, d, p, c)
        end
    end

    # assemble output
    o = Operator{keytype(d),valtype(d)}(collect(keys(d)), collect(values(d)))
    return (eltype(o.coeffs) == ComplexF64) ? cutoff(o, 1e-16) : o
end


function Base.:*(o1::Operator, o2::PauliString)
    return binary_kernel(prod, o1, o2)
end

Base.:*(o2::PauliString, o1::Operator) = o1 * o2

function commutator(o1::Operator, o2::PauliString)
    return binary_kernel(commutator, o1, o2)
end

commutator(o2::PauliString, o1::Operator) = -commutator(o1, o2)

function anticommutator(o1::Operator, o2::PauliString)
    return binary_kernel(anticommutator, o1, o2)
end

anticommutator(o2::PauliString, o1::Operator) = anticommutator(o1, o2)
