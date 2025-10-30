


function add(o::AbstractOperator, s::AbstractPauliString; sign=1)
    o2 = copy(o)
    c = sign * (1im)^ycount(s)
    if s in o2.strings
        i = findfirst(==(s), o2.strings)
        o2.coeffs[i] += c
    else
        push!(o2.strings, s)
        push!(o2.coeffs, c)
    end
    return o2
end

Base.:+(o::AbstractOperator, s::AbstractPauliString) = add(o, s; sign=1)
Base.:-(o::AbstractOperator, s::AbstractPauliString) = add(o, s; sign=-1)
Base.:+(s::AbstractPauliString, o::AbstractOperator) = o + s
Base.:-(s::AbstractPauliString, o::AbstractOperator) = o - s

Base.:+(s1::T, s2::T) where {T<:PauliString} = Operator(qubitlength(s1)) + s1 + s2
Base.:-(s1::T, s2::T) where {T<:PauliString} = Operator(qubitlength(s1)) + s1 - s2



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



function Base.:*(s1::T, s2::T) where {T<:PauliString}
    p, k = prod(s1, s2)
    return Operator([p], [k])
end


Base.:+(s1::T, s2::T) where {T<:PauliStringTS} = OperatorTS(s1) + s2
Base.:-(s1::T, s2::T) where {T<:PauliStringTS} = OperatorTS(s1) - s2


function binary_kernel(op, A::Operator{<:PauliStringTS}, B::PauliStringTS; maxlength=1000, epsilon=1e-16)
    # checklength(A, B)
    Ls = qubitsize(A)

    d = emptydict(A)
    p1s, c1s = A.strings, A.coeffs
    p2 = B
    c2 = (1im)^ycount(p2)

    # check lengths to safely use `@inbounds`
    length(p1s) == length(c1s) || throw(DimensionMismatch("strings and coefficients must have the same length"))

    # core kernel logic
    @inbounds for (p1, c1) in zip(p1s, c1s)
        rep1 = representative(p1)
        rep2 = representative(p2)
        for s in all_shifts(paulistringtype(A))
            p, k = op(rep1, shift(rep2, Ls, s))
            c = c1 * c2 * k
            if (k != 0) && (abs(c) > epsilon) && pauli_weight(p) < maxlength
                setwith!(+, d, PauliStringTS{Ls}(p), c)
            end
        end
    end

    o = typeof(A)(collect(keys(d)), collect(values(d)))
    return (eltype(o.coeffs) == ComplexF64) ? cutoff(o, epsilon) : o
end


function Base.:*(o1::Operator{<:PauliStringTS}, o2::PauliStringTS)
    return binary_kernel(prod, o1, o2)
end

Base.:*(o2::PauliStringTS, o1::Operator{<:PauliStringTS}) = o1 * o2

function commutator(o1::Operator{<:PauliStringTS}, o2::PauliStringTS)
    return binary_kernel(commutator, o1, o2)
end

commutator(o2::PauliStringTS, o1::Operator{<:PauliStringTS}) = -commutator(o1, o2)

function anticommutator(o1::Operator{<:PauliStringTS}, o2::PauliStringTS)
    return binary_kernel(anticommutator, o1, o2)
end

anticommutator(o2::PauliStringTS, o1::Operator{<:PauliStringTS}) = anticommutator(o1, o2)



function binary_kernel(op, A::PauliStringTS, B::PauliStringTS; maxlength=1000, epsilon=1e-16)
    # checklength(A, B)
    Ls = qubitsize(A)

    d = emptydict(A)

    p1 = A
    c1 = (1im)^ycount(p1)
    p2 = B
    c2 = (1im)^ycount(p2)
    # core kernel logic

    rep1 = representative(p1)
    rep2 = representative(p2)
    for s in all_shifts(paulistringtype(A))
        p, k = op(rep1, shift(rep2, Ls, s))
        c = c1 * c2 * k
        if (k != 0) && (abs(c) > epsilon) && pauli_weight(p) < maxlength
            setwith!(+, d, PauliStringTS{Ls}(p), c)
        end
    end

    o = typeof(A)(collect(keys(d)), collect(values(d)))
    return (eltype(o.coeffs) == ComplexF64) ? cutoff(o, epsilon) : o
end
