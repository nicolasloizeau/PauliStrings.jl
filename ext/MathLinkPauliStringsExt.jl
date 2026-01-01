module MathLinkPauliStringsExt
using LinearAlgebra
using PauliStrings
using MathLink
export OperatorMathLink, simplify_operator, simplify


# Define MathLinkNumber, a Number type that wraps MathLink expressions
# ---------------------------------------------------------------------


struct MathLinkNumber <: Number
    expression::Union{MathLink.WTypes, Number}
end



function Base.show(io::IO, x::MathLinkNumber)
    if x.expression isa MathLink.WExpr
        print(io, string(x.expression)[3:end-1])
    else
        print(io, string(x.expression))
    end
end


function LinearAlgebra.norm(x::Vector{MathLinkNumber})
    s = sum( xi->(W"Abs"(xi.expression))^2 , x)
    e = simplify(MathLinkNumber(W"Sqrt"(s)))
    return e
end


Base.:+(a::MathLinkNumber, b::MathLinkNumber) = MathLinkNumber(weval(a.expression + b.expression))
Base.:+(a::MathLinkNumber, b::Number) = MathLinkNumber(weval(a.expression + b))
Base.:+(a::Number, b::MathLinkNumber) = b+a

Base.:-(a::MathLinkNumber, b::MathLinkNumber) = MathLinkNumber(weval(a.expression - b.expression))
Base.:-(a::MathLinkNumber, b::Number) = MathLinkNumber(weval(a.expression - b))
Base.:-(a::Number, b::MathLinkNumber) = MathLinkNumber(weval(a - b.expression))
Base.:-(a::MathLinkNumber) = Number(weval(-a.expression))

Base.:*(a::MathLinkNumber, b::MathLinkNumber) = MathLinkNumber(weval(a.expression * b.expression))
Base.:*(a::MathLinkNumber, b::Number) = MathLinkNumber(weval(a.expression * b))
Base.:*(a::Number, b::MathLinkNumber) = b*a

Base.:/(a::MathLinkNumber, b::MathLinkNumber) = MathLinkNumber(weval(a.expression / b.expression))
Base.:/(a::MathLinkNumber, b::Number) = MathLinkNumber(weval(a.expression / b))
Base.:/(a::Number, b::MathLinkNumber) = MathLinkNumber(weval(a / b.expression))


Base.sqrt(a::MathLinkNumber) = MathLinkNumber(weval(W"Sqrt"(a.expression)))
Base.conj(a::MathLinkNumber) = MathLinkNumber(weval(W"Conjugate"(a.expression)))
Base.abs(a::MathLinkNumber) = MathLinkNumber(weval(W"Abs"(a.expression)))
Base.:^(a::MathLinkNumber, b::Integer) = MathLinkNumber(weval(a.expression ^ b))

function PauliStrings.simplify(expression::MathLink.WTypes; assumptions=nothing)
    if assumptions === nothing
        return weval(W"Simplify"(expression))
    else
        return weval(W"Simplify"(expression, assumptions))
    end
end

PauliStrings.simplify(a::Number; assumptions=nothing) = a


function PauliStrings.simplify(a::MathLinkNumber; assumptions=nothing)
    return MathLinkNumber(weval(simplify(a.expression, assumptions=assumptions)))
end

Base.Number(x::MathLink.WTypes) = MathLinkNumber(x)

# Define OperatorMathLink
# --------------------------

PauliStrings.OperatorMathLink(N::Int) = Operator{paulistringtype(N),MathLinkNumber}()

function Base.:+(o::Operator, args::Tuple{MathLink.WTypes,Vararg{Any}})
    args2 = (MathLinkNumber(args[1]), args[2:end]...)
    return o + args2
end

function PauliStrings.simplify_operator(o::Operator{P, MathLinkNumber}; assumptions=nothing) where {P}
    coeffs::Vector{MathLinkNumber} = [simplify(c, assumptions=assumptions) for c in o.coeffs]
    Operator{P, MathLinkNumber}(o.strings, coeffs)
end


Base.:+(o::Operator{P, MathLinkNumber}, a::MathLink.WTypes) where {P} = o + MathLinkNumber(a)
Base.:+(a::MathLink.WTypes, o::Operator{P, MathLinkNumber}) where {P} = o + a
Base.:-(o::Operator{P, MathLinkNumber}, a::MathLink.WTypes) where {P} = o - MathLinkNumber(a)
Base.:-(a::MathLink.WTypes, o::Operator{P, MathLinkNumber}) where {P} = -o + a
Base.:*(o::Operator{P, MathLinkNumber}, a::MathLink.WTypes) where {P} = o * MathLinkNumber(a)
Base.:*(a::MathLink.WTypes, o::Operator{P, MathLinkNumber}) where {P} = o * a
Base.:/(o::Operator{P, MathLinkNumber}, a::MathLink.WTypes) where {P} = o / MathLinkNumber(a)

function LinearAlgebra.norm(o::Operator{P, MathLinkNumber}; normalize=false) where {P}
    normalize ? norm(o.coeffs) : norm(o.coeffs) * MathLinkNumber(W"Sqrt"(2^(qubitlength(o))))
end

end
