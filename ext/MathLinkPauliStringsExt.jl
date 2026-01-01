module MathLinkPauliStringsExt
using LinearAlgebra
using PauliStrings
using MathLink
export OperatorMathLink, simplify_operator


# Define MathLinkNumber, a Number type that wraps MathLink expressions
# ---------------------------------------------------------------------


struct MathLinkNumber <: Number
    expression::Union{MathLink.WTypes, Number}
end



function Base.show(io::IO, x::MathLinkNumber)
    if x.expression isa MathLink.WSymbol
        print(io, string(x.expression))
    else
        print(io, string(x.expression)[3:end-1])
    end
end


LinearAlgebra.norm(x::Vector{MathLinkNumber}) = sqrt(sum( xi->abs(xi)^2 , x))


Base.:+(a::MathLinkNumber, b::MathLinkNumber) = MathLinkNumber(weval(a.expression + b.expression))
Base.:+(a::MathLinkNumber, b::Number) = MathLinkNumber(weval(a.expression + b))
Base.:+(a::Number, b::MathLinkNumber) = b+a

Base.:-(a::MathLinkNumber, b::MathLinkNumber) = MathLinkNumber(weval(a.expression - b.expression))
Base.:-(a::MathLinkNumber, b::Number) = MathLinkNumber(weval(a.expression - b))
Base.:-(a::Number, b::MathLinkNumber) = MathLinkNumber(weval(a - b.expression))
Base.:-(a::MathLinkNumber) = MathLinkNumber(weval(-a.expression))

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


simplify(a::MathLinkNumber) = MathLinkNumber(weval(W"Simplify"(a.expression)))

PauliStrings.simplify_operator(o::Operator{P, MathLinkNumber}) where {P} = Operator{P, MathLinkNumber}(o.strings, simplify.(o.coeffs))




# Define OperatorMathLink
# --------------------------

PauliStrings.OperatorMathLink(N::Int) = Operator{paulistringtype(N),MathLinkNumber}()

function Base.:+(o::Operator, args::Tuple{MathLink.WTypes,Vararg{Any}})
    args2 = (MathLinkNumber(args[1]), args[2:end]...)
    return o + args2
end

end
