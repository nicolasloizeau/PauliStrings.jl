module MathLinkPauliStringsExt
using LinearAlgebra
using PauliStrings
using MathLink
using ProgressBars: ProgressBar
export OperatorMathLink, simplify_operator, simplify, lanczos


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

function Base.inv(a::MathLinkNumber)
    return MathLinkNumber(weval(1 // a.expression))
end


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

Base.:/(o::Operator{P, MathLinkNumber}, a::Int) where {P} = o * W`1/$a`

function LinearAlgebra.norm(o::Operator{P, MathLinkNumber}; normalize=false) where {P}
    normalize ? norm(o.coeffs) : norm(o.coeffs) * MathLinkNumber(W"Sqrt"(2^(qubitlength(o))))
end



# Symbolic Lanczos algorithm
# ---------------------------


function PauliStrings.lanczos(H::Operator{P, MathLinkNumber}, O::Operator{P, MathLinkNumber}, steps::Int; assumptions=nothing, returnOn=false, observer=false, show_progress=true) where {P}
    @assert typeof(H) == typeof(O)
    @assert observer === false || returnOn === false
    O0 = deepcopy(O)
    O0 /= norm(O0, normalize=true)
    O0 = simplify_operator(O0, assumptions=assumptions)
    LHO = simplify_operator(commutator(H, O0), assumptions=assumptions)
    b = simplify(norm(LHO, normalize=true), assumptions=assumptions)
    O1 = simplify_operator(commutator(H, O0) / b, assumptions=assumptions)
    bs = [b]
    returnOn && (Ons = [O0, O1])
    (observer !== false) && (obs = [observer(O0), observer(O1)])
    progress = collect
    show_progress && (progress = ProgressBar)
    for n in progress(0:steps-2)
        LHO = simplify_operator(commutator(H, O1), assumptions=assumptions)
        O2 = simplify_operator(LHO - b * O0, assumptions=assumptions)
        b = simplify(norm(O2, normalize=true), assumptions=assumptions)
        O2 /= b
        returnOn && push!(Ons, O2)
        (observer !== false) && push!(obs, observer(O2))
        O0 = deepcopy(O1)
        O1 = deepcopy(O2)
        push!(bs, b)
    end
    (observer !== false) && return (bs, obs)
    returnOn && (return bs, Ons)
    return bs
end


end
