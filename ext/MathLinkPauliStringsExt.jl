module MathLinkPauliStringsExt
using LinearAlgebra
using PauliStrings
using MathLink
using ProgressBars: ProgressBar
export OperatorMathLink, OperatorMathLinkTS, MathLinkNumber, simplify_operator, simplify, lanczos

using PauliStrings: PauliStringTS, periodicpaulistringtype, qubitsize, periodicflags,
    paulistringtype, ycount, compress

# Define MathLinkNumber, a Number type that wraps MathLink expressions
# ---------------------------------------------------------------------



"""
    MathLinkNumber(expression::Union{MathLink.WTypes, Number})

A wrapper type for MathLink expressions that behaves like a Number.
Do `x.expression` to get the underlying MathLink expression.
"""
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
    (length(x) == 0) && (return 0)
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
    return MathLinkNumber(weval(W"Power"(a.expression, -1)))
end


Base.sqrt(a::MathLinkNumber) = MathLinkNumber(weval(W"Sqrt"(a.expression)))
Base.conj(a::MathLinkNumber) = MathLinkNumber(weval(W"Conjugate"(a.expression)))
Base.abs(a::MathLinkNumber) = MathLinkNumber(weval(W"Abs"(a.expression)))
Base.:^(a::MathLinkNumber, b::Integer) = MathLinkNumber(weval(a.expression ^ b))

Base.zero(::Type{MathLinkNumber}) = MathLinkNumber(0)
Base.one(::Type{MathLinkNumber}) = MathLinkNumber(1)
Base.zero(::MathLinkNumber) = MathLinkNumber(0)
Base.one(::MathLinkNumber) = MathLinkNumber(1)

Base.iszero(a::MathLinkNumber) = a.expression isa Number && iszero(a.expression)

function PauliStrings.simplify(expression::MathLink.WTypes; assumptions=nothing)
    if assumptions === nothing
        return weval(W"Simplify"(expression))
    else
        return weval(W"Simplify"(expression, assumptions))
    end
end

PauliStrings.simplify(a::Number; assumptions=nothing) = a

"""
    simplify(a::MathLinkNumber; assumptions=nothing)

Simplifies a `MathLinkNumber` using Mathematica's `Simplify` function.
Assumptions can be provided, for example as ``assumptions = W`Assumptions -> {a > 0, b > 2}` ``.
"""
function PauliStrings.simplify(a::MathLinkNumber; assumptions=nothing)
    return MathLinkNumber(weval(simplify(a.expression, assumptions=assumptions)))
end

Base.Number(x::MathLink.WTypes) = MathLinkNumber(x)

# Define OperatorMathLink
# --------------------------

"""
    OperatorMathLink(N::Int)

Creates a `PauliStrings.Operator` for `N` qubits with [`MathLinkNumber`](@ref) coefficients.
"""
PauliStrings.OperatorMathLink(N::Int) = Operator{paulistringtype(N),MathLinkNumber}()


function Base.:+(o::Operator, args::Tuple{MathLink.WTypes,Vararg{Any}})
    args2 = (MathLinkNumber(args[1]), args[2:end]...)
    return o + args2
end

"""
    simplify_operator(o::Operator{P, MathLinkNumber}; assumptions=nothing) where {P}

Simplifies a `Operator{P, MathLinkNumber}` using Mathematica's `Simplify` function.
Assumptions can be provided, for example as ``assumptions = W`Assumptions -> {a > 0, b > 2}` ``.
"""
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

using PauliStrings: qubitsize
function LinearAlgebra.norm(o::Operator{P, MathLinkNumber}; normalize=false) where P<:PauliStringTS
    scale = MathLinkNumber(W"Power"(2, qubitlength(o)))
    n = sqrt(trace_product(o', o; scale=scale))
    if normalize
        return n / MathLinkNumber(W"Sqrt"(2^(qubitlength(o))))
    else
        return n
    end
end


function _add_string_ts!(o::Operator{P, MathLinkNumber}, pauli::String, J::Number) where {P<:PauliStringTS}
    ps = paulistringtype(qubitlength(o))(pauli)
    Ls = qubitsize(P); Ps = periodicflags(P)
    p = PauliStringTS{Ls,Ps}(ps)
    c = (1im)^ycount(p) * J
    push!(o.strings, p)
    push!(o.coeffs, c)
    return o
end

function Base.:+(o::Operator{P, MathLinkNumber}, term::Tuple{Number,Char,Int,Char,Int}) where {P<:PauliStringTS}
    o1 = deepcopy(o)
    J, Pi, i, Pj, j = term
    pauli = fill('1', qubitlength(o1))
    pauli[i] = Pi
    pauli[j] = Pj
    _add_string_ts!(o1, join(pauli), J)
    return compress(o1)
end

function Base.:+(o::Operator{P, MathLinkNumber}, term::Tuple{Number,Char,Int}) where {P<:PauliStringTS}
    o1 = deepcopy(o)
    J, Pi, i = term
    pauli = fill('1', qubitlength(o1))
    pauli[i] = Pi
    _add_string_ts!(o1, join(pauli), J)
    return compress(o1)
end

function Base.:+(o::Operator{P, MathLinkNumber}, args::Tuple{Number,Vararg{Any}}) where {P<:PauliStringTS}
    pauli = fill('1', qubitlength(o))
    c = args[1]
    i = 2
    while i <= length(args)
        symbol = args[i]::String
        site = args[i+1]::Int
        if occursin(symbol, "XYZ")
            pauli[site] = only(symbol)
        else
            return invoke(Base.:+, Tuple{Operator,Tuple{Number,Vararg{Any}}}, o, args)
        end
        i += 2
    end
    o1 = deepcopy(o)
    _add_string_ts!(o1, join(pauli), c)
    return compress(o1)
end

Base.:+(o::Operator{P, MathLinkNumber}, args::Tuple{Vararg{Any}}) where {P<:PauliStringTS} = o + (1, args...)
Base.:-(o::Operator{P, MathLinkNumber}, args::Tuple{Number,Vararg{Any}}) where {P<:PauliStringTS} = o + (-args[1], args[2:end]...)
Base.:-(o::Operator{P, MathLinkNumber}, args::Tuple{Vararg{Any}}) where {P<:PauliStringTS} = o + (-1, args...)

function Base.:+(o::Operator{P, MathLinkNumber}, args::Tuple{MathLink.WTypes,Vararg{Any}}) where {P<:PauliStringTS}
    return o + (MathLinkNumber(args[1]), args[2:end]...)
end

"""
    OperatorMathLinkTS{Ls}(N::Int)

Create an empty `N`-qubit translation symmetric operator with `MathLinkNumber` coefficients.
"""
function (::Type{OperatorMathLinkTS{Ls}})(N::Int) where {Ls}
    Base.prod(Ls) == N || error("OperatorMathLinkTS{$Ls}(N): prod(Ls)=$(Base.prod(Ls)) must equal N=$N")
    P = periodicpaulistringtype(Ls)
    return Operator{P, MathLinkNumber}()
end

@doc raw"""
    OperatorMathLinkTS{Ls}(o::Operator{P, MathLinkNumber}; full=false)

Convert a MathLink operator `o` to a translation symmetric operator, mapping its
Pauli strings to orbit representatives.

Set `full=true` when `o` is supported on the whole chain (i.e. ``O=\sum_i T_i(O_0)``);
the coefficients are then divided by ``\prod Ls`` using an exact Mathematica rational.
Set `full=false` (default) when `o` holds only the local seed ``O_0``.
"""
function (::Type{OperatorMathLinkTS{Ls}})(o::Operator{P, MathLinkNumber}; full::Bool=false) where {Ls, P<:PauliString}
    Ps = ntuple(_ -> true, length(Ls))
    periodic_strings = PauliStringTS{Ls,Ps}.(o.strings)
    coeffs = copy(o.coeffs)
    o_ts = compress(Operator{eltype(periodic_strings), MathLinkNumber}(periodic_strings, coeffs))
    full && (o_ts = o_ts * W`1/$(Base.prod(Ls))`)
    return o_ts
end















# Symbolic Lanczos algorithm
# ---------------------------

"""
    lanczos(H::Operator{P, MathLinkNumber}, O::Operator{P, MathLinkNumber}, steps::Int; assumptions=nothing, returnOn=false, observer=false, show_progress=true) where {P}

Lanczos algorithm for symbolic `MathLink` operators. Assumptions can be provided to simplify the expressions during the algorithm (cf [`simplify`](@ref))
"""

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
