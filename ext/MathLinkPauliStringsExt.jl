module MathLinkPauliStringsExt
using LinearAlgebra
using PauliStrings
using MathLink
using ProgressBars: ProgressBar
export OperatorMathLink, OperatorMathLinkTS, MathLinkNumber, simplify_operator, simplify, lanczos

# ── MathLinkNumber ─────────────────────────────────────────────────────────────

struct MathLinkNumber <: PauliStrings.MathLinkNumber
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
    s = sum(xi -> (W"Abs"(xi.expression))^2, x)
    return simplify(MathLinkNumber(W"Sqrt"(s)))
end

Base.:+(a::MathLinkNumber, b::MathLinkNumber) = MathLinkNumber(weval(a.expression + b.expression))
Base.:+(a::MathLinkNumber, b::Number)         = MathLinkNumber(weval(a.expression + b))
Base.:+(a::Number, b::MathLinkNumber)         = b + a
Base.:-(a::MathLinkNumber, b::MathLinkNumber) = MathLinkNumber(weval(a.expression - b.expression))
Base.:-(a::MathLinkNumber, b::Number)         = MathLinkNumber(weval(a.expression - b))
Base.:-(a::Number, b::MathLinkNumber)         = MathLinkNumber(weval(a - b.expression))
Base.:-(a::MathLinkNumber)                    = Number(weval(-a.expression))
Base.:*(a::MathLinkNumber, b::MathLinkNumber) = MathLinkNumber(weval(a.expression * b.expression))
Base.:*(a::MathLinkNumber, b::Number)         = MathLinkNumber(weval(a.expression * b))
Base.:*(a::Number, b::MathLinkNumber)         = b * a
Base.:/(a::MathLinkNumber, b::MathLinkNumber) = MathLinkNumber(weval(a.expression / b.expression))
Base.:/(a::MathLinkNumber, b::Number)         = MathLinkNumber(weval(a.expression / b))
Base.:/(a::Number, b::MathLinkNumber)         = MathLinkNumber(weval(a / b.expression))
Base.inv(a::MathLinkNumber)                   = MathLinkNumber(weval(1 // a.expression))
Base.sqrt(a::MathLinkNumber)                  = MathLinkNumber(weval(W"Sqrt"(a.expression)))
Base.conj(a::MathLinkNumber)                  = MathLinkNumber(weval(W"Conjugate"(a.expression)))
Base.abs(a::MathLinkNumber)                   = MathLinkNumber(weval(W"Abs"(a.expression)))
Base.:^(a::MathLinkNumber, b::Integer)        = MathLinkNumber(weval(a.expression ^ b))

Base.zero(::Type{MathLinkNumber}) = MathLinkNumber(0)
Base.one(::Type{MathLinkNumber})  = MathLinkNumber(1)
Base.zero(::MathLinkNumber)       = MathLinkNumber(0)
Base.one(::MathLinkNumber)        = MathLinkNumber(1)
Base.iszero(::MathLinkNumber)     = false

Base.Number(x::MathLink.WTypes) = MathLinkNumber(x)

# ── simplify ───────────────────────────────────────────────────────────────────

function PauliStrings.simplify(expression::MathLink.WTypes; assumptions=nothing)
    assumptions === nothing ?
        weval(W"Simplify"(expression)) :
        weval(W"Simplify"(expression, assumptions))
end
PauliStrings.simplify(a::Number; assumptions=nothing) = a
function PauliStrings.simplify(a::MathLinkNumber; assumptions=nothing)
    MathLinkNumber(weval(simplify(a.expression; assumptions=assumptions)))
end

# ── OperatorMathLink ───────────────────────────────────────────────────────────

PauliStrings.OperatorMathLink(N::Int) = Operator{paulistringtype(N), MathLinkNumber}()

# ── add_string helper (internal) ───────────────────────────────────────────────
# PauliStrings.add_string is not exported, so we define our own that routes
# a full-operator String through PauliString{N} before converting to PauliStringTS.

function _add_string_ts!(o::Operator{P, MathLinkNumber}, pauli::String, J::Number) where {P <: PauliStringTS}
    N  = qubitlength(o)
    ps = paulistringtype(N)(pauli)  # PauliString{N,...} from String
    Ls = qubitsize(P); Ps = periodicflags(P)
    p  = PauliStringTS{Ls,Ps}(ps)   # convert to PauliStringTS
    c  = (1im)^ycount(p) * J
    push!(o.strings, p)
    push!(o.coeffs, c)
    return o
end

# ── += overloads for TS+MathLink ───────────────────────────────────────────────
# The generic io.jl vararg dispatch builds multi-site operators via TS binary_kernel
# multiplication (term *= o2 for each site), which expands all N orbit translations
# instead of storing a single representative. We override all relevant paths to
# build the full Pauli char array directly and call _add_string_ts! once.

# Two-site Char: (J, 'X', i, 'Z', j)
function Base.:+(o::Operator{P, MathLinkNumber}, term::Tuple{Number,Char,Int,Char,Int}) where {P <: PauliStringTS}
    o1 = deepcopy(o)
    J, Pi, i, Pj, j = term
    pauli = fill('1', qubitlength(o1))
    pauli[i] = Pi
    pauli[j] = Pj
    _add_string_ts!(o1, join(pauli), J)
    return compress(o1)
end

# One-site Char: (J, 'X', i)
function Base.:+(o::Operator{P, MathLinkNumber}, term::Tuple{Number,Char,Int}) where {P <: PauliStringTS}
    o1 = deepcopy(o)
    J, Pi, i = term
    pauli = fill('1', qubitlength(o1))
    pauli[i] = Pi
    _add_string_ts!(o1, join(pauli), J)
    return compress(o1)
end

# Vararg String path: (J, "X", i, "Z", j, ...) — only handles basic XYZ Paulis.
# Intercept for TS to avoid binary_kernel expanding all N orbits.
function Base.:+(o::Operator{P, MathLinkNumber}, args::Tuple{Number,Vararg{Any}}) where {P <: PauliStringTS}
    N = qubitlength(o)
    pauli = fill('1', N)
    c = args[1]
    i = 2
    while i <= length(args)
        symbol = args[i]::String
        site   = args[i+1]::Int
        if occursin(symbol, "XYZ")
            pauli[site] = only(symbol)
        else
            # fall back to generic path for compound operators (S+, S-, Pup, Pdown, Sx, Sy, Sz)
            return invoke(Base.:+, Tuple{Operator, Tuple{Number,Vararg{Any}}}, o, args)
        end
        i += 2
    end
    o1 = deepcopy(o)
    _add_string_ts!(o1, join(pauli), c)
    return compress(o1)
end

Base.:+(o::Operator{P,MathLinkNumber}, args::Tuple{Vararg{Any}}) where {P <: PauliStringTS} = o + (1, args...)
Base.:-(o::Operator{P,MathLinkNumber}, args::Tuple{Number,Vararg{Any}}) where {P <: PauliStringTS} = o + (-args[1], args[2:end]...)
Base.:-(o::Operator{P,MathLinkNumber}, args::Tuple{Vararg{Any}}) where {P <: PauliStringTS} = o + (-1, args...)

# ── WTypes tuple dispatch (converts W`expr` → MathLinkNumber) ─────────────────

# TS + WTypes: more specific, resolves ambiguity with the Vararg{Any} overload above
function Base.:+(o::Operator{P,MathLinkNumber}, args::Tuple{MathLink.WTypes, Vararg{Any}}) where {P <: PauliStringTS}
    o + (MathLinkNumber(args[1]), args[2:end]...)
end
# Non-TS: generic fallback
function Base.:+(o::Operator, args::Tuple{MathLink.WTypes, Vararg{Any}})
    o + (MathLinkNumber(args[1]), args[2:end]...)
end
Base.:+(o::Operator{P,MathLinkNumber}, a::MathLink.WTypes) where {P} = o + MathLinkNumber(a)
Base.:+(a::MathLink.WTypes, o::Operator{P,MathLinkNumber}) where {P} = o + a
Base.:-(o::Operator{P,MathLinkNumber}, a::MathLink.WTypes) where {P} = o - MathLinkNumber(a)
Base.:-(a::MathLink.WTypes, o::Operator{P,MathLinkNumber}) where {P} = -o + a
Base.:*(o::Operator{P,MathLinkNumber}, a::MathLink.WTypes) where {P} = o * MathLinkNumber(a)
Base.:*(a::MathLink.WTypes, o::Operator{P,MathLinkNumber}) where {P} = o * a
Base.:/(o::Operator{P,MathLinkNumber}, a::MathLink.WTypes) where {P} = o / MathLinkNumber(a)
Base.:/(o::Operator{P,MathLinkNumber}, a::Int)             where {P} = o * W`1/$a`

# ── simplify_operator ──────────────────────────────────────────────────────────

function PauliStrings.simplify_operator(o::Operator{P, MathLinkNumber}; assumptions=nothing) where {P}
    coeffs = MathLinkNumber[simplify(c; assumptions=assumptions) for c in o.coeffs]
    Operator{P, MathLinkNumber}(o.strings, coeffs)
end

# ── norm ───────────────────────────────────────────────────────────────────────

function LinearAlgebra.norm(o::Operator{P, MathLinkNumber}; normalize=false) where {P}
    normalize ? norm(o.coeffs) : norm(o.coeffs) * MathLinkNumber(W"Sqrt"(2^(qubitlength(o))))
end

using PauliStrings: qubitsize
function LinearAlgebra.norm(o::Operator{P, MathLinkNumber}; normalize=false) where {P <: PauliStringTS}
    Ls    = qubitsize(o)
    scale = MathLinkNumber(W"Sqrt"(2^(Base.prod(Ls))))
    n     = sqrt(trace_product(o', o; scale=scale))
    normalize ? n / MathLinkNumber(W"Sqrt"(2^(qubitlength(o)))) : n
end

# ── OperatorMathLinkTS constructors ────────────────────────────────────────────

function (::Type{PauliStrings.OperatorMathLinkTS{Ls}})(::Int) where {Ls}
    P = PauliStrings.periodicpaulistringtype(Ls)
    return Operator{P, MathLinkNumber}()
end

function (::Type{PauliStrings.OperatorMathLinkTS{Ls}})(
    o::Operator{P, T};
    full::Bool = false
) where {Ls, P <: PauliString, T}
    full && (o = o / Base.prod(Ls))
    Ps = ntuple(_ -> true, length(Ls))
    periodic_strings = PauliStringTS{Ls, Ps}.(o.strings)
    coeffs = MathLinkNumber[c isa MathLinkNumber ? c : MathLinkNumber(c) for c in o.coeffs]
    return compress(Operator{eltype(periodic_strings), MathLinkNumber}(periodic_strings, coeffs))
end

# ── lanczos ────────────────────────────────────────────────────────────────────

function PauliStrings.lanczos(
    H::Operator{P, MathLinkNumber},
    O::Operator{P, MathLinkNumber},
    steps::Int;
    assumptions=nothing, returnOn=false, observer=false, show_progress=true
) where {P}
    @assert typeof(H) == typeof(O)
    @assert observer === false || returnOn === false
    O0 = deepcopy(O)
    O0 /= norm(O0, normalize=true)
    O0 = simplify_operator(O0; assumptions=assumptions)
    LHO = simplify_operator(commutator(H, O0); assumptions=assumptions)
    b   = simplify(norm(LHO, normalize=true); assumptions=assumptions)
    O1  = simplify_operator(commutator(H, O0) / b; assumptions=assumptions)
    bs  = [b]
    returnOn             && (Ons = [O0, O1])
    (observer !== false) && (obs = [observer(O0), observer(O1)])
    progress = show_progress ? ProgressBar : collect
    for _ in progress(0:steps-2)
        LHO = simplify_operator(commutator(H, O1); assumptions=assumptions)
        O2  = simplify_operator(LHO - b * O0; assumptions=assumptions)
        b   = simplify(norm(O2, normalize=true); assumptions=assumptions)
        O2 /= b
        returnOn             && push!(Ons, O2)
        (observer !== false) && push!(obs, observer(O2))
        O0 = deepcopy(O1)
        O1 = deepcopy(O2)
        push!(bs, b)
    end
    (observer !== false) && return (bs, obs)
    returnOn && return (bs, Ons)
    return bs
end

# TS Hamiltonian: expand via resum then run dense lanczos
function PauliStrings.lanczos(
    H::Operator{P, MathLinkNumber},
    O::Operator{Q, MathLinkNumber},
    steps::Int;
    assumptions=nothing, returnOn=false, observer=false, show_progress=true
) where {P <: PauliStringTS, Q}
    H_dense = resum(H)
    return PauliStrings.lanczos(
        H_dense, O, steps;
        assumptions=assumptions, returnOn=returnOn,
        observer=observer, show_progress=show_progress,
    )
end

end 