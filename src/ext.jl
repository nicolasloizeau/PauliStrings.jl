# Symbolics
export OperatorSymbolics, simplify_operator, substitute_operator
function OperatorSymbolics() end
function simplify_operator() end
function substitute_operator() end

# MathLink
export OperatorMathLink, OperatorMathLinkTS, MathLinkNumber, simplify_operator, simplify
function OperatorMathLink() end
function simplify() end

"""
    OperatorMathLinkTS{Ls}

Parametric tag type whose constructors (defined in the MathLink extension) return
a translation-symmetric operator with symbolic `MathLinkNumber` coefficients.

    OperatorMathLinkTS{(N,)}(N)          # empty 1-D TS symbolic operator
    OperatorMathLinkTS{(N,)}(o::Operator) # convert existing operator to TS symbolic form

Requires `using MathLink` to activate.
"""
struct OperatorMathLinkTS{Ls} end

# Fallback constructors: give a clear error when MathLink is not loaded
(::Type{OperatorMathLinkTS{Ls}})(args...; kwargs...) where {Ls} =
    error("OperatorMathLinkTS requires MathLink to be loaded. Run `using MathLink` first.")

"""
    MathLinkNumber <: Number

Abstract stub type exported by `PauliStrings`. The concrete struct is defined
in the MathLink extension and subtypes this. Requires `using MathLink`.
"""
abstract type MathLinkNumber <: Number end
