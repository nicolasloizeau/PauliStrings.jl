

# Symbolics
export OperatorSymbolics, simplify_operator, substitute_operator
function OperatorSymbolics() end
function simplify_operator() end
function substitute_operator() end

# MathLink
export OperatorMathLink, OperatorMathLinkTS, simplify_operator, simplify
function OperatorMathLink() end
abstract type OperatorMathLinkTS{Ls} end
function simplify() end
