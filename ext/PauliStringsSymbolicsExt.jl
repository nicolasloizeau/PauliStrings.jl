module PauliStringsSymbolicsExt

using PauliStrings
import Symbolics
import Symbolics: simplify, Num, expand, simplify_fractions, factor
import PauliStrings: cutoff, binary_kernel, Operator, PauliString
import Base: show

function cutoff(o::Operator{P,Complex{Num}}, ::Real) where {P<:PauliString}
    return o
end

function binary_kernel(f, A::Operator{P,Complex{Num}}, B::Operator{P,Complex{Num}}; epsilon::Real=0, maxlength::Int=typemax(Int)) where {P<:PauliString}
    d = PauliStrings.emptydict(A)
    for (p1,c1) in zip(A.strings, A.coeffs), (p2,c2) in zip(B.strings, B.coeffs)
        p, k = f(p1, p2)
        c = c1 * c2 * k
        PauliStrings.setwith!(+, d, p, c)
    end
    o2 = Operator{P, typeof(first(values(d)))}(collect(keys(d)), collect(values(d)))
    return PauliStrings.compress(o2)
end

function show(io::IO, o::Operator{P,Complex{Num}}) where {P<:PauliString}
    for (p, c) in zip(o.strings, o.coeffs)
        print(io, c, ' ', string(p), '\n')
    end
end

"""Applies Symbolics.simplify to each coefficient and compresses like terms."""
function simplify(o::Operator{P,Complex{Num}}) where {P<:PauliString}
    newcoeffs = Complex{Num}[]
    for c in o.coeffs
        r = simplify_fractions(simplify(expand(real(c))))
        imc = simplify_fractions(simplify(expand(imag(c))))
        try
            r = factor(r)
            imc = factor(imc)
        catch
        end
        push!(newcoeffs, Complex{Num}(r, imc))
    end
    o2 = Operator{P, Complex{Num}}(copy(o.strings), newcoeffs)
end

end
