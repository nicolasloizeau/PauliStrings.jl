
import Symbolics as sym

function simplify(o::OperatorSymbolic)
    o2 = typeof(o)(o.N)
    for i in 1:length(o)
        c = sym.simplify(o.coef[i])
        if !iszero(c)
            push!(o2.coef, c)
            push!(o2.v, o.v[i])
            push!(o2.w, o.w[i])
        end
    end
    return o2
end
