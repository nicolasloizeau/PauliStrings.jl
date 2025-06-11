import .Symbolics: simplify, substitute, Num

function OperatorSymbolic(N::Int)
    return Operator{paulistringtype(N),Complex{Num}}()
end

function simplify_op(o::Operator)
    o2 = typeof(o)()
    for i in 1:length(o)
        c = simplify(o.coeffs[i])
        if !iszero(c)
            push!(o2.coeffs, c)
            push!(o2.strings, o.strings[i])
        end
    end
    return o2
end

function substitute_op(o::Operator, dict::Dict)
    o = simplify_op(o)
    ps, cs = o.strings, o.coeffs
    cs_expr = substitute.(o.coeffs, (dict,))
    cs_vals = ComplexF64[]

    # Attempt to convert all the coefficients to ComplexF64, not possible if one or more variables remained unassigned
    all_vals = true
    for c in cs_expr
        try
            push!(cs_vals, ComplexF64(Symbolics.value(c)))
        catch
            all_vals = false
            break
        end
    end

    if all_vals
        return Operator{paulistringtype(qubitlength(o)),ComplexF64}(copy(ps), cs_vals)
    else
        return Operator{paulistringtype(qubitlength(o)),Complex{Num}}(copy(ps), cs_expr)
    end
end