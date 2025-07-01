using Symbolics
using PauliStrings

"""
Initialize a zero operator on `N` qubits with symbolic coefficients. Uses `Symbolics.jl` to store this coefficients.
"""
OperatorSymbolic(N::Int) = Operator{paulistringtype(N),Complex{Num}}()

"""
Simplifies an Operator defined with symbolic coefficients. Uses `Symbolics.simplify` to simplify the symbolic 
expressions in each of the coefficients of `o`. Returns a new `Operator`.
"""
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

"""
Substitutes some or all of the variables in `o` according to the rule(s) in dict.
If all the substitutions are to concrete numeric values, then it will return an `Operator` with 
`Complex64` coefficients.
"""
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