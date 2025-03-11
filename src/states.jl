# functions that involve pure states

using PauliStrings.Circuits


"""
    trace_zpart(o::Operator)

Computes `<0|o|0>`.
"""
function trace_zpart(o::Operator)
    s = 0
    for i in 1:length(o)
        if xcount(o.v[i], o.w[i]) == 0 && ycount(o.v[i], o.w[i]) == 0
            s += o.coef[i]
        end
    end
    return s * 2^o.N
end


function get_ox(state)
    N = length(state)
        ox = eye(N)
    for i in 1:N
        if state[i] == '1'
            x = XGate(N, i)
            ox = ox*x
        end
    end
    return compress(ox)
end

"""
    expect(o::Operator, state::String)

Computes the expectation value `<state|o|state>`.
State is a single binary string that represents a pure state in the computational basis.
"""
function expect(o::Operator, state::String)
    @assert Set(state) ⊆ Set("01") "State must be a string of 0s and 1s"
    @assert o.N == length(state) "State length does not match operator size"
    ox = get_ox(state)
    o2 = ox*o*ox
    return trace_zpart(o2) / 2^o.N
end

"""
    expect_product(o1::Operator, o2::Operator, state::String)

Computes the expectation value `<state|o1*o2|state>`.
State is a single binary string that represents a pure state in the computational basis.
"""
function expect_product(o1::Operator, o2::Operator, state::String)
    @assert Set(state) ⊆ Set("01") "State must be a string of 0s and 1s"
    @assert o1.N == length(state) "State length does not match operator size"
    @assert o2.N == length(state) "State length does not match operator size"
    ox = get_ox(state)
    trace_product_z(ox*o1, o2*ox; scale=0) / 2^o1.N
end
