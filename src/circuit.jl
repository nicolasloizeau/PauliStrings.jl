

single_gates = ["X", "Y", "Z", "H", "S", "T", "Tdg", "Phase"]
two_gates = ["CNOT", "Swap", "CX", "CY", "CZ"]

allowed_gates = vcat(single_gates, two_gates, ["Noise"])



mutable struct Circuit
    N::Int
    gates::Vector{Tuple{String,Vector{Real}}}
    max_strings::Int
end

Circuit(N::Int; max_strings=2^30) = Circuit(N, [], max_strings)




function Base.push!(c::Circuit, gate::String, sites::Real...)
    if !(gate in allowed_gates)
        error("Unknown gate: $(gate)")
    end
    push!(c.gates, (gate, collect(sites)))
end



function XYZGate(N::Int, i::Int, type::String)
    O = Operator(N)
    O += type, i
    return O
end

XGate(N::Int, i::Int) = XYZGate(N, i, "X")
YGate(N::Int, i::Int) = XYZGate(N, i, "Y")
ZGate(N::Int, i::Int) = XYZGate(N, i, "Z")
HGate(N::Int, i::Int) = (XGate(N, i) + ZGate(N, i)) / sqrt(2)

PhaseGate(N::Int, i::Int, theta::Real) = (eye(N) + ZGate(N, i)) / 2 + (eye(N) - ZGate(N, i)) / 2 * exp(1im * theta)
SGate(N::Int, i::Int) = PhaseGate(N, i, pi / 2)
TGate(N::Int, i::Int) = PhaseGate(N, i, pi / 4)
TdgGate(N::Int, i::Int) = dagger(TGate(N::Int, i::Int))
SXGate(N::Int, i::Int) = ((1-1im)*XGate(N, i) + (1+1im)*eye(N))/2


function CXYZGate(N::Int, i::Int, j::Int, type::String)
    O = Operator(N)
    O -= "Z", i, type, j
    O += type, j
    O += "Z", i
    O += eye(N)
    return O / 2
end

CXGate(N::Int, i::Int, j::Int) = CXYZGate(N, i, j, "X")
CYGate(N::Int, i::Int, j::Int) = CXYZGate(N, i, j, "Y")
CZGate(N::Int, i::Int, j::Int) = CXYZGate(N, i, j, "Z")
CNOTGate(N::Int, i::Int, j::Int) = CXGate(N, i, j)

function SwapGate(N::Int, i::Int, j::Int)
    O = Operator(N)
    O += "X", i, "X", j
    O += "Y", i, "Y", j
    O += "Z", i, "Z", j
    O += eye(N)
    return O / 2
end

function CSXGate(N::Int, i::Int, j::Int)
    O = Operator(N)
    O += -1+1im, "Z", i, "X", j
    O += 1-1im, "Z", i
    O += 1-1im, "X", j
    O += eye(N)*(3+1im)
    return O / 4
end

CSXdgGate(N::Int, i::Int, j::Int) = dagger(CSXGate(N, i, j))

function compile(c::Circuit)
    U = eye(c.N)
    for (gate, args) in c.gates # args is a vector of the form [site1, site2, ... other args]
        if gate == "Noise"
            U = add_noise(U, args...)
        elseif gate in allowed_gates
            O = eval(Symbol(gate * "Gate"))(c.N, args...)
            U = O * U
        else
            error("Unknown gate: $gate")
        end
        U = compress(U)
        U = trim(U, c.max_strings)
    end
    return U
end


function CCXGate(N::Int, i::Int, j::Int, k::Int)
    O = Operator(N)
    O += "Z", i, "Z", j, "X", k
    O -= "Z", i, "X", k
    O -= "Z", j, "X", k
    O += "X", k
    O -= "Z", i, "Z", j
    O += "Z", i
    O += "Z", j
    O += eye(N) * 3
    return O / 4
end


# function MCXGate(N::int, sites::int...)
#     if length(sites) == 2
#         return CXGate(N, sites...)
#     elif length(sites) == 3
#         return CCXGate(N, sites...)
#     else
#         target = sites[end]
#         target2 = sites[end-1]
#         control = sites[1:end-2]


#     end
# end
