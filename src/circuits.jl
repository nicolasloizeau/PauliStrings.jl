

module Circuits
export CCXGate, TGate, TdgGate, HGate, SwapGate, PhaseGate, SXGate
export SGate, SXGate
export RXGate, RYGate, RZGate
export UGate
export CXGate, CYGate, CZGate, CNOTGate, CPhaseGate
export MCZGate
export CSXGate, CSXdgGate
export Circuit, compile, expect
export from_openqasm, from_openqasm_file
export grover_diffusion
export XGate, YGate, ZGate
export XXPlusYYGate
using PauliStrings

single_gates = ["X", "Y", "Z", "H", "RX", "RY", "RZ", "S", "SX", "T", "Tdg", "Phase", "U"]
two_gates = ["CNOT", "Swap", "CX", "CY", "CZ", "CSX", "CSXdg", "XXPlusYY", "CPhase"]
other = ["CCX", "Noise", "MCZ"]

allowed_gates = vcat(single_gates, two_gates, other)



mutable struct Circuit
    N::Int
    gates::Vector{Tuple{String,Vector{Int},Vector{Real}}}
    max_strings::Int
    noise_amplitude::Real
end


"""
    Circuit(N::Int; max_strings=2^30, noise_amplitude=0)

Creates an empty quantum circuit with `N` qubits. `max_strings` is the maximum number of strings to keep in the operator. `noise_amplitude` is the amplitude of the noise gate.
"""
Circuit(N::Int; max_strings=2^30, noise_amplitude=0) = Circuit(N, [], max_strings, noise_amplitude)


function get_site_pars(gate::String, site_pars)
    site_pars = collect(site_pars)
    if gate in single_gates
        sites = [Int(site_pars[1])]
        pars = site_pars[2:end]
    elseif gate in two_gates
        sites = [Int(site_pars[1]), Int(site_pars[2])]
        pars = site_pars[3:end]
    else
        gate in other
        sites = Int.(site_pars)
        pars = []
    end
    return sites, pars
end


"""
    push!(c::Circuit, gate::String, sites::Real...)

Adds a gate to the circuit `c`. The gate is specified by a string `gate` and a list of sites `sites`.
The gates have the same naming convention as in Qiskit.
Allowed gates are: "X", "Y", "Z", "H", "S", "T", "Tdg", "Phase", "CNOT", "Swap", "CX", "CY", "CZ", "CCX", "CSX", "CSXdg", "MCZ", "Noise".
"""
function Base.push!(c::Circuit, gate::String, site_pars::Real...)
    @assert gate in allowed_gates "Unknown gate: $(gate)"
    sites, pars = get_site_pars(gate, site_pars)
    push!(c.gates, (gate, sites, pars))
end

"""
    pushfirst!(c::Circuit, gate::String, sites::Real...)

Adds a gate to the beginning of the circuit `c`.
"""
function Base.pushfirst!(c::Circuit, gate::String, site_pars::Real...)
    @assert gate in allowed_gates "Unknown gate: $(gate)"
    sites, pars = get_site_pars(gate, site_pars)
    pushfirst!(c.gates, (gate, sites, pars))
end


const OPENQASM_GATE_ALIASES = Dict(
    "x" => "X",
    "y" => "Y",
    "z" => "Z",
    "h" => "H",
    "s" => "S",
    "sx" => "SX",
    "t" => "T",
    "tdg" => "Tdg",
    "sdg" => "Phase",
    "rx" => "RX",
    "ry" => "RY",
    "rz" => "RZ",
    "u1" => "Phase",
    "u2" => "U",
    "p" => "Phase",
    "phase" => "Phase",
    "u" => "U",
    "u3" => "U",
    "cx" => "CX",
    "cnot" => "CNOT",
    "cy" => "CY",
    "cz" => "CZ",
    "swap" => "Swap",
    "ccx" => "CCX",
    "cp" => "CPhase",
    "cu1" => "CPhase",
    "id" => "Id",
)


function _strip_openqasm_comments(src::AbstractString)
    src = replace(src, r"/\*.*?\*/"s => "")
    return join(first.(split.(split(src, '\n'), "//")), "\n")
end


function _openqasm_value(expr::AbstractString)
    expr = strip(expr)
    if !occursin(r"^[0-9eEpiPI\.\+\-\*\/\(\)\s]+$", expr)
        error("Unsupported OpenQASM parameter expression: $expr")
    end
    expr = replace(lowercase(expr), "pi" => string(pi))
    parsed = Meta.parse(expr)
    value = Core.eval(@__MODULE__, parsed)
    value isa Real || error("OpenQASM parameter is not real-valued: $expr")
    return value
end


function _openqasm_sites(args::AbstractString, registers)
    sites = Int[]
    for arg in split(args, ',')
        m = match(r"^\s*(\w+)\[(\d+)\]\s*$", arg)
        m === nothing && error("Unsupported OpenQASM qubit operand: $arg")
        name = m.captures[1]
        idx = parse(Int, m.captures[2])
        haskey(registers, name) || error("Unknown OpenQASM register: $name")
        offset, len = registers[name]
        0 <= idx < len || error("Qubit index $idx out of bounds for $name")
        push!(sites, offset + idx + 1)
    end
    return sites
end


"""
    from_openqasm(src::AbstractString; max_strings=2^30, noise_amplitude=0)

Parse a basic OpenQASM 2.0 circuit into a [`Circuit`](@ref). The parser
supports `qreg` declarations and common gates from `qelib1.inc`. Measurements,
classical registers, and barriers are ignored because `Circuit` represents the
unitary part of a circuit.
"""
function from_openqasm(src::AbstractString; max_strings=2^30, noise_amplitude=0)
    registers = Dict{String,Tuple{Int,Int}}()
    c = nothing
    nqubits = 0

    src = _strip_openqasm_comments(src)
    for raw_statement in split(src, ';')
        line = strip(raw_statement)
        isempty(line) && continue

        if startswith(line, "OPENQASM") || startswith(line, "include")
            continue
        end

        m = match(r"^qreg\s+(\w+)\[(\d+)\]$", line)
        if m !== nothing
            name = m.captures[1]
            len = parse(Int, m.captures[2])
            len > 0 || error("Invalid OpenQASM qreg size: $len")
            !haskey(registers, name) ||
                error("Duplicate OpenQASM qreg declaration: $name")
            offset = nqubits
            registers[name] = (offset, len)
            nqubits += len
            if c === nothing
                c = Circuit(
                    nqubits;
                    max_strings=max_strings,
                    noise_amplitude=noise_amplitude,
                )
            else
                c.N = nqubits
            end
            continue
        end

        if startswith(line, "creg") || startswith(line, "barrier") ||
           startswith(line, "measure")
            continue
        end

        startswith(line, "gate ") &&
            error("Custom OpenQASM gate definitions are not supported")
        startswith(line, "opaque ") &&
            error("Opaque OpenQASM gates are not supported")
        startswith(line, "if") &&
            error("Classically controlled OpenQASM gates are not supported")

        c === nothing && error("OpenQASM circuit has no qreg declaration")

        m = match(r"^(\w+)(?:\((.*)\))?\s+(.+)$", line)
        m === nothing && error("Could not parse OpenQASM statement: $line")
        qasm_gate, par_src, arg_src = m.captures
        key = lowercase(qasm_gate)
        haskey(OPENQASM_GATE_ALIASES, key) ||
            error("Unsupported OpenQASM gate: $qasm_gate")
        gate = OPENQASM_GATE_ALIASES[key]
        sites = _openqasm_sites(arg_src, registers)
        pars = isnothing(par_src) ? Real[] : Real[
            _openqasm_value(par) for par in split(par_src, ',')
        ]
        if key == "id"
            isempty(pars) || error("OpenQASM id gate does not take parameters")
            continue
        elseif key == "sdg"
            isempty(pars) || error("OpenQASM sdg gate does not take parameters")
            push!(c, gate, sites..., -pi / 2)
            continue
        elseif key == "u2"
            length(pars) == 2 || error("OpenQASM u2 gate expects 2 parameters")
            push!(c, gate, sites..., pi / 2, pars...)
            continue
        end
        push!(c, gate, sites..., pars...)
    end

    c === nothing && error("OpenQASM circuit has no qreg declaration")
    return c
end


"""
    from_openqasm_file(path::AbstractString; kwargs...)

Load a basic OpenQASM 2.0 circuit file into a [`Circuit`](@ref).
"""
function from_openqasm_file(path::AbstractString; kwargs...)
    return from_openqasm(read(path, String); kwargs...)
end



function XYZGate(N::Int, i::Int, type::String)
    O = Operator(N)
    O += type, i
    return O
end



"""
    XGate(N::Int, i::Int)
    YGate(N::Int, i::Int)
    ZGate(N::Int, i::Int)
    HGate(N::Int, i::Int)
    SGate(N::Int, i::Int)
    TGate(N::Int, i::Int)
    TdgGate(N::Int, i::Int)
    SXGate(N::Int, i::Int)

Creates a single qubit gate acting on qubit `i` of a `N` qubit system.
"""
XGate(N::Int, i::Int) = XYZGate(N, i, "X")
YGate(N::Int, i::Int) = XYZGate(N, i, "Y")
ZGate(N::Int, i::Int) = XYZGate(N, i, "Z")
HGate(N::Int, i::Int) = (XGate(N, i) + ZGate(N, i)) / sqrt(2)
SGate(N::Int, i::Int) = PhaseGate(N, i, pi / 2)
TGate(N::Int, i::Int) = PhaseGate(N, i, pi / 4)
TdgGate(N::Int, i::Int) = TGate(N::Int, i::Int)'
SXGate(N::Int, i::Int) = ((1 - 1im) * XGate(N, i) + (1 + 1im) * eye(N)) / 2

"""
    PhaseGate(N::Int, i::Int, theta::Real)

Creates a phase gate acting on qubit `i` of a `N` qubit system with phase `theta`.
"""
PhaseGate(N::Int, i::Int, theta::Real) = (eye(N) + ZGate(N, i)) / 2 + (eye(N) - ZGate(N, i)) / 2 * exp(1im * theta)

"""
    UGate(N::Int, i::Int, theta::Real, phi::Real, lam::Real)

General 1-qubit rotation of qubit `i` of a `N` qubit system with Euler angles `theta`, `phi`, `lam` of form

``U(\\theta, \\phi, \\lambda) =\\begin{pmatrix}\\cos\\left(\\theta/2\\right) & -e^{i\\lambda}\\sin\\left(\\theta/2\\right) \\\\e^{i\\phi}\\sin\\left(\\theta/2\\right) & e^{i(\\phi+\\lambda)}\\cos\\left(\\theta/2\\right)\\end{pmatrix}``

"""
UGate(N::Int, i::Int, theta::Real, phi::Real, lam::Real) =
	(cos(theta / 2) * (1 + exp(1im * (lam + phi))) / 2 * eye(N) -
	sin(theta / 2) * (exp(1im * lam) - exp(1im * phi)) / 2 * XGate(N, i) -
	1im * sin(theta / 2) * (exp(1im * lam) + exp(1im * phi)) / 2 * YGate(N, i) -
	cos(theta / 2) * (-1 + exp(1im * (lam + phi))) / 2 * ZGate(N, i))

"""
    RXGate(N::Int, i::Int, phi::Real)
    RYGate(N::Int, i::Int, theta::Real)
    RZGate(N::Int, i::Int, phi::Real)

1-qubit rotation of qubit `i` of a `N` qubit system around specific axis.
"""
RXGate(N::Int, i::Int, theta::Real) = UGate(N, i, theta, pi/2, -pi/2)
RYGate(N::Int, i::Int, theta::Real) = UGate(N, i, theta, 0, 0)
RZGate(N::Int, i::Int, phi::Real) = exp(1im * phi / 2) * UGate(N, i, 0, 0, phi)

"""
    CPhaseGate(N::Int, i::Int, j::Int, theta::Real)

Controled phase gate with control qubit `i` and target qubit `j` of a `N` qubit system.
"""
function CPhaseGate(N::Int, i::Int, j::Int, theta::Real)
    O = Operator(N)
    c = (exp(1im * theta) - 1) / 4
    O += c, "Z", i, "Z", j
    O -= c, "Z", i
    O -= c, "Z", j
    O += (exp(1im * theta) + 3) / 4
    return O
end


function CXYZGate(N::Int, i::Int, j::Int, type::String)
    O = Operator(N)
    O -= "Z", i, type, j
    O += type, j
    O += "Z", i
    O += eye(N)
    return O / 2
end


"""
    CXGate(N::Int, i::Int, j::Int)
    CYGate(N::Int, i::Int, j::Int)
    CZGate(N::Int, i::Int, j::Int)
    CNOTGate(N::Int, i::Int, j::Int)

Creates a two qubit gate with control qubit `i` and target qubit `j` of a `N` qubit system.
"""
CXGate(N::Int, i::Int, j::Int) = CXYZGate(N, i, j, "X")
CYGate(N::Int, i::Int, j::Int) = CXYZGate(N, i, j, "Y")
CZGate(N::Int, i::Int, j::Int) = CXYZGate(N, i, j, "Z")
CNOTGate(N::Int, i::Int, j::Int) = CXGate(N, i, j)


"""
    SwapGate(N::Int, i::Int, j::Int)

Creates a swap gate between qubits `i` and `j` of a `N` qubit system.
"""
function SwapGate(N::Int, i::Int, j::Int)
    O = Operator(N)
    O += "X", i, "X", j
    O += "Y", i, "Y", j
    O += "Z", i, "Z", j
    O += eye(N)
    return O / 2
end


"""
    XXPlusYYGate(N::Int, i::Int, j::Int, theta::Real, beta::Real)

XX+YY gate between qubits `i` and `j` of a `N` qubit system.
"""
function XXPlusYYGate(N::Int, i::Int, j::Int, theta::Real, beta::Real)
    O = Operator(N)
    O += 0.5, "Z", i, "Z", j
    O += 0.5
    O += cos(theta / 2) / 2
    O -= cos(theta / 2) / 2, "Z", i, "Z", j
    c = -1im * sin(theta / 2) * exp(1im * beta) / 4
    O += c, "X", i, "X", j
    O -= 1im * c, "Z", i, "Y", j
    O += 1im * c, "Y", i, "Z", j
    O += c, "Y", i, "Y", j
    c = -1im * sin(theta / 2) * exp(-1im * beta) / 4
    O += c, "X", i, "X", j
    O += 1im * c, "Z", i, "Y", j
    O -= 1im * c, "Y", i, "Z", j
    O += c, "Y", i, "Y", j
    return O
end

"""
    CSXGate(N::Int, i::Int, j::Int)
    CSXdgGate(N::Int, i::Int, j::Int)

Controlled sqrt X gate and its dagger
"""
function CSXGate(N::Int, i::Int, j::Int)
    O = Operator(N)
    O += -1 + 1im, "Z", i, "X", j
    O += 1 - 1im, "Z", i
    O += 1 - 1im, "X", j
    O += eye(N) * (3 + 1im)
    return O / 4
end
CSXdgGate(N::Int, i::Int, j::Int) = CSXGate(N, i, j)'

"""
    CCXGate(N::Int, i::Int, j::Int, k::Int)

Tofolli gate with control qubits `i` and `j` and target qubit `k` of a `N` qubit system.
"""
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

"""
    MCZGate(N::Int, sites::Int...)

Creates a multi-controlled Z gate acting on `sites` qubits of a `N` qubit system.
"""
function MCZGate(N::Int, sites::Int...)
    sites = collect(sites)
    U = eye(N) - 2 * all_z(N, sites) / 2^N
    for i in sites
        U = XGate(N, i) * U * XGate(N, i)
    end
    return compress(U)
end

MCZGate(N::Int) = MCZGate(N, 1:N...)

"""
    grover_diffusion(N::Int, sites::Int...)

Creates the Grover diffusion operator acting on `sites` qubits of a `N` qubit system.
"""
function grover_diffusion(N::Int, sites::Int...)
    U = MCZGate(N, sites...)
    for i in sites
        U = HGate(N, i) * U * HGate(N, i)
    end
    return compress(U)
end



"""
    compile(c::Circuit)

Compiles the quantum circuit `c` into a unitary operator. Applies the gates in the order they were added.
Applies noise gates if present and trim the operator to `c.max_strings` strings at each step.
"""
function compile(c::Circuit)
    U = eye(c.N)
    for (gate, sites, args) in c.gates
        if gate == "Noise"
            U = add_noise(U, c.noise_amplitude)
        elseif gate in allowed_gates
            O = eval(Symbol(gate * "Gate"))(c.N, sites..., args...)
            U = O * U
        else
            error("Unknown gate: $gate")
        end
        U = compress(U)
        U = trim(U, c.max_strings)
    end
    return U
end




end
