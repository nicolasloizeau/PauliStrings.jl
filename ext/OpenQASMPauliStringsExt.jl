module OpenQASMPauliStringsExt

using PauliStrings
using PauliStrings.Circuits: Circuit
import PauliStrings.Circuits: parse_qasm, load_qasm
using OpenQASM
using OpenQASM.Types: MainProgram, Instruction, RegDecl, Bit, Neg, Call

# Evaluate an OpenQASM gate-parameter expression directly from the parsed AST. A parameter
# node is one of: a lexer token (number, or the reserved word `pi`), a `Neg` (unary minus),
# an infix `Tuple` `(left, op_token, right)`, or a `Call` (qelib1 function). We evaluate
# over this fixed grammar only, never arbitrary Julia.
function _eval_param(node)
    if node isa Number
        return Float64(node)
    elseif node isa Tuple
        length(node) == 3 || throw(ArgumentError("Unsupported OpenQASM expression: $node"))
        left = _eval_param(node[1])
        right = _eval_param(node[3])
        op = node[2].str
        op == "+" && return left + right
        op == "-" && return left - right
        op == "*" && return left * right
        op == "/" && return left / right
        op == "^" && return left^right
        throw(ArgumentError("Unsupported OpenQASM operator `$op`"))
    elseif node isa Neg
        return -_eval_param(node.val)
    elseif node isa Call
        x = _eval_param(node.args)
        fn = node.name
        fn === :sin && return sin(x)
        fn === :cos && return cos(x)
        fn === :tan && return tan(x)
        fn === :exp && return exp(x)
        fn === :ln && return log(x)
        fn === :sqrt && return sqrt(x)
        throw(ArgumentError("Unsupported OpenQASM function `$fn`"))
    else
        s = node.str
        s == "pi" && return Float64(pi)
        return parse(Float64, s)
    end
end

# qelib1 gates that map directly onto an existing Circuit gate (name + qubit count + params).
const _DIRECT = Dict{String,String}(
    "x" => "X", "y" => "Y", "z" => "Z", "h" => "H",
    "s" => "S", "t" => "T", "sx" => "SX",
    "rx" => "RX", "ry" => "RY", "rz" => "RZ",
    "u" => "U", "u3" => "U",
    "cx" => "CNOT", "cy" => "CY", "cz" => "CZ", "swap" => "Swap",
    "ccx" => "CCX",
)

function _reg_qubits(bit::Bit, offsets::Dict{String,Int}, sizes::Dict{String,Int})
    reg = bit.name.str
    haskey(offsets, reg) || throw(ArgumentError("Unknown OpenQASM qreg `$reg`"))
    if bit.address === nothing
        # a whole-register argument, e.g. `h q;`
        return collect(offsets[reg]+1:offsets[reg]+sizes[reg])
    else
        return [offsets[reg] + parse(Int, bit.address.str) + 1]   # Circuit sites are 1-based
    end
end

function _push_instruction!(c::Circuit, inst::Instruction, offsets::Dict{String,Int}, sizes::Dict{String,Int})
    name = lowercase(inst.name)
    params = Float64[_eval_param(carg) for carg in inst.cargs]
    operands = [_reg_qubits(q, offsets, sizes) for q in inst.qargs]

    # OpenQASM broadcasts a gate over whole-register arguments. Single qubits repeat;
    # registers must share the broadcast length.
    width = maximum(length, operands)
    all(o -> length(o) == 1 || length(o) == width, operands) ||
        throw(ArgumentError("OpenQASM register sizes do not match for broadcasting `$(inst.name)`"))

    for k in 1:width
        qubits = [length(o) == 1 ? o[1] : o[k] for o in operands]
        _apply_gate!(c, name, inst.name, qubits, params)
    end
    return c
end

function _apply_gate!(c::Circuit, name::String, raw::AbstractString, qubits::Vector{Int}, params::Vector{Float64})
    # gates needing a small rewrite onto the Circuit's basis set
    if name in ("tdg", "sdg")
        push!(c, "Phase", qubits[1], name == "sdg" ? -pi / 2 : -pi / 4)
    elseif name in ("p", "u1")
        push!(c, "Phase", qubits[1], params[1])
    elseif name == "u2"
        push!(c, "U", qubits[1], pi / 2, params[1], params[2])
    elseif haskey(_DIRECT, name)
        push!(c, _DIRECT[name], qubits..., params...)
    else
        throw(ArgumentError("Unsupported OpenQASM gate `$raw`"))
    end
    return c
end

function parse_qasm(source::AbstractString; max_strings = 2^30, noise_amplitude = 0)
    program = OpenQASM.parse(String(source))::MainProgram

    # first pass: size the circuit and lay out qreg -> contiguous site offsets
    offsets = Dict{String,Int}()
    sizes = Dict{String,Int}()
    n = 0
    for stmt in program.prog
        if stmt isa RegDecl && stmt.type.str == "qreg"
            sz = parse(Int, stmt.size.str)
            offsets[stmt.name.str] = n
            sizes[stmt.name.str] = sz
            n += sz
        end
    end
    n > 0 || throw(ArgumentError("OpenQASM program declares no qreg"))

    c = Circuit(n; max_strings = max_strings, noise_amplitude = noise_amplitude)
    for stmt in program.prog
        if stmt isa Instruction
            _push_instruction!(c, stmt, offsets, sizes)
        end
        # RegDecl (creg/qreg), Include, Measure, Barrier are ignored for a unitary import
    end
    return c
end

load_qasm(path::AbstractString; kwargs...) = parse_qasm(read(path, String); kwargs...)

end # module
