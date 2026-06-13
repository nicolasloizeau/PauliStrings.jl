struct OpenQASMError <: Exception
    message::String
end

Base.showerror(io::IO, err::OpenQASMError) = print(io, err.message)

function _qasm_error(message, statement = nothing)
    suffix = isnothing(statement) ? "" : ": `$(strip(statement))`"
    throw(OpenQASMError("OpenQASM import error: $message$suffix"))
end

function _qasm_strip_comments(source::AbstractString)
    without_blocks = replace(source, r"/\*.*?\*/"s => " ")
    return replace(without_blocks, r"//[^\n\r]*" => " ")
end

function _qasm_split_arguments(value::AbstractString)
    arguments = String[]
    depth = 0
    start = firstindex(value)

    for index in eachindex(value)
        character = value[index]
        if character == '('
            depth += 1
        elseif character == ')'
            depth -= 1
            depth < 0 && _qasm_error("unbalanced parentheses", value)
        elseif character == ',' && depth == 0
            push!(arguments, strip(value[start:prevind(value, index)]))
            start = nextind(value, index)
        end
    end

    depth == 0 || _qasm_error("unbalanced parentheses", value)
    start <= lastindex(value) && push!(arguments, strip(value[start:end]))
    return filter(!isempty, arguments)
end

function _qasm_eval_expression(expression::AbstractString)
    parsed = try
        Meta.parse(strip(expression))
    catch
        _qasm_error("invalid parameter expression", expression)
    end
    return Float64(_qasm_eval_expression(parsed, expression))
end

function _qasm_eval_expression(value::Real, source)
    return value
end

function _qasm_eval_expression(symbol::Symbol, source)
    symbol == :pi && return pi
    _qasm_error("unsupported symbol `$symbol` in parameter expression", source)
end

function _qasm_eval_expression(expression::Expr, source)
    expression.head == :call || _qasm_error("unsupported parameter expression", source)
    operation = expression.args[1]
    arguments =
        [_qasm_eval_expression(argument, source) for argument in expression.args[2:end]]

    operation == :+ && return +(arguments...)
    operation == :- && return -(arguments...)
    operation == :* && return *(arguments...)
    operation == :/ && length(arguments) == 2 && return arguments[1] / arguments[2]
    operation == :^ && length(arguments) == 2 && return arguments[1] ^ arguments[2]
    operation == :sin && length(arguments) == 1 && return sin(arguments[1])
    operation == :cos && length(arguments) == 1 && return cos(arguments[1])
    operation == :tan && length(arguments) == 1 && return tan(arguments[1])
    operation == :exp && length(arguments) == 1 && return exp(arguments[1])
    operation == :ln && length(arguments) == 1 && return log(arguments[1])
    operation == :sqrt && length(arguments) == 1 && return sqrt(arguments[1])
    _qasm_error("unsupported operation `$operation` in parameter expression", source)
end

function _qasm_eval_expression(value, source)
    _qasm_error("unsupported parameter expression", source)
end

function _qasm_split_invocation(statement::AbstractString)
    match_result = match(r"^([A-Za-z_][A-Za-z0-9_]*)", statement)
    isnothing(match_result) && _qasm_error("invalid gate invocation", statement)

    gate = lowercase(match_result.captures[1])
    position = nextind(statement, match_result.offset, length(match_result.match))
    while position <= lastindex(statement) && isspace(statement[position])
        position = nextind(statement, position)
    end

    parameters = String[]
    if position <= lastindex(statement) && statement[position] == '('
        depth = 1
        parameter_start = nextind(statement, position)
        position = parameter_start
        while position <= lastindex(statement) && depth > 0
            character = statement[position]
            character == '(' && (depth += 1)
            character == ')' && (depth -= 1)
            depth == 0 && break
            position = nextind(statement, position)
        end
        depth == 0 || _qasm_error("unbalanced gate parameters", statement)
        parameter_end = prevind(statement, position)
        if parameter_start <= parameter_end
            parameters = _qasm_split_arguments(statement[parameter_start:parameter_end])
        end
        position = nextind(statement, position)
    end

    while position <= lastindex(statement) && isspace(statement[position])
        position = nextind(statement, position)
    end
    position <= lastindex(statement) || _qasm_error("missing gate operands", statement)
    operands = _qasm_split_arguments(statement[position:end])
    return gate, parameters, operands
end

function _qasm_resolve_operand(
    operand::AbstractString,
    registers::Dict{String,UnitRange{Int}},
)
    indexed = match(r"^([A-Za-z_][A-Za-z0-9_]*)\s*\[\s*(\d+)\s*\]$", operand)
    if !isnothing(indexed)
        name = indexed.captures[1]
        haskey(registers, name) || _qasm_error("unknown quantum register `$name`", operand)
        index = parse(Int, indexed.captures[2])
        register = registers[name]
        0 <= index < length(register) ||
            _qasm_error("qubit index $index is out of bounds", operand)
        return [first(register) + index]
    end

    name_match = match(r"^([A-Za-z_][A-Za-z0-9_]*)$", strip(operand))
    isnothing(name_match) && _qasm_error("invalid quantum operand", operand)
    name = name_match.captures[1]
    haskey(registers, name) || _qasm_error("unknown quantum register `$name`", operand)
    return collect(registers[name])
end

function _qasm_require_arity(gate, sites, parameters, site_count, parameter_count)
    length(sites) == site_count ||
        _qasm_error("gate `$gate` expects $site_count qubit operand(s)")
    length(parameters) == parameter_count ||
        _qasm_error("gate `$gate` expects $parameter_count parameter(s)")
end

function _qasm_append_gate!(circuit::Circuit, gate, sites, parameters)
    if gate == "id"
        _qasm_require_arity(gate, sites, parameters, 1, 0)
    elseif gate in ("x", "y", "z", "h", "s", "t", "sx")
        _qasm_require_arity(gate, sites, parameters, 1, 0)
        push!(circuit, uppercase(gate), sites[1])
    elseif gate == "sdg"
        _qasm_require_arity(gate, sites, parameters, 1, 0)
        push!(circuit, "Phase", sites[1], -pi / 2)
    elseif gate == "tdg"
        _qasm_require_arity(gate, sites, parameters, 1, 0)
        push!(circuit, "Tdg", sites[1])
    elseif gate == "sxdg"
        _qasm_require_arity(gate, sites, parameters, 1, 0)
        for _ = 1:3
            push!(circuit, "SX", sites[1])
        end
    elseif gate in ("rx", "ry", "rz")
        _qasm_require_arity(gate, sites, parameters, 1, 1)
        push!(circuit, uppercase(gate), sites[1], parameters[1])
    elseif gate in ("p", "phase", "u1")
        _qasm_require_arity(gate, sites, parameters, 1, 1)
        push!(circuit, "Phase", sites[1], parameters[1])
    elseif gate == "u2"
        _qasm_require_arity(gate, sites, parameters, 1, 2)
        push!(circuit, "U", sites[1], pi / 2, parameters...)
    elseif gate in ("u", "u3")
        _qasm_require_arity(gate, sites, parameters, 1, 3)
        push!(circuit, "U", sites[1], parameters...)
    elseif gate in ("cx", "cnot")
        _qasm_require_arity(gate, sites, parameters, 2, 0)
        push!(circuit, "CNOT", sites...)
    elseif gate in ("cy", "cz", "swap", "csx", "csxdg")
        _qasm_require_arity(gate, sites, parameters, 2, 0)
        mapped = Dict(
            "cy" => "CY",
            "cz" => "CZ",
            "swap" => "Swap",
            "csx" => "CSX",
            "csxdg" => "CSXdg",
        )
        push!(circuit, mapped[gate], sites...)
    elseif gate in ("cp", "cu1")
        _qasm_require_arity(gate, sites, parameters, 2, 1)
        push!(circuit, "CPhase", sites..., parameters[1])
    elseif gate == "ccx"
        _qasm_require_arity(gate, sites, parameters, 3, 0)
        push!(circuit, "CCX", sites...)
    else
        _qasm_error("unsupported gate `$gate`")
    end
    return circuit
end

"""
    from_openqasm(source; ignore_measurements=false, max_strings=2^30,
                  noise_amplitude=0)

Import an OpenQASM 2 program into a [`Circuit`](@ref). Qubit indices are converted
from OpenQASM's zero-based convention to the one-based convention used by
`PauliStrings.Circuits`. Register-wide gate applications and multiple quantum
registers are supported.

The `Circuit` type represents unitary circuits only. Set `ignore_measurements=true`
to load benchmark programs that end in measurement statements; otherwise
`measure` and `reset` statements raise an `OpenQASMError`.
"""
function from_openqasm(
    source::AbstractString;
    ignore_measurements::Bool = false,
    max_strings = 2^30,
    noise_amplitude = 0,
)
    cleaned = _qasm_strip_comments(source)
    occursin(r"(?i)\bgate\s+[A-Za-z_]", cleaned) &&
        _qasm_error("custom gate declarations are not supported")
    statements = filter(!isempty, strip.(split(cleaned, ';')))

    registers = Dict{String,UnitRange{Int}}()
    gate_statements = String[]
    qubit_count = 0
    saw_header = false

    for statement in statements
        if occursin(r"(?i)^OPENQASM\s+2(?:\.0)?$", statement)
            saw_header = true
        elseif occursin(r"(?i)^include\s+\"qelib1\.inc\"$", statement)
            continue
        elseif (
            register_match = match(
                r"(?i)^qreg\s+([A-Za-z_][A-Za-z0-9_]*)\s*\[\s*(\d+)\s*\]$",
                statement,
            )
        ) !== nothing
            name = register_match.captures[1]
            haskey(registers, name) && _qasm_error("duplicate quantum register `$name`")
            size = parse(Int, register_match.captures[2])
            size > 0 || _qasm_error("quantum registers must contain at least one qubit")
            registers[name] = (qubit_count+1):(qubit_count+size)
            qubit_count += size
        elseif occursin(r"(?i)^creg\s+", statement) ||
               occursin(r"(?i)^barrier\b", statement)
            continue
        elseif occursin(r"(?i)^(measure|reset)\b", statement)
            ignore_measurements || _qasm_error(
                "non-unitary statements require `ignore_measurements=true`",
                statement,
            )
        elseif occursin(r"(?i)^(opaque|if)\b", statement)
            _qasm_error(
                "classical control and opaque declarations are not supported",
                statement,
            )
        else
            push!(gate_statements, statement)
        end
    end

    saw_header || _qasm_error("missing `OPENQASM 2.0` header")
    qubit_count > 0 || _qasm_error("program does not declare a quantum register")
    circuit =
        Circuit(qubit_count; max_strings = max_strings, noise_amplitude = noise_amplitude)

    for statement in gate_statements
        gate, parameter_sources, operand_sources = _qasm_split_invocation(statement)
        parameters = _qasm_eval_expression.(parameter_sources)
        operands =
            [_qasm_resolve_operand(operand, registers) for operand in operand_sources]
        width = maximum(length, operands)
        all(length(operand) in (1, width) for operand in operands) ||
            _qasm_error("register operands must have matching sizes", statement)

        for index = 1:width
            sites =
                [length(operand) == 1 ? operand[1] : operand[index] for operand in operands]
            _qasm_append_gate!(circuit, gate, sites, parameters)
        end
    end
    return circuit
end

"""
    from_openqasm_file(path; kwargs...)

Read an OpenQASM 2 program from `path` and import it with [`from_openqasm`](@ref).
"""
function from_openqasm_file(path::AbstractString; kwargs...)
    return from_openqasm(read(path, String); kwargs...)
end
