# Quantum circuits

The module `Circuits` provides a way to construct and simulate circuits.
Initialise an empty circuit with `Circuit(n)` where `n` is the number of qubits.
Then, add gates to the circuit with `push!(circuit, gate, qubits...)` where `gate` is a string representing the gate and `qubits` are the qubits the gate acts on.
For example, a CNOT gate with control qubit 1 and target qubit 2 is added to the circuit `c` with `push!(c, "CNOT", 1, 2)`.
Depolarizing noise can be added to the circuit with `push!(c, "Noise")`. The noise amplitude can be set with `c.noise_amplitude=p`.


Lets construct a CCX gate (Toffoli gate) out of CNOT and single qubit gates.


```@example circuits
using PauliStrings
using PauliStrings.Circuits

function noisy_toffoli()
    c = Circuit(3)
    push!(c, "H", 3)
    push!(c, "CNOT", 2, 3); push!(c, "Noise")
    push!(c, "Tdg", 3)
    push!(c, "CNOT", 1, 3); push!(c, "Noise")
    push!(c, "T", 3)
    push!(c, "CNOT", 2, 3); push!(c, "Noise")
    push!(c, "Tdg", 3)
    push!(c, "CNOT", 1, 3); push!(c, "Noise")
    push!(c, "T", 2)
    push!(c, "T", 3)
    push!(c, "CNOT", 1, 2); push!(c, "Noise")
    push!(c, "H", 3)
    push!(c, "T", 1)
    push!(c, "Tdg", 2)
    push!(c, "CNOT", 1, 2); push!(c, "Noise")
    return c
end
```

We can plot the circuit using [QuantumCircuitDraw.jl](https://github.com/nicolasloizeau/QuantumCircuitDraw.jl):
```julia
using QuantumCircuitDraw
c = noisy_toffoli()
paulistrings_plot(c)
```

![plot](./assets/toffoli.png)


The circuit can be compiled to a unitary matrix with `compile(c)`.
Before compiling, set `c.noise_amplitude` to the desired noise amplitude and `c.max_strings` to the number of strings to keep at each step.
Compile will multiply the gates in the circuit from left to right, apply the noise using [`add_noise`](@ref) and trim the operator at each step using [`trim`](@ref).

Lets check our Toffoli gate against the built-in `CCX` gate in `Circuits`:

```@example circuits
using PauliStrings.Circuits
using LinearAlgebra: norm
c = noisy_toffoli()
c.noise_amplitude = 0
U1 = compile(c)
U2 = CCXGate(3,1,2,3)
println(norm(U1-U2))
```

We can also compute expectation values of states in the computational basis with `expect(c, state_out)` and `expect(c, state_in, state_out)` .
Let's compute the expectation values of the output states of the Toffoli gate when the input state is $|111\rangle$:
```@example circuits
in_state = "111"
out_states = ["000", "001", "010", "011", "100", "101", "110", "111"]
p = [real(expect(c, in_state, out_state)) for out_state in out_states]
using Plots
bar(out_states, p, legend=false, xlabel="out state", ylabel="<out|U|in>")
```

Same as above but seting the noise amplitude to 0.05:
```@example circuits
c.noise_amplitude = 0.05
p = [real(expect(c, in_state, out_state)) for out_state in out_states]
bar(out_states, p, legend=false, xlabel="out state", ylabel="<out|U|in>")
```

## Importing OpenQASM circuits

OpenQASM 2 programs can be imported from a string with [`from_openqasm`](@ref) or
from a file with [`from_openqasm_file`](@ref). OpenQASM qubits are zero-indexed;
the importer maps them to the one-indexed qubits used by `PauliStrings.Circuits`.

For example, the small `cat_state_n4.qasm` circuit from
[QASMBench](https://github.com/pnnl/QASMBench/tree/master/small/cat_state_n4)
can be loaded and compiled as follows:

```julia
using PauliStrings.Circuits

circuit = from_openqasm_file(
    "path/to/QASMBench/small/cat_state_n4/cat_state_n4.qasm";
    ignore_measurements=true,
)
unitary = compile(circuit)
```

`Circuit` represents unitary evolution, so measurements and resets are rejected
by default. Pass `ignore_measurements=true` when a benchmark ends with classical
readout and only its unitary portion should be imported. Multiple quantum
registers, register-wide gates, comments, and standard `qelib1.inc` gates are
supported. Custom gate declarations and classical control are reported as
unsupported instead of being silently ignored.
