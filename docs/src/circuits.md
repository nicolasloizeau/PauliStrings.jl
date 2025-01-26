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

![plot](./toffoli.png)


The circuit can be compiled to a unitary matrix with `compile(c)`.
Before compiling, set `c.noise_amplitude` to the desired noise amplitude and `c.max_strings` to the number of strings to keep at each step.
Compile will multiply the gates in the circuit from left to right, apply the noise using [`add_noise`](@ref) and trim the operator at each step using [`trim`](@ref).

Lets check our Toffoli gate against the built-in `CCX` gate in `Circuits`:

```@example circuits
using PauliStrings.Circuits
c = noisy_toffoli()
c.noise_amplitude = 0
U1 = compile(c)
U2 = CCXGate(3,1,2,3)
println(opnorm(U1-U2))
```

We can also compute expectation values of states in the computational basis with `expect(c, state_out)` and `expect(c, state_in, state_out)` .
Let's compute the expectation values of the output states of the Toffoli gate when the input state is $|111\rangle$:
```@example circuits
in_state = "111"
out_states = ["000", "001", "010", "011", "100", "101", "110", "111"]
p = [expect(c, in_state, out_state) for out_state in out_states]
using Plots
bar(out_states, p, legend=false, xlabel="out state", ylabel="<out|U|in>")
```

Same as above but seting the noise amplitude to 0.05:
```@example circuits
c.noise_amplitude = 0.05
p = [expect(c, in_state, out_state) for out_state in out_states]
bar(out_states, p, legend=false, xlabel="out state", ylabel="<out|U|in>")
```
