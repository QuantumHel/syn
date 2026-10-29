# synIR

An intermediate representation for quantum circuits, together with architecture-aware passes
that synthesize it back into gates.

A circuit is rewritten into the normal form `U = T * U_pp` (matrix order, so `U_pp` runs first):
a run of Pauli gadgets followed by a single Clifford operator. Every Clifford in the input is
absorbed into `T`, leaving only genuinely non-Clifford rotations in `U_pp`; synthesis is then free
to choose new Clifford structure for both halves, subject to a hardware coupling graph.

The workspace has two crates:

- **`synir`** — the IR and the synthesis passes.
- **`synpy`** — Python bindings (the `synpy_rust` extension module), plus the `synpy` Python
  package with Qiskit and pauliopt integrations.

## Constraining synthesis to a coupling graph

Synthesis can target a hardware graph directly, rather than emitting an all-to-all circuit and
paying a routing pass afterwards:

```python
from synpy.qiskit.connectivity import target_to_connectivity
from synpy.qiskit.plugin import qiskit_to_synir
from synpy.synpy_rust import PyConnectivity, QiskitSynIR

pe = qiskit_to_synir(circuit)
pe.set_pauli_strategy("PSGS")          # the connectivity-aware gadget synthesizer
pe.set_tableau_strategy("PermRowCol")  # ... and the connectivity-aware Clifford one

sink = QiskitSynIR(circuit.copy_empty_like())
pe.synthesize_to_qiskit(sink, target_to_connectivity(backend.target))
result = sink.get_circuit()
```

Every two-qubit gate in `result` then lies on an edge of the graph. Omit the second argument for
all-to-all. `PyConnectivity` also has `line`, `grid`, `complete` and a raw `(num_qubits, edges)`
constructor, so no Qiskit backend is needed.

Conformance takes *both* architecture-aware strategies: a `"Naive"` half routes wherever it likes.
Pairing a graph with one `"Naive"` half warns, with two it raises.

## Documentation

Build and open the Rust API docs:

```
cargo doc --open -p synir
```

Start at the `synir` crate root, which covers the normal form, the circuit-order versus
matrix-order conventions, and the module map.

Python users get the same documentation through `help()` — pyo3 installs the Rust doc comments as
runtime docstrings — and type information from `synpy/python/synpy/synpy_rust.pyi`.

## Installation

### Python

To run the basic python setup, run:

```
make python_project_setup
```

This will install the required dependencies.
*It is recommended to do this in a virtual environment*.

### Rust

To run the basic rust setup, run:

```
make rust_project_setup
```

Build the project:

```
make build
```

Run tests:

```
make test
```

Format-Check:

```
make format-check
```

(You can run format and check individually)

Clean:

```
make clean
```