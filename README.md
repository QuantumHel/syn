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