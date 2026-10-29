//! Python bindings for [`synir`], exported as the `synpy_rust` extension module.
//!
//! The Rust side stays thin on purpose: import and synthesis logic lives in `synir` so that its
//! own test suite exercises the exact code path Python drives, rather than a lookalike that can
//! drift. What this crate adds is the pyo3 surface, the Qiskit circuit sink, and conversion
//! between Python values and `synir` types.
//!
//! # What is exported
//!
//! | Python name | Rust type | Purpose |
//! |---|---|---|
//! | `PyPauliExponential` | [`PyPauliExponential`] | Build a circuit gate by gate, then synthesize it |
//! | `PyConnectivity` | [`PyConnectivity`] | A hardware coupling graph, to constrain synthesis |
//! | `PyCliffordTableau` | [`PyCliffordTableau`] | A Clifford operator; synthesizes on its own |
//! | `QiskitSynIR` | [`QiskitSynIR`] | Gate sink writing into a Qiskit `QuantumCircuit` |
//! | `PyCommand` | [`PyCommand`] | One emitted gate, as a Python enum |
//! | `PyPauliString` | [`PyPauliString`] | A Pauli string with a phase |
//! | `synthesize_pauli_exponential` | [`synthesize_pauli_exponential`] | One-shot synthesis from a Hamiltonian |
//!
//! The Python package that wraps this module — including the Qiskit transpiler plugin — lives in
//! `synpy/python/synpy`.

pub mod connectivity;
pub mod synthesis;
pub mod tableau;
pub mod validation;
pub mod wrapper;

use crate::connectivity::PyConnectivity;
use crate::synthesis::synthesize_pauli_exponential;
use crate::synthesis::{PyCommand, PyPauliString};
use crate::tableau::PyCliffordTableau;
use crate::wrapper::qiskit::QiskitSynIR;
use crate::wrapper::PyPauliExponential;
use pyo3::prelude::{PyModule, PyModuleMethods};
use pyo3::{pymodule, wrap_pyfunction, Bound, PyResult};
use synir::ir::{CliffordGates, Gates};

/// Registers every exported class and function on the `synpy_rust` module.
#[pymodule]
#[pyo3(name = "synpy_rust")]
fn synpy_rust(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyPauliString>()?;
    m.add_class::<PyCommand>()?;
    let _ = m.add_function(wrap_pyfunction!(synthesize_pauli_exponential, m)?);
    m.add_class::<PyCliffordTableau>()?;
    m.add_class::<PyPauliExponential>()?;
    m.add_class::<PyConnectivity>()?;
    m.add_class::<QiskitSynIR>()?;
    Ok(())
}
