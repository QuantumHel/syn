mod synthesis;
mod tableau;
mod validation;

use crate::synthesis::synthesize_pauli_exponential;
use crate::synthesis::{PyCommand, PyPauliString};
use crate::tableau::PyCliffordTableau;
use pyo3::prelude::{PyModule, PyModuleMethods};
use pyo3::{pymodule, wrap_pyfunction, Bound, FromPyObject, PyRef, PyResult, Python};
use std::ops::Deref;
use syn::ir::Synthesizer;

#[pymodule]
#[pyo3(name = "synpy_rust")]
fn synpy_rust(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyPauliString>()?;
    m.add_class::<PyCommand>()?;
    let _ = m.add_function(wrap_pyfunction!(synthesize_pauli_exponential, m)?);
    m.add_class::<PyCliffordTableau>()?;

    Ok(())
}