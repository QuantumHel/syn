mod synthesis;
mod tableau;
mod validation;
mod qiskit;

use crate::synthesis::synthesize_pauli_exponential;
use crate::synthesis::{PyCommand, PyPauliString};
use crate::tableau::PyCliffordTableau;
use crate::qiskit::{PauliExponentialWrap, QiskitSynIR, synthesize_to_qiskit};
use pyo3::prelude::{PyModule, PyModuleMethods};
use pyo3::{pymodule, wrap_pyfunction, Bound, PyResult};

#[pymodule]
#[pyo3(name = "synpy_rust")]
fn synpy_rust(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyPauliString>()?;
    m.add_class::<PyCommand>()?;
    m.add_function(wrap_pyfunction!(synthesize_pauli_exponential, m)?)?;
    m.add_class::<PyCliffordTableau>()?;
    m.add_class::<QiskitSynIR>()?;
    m.add_class::<PauliExponentialWrap>()?;
    m.add_function(wrap_pyfunction!(synthesize_to_qiskit, m)?)?;

    Ok(())
}
