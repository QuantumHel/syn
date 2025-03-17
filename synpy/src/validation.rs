use crate::synthesis::{PyCommand, PyPauliString};
use pyo3::exceptions::PyException;
use pyo3::PyErr;

pub fn validate(
    converted_hamiltonian: &Vec<Vec<PyPauliString>>,
    clifford_gates: &Vec<PyCommand>,
    nr_qubits: usize,
) -> Result<(), PyErr> {
    let qubits: Vec<usize> = converted_hamiltonian
        .iter()
        .map(|inner| {
            inner
                .iter()
                .map(PyPauliString::get_qubits)
                .collect::<Vec<usize>>()
        })
        .flatten()
        .collect();

    if let Some(&first) = qubits.first() {
        if !qubits.iter().all(|&x| x == first) {
            return Err(PyException::new_err(
                "All Paulistrings provided to the hamiltonian must be of the same length.",
            ));
        }

        if first != nr_qubits {
            return Err(PyException::new_err(
                "Hamiltonian has invalid number of qubits.",
            ));
        }
    } else {
        return Err(PyException::new_err("Hamiltonian has no pauli strings."));
    }
    if let Some(nr_qubits_clifford) = clifford_gates
        .iter()
        .map(PyCommand::get_max_nr_qubits)
        .max()
    {
        if nr_qubits_clifford >= nr_qubits {
            return Err(PyException::new_err("Clifford have too many qubits."));
        }
    }
    Ok(())
}
