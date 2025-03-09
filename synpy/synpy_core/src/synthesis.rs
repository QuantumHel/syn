use std::ops::Deref;
use pyo3::{pyfunction, Py};
use crate::PyRef;
use crate::Python;
use pyo3::exceptions::PyException;
use pyo3::{pyclass, pymethods, PyErr};
use syn::data_structures::{CliffordTableau, PauliPolynomial};
use syn::data_structures::PropagateClifford;
use syn::ir::clifford_tableau::CliffordTableauSynthStrategy;
use syn::ir::CliffordGates;
use syn::ir::Gates;
use syn::ir::pauli_exponential::{PauliExponential, PauliExponentialSynthesizer};
use syn::ir::pauli_polynomial::PauliPolynomialSynthStrategy;
use crate::PyResult;
use crate::validation::validate;
use syn::ir::Synthesizer;

#[pyclass]
#[derive(Debug, Copy, Clone)]
pub enum PyCommand {
    CX(usize, usize),
    CZ(usize, usize),
    X(usize),
    Y(usize),
    Z(usize),
    H(usize),
    S(usize),
    V(usize),
    SDgr(usize),
    VDgr(usize),
    Rx(usize, f64),
    Ry(usize, f64),
    Rz(usize, f64),
}

impl PyCommand {
    pub(crate) fn get_max_nr_qubits(&self) -> usize {
        match *self {
            PyCommand::CX(a, b) | PyCommand::CZ(a, b) => std::cmp::max(a, b) + 1,
            PyCommand::X(q)
            | PyCommand::Y(q)
            | PyCommand::Z(q)
            | PyCommand::H(q)
            | PyCommand::S(q)
            | PyCommand::V(q)
            | PyCommand::SDgr(q)
            | PyCommand::VDgr(q)
            | PyCommand::Rx(q, _)
            | PyCommand::Ry(q, _)
            | PyCommand::Rz(q, _) => q + 1,
        }
    }
}

pub fn parse_clifford_commands(
    size: usize,
    commands: &[PyCommand],
) -> Result<CliffordTableau, PyErr> {
    let mut tableau = CliffordTableau::new(size);
    for command in commands.iter() {
        match command {
            PyCommand::H(target) => {
                tableau.h(*target);
            }
            PyCommand::S(target) => {
                tableau.s(*target);
            }
            PyCommand::V(target) => {
                tableau.v(*target);
            }
            PyCommand::CX(control, target) => {
                tableau.cx(*control, *target);
            }
            PyCommand::CZ(control, target) => {
                tableau.cz(*control, *target);
            }
            PyCommand::X(target) => {
                tableau.x(*target);
            }
            PyCommand::Z(target) => {
                tableau.z(*target);
            }
            PyCommand::Y(target) => {
                tableau.y(*target);
            }
            PyCommand::VDgr(target) => {
                tableau.v_dgr(*target);
            }
            PyCommand::SDgr(target) => {
                tableau.s_dgr(*target);
            }
            _ => {
                return Err(PyException::new_err(
                    "Clifford Circuit contained non clifford command",
                ));
            }
        }
    }
    Ok(tableau)
}

#[pyclass]
#[derive(Clone)]
pub struct PyPauliString {
    pub pauli_string: String,
    pub phase: f64,
}

#[pymethods]
impl PyPauliString {
    #[new]
    fn new(pauli_string: String, phase: f64) -> Self {
        PyPauliString {
            pauli_string,
            phase,
        }
    }

    pub fn get_qubits(&self) -> usize {
        self.pauli_string.len()
    }
    pub fn as_tuple(&self) -> (&str, f64) {
        (self.pauli_string.as_str(), self.phase)
    }
}

#[derive(Debug, Default)]
pub struct CommandCollector {
    commands: Vec<PyCommand>,
}

impl CommandCollector {
    pub fn new() -> Self {
        Self {
            commands: Vec::new(),
        }
    }

    pub fn commands(&self) -> Vec<PyCommand>  {
        self.commands.clone()
    }
}

impl CliffordGates for CommandCollector {
    fn s(&mut self, target: syn::IndexType) {
        self.commands.push(PyCommand::S(target));
    }

    fn v(&mut self, target: syn::IndexType) {
        self.commands.push(PyCommand::V(target));
    }

    fn s_dgr(&mut self, target: syn::IndexType) {
        self.commands.push(PyCommand::SDgr(target));
    }

    fn v_dgr(&mut self, target: syn::IndexType) {
        self.commands.push(PyCommand::VDgr(target));
    }

    fn x(&mut self, target: syn::IndexType) {
        self.commands.push(PyCommand::X(target));
    }

    fn y(&mut self, target: syn::IndexType) {
        self.commands.push(PyCommand::Y(target));
    }

    fn z(&mut self, target: syn::IndexType) {
        self.commands.push(PyCommand::Z(target));
    }

    fn h(&mut self, target: syn::IndexType) {
        self.commands.push(PyCommand::H(target));
    }

    fn cx(&mut self, control: syn::IndexType, target: syn::IndexType) {
        self.commands.push(PyCommand::CX(control, target));
    }

    fn cz(&mut self, control: syn::IndexType, target: syn::IndexType) {
        self.commands.push(PyCommand::CZ(control, target));
    }
}

impl Gates for CommandCollector {
    fn rx(&mut self, target: syn::IndexType, angle: f64) {
        self.commands.push(PyCommand::Rx(target, angle));
    }

    fn ry(&mut self, target: syn::IndexType, angle: f64) {
        self.commands.push(PyCommand::Ry(target, angle));
    }

    fn rz(&mut self, target: syn::IndexType, angle: f64) {
        self.commands.push(PyCommand::Rz(target, angle));
    }
}

#[pyfunction]
#[pyo3(signature = (hamiltonian, clifford_gates, nr_qubits), text_signature = "(hamiltonian: list[list[PyPauliString]], clifford_gates: list[PyCommand], nr_qubits: int)")]
pub fn synthesize_pauli_exponential(
    hamiltonian: Vec<Vec<PyPauliString>>,
    clifford_gates: Vec<PyRef<PyCommand>>,
    nr_qubits: usize,
) -> PyResult<Vec<PyCommand>> {
    let converted_hamiltonian = hamiltonian
        .iter()
        .map(|inner| inner.iter().map(PyPauliString::as_tuple).collect())
        .map(|inner: Vec<(&str, f64)>| PauliPolynomial::from_hamiltonian(inner))
        .collect();
    let clifford_gates: Vec<PyCommand> = clifford_gates.iter().map(|cmd| *(cmd.deref())).collect();

    validate(&hamiltonian, &clifford_gates, nr_qubits)?;

    let clifford_tableau = parse_clifford_commands(nr_qubits, &clifford_gates)?;
    let pauli_exponential = PauliExponential::new(converted_hamiltonian, clifford_tableau);

    let mut command_collector = CommandCollector::new();
    let mut synthesizer = PauliExponentialSynthesizer::from_strategy(
        PauliPolynomialSynthStrategy::Naive,
        CliffordTableauSynthStrategy::Naive,
    );
    synthesizer.synthesize(pauli_exponential, &mut command_collector);
    Ok(command_collector.commands())
}
