//! Gate commands, the collector sink, and the one-shot synthesis entry point.

use pyo3::exceptions::PyException;
use pyo3::{pyclass, pymethods, PyErr};
use std::ops::Deref;
use synir::data_structures::Angle;

use pyo3::{pyfunction, PyRef, PyResult};
use synir::{
    data_structures::{CliffordTableau, PauliExponential, PauliPolynomial, PropagateClifford},
    ir::{
        clifford_tableau::CliffordTableauSynthStrategy,
        pauli_exponential::PauliExponentialSynthesizer,
        pauli_polynomial::PauliPolynomialSynthStrategy, CliffordGates, Gates, Synthesizer,
    },
    IndexType,
};

use crate::validation::validate;

/// A single emitted gate, carrying its qubit indices and — for rotations — its angle in radians.
///
/// This is what synthesis hands back to Python: a flat list in circuit order. Convert it with
/// `synpy.utils.pycommand_to_tuple` or `synpy.utils.pycommand_to_qasm`.
#[pyclass]
#[derive(Debug, Copy, Clone)]
pub enum PyCommand {
    /// Controlled-X, `(control, target)`.
    CX(usize, usize),
    /// Controlled-Z, `(control, target)`.
    CZ(usize, usize),
    /// Pauli X.
    X(usize),
    /// Pauli Y.
    Y(usize),
    /// Pauli Z.
    Z(usize),
    /// Hadamard.
    H(usize),
    /// Phase gate `S`.
    S(usize),
    /// `V = sqrt(X)`.
    V(usize),
    /// `S^dag`.
    SDgr(usize),
    /// `V^dag`.
    VDgr(usize),
    /// X rotation, `(target, angle_in_radians)`.
    Rx(usize, f64),
    /// Y rotation, `(target, angle_in_radians)`.
    Ry(usize, f64),
    /// Z rotation, `(target, angle_in_radians)`.
    Rz(usize, f64),
}

impl PyCommand {
    /// The smallest register this command fits in — its highest qubit index plus one.
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

/// Replays a list of Clifford commands into a tableau.
///
/// # Errors
///
/// Returns a Python `Exception` if the list contains a rotation, which has no tableau
/// representation.
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

/// A wrapper type that can synthesize itself into a command list.
pub trait Synthesize {
    /// The gates realizing this object, in circuit order.
    fn synthesize(&self) -> Vec<PyCommand>;
}

/// A Pauli string together with a rotation angle — one gadget of a Hamiltonian.
///
/// `pauli_string` is one letter per qubit, e.g. `"ZZII"`; `phase` is the angle in radians.
#[pyclass]
#[derive(Clone)]
pub struct PyPauliString {
    /// One letter from `IXYZ` per qubit.
    pub pauli_string: String,
    /// The rotation angle, in radians.
    pub phase: f64,
}

#[pymethods]
impl PyPauliString {
    /// Builds a gadget from its string and angle.
    #[new]
    fn new(pauli_string: String, phase: f64) -> Self {
        PyPauliString {
            pauli_string,
            phase,
        }
    }

    /// The number of qubits, i.e. the length of the string.
    pub fn get_qubits(&self) -> usize {
        self.pauli_string.len()
    }
    /// This gadget as a `(pauli_string, phase)` tuple.
    pub fn as_tuple(&self) -> (&str, f64) {
        (self.pauli_string.as_str(), self.phase)
    }
}

/// A gate sink that records everything emitted into it as [`PyCommand`]s.
///
/// This is the bridge from `synir`'s [`CliffordGates`]/[`Gates`] traits to a Python-visible list.
/// Unlike [`QiskitSynIR`](crate::wrapper::qiskit::QiskitSynIR), it has no opinion about the target
/// framework — it just buffers.
///
/// It inherits the default `add_final_permutation`, so a leftover permutation arrives expanded
/// into `CX` triples rather than reported separately.
#[derive(Debug, Default)]
pub struct CommandCollector {
    commands: Vec<PyCommand>,
}

impl CommandCollector {
    /// An empty collector.
    pub fn new() -> Self {
        Self {
            commands: Vec::new(),
        }
    }

    /// A copy of everything collected so far, in emission order.
    pub fn commands(&self) -> Vec<PyCommand> {
        self.commands.clone()
    }
}

impl CliffordGates for CommandCollector {
    fn s(&mut self, target: IndexType) {
        self.commands.push(PyCommand::S(target));
    }

    fn v(&mut self, target: IndexType) {
        self.commands.push(PyCommand::V(target));
    }

    fn s_dgr(&mut self, target: IndexType) {
        self.commands.push(PyCommand::SDgr(target));
    }

    fn v_dgr(&mut self, target: IndexType) {
        self.commands.push(PyCommand::VDgr(target));
    }

    fn x(&mut self, target: IndexType) {
        self.commands.push(PyCommand::X(target));
    }

    fn y(&mut self, target: IndexType) {
        self.commands.push(PyCommand::Y(target));
    }

    fn z(&mut self, target: IndexType) {
        self.commands.push(PyCommand::Z(target));
    }

    fn h(&mut self, target: IndexType) {
        self.commands.push(PyCommand::H(target));
    }

    fn cx(&mut self, control: IndexType, target: IndexType) {
        self.commands.push(PyCommand::CX(control, target));
    }

    fn cz(&mut self, control: IndexType, target: IndexType) {
        self.commands.push(PyCommand::CZ(control, target));
    }
}

impl Gates for CommandCollector {
    fn rx(&mut self, target: IndexType, angle: f64) {
        self.commands.push(PyCommand::Rx(target, angle));
    }

    fn ry(&mut self, target: IndexType, angle: f64) {
        self.commands.push(PyCommand::Ry(target, angle));
    }

    fn rz(&mut self, target: IndexType, angle: f64) {
        self.commands.push(PyCommand::Rz(target, angle));
    }
}

/// Synthesizes a Hamiltonian plus a trailing Clifford into a flat command list.
///
/// The one-shot alternative to building up a
/// [`PyPauliExponential`](crate::wrapper::PyPauliExponential) gate by gate — it takes the IR's two
/// halves directly. `hamiltonian` is a list of polynomials, each a list of gadgets;
/// `clifford_gates` is replayed into the trailing tableau.
///
/// Uses the naive strategy for both halves, so the result assumes all-to-all connectivity.
///
/// # Errors
///
/// Returns a Python `Exception` if the arguments fail [`validate`], or if `clifford_gates`
/// contains a rotation.
#[pyfunction]
#[pyo3(signature = (hamiltonian, clifford_gates, nr_qubits), text_signature = "(hamiltonian: list[list[PyPauliString]], clifford_gates: list[PyCommand], nr_qubits: int)")]
pub fn synthesize_pauli_exponential(
    hamiltonian: Vec<Vec<PyPauliString>>,
    clifford_gates: Vec<PyRef<PyCommand>>,
    nr_qubits: usize,
) -> PyResult<Vec<PyCommand>> {
    let converted_hamiltonian = hamiltonian
        .iter()
        .map(|inner| {
            inner
                .iter()
                .map(PyPauliString::as_tuple)
                .map(|(pauli_str, phase)| (pauli_str, Angle::from_angle(phase)))
                .collect::<Vec<(&str, Angle)>>()
        })
        .map(|inner: Vec<(&str, Angle)>| PauliPolynomial::from_hamiltonian(inner))
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
