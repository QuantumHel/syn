//! The `PyPauliExponential` builder and the synthesis entry point behind it.

pub mod qiskit;

extern crate pyo3;
extern crate pyo3_ffi;


use pyo3::prelude::*;
use synir::{
    data_structures::PauliExponential,
    import::CircuitImporter,
    ir::{
        clifford_tableau::CliffordTableauSynthStrategy,
        pauli_exponential::PauliExponentialSynthesizer,
        pauli_polynomial::PauliPolynomialSynthStrategy, CliffordGates, Gates, Synthesizer,
    },
};

use crate::wrapper::qiskit::QiskitSynIR;

/// Builds a circuit gate by gate, then synthesizes it.
///
/// Gates must be added in **reverse circuit order** — last executed first. That is what
/// `synpy.qiskit.plugin.qiskit_to_synir` does with `reversed(circuit.data)`, and it is what lets
/// each Clifford conjugate the gadgets that come after it as it arrives.
///
/// Clifford gates never become gadgets: they are folded into the tableau immediately. `add_rz`
/// checks its angle and does the same whenever the rotation turns out to be Clifford, so only
/// genuinely non-Clifford rotations survive as gadgets.
///
/// Calling `synthesize_to_qiskit` **consumes** the accumulated gates; the object is left empty and
/// ready to be built up again.
///
/// ```python
/// from synpy.synpy_rust import PyPauliExponential, QiskitSynIR
///
/// pe = PyPauliExponential(2)
/// pe.set_pauli_strategy("PSGS")
/// for gate in reversed(circuit.data):
///     ...          # pe.add_cx(...), pe.add_h(...), pe.add_rz(...)
///
/// out = QiskitSynIR(2)
/// pe.synthesize_to_qiskit(out)
/// result = out.get_circuit()
/// ```
#[pyclass]
pub struct PyPauliExponential {
    /// The gate-by-gate import lives in `synir::import` so that `synir`'s own tests exercise the
    /// same code path this wrapper does. `qiskit_to_synir` feeds gates in reverse circuit order.
    importer: CircuitImporter,
    pauli_strategy: PauliPolynomialSynthStrategy,
    tableau_strategy: CliffordTableauSynthStrategy,
}

#[pymethods]
impl PyPauliExponential {
    /// An empty builder over `num_qubits` qubits.
    ///
    /// Defaults to the `Naive` polynomial strategy and the `PermRowCol` tableau strategy.
    #[new]
    pub fn new(num_qubits: usize) -> Self {
        Self {
            importer: CircuitImporter::new(num_qubits),
            pauli_strategy: PauliPolynomialSynthStrategy::Naive,
            tableau_strategy: CliffordTableauSynthStrategy::PermRowCol,
        }
    }

    /// Synthesizes everything accumulated so far into `circuit`, and resets this builder.
    pub fn synthesize_to_qiskit(&mut self, circuit: &mut QiskitSynIR) {
        synthesize(self, circuit);
    }

    /// Chooses the Pauli polynomial algorithm: `"Naive"` or `"PSGS"`.
    ///
    /// Raises `pyo3_runtime.PanicException` on any other name.
    pub fn set_pauli_strategy(&mut self, strategy: String) {
        match strategy.as_str() {
            "Naive" => self.pauli_strategy = PauliPolynomialSynthStrategy::Naive,
            "PSGS" => self.pauli_strategy = PauliPolynomialSynthStrategy::PSGS,
            _ => panic!("Unknown Pauli polynomial synthesis strategy: {}", strategy),
        }
    }
    /// Chooses the Clifford tableau algorithm: `"Naive"` or `"PermRowCol"`.
    ///
    /// Raises `pyo3_runtime.PanicException` on any other name.
    pub fn set_tableau_strategy(&mut self, strategy: String) {
        match strategy.as_str() {
            "Naive" => self.tableau_strategy = CliffordTableauSynthStrategy::Naive,
            "PermRowCol" => self.tableau_strategy = CliffordTableauSynthStrategy::PermRowCol,
            _ => panic!("Unknown Clifford tableau synthesis strategy: {}", strategy),
        }
    }

    /// Adds a Hadamard on `target`.
    pub fn add_h(&mut self, target: usize) {
        self.importer.add_h(target);
    }

    /// Adds `V = sqrt(X)` on `target`.
    pub fn add_v(&mut self, target: usize) {
        self.importer.add_v(target);
    }

    /// Adds the phase gate `S` on `target`.
    pub fn add_s(&mut self, target: usize) {
        self.importer.add_s(target);
    }

    /// Adds `S^dag` on `target`.
    pub fn add_s_dgr(&mut self, target: usize) {
        self.importer.add_s_dgr(target);
    }

    /// Adds Pauli `X` on `target`.
    pub fn add_x(&mut self, target: usize) {
        self.importer.add_x(target);
    }

    /// Adds Pauli `Y` on `target`.
    pub fn add_y(&mut self, target: usize) {
        self.importer.add_y(target);
    }

    /// Adds Pauli `Z` on `target`.
    pub fn add_z(&mut self, target: usize) {
        self.importer.add_z(target);
    }

    /// Adds a `CX` from `control` onto `target`.
    pub fn add_cx(&mut self, control: usize, target: usize) {
        self.importer.add_cx(control, target);
    }

    /// Adds a `Z` rotation of `angle` radians on `target`.
    ///
    /// A Clifford angle (a multiple of `pi/2`) is folded into the tableau as the equivalent
    /// Clifford gate rather than kept as a gadget; an angle of 0 adds nothing at all.
    pub fn add_rz(&mut self, target: usize, angle: f64) {
        self.importer.add_rz(target, angle);
    }

    /// Prints the gadgets and tableau accumulated so far, for debugging.
    ///
    /// The tableau shown is still in the inverse form the importer keeps it in; it is only
    /// inverted when synthesis consumes the builder.
    pub fn print(&self) {
        pe_print_helper(self.importer.raw());
    }
}

/// Runs the preprocessing passes and then synthesis, emitting into `circuit`.
///
/// Takes the accumulated gates out of `pe`, leaving it empty over the same number of qubits.
pub fn synthesize<G>(pe: &mut PyPauliExponential, circuit: &mut G)
where
    G: CliffordGates + Gates,
{
    let mut synth = PauliExponentialSynthesizer::from_strategy(
        pe.pauli_strategy.clone(),
        pe.tableau_strategy.clone(),
    );
    let num_qubits = pe.importer.num_qubits();
    let importer = std::mem::replace(&mut pe.importer, CircuitImporter::new(num_qubits));
    let pe = importer.finish();
    synth.synthesize(pe, circuit);
}

/// Prints a Pauli exponential's angles, chains and tableau to stdout.
fn pe_print_helper(pe: &PauliExponential) {
    println!("PE:");
    for pp in pe.chains() {
        for a in pp.angles() {
            print!("{} ", *a);
        }
        println!("");
        for chain in pp.chains() {
            println!("{}", chain);
        }
        println!("---");
    }
    println!("{}", pe.clifford_tableau());
}
