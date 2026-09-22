pub(crate) mod qiskit;

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
    #[new]
    pub fn new(num_qubits: usize) -> Self {
        Self {
            importer: CircuitImporter::new(num_qubits),
            pauli_strategy: PauliPolynomialSynthStrategy::Naive,
            tableau_strategy: CliffordTableauSynthStrategy::PermRowCol,
        }
    }

    pub fn synthesize_to_qiskit(&mut self, circuit: &mut QiskitSynIR) {
        synthesize(self, circuit);
    }

    pub fn set_pauli_strategy(&mut self, strategy: String) {
        match strategy.as_str() {
            "Naive" => self.pauli_strategy = PauliPolynomialSynthStrategy::Naive,
            "PSGS" => self.pauli_strategy = PauliPolynomialSynthStrategy::PSGS,
            _ => panic!("Unknown Pauli polynomial synthesis strategy: {}", strategy),
        }
    }
    pub fn set_tableau_strategy(&mut self, strategy: String) {
        match strategy.as_str() {
            "Naive" => self.tableau_strategy = CliffordTableauSynthStrategy::Naive,
            "PermRowCol" => self.tableau_strategy = CliffordTableauSynthStrategy::PermRowCol,
            _ => panic!("Unknown Clifford tableau synthesis strategy: {}", strategy),
        }
    }

    pub fn add_h(&mut self, target: usize) {
        self.importer.add_h(target);
    }

    pub fn add_v(&mut self, target: usize) {
        self.importer.add_v(target);
    }

    pub fn add_s(&mut self, target: usize) {
        self.importer.add_s(target);
    }

    pub fn add_s_dgr(&mut self, target: usize) {
        self.importer.add_s_dgr(target);
    }

    pub fn add_x(&mut self, target: usize) {
        self.importer.add_x(target);
    }

    pub fn add_y(&mut self, target: usize) {
        self.importer.add_y(target);
    }

    pub fn add_z(&mut self, target: usize) {
        self.importer.add_z(target);
    }

    pub fn add_cx(&mut self, control: usize, target: usize) {
        self.importer.add_cx(control, target);
    }

    pub fn add_rz(&mut self, target: usize, angle: f64) {
        self.importer.add_rz(target, angle);
    }

    pub fn print(&self) {
        pe_print_helper(self.importer.raw());
    }
}

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
