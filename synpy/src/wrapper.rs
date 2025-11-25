pub(crate) mod qiskit;

extern crate pyo3;
extern crate pyo3_ffi;

use std::collections::VecDeque;

use pyo3::prelude::*;
use synir::{
    data_structures::{CliffordTableau, PauliExponential},
    ir::{
        CliffordGates, Gates, Synthesizer, clifford_tableau::CliffordTableauSynthStrategy, pauli_exponential::PauliExponentialSynthesizer, pauli_polynomial::PauliPolynomialSynthStrategy
    },
};

use crate::wrapper::qiskit::QiskitSynIR;

#[pyclass]
pub struct PyPauliExponential {
    pe: PauliExponential,
    pauli_strategy: PauliPolynomialSynthStrategy,
    tableau_strategy: CliffordTableauSynthStrategy,
}

#[pymethods]
impl PyPauliExponential {
    #[new]
    pub fn new(num_qubits: usize) -> Self {
        let pe = PauliExponential::new(VecDeque::from(vec![]), CliffordTableau::new(num_qubits));
        Self {
            pe,
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
}

pub fn synthesize<G>(pe: &mut PyPauliExponential, circuit: &mut G)
where
    G: CliffordGates + Gates,
{
    let mut synth = PauliExponentialSynthesizer::from_strategy(
        pe.pauli_strategy.clone(),
        pe.tableau_strategy.clone(),
    );
    let pe = std::mem::take(&mut pe.pe);
    synth.synthesize(pe, circuit)
}
