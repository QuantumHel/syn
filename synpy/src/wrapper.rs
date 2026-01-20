pub(crate) mod qiskit;

extern crate pyo3;
extern crate pyo3_ffi;

use std::collections::VecDeque;

use pyo3::prelude::*;
use synir::{
    data_structures::{
        Angle, CliffordTableau, PauliExponential, PauliPolynomial, PropagateClifford,
    },
    ir::{
        clifford_tableau::CliffordTableauSynthStrategy,
        pauli_exponential::PauliExponentialSynthesizer,
        pauli_polynomial::PauliPolynomialSynthStrategy, CliffordGates, Gates, Synthesizer,
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

    pub fn add_h(&mut self, target: usize) {
        self.pe.h(target);
    }

    pub fn add_s(&mut self, target: usize) {
        self.pe.s(target);
    }

    pub fn add_s_dgr(&mut self, target: usize) {
        self.pe.s_dgr(target);
    }

    pub fn add_x(&mut self, target: usize) {
        self.pe.x(target);
    }

    pub fn add_y(&mut self, target: usize) {
        self.pe.y(target);
    }

    pub fn add_z(&mut self, target: usize) {
        self.pe.z(target);
    }

    pub fn add_cx(&mut self, control: usize, target: usize) {
        self.pe.cx(control, target);
    }

    pub fn add_rz(&mut self, target: usize, angle: f64) {
        let size = self.pe.size();
        let ppvec = self.pe.chains();
        
        let newpp = PauliPolynomial::from_hamiltonian(vec![(
            &to_pauli_component(size, &target, 'Z'),
            Angle::Arbitrary(angle),
        )]);
        ppvec.push_front(newpp);
    }
}

fn to_pauli_component(size: usize, target: &usize, pauli: char) -> String {
    let mut term = String::new();
    for i in 0..size {
        if i == *target {
            term.push(pauli);
        } else {
            term.push('I');
        }
    }
    term
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
    synth.synthesize(pe, circuit);
}
