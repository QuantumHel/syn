use std::collections::VecDeque;

use crate::data_structures::{CliffordTableau, PauliPolynomial};

pub mod naive;

#[derive(Default)]
pub struct PauliExponential {
    pauli_polynomials: VecDeque<PauliPolynomial>,
    clifford_tableau: CliffordTableau,
}

impl PauliExponential {
    pub fn new(
        pauli_polynomials: VecDeque<PauliPolynomial>,
        clifford_tableau: CliffordTableau,
    ) -> Self {
        PauliExponential {
            pauli_polynomials,
            clifford_tableau,
        }
    }
}

pub trait PauliExponentialSynthesizer<G> {
    fn synthesize(&mut self, external_repr: &mut G);
}
