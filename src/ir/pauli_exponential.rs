use crate::data_structures::{CliffordTableau, PauliPolynomial};

pub mod naive;

#[derive(Default)]
pub struct PauliExponential {
    pauli_polynomial: PauliPolynomial,
    clifford_tableau: CliffordTableau,
}

impl PauliExponential {
    pub fn new(pauli_polynomial: PauliPolynomial, clifford_tableau: CliffordTableau) -> Self {
        PauliExponential {
            pauli_polynomial,
            clifford_tableau,
        }
    }
}

pub trait PauliExponentialSynthesizer<G> {
    fn synthesize(&mut self, external_repr: &mut G);
}
