use std::collections::VecDeque;

use crate::data_structures::{CliffordTableau, PauliPolynomial};

#[derive(Default)]
pub struct PauliExponential {
    pub(crate) pauli_polynomials: VecDeque<PauliPolynomial>,
    pub(crate) clifford_tableau: CliffordTableau,
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
