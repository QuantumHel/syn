use std::collections::VecDeque;

use crate::{
    data_structures::{CliffordTableau, PauliPolynomial},
    ir::{CliffordGates, Gates, Synthesizer},
};
use bitvec::{bitvec, order::Lsb0};

use super::helper::push_down_pauli_polynomial_update;

#[derive(Default)]
pub struct NaivePauliPolynomialSynthesizer {
    clifford_tableau: CliffordTableau,
}

impl NaivePauliPolynomialSynthesizer {
    pub fn set_clifford_tableau(&mut self, clifford_tableau: CliffordTableau) -> &mut Self {
        self.clifford_tableau = clifford_tableau;
        self
    }
}

impl<G> Synthesizer<VecDeque<PauliPolynomial>, G, CliffordTableau>
    for NaivePauliPolynomialSynthesizer
where
    G: CliffordGates + Gates,
{
    fn synthesize(
        &mut self,
        mut pauli_polynomials: VecDeque<PauliPolynomial>,
        repr: &mut G,
    ) -> CliffordTableau {
        let mut clifford_tableau = std::mem::take(&mut self.clifford_tableau);
        while !pauli_polynomials.is_empty() {
            let pauli_polynomial = pauli_polynomials.pop_front().unwrap();
            let num_gadgets: usize = pauli_polynomial.length();
            let mask = bitvec![usize, Lsb0; 1; num_gadgets];
            push_down_pauli_polynomial_update(
                &mut pauli_polynomials,
                repr,
                &mut clifford_tableau,
                pauli_polynomial,
                num_gadgets,
                mask,
            );
        }

        clifford_tableau
    }
}
