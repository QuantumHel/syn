use std::collections::VecDeque;

use crate::{
    architecture::connectivity::Connectivity,
    data_structures::{CliffordTableau, PauliPolynomial},
    ir::{
        pauli_polynomial::helper::{check_columns, identity_recurse},
        CliffordGates, Gates, Synthesizer,
    },
};
use bitvec::{bitvec, order::Lsb0};

#[derive(Default)]
pub struct PSGSPauliPolynomialSynthesizer {
    connectivity: Connectivity,
}

impl PSGSPauliPolynomialSynthesizer {
    pub fn set_connectivity(&mut self, connectivity: Connectivity) -> &mut Self {
        self.connectivity = connectivity;
        self
    }
}

impl<G> Synthesizer<VecDeque<PauliPolynomial>, G, CliffordTableau>
    for PSGSPauliPolynomialSynthesizer
where
    G: CliffordGates + Gates,
{
    fn synthesize(
        &mut self,
        mut pauli_polynomials: VecDeque<PauliPolynomial>,
        repr: &mut G,
    ) -> CliffordTableau {
        if pauli_polynomials.is_empty() {
            panic!("You are trying to synthesize an empty PauliPolynomial.")
        }
        let mut clifford_tableau = CliffordTableau::new(pauli_polynomials[0].size());
        while !pauli_polynomials.is_empty() {
            let mut pauli_polynomial = pauli_polynomials.pop_front().unwrap();
            let num_gadgets: usize = pauli_polynomial.length();
            let mut polynomial_mask = bitvec![usize, Lsb0; 1; num_gadgets];
            check_columns(repr, &mut pauli_polynomial, &mut polynomial_mask);
            identity_recurse(
                &mut pauli_polynomial,
                &mut clifford_tableau,
                &self.connectivity,
                polynomial_mask,
                repr,
            );
        }
        clifford_tableau
    }
}
