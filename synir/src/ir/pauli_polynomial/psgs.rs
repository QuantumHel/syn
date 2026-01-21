use std::{collections::VecDeque, ops::Deref};

use crate::{
    architecture::connectivity::Connectivity,
    data_structures::{CliffordTableau, PauliExponential, PauliPolynomial},
    ir::{
        pauli_polynomial::{
            self,
            helper::{check_columns, identity_recurse},
        },
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
        let num_qubits = &pauli_polynomials[0].size();
        let mut clifford_tableau = CliffordTableau::new(*num_qubits);
        let mut maybe_pp = pauli_polynomials.pop_front();
        while maybe_pp.as_ref().is_some() {
            let mut pauli_polynomial = maybe_pp.unwrap();
            let num_gadgets: usize = pauli_polynomial.length();
            let mut polynomial_mask = bitvec![usize, Lsb0; 1; num_gadgets];
            let mut remainder_pe = PauliExponential {
                pauli_polynomials: pauli_polynomials.to_owned(),
                clifford_tableau: clifford_tableau.to_owned(),
            };
            check_columns(repr, &mut pauli_polynomial, &mut polynomial_mask);
            identity_recurse(
                &mut pauli_polynomial,
                &mut remainder_pe,
                &self.connectivity,
                polynomial_mask,
                repr,
            );
            pauli_polynomials = remainder_pe.pauli_polynomials;
            clifford_tableau = remainder_pe.clifford_tableau;
            maybe_pp = pauli_polynomials.pop_front();
        }
        clifford_tableau
    }
}
