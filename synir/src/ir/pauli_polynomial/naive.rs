use std::collections::VecDeque;

use crate::{
    data_structures::{CliffordTableau, PauliPolynomial},
    ir::{CliffordGates, Gates, Synthesizer},
};
use bitvec::{bitvec, order::Lsb0};

use super::helper::push_down_pauli_polynomial_update;

/// Synthesizes each polynomial by pushing its gadgets down one at a time.
///
/// No structure is shared between gadgets: each is diagonalized, emitted, and propagated through
/// whatever remains. Simple and connectivity-agnostic — the baseline that
/// [`PSGSPauliPolynomialSynthesizer`](super::psgs::PSGSPauliPolynomialSynthesizer) is measured
/// against.
///
/// Returns the tableau of basis changes left in place; see the module documentation for why they
/// are not un-computed here.
#[derive(Default)]
pub struct NaivePauliPolynomialSynthesizer {}

impl NaivePauliPolynomialSynthesizer {}

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
        let mut clifford_tableau = CliffordTableau::new(pauli_polynomials[0].num_qubits());
        while !pauli_polynomials.is_empty() {
            let pauli_polynomial = pauli_polynomials.pop_front().unwrap();
            let num_gadgets = pauli_polynomial.length();
            let mask = bitvec![usize, Lsb0; 1; num_gadgets];
            push_down_pauli_polynomial_update(
                &mut pauli_polynomials,
                repr,
                &mut clifford_tableau,
                pauli_polynomial,
                num_gadgets,
                mask,
                false,
            );
        }

        clifford_tableau
    }
}
