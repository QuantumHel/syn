use std::{cell::Cell, collections::VecDeque, ops::Deref};

use crate::{
    architecture::connectivity::Connectivity,
    data_structures::{CliffordTableau, PauliExponential, PauliPolynomial},
    ir::{
        pauli_polynomial::{self, helper::identity_recurse},
        CliffordGates, Gates, Synthesizer,
    },
};
use bitvec::{bitvec, order::Lsb0};

/// Pauli String Gadget Synthesis: recursive partitioning that shares CNOT structure.
///
/// Instead of handling one gadget at a time, PSGS looks for a qubit whose column splits the
/// remaining gadgets into two sizeable groups, diagonalizes on that qubit, and recurses into each
/// group with a mask restricting it to those gadgets. Gadgets that share legs therefore share the
/// CNOT ladder that folds those legs together, which is where the saving over
/// [`NaivePauliPolynomialSynthesizer`](super::NaivePauliPolynomialSynthesizer) comes from.
///
/// The recursion consumes gadgets by *blanking* them
/// ([`PauliPolynomial::clear_gadget`](crate::data_structures::PauliPolynomial::clear_gadget))
/// rather than removing them, so the positional masks held across recursive calls stay valid.
///
/// Set the target hardware graph with [`set_connectivity`](Self::set_connectivity); the default is
/// an empty graph, so callers that do not set one get a single-qubit architecture.
///
/// # Panics
///
/// Panics if handed an empty deque, or if the recursion finishes with a gadget left
/// unsynthesized — a bug rather than a usage error, but one worth failing loudly on, since
/// otherwise it surfaces as a silently wrong unitary.
#[derive(Default)]
pub struct PSGSPauliPolynomialSynthesizer {
    connectivity: Connectivity,
}

impl PSGSPauliPolynomialSynthesizer {
    /// Sets the hardware graph the recursion routes against.
    pub fn set_connectivity(&mut self, connectivity: Connectivity) -> &mut Self {
        self.connectivity = connectivity;
        self
    }

    /// The hardware graph this synthesizer targets.
    pub fn get_connectivity(&self) -> &Connectivity {
        &self.connectivity
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
        let num_qubits = &pauli_polynomials[0].num_qubits();
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
            identity_recurse(
                &mut pauli_polynomial,
                &mut remainder_pe,
                &self.connectivity,
                polynomial_mask,
                repr,
            );
            // `identity_recurse` blanks each gadget as it emits it, so anything left non-identity
            // is a gadget that fell out of the recursion's masks and was never synthesized. That
            // used to surface only as a wrong unitary, many gates downstream.
            let unsynthesized = (0..pauli_polynomial.length())
                .filter(|i| pauli_polynomial.pauli_gadget(*i).weight() > 0)
                .map(|i| format!("{} @ {}", pauli_polynomial.pauli_gadget(i), pauli_polynomial.angle(i)))
                .collect::<Vec<_>>();
            assert!(
                unsynthesized.is_empty(),
                "PSGS left {} gadget(s) unsynthesized: [{}]",
                unsynthesized.len(),
                unsynthesized.join(", ")
            );
            pauli_polynomials = remainder_pe.pauli_polynomials;
            clifford_tableau = remainder_pe.clifford_tableau;
            maybe_pp = pauli_polynomials.pop_front();
        }
        clifford_tableau
    }
}
