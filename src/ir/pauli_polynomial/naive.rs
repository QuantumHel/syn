use crate::{
    data_structures::{
        CliffordTableau, MaskedPropagateClifford, PauliLetter, PauliPolynomial, PropagateClifford,
    },
    ir::{CliffordGates, Gates},
};
use bitvec::{bitvec, order::Lsb0};
use itertools::Itertools;

use super::PauliPolynomialSynthesizer;
pub struct NaivePauliPolynomialSynthesizer {
    pauli_polynomial: PauliPolynomial,
    clifford_tableau: CliffordTableau,
}

impl NaivePauliPolynomialSynthesizer {
    pub fn new(pauli_polynomial: PauliPolynomial, clifford_tableau: CliffordTableau) -> Self {
        assert!(pauli_polynomial.size() == clifford_tableau.size());
        Self {
            pauli_polynomial,
            clifford_tableau,
        }
    }
}

impl<G> PauliPolynomialSynthesizer<G> for NaivePauliPolynomialSynthesizer
where
    G: CliffordGates + Gates,
{
    fn synthesize(&mut self, repr: &mut G) -> CliffordTableau {
        let pauli_polynomial = &mut self.pauli_polynomial;
        let mut clifford_tableau = std::mem::take(&mut self.clifford_tableau);

        let num_gadgets = pauli_polynomial.length();
        let gadget_length = pauli_polynomial.length();
        let mut mask = bitvec![usize, Lsb0; 1; gadget_length];
        for col in 0..num_gadgets {
            let mut affected_qubits = Vec::new();
            for (i, row) in pauli_polynomial.chains().iter().enumerate() {
                match row.pauli(col) {
                    PauliLetter::I => {}
                    PauliLetter::X => {
                        affected_qubits.push(i);
                        pauli_polynomial.masked_h(i, &mask);
                        clifford_tableau.h(i);
                        repr.h(i);
                    }
                    PauliLetter::Y => {
                        affected_qubits.push(i);
                        pauli_polynomial.masked_s(i, &mask);
                        clifford_tableau.s(i);
                        repr.s(i);
                    }
                    PauliLetter::Z => {
                        affected_qubits.push(i);
                    }
                }
            }
            if affected_qubits.len() > 1 {
                for (&control, &target) in affected_qubits.iter().tuple_windows() {
                    pauli_polynomial.masked_cx(control, target, &mask);
                    clifford_tableau.cx(control, target);
                    repr.cx(control, target);
                }
            }
            let last_qubit = *affected_qubits.last().unwrap();
            repr.rz(last_qubit, pauli_polynomial.angle(col));
            mask.replace(col, false);
        }
        clifford_tableau
    }
}
