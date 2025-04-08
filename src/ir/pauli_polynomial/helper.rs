use std::collections::VecDeque;

use itertools::Itertools;

use crate::{
    data_structures::{
        CliffordTableau, MaskedPropagateClifford, PauliLetter, PauliPolynomial, PropagateClifford,
    },
    ir::{CliffordGates, Gates},
};

pub(super) fn push_down_pauli_polynomial_update<G>(
    pauli_polynomials: &VecDeque<PauliPolynomial>,
    repr: &mut G,
    clifford_tableau: &mut CliffordTableau,
    pauli_polynomial: PauliPolynomial,
    num_gadgets: usize,
    mut mask: bitvec::prelude::BitVec,
) where
    G: CliffordGates + Gates,
{
    for col in 0..num_gadgets {
        let mut affected_qubits = Vec::new();
        for (i, row) in pauli_polynomial.chains().iter().enumerate() {
            match row.pauli(col) {
                PauliLetter::I => {}
                PauliLetter::X => {
                    affected_qubits.push(i);
                    pauli_polynomial.masked_h(i, &mask);
                    pauli_polynomials.masked_h(i, &mask);
                    clifford_tableau.h(i);
                    repr.h(i);
                }
                PauliLetter::Y => {
                    affected_qubits.push(i);
                    pauli_polynomial.masked_s(i, &mask);
                    pauli_polynomials.masked_s(i, &mask);
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
}
