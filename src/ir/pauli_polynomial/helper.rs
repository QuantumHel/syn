use std::collections::VecDeque;

use bitvec::{bitvec, order::Lsb0};
use itertools::Itertools;

use crate::{
    data_structures::{
        CliffordTableau, MaskedPropagateClifford, PauliLetter, PauliPolynomial, PropagateClifford,
    },
    ir::{CliffordGates, Gates},
};
use bitvec::prelude::BitVec;

impl PropagateClifford for VecDeque<PauliPolynomial> {
    fn cx(&mut self, control: usize, target: usize) -> &mut Self {
        for pauli_polynomial in self.iter_mut() {
            pauli_polynomial.cx(control, target);
        }
        self
    }

    fn s(&mut self, target: usize) -> &mut Self {
        for pauli_polynomial in self.iter_mut() {
            pauli_polynomial.s(target);
        }
        self
    }

    fn v(&mut self, target: usize) -> &mut Self {
        for pauli_polynomial in self.iter_mut() {
            pauli_polynomial.v(target);
        }
        self
    }
}

impl MaskedPropagateClifford for VecDeque<PauliPolynomial> {
    fn masked_cx(&self, control: usize, target: usize, mask: &BitVec) -> &Self {
        self[0].masked_cx(control, target, mask);
        for pauli_polynomial in self.iter().skip(1) {
            pauli_polynomial.masked_cx(
                control,
                target,
                &bitvec![usize, Lsb0; 1; pauli_polynomial.length()],
            );
        }
        self
    }

    fn masked_s(&self, target: usize, mask: &BitVec) -> &Self {
        self[0].masked_s(target, mask);
        for pauli_polynomial in self.iter().skip(1) {
            pauli_polynomial.masked_s(target, &bitvec![usize, Lsb0; 1; pauli_polynomial.length()]);
        }
        self
    }

    fn masked_v(&self, target: usize, mask: &BitVec) -> &Self {
        self[0].masked_v(target, mask);
        for pauli_polynomial in self.iter().skip(1) {
            pauli_polynomial.masked_v(target, &bitvec![usize, Lsb0; 1; pauli_polynomial.length()]);
        }
        self
    }
}

pub(super) fn push_down_pauli_polynomial_update<G>(
    pauli_polynomials: &mut VecDeque<PauliPolynomial>,
    repr: &mut G,
    clifford_tableau: &mut CliffordTableau,
    pauli_polynomial: PauliPolynomial,
    num_gadgets: usize,
    mut mask: BitVec,
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
                    pauli_polynomials.h(i);
                    clifford_tableau.h(i);
                    repr.h(i);
                }
                PauliLetter::Y => {
                    affected_qubits.push(i);
                    pauli_polynomial.masked_s(i, &mask);
                    pauli_polynomials.s(i);
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
