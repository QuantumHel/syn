use std::{collections::VecDeque, iter::zip};

use bitvec::{bitvec, order::Lsb0};
use itertools::Itertools;

use crate::{
    architecture::{connectivity::Connectivity, Architecture},
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
    fn masked_cx(&mut self, control: usize, target: usize, mask: &BitVec) -> &mut Self {
        self[0].masked_cx(control, target, mask);
        for pauli_polynomial in self.iter_mut().skip(1) {
            pauli_polynomial.masked_cx(
                control,
                target,
                &bitvec![usize, Lsb0; 1; pauli_polynomial.length()],
            );
        }
        self
    }

    fn masked_s(&mut self, target: usize, mask: &BitVec) -> &mut Self {
        self[0].masked_s(target, mask);
        for pauli_polynomial in self.iter_mut().skip(1) {
            pauli_polynomial.masked_s(target, &bitvec![usize, Lsb0; 1; pauli_polynomial.length()]);
        }
        self
    }

    fn masked_v(&mut self, target: usize, mask: &BitVec) -> &mut Self {
        self[0].masked_v(target, mask);
        for pauli_polynomial in self.iter_mut().skip(1) {
            pauli_polynomial.masked_v(target, &bitvec![usize, Lsb0; 1; pauli_polynomial.length()]);
        }
        self
    }
}

pub(super) fn push_down_pauli_polynomial_update<G>(
    pauli_polynomials: &mut VecDeque<PauliPolynomial>,
    repr: &mut G,
    clifford_tableau: &mut CliffordTableau,
    mut pauli_polynomial: PauliPolynomial,
    num_gadgets: usize,
    mut mask: BitVec,
) where
    G: CliffordGates + Gates,
{
    for col in 0..num_gadgets {
        let mut affected_qubits = Vec::new();
        for i in 0..pauli_polynomial.size() {
            let row = pauli_polynomial.chain(i);
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
                    pauli_polynomial.masked_v(i, &mask);
                    pauli_polynomials.v(i);
                    clifford_tableau.v(i);
                    repr.v(i);
                }
                PauliLetter::Z => {
                    affected_qubits.push(i);
                }
            }
        }
        if affected_qubits.len() > 1 {
            for (&control, &target) in affected_qubits.iter().tuple_windows() {
                pauli_polynomial.masked_cx(control, target, &mask);
                pauli_polynomials.cx(control, target);
                clifford_tableau.cx(control, target);
                repr.cx(control, target);
            }
        }
        let last_qubit = *affected_qubits.last().unwrap();
        repr.rz(last_qubit, pauli_polynomial.angle(col));
        mask.replace(col, false);
    }
}

pub(super) fn pick_qubit(
    pauli_polynomial: &PauliPolynomial,
    polynomial_mask: &BitVec,
    selected_qubits: &[usize],
) -> usize {
    let weight_i = 10;
    let mut costs = vec![0usize; selected_qubits.len()];
    for (index, qubit) in selected_qubits.iter().enumerate() {
        let chain = pauli_polynomial.chain(*qubit);
        for (bit, mask) in zip(chain.combine(), polynomial_mask.iter()) {
            if !mask {
                continue;
            }
            let weight = if bit { 1 } else { weight_i };
            costs[index] += weight;
        }
    }
    // If all costs are zero, return the first qubit
    // Find the qubit with the maximum cost
    selected_qubits[costs.iter().position_max().unwrap()]
}

/// Checks if any Pauli Gadgets have no or 1 identity legs.
/// The gadget is removed from scope and the bit in polynomial mask is set to 0
pub(super) fn check_columns<G>(
    repr: &mut G,
    pauli_polynomial: &mut PauliPolynomial,
    polynomial_mask: &mut BitVec,
) where
    G: CliffordGates + Gates,
{
    let invalid = {
        let length = pauli_polynomial.length();

        let mut invalid = BitVec::repeat(false, length);
        let mut seen_one = BitVec::repeat(false, length);

        for chain in pauli_polynomial.chains().iter() {
            let combination = chain.combine();
            invalid |= seen_one.clone() & &combination;
            seen_one |= combination;
        }
        invalid
    };
    let length = pauli_polynomial.length();
    let PauliPolynomial {
        ref mut chains,
        angles,
        ..
    } = pauli_polynomial;
    for index in (0..length).rev() {
        if !invalid[index] && polynomial_mask[index] {
            polynomial_mask.swap_remove(index);
            let angle = angles.swap_remove(index);
            for (qubit, chain) in chains.iter_mut().enumerate() {
                let letter = chain.swap_remove(index);
                match letter {
                    PauliLetter::X => {
                        repr.rx(qubit, angle);
                    }
                    PauliLetter::Y => {
                        repr.ry(qubit, angle);
                    }
                    PauliLetter::Z => {
                        repr.rz(qubit, angle);
                    }
                    PauliLetter::I => {}
                }
            }
        }
    }
}

/// Partition the input polynomial mask on the selected qubit. If the gadget has an identity leg, it is assigned to the identity mask. Otherwise, it is assigned to the other mask.
pub(super) fn identity_partition(
    pauli_polynomial: &PauliPolynomial,
    mut polynomial_mask: BitVec,
    selected_qubit: usize,
) -> (BitVec, BitVec) {
    let pauli_polynomial_length = pauli_polynomial.length();
    let mut identity_mask = BitVec::with_capacity(pauli_polynomial_length);
    let mut other_mask = BitVec::with_capacity(pauli_polynomial_length);

    let pauli_chain = pauli_polynomial.chain(selected_qubit).iter();

    for pauli in pauli_chain.iter() {
        match pauli {
            PauliLetter::I => {
                identity_mask.push(polynomial_mask.pop().unwrap());
                other_mask.push(false);
            }
            _ => {
                other_mask.push(polynomial_mask.pop().unwrap());
                identity_mask.push(false);
            }
        }
    }

    (identity_mask, other_mask)
}

/// Partition the input polynomial mask on the selected qubit. If the gadget has an identity leg, it is assigned to the identity mask. Otherwise, it is assigned to the other mask.
pub(super) fn max_partition(
    pauli_polynomial: &PauliPolynomial,
    mut polynomial_mask: BitVec,
    selected_qubit: usize,
) -> (
    BitVec,
    PauliLetter,
    BitVec,
    PauliLetter,
    BitVec,
    PauliLetter,
) {
    let pauli_polynomial_length = pauli_polynomial.length();

    let mut x_mask = BitVec::with_capacity(pauli_polynomial_length);
    let mut y_mask = BitVec::with_capacity(pauli_polynomial_length);
    let mut z_mask = BitVec::with_capacity(pauli_polynomial_length);

    let pauli_chain = pauli_polynomial.chain(selected_qubit).iter();

    for pauli in pauli_chain.iter() {
        match pauli {
            PauliLetter::I => {
                panic!("Cannot partition polynomial with identity leg on selected qubit");
            }
            PauliLetter::X => {
                x_mask.push(polynomial_mask.pop().unwrap());
                y_mask.push(false);
                z_mask.push(false);
            }
            PauliLetter::Y => {
                x_mask.push(false);
                y_mask.push(polynomial_mask.pop().unwrap());
                z_mask.push(false);
            }
            PauliLetter::Z => {
                x_mask.push(false);
                y_mask.push(false);
                z_mask.push(polynomial_mask.pop().unwrap());
            }
        }
    }
    let mut polynomial_parts = vec![
        (x_mask, PauliLetter::X),
        (y_mask, PauliLetter::Y),
        (z_mask, PauliLetter::Z),
    ];
    polynomial_parts.sort_by(|a, b| a.0.count_ones().cmp(&b.0.count_ones()));

    let (largest_mask, largest_pauli) = polynomial_parts.pop().unwrap();
    let (second_mask, second_pauli) = polynomial_parts.pop().unwrap();
    let (third_mask, third_pauli) = polynomial_parts.pop().unwrap();

    (
        largest_mask,
        largest_pauli,
        second_mask,
        second_pauli,
        third_mask,
        third_pauli,
    )
}

pub(super) fn identity_recurse<G>(
    pauli_polynomial: &mut PauliPolynomial,
    clifford_tableau: &mut CliffordTableau,
    connectivity: &Connectivity,
    mut polynomial_mask: BitVec,
    // selected_qubits: &[usize],
    repr: &mut G,
) where
    G: CliffordGates + Gates,
{
    // Remove gadgets with single rotation.
    check_columns(repr, pauli_polynomial, &mut polynomial_mask);
    if polynomial_mask.count_ones() == 0 {
        return;
    }
    let selected_qubits = connectivity.non_cutting();
    let selected_qubit = pick_qubit(pauli_polynomial, &polynomial_mask, selected_qubits);
    // Create new connectivity without the selected qubit
    let reduced_connectivity = connectivity.disconnect(selected_qubit);

    let (identity_mask, other_mask) =
        identity_partition(pauli_polynomial, polynomial_mask, selected_qubit);

    if identity_mask.count_ones() > 0 {
        // recurse down identity mask
        identity_recurse(
            pauli_polynomial,
            clifford_tableau,
            &reduced_connectivity,
            identity_mask,
            repr,
        );
        // ensure remainder is synthesized
        identity_recurse(
            pauli_polynomial,
            clifford_tableau,
            connectivity,
            other_mask,
            repr,
        )
    } else {
        // `identity_mask` is empty, we do not process it
        let (largest_mask, largest_pauli, mut remaining_mask, _, third_mask, _) =
            max_partition(pauli_polynomial, other_mask, selected_qubit);

        remaining_mask |= third_mask.as_bitslice();

        // Ensure `next_qubit` is always a neighbor of `selected_qubit`
        let next_qubit = pick_qubit(
            pauli_polynomial,
            &largest_mask,
            &connectivity.neighbors(selected_qubit),
        );

        // Check if there are identities on `next_qubit`
        let (next_identity_mask, next_other_mask) =
            identity_partition(pauli_polynomial, largest_mask, next_qubit);

        // Ensure that selected qubit is always Pauli::Z
        diagonalize_qubit(
            pauli_polynomial,
            clifford_tableau,
            repr,
            selected_qubit,
            largest_pauli,
        );

        if next_identity_mask.count_ones() > 0 {
            disconnect_i(
                pauli_polynomial,
                clifford_tableau,
                repr,
                selected_qubit,
                next_qubit,
            );
            remaining_mask |= next_other_mask.as_bitslice();

            let (identity_mask, other_mask) =
                identity_partition(pauli_polynomial, remaining_mask, selected_qubit);

            identity_recurse(
                pauli_polynomial,
                clifford_tableau,
                &reduced_connectivity,
                identity_mask,
                repr,
            );

            identity_recurse(
                pauli_polynomial,
                clifford_tableau,
                connectivity,
                other_mask,
                repr,
            )
        } else {
            let (
                mut largest_mask,
                largest_next_pauli,
                second_mask,
                second_next_pauli,
                third_mask,
                _,
            ) = max_partition(pauli_polynomial, next_other_mask, next_qubit);

            largest_mask |= second_mask.as_bitslice();
            let is_x = largest_next_pauli == PauliLetter::X || second_next_pauli == PauliLetter::X;
            let is_y = largest_next_pauli == PauliLetter::Y || second_next_pauli == PauliLetter::Y;

            disconnect(
                pauli_polynomial,
                clifford_tableau,
                repr,
                selected_qubit,
                next_qubit,
                is_x,
                is_y,
            );

            identity_recurse(
                pauli_polynomial,
                clifford_tableau,
                &reduced_connectivity,
                largest_mask,
                repr,
            );
            identity_recurse(
                pauli_polynomial,
                clifford_tableau,
                connectivity,
                third_mask,
                repr,
            );
        }
    }
}

fn disconnect_i<G>(
    pauli_polynomial: &mut PauliPolynomial,
    clifford_tableau: &mut CliffordTableau,
    repr: &mut G,
    selected_qubit: usize,
    next_qubit: usize,
) where
    G: CliffordGates + Gates,
{
    pauli_polynomial.cx(selected_qubit, next_qubit);
    pauli_polynomial.cx(next_qubit, selected_qubit);
    clifford_tableau.cx(selected_qubit, next_qubit);
    clifford_tableau.cx(next_qubit, selected_qubit);
    repr.cx(selected_qubit, next_qubit);
    repr.cx(next_qubit, selected_qubit);
}

fn disconnect<G>(
    pauli_polynomial: &mut PauliPolynomial,
    clifford_tableau: &mut CliffordTableau,
    repr: &mut G,
    selected_qubit: usize,
    next_qubit: usize,
    is_x: bool,
    is_y: bool,
) where
    G: CliffordGates + Gates,
{
    match (is_x, is_y) {
        (true, true) => {
            pauli_polynomial.h(next_qubit);
            clifford_tableau.h(next_qubit);
            repr.h(next_qubit);
        }
        (true, false) => {
            pauli_polynomial.s(next_qubit);
            clifford_tableau.s(next_qubit);
            repr.s(next_qubit);
        }
        _ => {}
    }
    pauli_polynomial.cx(selected_qubit, next_qubit);
    clifford_tableau.cx(selected_qubit, next_qubit);
    repr.cx(selected_qubit, next_qubit);
}

fn diagonalize_qubit<G>(
    pauli_polynomial: &mut PauliPolynomial,
    clifford_tableau: &mut CliffordTableau,
    repr: &mut G,
    selected_qubit: usize,
    largest_pauli: PauliLetter,
) where
    G: CliffordGates + Gates,
{
    match largest_pauli {
        PauliLetter::I => panic!("Should not have Pauli::I here"),
        PauliLetter::X => {
            pauli_polynomial.h(selected_qubit);
            clifford_tableau.h(selected_qubit);
            repr.h(selected_qubit);
        }
        PauliLetter::Y => {
            pauli_polynomial.v(selected_qubit);
            clifford_tableau.v(selected_qubit);
            repr.v(selected_qubit);
        }
        PauliLetter::Z => {}
    }
}
