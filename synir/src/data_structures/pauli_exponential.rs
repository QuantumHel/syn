use std::borrow::{Borrow, BorrowMut};
use std::collections::VecDeque;
use std::iter::zip;

use crate::data_structures::pauli_polynomial::simplify::{
    merge_repeats as merge_repeats_pp, split_off_cliffords,
};
use crate::data_structures::{CliffordTableau, PauliPolynomial, PropagateClifford};

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

    pub fn chains(&self) -> &VecDeque<PauliPolynomial> {
        &self.pauli_polynomials
    }

    pub fn mut_chains(&mut self) -> &mut VecDeque<PauliPolynomial> {
        &mut self.pauli_polynomials
    }

    pub fn size(&self) -> usize {
        self.clifford_tableau.size()
    }

    pub fn clifford_tableau(&self) -> &CliffordTableau {
        &self.clifford_tableau
    }

    pub fn mut_clifford_tableau(&mut self) -> &mut CliffordTableau {
        &mut self.clifford_tableau
    }
}

impl PropagateClifford for PauliExponential {
    fn cx(&mut self, control: crate::IndexType, target: crate::IndexType) -> &mut Self {
        self.pauli_polynomials.cx(control, target);
        self.clifford_tableau.cx(control, target);
        self
    }

    fn s(&mut self, target: crate::IndexType) -> &mut Self {
        self.pauli_polynomials.s(target);
        self.clifford_tableau.s(target);
        self
    }

    fn v(&mut self, target: crate::IndexType) -> &mut Self {
        self.pauli_polynomials.v(target);
        self.clifford_tableau.v(target);
        self
    }
}

pub fn merge_repeats(pe: &mut PauliExponential) {
    for pp in pe.mut_chains() {
        merge_repeats_pp(pp);
    }
}

pub fn merge_commuting_pp(pe: &mut PauliExponential) {
    let mut new_pps = vec![];
    pe.pauli_polynomials.make_contiguous(); // Stores all elements in as_mut_slices.0
    let (mut pp, mut old_pps) = match pe.pauli_polynomials.as_mut_slices().0.split_first_mut() {
        Some(split) => split,
        None => return,
    };
    while !old_pps.is_empty() {
        // Get first
        let (next_pp, tmp_pps) = old_pps.split_first_mut().unwrap();
        // Partition commuting
        let (commuting_idx, other_idx): (Vec<_>, Vec<_>) =
            (0..next_pp.length()).partition(|i| pp.commutes_with_gadget(next_pp.pauli_gadget(*i)));
        // Add commuting to pp
        if !commuting_idx.is_empty() {
            for i in commuting_idx {
                let (gadget, angle) = next_pp.remove_gadget(i);
                pp.append_gadget(gadget, angle);
            }
        }
        // Start new pp if next_pp did not commute
        if !other_idx.is_empty() {
            new_pps.push(pp.to_owned());
            pp = next_pp;
        }
        // Update old_pps to no longer contain the first.
        old_pps = tmp_pps;
    }
    // Add the last pp
    new_pps.push(pp.to_owned());
    // Update pauli_polynomials
    pe.pauli_polynomials = VecDeque::from(new_pps);
}

pub fn push_clifford_angles(pe: &mut PauliExponential) {
    let mut new_pps = VecDeque::new();
    pe.pauli_polynomials.make_contiguous();
    let mut maybe_fst = pe.pauli_polynomials.as_mut_slices().0.split_first_mut();
    while maybe_fst.is_some() {
        let (pp, others) = maybe_fst.unwrap();
        let (non_clifford, clifford) = split_off_cliffords(pp.to_owned());
        if non_clifford.length() > 0 {
            new_pps.push_back(non_clifford);
        }
        if clifford.length() > 0 {
            println!("Found {} clifford gadgets", clifford.length());
            for other_pp in others.iter_mut() {
                other_pp.push_clifford(&clifford);
            }
            for (gadget, angle) in zip(clifford.pauli_gadgets(), clifford.angles()) {
                pe.clifford_tableau
                    .compose_gadget((gadget, *angle))
                    .unwrap();
            }
        }
        maybe_fst = others.split_first_mut();
    }
    pe.pauli_polynomials = new_pps;
}
