//! The top-level normal form, plus the preprocessing passes that tidy it.
//!
//! See [`PauliExponential`].

use std::borrow::{Borrow, BorrowMut};
use std::collections::VecDeque;
use std::iter::zip;

use crate::data_structures::pauli_polynomial::simplify::{
    merge_repeats as merge_repeats_pp, split_off_cliffords,
};
use crate::data_structures::{CliffordTableau, PauliPolynomial, PropagateClifford};

/// A whole circuit in synIR's normal form: gadgets first, then one Clifford.
///
/// ```text
/// U = T * U_pp
/// ```
///
/// in matrix order, so `U_pp` runs first. `U_pp` is the deque of
/// [`PauliPolynomial`]s in circuit order, and `T` is the trailing [`CliffordTableau`].
///
/// Why a deque of polynomials rather than one flat list of gadgets: the gadgets *within* a single
/// polynomial pairwise commute, so a synthesizer may reorder them freely and pay for their shared
/// structure once. A gadget that fails to commute starts a new polynomial. That grouping is what
/// [`merge_commuting_pp`] computes.
///
/// Build one with [`crate::import`] rather than by hand; consume it with
/// [`PauliExponentialSynthesizer`](crate::ir::pauli_exponential::PauliExponentialSynthesizer).
#[derive(Default)]
pub struct PauliExponential {
    pub(crate) pauli_polynomials: VecDeque<PauliPolynomial>,
    pub(crate) clifford_tableau: CliffordTableau,
}

impl PauliExponential {
    /// Pairs a run of polynomials with the Clifford that follows them.
    pub fn new(
        pauli_polynomials: VecDeque<PauliPolynomial>,
        clifford_tableau: CliffordTableau,
    ) -> Self {
        PauliExponential {
            pauli_polynomials,
            clifford_tableau,
        }
    }

    /// The polynomials, in circuit order.
    pub fn chains(&self) -> &VecDeque<PauliPolynomial> {
        &self.pauli_polynomials
    }

    /// Mutable access to the polynomials.
    pub fn mut_chains(&mut self) -> &mut VecDeque<PauliPolynomial> {
        &mut self.pauli_polynomials
    }

    /// The number of qubits, taken from the tableau.
    pub fn size(&self) -> usize {
        self.clifford_tableau.size()
    }

    /// The trailing Clifford `T`.
    pub fn clifford_tableau(&self) -> &CliffordTableau {
        &self.clifford_tableau
    }

    /// Mutable access to the trailing Clifford.
    pub fn mut_clifford_tableau(&mut self) -> &mut CliffordTableau {
        &mut self.clifford_tableau
    }
}

/// Records a gate that a synthesizer has just **emitted**.
///
/// With the outstanding work written as `T^dag * U_rest`, emitting `G` means `T <- G * T` (append)
/// and `U_rest <- G U_rest G^dag`. `PropagateClifford for PauliPolynomial` conjugates the other
/// way (`Q -> G^dag Q G`), so the polynomial side takes the **adjoint** gate while the tableau
/// takes the gate itself. The two sides therefore disagree by a dagger on purpose; `h`, `x`, `y`,
/// `z` and `cx` are self-adjoint, which is why the mismatch stayed hidden.
///
/// This is *not* the operation an importer wants -- folding a circuit into the IR needs the
/// tableau prepended instead. See [`crate::import::CircuitImporter`].
impl PropagateClifford for PauliExponential {
    fn cx(&mut self, control: crate::IndexType, target: crate::IndexType) -> &mut Self {
        self.pauli_polynomials.cx(control, target);
        self.clifford_tableau.cx(control, target);
        self
    }

    fn s(&mut self, target: crate::IndexType) -> &mut Self {
        self.pauli_polynomials.s_dgr(target);
        self.clifford_tableau.s(target);
        self
    }

    fn v(&mut self, target: crate::IndexType) -> &mut Self {
        self.pauli_polynomials.v_dgr(target);
        self.clifford_tableau.v(target);
        self
    }
}

/// Merges gadgets that share a Pauli string, within each polynomial.
///
/// Two gadgets on the same string are one rotation by the sum of their angles, so this strictly
/// shrinks the polynomial. Runs the `simplify::merge_repeats` pass on each member.
pub fn merge_repeats(pe: &mut PauliExponential) {
    for pp in pe.mut_chains() {
        merge_repeats_pp(pp);
    }
}

/// Regroups the polynomials so each one holds only mutually commuting gadgets, greedily.
///
/// Walks the deque in order, pulling every gadget that commutes with the polynomial being built
/// into it and starting a fresh polynomial at the first one that does not. Larger commuting blocks
/// give the synthesizers more freedom to reorder, so this is worth running before synthesis.
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

/// Splits Clifford-angle gadgets out of the polynomials and folds them into the tableau.
///
/// Only valid while the tableau is held in **inverse form**, which is the case inside
/// `CircuitImporter::finish` -- its sole caller. Splitting `P_1 = D * P_1'` and pushing `D` left
/// past the later polynomials gives `T <- T * D`, so the inverse-form tableau must absorb `D^dag`,
/// i.e. the gadget with the angle negated.
pub fn push_clifford_angles_inverse(pe: &mut PauliExponential) {
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
            for other_pp in others.iter_mut() {
                other_pp.push_clifford(&clifford);
            }
            for (gadget, angle) in zip(clifford.pauli_gadgets(), clifford.angles()) {
                // `-angle` because the tableau is in inverse form here: see the doc comment.
                pe.clifford_tableau
                    .compose_gadget((gadget, -*angle))
                    .unwrap();
            }
        }
        maybe_fst = others.split_first_mut();
    }
    pe.pauli_polynomials = new_pps;
}
