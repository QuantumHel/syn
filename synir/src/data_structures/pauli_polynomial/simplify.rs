//! Rewrites that shrink a [`PauliPolynomial`] without changing the unitary it denotes.
//!
//! Two passes, both run by the importer before synthesis:
//!
//! * [`merge_repeats`] — gadgets on the same Pauli string collapse into one, adding their angles.
//! * [`split_off_cliffords`] — gadgets whose angle is Clifford are separated out, so the caller
//!   can fold them into a tableau instead of emitting rotations for them.

use crate::data_structures::{PauliLetter, PauliPolynomial, PauliString};
use itertools::{Either, Itertools};
use std::collections::HashMap;

/// Finds gadgets that share a Pauli string.
///
/// Returns `(key, indices)` pairs for each string appearing more than once, sorted by key. The key
/// is the string packed into a `usize`, two bits per qubit — which caps this at 32 qubits on a
/// 64-bit target before distinct strings start colliding.
pub(crate) fn check_repeats(pp: &PauliPolynomial) -> Vec<(usize, Vec<usize>)> {
    let size = pp.num_qubits();
    let length = pp.length();
    let mut repeats = HashMap::<usize, Vec<usize>>::new();
    for index in 0..length {
        let mut num = 0;
        for letter in 0..size {
            num += (pp.chain(letter).x(index) as usize) << 2 * letter;
            num += (pp.chain(letter).z(index) as usize) << 2 * letter + 1;
        }
        repeats
            .entry(num)
            .and_modify(|e: &mut Vec<usize>| e.push(index))
            .or_insert(vec![index]);
    }
    repeats
        .into_iter()
        .filter(|(_, v)| v.len() > 1)
        .sorted()
        .collect_vec()
}

/// Applies the grouping from [`check_repeats`], summing angles and deleting the duplicates.
///
/// Each group's angle lands on its lowest index; the rest are removed. Removal preserves the order
/// of the surviving gadgets, unlike
/// [`PauliPolynomial::remove_gadget`](crate::data_structures::PauliPolynomial::remove_gadget).
pub(crate) fn _merge_repeats(
    pp: &mut PauliPolynomial,
    merge_list: Vec<(usize, Vec<usize>)>,
) -> &PauliPolynomial {
    let mut pp_merge_list = Vec::<usize>::new();
    // merge all the angles first
    for (_, angle_merge_list) in merge_list {
        let merge_index = angle_merge_list[0];
        let mut angle = pp.angle(merge_index);
        for angle_index in angle_merge_list.iter().skip(1) {
            angle += pp.angle(*angle_index);
        }
        pp.angles[merge_index] = angle;
        pp_merge_list.extend_from_slice(&angle_merge_list[1..]);
    }
    // remove duplicate entries
    pp_merge_list.sort_by(|a, b| b.cmp(a));
    for remove_index in pp_merge_list {
        pp.angles.remove(remove_index);
        for chain_index in 0..pp.chains().len() {
            pp.chains[chain_index].x.remove(remove_index);
            pp.chains[chain_index].z.remove(remove_index);
        }
    }
    pp
}

/// Collapses every repeated Pauli string into a single gadget carrying the summed angle.
///
/// Sound because gadgets on the same string commute and compose by adding angles.
pub fn merge_repeats(pp: &mut PauliPolynomial) -> &PauliPolynomial {
    let repeats = check_repeats(pp);
    _merge_repeats(pp, repeats)
}

/// Partitions a polynomial by whether each gadget's angle is Clifford.
///
/// Returns `(non_clifford, clifford)`. The Clifford half can be absorbed into a
/// [`CliffordTableau`](crate::data_structures::CliffordTableau) and stop costing rotation gates;
/// see [`push_clifford_angles_inverse`](crate::data_structures::pauli_exponential::push_clifford_angles_inverse),
/// which is careful about the direction that folding has to go.
///
/// Relative order is preserved within each half.
pub fn split_off_cliffords(pp: PauliPolynomial) -> (PauliPolynomial, PauliPolynomial) {
    // Removes clifford gates from this PauliPolynomial and returns them as a separate PP
    let mut to_remove = vec![];
    for (i, angle) in pp.angles.iter().enumerate() {
        if angle.is_clifford() {
            to_remove.push(i);
        }
    }
    let (filtered_chains, removed_chains) = pp
        .chains
        .iter()
        .map(|ps| {
            let (left, right): (Vec<PauliLetter>, Vec<PauliLetter>) =
                (0..ps.len()).into_iter().partition_map(|i| {
                    if !to_remove.contains(&i) {
                        Either::Left(ps.pauli(i))
                    } else {
                        Either::Right(ps.pauli(i))
                    }
                });
            (
                PauliString::from_letters(&left),
                PauliString::from_letters(&right),
            )
        })
        .unzip();
    let (filtered_angles, removed_angles) = pp.angles.iter().enumerate().partition_map(|(j, p)| {
        if !to_remove.contains(&j) {
            Either::Left(p)
        } else {
            Either::Right(p)
        }
    });
    (
        PauliPolynomial {
            chains: filtered_chains,
            angles: filtered_angles,
        },
        PauliPolynomial {
            chains: removed_chains,
            angles: removed_angles,
        },
    )
}

#[cfg(test)]
mod tests {
    use std::borrow::BorrowMut;

    use crate::data_structures::Angle;
    use crate::data_structures::PauliString;

    use super::*;

    #[test]
    fn test_simple_check_repeats() {
        // Combined reading from back -> 01 = 1
        let pp = PauliPolynomial::from_hamiltonian(vec![
            ("I", Angle::from_angle(1.0)),
            ("X", Angle::from_angle(2.0)),
            ("Z", Angle::from_angle(3.0)),
            ("X", Angle::from_angle(4.0)),
            ("X", Angle::from_angle(5.0)),
        ]);
        let repeats = check_repeats(&pp);
        assert!(repeats.len() == 1);
        assert_eq!(repeats, vec![(1, vec![1, 3, 4])]);
    }

    #[test]
    fn test_check_repeats() {
        // XIZY appears twice
        // Z string -> 0011
        // X string -> 1001
        // Combined reading from back -> 11 10 00 01 = 225
        let pp = PauliPolynomial::from_hamiltonian(vec![
            ("XIZY", Angle::from_angle(1.0)),
            ("XIZY", Angle::from_angle(2.0)),
            ("YZZI", Angle::from_angle(3.0)),
        ]);
        let repeats = check_repeats(&pp);
        assert!(repeats.len() == 1);
        assert_eq!(repeats, vec![(225, vec![0, 1])]);
    }

    #[test]
    fn test_multiple_repeats() {
        // Combined reading from back -> 01 = 1
        let pp = PauliPolynomial::from_hamiltonian(vec![
            ("II", Angle::from_angle(1.0)),
            ("IX", Angle::from_angle(2.0)),
            ("ZZ", Angle::from_angle(3.0)),
            ("IX", Angle::from_angle(4.0)),
            ("ZZ", Angle::from_angle(5.0)),
        ]);
        let repeats = check_repeats(&pp);
        assert!(repeats.len() == 2);
        assert_eq!(repeats, vec![(4, vec![1, 3]), (10, vec![2, 4])]);
    }

    #[test]
    fn test_simple_merge_repeats() {
        // Combined reading from back -> 01 = 1
        let mut pp = PauliPolynomial::from_hamiltonian(vec![
            ("I", Angle::from_angle(1.0)),
            ("X", Angle::from_angle(2.0)),
            ("Z", Angle::from_angle(3.0)),
            ("X", Angle::from_angle(4.0)),
            ("X", Angle::from_angle(5.0)),
        ]);
        let repeats = check_repeats(&pp);

        _merge_repeats(pp.borrow_mut(), repeats);

        assert!(pp.chain(0).len() == 3);
        assert_eq!(pp.chain(0), &PauliString::from_text("IXZ"));
        assert_eq!(pp.angles, Angle::from_angles(&[1.0, 11.0, 3.0]));
    }

    #[test]
    fn test_merge_repeats() {
        let mut pp = PauliPolynomial::from_hamiltonian(vec![
            ("XIZY", Angle::from_angle(1.0)),
            ("XIZY", Angle::from_angle(2.0)),
            ("YZZI", Angle::from_angle(3.0)),
        ]);
        let repeats = check_repeats(&pp);
        let pp = _merge_repeats(pp.borrow_mut(), repeats);

        assert!(pp.chain(0).len() == 2);
        assert_eq!(pp.chain(0), &PauliString::from_text("XY"));
        assert_eq!(pp.chain(1), &PauliString::from_text("IZ"));
        assert_eq!(pp.chain(2), &PauliString::from_text("ZZ"));
        assert_eq!(pp.chain(3), &PauliString::from_text("YI"));
        assert_eq!(pp.angles, Angle::from_angles(&[3.0, 3.0]));
    }

    #[test]
    fn test_multiple_merge_repeats() {
        // Combined reading from back -> 01 = 1
        let mut original_pp: PauliPolynomial = PauliPolynomial::from_hamiltonian(vec![
            ("II", Angle::from_angle(1.0)),
            ("IX", Angle::from_angle(2.0)),
            ("ZZ", Angle::from_angle(3.0)),
            ("IX", Angle::from_angle(4.0)),
            ("ZZ", Angle::from_angle(5.0)),
        ]);
        let repeats: Vec<(usize, Vec<usize>)> = check_repeats(&original_pp);
        let pp = _merge_repeats(original_pp.borrow_mut(), repeats);

        assert!(pp.chain(0).len() == 3);
        assert_eq!(pp.chain(0), &PauliString::from_text("IIZ"));
        assert_eq!(pp.chain(1), &PauliString::from_text("IXZ"));
        assert_eq!(pp.angles, Angle::from_angles(&[1.0, 6.0, 8.0]));
    }
}
