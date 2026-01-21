use crate::data_structures::{PauliLetter, PauliPolynomial, PauliString};
use itertools::{Either, Itertools};
use std::collections::HashMap;

/// Check for items in PauliPolynomial that has the same Pauli string
/// Returns a vector of a tuple of Pauli string identifier and the locations of this repeated string
pub(crate) fn check_repeats(pp: &PauliPolynomial) -> Vec<(usize, Vec<usize>)> {
    let size = pp.size();
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

/// Takes output of check_repeats and merges items with the same Pauli string
/// Assumes that all angles are of the same type
pub(crate) fn _merge_repeats(
    mut pp: PauliPolynomial,
    merge_list: Vec<(usize, Vec<usize>)>,
) -> PauliPolynomial {
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

/// Merges all items in PauliPolynomial that share the same Pauli string
pub fn merge_repeats(mut pp: PauliPolynomial) -> PauliPolynomial {
    let repeats = check_repeats(&pp);
    _merge_repeats(pp, repeats)
}

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
                    if to_remove.contains(&i) {
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
        if to_remove.contains(&j) {
            Either::Left(p)
        } else {
            Either::Right(p)
        }
    });
    (
        PauliPolynomial {
            chains: filtered_chains,
            angles: filtered_angles,
            size: pp.size,
        },
        PauliPolynomial {
            chains: removed_chains,
            angles: removed_angles,
            size: pp.size,
        },
    )
}

#[cfg(test)]
mod tests {
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
        let pp = PauliPolynomial::from_hamiltonian(vec![
            ("I", Angle::from_angle(1.0)),
            ("X", Angle::from_angle(2.0)),
            ("Z", Angle::from_angle(3.0)),
            ("X", Angle::from_angle(4.0)),
            ("X", Angle::from_angle(5.0)),
        ]);
        let repeats = check_repeats(&pp);

        let pp = _merge_repeats(pp, repeats);

        assert!(pp.chain(0).len() == 3);
        assert_eq!(pp.chain(0), &PauliString::from_text("IXZ"));
        assert_eq!(pp.angles, Angle::from_angles(&[1.0, 11.0, 3.0]));
    }

    #[test]
    fn test_merge_repeats() {
        let pp = PauliPolynomial::from_hamiltonian(vec![
            ("XIZY", Angle::from_angle(1.0)),
            ("XIZY", Angle::from_angle(2.0)),
            ("YZZI", Angle::from_angle(3.0)),
        ]);
        let repeats = check_repeats(&pp);
        let pp = _merge_repeats(pp, repeats);

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
        let pp = PauliPolynomial::from_hamiltonian(vec![
            ("II", Angle::from_angle(1.0)),
            ("IX", Angle::from_angle(2.0)),
            ("ZZ", Angle::from_angle(3.0)),
            ("IX", Angle::from_angle(4.0)),
            ("ZZ", Angle::from_angle(5.0)),
        ]);

        let repeats = check_repeats(&pp);
        let pp = _merge_repeats(pp, repeats);

        assert!(pp.chain(0).len() == 3);
        assert_eq!(pp.chain(0), &PauliString::from_text("IIZ"));
        assert_eq!(pp.chain(1), &PauliString::from_text("IXZ"));
        assert_eq!(pp.angles, Angle::from_angles(&[1.0, 6.0, 8.0]));
    }
}
