use crate::data_structures::PauliPolynomial;
use itertools::Itertools;
use std::collections::HashMap;

pub fn check_repeats(pp: &PauliPolynomial) -> Vec<(usize, Vec<usize>)> {
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

pub fn merge_repeats(
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

#[cfg(test)]
mod tests {
    use crate::data_structures::PauliString;

    use super::*;

    #[test]
    fn test_simple_check_repeats() {
        // Combined reading from back -> 01 = 1
        let pp = PauliPolynomial::from_hamiltonian(vec![
            ("I", 1.0),
            ("X", 2.0),
            ("Z", 3.0),
            ("X", 4.0),
            ("X", 5.0),
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
        let pp =
            PauliPolynomial::from_hamiltonian(vec![("XIZY", 1.0), ("XIZY", 2.0), ("YZZI", 3.0)]);
        let repeats = check_repeats(&pp);
        assert!(repeats.len() == 1);
        assert_eq!(repeats, vec![(225, vec![0, 1])]);
    }

    #[test]
    fn test_multiple_repeats() {
        // Combined reading from back -> 01 = 1
        let pp = PauliPolynomial::from_hamiltonian(vec![
            ("II", 1.0),
            ("IX", 2.0),
            ("ZZ", 3.0),
            ("IX", 4.0),
            ("ZZ", 5.0),
        ]);
        let repeats = check_repeats(&pp);
        assert!(repeats.len() == 2);
        assert_eq!(repeats, vec![(4, vec![1, 3]), (10, vec![2, 4])]);
    }

    #[test]
    fn test_simple_merge_repeats() {
        // Combined reading from back -> 01 = 1
        let pp = PauliPolynomial::from_hamiltonian(vec![
            ("I", 1.0),
            ("X", 2.0),
            ("Z", 3.0),
            ("X", 4.0),
            ("X", 5.0),
        ]);
        let repeats = check_repeats(&pp);

        let pp = merge_repeats(pp, repeats);

        assert!(pp.chain(0).len() == 3);
        assert_eq!(pp.chain(0), &PauliString::from_text("IXZ"));
        assert_eq!(pp.angles, &[1.0, 11.0, 3.0]);
    }

    #[test]
    fn test_merge_repeats() {
        let pp =
            PauliPolynomial::from_hamiltonian(vec![("XIZY", 1.0), ("XIZY", 2.0), ("YZZI", 3.0)]);
        let repeats = check_repeats(&pp);
        let pp = merge_repeats(pp, repeats);

        assert!(pp.chain(0).len() == 2);
        assert_eq!(pp.chain(0), &PauliString::from_text("XY"));
        assert_eq!(pp.chain(1), &PauliString::from_text("IZ"));
        assert_eq!(pp.chain(2), &PauliString::from_text("ZZ"));
        assert_eq!(pp.chain(3), &PauliString::from_text("YI"));
        assert_eq!(pp.angles, &[3.0, 3.0]);
    }

    #[test]
    fn test_multiple_merge_repeats() {
        // Combined reading from back -> 01 = 1
        let pp = PauliPolynomial::from_hamiltonian(vec![
            ("II", 1.0),
            ("IX", 2.0),
            ("ZZ", 3.0),
            ("IX", 4.0),
            ("ZZ", 5.0),
        ]);

        let repeats = check_repeats(&pp);
        let pp = merge_repeats(pp, repeats);

        assert!(pp.chain(0).len() == 3);
        assert_eq!(pp.chain(0), &PauliString::from_text("IIZ"));
        assert_eq!(pp.chain(1), &PauliString::from_text("IXZ"));
        assert_eq!(pp.angles, &[1.0, 6.0, 8.0]);
    }
}
