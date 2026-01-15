use std::iter::zip;

use super::{pauli_string::PauliString, IndexType, MaskedPropagateClifford, PropagateClifford};
use crate::data_structures::{Angle, PauliLetter};
use bitvec::vec::BitVec;
use itertools::zip_eq;

mod simplify;

#[derive(Debug, Clone, Default)]
pub struct PauliPolynomial {
    chains: Vec<PauliString>,
    angles: Vec<Angle>,
    size: usize,
}

impl PauliPolynomial {
    pub fn from_hamiltonian(hamiltonian_representation: Vec<(&str, Angle)>) -> Self {
        assert!(!hamiltonian_representation.is_empty());
        let terms = hamiltonian_representation.len();
        let num_qubits = hamiltonian_representation[0].0.len();
        let mut angles = Vec::<Angle>::with_capacity(terms);
        let mut chain_strings = vec![String::with_capacity(terms); num_qubits];
        //let chains = vec![PauliString::new(); num_qubits];
        for (pauli_string, angle) in hamiltonian_representation {
            zip_eq(chain_strings.iter_mut(), pauli_string.chars()).for_each(
                |(chain, pauli_letter)| {
                    chain.push(pauli_letter);
                },
            );
            angles.push(angle);
        }
        let chains = chain_strings
            .iter()
            .map(|gadget| PauliString::from_text(gadget))
            .collect::<Vec<_>>();

        PauliPolynomial {
            chains,
            angles,
            size: num_qubits,
        }
    }

    pub fn size(&self) -> usize {
        self.size
    }

    pub fn length(&self) -> usize {
        self.angles.len()
    }

    pub fn chain(&self, i: usize) -> &PauliString {
        &self.chains[i]
    }

    pub fn chains(&self) -> &Vec<PauliString> {
        &self.chains
    }

    pub fn angle(&self, i: usize) -> Angle {
        self.angles[i]
    }

    pub fn extend_z(&mut self, target: usize, angle: f64) {
        for (i, chain) in self.chains.iter_mut().enumerate() {
            if i == target {
                chain.z.push(true);
            } else {
                chain.z.push(false);
            }
            chain.x.push(false);
        }
        self.angles.push(Angle::Arbitrary(angle));
    }

    pub fn commutes_with(&self, other: &PauliPolynomial) -> bool {
        let size = self.size();
        assert_eq!(size, other.size());

        let self_length = self.length();
        let other_length = other.length();

        for index_1 in 0..self_length {
            let mut pauli_string = Vec::with_capacity(size);
            for q1 in 0..size {
                pauli_string.push(self.chain(q1).pauli(index_1));
            }
            for index_2 in 0..other_length {
                let other_pauli_string = (0..size).map(|q2| other.chain(q2).pauli(index_2));
                let mut commutes = true;
                for (p1, p2) in zip(&pauli_string, other_pauli_string) {
                    if *p1 == PauliLetter::I || p2 == PauliLetter::I || p1 == &p2 {
                        continue;
                    }
                    commutes = !commutes;
                }
                if !commutes {
                    return false;
                }
            }
        }
        true
    }
}

impl PropagateClifford for PauliPolynomial {
    fn cx(&mut self, control: IndexType, target: IndexType) -> &mut Self {
        let mut bit_mask: BitVec = BitVec::repeat(true, self.length());

        let [control, target] = self.chains.get_disjoint_mut([control, target]).unwrap();

        bit_mask ^= &control.z;
        bit_mask ^= &target.x;
        bit_mask &= &control.x;
        bit_mask &= &target.z;

        super::pauli_string::cx(control, target);
        for (angle, flip) in zip(self.angles.iter_mut(), bit_mask.iter()) {
            if *flip {
                angle.flip();
            }
        }

        self
    }

    fn s(&mut self, target: IndexType) -> &mut Self {
        let chains_target = self.chains.get_mut(target).unwrap();
        // Update angles
        let y_vec = chains_target.y_bitmask();
        for (angle, flip) in zip(self.angles.iter_mut(), y_vec.iter()) {
            if *flip {
                angle.flip();
            }
        }
        chains_target.s();
        self
    }

    fn v(&mut self, target: IndexType) -> &mut Self {
        let chains_target = self.chains.get_mut(target).unwrap();
        chains_target.v();
        // Update angles
        let y_vec = chains_target.y_bitmask();
        for (angle, flip) in zip(self.angles.iter_mut(), y_vec.iter()) {
            if *flip {
                angle.flip();
            }
        }
        self
    }
}

impl MaskedPropagateClifford for PauliPolynomial {
    fn masked_cx(&mut self, control: IndexType, target: IndexType, mask: &BitVec) -> &mut Self {
        let mut bit_mask = BitVec::repeat(true, self.length());
        let [control, target] = self.chains.get_disjoint_mut([control, target]).unwrap();

        bit_mask ^= &control.z;
        bit_mask ^= &target.x;
        bit_mask &= &control.x;
        bit_mask &= &target.z;
        bit_mask &= mask;

        super::pauli_string::masked_cx(control, target, mask);
        for (angle, flip) in zip(self.angles.iter_mut(), bit_mask.iter()) {
            if *flip {
                angle.flip();
            }
        }

        self
    }

    fn masked_s(&mut self, target: IndexType, mask: &BitVec) -> &mut Self {
        let chains_target = &mut self.chains[target];

        // Update angles
        let y_vec = chains_target.masked_y_bitmask(mask);
        for (angle, flip) in zip(self.angles.iter_mut(), y_vec.iter()) {
            if *flip {
                angle.flip();
            }
        }
        chains_target.masked_s(mask);
        self
    }

    fn masked_v(&mut self, target: IndexType, mask: &BitVec) -> &mut Self {
        let chains_target = &mut self.chains[target];
        chains_target.masked_v(mask);
        // Update angles
        let y_vec = chains_target.masked_y_bitmask(mask);
        for (angle, flip) in zip(self.angles.iter_mut(), y_vec.iter()) {
            if *flip {
                angle.flip();
            }
        }
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use itertools::Itertools;

    impl PartialEq for PauliPolynomial {
        fn eq(&self, other: &Self) -> bool {
            self.chains == other.chains && self.angles == other.angles
        }
    }

    #[test]
    fn test_pauli_polynomial_constructor() {
        let size = 3;
        let ham = vec![
            ("IXYZ", Angle::from_angle(0.3)),
            ("XXII", Angle::from_angle(0.7)),
            ("YYII", Angle::from_angle(0.12)),
        ];
        let pp = PauliPolynomial::from_hamiltonian(ham);

        let pg1_ref = PauliString::from_text("IXY");
        let pg2_ref = PauliString::from_text("XXY");
        let pg3_ref = PauliString::from_text("YII");
        let pg4_ref = PauliString::from_text("ZII");

        let angles_ref = Angle::from_angles(&[0.3, 0.7, 0.12]);

        let pp_ref = PauliPolynomial {
            chains: vec![pg1_ref, pg2_ref, pg3_ref, pg4_ref],
            angles: angles_ref,
            size,
        };
        assert_eq!(pp, pp_ref);
    }

    #[test]
    #[should_panic]
    fn test_pauli_polynomial_constructor_empty_hamiltonian() {
        let ham = vec![];
        let _ = PauliPolynomial::from_hamiltonian(ham);
    }

    #[test]
    #[should_panic]
    fn test_pauli_polynomial_constructor_unequal_strings() {
        let ham = vec![
            ("IXYZ", Angle::from_angle(0.3)),
            ("XXI", Angle::from_angle(0.7)),
            ("YYII", Angle::from_angle(0.12)),
        ];
        let _ = PauliPolynomial::from_hamiltonian(ham);
    }

    fn setup_sample_pp() -> PauliPolynomial {
        let size = 3;
        let pg1_ref = PauliString::from_text("IXY");
        let pg2_ref = PauliString::from_text("ZYX");
        let pg3_ref = PauliString::from_text("YIX");
        let angles_ref = Angle::from_angles(&[0.3, 0.7, 0.12]);
        PauliPolynomial {
            chains: vec![pg1_ref, pg2_ref, pg3_ref],
            angles: angles_ref,
            size,
        }
    }

    #[test]
    fn test_pauli_polynomial_s() {
        // Polynomials: IZY, XYI, YXX
        let size = 3;
        let mut pp = setup_sample_pp();

        // Apply S to qubits 0 and 1.
        pp.s(0);
        pp.s(1);

        // Polynomials: IZY, -YXI, -XYX

        // IXY -> IY(-X)
        let pg1_ref = PauliString::from_text("IYX");
        // ZYX -> Z(-X)Y
        let pg2_ref = PauliString::from_text("ZXY");
        // YIX
        let pg3_ref = PauliString::from_text("YIX");
        let angles_ref = Angle::from_angles(&[0.3, -0.7, -0.12]);
        let pp_ref = PauliPolynomial {
            chains: vec![pg1_ref, pg2_ref, pg3_ref],
            angles: angles_ref,
            size,
        };
        assert_eq!(pp, pp_ref);
    }

    #[test]
    fn test_pauli_polynomial_v() {
        // Polynomials: IZY, XYI, YXX
        let size = 3;
        let mut pp = setup_sample_pp();

        // Apply V to qubits 0 and 1.
        pp.v(1);
        pp.v(2);

        // IXY
        let pg1_ref = PauliString::from_text("IXY");
        // ZYX -> (-Y)ZX
        let pg2_ref = PauliString::from_text("YZX");
        // YIX -> ZIX
        let pg3_ref = PauliString::from_text("ZIX");
        let angles_ref = Angle::from_angles(&[-0.3, 0.7, 0.12]);
        let pp_ref = PauliPolynomial {
            chains: vec![pg1_ref, pg2_ref, pg3_ref],
            angles: angles_ref,
            size,
        };
        assert_eq!(pp, pp_ref);
    }

    #[test]
    fn test_pauli_polynomial_s_dgr() {
        // Polynomials: IZY, XYI, YXX
        let size = 3;
        let mut pp = setup_sample_pp();

        // Apply Sdgr to qubits 1 and 2.
        pp.s_dgr(1);
        pp.s_dgr(2);

        // IXY
        let pg1_ref = PauliString::from_text("IXY");
        // ZYX -> ZX(-Y)
        let pg2_ref = PauliString::from_text("ZXY");
        // YIX -> XI(-Y)
        let pg3_ref = PauliString::from_text("XIY");
        let angles_ref = Angle::from_angles(&[0.3, 0.7, 0.12]);
        let pp_ref = PauliPolynomial {
            chains: vec![pg1_ref, pg2_ref, pg3_ref],
            angles: angles_ref,
            size,
        };
        assert_eq!(pp, pp_ref);
    }

    #[test]
    fn test_pauli_polynomial_v_dgr() {
        // Polynomials: IZY, XYI, YXX
        let size = 3;
        let mut pp = setup_sample_pp();

        // Apply Vdgr to qubits 1 and 2.
        pp.v_dgr(1);
        pp.v_dgr(2);

        // IXY
        let pg1_ref = PauliString::from_text("IXY");
        // ZYX -> Y(-Z)X
        let pg2_ref = PauliString::from_text("YZX");
        // YIX -> (-Z)IX
        let pg3_ref = PauliString::from_text("ZIX");
        let angles_ref = Angle::from_angles(&[-0.3, -0.7, 0.12]);
        let pp_ref = PauliPolynomial {
            chains: vec![pg1_ref, pg2_ref, pg3_ref],
            angles: angles_ref,
            size,
        };
        assert_eq!(pp, pp_ref);
    }

    #[test]
    fn test_pauli_polynomial_h() {
        // Polynomials: IZY, XYI, YXX
        let size = 3;
        let mut pp = setup_sample_pp();

        // Apply H to qubits 0 and 1.
        pp.h(0);
        pp.h(1);

        // IXY -> IZ(-Y)
        let pg1_ref = PauliString::from_text("IZY");
        // ZYX -> X(-Y)Z
        let pg2_ref = PauliString::from_text("XYZ");
        // YIX -
        let pg3_ref = PauliString::from_text("YIX");
        let angles_ref = Angle::from_angles(&[0.3, -0.7, -0.12]);
        let pp_ref = PauliPolynomial {
            chains: vec![pg1_ref, pg2_ref, pg3_ref],
            angles: angles_ref,
            size,
        };
        assert_eq!(pp, pp_ref);
    }

    fn setup_sample_two_qubit_pp(pauli_letter: char) -> PauliPolynomial {
        let size = 3;
        let pg1_ref = match pauli_letter {
            'i' => PauliString::from_text("IIII"),
            'x' => PauliString::from_text("XXXX"),
            'y' => PauliString::from_text("YYYY"),
            'z' => PauliString::from_text("ZZZZ"),
            _ => panic!("Pauli letter not recognized"),
        };

        let pg2_ref = PauliString::from_text("IXYZ");
        let pg3_ref = PauliString::from_text("YIXZ");

        let angles_ref = Angle::from_angles(&[0.3, 0.7, 0.12, 0.15]);

        PauliPolynomial {
            chains: vec![pg1_ref, pg2_ref, pg3_ref],
            angles: angles_ref,
            size,
        }
    }

    #[test]
    fn test_pauli_polynomial_cx_i() {
        let size = 3;
        // Setup:
        // q0 -> IIII
        // q1 -> IXYZ
        let mut pp = setup_sample_two_qubit_pp('i');

        pp.cx(0, 1);

        // IIII -> IIZZ
        // IXYZ -> IXYZ

        let pg1_ref = PauliString::from_text("IIZZ");
        let pg2_ref = PauliString::from_text("IXYZ");
        // YIXZ
        let pg3_ref = PauliString::from_text("YIXZ");
        // [1, 1, 1, 1]
        let angles_ref = Angle::from_angles(&[0.3, 0.7, 0.12, 0.15]);
        let pp_ref = PauliPolynomial {
            chains: vec![pg1_ref, pg2_ref, pg3_ref],
            angles: angles_ref,
            size,
        };
        assert_eq!(pp, pp_ref);
    }

    #[test]
    fn test_pauli_polynomial_cx_x() {
        let size = 3;
        // Setup:
        // q0 -> XXXX
        // q1 -> IXYZ
        let mut pp = setup_sample_two_qubit_pp('x');

        pp.cx(0, 1);

        // XXXX -> XXYY
        // IXYZ -> XIZY
        let pg1_ref = PauliString::from_text("XXYY");
        let pg2_ref = PauliString::from_text("XIZY");
        // YIXZ
        let pg3_ref = PauliString::from_text("YIXZ");
        // [1, 1, 1, -1]
        let angles_ref = Angle::from_angles(&[0.3, 0.7, 0.12, -0.15]);
        let pp_ref = PauliPolynomial {
            chains: vec![pg1_ref, pg2_ref, pg3_ref],
            angles: angles_ref,
            size,
        };
        assert_eq!(pp, pp_ref);
    }

    #[test]
    fn test_pauli_polynomial_cx_y() {
        let size = 3;
        // Setup:
        // q0 -> YYYY
        // q1 -> IXYZ
        let mut pp = setup_sample_two_qubit_pp('y');

        pp.cx(0, 1);

        // YYYY -> YYXX
        // IXYZ -> XIZY
        let pg1_ref = PauliString::from_text("YYXX");
        let pg2_ref = PauliString::from_text("XIZY");
        // YIXZ
        let pg3_ref = PauliString::from_text("YIXZ");
        // [1, 1, -1, 1]
        let angles_ref = Angle::from_angles(&[0.3, 0.7, -0.12, 0.15]);
        let pp_ref = PauliPolynomial {
            chains: vec![pg1_ref, pg2_ref, pg3_ref],
            angles: angles_ref,
            size,
        };
        assert_eq!(pp, pp_ref);
    }

    #[test]
    fn test_pauli_polynomial_cx_z() {
        let size = 3;
        // Setup:
        // q0 -> ZZZZ
        // q1 -> IXYZ
        let mut pp = setup_sample_two_qubit_pp('z');

        pp.cx(0, 1);

        // ZZZZ -> ZZII
        // IXYZ -> IXYZ
        let pg1_ref = PauliString::from_text("ZZII");
        let pg2_ref = PauliString::from_text("IXYZ");
        // YIXZ
        let pg3_ref = PauliString::from_text("YIXZ");
        // [1, 1, 1, 1]
        let angles_ref = Angle::from_angles(&[0.3, 0.7, 0.12, 0.15]);
        let pp_ref = PauliPolynomial {
            chains: vec![pg1_ref, pg2_ref, pg3_ref],
            angles: angles_ref,
            size,
        };
        assert_eq!(pp, pp_ref);
    }

    #[test]
    fn test_pauli_polynomial_cz_i() {
        let size = 3;
        // Setup:
        // q0 -> IIII
        // q1 -> IXYZ
        let mut pp = setup_sample_two_qubit_pp('i');

        pp.cz(0, 1);

        // IIII -> IZZI
        // IXYZ -> IXYZ

        let pg1_ref = PauliString::from_text("IZZI");
        let pg2_ref = PauliString::from_text("IXYZ");
        // YIXZ
        let pg3_ref = PauliString::from_text("YIXZ");
        // [1, 1, 1, 1]
        let angles_ref = Angle::from_angles(&[0.3, 0.7, 0.12, 0.15]);
        let pp_ref = PauliPolynomial {
            chains: vec![pg1_ref, pg2_ref, pg3_ref],
            angles: angles_ref,
            size,
        };
        assert_eq!(pp, pp_ref);
    }

    #[test]
    fn test_pauli_polynomial_cz_x() {
        let size = 3;
        // Setup:
        // q0 -> XXXX
        // q1 -> IXYZ
        let mut pp = setup_sample_two_qubit_pp('x');

        pp.cz(0, 1);

        // XXXX -> XYYX
        // IXYZ -> ZYXI
        let pg1_ref = PauliString::from_text("XYYX");
        let pg2_ref = PauliString::from_text("ZYXI");
        // YIXZ
        let pg3_ref = PauliString::from_text("YIXZ");
        // [1, 1, -1, 1]
        let angles_ref = Angle::from_angles(&[0.3, 0.7, -0.12, 0.15]);
        let pp_ref = PauliPolynomial {
            chains: vec![pg1_ref, pg2_ref, pg3_ref],
            angles: angles_ref,
            size,
        };
        assert_eq!(pp, pp_ref);
    }

    #[test]
    fn test_pauli_polynomial_cz_y() {
        let size = 3;
        // Setup:
        // q0 -> YYYY
        // q1 -> IXYZ
        let mut pp = setup_sample_two_qubit_pp('y');

        pp.cz(0, 1);

        // YYYY -> YXXY
        // IXYZ -> ZYXI
        let pg1_ref = PauliString::from_text("YXXY");
        let pg2_ref = PauliString::from_text("ZYXI");
        // YIXZ
        let pg3_ref = PauliString::from_text("YIXZ");
        // [1, -1, 1, 1]
        let angles_ref = Angle::from_angles(&[0.3, -0.7, 0.12, 0.15]);
        let pp_ref = PauliPolynomial {
            chains: vec![pg1_ref, pg2_ref, pg3_ref],
            angles: angles_ref,
            size,
        };
        assert_eq!(pp, pp_ref);
    }

    #[test]
    fn test_pauli_polynomial_cz_z() {
        let size = 3;
        // Setup:
        // q0 -> ZZZZ
        // q1 -> IXYZ
        let mut pp = setup_sample_two_qubit_pp('z');

        pp.cz(0, 1);

        // ZZZZ -> ZIIZ
        // IXYZ -> IXYZ
        let pg1_ref = PauliString::from_text("ZIIZ");
        let pg2_ref = PauliString::from_text("IXYZ");
        // YIXZ
        let pg3_ref = PauliString::from_text("YIXZ");
        // [1, 1, 1, 1]
        let angles_ref = Angle::from_angles(&[0.3, 0.7, 0.12, 0.15]);
        let pp_ref = PauliPolynomial {
            chains: vec![pg1_ref, pg2_ref, pg3_ref],
            angles: angles_ref,
            size,
        };
        assert_eq!(pp, pp_ref);
    }

    #[test]
    fn test_commutes_with_simple() {
        let pp1s = vec![
            vec![("I", Angle::from_angle(0.3))],
            vec![("X", Angle::from_angle(0.5))],
            vec![("Y", Angle::from_angle(0.7))],
            vec![("Z", Angle::from_angle(0.9))],
        ]
        .into_iter()
        .map(|ham| PauliPolynomial::from_hamiltonian(ham))
        .collect::<Vec<_>>();

        let pp2s = pp1s.clone();

        for (i, (pp1, pp2)) in pp1s.iter().cartesian_product(pp2s.iter()).enumerate() {
            if i <= 5 || i == 8 || i == 10 || i == 12 || i == 15 {
                assert!(pp1.commutes_with(pp2));
            } else {
                assert!(!pp1.commutes_with(pp2));
            }
        }
    }

    #[test]
    fn test_commutes_with() {
        let pp1 = PauliPolynomial::from_hamiltonian(vec![
            ("IYYX", Angle::from_angle(0.3)),
            ("XXXI", Angle::from_angle(0.5)),
        ]);

        let pp2 = PauliPolynomial::from_hamiltonian(vec![
            ("IYZZ", Angle::from_angle(0.7)),
            ("ZZXI", Angle::from_angle(0.9)),
        ]);

        assert!(pp1.commutes_with(&pp2));
    }

    #[test]
    fn test_not_commutes_with() {
        let pp1 = PauliPolynomial::from_hamiltonian(vec![
            ("IYYX", Angle::from_angle(0.3)),
            ("XXXI", Angle::from_angle(0.5)),
        ]);

        let pp2 = PauliPolynomial::from_hamiltonian(vec![
            ("IYZZ", Angle::from_angle(0.7)),
            ("ZZXY", Angle::from_angle(0.9)),
        ]);

        assert!(!pp1.commutes_with(&pp2));
    }
}
