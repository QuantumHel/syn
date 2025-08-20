use std::{iter::zip, sync::RwLock};

use bitvec::vec::BitVec;
use itertools::zip_eq;

use super::{pauli_string::PauliString, IndexType, MaskedPropagateClifford, PropagateClifford};

// todo: Make this into a union / type Angle
type Angle = f64;

#[derive(Debug, Default)]
pub struct PauliPolynomial {
    chains: Vec<PauliString>,
    angles: RwLock<Vec<Angle>>,
    size: usize,
}

impl Clone for PauliPolynomial {
    fn clone(&self) -> Self {
        PauliPolynomial {
            chains: self.chains.clone(),
            angles: RwLock::new(self.angles.read().unwrap().clone()),
            size: self.size,
        }
    }
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
            .map(|gadget| (PauliString::from_text(gadget)))
            .collect::<Vec<_>>();

        PauliPolynomial {
            chains,
            angles: RwLock::new(angles),
            size: num_qubits,
        }
    }

    pub fn size(&self) -> usize {
        self.size
    }

    pub fn length(&self) -> usize {
        self.angles.read().unwrap().len()
    }

    pub fn chains(&self) -> &Vec<PauliString> {
        &self.chains
    }

    pub fn angle(&self, i: usize) -> Angle {
        self.angles.read().unwrap()[i]
    }
}

impl PropagateClifford for PauliPolynomial {
    fn cx(&mut self, control: IndexType, target: IndexType) -> &mut Self {
        let mut bit_mask: BitVec = BitVec::repeat(true, self.length());

        let [control, target] = self.chains.get_disjoint_mut([control, target]).unwrap();

        bit_mask ^= control.z.read().unwrap().as_bitslice();
        bit_mask ^= target.x.read().unwrap().as_bitslice();
        bit_mask &= control.x.read().unwrap().as_bitslice();
        bit_mask &= target.z.read().unwrap().as_bitslice();

        super::pauli_string::cx(control, target);
        for (angle, flip) in zip(self.angles.write().unwrap().iter_mut(), bit_mask.iter()) {
            if *flip {
                *angle *= -1.0;
            }
        }

        self
    }

    fn s(&mut self, target: IndexType) -> &mut Self {
        let chains_target = self.chains.get_mut(target).unwrap();
        // Update angles
        let y_vec = chains_target.y_bitmask();
        for (angle, flip) in zip(self.angles.write().unwrap().iter_mut(), y_vec.iter()) {
            if *flip {
                *angle *= -1.0;
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
        for (angle, flip) in zip(self.angles.write().unwrap().iter_mut(), y_vec.iter()) {
            if *flip {
                *angle *= -1.0;
            }
        }
        self
    }
}

impl MaskedPropagateClifford for PauliPolynomial {
    fn masked_cx(&self, control: IndexType, target: IndexType, mask: &BitVec) -> &Self {
        let mut bit_mask = BitVec::repeat(true, self.angles.read().unwrap().len());
        // let [control, target] = self.chains.get_many([control, target]).unwrap();
        let control = self.chains.get(control).unwrap();
        let target = self.chains.get(target).unwrap();

        bit_mask ^= control.z.read().unwrap().as_bitslice();
        bit_mask ^= target.x.read().unwrap().as_bitslice();
        bit_mask &= control.x.read().unwrap().as_bitslice();
        bit_mask &= target.z.read().unwrap().as_bitslice();
        bit_mask &= mask;

        super::pauli_string::masked_cx(control, target, mask);
        for (angle, flip) in zip(self.angles.write().unwrap().iter_mut(), bit_mask.iter()) {
            if *flip {
                *angle *= -1.0;
            }
        }

        self
    }

    fn masked_s(&self, target: IndexType, mask: &BitVec) -> &Self {
        let chains_target = self.chains.get(target).unwrap();

        // Update angles
        let y_vec = chains_target.masked_y_bitmask(mask);
        for (angle, flip) in zip(self.angles.write().unwrap().iter_mut(), y_vec.iter()) {
            if *flip {
                *angle *= -1.0;
            }
        }
        chains_target.masked_s(mask);
        self
    }

    fn masked_v(&self, target: IndexType, mask: &BitVec) -> &Self {
        let chains_target = self.chains.get(target).unwrap();
        chains_target.masked_v(mask);
        // Update angles
        let y_vec = chains_target.masked_y_bitmask(mask);
        for (angle, flip) in zip(self.angles.write().unwrap().iter_mut(), y_vec.iter()) {
            if *flip {
                *angle *= -1.0;
            }
        }
        self
    }
}
#[cfg(test)]
mod tests {
    use super::*;

    impl PartialEq for PauliPolynomial {
        fn eq(&self, other: &Self) -> bool {
            self.chains == other.chains
                && *self.angles.read().unwrap() == *other.angles.read().unwrap()
        }
    }

    #[test]
    fn from_hamiltonian() {
        let size = 3;
        let ham = vec![("IXYZ", 0.3), ("XXII", 0.7), ("YYII", 0.12)];
        let pp = PauliPolynomial::from_hamiltonian(ham);

        let pg1_ref = PauliString::from_text("IXY");
        let pg2_ref = PauliString::from_text("XXY");
        let pg3_ref = PauliString::from_text("YII");
        let pg4_ref = PauliString::from_text("ZII");

        let angles_ref = vec![0.3, 0.7, 0.12];

        let pp_ref = PauliPolynomial {
            chains: vec![pg1_ref, pg2_ref, pg3_ref, pg4_ref],
            angles: RwLock::new(angles_ref),
            size,
        };
        assert_eq!(pp, pp_ref);
    }

    #[test]
    #[should_panic]
    fn from_hamiltonian_empty_hamiltonian() {
        let ham = vec![];
        let _ = PauliPolynomial::from_hamiltonian(ham);
    }

    #[test]
    #[should_panic]
    fn from_hamiltonian_unequal_strings() {
        let ham = vec![("IXYZ", 0.3), ("XXI", 0.7), ("YYII", 0.12)];
        let _ = PauliPolynomial::from_hamiltonian(ham);
    }

    fn setup_sample_pp() -> PauliPolynomial {
        let size = 3;
        let pg1_ref = PauliString::from_text("IXY");
        let pg2_ref = PauliString::from_text("ZYX");
        let pg3_ref = PauliString::from_text("YIX");
        let angles_ref = vec![0.3, 0.7, 0.12];
        PauliPolynomial {
            chains: vec![pg1_ref, pg2_ref, pg3_ref],
            angles: RwLock::new(angles_ref),
            size,
        }
    }

    #[test]
    fn propagate_s() {
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
        let angles_ref = vec![0.3, -0.7, -0.12];
        let pp_ref = PauliPolynomial {
            chains: vec![pg1_ref, pg2_ref, pg3_ref],
            angles: RwLock::new(angles_ref),
            size,
        };
        assert_eq!(pp, pp_ref);
    }

    #[test]
    fn propagate_v() {
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
        let angles_ref = vec![-0.3, 0.7, 0.12];
        let pp_ref = PauliPolynomial {
            chains: vec![pg1_ref, pg2_ref, pg3_ref],
            angles: RwLock::new(angles_ref),
            size,
        };
        assert_eq!(pp, pp_ref);
    }

    #[test]
    fn propagate_s_dgr() {
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
        let angles_ref = vec![0.3, 0.7, 0.12];
        let pp_ref = PauliPolynomial {
            chains: vec![pg1_ref, pg2_ref, pg3_ref],
            angles: RwLock::new(angles_ref),
            size,
        };
        assert_eq!(pp, pp_ref);
    }

    #[test]
    fn propagate_v_dgr() {
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
        let angles_ref = vec![-0.3, -0.7, 0.12];
        let pp_ref = PauliPolynomial {
            chains: vec![pg1_ref, pg2_ref, pg3_ref],
            angles: RwLock::new(angles_ref),
            size,
        };
        assert_eq!(pp, pp_ref);
    }

    #[test]
    fn propagate_h() {
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
        let angles_ref = vec![0.3, -0.7, -0.12];
        let pp_ref = PauliPolynomial {
            chains: vec![pg1_ref, pg2_ref, pg3_ref],
            angles: RwLock::new(angles_ref),
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

        let angles_ref = vec![0.3, 0.7, 0.12, 0.15];

        PauliPolynomial {
            chains: vec![pg1_ref, pg2_ref, pg3_ref],
            angles: RwLock::new(angles_ref),
            size,
        }
    }

    #[test]
    fn propagate_cx_i() {
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
        let angles_ref = vec![0.3, 0.7, 0.12, 0.15];
        let pp_ref = PauliPolynomial {
            chains: vec![pg1_ref, pg2_ref, pg3_ref],
            angles: RwLock::new(angles_ref),
            size,
        };
        assert_eq!(pp, pp_ref);
    }

    #[test]
    fn propagate_cx_x() {
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
        let angles_ref = vec![0.3, 0.7, 0.12, -0.15];
        let pp_ref = PauliPolynomial {
            chains: vec![pg1_ref, pg2_ref, pg3_ref],
            angles: RwLock::new(angles_ref),
            size,
        };
        assert_eq!(pp, pp_ref);
    }

    #[test]
    fn propagate_cx_y() {
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
        let angles_ref = vec![0.3, 0.7, -0.12, 0.15];
        let pp_ref = PauliPolynomial {
            chains: vec![pg1_ref, pg2_ref, pg3_ref],
            angles: RwLock::new(angles_ref),
            size,
        };
        assert_eq!(pp, pp_ref);
    }

    #[test]
    fn propagate_cx_z() {
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
        let angles_ref = vec![0.3, 0.7, 0.12, 0.15];
        let pp_ref = PauliPolynomial {
            chains: vec![pg1_ref, pg2_ref, pg3_ref],
            angles: RwLock::new(angles_ref),
            size,
        };
        assert_eq!(pp, pp_ref);
    }

    #[test]
    fn propagate_cz_i() {
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
        let angles_ref = vec![0.3, 0.7, 0.12, 0.15];
        let pp_ref = PauliPolynomial {
            chains: vec![pg1_ref, pg2_ref, pg3_ref],
            angles: RwLock::new(angles_ref),
            size,
        };
        assert_eq!(pp, pp_ref);
    }

    #[test]
    fn propagate_cz_x() {
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
        let angles_ref = vec![0.3, 0.7, -0.12, 0.15];
        let pp_ref = PauliPolynomial {
            chains: vec![pg1_ref, pg2_ref, pg3_ref],
            angles: RwLock::new(angles_ref),
            size,
        };
        assert_eq!(pp, pp_ref);
    }

    #[test]
    fn propagate_cz_y() {
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
        let angles_ref = vec![0.3, -0.7, 0.12, 0.15];
        let pp_ref = PauliPolynomial {
            chains: vec![pg1_ref, pg2_ref, pg3_ref],
            angles: RwLock::new(angles_ref),
            size,
        };
        assert_eq!(pp, pp_ref);
    }

    #[test]
    fn propagate_cz_z() {
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
        // [1, 1, -1, 1]
        let angles_ref = vec![0.3, 0.7, 0.12, 0.15];
        let pp_ref = PauliPolynomial {
            chains: vec![pg1_ref, pg2_ref, pg3_ref],
            angles: RwLock::new(angles_ref),
            size,
        };
        assert_eq!(pp, pp_ref);
    }
}
