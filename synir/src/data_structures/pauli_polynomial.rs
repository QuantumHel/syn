// use std::iter::zip;

use bitvec::vec::BitVec;
use itertools::zip_eq;
// use itertools::Itertools;
use std::fmt;
use std::{iter::zip, sync::RwLock};

use super::{pauli_string::PauliString, IndexType, MaskedPropagateClifford, PropagateClifford};

// todo: Make this into a union / type Angle
type Angle = f64;

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

    pub fn get_line_string(&self, i: usize) -> String {
        let mut out = String::new();
        let chain_str = self.chains[i].to_string();
        for ch in chain_str.chars() {
            out.push(ch);
            if !ch.is_whitespace() {
                out.push_str("     |");
            }
        }
        out
    }

    pub fn get_first_line_string(&self) -> String {
        let mut out = String::new();
        for angle in self.angles.iter() {
            out.push_str(&format!(" {:.3}", angle)); //force 3 decimal place for formatting
            out.push_str(" |");
        }
        out
    }

    pub fn empty(i: usize) -> Self {
        PauliPolynomial {
            chains: vec![],
            angles: vec![],
            size: i,
        }
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
                *angle *= -1.0;
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
        for (angle, flip) in zip(self.angles.iter_mut(), y_vec.iter()) {
            if *flip {
                *angle *= -1.0;
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
                *angle *= -1.0;
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
                *angle *= -1.0;
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
                *angle *= -1.0;
            }
        }
        self
    }
}

impl fmt::Display for PauliPolynomial {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut out = String::new();
        if self.angles.is_empty() {
            out.push_str("Angles || None|\n");
            for i in 0..self.size() {
                out.push_str("QB");
                out.push_str(&i.to_string());
                out.push_str("    || ");
                out.push_str("_   |\n");
            }
            writeln!(f, "{}", out)?;
        } else {
            // write first line
            out.push_str("Angles ||"); // I take this out from get_first_line_string because I want to reuse that function for pauli exponential
            out.push_str(&self.get_first_line_string());
            writeln!(f, "{}", out)?;

            // write subsequent lines
            let chains = self.chains();
            for (i, _) in chains.iter().enumerate() {
                let mut out = String::new();
                out.push_str("QB");
                out.push_str(&i.to_string());
                out.push_str("    || ");
                out.push_str(&self.get_line_string(i));
                writeln!(f, "{}", out)?;
            }
        }

        writeln!(f)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    impl PartialEq for PauliPolynomial {
        fn eq(&self, other: &Self) -> bool {
            self.chains == other.chains && self.angles == other.angles
        }
    }

    #[test]
    fn test_pauli_polynomial_constructor() {
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
        let angles_ref = vec![0.3, -0.7, -0.12];
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
        let angles_ref = vec![-0.3, 0.7, 0.12];
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
        let angles_ref = vec![0.3, 0.7, 0.12];
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
        let angles_ref = vec![-0.3, -0.7, 0.12];
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
        let angles_ref = vec![0.3, -0.7, -0.12];
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

        let angles_ref = vec![0.3, 0.7, 0.12, 0.15];

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
        let angles_ref = vec![0.3, 0.7, 0.12, 0.15];
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
        let angles_ref = vec![0.3, 0.7, 0.12, -0.15];
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
        let angles_ref = vec![0.3, 0.7, -0.12, 0.15];
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
        let angles_ref = vec![0.3, 0.7, 0.12, 0.15];
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
        let angles_ref = vec![0.3, 0.7, 0.12, 0.15];
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
        let angles_ref = vec![0.3, 0.7, -0.12, 0.15];
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
        let angles_ref = vec![0.3, -0.7, 0.12, 0.15];
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
        // [1, 1, -1, 1]
        let angles_ref = vec![0.3, 0.7, 0.12, 0.15];
        let pp_ref = PauliPolynomial {
            chains: vec![pg1_ref, pg2_ref, pg3_ref],
            angles: angles_ref,
            size,
        };
        assert_eq!(pp, pp_ref);
    }
    #[test]
    fn test_pauli_polynomial_display() {
        let pp = setup_sample_pp();
        assert_eq!(
            pp.to_string(),
            "Angles || 0.300 | 0.700 | 0.120 |\nQB0    || I     | X     | Y     |\nQB1    || Z     | Y     | X     |\nQB2    || Y     | I     | X     |\n\n"
        );
    }
}
