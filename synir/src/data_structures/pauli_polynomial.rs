use std::borrow::{Borrow, BorrowMut};
use std::iter::zip;

use super::{pauli_string::PauliString, IndexType, MaskedPropagateClifford, PropagateClifford};
use crate::data_structures::{Angle, PauliLetter};
use bitvec::vec::BitVec;
use itertools::{zip_eq, Itertools};

pub mod simplify;

#[derive(Debug, Clone, Default)]
pub struct PauliPolynomial {
    chains: Vec<PauliString>,
    angles: Vec<Angle>,
}

impl PauliPolynomial {
    pub fn new(num_qubits: usize) -> Self {
        PauliPolynomial::from_components(
            (0..num_qubits)
                .map(|_| PauliString::from_text(""))
                .collect_vec(),
            vec![],
        )
    }

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

        Self { chains, angles }
    }

    pub fn from_components(chains: Vec<PauliString>, angles: Vec<Angle>) -> Self {
        assert!(
            chains.len() > 0,
            "Cannot construct PauliPolynomial from empty chains; unknown number of qubits."
        );
        for i in 0..chains.len() {
            assert!(
                chains[0].len() == angles.len(),
                "Each PauliPolynomial chain need to have the same length als angles."
            );
        }
        Self { chains, angles }
    }

    pub fn from_pauli_gadgets(gadgets: Vec<PauliString>, angles: Vec<Angle>) -> Self {
        assert!(
            gadgets.len() > 0,
            "Cannot construct PauliPolynomial from empty gadgets; unknown number of qubits."
        );
        assert!(
            gadgets.len() == angles.len(),
            "PauliPolynomials need the same number of gadgets as angles"
        );
        Self::from_components(
            (0..gadgets.len())
                .map(|i| {
                    PauliString::from_letters(
                        gadgets.iter().map(|g| g.pauli(i)).collect_vec().borrow(),
                    )
                })
                .collect_vec(),
            angles,
        )
    }

    pub fn num_qubits(&self) -> usize {
        self.chains.len()
    }

    pub fn length(&self) -> usize {
        self.angles.len()
    }

    pub fn pauli_gadget(&self, i: usize) -> PauliString {
        assert!(i < self.length());
        PauliString::from_letters(&self.chains.iter().map(|ps| ps.pauli(i)).collect_vec())
    }

    pub fn pauli_gadgets(&self) -> Vec<PauliString> {
        (0..self.length())
            .map(|i| self.pauli_gadget(i))
            .collect_vec()
    }

    pub fn chain(&self, i: usize) -> &PauliString {
        assert!(i < self.num_qubits());
        &self.chains[i]
    }

    pub fn chains(&self) -> &Vec<PauliString> {
        &self.chains
    }

    pub fn angles(&self) -> &Vec<Angle> {
        &self.angles
    }

    pub fn mut_angles(&mut self) -> &mut Vec<Angle> {
        &mut self.angles
    }

    pub fn angle(&self, i: usize) -> Angle {
        assert!(i < self.length());
        self.angles[i]
    }

    pub fn mut_chains(&mut self) -> &mut Vec<PauliString> {
        &mut self.chains
    }

    pub fn mut_chains_and_angles(&mut self) -> (&mut Vec<PauliString>, &mut Vec<Angle>) {
        assert!(self.angles.len() == self.chains[0].len(), "PauliPolynomials should always have the same amount of Chains and Angles, but they are {} {}.", self.chains[0].len(), self.angles.len());
        (&mut self.chains, &mut self.angles)
    }

    pub fn extend_z(&mut self, target: usize, angle: f64) {
        assert!(target < self.num_qubits());
        for (i, chain) in self.chains.iter_mut().enumerate() {
            if i == target {
                chain.z.push(true);
            } else {
                chain.z.push(false);
            }
            chain.x.push(false);
        }
        self.angles.push(Angle::from_angle(angle));
    }

    pub fn append_other(&mut self, other: &PauliPolynomial) {
        for (gadget, angle) in zip(other.pauli_gadgets(), other.angles()) {
            self.append_gadget(gadget, *angle);
        }
    }

    pub fn append_gadget(&mut self, gadget: PauliString, angle: Angle) {
        if self.num_qubits() != gadget.len() {
            panic!(
                "Appending Gadget to PauliPolynomial of different size {} {}",
                gadget.len(),
                self.num_qubits()
            );
        }
        self.angles.push(angle);
        for i in 0..self.num_qubits() {
            self.chains[i].x.push(gadget.x(i));
            self.chains[i].z.push(gadget.z(i));
        }
    }

    pub fn remove_gadget(&mut self, i: usize) -> (PauliString, Angle) {
        assert!(i < self.length());
        let angle = self.angles.swap_remove(i);
        let gadget = PauliString::from_letters(
            (0..self.num_qubits())
                .map(|j| self.chains[j].swap_remove(i))
                .collect_vec()
                .borrow(),
        );
        (gadget, angle)
    }

    pub fn commutes_with(&self, other: &PauliPolynomial) -> bool {
        if self.num_qubits() != other.num_qubits() {
            panic!("Commutation checking only works with PauliPolynomials of the same size, but found {} and {}", self.num_qubits(), other.num_qubits());
        }

        for gadget1 in (0..self.length()).map(|i| self.pauli_gadget(i)) {
            if !other.commutes_with_gadget(gadget1) {
                return false;
            }
        }
        true
    }

    pub fn commutes_with_gadget(&self, gadget: PauliString) -> bool {
        if self.num_qubits() != gadget.len() {
            panic!("Commutation checking only works with PauliPolynomials of the same size, but found {} and {}", self.num_qubits(), gadget.len());
        }
        for self_gadget in (0..self.length()).map(|i| self.pauli_gadget(i)) {
            if !self_gadget.commutes_with(&gadget) {
                return false;
            }
        }
        true
    }

    pub fn push_clifford_gadget(&mut self, gadget: PauliString, angle: Angle) {
        self.propate_gadget(gadget, angle).unwrap();
    }

    pub fn push_clifford(&mut self, other: &PauliPolynomial) {
        for (gadget, angle) in zip(other.pauli_gadgets(), other.angles()) {
            self.push_clifford_gadget(gadget, *angle);
        }
    }
}

impl PropagateClifford for PauliPolynomial {
    fn cx(&mut self, control: IndexType, target: IndexType) -> &mut Self {
        self.masked_cx(control, target, &BitVec::repeat(true, self.length()))
    }

    fn s(&mut self, target: IndexType) -> &mut Self {
        self.masked_s(target, &BitVec::repeat(true, self.length()))
    }

    fn v(&mut self, target: IndexType) -> &mut Self {
        self.masked_v(target, &BitVec::repeat(true, self.length()))
    }
}

impl MaskedPropagateClifford for PauliPolynomial {
    fn masked_cx(&mut self, control: IndexType, target: IndexType, mask: &BitVec) -> &mut Self {
        //let mut bit_mask = BitVec::repeat(true, self.length());
        let [control, target] = self.chains.get_disjoint_mut([control, target]).unwrap();
        self.angles = self
            .angles
            .iter_mut()
            .enumerate()
            .map(|(i, angle)| {
                match (control.pauli(i), target.pauli(i)) {
                    (PauliLetter::X, PauliLetter::Z) => angle.flip(),
                    (PauliLetter::Y, PauliLetter::Y) => angle.flip(),
                    _ => (),
                };
                //println!("{} {} {}", control.pauli(i), target.pauli(i), angle);
                *angle
            })
            .collect();
        super::pauli_string::masked_cx(control, target, mask);
        self
    }

    fn masked_s(&mut self, target: IndexType, mask: &BitVec) -> &mut Self {
        let chains_target = &mut self.chains[target];
        self.angles = self
            .angles
            .iter_mut()
            .enumerate()
            .map(|(i, angle)| {
                match chains_target.pauli(i) {
                    PauliLetter::X => angle.flip(),
                    _ => (),
                };
                *angle
            })
            .collect();
        chains_target.masked_s(mask);
        self
    }

    fn masked_v(&mut self, target: IndexType, mask: &BitVec) -> &mut Self {
        let chains_target = &mut self.chains[target];
        self.angles = self
            .angles
            .iter_mut()
            .enumerate()
            .map(|(i, angle)| {
                match chains_target.pauli(i) {
                    PauliLetter::Y => angle.flip(),
                    _ => (),
                };
                *angle
            })
            .collect();
        chains_target.masked_v(mask);
        self
    }
}

#[cfg(test)]
mod tests {
    use std::collections::VecDeque;

    use crate::ir::pauli_exponential::print_pp_help;

    use super::*;
    use itertools::Itertools;

    impl PartialEq for PauliPolynomial {
        fn eq(&self, other: &Self) -> bool {
            self.chains == other.chains && self.angles == other.angles
        }
    }

    #[test]
    fn test_pauli_polynomial_constructor() {
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
        let pg1_ref = PauliString::from_text("IXY");
        let pg2_ref = PauliString::from_text("ZYX");
        let pg3_ref = PauliString::from_text("YIX");
        let angles_ref = vec![
            Angle::from_angle(0.3),
            Angle::from_pi4_rotation(2),
            Angle::from_angle(0.12),
        ];
        let pp = PauliPolynomial {
            chains: vec![pg1_ref, pg2_ref, pg3_ref],
            angles: angles_ref,
        };
        println!("Original");
        print_pp_help(&VecDeque::from(vec![pp.to_owned()]));
        pp
    }

    #[test]
    fn test_pauli_polynomial_s() {
        // Polynomials: IZY, XYI, YIX
        let mut pp = setup_sample_pp();
        pp.s(1);
        pp.s(2);

        // IXY
        let pg1_ref = PauliString::from_text("IXY");
        // ZYX -> ZX(-Y)
        let pg2_ref = PauliString::from_text("ZXY");
        // YIX -> XI(-Y)
        let pg3_ref = PauliString::from_text("XIY");
        let angles_ref = vec![
            Angle::from_angle(0.3),
            Angle::from_pi4_rotation(2),
            Angle::from_angle(0.12),
        ];
        let pp_ref = PauliPolynomial {
            chains: vec![pg1_ref, pg2_ref, pg3_ref],
            angles: angles_ref,
        };
        print_pp_help(&VecDeque::from(vec![pp.to_owned()]));
        print_pp_help(&VecDeque::from(vec![pp_ref.to_owned()]));
        assert_eq!(pp, pp_ref);
    }

    #[test]
    fn test_pauli_polynomial_v() {
        // Polynomials: IZY, XYI, YXX
        let mut pp = setup_sample_pp();

        // Apply V to qubits 0 and 1.
        pp.v(1);
        pp.v(2);

        // IXY
        let pg1_ref = PauliString::from_text("IXY");
        // ZYX -> Y(-Z)X
        let pg2_ref = PauliString::from_text("YZX");
        // YIX -> (-Z)IX
        let pg3_ref = PauliString::from_text("ZIX");
        let angles_ref = vec![
            Angle::from_angle(-0.3),
            Angle::from_pi4_rotation(6),
            Angle::from_angle(0.12),
        ];
        let pp_ref = PauliPolynomial {
            chains: vec![pg1_ref, pg2_ref, pg3_ref],
            angles: angles_ref,
        };
        assert_eq!(pp, pp_ref);
    }

    #[test]
    fn test_pauli_polynomial_s_dgr() {
        // Polynomials: IZY, XYI, YXX
        let mut pp = setup_sample_pp();

        // Apply Sdgr to qubits 1 and 2.
        pp.s_dgr(1);
        pp.s_dgr(2);

        // IXY
        let pg1_ref = PauliString::from_text("IXY");
        // ZYX -> Z(-X)Y
        let pg2_ref = PauliString::from_text("ZXY");
        // YIX -> (-X)IY
        let pg3_ref = PauliString::from_text("XIY");
        let angles_ref = vec![
            Angle::from_angle(-0.3),
            Angle::from_pi4_rotation(6),
            Angle::from_angle(0.12),
        ];
        let pp_ref = PauliPolynomial {
            chains: vec![pg1_ref, pg2_ref, pg3_ref],
            angles: angles_ref,
        };
        assert_eq!(pp, pp_ref);
    }

    #[test]
    fn test_pauli_polynomial_v_dgr() {
        // Polynomials: IZY, XYI, YXX
        let mut pp = setup_sample_pp();

        // Apply Vdgr to qubits 1 and 2.
        pp.v_dgr(1);
        pp.v_dgr(2);

        // IXY
        let pg1_ref = PauliString::from_text("IXY");
        // ZYX -> (-Y)ZX
        let pg2_ref = PauliString::from_text("YZX");
        // YIX -> ZIX
        let pg3_ref = PauliString::from_text("ZIX");
        let angles_ref = vec![
            Angle::from_angle(-0.3),
            Angle::from_pi4_rotation(2),
            Angle::from_angle(0.12),
        ];
        let pp_ref = PauliPolynomial {
            chains: vec![pg1_ref, pg2_ref, pg3_ref],
            angles: angles_ref,
        };
        assert_eq!(pp, pp_ref);
    }

    #[test]
    fn test_pauli_polynomial_h() {
        // Polynomials: IZY, XYI, YXX
        let mut pp = setup_sample_pp();

        // Apply H to qubits 0 and 1.
        //pp.h(0);
        pp.h(1);

        // IXY -> IZ(-Y)
        //let pg1_ref = PauliString::from_text("IZY");
        let pg1_ref = PauliString::from_text("IXY");
        // ZYX -> X(-Y)Z
        let pg2_ref = PauliString::from_text("XYZ");
        // YIX -
        let pg3_ref = PauliString::from_text("YIX");
        let angles_ref = vec![
            Angle::from_angle(0.3),
            Angle::from_pi4_rotation(6),
            Angle::from_angle(0.12),
            //Angle::from_angle(-0.12),
        ];
        let pp_ref = PauliPolynomial {
            chains: vec![pg1_ref, pg2_ref, pg3_ref],
            angles: angles_ref,
        };
        print_pp_help(&VecDeque::from(vec![pp.to_owned()]));
        print_pp_help(&VecDeque::from(vec![pp_ref.to_owned()]));
        assert_eq!(pp, pp_ref);
    }

    fn setup_sample_two_qubit_pp(pauli_letter: char) -> PauliPolynomial {
        let pg1_ref = match pauli_letter {
            'i' => PauliString::from_text("IIII"),
            'x' => PauliString::from_text("XXXX"),
            'y' => PauliString::from_text("YYYY"),
            'z' => PauliString::from_text("ZZZZ"),
            _ => panic!("Pauli letter not recognized"),
        };

        let pg2_ref = PauliString::from_text("IXYZ");
        let pg3_ref = PauliString::from_text("YIXZ");

        let angles_ref = vec![
            Angle::from_angle(0.3),
            Angle::from_pi4_rotation(3),
            Angle::from_pi4_rotation(2),
            Angle::from_angle(0.15),
        ];

        PauliPolynomial {
            chains: vec![pg1_ref, pg2_ref, pg3_ref],
            angles: angles_ref,
        }
    }

    #[test]
    fn test_pauli_polynomial_cx_i() {
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
        let angles_ref = vec![
            Angle::from_angle(0.3),
            Angle::from_pi4_rotation(3),
            Angle::from_pi4_rotation(2),
            Angle::from_angle(0.15),
        ];
        let pp_ref = PauliPolynomial {
            chains: vec![pg1_ref, pg2_ref, pg3_ref],
            angles: angles_ref,
        };
        print_pp_help(&VecDeque::from(vec![pp.to_owned()]));
        print_pp_help(&VecDeque::from(vec![pp_ref.to_owned()]));
        assert_eq!(pp, pp_ref);
    }

    #[test]
    fn test_pauli_polynomial_cx_x() {
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
        let angles_ref = vec![
            Angle::from_angle(0.3),
            Angle::from_pi4_rotation(3),
            Angle::from_pi4_rotation(2),
            Angle::from_angle(-0.15),
        ];
        let pp_ref = PauliPolynomial {
            chains: vec![pg1_ref, pg2_ref, pg3_ref],
            angles: angles_ref,
        };
        print_pp_help(&VecDeque::from(vec![pp.to_owned()]));
        print_pp_help(&VecDeque::from(vec![pp_ref.to_owned()]));
        assert_eq!(pp, pp_ref);
    }

    #[test]
    fn test_pauli_polynomial_cx_y() {
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
        let angles_ref = vec![
            Angle::from_angle(0.3),
            Angle::from_pi4_rotation(3),
            Angle::from_pi4_rotation(6),
            Angle::from_angle(0.15),
        ];
        let pp_ref = PauliPolynomial {
            chains: vec![pg1_ref, pg2_ref, pg3_ref],
            angles: angles_ref,
        };
        print_pp_help(&VecDeque::from(vec![pp.to_owned()]));
        print_pp_help(&VecDeque::from(vec![pp_ref.to_owned()]));
        assert_eq!(pp, pp_ref);
    }

    #[test]
    fn test_pauli_polynomial_cx_z() {
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
        let angles_ref = vec![
            Angle::from_angle(0.3),
            Angle::from_pi4_rotation(3),
            Angle::from_pi4_rotation(2),
            Angle::from_angle(0.15),
        ];
        let pp_ref = PauliPolynomial {
            chains: vec![pg1_ref, pg2_ref, pg3_ref],
            angles: angles_ref,
        };
        print_pp_help(&VecDeque::from(vec![pp.to_owned()]));
        print_pp_help(&VecDeque::from(vec![pp_ref.to_owned()]));
        assert_eq!(pp, pp_ref);
    }

    #[test]
    fn test_pauli_polynomial_cz_i() {
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
        let angles_ref = vec![
            Angle::from_angle(0.3),
            Angle::from_pi4_rotation(3),
            Angle::from_pi4_rotation(2),
            Angle::from_angle(0.15),
        ];
        let pp_ref = PauliPolynomial {
            chains: vec![pg1_ref, pg2_ref, pg3_ref],
            angles: angles_ref,
        };
        print_pp_help(&VecDeque::from(vec![pp.to_owned()]));
        print_pp_help(&VecDeque::from(vec![pp_ref.to_owned()]));
        assert_eq!(pp, pp_ref);
    }

    #[test]
    fn test_pauli_polynomial_cz_x() {
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
        let angles_ref = vec![
            Angle::from_angle(0.3),
            Angle::from_pi4_rotation(3),
            Angle::from_pi4_rotation(6),
            Angle::from_angle(0.15),
        ];
        let pp_ref = PauliPolynomial {
            chains: vec![pg1_ref, pg2_ref, pg3_ref],
            angles: angles_ref,
        };
        print_pp_help(&VecDeque::from(vec![pp.to_owned()]));
        print_pp_help(&VecDeque::from(vec![pp_ref.to_owned()]));
        assert_eq!(pp, pp_ref);
    }

    #[test]
    fn test_pauli_polynomial_cz_y() {
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
        let angles_ref = vec![
            Angle::from_angle(0.3),
            Angle::from_pi4_rotation(5),
            Angle::from_pi4_rotation(2),
            Angle::from_angle(0.15),
        ];
        let pp_ref = PauliPolynomial {
            chains: vec![pg1_ref, pg2_ref, pg3_ref],
            angles: angles_ref,
        };
        print_pp_help(&VecDeque::from(vec![pp.to_owned()]));
        print_pp_help(&VecDeque::from(vec![pp_ref.to_owned()]));
        assert_eq!(pp, pp_ref);
    }

    #[test]
    fn test_pauli_polynomial_cz_z() {
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
        let angles_ref = vec![
            Angle::from_angle(0.3),
            Angle::from_pi4_rotation(3),
            Angle::from_pi4_rotation(2),
            Angle::from_angle(0.15),
        ];
        let pp_ref = PauliPolynomial {
            chains: vec![pg1_ref, pg2_ref, pg3_ref],
            angles: angles_ref,
        };
        print_pp_help(&VecDeque::from(vec![pp.to_owned()]));
        print_pp_help(&VecDeque::from(vec![pp_ref.to_owned()]));
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

    #[test]
    fn test_cliffords() {
        // Check whether the gadgets are properly updated according to https://arxiv.org/pdf/2007.10515
        let mut pp = PauliPolynomial::from_hamiltonian(vec![
            ("IXII", Angle::from_angle(0.234)),
            ("IIYI", Angle::from_angle(0.789)),
            ("IIIZ", Angle::from_angle(0.649)),
            ("XXXX", Angle::from_pi4_rotation(1)),
            ("YYYY", Angle::from_pi4_rotation(3)),
            ("ZZZZ", Angle::from_pi4_rotation(5)),
        ]);
        pp.h(0);
        pp.h(1);
        pp.h(2);
        pp.h(3);
        let pp_ref_h = PauliPolynomial::from_hamiltonian(vec![
            ("IZII", Angle::from_angle(0.234)),
            ("IIYI", Angle::from_angle(-0.789)),
            ("IIIX", Angle::from_angle(0.649)),
            ("ZZZZ", Angle::from_pi4_rotation(1)),
            ("YYYY", Angle::from_pi4_rotation(3)),
            ("XXXX", Angle::from_pi4_rotation(5)),
        ]);
        assert_eq!(pp, pp_ref_h);
        pp.s(0);
        pp.s(1);
        pp.s(2);
        pp.s(3);
        let pp_ref_s = PauliPolynomial::from_hamiltonian(vec![
            ("IZII", Angle::from_angle(0.234)),
            ("IIXI", Angle::from_angle(-0.789)),
            ("IIIY", Angle::from_angle(-0.649)),
            ("ZZZZ", Angle::from_pi4_rotation(1)),
            ("XXXX", Angle::from_pi4_rotation(3)),
            ("YYYY", Angle::from_pi4_rotation(5)),
        ]);
        assert_eq!(pp, pp_ref_s);
        pp.v(0);
        pp.v(1);
        pp.v(2);
        pp.v(3);
        let pp_ref_v = PauliPolynomial::from_hamiltonian(vec![
            ("IYII", Angle::from_angle(0.234)),
            ("IIXI", Angle::from_angle(-0.789)),
            ("IIIZ", Angle::from_angle(0.649)),
            ("YYYY", Angle::from_pi4_rotation(1)),
            ("XXXX", Angle::from_pi4_rotation(3)),
            ("ZZZZ", Angle::from_pi4_rotation(5)),
        ]);
        assert_eq!(pp, pp_ref_v);
        pp.cx(0, 3);
        let pp_ref_cnot1 = PauliPolynomial::from_hamiltonian(vec![
            ("IYII", Angle::from_angle(0.234)),
            ("IIXI", Angle::from_angle(-0.789)),
            ("ZIIZ", Angle::from_angle(0.649)),
            ("XYYZ", Angle::from_pi4_rotation(7)),
            ("XXXI", Angle::from_pi4_rotation(3)),
            ("IZZZ", Angle::from_pi4_rotation(5)),
        ]);
        print_pp_help(&VecDeque::from(vec![pp.to_owned()]));
        print_pp_help(&VecDeque::from(vec![pp_ref_cnot1.to_owned()]));
        assert_eq!(pp, pp_ref_cnot1);
        pp.cx(2, 0);
        let pp_ref_cnot2 = PauliPolynomial::from_hamiltonian(vec![
            ("IYII", Angle::from_angle(0.234)),
            ("XIXI", Angle::from_angle(-0.789)),
            ("ZIZZ", Angle::from_angle(0.649)),
            ("IYYZ", Angle::from_pi4_rotation(7)),
            ("IXXI", Angle::from_pi4_rotation(3)),
            ("IZZZ", Angle::from_pi4_rotation(5)),
        ]);
        assert_eq!(pp, pp_ref_cnot2);
        pp.cx(0, 1);
        let pp_ref_cnot3 = PauliPolynomial::from_hamiltonian(vec![
            ("ZYII", Angle::from_angle(0.234)),
            ("XXXI", Angle::from_angle(-0.789)),
            ("ZIZZ", Angle::from_angle(0.649)),
            ("ZYYZ", Angle::from_pi4_rotation(7)),
            ("IXXI", Angle::from_pi4_rotation(3)),
            ("ZZZZ", Angle::from_pi4_rotation(5)),
        ]);
        assert_eq!(pp, pp_ref_cnot3);
    }

    #[test]
    fn test_cnot_prop() {
        // I in control or target
        let mut pp_i = PauliPolynomial::from_hamiltonian(vec![
            ("IXII", Angle::from_angle(1.)),
            ("IYII", Angle::from_angle(2.)),
            ("IZII", Angle::from_angle(3.)),
            ("IIXI", Angle::from_angle(4.)),
            ("IIYI", Angle::from_angle(5.)),
            ("IIZI", Angle::from_angle(6.)),
        ]);
        pp_i.cx(0, 1);
        pp_i.cx(2, 3);
        let pp_ref_i = PauliPolynomial::from_hamiltonian(vec![
            ("IXII", Angle::from_angle(1.)),
            ("ZYII", Angle::from_angle(2.)),
            ("ZZII", Angle::from_angle(3.)),
            ("IIXX", Angle::from_angle(4.)),
            ("IIYX", Angle::from_angle(5.)),
            ("IIZI", Angle::from_angle(6.)),
        ]);
        assert_eq!(pp_i, pp_ref_i);
        // X in control or target
        let mut pp_x = PauliPolynomial::from_hamiltonian(vec![
            ("XXII", Angle::from_angle(1.)),
            ("XYII", Angle::from_angle(2.)),
            ("XZII", Angle::from_angle(3.)),
            ("IIXX", Angle::from_angle(4.)),
            ("IIYX", Angle::from_angle(5.)),
            ("IIZX", Angle::from_angle(6.)),
        ]);
        pp_x.cx(0, 1);
        pp_x.cx(2, 3);
        let pp_ref_x = PauliPolynomial::from_hamiltonian(vec![
            ("XIII", Angle::from_angle(1.)),
            ("YZII", Angle::from_angle(2.)),
            ("YYII", Angle::from_angle(-3.)),
            ("IIXI", Angle::from_angle(4.)),
            ("IIYI", Angle::from_angle(5.)),
            ("IIZX", Angle::from_angle(6.)),
        ]);
        assert_eq!(pp_x, pp_ref_x);
        // Y in control or target
        let mut pp_y = PauliPolynomial::from_hamiltonian(vec![
            ("YXII", Angle::from_angle(1.)),
            ("YYII", Angle::from_angle(2.)),
            ("YZII", Angle::from_angle(3.)),
            ("IIXY", Angle::from_angle(4.)),
            ("IIYY", Angle::from_angle(5.)),
            ("IIZY", Angle::from_angle(6.)),
        ]);
        pp_y.cx(0, 1);
        pp_y.cx(2, 3);
        let pp_ref_y = PauliPolynomial::from_hamiltonian(vec![
            ("YIII", Angle::from_angle(1.)),
            ("XZII", Angle::from_angle(-2.)),
            ("XYII", Angle::from_angle(3.)),
            ("IIYZ", Angle::from_angle(4.)),
            ("IIXZ", Angle::from_angle(-5.)),
            ("IIIY", Angle::from_angle(6.)),
        ]);
        assert_eq!(pp_y, pp_ref_y);
        // Z in control or target
        let mut pp_z = PauliPolynomial::from_hamiltonian(vec![
            ("ZXII", Angle::from_angle(1.)),
            ("ZYII", Angle::from_angle(2.)),
            ("ZZII", Angle::from_angle(3.)),
            ("IIXZ", Angle::from_angle(4.)),
            ("IIYZ", Angle::from_angle(5.)),
            ("IIZZ", Angle::from_angle(6.)),
        ]);
        pp_z.cx(0, 1);
        pp_z.cx(2, 3);
        let pp_ref_z = PauliPolynomial::from_hamiltonian(vec![
            ("ZXII", Angle::from_angle(1.)),
            ("IYII", Angle::from_angle(2.)),
            ("IZII", Angle::from_angle(3.)),
            ("IIYY", Angle::from_angle(-4.)),
            ("IIXY", Angle::from_angle(5.)),
            ("IIIZ", Angle::from_angle(6.)),
        ]);
        assert_eq!(pp_z, pp_ref_z);
    }
}
