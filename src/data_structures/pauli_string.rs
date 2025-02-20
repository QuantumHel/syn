use bitvec::{prelude::BitVec, slice::BitSlice};
use std::cell::RefCell;
use std::fmt;
use std::iter::zip;

use super::PauliLetter;

#[derive(Clone, PartialEq, Eq, Debug)]
pub struct PauliString {
    pub(super) x: RefCell<BitVec>,
    pub(super) z: RefCell<BitVec>,
}

impl PauliString {
    /// Constructor for PauliString
    pub fn new(pauli_x: BitVec, pauli_z: BitVec) -> Self {
        assert!(pauli_x.len() == pauli_z.len());
        PauliString {
            x: RefCell::new(pauli_x),
            z: RefCell::new(pauli_z),
        }
    }

    /// Construct identity PauliString for position `i`
    pub fn from_basis_int(i: usize, length: usize) -> Self {
        assert!(length > i);
        let mut x = BitVec::repeat(false, 2 * length);
        let mut z = BitVec::repeat(false, 2 * length);
        x.set(i, true);
        z.set(i + length, true);
        PauliString::new(x, z)
    }

    /// Takes in a String containing "I"
    pub fn from_text(pauli: &str) -> Self {
        let (x, z): (BitVec, BitVec) = pauli
            .chars()
            .map(|pauli_char| {
                let (x, z) = match pauli_char {
                    'I' => (false, false),
                    'X' => (true, false),
                    'Y' => (true, true),
                    'Z' => (false, true),
                    _ => panic!("Letter not recognized"),
                };
                (x, z)
            })
            .collect();

        PauliString::new(x, z)
    }

    pub fn x(&self, i: usize) -> bool {
        self.x.borrow()[i]
    }

    pub fn z(&self, i: usize) -> bool {
        self.z.borrow()[i]
    }

    pub fn pauli(&self, i: usize) -> PauliLetter {
        PauliLetter::new(self.x.borrow()[i], self.z.borrow()[i])
    }

    pub fn len(&self) -> usize {
        self.x.borrow().len()
    }

    pub fn is_empty(&self) -> bool {
        self.x.borrow().is_empty()
    }

    pub(crate) fn s(&self) {
        *self.z.borrow_mut() ^= self.x.borrow().as_bitslice();
    }

    pub(crate) fn masked_s(&self, mask: &BitSlice) {
        let mut mask = mask.to_owned();
        mask &= self.x.borrow().clone();
        *self.z.borrow_mut() ^= &mask;
    }

    pub(crate) fn v(&self) {
        *self.x.borrow_mut() ^= self.z.borrow().as_bitslice();
    }

    pub(crate) fn masked_v(&self, mask: &BitSlice) {
        let mut mask = mask.to_owned();
        mask &= self.z.borrow().as_bitslice();
        *self.x.borrow_mut() ^= &mask;
    }

    pub(crate) fn h(&self) {
        let tmp = self.z.borrow().clone();
        *self.z.borrow_mut() = self.x.borrow().clone();
        *self.x.borrow_mut() = tmp;
    }

    pub(crate) fn masked_h(&self, mask: &BitSlice) {
        let mut mask = mask.to_owned();
        *self.x.borrow_mut() ^= self.z.borrow().as_bitslice();
        mask &= self.x.borrow().as_bitslice();
        *self.z.borrow_mut() ^= &mask;
        *self.x.borrow_mut() ^= self.z.borrow().as_bitslice();
    }

    pub(crate) fn y_bitmask(&self) -> BitVec {
        let mut mask = self.x.borrow().clone();
        mask &= self.z.borrow().as_bitslice();
        mask
    }

    pub(crate) fn masked_y_bitmask(&self, mask: &BitSlice) -> BitVec {
        let mut mask = mask.to_owned();
        mask &= self.x.borrow().as_bitslice();
        mask &= self.z.borrow().as_bitslice();
        mask
    }
}

pub(crate) fn cx(control: &PauliString, target: &PauliString) {
    assert!(control.len() == target.len());
    *target.x.borrow_mut() ^= control.x.borrow().as_bitslice();
    *control.z.borrow_mut() ^= target.z.borrow().as_bitslice();
}

pub(crate) fn masked_cx(control: &PauliString, target: &PauliString, mask: &BitSlice) {
    assert!(control.len() == target.len());
    let mut x_mask = mask.to_owned();
    let mut z_mask = mask.to_owned();
    x_mask &= control.x.borrow().as_bitslice();
    z_mask &= target.z.borrow().as_bitslice();
    *target.x.borrow_mut() ^= &x_mask;
    *control.z.borrow_mut() ^= &z_mask;
}

impl fmt::Display for PauliString {
    // This trait requires `fmt` with this exact signature.
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut pauli_str = String::new();
        for (x, z) in zip(self.x.borrow().iter(), self.z.borrow().iter()) {
            match (*x, *z) {
                (false, false) => pauli_str.push('I'),
                (false, true) => pauli_str.push('Z'),
                (true, false) => pauli_str.push('X'),
                (true, true) => pauli_str.push('Y'),
            }
            pauli_str.push(' ');
        }
        pauli_str.pop();
        // pauli_str.push_str("]");
        write!(f, "{}", pauli_str)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bitvec::prelude::Lsb0;
    use bitvec::{bits, bitvec};

    #[test]
    fn test_from_basis_int() {
        let i = 3;
        let length = 5;

        let paulivec = PauliString::from_basis_int(i, length);
        assert!(paulivec.x.borrow().get(i).unwrap());
        assert!(paulivec.z.borrow().get(i + length).unwrap());
    }

    #[test]
    #[should_panic]
    fn test_from_basis_int_oversized_i() {
        let i = 5;
        let length = 3;
        PauliString::from_basis_int(i, length);
    }

    #[test]
    fn test_from_text_string() {
        let pauli_string = "IXYZ";
        let paulivec = PauliString::from_text(pauli_string);
        let x_ref = bitvec![0, 1, 1, 0];
        let z_ref = bitvec![0, 0, 1, 1];
        assert_eq!(paulivec.x.borrow().as_bitslice(), &x_ref);
        assert_eq!(paulivec.z.borrow().as_bitslice(), &z_ref);
    }

    #[test]
    fn test_xz_access() {
        let pauli_string = "IXYZ";
        let paulivec = PauliString::from_text(pauli_string);

        assert!(paulivec.x(1));
        assert!(!paulivec.z(1));

        assert!(!paulivec.x(3));
        assert!(paulivec.z(3));
    }

    #[test]
    fn test_pauli_string_s() {
        let paulivec = PauliString::from_text("IXYZ");
        paulivec.s();
        let paulivec_ref = PauliString::from_text("IYXZ");
        assert!(paulivec.x == paulivec_ref.x);
        assert!(paulivec.z == paulivec_ref.z);
    }

    #[test]
    fn test_pauli_string_v() {
        let paulivec = PauliString::from_text("IXYZ");
        paulivec.v();
        let paulivec_ref = PauliString::from_text("IXZY");
        assert!(paulivec.x == paulivec_ref.x);
        assert!(paulivec.z == paulivec_ref.z);
    }

    #[test]
    fn test_pauli_string_h() {
        let paulivec = PauliString::from_text("IXYZ");
        paulivec.h();
        let paulivec_ref = PauliString::from_text("IZYX");
        assert!(paulivec.x == paulivec_ref.x);
        assert!(paulivec.z == paulivec_ref.z);
    }

    #[test]
    fn test_pauli_string_masked_s() {
        let paulivec = PauliString::from_text("IXYZIXYZ");
        let mask = bits![usize, Lsb0; 0, 0, 0, 0, 1, 1, 1, 1];
        paulivec.masked_s(mask);
        let paulivec_ref = PauliString::from_text("IXYZIYXZ");

        assert!(paulivec.x == paulivec_ref.x);
        assert!(paulivec.z == paulivec_ref.z);
    }

    #[test]
    fn test_pauli_string_masked_v() {
        let paulivec = PauliString::from_text("IXYZIXYZ");
        let mask = bits![usize, Lsb0; 0, 0, 0, 0, 1, 1, 1, 1];
        paulivec.masked_v(mask);
        let paulivec_ref = PauliString::from_text("IXYZIXZY");
        assert!(paulivec.x == paulivec_ref.x);
        assert!(paulivec.z == paulivec_ref.z);
    }

    #[test]
    fn test_pauli_string_masked_h() {
        let paulivec = PauliString::from_text("IXYZIXYZ");
        let mask = bits![usize, Lsb0; 0, 0, 0, 0, 1, 1, 1, 1];
        paulivec.masked_h(mask);
        let paulivec_ref = PauliString::from_text("IXYZIZYX");
        assert!(paulivec.x == paulivec_ref.x);
        assert!(paulivec.z == paulivec_ref.z);
    }

    #[test]
    fn test_pauli_string_cx() {
        let control = PauliString::from_text("IIIIXXXXYYYYZZZZ");
        let target = PauliString::from_text("IXYZIXYZIXYZIXYZ");
        cx(&control, &target);
        let control_ref = PauliString::from_text("IIZZXXYYYYXXZZII");
        let target_ref = PauliString::from_text("IXYZXIZYXIZYIXYZ");

        assert!(control.x == control_ref.x);
        assert!(control.z == control_ref.z);

        assert!(target.x == target_ref.x);
        assert!(target.z == target_ref.z);
    }

    #[test]
    fn test_pauli_string_masked_cx() {
        let control = PauliString::from_text("IIIIXXXXYYYYZZZZIIIIXXXXYYYYZZZZ");
        let target = PauliString::from_text("IXYZIXYZIXYZIXYZIXYZIXYZIXYZIXYZ");
        let mask = bits![usize, Lsb0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1];
        masked_cx(&control, &target, mask);
        let control_ref = PauliString::from_text("IIIIXXXXYYYYZZZZIIZZXXYYYYXXZZII");
        let target_ref = PauliString::from_text("IXYZIXYZIXYZIXYZIXYZXIZYXIZYIXYZ");

        assert!(control.x == control_ref.x);
        assert!(control.z == control_ref.z);

        assert!(target.x == target_ref.x);
        assert!(target.z == target_ref.z);
    }

    #[test]
    fn test_y_bitmask() {
        let paulivec = PauliString::from_text("IYXYZY");
        let y_bitmask = paulivec.y_bitmask();
        let y_bitmask_ref = bitvec![0, 1, 0, 1, 0, 1];
        assert_eq!(y_bitmask, y_bitmask_ref);
    }

    #[test]
    fn test_pauli_string_display() {
        let pauli_string = PauliString::from_text("IXYZI");
        assert_eq!(String::from("I X Y Z I"), pauli_string.to_string());
    }
}
