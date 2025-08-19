use bitvec::{prelude::BitVec, slice::BitSlice};
use std::fmt;
use std::iter::zip;
use std::sync::RwLock;

use super::PauliLetter;

#[derive(Debug)]
pub struct PauliString {
    pub(super) x: RwLock<BitVec>,
    pub(super) z: RwLock<BitVec>,
}

impl PartialEq for PauliString {
    fn eq(&self, other: &Self) -> bool {
        *self.x.read().unwrap() == *other.x.read().unwrap()
            && *self.z.read().unwrap() == *other.z.read().unwrap()
    }
}

impl Eq for PauliString {}

impl Clone for PauliString {
    fn clone(&self) -> Self {
        PauliString {
            x: RwLock::new(self.x.read().unwrap().clone()),
            z: RwLock::new(self.z.read().unwrap().clone()),
        }
    }
}

impl PauliString {
    /// Constructor for PauliString
    pub fn new(pauli_x: BitVec, pauli_z: BitVec) -> Self {
        assert_eq!(pauli_x.len(), pauli_z.len());
        PauliString {
            x: RwLock::new(pauli_x),
            z: RwLock::new(pauli_z),
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
        let (x, z) = pauli
            .chars()
            .filter_map(|pauli_char| match pauli_char {
                'I' => Some((false, false)),
                'X' => Some((true, false)),
                'Y' => Some((true, true)),
                'Z' => Some((false, true)),
                ' ' => None,
                _ => panic!("Letter not recognized"),
            })
            .collect();

        PauliString::new(x, z)
    }

    pub fn x(&self, i: usize) -> bool {
        self.x.read().unwrap()[i]
    }

    pub fn z(&self, i: usize) -> bool {
        self.z.read().unwrap()[i]
    }

    pub fn pauli(&self, i: usize) -> PauliLetter {
        PauliLetter::new(self.x(i), self.z(i))
    }

    pub fn len(&self) -> usize {
        self.x.read().unwrap().len()
    }

    pub fn is_empty(&self) -> bool {
        self.x.read().unwrap().is_empty()
    }

    pub(crate) fn s(&self) {
        *self.z.write().unwrap() ^= self.x.read().unwrap().as_bitslice();
    }

    pub(crate) fn masked_s(&self, mask: &BitSlice) {
        let mut mask = mask.to_owned();
        mask &= self.x.read().unwrap().as_bitslice();
        *self.z.write().unwrap() ^= &mask;
    }

    pub(crate) fn v(&self) {
        *self.x.write().unwrap() ^= self.z.read().unwrap().as_bitslice();
    }

    pub(crate) fn masked_v(&self, mask: &BitSlice) {
        let mut mask = mask.to_owned();
        mask &= self.z.read().unwrap().as_bitslice();
        *self.x.write().unwrap() ^= &mask;
    }

    #[allow(dead_code)]
    pub(crate) fn h(&self) {
        std::mem::swap(&mut *self.x.write().unwrap(), &mut *self.z.write().unwrap());
    }

    #[allow(dead_code)]
    pub(crate) fn masked_h(&self, mask: &BitSlice) {
        let mut mask = mask.to_owned();
        *self.x.write().unwrap() ^= self.z.read().unwrap().as_bitslice();
        mask &= self.x.read().unwrap().as_bitslice();
        *self.z.write().unwrap() ^= &mask;
        *self.x.write().unwrap() ^= self.z.read().unwrap().as_bitslice();
    }

    pub(crate) fn y_bitmask(&self) -> BitVec {
        let mut mask = self.x.read().unwrap().clone();
        mask &= self.z.read().unwrap().as_bitslice();
        mask
    }

    pub(crate) fn masked_y_bitmask(&self, mask: &BitSlice) -> BitVec {
        let mut mask = mask.to_owned();
        mask &= self.x.read().unwrap().as_bitslice();
        mask &= self.z.read().unwrap().as_bitslice();
        mask
    }
}

pub(crate) fn cx(control: &PauliString, target: &PauliString) {
    assert_eq!(control.len(), target.len());
    *target.x.write().unwrap() ^= control.x.read().unwrap().as_bitslice();
    *control.z.write().unwrap() ^= target.z.read().unwrap().as_bitslice();
}

pub(crate) fn masked_cx(control: &PauliString, target: &PauliString, mask: &BitSlice) {
    assert_eq!(control.len(), target.len());
    let mut x_mask = mask.to_owned();
    let mut z_mask = mask.to_owned();
    x_mask &= control.x.read().unwrap().as_bitslice();
    z_mask &= target.z.read().unwrap().as_bitslice();
    *target.x.write().unwrap() ^= &x_mask;
    *control.z.write().unwrap() ^= &z_mask;
}

impl fmt::Display for PauliString {
    // This trait requires `fmt` with this exact signature.
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut pauli_str = String::new();
        for (x, z) in zip(self.x.read().unwrap().iter(), self.z.read().unwrap().iter()) {
            match (*x, *z) {
                (false, false) => pauli_str.push('I'),
                (false, true) => pauli_str.push('Z'),
                (true, false) => pauli_str.push('X'),
                (true, true) => pauli_str.push('Y'),
            }
            pauli_str.push(' ');
        }
        pauli_str.pop();
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
        assert!(paulivec.x.read().unwrap().get(i).unwrap());
        assert!(paulivec.z.read().unwrap().get(i + length).unwrap());
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
        assert_eq!(paulivec.x.read().unwrap().as_bitslice(), &x_ref);
        assert_eq!(paulivec.z.read().unwrap().as_bitslice(), &z_ref);
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

        assert_eq!(paulivec, paulivec_ref);
    }

    #[test]
    fn test_pauli_string_v() {
        let paulivec = PauliString::from_text("IXYZ");
        paulivec.v();
        let paulivec_ref = PauliString::from_text("IXZY");

        assert_eq!(paulivec, paulivec_ref);
    }

    #[test]
    fn test_pauli_string_h() {
        let paulivec = PauliString::from_text("IXYZ");
        paulivec.h();
        let paulivec_ref = PauliString::from_text("IZYX");

        assert_eq!(paulivec, paulivec_ref);
    }

    #[test]
    fn test_pauli_string_masked_s() {
        let paulivec = PauliString::from_text("IXYZIXYZ");
        let mask = bits![usize, Lsb0; 0, 0, 0, 0, 1, 1, 1, 1];
        paulivec.masked_s(mask);
        let paulivec_ref = PauliString::from_text("IXYZIYXZ");

        assert_eq!(paulivec, paulivec_ref);
    }

    #[test]
    fn test_pauli_string_masked_v() {
        let paulivec = PauliString::from_text("IXYZIXYZ");
        let mask = bits![usize, Lsb0; 0, 0, 0, 0, 1, 1, 1, 1];
        paulivec.masked_v(mask);
        let paulivec_ref = PauliString::from_text("IXYZIXZY");

        assert_eq!(paulivec, paulivec_ref);
    }

    #[test]
    fn test_pauli_string_masked_h() {
        let paulivec = PauliString::from_text("IXYZIXYZ");
        let mask = bits![usize, Lsb0; 0, 0, 0, 0, 1, 1, 1, 1];
        paulivec.masked_h(mask);
        let paulivec_ref = PauliString::from_text("IXYZIZYX");

        assert_eq!(paulivec, paulivec_ref);
    }

    #[test]
    fn test_pauli_string_cx() {
        let control = PauliString::from_text("IIIIXXXXYYYYZZZZ");
        let target = PauliString::from_text("IXYZIXYZIXYZIXYZ");
        cx(&control, &target);
        let control_ref = PauliString::from_text("IIZZXXYYYYXXZZII");
        let target_ref = PauliString::from_text("IXYZXIZYXIZYIXYZ");

        assert_eq!(control, control_ref);
        assert_eq!(target, target_ref);
    }

    #[test]
    fn test_pauli_string_masked_cx() {
        let control = PauliString::from_text("IIIIXXXXYYYYZZZZIIIIXXXXYYYYZZZZ");
        let target = PauliString::from_text("IXYZIXYZIXYZIXYZIXYZIXYZIXYZIXYZ");
        let mask = bits![usize, Lsb0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1];
        masked_cx(&control, &target, mask);
        let control_ref = PauliString::from_text("IIIIXXXXYYYYZZZZIIZZXXYYYYXXZZII");
        let target_ref = PauliString::from_text("IXYZIXYZIXYZIXYZIXYZXIZYXIZYIXYZ");

        assert_eq!(control, control_ref);
        assert_eq!(target, target_ref);
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
        assert_eq!(pauli_string.to_string(), String::from("I X Y Z I"));
    }
}
