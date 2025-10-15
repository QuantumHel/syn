use bitvec::{prelude::BitVec, slice::BitSlice};
use std::fmt;
use std::iter::zip;

use super::PauliLetter;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PauliString {
    pub(super) x: BitVec,
    pub(super) z: BitVec,
}

impl PauliString {
    /// Constructor for PauliString
    pub fn new(pauli_x: BitVec, pauli_z: BitVec) -> Self {
        assert_eq!(pauli_x.len(), pauli_z.len());
        PauliString {
            x: pauli_x,
            z: pauli_z,
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
        self.x[i]
    }

    pub fn x_weight(&self) -> usize {
        self.x.count_ones()
    }

    pub fn z_weight(&self) -> usize {
        self.z.count_ones()
    }

    pub fn z(&self, i: usize) -> bool {
        self.z[i]
    }

    pub fn pauli(&self, i: usize) -> PauliLetter {
        PauliLetter::new(self.x(i), self.z(i))
    }

    pub fn len(&self) -> usize {
        self.x.len()
    }

    pub fn is_empty(&self) -> bool {
        self.x.is_empty()
    }

    pub(crate) fn s(&mut self) {
        self.z ^= &self.x;
    }

    pub(crate) fn masked_s(&mut self, mask: &BitSlice) {
        let mut mask = mask.to_owned();
        mask &= &self.x;
        self.z ^= &mask;
    }

    pub(crate) fn v(&mut self) {
        self.x ^= &self.z;
    }

    pub(crate) fn masked_v(&mut self, mask: &BitSlice) {
        let mut mask = mask.to_owned();
        mask &= &self.z;
        self.x ^= &mask;
    }

    #[allow(dead_code)]
    pub(crate) fn h(&mut self) {
        std::mem::swap(&mut self.x, &mut self.z);
    }

    #[allow(dead_code)]
    pub(crate) fn masked_h(&mut self, mask: &BitSlice) {
        let mut mask = mask.to_owned();
        self.x ^= &self.z;
        mask &= &self.x;
        self.z ^= &mask;
        self.x ^= &self.z;
    }

    pub(crate) fn y_bitmask(&self) -> BitVec {
        let mut mask = self.x.clone();
        mask &= &self.z;
        mask
    }

    pub(crate) fn masked_y_bitmask(&self, mask: &BitSlice) -> BitVec {
        let mut mask = mask.to_owned();
        mask &= &self.x;
        mask &= &self.z;
        mask
    }
}

pub(crate) fn cx(control: &mut PauliString, target: &mut PauliString) {
    assert_eq!(control.len(), target.len());
    target.x ^= &control.x;
    control.z ^= &target.z;
}

pub(crate) fn masked_cx(control: &mut PauliString, target: &mut PauliString, mask: &BitSlice) {
    assert_eq!(control.len(), target.len());
    let mut x_mask = mask.to_owned();
    let mut z_mask = mask.to_owned();
    x_mask &= &control.x;
    z_mask &= &target.z;
    target.x ^= &x_mask;
    control.z ^= &z_mask;
}

impl fmt::Display for PauliString {
    // This trait requires `fmt` with this exact signature.
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut pauli_str = String::new();
        for (x, z) in zip(&self.x, &self.z) {
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
        assert!(paulivec.x(i));
        assert!(paulivec.z(i + length));
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
        assert_eq!(paulivec.x, x_ref);
        assert_eq!(paulivec.z, z_ref);
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
        let mut paulivec = PauliString::from_text("IXYZ");
        paulivec.s();
        let paulivec_ref = PauliString::from_text("IYXZ");

        assert_eq!(paulivec, paulivec_ref);
    }

    #[test]
    fn test_pauli_string_v() {
        let mut paulivec = PauliString::from_text("IXYZ");
        paulivec.v();
        let paulivec_ref = PauliString::from_text("IXZY");

        assert_eq!(paulivec, paulivec_ref);
    }

    #[test]
    fn test_pauli_string_h() {
        let mut paulivec = PauliString::from_text("IXYZ");
        paulivec.h();
        let paulivec_ref = PauliString::from_text("IZYX");

        assert_eq!(paulivec, paulivec_ref);
    }

    #[test]
    fn test_pauli_string_masked_s() {
        let mut paulivec = PauliString::from_text("IXYZIXYZ");
        let mask = bits![usize, Lsb0; 0, 0, 0, 0, 1, 1, 1, 1];
        paulivec.masked_s(mask);
        let paulivec_ref = PauliString::from_text("IXYZIYXZ");

        assert_eq!(paulivec, paulivec_ref);
    }

    #[test]
    fn test_pauli_string_masked_v() {
        let mut paulivec = PauliString::from_text("IXYZIXYZ");
        let mask = bits![usize, Lsb0; 0, 0, 0, 0, 1, 1, 1, 1];
        paulivec.masked_v(mask);
        let paulivec_ref = PauliString::from_text("IXYZIXZY");

        assert_eq!(paulivec, paulivec_ref);
    }

    #[test]
    fn test_pauli_string_masked_h() {
        let mut paulivec = PauliString::from_text("IXYZIXYZ");
        let mask = bits![usize, Lsb0; 0, 0, 0, 0, 1, 1, 1, 1];
        paulivec.masked_h(mask);
        let paulivec_ref = PauliString::from_text("IXYZIZYX");

        assert_eq!(paulivec, paulivec_ref);
    }

    #[test]
    fn test_pauli_string_cx() {
        let mut control = PauliString::from_text("IIIIXXXXYYYYZZZZ");
        let mut target = PauliString::from_text("IXYZIXYZIXYZIXYZ");
        cx(&mut control, &mut target);
        let control_ref = PauliString::from_text("IIZZXXYYYYXXZZII");
        let target_ref = PauliString::from_text("IXYZXIZYXIZYIXYZ");

        assert_eq!(control, control_ref);
        assert_eq!(target, target_ref);
    }

    #[test]
    fn test_pauli_string_masked_cx() {
        let mut control = PauliString::from_text("IIIIXXXXYYYYZZZZIIIIXXXXYYYYZZZZ");
        let mut target = PauliString::from_text("IXYZIXYZIXYZIXYZIXYZIXYZIXYZIXYZ");
        let mask = bits![usize, Lsb0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1];
        masked_cx(&mut control, &mut target, mask);
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
