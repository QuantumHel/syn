use bitvec::{prelude::BitVec, slice::BitSlice};
use std::fmt;
use std::iter::zip;
use std::sync::RwLock;

use super::PauliLetter;

/// A vector of Pauli operators (`I`, `X`, `Y`, `Z`)-
#[derive(Debug)]
pub struct PauliVec {
    /// The `x` bits. There is one bit for every operator in this vec. When it is
    /// set, the corresponding operator is either `X` or `Y` (depending on the `z`
    /// vec).
    pub(super) x: RwLock<BitVec>,
    /// The `z` bits. There is one bit for every operator in this vec. When it is
    /// set, the corresponding operator is either `Z` or `Y` (depending on the `x`
    /// vec).
    pub(super) z: RwLock<BitVec>,
}

impl PartialEq for PauliVec {
    fn eq(&self, other: &Self) -> bool {
        *self.x.read().unwrap() == *other.x.read().unwrap()
            && *self.z.read().unwrap() == *other.z.read().unwrap()
    }
}

impl Eq for PauliVec {}

impl Clone for PauliVec {
    fn clone(&self) -> Self {
        PauliVec {
            x: RwLock::new(self.x.read().unwrap().clone()),
            z: RwLock::new(self.z.read().unwrap().clone()),
        }
    }
}

impl PauliVec {
    /// Constructs a new Pauli vector from separate `pauli_x` and `pauli_z` vectors.
    /// Both must have same length. The letters are then represented in the following
    /// way:
    ///
    ///  ```text
    /// | x     | z     | letter |
    /// |-------|-------|--------|
    /// | false | false | I      |
    /// | true  | false | X      |
    /// | true  | true  | Y      |
    /// | false | true  | Z      |
    /// ```
    ///
    /// # Panics
    /// Panics if `pauli_x` and `pauli_z` are not of the same length.
    pub fn new(pauli_x: BitVec, pauli_z: BitVec) -> Self {
        assert_eq!(pauli_x.len(), pauli_z.len());
        PauliVec {
            x: RwLock::new(pauli_x),
            z: RwLock::new(pauli_z),
        }
    }

    /// Constructs a new Pauli vector from a string of Pauli letters.
    ///
    /// The letters must be upper case (to avoid confusion with the complex `i`).
    /// Spaces are ignored. The valid letters are `I`, `X`, `Y` and `Z`.
    ///
    /// # Panics
    /// Panics if an unknown letter is encountered.
    ///
    /// # Examples
    /// ```
    /// # use syn::data_structures::PauliVec;
    /// let p1 = PauliVec::from_text("XZZ");
    /// let p2 = PauliVec::from_text("X Z   Z");
    /// assert_eq!(p1, p2);
    /// ```
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

        PauliVec::new(x, z)
    }

    /// Returns whether the `i`th operator in the Pauli vector is `X` or `Y`.
    ///
    /// # Panics
    /// Panics if `i` is out of bounds.
    pub fn x(&self, i: usize) -> bool {
        self.x.read().unwrap()[i]
    }

    /// Returns whether the `i`th operator in the Pauli vector is `Z` or `Y`.
    ///
    /// # Panics
    /// Panics if `i` is out of bounds.
    pub fn z(&self, i: usize) -> bool {
        self.z.read().unwrap()[i]
    }

    /// Returns the pauli letter at the `i`th position of the Pauli vector.
    ///
    /// # Panics
    /// Panics if `i` is out of bounds.
    pub fn pauli(&self, i: usize) -> PauliLetter {
        PauliLetter::new(self.x(i), self.z(i))
    }

    /// Returns the length of the Pauli vector.
    pub fn len(&self) -> usize {
        self.x.read().unwrap().len()
    }

    /// Returns whether the Pauli vector is empty.
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

pub(crate) fn cx(control: &PauliVec, target: &PauliVec) {
    assert_eq!(control.len(), target.len());
    *target.x.write().unwrap() ^= control.x.read().unwrap().as_bitslice();
    *control.z.write().unwrap() ^= target.z.read().unwrap().as_bitslice();
}

pub(crate) fn masked_cx(control: &PauliVec, target: &PauliVec, mask: &BitSlice) {
    assert_eq!(control.len(), target.len());
    let mut x_mask = mask.to_owned();
    let mut z_mask = mask.to_owned();
    x_mask &= control.x.read().unwrap().as_bitslice();
    z_mask &= target.z.read().unwrap().as_bitslice();
    *target.x.write().unwrap() ^= &x_mask;
    *control.z.write().unwrap() ^= &z_mask;
}

impl fmt::Display for PauliVec {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut out = String::new();
        for (x, z) in zip(self.x.read().unwrap().iter(), self.z.read().unwrap().iter()) {
            match (*x, *z) {
                (false, false) => out.push('I'),
                (false, true) => out.push('Z'),
                (true, false) => out.push('X'),
                (true, true) => out.push('Y'),
            }
            out.push(' ');
        }
        out.pop();
        write!(f, "{}", out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bitvec::prelude::Lsb0;
    use bitvec::{bits, bitvec};

    #[test]
    fn from_text_string() {
        let pauli_string = "IXYZ";
        let paulivec = PauliVec::from_text(pauli_string);
        let x_ref = bitvec![0, 1, 1, 0];
        let z_ref = bitvec![0, 0, 1, 1];
        assert_eq!(paulivec.x.read().unwrap().as_bitslice(), &x_ref);
        assert_eq!(paulivec.z.read().unwrap().as_bitslice(), &z_ref);
    }

    #[test]
    fn xz_access() {
        let pauli_string = "IXYZ";
        let paulivec = PauliVec::from_text(pauli_string);

        assert!(paulivec.x(1));
        assert!(!paulivec.z(1));

        assert!(!paulivec.x(3));
        assert!(paulivec.z(3));
    }

    #[test]
    fn apply_s() {
        let paulivec = PauliVec::from_text("IXYZ");
        paulivec.s();
        let paulivec_ref = PauliVec::from_text("IYXZ");

        assert_eq!(paulivec, paulivec_ref);
    }

    #[test]
    fn apply_v() {
        let paulivec = PauliVec::from_text("IXYZ");
        paulivec.v();
        let paulivec_ref = PauliVec::from_text("IXZY");

        assert_eq!(paulivec, paulivec_ref);
    }

    #[test]
    fn apply_h() {
        let paulivec = PauliVec::from_text("IXYZ");
        paulivec.h();
        let paulivec_ref = PauliVec::from_text("IZYX");

        assert_eq!(paulivec, paulivec_ref);
    }

    #[test]
    fn apply_masked_s() {
        let paulivec = PauliVec::from_text("IXYZIXYZ");
        let mask = bits![usize, Lsb0; 0, 0, 0, 0, 1, 1, 1, 1];
        paulivec.masked_s(mask);
        let paulivec_ref = PauliVec::from_text("IXYZIYXZ");

        assert_eq!(paulivec, paulivec_ref);
    }

    #[test]
    fn apply_masked_v() {
        let paulivec = PauliVec::from_text("IXYZIXYZ");
        let mask = bits![usize, Lsb0; 0, 0, 0, 0, 1, 1, 1, 1];
        paulivec.masked_v(mask);
        let paulivec_ref = PauliVec::from_text("IXYZIXZY");

        assert_eq!(paulivec, paulivec_ref);
    }

    #[test]
    fn apply_masked_h() {
        let paulivec = PauliVec::from_text("IXYZIXYZ");
        let mask = bits![usize, Lsb0; 0, 0, 0, 0, 1, 1, 1, 1];
        paulivec.masked_h(mask);
        let paulivec_ref = PauliVec::from_text("IXYZIZYX");

        assert_eq!(paulivec, paulivec_ref);
    }

    #[test]
    fn apply_cx() {
        let control = PauliVec::from_text("IIIIXXXXYYYYZZZZ");
        let target = PauliVec::from_text("IXYZIXYZIXYZIXYZ");
        cx(&control, &target);
        let control_ref = PauliVec::from_text("IIZZXXYYYYXXZZII");
        let target_ref = PauliVec::from_text("IXYZXIZYXIZYIXYZ");

        assert_eq!(control, control_ref);
        assert_eq!(target, target_ref);
    }

    #[test]
    fn apply_masked_cx() {
        let control = PauliVec::from_text("IIIIXXXXYYYYZZZZIIIIXXXXYYYYZZZZ");
        let target = PauliVec::from_text("IXYZIXYZIXYZIXYZIXYZIXYZIXYZIXYZ");
        let mask = bits![usize, Lsb0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1];
        masked_cx(&control, &target, mask);
        let control_ref = PauliVec::from_text("IIIIXXXXYYYYZZZZIIZZXXYYYYXXZZII");
        let target_ref = PauliVec::from_text("IXYZIXYZIXYZIXYZIXYZXIZYXIZYIXYZ");

        assert_eq!(control, control_ref);
        assert_eq!(target, target_ref);
    }

    #[test]
    fn y_bitmask() {
        let paulivec = PauliVec::from_text("IYXYZY");
        let y_bitmask = paulivec.y_bitmask();
        let y_bitmask_ref = bitvec![0, 1, 0, 1, 0, 1];
        assert_eq!(y_bitmask, y_bitmask_ref);
    }

    #[test]
    fn string_display() {
        let paulivec = PauliVec::from_text("IXYZI");
        assert_eq!(paulivec.to_string(), String::from("I X Y Z I"));
    }
}
