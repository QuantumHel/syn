use bitvec::prelude::BitVec;

#[derive(Clone, PartialEq, Eq, Debug)]
pub(super) struct PauliString {
    pub(super) x: BitVec,
    pub(super) z: BitVec,
}

impl PauliString {
    /// Constructor for PauliString
    pub fn new(pauli_x: BitVec, pauli_z: BitVec) -> Self {
        assert!(pauli_x.len() == pauli_z.len());
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
    pub fn from_text_string(mut pauli: String) -> Self {
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

    pub(super) fn s(&mut self) {
        self.x ^= &self.z;
    }

    pub(super) fn v(&mut self) {
        self.z ^= &self.x;
    }

    pub(super) fn y_bitmask(&self) -> BitVec {
        let mut mask = self.x.clone();
        mask &= &self.z;
        mask
    }
}

pub(super) fn cx(control: &mut PauliString, target: &mut PauliString) {
    assert!(control.len() == target.len());
    target.x ^= &control.x;
    control.z ^= &target.z;
}

#[cfg(test)]
mod tests {
    use super::*;
    use bitvec::bitvec;
    use bitvec::prelude::Lsb0;

    #[test]
    fn test_from_basis_int() {
        let i = 3;
        let length = 5;

        let paulivec = PauliString::from_basis_int(i, length);
        assert!(paulivec.x.get(i).unwrap());
        assert!(paulivec.z.get(i + length).unwrap());
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
        let pauli_string = String::from("IXYZ");
        let paulivec = PauliString::from_text_string(pauli_string);
        let x_ref = bitvec![0, 1, 1, 0];
        let z_ref = bitvec![0, 0, 1, 1];
        assert!(paulivec.x == x_ref);
        assert!(paulivec.z == z_ref);
    }

    #[test]
    fn test_s() {
        let mut paulivec = PauliString::from_text_string(String::from("IXYZ"));
        paulivec.s();
        let paulivec_ref = PauliString::from_text_string(String::from("IXZY"));
        assert!(paulivec.x == paulivec_ref.x);
        assert!(paulivec.z == paulivec_ref.z);
    }

    #[test]
    fn test_v() {
        let mut paulivec = PauliString::from_text_string(String::from("IXYZ"));
        paulivec.v();
        let paulivec_ref = PauliString::from_text_string(String::from("IYXZ"));
        assert!(paulivec.x == paulivec_ref.x);
        assert!(paulivec.z == paulivec_ref.z);
    }

    #[test]
    fn test_y_bitmask() {
        let paulivec = PauliString::from_text_string(String::from("IYXYZY"));
        let y_bitmask = paulivec.y_bitmask();
        let y_bitmask_ref = bitvec![0, 1, 0, 1, 0, 1];
        assert_eq!(y_bitmask, y_bitmask_ref);
    }
}
