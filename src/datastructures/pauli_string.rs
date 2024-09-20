use bitvec::prelude::BitVec;

#[derive(Clone)]
pub(super) struct PauliString{
    pub(super) x: BitVec,
    pub(super) z: BitVec,
}

// type PauliString = Vec<PauliBit>;
impl PauliString{
    pub fn new(pauli_x: BitVec, pauli_z: BitVec) -> Self{
        assert!(pauli_x.len() == pauli_z.len());
        PauliString{ 
            x: pauli_x,
            z: pauli_z
        }
    }

    pub fn from_basis_int(i: usize, len: usize) ->Self{
        assert!(len>i);
        let mut x = BitVec::new();
        let mut z = BitVec::new();
        for j in 0..len{
            x.push(i==j);
            z.push(false);
        }
        for j in 0..len{
            x.push(false);
            z.push(i==j);
        }
        PauliString::new(x, z)
    }

    pub fn len(&self) -> usize{
        self.x.len()
    }

    pub fn from_text_string(mut pauli: String) -> Self{
        let mut x = BitVec::new();
        let mut z = BitVec::new();
        while !pauli.is_empty(){
            let (x_bit, z_bit) = match pauli.pop().unwrap(){
                'I' => (false, false),
                'X' => (true, false),
                'Y' => (true, true),
                'Z' => (false, true),
                _ => panic!("Letter not recognized"),
            };
            x.push(x_bit);
            z.push(z_bit);
        }
        PauliString::new(x, z)
    }

    pub(super) fn s(&mut self) {
        self.x ^= &self.z;
    }
    
    pub(super) fn v(&mut self) {
        self.z ^= &self.x;
    }

    pub(super) fn y_bitmask(self) -> BitVec{
        self.x & self.z
    }

}


pub(super) fn cx(control: &mut PauliString, target: &mut PauliString){
    assert!(control.len() == target.len());
    target.x ^= &control.x;
    control.z ^= &target.z;
}
