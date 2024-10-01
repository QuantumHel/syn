use bitvec::prelude::BitVec;
use std::iter;

use super::{
    pauli_string::{cx, PauliString},
    IndexType, PropagateClifford,
};

#[derive(PartialEq, Eq, Debug)]
pub struct CliffordTableau {
    stabilizers: Vec<PauliString>,
    x_signs: BitVec,
    z_signs: BitVec,
    // https://quantumcomputing.stackexchange.com/questions/28740/tracking-the-signs-of-the-inverse-tableau
}

impl CliffordTableau {
    /// Constructs a Clifford Tableau of `n` qubits initialized to the identity operation
    pub fn new(n: usize) -> Self {
        CliffordTableau {
            stabilizers: { (0..n).map(|i| PauliString::from_basis_int(i, n)).collect() },
            x_signs: BitVec::from_iter(iter::repeat(false).take(n)),
            z_signs: BitVec::from_iter(iter::repeat(false).take(n)),
        }
    }
}

impl PropagateClifford for CliffordTableau {
    fn cx(&mut self, control: IndexType, target: IndexType) -> &mut Self {
        let [control, target] = self.stabilizers.get_many_mut([control, target]).unwrap();
        cx(control, target);
        self
    }

    fn s(&mut self, target: IndexType) -> &mut Self {
        let chains_target = self.stabilizers.get_mut(target).unwrap();
        chains_target.s();
        // Defined for Phase gate in https://arxiv.org/pdf/quant-ph/0406196
        self.x_signs ^= chains_target.y_bitmask();
        self
    }

    fn v(&mut self, target: IndexType) -> &mut Self {
        let chains_target = self.stabilizers.get_mut(target).unwrap();
        self.z_signs ^= chains_target.y_bitmask();
        // TODO Double check if this works as intended.
        chains_target.s();
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bitvec::bitvec;
    use bitvec::prelude::Lsb0;

    #[test]
    fn test_clifford_tableau_constructor() {
        let ct_size = 3;
        let ct = CliffordTableau::new(ct_size);
        let x_1 = bitvec![1, 0, 0, 0, 0, 0];
        let z_1 = bitvec![0, 0, 0, 1, 0, 0];
        let pauli_1 = PauliString { x: x_1, z: z_1 };
        let x_2 = bitvec![0, 1, 0, 0, 0, 0];
        let z_2 = bitvec![0, 0, 0, 0, 1, 0];
        let pauli_2 = PauliString { x: x_2, z: z_2 };
        let x_3 = bitvec![0, 0, 1, 0, 0, 0];
        let z_3 = bitvec![0, 0, 0, 0, 0, 1];
        let pauli_3 = PauliString { x: x_3, z: z_3 };
        let x_signs = bitvec![0, 0, 0,];
        let z_signs = bitvec![0, 0, 0,];
        let clifford_tableau_ref = CliffordTableau {
            stabilizers: vec![pauli_1, pauli_2, pauli_3],
            x_signs,
            z_signs,
        };
        assert_eq!(ct, clifford_tableau_ref);
    }

    #[test]
    fn test_clifford_tableau_cx() {
        // implement a swap on qubits 0 and 1
        let mut ct = CliffordTableau::new(3);
        ct.cx(0, 1);
        ct.cx(1, 0);
        ct.cx(0, 1);

        let x_1 = bitvec![1, 0, 0, 0, 0, 0];
        let z_1 = bitvec![0, 0, 0, 1, 0, 0];
        let pauli_1 = PauliString { x: x_1, z: z_1 };
        let x_2 = bitvec![0, 1, 0, 0, 0, 0];
        let z_2 = bitvec![0, 0, 0, 0, 1, 0];
        let pauli_2 = PauliString { x: x_2, z: z_2 };
        let x_3 = bitvec![0, 0, 1, 0, 0, 0];
        let z_3 = bitvec![0, 0, 0, 0, 0, 1];
        let pauli_3 = PauliString { x: x_3, z: z_3 };
        let x_signs = bitvec![0, 0, 0,];
        let z_signs = bitvec![0, 0, 0,];
        let clifford_tableau_ref = CliffordTableau {
            stabilizers: vec![pauli_2, pauli_1, pauli_3],
            x_signs,
            z_signs,
        };
        assert_eq!(ct, clifford_tableau_ref);
    }

    #[test]
    fn test_clifford_tableau_s() {
        let mut ct = CliffordTableau::new(3);
        ct.s(2);

        let x_1 = bitvec![1, 0, 0, 0, 0, 0];
        let z_1 = bitvec![0, 0, 0, 1, 0, 0];
        let pauli_1 = PauliString { x: x_1, z: z_1 };
        let x_2 = bitvec![0, 1, 0, 0, 0, 0];
        let z_2 = bitvec![0, 0, 0, 0, 1, 0];
        let pauli_2 = PauliString { x: x_2, z: z_2 };
        let x_3 = bitvec![0, 0, 1, 0, 0, 0];
        let z_3 = bitvec![0, 0, 0, 0, 0, 1];
        let pauli_3 = PauliString { x: x_3, z: z_3 };
        let x_signs = bitvec![0, 0, 0,];
        let z_signs = bitvec![0, 0, 0,];
        let clifford_tableau_ref = CliffordTableau {
            stabilizers: vec![pauli_2, pauli_1, pauli_3],
            x_signs,
            z_signs,
        };
        assert_eq!(ct, clifford_tableau_ref);
    }
}
