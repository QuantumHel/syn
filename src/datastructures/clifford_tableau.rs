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