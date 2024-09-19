use std::{cell::RefCell, iter, rc::Rc};

use bitvec::prelude::BitVec;

use super::{pauli_string::{cx, PauliString}, PropagateClifford};

pub struct CliffordTableau{
    stabilizers: Vec<PauliString>,
    x_signs: BitVec,
    z_signs: BitVec,
    // https://quantumcomputing.stackexchange.com/questions/28740/tracking-the-signs-of-the-inverse-tableau
}

impl CliffordTableau {
    pub fn new(n: usize) -> Self{
        CliffordTableau{
            stabilizers: {
                let mut stab = Vec::new();
                for i in 0..n{
                    stab.push(PauliString::from_basis_int(i, n));
                }
                stab
            },
            x_signs: BitVec::from_iter(iter::repeat(false).take(n)),
            z_signs: BitVec::from_iter(iter::repeat(false).take(n)),
        }
    }
}

impl PropagateClifford for CliffordTableau{
    fn cx(&mut self, control: super::IndexType, target: super::IndexType) -> &mut Self {
        match control < target {
            true => {
                let split = self.stabilizers.split_at_mut(target); 
                cx(split.1.get_mut(0).unwrap(), split.0.get_mut(control).unwrap())
            },
            false => {
                let split = self.stabilizers.split_at_mut(control);
                cx(split.0.get_mut(target).unwrap(), split.1.get_mut(0).unwrap())
            },
        };
        self
    }

    fn s(&mut self, target: super::IndexType) -> &mut Self {
        let chains_target = self.stabilizers.get_mut(target).unwrap();
        chains_target.s();
        // Defined for Phase gate in https://arxiv.org/pdf/quant-ph/0406196
        self.x_signs ^= chains_target.to_owned().y_bitmask();
        self
    }

    fn v(&mut self, target: super::IndexType) -> &mut Self {
        let chains_target = self.stabilizers.get_mut(target).unwrap();
        self.z_signs ^= chains_target.to_owned().y_bitmask();
        // TODO Double check if this works as intended.
        chains_target.s();
        self
    }
}