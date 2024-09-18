use std::{cell::RefCell, rc::Rc};

use bitvec::prelude::BitVec;

use super::{pauli_string::PauliString, PropagateClifford};

pub struct CliffordTableau{
    stabilizers: Vec<Rc<RefCell<PauliString>>>,
    // todo: Make this into a union / type Angle
    phase_flips: BitVec
}

impl CliffordTableau {
    pub fn new(n: usize) -> Self{
        todo!("impl")
    }
}

impl PropagateClifford for CliffordTableau{
    fn cx(&mut self, control: super::IndexType, target: super::IndexType) -> &mut Self {
        todo!()
    }

    fn s(&mut self, target: super::IndexType) -> &mut Self {
        todo!()
    }

    fn v(&mut self, target: super::IndexType) -> &mut Self {
        todo!()
    }
}