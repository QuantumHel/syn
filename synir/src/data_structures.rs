use crate::IndexType;

mod clifford_tableau;
mod pauli_polynomial;
mod pauli_string;

use bitvec::vec::BitVec;
pub use clifford_tableau::CliffordTableau;
pub use pauli_polynomial::PauliPolynomial;
pub use pauli_string::PauliString;

pub type Angle = f64;

pub trait HasAdjoint {
    fn adjoint(&self) -> Self;
}
pub trait PropagateClifford
where
    Self: Sized,
{
    fn cx(&mut self, control: IndexType, target: IndexType) -> &mut Self;
    fn s(&mut self, target: IndexType) -> &mut Self;
    fn v(&mut self, target: IndexType) -> &mut Self;

    fn s_dgr(&mut self, target: IndexType) -> &mut Self {
        self.z(target).s(target)
    }

    fn v_dgr(&mut self, target: IndexType) -> &mut Self {
        self.x(target).v(target)
    }

    fn x(&mut self, target: IndexType) -> &mut Self {
        self.v(target).v(target)
    }

    fn y(&mut self, target: IndexType) -> &mut Self {
        self.s_dgr(target).x(target).s(target)
    }

    fn z(&mut self, target: IndexType) -> &mut Self {
        self.s(target).s(target)
    }

    fn h(&mut self, target: IndexType) -> &mut Self {
        self.s(target).v(target).s(target)
    }

    fn cz(&mut self, control: IndexType, target: IndexType) -> &mut Self {
        self.h(target);
        self.cx(control, target);
        self.h(target)
    }
}

pub trait MaskedPropagateClifford
where
    Self: Sized,
{
    fn masked_cx(&mut self, control: IndexType, target: IndexType, mask: &BitVec) -> &mut Self;
    fn masked_s(&mut self, target: IndexType, mask: &BitVec) -> &mut Self;
    fn masked_v(&mut self, target: IndexType, mask: &BitVec) -> &mut Self;

    fn masked_s_dgr(&mut self, target: IndexType, mask: &BitVec) -> &mut Self {
        self.masked_z(target, mask).masked_s(target, mask)
    }

    fn masked_v_dgr(&mut self, target: IndexType, mask: &BitVec) -> &mut Self {
        self.masked_x(target, mask).masked_v(target, mask)
    }

    fn masked_x(&mut self, target: IndexType, mask: &BitVec) -> &mut Self {
        self.masked_v(target, mask).masked_v(target, mask)
    }

    fn masked_y(&mut self, target: IndexType, mask: &BitVec) -> &mut Self {
        self.masked_s_dgr(target, mask)
            .masked_x(target, mask)
            .masked_s(target, mask)
    }

    fn masked_z(&mut self, target: IndexType, mask: &BitVec) -> &mut Self {
        self.masked_s(target, mask).masked_s(target, mask)
    }

    fn masked_h(&mut self, target: IndexType, mask: &BitVec) -> &mut Self {
        self.masked_s(target, mask)
            .masked_v(target, mask)
            .masked_s(target, mask)
    }

    fn masked_cz(&mut self, control: IndexType, target: IndexType, mask: &BitVec) -> &mut Self {
        self.masked_h(target, mask);
        self.masked_cx(control, target, mask);
        self.masked_h(target, mask)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PauliLetter {
    I,
    X,
    Y,
    Z,
}

impl PauliLetter {
    pub fn new(x: bool, z: bool) -> Self {
        match (x, z) {
            (false, false) => PauliLetter::I,
            (true, false) => PauliLetter::X,
            (true, true) => PauliLetter::Y,
            (false, true) => PauliLetter::Z,
        }
    }
}
