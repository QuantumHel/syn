use crate::{data_structures::HasAdjoint, IndexType};

pub mod clifford_tableau;
pub mod pauli_exponential;
pub mod pauli_polynomial;

pub trait CliffordGates {
    fn s(&mut self, target: IndexType);
    fn v(&mut self, target: IndexType);
    fn s_dgr(&mut self, target: IndexType);
    fn v_dgr(&mut self, target: IndexType);
    fn x(&mut self, target: IndexType);
    fn y(&mut self, target: IndexType);
    fn z(&mut self, target: IndexType);
    fn h(&mut self, target: IndexType);
    fn cx(&mut self, control: IndexType, target: IndexType);
    fn cz(&mut self, control: IndexType, target: IndexType);
}

pub trait Gates {
    fn rx(&mut self, target: IndexType, angle: f64);
    fn ry(&mut self, target: IndexType, angle: f64);
    fn rz(&mut self, target: IndexType, angle: f64);
}

pub trait AdjointSynthesizer<From, To, Returns = ()>
where
    From: HasAdjoint,
{
    fn synthesize_adjoint(&mut self, ir: From, external_repr: &mut To) -> Returns;
}

pub trait Synthesizer<From, To, Returns = ()> {
    fn synthesize(&mut self, ir: From, repr: &mut To) -> Returns;
}
