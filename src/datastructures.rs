use crate::IndexType;

mod clifford_tableau;
mod pauli_polynomial;
mod pauli_string;

pub use clifford_tableau::CliffordTableau;
pub use pauli_polynomial::PauliPolynomial;
pub(crate) use pauli_string::PauliString;
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
