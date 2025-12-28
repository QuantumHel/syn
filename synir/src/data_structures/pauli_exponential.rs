use std::collections::VecDeque;

use crate::data_structures::{CliffordTableau, PauliPolynomial, PropagateClifford};

#[derive(Default)]
pub struct PauliExponential {
    pub(crate) pauli_polynomials: VecDeque<PauliPolynomial>,
    pub(crate) clifford_tableau: CliffordTableau,
}

impl PauliExponential {
    pub fn new(
        pauli_polynomials: VecDeque<PauliPolynomial>,
        clifford_tableau: CliffordTableau,
    ) -> Self {
        PauliExponential {
            pauli_polynomials,
            clifford_tableau,
        }
    }

    pub fn chains(&mut self) -> &mut VecDeque<PauliPolynomial>{
        &mut self.pauli_polynomials
    }

    pub fn size(&self) -> usize{
        self.clifford_tableau.size()
    }
}

impl PropagateClifford for PauliExponential{
    fn cx(&mut self, control: crate::IndexType, target: crate::IndexType) -> &mut Self {
        self.pauli_polynomials.cx(control, target);
        self.clifford_tableau.cx(control, target);
        self
    }

    fn s(&mut self, target: crate::IndexType) -> &mut Self {
        self.pauli_polynomials.s(target);
        self.clifford_tableau.s(target);
        self
    }

    fn v(&mut self, target: crate::IndexType) -> &mut Self {
        self.pauli_polynomials.v(target);
        self.clifford_tableau.v(target);
        self
    }
}