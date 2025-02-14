use std::collections::VecDeque;

use crate::data_structures::{CliffordTableau, PauliPolynomial};

pub mod naive;

pub trait PauliPolynomialSynthesizer<G> {
    fn synthesize(
        &mut self,
        pauli_polynomial: VecDeque<PauliPolynomial>,
        external_repr: &mut G,
    ) -> CliffordTableau;
}

#[derive(Default)]
pub enum PauliPolynomialSynthStrategy {
    #[default]
    Naive,
}
