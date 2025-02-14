use crate::data_structures::CliffordTableau;

pub mod naive;

pub trait PauliPolynomialSynthesizer<G> {
    fn synthesize(&mut self, external_repr: &mut G) -> CliffordTableau;
}

pub enum PauliPolynomialSynthStrategy {
    Naive,
}
