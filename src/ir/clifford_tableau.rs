use crate::data_structures::CliffordTableau;

pub mod custom_pivots;
mod helper;
pub mod naive;

pub trait CliffordTableauSynthesizer<G> {
    fn synthesize(&mut self, clifford_tableau: CliffordTableau, external_repr: &mut G);
    fn synthesize_adjoint(&mut self, clifford_tableau: CliffordTableau, external_repr: &mut G);
}

#[derive(Default)]
pub enum CliffordTableauSynthStrategy {
    #[default]
    Naive,
    Custom(Vec<usize>, Vec<usize>),
}
