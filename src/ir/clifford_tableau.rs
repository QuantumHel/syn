use crate::data_structures::CliffordTableau;

pub mod custom_pivots;
mod helper;
pub mod naive;

pub trait CliffordTableauSynthesizer<G> {
    fn synthesize(&mut self, mut clifford_tableau: CliffordTableau, repr: &mut G) {
        clifford_tableau = clifford_tableau.adjoint();
        CliffordTableauSynthesizer::<G>::synthesize_adjoint(self, clifford_tableau, repr)
    }
    fn synthesize_adjoint(&mut self, clifford_tableau: CliffordTableau, external_repr: &mut G);
}

#[derive(Default)]
pub enum CliffordTableauSynthStrategy {
    #[default]
    Naive,
    Custom(Vec<usize>, Vec<usize>),
}
