pub mod custom_pivots;
mod helper;
pub mod naive;

pub trait CliffordTableauSynthesizer<G> {
    fn synthesize(&mut self, external_repr: &mut G);
    fn synthesize_adjoint(&mut self, external_repr: &mut G);
}

pub enum CliffordTableauSynthStrategy {
    Naive,
    Custom(Vec<usize>, Vec<usize>),
}
