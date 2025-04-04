pub mod custom_callback;
mod helper;
pub mod naive;

use crate::data_structures::CliffordTableau;
pub use custom_callback::CallbackCliffordSynthesizer;
pub use naive::NaiveCliffordSynthesizer;

use super::{AdjointSynthesizer, Synthesizer};

#[derive(Default)]
pub enum CliffordTableauSynthStrategy {
    #[default]
    Naive,
    Custom(Vec<usize>, Vec<usize>),
}

impl<T: AdjointSynthesizer<CliffordTableau, To>, To> Synthesizer<CliffordTableau, To> for T {
    fn synthesize(&mut self, ir: CliffordTableau, repr: &mut To) {
        self.synthesize_adjoint(ir, repr)
    }
}
