pub mod custom_callback;
pub mod custom_pivots;
mod helper;
pub mod naive;

pub use custom_pivots::CustomPivotCliffordSynthesizer;
pub use naive::NaiveCliffordSynthesizer;

use crate::data_structures::CliffordTableau;

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
