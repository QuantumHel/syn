use super::{AdjointSynthesizer, Synthesizer};
use crate::data_structures::{CliffordTableau, HasAdjoint};

pub use custom_callback::CallbackCliffordSynthesizer;
pub use naive::NaiveCliffordSynthesizer;
pub use permrowcol::PermRowColCliffordSynthesizer;

mod custom_callback;
mod helper;
mod naive;
mod permrowcol;

#[derive(Default)]
pub enum CliffordTableauSynthStrategy {
    #[default]
    Naive,
    PermRowCol,
    Custom(Vec<usize>, Vec<usize>),
}

impl<T: AdjointSynthesizer<CliffordTableau, To>, To> Synthesizer<CliffordTableau, To> for T {
    fn synthesize(&mut self, ir: CliffordTableau, repr: &mut To) {
        let ir = ir.adjoint();
        self.synthesize_adjoint(ir, repr)
    }
}
