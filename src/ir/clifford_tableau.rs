pub mod custom_callback;
pub mod custom_pivots;
mod helper;
pub mod naive;

pub use custom_pivots::CustomPivotCliffordSynthesizer;
pub use naive::NaiveCliffordSynthesizer;

use crate::data_structures::CliffordTableau;

use super::{HasAdjoint, Synthesizer};

#[derive(Default)]
pub enum CliffordTableauSynthStrategy {
    #[default]
    Naive,
    Custom(Vec<usize>, Vec<usize>),
}

impl<T: HasAdjoint<CliffordTableau, To, Return>, To, Return>
    Synthesizer<CliffordTableau, To, Return> for T
{
    fn synthesize(&mut self, ir: CliffordTableau, repr: &mut To) -> Return {
        self.synthesize_adjoint(ir, repr)
    }
}
