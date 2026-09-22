//! Clifford tableau synthesis.
//!
//! Every algorithm here reduces a tableau to the identity by emitting gates, so they all implement
//! [`AdjointSynthesizer`] rather than [`Synthesizer`]; the blanket impl at the bottom of this
//! module supplies the forward direction. Each returns the residual tableau, which is a bare qubit
//! permutation once reduction has finished.

use super::{AdjointSynthesizer, Synthesizer};
use crate::data_structures::{CliffordTableau, HasAdjoint};

pub use custom_callback::CallbackCliffordSynthesizer;
pub use naive::NaiveCliffordSynthesizer;
pub use permrowcol::PermRowColCliffordSynthesizer;

mod custom_callback;
mod helper;
mod naive;
mod permrowcol;

/// Which tableau synthesis algorithm
/// [`PauliExponentialSynthesizer`](crate::ir::pauli_exponential::PauliExponentialSynthesizer)
/// should use.
#[derive(Default, Clone)]
pub enum CliffordTableauSynthStrategy {
    /// [`NaiveCliffordSynthesizer`] — reduce qubit by qubit, ignoring connectivity.
    #[default]
    Naive,
    /// [`PermRowColCliffordSynthesizer`] — connectivity-aware, choosing its own pivot order.
    PermRowCol,
    /// [`CallbackCliffordSynthesizer`] driven by a fixed pivot schedule, as `(rows, columns)`.
    ///
    /// Mostly a testing and benchmarking hook: it makes a synthesis run reproducible by pinning
    /// the choices the other strategies make heuristically.
    Custom(Vec<usize>, Vec<usize>),
}

/// Every tableau [`AdjointSynthesizer`] is also a [`Synthesizer`]: invert first, then reduce.
impl<T: AdjointSynthesizer<CliffordTableau, To, Returns>, To, Returns>
    Synthesizer<CliffordTableau, To, Returns> for T
{
    fn synthesize(&mut self, ir: CliffordTableau, repr: &mut To) -> Returns {
        let ir = ir.adjoint();
        return self.synthesize_adjoint(ir, repr);
    }
}
