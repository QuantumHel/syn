//! Pauli polynomial synthesis.
//!
//! Both algorithms consume a `VecDeque<PauliPolynomial>` and return the
//! [`CliffordTableau`](crate::data_structures::CliffordTableau) describing the basis changes they
//! left in place. They deliberately do *not* un-compute those basis changes: the caller folds the
//! residual tableau into the Clifford half of the normal form and synthesizes it once, which is
//! cheaper than undoing the changes gadget by gadget. See
//! [`PauliExponentialSynthesizer`](crate::ir::pauli_exponential::PauliExponentialSynthesizer).

mod helper;
pub mod naive;
pub mod psgs;

pub use naive::NaivePauliPolynomialSynthesizer;

/// Which polynomial synthesis algorithm
/// [`PauliExponentialSynthesizer`](crate::ir::pauli_exponential::PauliExponentialSynthesizer)
/// should use.
#[derive(Default, Clone)]
pub enum PauliPolynomialSynthStrategy {
    /// [`NaivePauliPolynomialSynthesizer`] — one gadget at a time, no sharing.
    #[default]
    Naive,
    /// [`PSGSPauliPolynomialSynthesizer`](psgs::PSGSPauliPolynomialSynthesizer) — recursive
    /// partitioning that shares CNOT structure between gadgets.
    PSGS,
}
