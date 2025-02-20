pub mod custom_callback;
pub mod custom_pivots;
mod helper;
pub mod naive;

pub use custom_pivots::CustomPivotCliffordSynthesizer;
pub use naive::NaiveCliffordSynthesizer;

#[derive(Default, PartialEq, Eq)]
pub enum CliffordTableauSynthStrategy {
    #[default]
    Naive,
    Custom(Vec<usize>, Vec<usize>),
}
