mod helper;
pub mod naive;
pub mod psgs;

pub use naive::NaivePauliPolynomialSynthesizer;

#[derive(Default)]
pub enum PauliPolynomialSynthStrategy {
    #[default]
    Naive,
    PSGS,
}
