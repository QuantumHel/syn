mod helper;
pub mod naive;
pub mod psgs;

pub use naive::NaivePauliPolynomialSynthesizer;

#[derive(Default, Clone)]
pub enum PauliPolynomialSynthStrategy {
    #[default]
    Naive,
    PSGS,
}
