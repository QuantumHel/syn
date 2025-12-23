mod helper;
pub mod naive;

pub use naive::NaivePauliPolynomialSynthesizer;

#[derive(Default, Clone)]
pub enum PauliPolynomialSynthStrategy {
    #[default]
    Naive,
}
