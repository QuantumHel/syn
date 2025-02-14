pub mod naive;

pub use naive::NaivePauliPolynomialSynthesizer;

#[derive(Default)]
pub enum PauliPolynomialSynthStrategy {
    #[default]
    Naive,
}
