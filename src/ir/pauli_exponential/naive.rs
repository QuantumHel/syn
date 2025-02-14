use crate::ir::{
    clifford_tableau::{
        custom_pivots::CustomPivotCliffordSynthesizer, naive::NaiveCliffordSynthesizer,
        CliffordTableauSynthStrategy, CliffordTableauSynthesizer,
    },
    pauli_polynomial::{
        naive::NaivePauliPolynomialSynthesizer, PauliPolynomialSynthStrategy,
        PauliPolynomialSynthesizer,
    },
    CliffordGates, Gates,
};

use super::{PauliExponential, PauliExponentialSynthesizer};

pub struct NaivePauliExponentialSynthesizer {
    pauli_exponential: PauliExponential,
    pauli_strategy: PauliPolynomialSynthStrategy,
    clifford_strategy: CliffordTableauSynthStrategy,
}

impl NaivePauliExponentialSynthesizer {
    pub fn new(
        pauli_exponential: PauliExponential,
        pauli_strategy: PauliPolynomialSynthStrategy,
        clifford_strategy: CliffordTableauSynthStrategy,
    ) -> Self {
        Self {
            pauli_exponential,
            pauli_strategy,
            clifford_strategy,
        }
    }
}

impl<G> PauliExponentialSynthesizer<G> for NaivePauliExponentialSynthesizer
where
    G: CliffordGates + Gates,
{
    fn synthesize(&mut self, repr: &mut G) {
        let PauliExponential {
            pauli_polynomial,
            clifford_tableau,
        } = std::mem::take(&mut self.pauli_exponential);

        let clifford_tableau = match self.pauli_strategy {
            PauliPolynomialSynthStrategy::Naive => {
                let mut pauli_synthesizer =
                    NaivePauliPolynomialSynthesizer::new(pauli_polynomial, clifford_tableau);
                pauli_synthesizer.synthesize(repr)
            }
        };

        match &self.clifford_strategy {
            CliffordTableauSynthStrategy::Naive => {
                let mut clifford_synthesizer = NaiveCliffordSynthesizer::new(clifford_tableau);
                clifford_synthesizer.synthesize(repr);
            }
            CliffordTableauSynthStrategy::Custom(custom_rows, custom_columns) => {
                let mut clifford_synthesizer = CustomPivotCliffordSynthesizer::new(
                    clifford_tableau,
                    custom_rows.to_owned(),
                    custom_columns.to_owned(),
                );
                clifford_synthesizer.synthesize(repr);
            }
        };
    }
}
