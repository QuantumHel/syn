use std::collections::VecDeque;

use crate::data_structures::{CliffordTableau, PauliPolynomial};

use crate::ir::{CliffordGates, Gates, Synthesizer};

use crate::ir::{
    clifford_tableau::{
        custom_pivots::CustomPivotCliffordSynthesizer, naive::NaiveCliffordSynthesizer,
        CliffordTableauSynthStrategy,
    },
    pauli_polynomial::{naive::NaivePauliPolynomialSynthesizer, PauliPolynomialSynthStrategy},
};

#[derive(Default)]
pub struct PauliExponential {
    pauli_polynomials: VecDeque<PauliPolynomial>,
    clifford_tableau: CliffordTableau,
}

impl PauliExponential {
    pub fn new(
        pauli_polynomials: VecDeque<PauliPolynomial>,
        clifford_tableau: CliffordTableau,
    ) -> Self {
        PauliExponential {
            pauli_polynomials,
            clifford_tableau,
        }
    }
}

#[derive(Default)]
pub struct PauliExponentialSynthesizer {
    pauli_strategy: PauliPolynomialSynthStrategy,
    clifford_strategy: CliffordTableauSynthStrategy,
}

impl PauliExponentialSynthesizer {
    pub fn from_strategy(
        pauli_strategy: PauliPolynomialSynthStrategy,
        clifford_strategy: CliffordTableauSynthStrategy,
    ) -> Self {
        Self {
            pauli_strategy,
            clifford_strategy,
        }
    }

    pub fn set_pauli_strategy(
        &mut self,
        pauli_strategy: PauliPolynomialSynthStrategy,
    ) -> &mut Self {
        self.pauli_strategy = pauli_strategy;
        self
    }

    pub fn set_clifford_strategy(
        &mut self,
        clifford_strategy: CliffordTableauSynthStrategy,
    ) -> &mut Self {
        self.clifford_strategy = clifford_strategy;
        self
    }
}

impl<G> Synthesizer<PauliExponential, G, ()> for PauliExponentialSynthesizer
where
    G: CliffordGates + Gates,
{
    fn synthesize(&mut self, pauli_exponential: PauliExponential, repr: &mut G) {
        let PauliExponential {
            pauli_polynomials,
            clifford_tableau,
        } = pauli_exponential;

        let clifford_tableau = match self.pauli_strategy {
            PauliPolynomialSynthStrategy::Naive => {
                let mut pauli_synthesizer = NaivePauliPolynomialSynthesizer::default();
                pauli_synthesizer.set_clifford_tableau(clifford_tableau);
                pauli_synthesizer.synthesize(pauli_polynomials, repr)
            }
        };

        match &self.clifford_strategy {
            CliffordTableauSynthStrategy::Naive => {
                let mut clifford_synthesizer = NaiveCliffordSynthesizer {};
                clifford_synthesizer.synthesize_adjoint(clifford_tableau, repr);
            }
            CliffordTableauSynthStrategy::Custom(custom_rows, custom_columns) => {
                let mut clifford_synthesizer = CustomPivotCliffordSynthesizer::default();
                clifford_synthesizer
                    .set_custom_columns(custom_columns.to_owned())
                    .set_custom_rows(custom_rows.to_owned());
                clifford_synthesizer.synthesize_adjoint(clifford_tableau, repr);
            }
        };
    }

    fn synthesize_adjoint(&mut self, _: PauliExponential, _: &mut G) {
        unimplemented!("Not required for PauliExponential")
    }
}
