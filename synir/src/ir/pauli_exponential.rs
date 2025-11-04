use std::collections::VecDeque;

use crate::architecture::connectivity::Connectivity;
use crate::data_structures::{CliffordTableau, HasAdjoint, PauliPolynomial};

use crate::ir::pauli_polynomial::psgs::PSGSPauliPolynomialSynthesizer;
use crate::ir::{CliffordGates, Gates, Synthesizer};

use crate::ir::clifford_tableau::NaiveCliffordSynthesizer;
use crate::ir::clifford_tableau::{CallbackCliffordSynthesizer, PermRowColCliffordSynthesizer};
use crate::ir::{
    clifford_tableau::CliffordTableauSynthStrategy,
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

impl<G> Synthesizer<PauliExponential, G> for PauliExponentialSynthesizer
where
    G: CliffordGates + Gates,
{
    fn synthesize(&mut self, pauli_exponential: PauliExponential, repr: &mut G) {
        let PauliExponential {
            pauli_polynomials,
            clifford_tableau,
        } = pauli_exponential;
        let num_qubits = clifford_tableau.size();
        let clifford_tableau = match self.pauli_strategy {
            PauliPolynomialSynthStrategy::Naive => {
                let mut pauli_synthesizer = NaivePauliPolynomialSynthesizer::default();
                pauli_synthesizer.synthesize(pauli_polynomials, repr)
            }
            PauliPolynomialSynthStrategy::PSGS => {
                let mut pauli_synthesizer = PSGSPauliPolynomialSynthesizer::default();
                pauli_synthesizer.set_connectivity(Connectivity::complete(num_qubits));
                pauli_synthesizer.synthesize(pauli_polynomials, repr)
            }
        };

        match &self.clifford_strategy {
            CliffordTableauSynthStrategy::Naive => {
                let mut clifford_synthesizer = NaiveCliffordSynthesizer::default();
                clifford_synthesizer.synthesize(clifford_tableau.adjoint(), repr);
            }
            CliffordTableauSynthStrategy::Custom(custom_rows, custom_columns) => {
                let mut clifford_synthesizer = CallbackCliffordSynthesizer::custom_pivot(
                    custom_columns.to_owned(),
                    custom_rows.to_owned(),
                );
                clifford_synthesizer.synthesize(clifford_tableau.adjoint(), repr);
            }
            CliffordTableauSynthStrategy::PermRowCol => {
                let size = clifford_tableau.size();
                let mut clifford_synthesizer =
                    PermRowColCliffordSynthesizer::new(Connectivity::complete(size));
                clifford_synthesizer.synthesize(clifford_tableau.adjoint(), repr);
            }
        };
    }
}
