use std::collections::VecDeque;
use std::fmt;

use crate::data_structures::{CliffordTableau, HasAdjoint, PauliPolynomial};

use crate::ir::{CliffordGates, Gates, Synthesizer};

use crate::ir::clifford_tableau::CallbackCliffordSynthesizer;
use crate::ir::clifford_tableau::NaiveCliffordSynthesizer;
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

impl fmt::Display for PauliExponential {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut out: String = String::new();
        if self.pauli_polynomials.is_empty() {
            out.push_str("Angles || No poly ||");
            out.push_str(&self.clifford_tableau.get_first_line_string());
            for i in 0..&self.clifford_tableau.column(0).len() / 2 {
                out.push_str("QB");
                out.push_str(&i.to_string());
                if i < 10 {
                    out.push(' ');
                }
                out.push_str("   || _______ || ");
                out.push_str(&self.clifford_tableau.get_line_string(i).as_str());
                out.push_str("\n");
            }
        } else {
            out.push_str("Angles ||");
            for pp in &self.pauli_polynomials {
                out.push_str(pp.get_first_line_string().as_str());
                out.push('|');
            }
            out.push_str(&self.clifford_tableau.get_first_line_string());

            for i in 0..&self.clifford_tableau.column(0).len() / 2 {
                out.push_str("QB");
                out.push_str(&i.to_string());
                if i < 10 {
                    out.push(' ');
                }
                out.push_str("   || ");
                for pp in &self.pauli_polynomials {
                    out.push_str(pp.get_line_string(i).as_str());
                    out.push_str("| ");
                }
                out.push_str(&self.clifford_tableau.get_line_string(i).as_str());
                out.push_str("\n");
            }
        }

        write!(f, "{}", out)?;
        writeln!(f)
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

        let clifford_tableau = match self.pauli_strategy {
            PauliPolynomialSynthStrategy::Naive => {
                let mut pauli_synthesizer = NaivePauliPolynomialSynthesizer::default();
                pauli_synthesizer.set_clifford_tableau(clifford_tableau);
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
        };
    }
}
