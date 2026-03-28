use core::num;
use std::collections::VecDeque;

use crate::{
    architecture::connectivity::Connectivity,
    data_structures::{CliffordTableau, HasAdjoint, PauliExponential, PauliPolynomial}, ir::AdjointSynthesizer,
};

use crate::ir::pauli_polynomial::psgs::PSGSPauliPolynomialSynthesizer;
use crate::ir::{CliffordGates, Gates, Synthesizer};

use crate::ir::clifford_tableau::NaiveCliffordSynthesizer;
use crate::ir::clifford_tableau::{CallbackCliffordSynthesizer, PermRowColCliffordSynthesizer};
use crate::ir::{
    clifford_tableau::CliffordTableauSynthStrategy,
    pauli_polynomial::{naive::NaivePauliPolynomialSynthesizer, PauliPolynomialSynthStrategy},
};

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
        let ct = match pauli_polynomials.is_empty() {
            true => CliffordTableau::new(num_qubits), // Skip PauliPolynomial synthesis
            false => match self.pauli_strategy {
                PauliPolynomialSynthStrategy::Naive => {
                    let mut pauli_synthesizer = NaivePauliPolynomialSynthesizer::default();
                    pauli_synthesizer.synthesize(pauli_polynomials, repr)
                }
                PauliPolynomialSynthStrategy::PSGS => {
                    let mut pauli_synthesizer = PSGSPauliPolynomialSynthesizer::default();
                    pauli_synthesizer.set_connectivity(Connectivity::complete(num_qubits));
                    pauli_synthesizer.synthesize(pauli_polynomials, repr)
                }
            },
        };
        // NOTE: The `ct` contains the current transformation of the qubits
        //       i.e. the Cliffords as synthesized - from left to right in the circuit
        //       i.e. the Cliffords are appended to the ct
        //       i.e. `ct.adjoint()` needs to be applied before `clifford_tableau`. 
        
        /*
        let (mut lhs, mut rhs) = (ct, clifford_tableau);
        println!("lhs.compose(rhs) {}", lhs.compose(&rhs));
        println!("lhs.compose(rhs.adjoint() {}", lhs.compose(&rhs.adjoint()));
        println!("lhs.compose(rhs).adjoint() {}", lhs.compose(&rhs).adjoint());
        println!("lhs.compose(rhs.adjoint()).adjoint() {}", lhs.compose(&rhs.adjoint()).adjoint());
        println!("lhs_adj.compose(rhs) {}", lhs.adjoint().compose(&rhs));
        println!("lhs_adj.compose(rhs.adjoint() {}", lhs.adjoint().compose(&rhs.adjoint()));
        println!("lhs_adj.compose(rhs).adjoint() {}", lhs.adjoint().compose(&rhs).adjoint());
        println!("lhs_adj.compose(rhs.adjoint()).adjoint() {}", lhs.adjoint().compose(&rhs.adjoint()).adjoint());
        (rhs, lhs) = (lhs, rhs);
        
        println!("rhs.compose(rhs) {}", lhs.compose(&rhs));
        println!("rhs.compose(rhs.adjoint() {}", lhs.compose(&rhs.adjoint()));
        println!("rhs.compose(rhs).adjoint() {}", lhs.compose(&rhs).adjoint());
        println!("rhs.compose(rhs.adjoint()).adjoint() {}", lhs.compose(&rhs.adjoint()).adjoint());
        println!("rhs_adj.compose(lhs) {}", lhs.adjoint().compose(&rhs));
        println!("rhs_adj.compose(lhs.adjoint() {}", lhs.adjoint().compose(&rhs.adjoint()));
        println!("rhs_adj.compose(lhs).adjoint() {}", lhs.adjoint().compose(&rhs).adjoint());
        println!("rhs_adj.compose(lhs.adjoint()).adjoint() {}", lhs.adjoint().compose(&rhs.adjoint()).adjoint());
        (rhs, lhs) = (lhs, rhs);
        //let combined_ct = lhs.adjoint().compose(&rhs.adjoint()).adjoint();
        let combined_ct = lhs.compose(&rhs.adjoint());
        */
        let combined_ct = ct.compose(&clifford_tableau.adjoint());
        //let combined_ct = ct.adjoint().compose(&clifford_tableau.adjoint()).adjoint();
        let final_ct = match &self.clifford_strategy {
            CliffordTableauSynthStrategy::Naive => {
                let mut clifford_synthesizer = NaiveCliffordSynthesizer::default();
                clifford_synthesizer.synthesize_adjoint(combined_ct, repr)
            }
            CliffordTableauSynthStrategy::Custom(custom_rows, custom_columns) => {
                let mut clifford_synthesizer = CallbackCliffordSynthesizer::custom_pivot(
                    custom_columns.to_owned(),
                    custom_rows.to_owned(),
                );
                clifford_synthesizer.synthesize_adjoint(combined_ct, repr)
            }
            CliffordTableauSynthStrategy::PermRowCol => {
                let mut clifford_synthesizer =
                    PermRowColCliffordSynthesizer::new(Connectivity::complete(num_qubits));
                clifford_synthesizer.synthesize_adjoint(combined_ct, repr)
            }
        };
        let final_perm = final_ct.get_permutation();
        match final_perm {
            Some(perm) => repr.add_final_permutation(perm),
            None => panic!("Final state was not a permutation: {final_ct}"),
        }
    }
}

pub fn print_pp_help(pauli_polynomials: &VecDeque<PauliPolynomial>) {
    for pp in pauli_polynomials {
        for i in 0..pp.length() {
            println!("{}, {}", pp.pauli_gadget(i), pp.angle(i));
        }
        println!("--");
    }
}
