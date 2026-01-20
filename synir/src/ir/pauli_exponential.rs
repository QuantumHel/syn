use core::num;
use std::collections::VecDeque;

use crate::{
    architecture::connectivity::Connectivity,
    data_structures::{CliffordTableau, HasAdjoint, PauliExponential, PauliPolynomial},
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
        print_pp_help(&pauli_polynomials);
        println!("Before Synth {}", clifford_tableau);
        let num_qubits = clifford_tableau.size();
        let ct = match pauli_polynomials.is_empty(){
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
            }
        };
        println!("Before combine: {ct}");
        //let combined_ct = ct.adjoint().compose(&clifford_tableau.adjoint());
        let combined_ct = clifford_tableau.compose(&ct).adjoint();
        println!("After synth: {}", combined_ct);
        
        let final_ct = match &self.clifford_strategy {
            CliffordTableauSynthStrategy::Naive => {
                let mut clifford_synthesizer = NaiveCliffordSynthesizer::default();
                clifford_synthesizer.synthesize(combined_ct, repr)
            }
            CliffordTableauSynthStrategy::Custom(custom_rows, custom_columns) => {
                let mut clifford_synthesizer = CallbackCliffordSynthesizer::custom_pivot(
                    custom_columns.to_owned(),
                    custom_rows.to_owned(),
                );
                clifford_synthesizer.synthesize(combined_ct, repr)
            }
            CliffordTableauSynthStrategy::PermRowCol => {
                let mut clifford_synthesizer =
                    PermRowColCliffordSynthesizer::new(Connectivity::complete(num_qubits));
                clifford_synthesizer.synthesize(combined_ct, repr)
            }
        };
        println!("Final perm: {}", final_ct);
        let final_perm = final_ct.get_permutation();
        match final_perm {
            Some(perm) => repr.add_final_permutation(perm),
            None => panic!("Final state was not a permutation: {final_ct}")
        }
    }
}

fn print_pp_help(pauli_polynomials: &VecDeque<PauliPolynomial>){
    for pp in pauli_polynomials {
        for i in 0..pp.length(){
            println!("{}, {}", pp.pauli_gadget(i), pp.angle(i));
        }
        println!("--");
    }

}