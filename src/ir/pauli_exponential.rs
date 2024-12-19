use bitvec::{bitvec, order::Lsb0};
use itertools::Itertools;

use crate::{
    data_structures::{
        CliffordTableau, MaskedPropagateClifford, PauliLetter, PauliPolynomial, PropagateClifford,
    },
    synthesis_methods::naive::{Naive, NaiveAdjoint},
};

use super::{clifford_tableau::CliffordTableauSynthesizer, CliffordGates, Gates};

pub struct PauliExponential {
    pauli_polynomial: PauliPolynomial,
    clifford_tableau: CliffordTableau,
}

impl PauliExponential {
    pub fn new(pauli_polynomial: PauliPolynomial, clifford_tableau: CliffordTableau) -> Self {
        PauliExponential {
            pauli_polynomial,
            clifford_tableau,
        }
    }
}

pub struct PauliExponentialSynthesizer;
impl<G> Naive<PauliExponential, G> for PauliExponentialSynthesizer
where
    G: CliffordGates + Gates,
{
    fn run_naive(pauli_exponential: PauliExponential, repr: &mut G) {
        let PauliExponential {
            pauli_polynomial,
            mut clifford_tableau,
        } = pauli_exponential;
        let num_gadgets = pauli_polynomial.length();
        let gadget_length = pauli_polynomial.length();
        let mut mask = bitvec![usize, Lsb0; 1; gadget_length];
        for col in 0..num_gadgets {
            let mut affected_qubits = Vec::new();
            for (i, row) in pauli_polynomial.chains().iter().enumerate() {
                match row.pauli(col) {
                    PauliLetter::I => {}
                    PauliLetter::X => {
                        affected_qubits.push(i);
                        pauli_polynomial.masked_h(i, &mask);
                        clifford_tableau.h(i);
                        repr.h(i);
                    }
                    PauliLetter::Y => {
                        affected_qubits.push(i);
                        pauli_polynomial.masked_s(i, &mask);
                        clifford_tableau.s(i);
                        repr.s(i);
                    }
                    PauliLetter::Z => {
                        affected_qubits.push(i);
                    }
                }
            }
            if affected_qubits.len() > 1 {
                for (&control, &target) in affected_qubits.iter().tuple_windows() {
                    pauli_polynomial.masked_cx(control, target, &mask);
                    clifford_tableau.cx(control, target);
                    repr.cx(control, target);
                }
            }
            let last_qubit = *affected_qubits.last().unwrap();
            repr.rz(last_qubit, pauli_polynomial.angle(col));
            mask.replace(col, false);
        }
        CliffordTableauSynthesizer::run_naive_adjoint(&clifford_tableau, repr);
    }
}
