use crate::{
    data_structures::{CliffordTableau, PauliPolynomial},
    synthesis_methods::naive::Naive,
};

use super::{CliffordGates, Gates};

struct PauliExponential {
    pauli_polynomial: PauliPolynomial,
    clifford_tableau: CliffordTableau,
}

pub struct PauliExponentialSynthesizer;
impl<G> Naive<PauliExponential, G> for PauliExponentialSynthesizer
where
    G: CliffordGates + Gates,
{
    fn run_naive(program: PauliExponential, external_rep: &mut G) {
        todo!()
    }
}
