mod common;

use std::collections::VecDeque;

use synir::{
    data_structures::{CliffordTableau, PauliExponential, PauliPolynomial},
    ir::{
        clifford_tableau::CliffordTableauSynthStrategy,
        pauli_exponential::PauliExponentialSynthesizer,
        pauli_polynomial::PauliPolynomialSynthStrategy, Synthesizer,
    },
};

use crate::common::mock_circuit::MockCircuit;

fn run_synthesizer(pe: PauliExponential) -> MockCircuit {
    let mut mock: MockCircuit = MockCircuit::new();
    let mut synthesizer = PauliExponentialSynthesizer::from_strategy(
        PauliPolynomialSynthStrategy::Naive,
        CliffordTableauSynthStrategy::Naive,
    );
    synthesizer.synthesize(pe, &mut mock);
    return mock;
}

#[test]
fn test_empty_pe() {
    let pe = PauliExponential::new(VecDeque::new(), CliffordTableau::new(5));
    let mock = run_synthesizer(pe);
    let ref_commands = [];
    assert_eq!(mock.commands(), &ref_commands);
}

#[test]
fn test_empty_pp_in_pe() {
    let pe = PauliExponential::new(
        VecDeque::from([PauliPolynomial::from_components(vec![], vec![], 5)]),
        CliffordTableau::new(5),
    );
    let mock = run_synthesizer(pe);
    let ref_commands = [];
    assert_eq!(mock.commands(), &ref_commands);
}
