mod common;

use std::collections::VecDeque;

use common::{parse_clifford_commands, MockCircuit, MockCommand};
use syn::data_structures::{CliffordTableau, PauliPolynomial};
use syn::ir::clifford_tableau::naive::NaiveCliffordSynthesizer;
use syn::ir::clifford_tableau::{CliffordTableauSynthStrategy, CliffordTableauSynthesizer};
use syn::ir::pauli_exponential::naive::NaivePauliExponentialSynthesizer;
use syn::ir::pauli_exponential::{PauliExponential, PauliExponentialSynthesizer};
use syn::ir::pauli_polynomial::PauliPolynomialSynthStrategy;

fn setup_simple_pe() -> PauliExponential {
    let ham = vec![("IZZZ", 0.3)];

    let pauli_polynomial = PauliPolynomial::from_hamiltonian(ham);
    let clifford_tableau = CliffordTableau::new(4);
    PauliExponential::new(VecDeque::from([pauli_polynomial]), clifford_tableau)
}

fn setup_complex_pe() -> PauliExponential {
    let ham = vec![("IXYZ", 0.3), ("XXII", 0.7), ("YYII", 0.12)];

    let pauli_polynomial = PauliPolynomial::from_hamiltonian(ham);
    let clifford_tableau = CliffordTableau::new(4);
    PauliExponential::new(VecDeque::from([pauli_polynomial]), clifford_tableau)
}

#[test]
fn test_naive_pauli_exponential_synthesis() {
    let pe = setup_simple_pe();
    let mut mock = MockCircuit::new();
    let mut synthesizer = NaivePauliExponentialSynthesizer::from_strategy(
        PauliPolynomialSynthStrategy::Naive,
        CliffordTableauSynthStrategy::Naive,
    );
    synthesizer.synthesize(pe, &mut mock);

    let ref_commands = [
        MockCommand::CX(1, 2),
        MockCommand::CX(2, 3),
        MockCommand::Rz(3, 0.3),
        MockCommand::CX(1, 2),
        MockCommand::CX(2, 3),
    ];

    assert_eq!(mock.commands(), &ref_commands);
}

#[test]
fn test_naive_pauli_exponential_complex() {
    let pe = setup_complex_pe();
    let mut mock = MockCircuit::new();

    let mut synthesizer = NaivePauliExponentialSynthesizer::from_strategy(
        PauliPolynomialSynthStrategy::Naive,
        CliffordTableauSynthStrategy::Naive,
    );
    synthesizer.synthesize(pe, &mut mock);

    let input = [
        MockCommand::H(1),
        MockCommand::S(2),
        MockCommand::CX(1, 2),
        MockCommand::CX(2, 3),
        MockCommand::H(0),
        MockCommand::CX(0, 1),
        MockCommand::H(0),
        MockCommand::H(2),
        MockCommand::H(3),
        MockCommand::CX(0, 1),
        MockCommand::CX(1, 2),
        MockCommand::CX(2, 3),
    ];

    let ref_ct = parse_clifford_commands(4, &input);
    let mut mock_ct = MockCircuit::new();

    let mut cliff_synthesizer = NaiveCliffordSynthesizer::default();

    cliff_synthesizer.synthesize(ref_ct.clone(), &mut mock_ct);

    let mock_ct_ref_commands = [
        MockCommand::S(2),
        MockCommand::H(1),
        MockCommand::CX(0, 1),
        MockCommand::CX(0, 2),
        MockCommand::H(1),
        MockCommand::H(2),
        MockCommand::CX(1, 0),
        MockCommand::CX(2, 0),
        MockCommand::H(1),
        MockCommand::S(1),
        MockCommand::S(2),
        MockCommand::CX(1, 2),
        MockCommand::CX(3, 2),
        MockCommand::CX(2, 3),
        MockCommand::CX(3, 2),
        MockCommand::H(2),
        MockCommand::CX(3, 2),
        MockCommand::S(3),
        MockCommand::X(1),
        MockCommand::X(2),
        MockCommand::X(3),
        MockCommand::Z(3),
    ];

    assert_eq!(mock_ct.commands(), &mock_ct_ref_commands);

    let ref_commands = [
        MockCommand::H(1),
        MockCommand::S(2),
        MockCommand::CX(1, 2),
        MockCommand::CX(2, 3),
        MockCommand::Rz(3, 0.3),
        MockCommand::H(0),
        MockCommand::CX(0, 1),
        MockCommand::Rz(1, 0.7),
        MockCommand::H(0),
        MockCommand::H(2),
        MockCommand::H(3),
        MockCommand::CX(0, 1),
        MockCommand::CX(1, 2),
        MockCommand::CX(2, 3),
        MockCommand::Rz(3, -0.12),
        MockCommand::S(2),
        MockCommand::H(1),
        MockCommand::CX(0, 1),
        MockCommand::CX(0, 2),
        MockCommand::H(1),
        MockCommand::H(2),
        MockCommand::CX(1, 0),
        MockCommand::CX(2, 0),
        MockCommand::H(1),
        MockCommand::S(1),
        MockCommand::S(2),
        MockCommand::CX(1, 2),
        MockCommand::CX(3, 2),
        MockCommand::CX(2, 3),
        MockCommand::CX(3, 2),
        MockCommand::H(2),
        MockCommand::CX(3, 2),
        MockCommand::S(3),
        MockCommand::X(1),
        MockCommand::X(2),
        MockCommand::X(3),
        MockCommand::Z(3),
    ];

    assert_eq!(mock.commands(), &ref_commands);
}
