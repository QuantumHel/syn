mod common;

use std::collections::VecDeque;

use common::{parse_clifford_commands, MockCircuit, MockCommand};
use syn::data_structures::{CliffordTableau, PauliPolynomial};
use syn::ir::pauli_polynomial::NaivePauliPolynomialSynthesizer;
use syn::ir::Synthesizer;

fn setup_complex_pp() -> VecDeque<PauliPolynomial> {
    let ham_1 = vec![("IZZZ", 0.3)];
    let ham_2 = vec![("XXII", 0.7)];

    let pp_1 = PauliPolynomial::from_hamiltonian(ham_1);
    let pp_2 = PauliPolynomial::from_hamiltonian(ham_2);
    VecDeque::from([pp_1, pp_2])
}

fn setup_simple_pp() -> VecDeque<PauliPolynomial> {
    let ham = vec![("IXYZ", 0.3)];

    let pauli_polynomial = PauliPolynomial::from_hamiltonian(ham);

    VecDeque::from([pauli_polynomial])
}

#[test]
fn test_naive_pauli_exponential_synthesis() {
    let pp = setup_simple_pp();
    let mut mock = MockCircuit::new();
    let mut synthesizer = NaivePauliPolynomialSynthesizer::default();
    synthesizer.set_clifford_tableau(CliffordTableau::new(4));
    let ct = synthesizer.synthesize(pp, &mut mock);

    let ref_commands = [
        MockCommand::H(1),
        MockCommand::S(2),
        MockCommand::CX(1, 2),
        MockCommand::CX(2, 3),
        MockCommand::Rz(3, 0.3),
    ];

    let ref_clifford_commands = [
        MockCommand::H(1),
        MockCommand::S(2),
        MockCommand::CX(1, 2),
        MockCommand::CX(2, 3),
    ];

    assert_eq!(mock.commands(), &ref_commands);
    assert_eq!(ct, parse_clifford_commands(4, &ref_clifford_commands));
}

#[test]
fn test_naive_pauli_exponential_synthesis_complex() {
    let pp = setup_complex_pp();
    let mut mock = MockCircuit::new();
    let mut synthesizer = NaivePauliPolynomialSynthesizer::default();
    synthesizer.set_clifford_tableau(CliffordTableau::new(4));
    let ct = synthesizer.synthesize(pp, &mut mock);

    let ref_commands = [
        MockCommand::CX(1, 2),
        MockCommand::CX(2, 3),
        MockCommand::Rz(3, 0.3),
        MockCommand::H(0),
        MockCommand::H(1),
        MockCommand::H(2),
        MockCommand::H(3),
        MockCommand::CX(0, 1),
        MockCommand::CX(1, 2),
        MockCommand::CX(2, 3),
        MockCommand::Rz(3, 0.7),
    ];

    let ref_clifford_commands = [
        MockCommand::CX(1, 2),
        MockCommand::CX(2, 3),
        MockCommand::H(0),
        MockCommand::H(1),
        MockCommand::H(2),
        MockCommand::H(3),
        MockCommand::CX(0, 1),
        MockCommand::CX(1, 2),
        MockCommand::CX(2, 3),
    ];

    assert_eq!(mock.commands(), &ref_commands);
    assert_eq!(ct, parse_clifford_commands(4, &ref_clifford_commands));
}
