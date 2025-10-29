use std::collections::VecDeque;

use crate::common::mock_circuit::{parse_clifford_commands, MockCircuit, MockCommand};
use crate::common::sample_pauli_poly::{setup_complex_pp, setup_simple_pp};
use synir::architecture::connectivity::Connectivity;
use synir::data_structures::{CliffordTableau, PauliPolynomial};
use synir::ir::pauli_polynomial::psgs::PSGSPauliPolynomialSynthesizer;
use synir::ir::Synthesizer;

fn run_synthesizer(pp: VecDeque<PauliPolynomial>) -> (MockCircuit, CliffordTableau) {
    let mut mock: MockCircuit = MockCircuit::new();
    let mut synthesizer = PSGSPauliPolynomialSynthesizer::default();
    synthesizer.set_clifford_tableau(CliffordTableau::new(4));
    let ct = synthesizer.synthesize(pp, &mut mock);
    return (mock, ct);
}

#[test]
fn test_psgs_pauli_exponential_synthesis_simple() {
    let pp = setup_simple_pp();
    let mut mock = MockCircuit::new();
    let mut synthesizer = PSGSPauliPolynomialSynthesizer::default();
    synthesizer.set_clifford_tableau(CliffordTableau::new(4));
    synthesizer.set_connectivity(Connectivity::complete(4));
    let ct = synthesizer.synthesize(pp, &mut mock);

    let ref_commands = [
        MockCommand::S(1),
        MockCommand::CX(3, 1),
        MockCommand::V(2),
        MockCommand::CX(2, 1),
        MockCommand::Ry(1, 0.3),
    ];

    let ref_clifford_commands = [
        MockCommand::S(1),
        MockCommand::CX(3, 1),
        MockCommand::V(2),
        MockCommand::CX(2, 1),
    ];

    assert_eq!(mock.commands(), &ref_commands);
    assert_eq!(ct, parse_clifford_commands(4, &ref_clifford_commands));
}

#[test]
fn test_psgs_pauli_exponential_synthesis_complex() {
    let pp = setup_complex_pp();
    let mut mock = MockCircuit::new();
    let mut synthesizer = PSGSPauliPolynomialSynthesizer::default();
    synthesizer.set_clifford_tableau(CliffordTableau::new(4));
    synthesizer.set_connectivity(Connectivity::complete(4));
    let ct = synthesizer.synthesize(pp, &mut mock);

    let ref_commands = [
        MockCommand::CX(3, 1),
        MockCommand::CX(2, 1),
        MockCommand::Rz(1, 0.3),
        MockCommand::H(1),
        MockCommand::S(0),
        MockCommand::CX(1, 0),
        MockCommand::Ry(0, 0.7),
    ];

    let ref_clifford_commands = [
        MockCommand::CX(3, 1),
        MockCommand::CX(2, 1),
        MockCommand::H(1),
        MockCommand::S(0),
        MockCommand::CX(1, 0),
    ];

    assert_eq!(mock.commands(), &ref_commands);
    assert_eq!(ct, parse_clifford_commands(4, &ref_clifford_commands));
}
