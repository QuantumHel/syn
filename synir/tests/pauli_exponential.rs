mod common;

use std::collections::VecDeque;

use common::mock_circuit::{MockCircuit, MockCommand};
use synir::data_structures::{
    Angle, CliffordTableau, PropagateClifford, PauliExponential, PauliPolynomial,
};
use synir::ir::clifford_tableau::{self, CliffordTableauSynthStrategy, NaiveCliffordSynthesizer};
use synir::ir::pauli_exponential::PauliExponentialSynthesizer;
use synir::ir::pauli_polynomial::PauliPolynomialSynthStrategy;
use synir::ir::{AdjointSynthesizer, Synthesizer};

fn setup_simple_pe() -> PauliExponential {
    let ham = vec![("IZZZ", Angle::from_angle(0.3))];

    let pauli_polynomial = PauliPolynomial::from_hamiltonian(ham);
    let clifford_tableau = CliffordTableau::new(4);
    PauliExponential::new(VecDeque::from([pauli_polynomial]), clifford_tableau)
}

fn setup_complex_pe() -> PauliExponential {
    let ham = vec![
        ("IXYZ", Angle::from_angle(0.3)),
        ("XXII", Angle::from_angle(0.7)),
        ("YYII", Angle::from_angle(0.12)),
    ];

    let pauli_polynomial = PauliPolynomial::from_hamiltonian(ham);
    let clifford_tableau = CliffordTableau::new(4);
    PauliExponential::new(VecDeque::from([pauli_polynomial]), clifford_tableau)
}

#[test]
fn test_naive_pauli_exponential_synthesis() {
    let pe = setup_simple_pe();
    let mut mock = MockCircuit::new();
    let mut synthesizer = PauliExponentialSynthesizer::from_strategy(
        PauliPolynomialSynthStrategy::Naive,
        CliffordTableauSynthStrategy::Naive,
    );
    synthesizer.synthesize(pe, &mut mock);
    for c in mock.commands() {
        println!("{:?}", c)
    }
    let final_stab_transform = mock.cliffords_only().to_clifford_tableau(4);
    println!("{}", final_stab_transform);
    assert_eq!(
        final_stab_transform,
        CliffordTableau::new(4)
    );
    let ref_commands = [
        MockCommand::CX(1, 2),
        MockCommand::CX(2, 3),
        MockCommand::Rz(3, 0.3),
        MockCommand::CX(1, 2),
        MockCommand::CX(1, 3),
        MockCommand::CX(2, 3),
    ];

    assert_eq!(mock.commands(), &ref_commands);
}

#[test]
fn test_naive_pauli_exponential_complex() {
    return;
    let pe = setup_complex_pe();
    let mut mock = MockCircuit::new();

    let mut synthesizer = PauliExponentialSynthesizer::from_strategy(
        PauliPolynomialSynthStrategy::Naive,
        CliffordTableauSynthStrategy::Naive,
    );
    synthesizer.synthesize(pe, &mut mock);
    let compute_uncompute_ct = mock.cliffords_only().to_clifford_tableau(4);
    println!("{}", compute_uncompute_ct);
    for c in mock.commands() {
        println!("{:?}", c);
    }
    println!("{}", compute_uncompute_ct);
    assert_eq!(compute_uncompute_ct, CliffordTableau::new(4));

    let input = [
        MockCommand::H(1),
        MockCommand::V(2),
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

    let ref_ct = MockCircuit::from_vec(input.into_iter().collect()).to_clifford_tableau(4);
    let mut mock_ct = MockCircuit::new();

    let mut cliff_synthesizer = NaiveCliffordSynthesizer::default();

    cliff_synthesizer.synthesize_adjoint(ref_ct.clone(), &mut mock_ct);

    let mock_ct_ref_commands = [
        MockCommand::CX(0, 1),
        MockCommand::CX(0, 2),
        MockCommand::CX(0, 3),
        MockCommand::H(1),
        MockCommand::H(2),
        MockCommand::H(3),
        MockCommand::CX(1, 0),
        MockCommand::CX(2, 0),
        MockCommand::CX(3, 0),
        MockCommand::V(1),
        MockCommand::V(3),
        MockCommand::CX(2, 1),
        MockCommand::CX(3, 1),
        MockCommand::CX(3, 2),
        MockCommand::CX(2, 3),
        MockCommand::CX(3, 2),
        MockCommand::CX(3, 2),
        MockCommand::X(1),
        MockCommand::X(2),
    ];

    assert_eq!(mock_ct.commands(), &mock_ct_ref_commands);
    /*
       ("IXYZ", Angle::from_angle(0.3)),
       ("XXII", Angle::from_angle(0.7)),
       ("YYII", Angle::from_angle(0.12)),
    */
    let ref_commands = [
        MockCommand::H(1),
        MockCommand::V(2),
        MockCommand::CX(1, 2),
        MockCommand::CX(2, 3),
        MockCommand::Rz(3, -0.3),
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
        MockCommand::CX(0, 1),
        MockCommand::CX(0, 2),
        MockCommand::CX(0, 3),
        MockCommand::H(1),
        MockCommand::H(2),
        MockCommand::H(3),
        MockCommand::CX(1, 0),
        MockCommand::CX(2, 0),
        MockCommand::CX(3, 0),
        MockCommand::V(1),
        MockCommand::V(3),
        MockCommand::CX(2, 1),
        MockCommand::CX(3, 1),
        MockCommand::CX(3, 2),
        MockCommand::CX(2, 3),
        MockCommand::CX(3, 2),
        MockCommand::CX(3, 2),
        MockCommand::X(1),
        MockCommand::X(2),
    ];

    assert_eq!(mock.commands(), &ref_commands);
}

#[test]
fn test_naive_pe_with_ct(){

    let ham = vec![("IZZZ", Angle::from_angle(0.3))];

    let pauli_polynomial = PauliPolynomial::from_hamiltonian(ham);
    let mut clifford_tableau = CliffordTableau::new(4);
    clifford_tableau.cx(2,3);
    clifford_tableau.cx(1,2);
    clifford_tableau.s(0);
    clifford_tableau.v(0);
    let mut ct_ref = CliffordTableau::new(4);
    // Add PP CNOTs
    ct_ref.cx(1,2);
    ct_ref.cx(2,3);
    // Add CT single qubit gates
    ct_ref.s(0);
    ct_ref.v(0);
    let pe = PauliExponential::new(VecDeque::from([pauli_polynomial]), clifford_tableau);
    let mut mock = MockCircuit::new();
    let mut synthesizer = PauliExponentialSynthesizer::from_strategy(
        PauliPolynomialSynthStrategy::Naive,
        CliffordTableauSynthStrategy::Naive,
    );
    synthesizer.synthesize(pe, &mut mock);
    for c in mock.commands() {
        println!("{:?}", c)
    }
    let final_stab_transform = mock.cliffords_only().to_clifford_tableau(4);
    println!("Final {}", final_stab_transform);
    println!("Ref {}", ct_ref);
    assert_eq!(
        final_stab_transform,
        ct_ref
    );
    let ref_commands = [
        MockCommand::CX(1, 2),
        MockCommand::CX(2, 3),
        MockCommand::Rz(3, 0.3),
        MockCommand::S(0), // TODO Generates H here, which causes a sign error.
        MockCommand::V(0) // Potential sign error in adjoint() or compose().
    ];

    assert_eq!(mock.commands(), &ref_commands);
}