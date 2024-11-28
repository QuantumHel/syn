// mod common;

// use common::{MockCircuit, MockCommand};
// use syn::data_structures::{CliffordTableau, PauliPolynomial};
// use syn::ir::pauli_exponential::{PauliExponential, PauliExponentialSynthesizer};
// use syn::synthesis_methods::naive::Naive;

// fn setup_simple_pe() -> PauliExponential {
//     let ham = vec![("IZZZ", 0.3)];

//     let pauli_polynomial = PauliPolynomial::from_hamiltonian(ham);
//     let clifford_tableau = CliffordTableau::new(4);
//     PauliExponential::new(pauli_polynomial, clifford_tableau)
// }

// fn setup_complex_pe() -> PauliExponential {
//     let ham = vec![("IXYZ", 0.3), ("XXII", 0.7), ("YYII", 0.12)];

//     let pauli_polynomial = PauliPolynomial::from_hamiltonian(ham);
//     let clifford_tableau = CliffordTableau::new(4);
//     PauliExponential::new(pauli_polynomial, clifford_tableau)
// }

// #[test]
// fn test_naive_clifford_synthesis() {
//     let pe = setup_simple_pe();
//     let mut mock = MockCircuit::new();
//     PauliExponentialSynthesizer::run(pe.clone(), &mut mock);

//     let ref_commands = [
//         MockCommand::CX(1, 2),
//         MockCommand::CX(2, 3),
//         MockCommand::Rz(3, 0.3),
//         MockCommand::CX(1, 2),
//         MockCommand::CX(1, 3),
//         MockCommand::CX(2, 3),
//     ];

//     assert_eq!(mock.commands(), &ref_commands);
// }

// #[test]
// fn test_naive_clifford_synthesis_complex() {
//     let pe = setup_complex_pe();
//     let mut mock = MockCircuit::new();
//     PauliExponentialSynthesizer::run(pe.clone(), &mut mock);

//     let ref_commands = [
//         MockCommand::CX(1, 2),
//         MockCommand::CX(2, 3),
//         MockCommand::Rz(3, 0.3),
//         MockCommand::CX(1, 2),
//         MockCommand::CX(1, 3),
//         MockCommand::CX(2, 3),
//     ];

//     assert_eq!(mock.commands(), &ref_commands);
// }
