mod common;

use bitvec::bitvec;
use bitvec::prelude::Lsb0;
use common::{parse_clifford_commands, MockCircuit};
use synir::architecture::connectivity::Connectivity;
use synir::data_structures::{CliffordTableau, PauliString, PropagateClifford};
use synir::ir::clifford_tableau::PermRowColCliffordSynthesizer;
use synir::ir::Synthesizer;

fn setup_sample_inverse_ct() -> CliffordTableau {
    // Stab: -ZIYZ, -ZZYZ, -XZXI, IZXX
    // Destab: -YYIZ, -YYXZ, ZIXX, -XZXZ
    // qubit 1x: ZZXI
    // qubit 1z: YYZX
    let pauli_1 = PauliString::from_text("ZZXIYYZX");

    // qubit 2x: IZZZ
    // qubit 2z: YYIZ
    let pauli_2 = PauliString::from_text("IZZZYYIZ");

    // qubit 3x: YYXX
    // qubit 3z: IXXX
    let pauli_3 = PauliString::from_text("YYXXIXXX");

    // qubit 3x: ZZIX
    // qubit 3z: ZZXZ
    let pauli_4 = PauliString::from_text("ZZIXZZXZ");

    let signs = bitvec![1, 1, 1, 0, 1, 1, 0, 1];
    CliffordTableau::from_parts(vec![pauli_1, pauli_2, pauli_3, pauli_4], signs)
}

#[test]
fn test_prc_clifford_synthesis_large() {
    let mut clifford_tableau = setup_sample_inverse_ct();
    let mut mock = MockCircuit::new();

    let connectivity = Connectivity::grid(2, 2);
    let mut synthesizer = PermRowColCliffordSynthesizer::new(connectivity);

    synthesizer.synthesize(clifford_tableau.clone(), &mut mock);

    let ref_ct = parse_clifford_commands(4, mock.commands());
    clifford_tableau.permute(synthesizer.permutation());

    assert_eq!(clifford_tableau, ref_ct);
}

#[test]
fn test_prc_clifford_synthesis_simple() {
    let num_qubits = 3;
    let mut clifford_tableau = CliffordTableau::new(num_qubits);

    clifford_tableau.cx(2, 1);
    clifford_tableau.cx(1, 2);
    clifford_tableau.cx(0, 2);
    let mut mock = MockCircuit::new();

    let connectivity = Connectivity::line(num_qubits);

    let mut synthesizer = PermRowColCliffordSynthesizer::new(connectivity);
    synthesizer.synthesize(clifford_tableau.clone(), &mut mock);

    let ref_ct = parse_clifford_commands(3, mock.commands());

    clifford_tableau.permute(synthesizer.permutation());
    assert_eq!(clifford_tableau, ref_ct);
}

#[test]
fn test_prc_swap_to_identity() {
    let num_qubits = 2;
    let mut clifford_tableau = CliffordTableau::new(num_qubits);

    clifford_tableau.cx(0, 1);
    clifford_tableau.cx(1, 0);
    clifford_tableau.cx(0, 1);
    println!("clifford_tableau: {}", clifford_tableau);
    let mut mock = MockCircuit::new();

    let connectivity = Connectivity::line(num_qubits);

    let mut synthesizer = PermRowColCliffordSynthesizer::new(connectivity);
    synthesizer.synthesize(clifford_tableau.clone(), &mut mock);

    let ref_ct = parse_clifford_commands(2, mock.commands());

    clifford_tableau.permute(synthesizer.permutation());
    // Check that the synthesized circuit and original are the same
    assert_eq!(clifford_tableau, ref_ct);
    // Check that the resulting circuit is empty
    assert_eq!(mock.commands().len(), 0);
}
