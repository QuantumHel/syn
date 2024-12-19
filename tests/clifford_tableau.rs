mod common;

use bitvec::bitvec;
use bitvec::prelude::Lsb0;
use common::{parse_clifford_commands, MockCircuit, MockCommand};
use syn::data_structures::{CliffordTableau, PauliString, PropagateClifford};
use syn::ir::clifford_tableau::CliffordTableauSynthesizer;
use syn::synthesis_methods::{custom::Custom, naive::Naive};

fn setup_sample_ct() -> CliffordTableau {
    // Stab: ZZZ, -YIY, XIX
    // Destab: -IXI, XXI, IYY
    let ct_size = 3;
    // qubit 1x: ZYI
    // qubit 1z: IZZ
    let pauli_1 = PauliString::from_text("ZYIIZZ");

    // qubit 2x: ZIX
    // qubit 2z: XII
    let pauli_2 = PauliString::from_text("ZIXXII");

    // qubit 3x: ZYY
    // qubit 3z: IIZ
    let pauli_3 = PauliString::from_text("ZYYIIZ");

    let signs = bitvec![0, 1, 0, 1, 0, 0];
    CliffordTableau::from_parts(vec![pauli_1, pauli_2, pauli_3], signs, ct_size)
}

fn setup_sample_inverse_ct() -> CliffordTableau {
    // Stab: -ZIYZ, -ZZYZ, -XZXI, IZXX
    // Destab: -YYIZ, -YYXZ, ZIXX, -XZXZ
    let ct_size = 4;
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
    CliffordTableau::from_parts(vec![pauli_1, pauli_2, pauli_3, pauli_4], signs, ct_size)
}

#[test]
fn test_clifford_synthesis() {
    let clifford = setup_sample_ct();
    let mut mock = MockCircuit::new();
    CliffordTableauSynthesizer::run_naive(&clifford, &mut mock);

    let ref_ct = parse_clifford_commands(3, mock.commands());
    assert_eq!(clifford, ref_ct);
}

#[test]
fn test_clifford_synthesis_large() {
    let clifford = setup_sample_inverse_ct();
    let mut mock = MockCircuit::new();
    CliffordTableauSynthesizer::run_naive(&clifford, &mut mock);

    let ref_ct = parse_clifford_commands(4, mock.commands());

    assert_eq!(clifford, ref_ct);
}

#[test]
fn test_clifford_synthesis_simple() {
    let mut clifford = CliffordTableau::new(3);
    clifford.cx(0, 1);
    clifford.cx(1, 2);
    let mut mock = MockCircuit::new();
    CliffordTableauSynthesizer::run_naive(&clifford, &mut mock);

    let ref_ct = parse_clifford_commands(3, mock.commands());
    assert_eq!(clifford, ref_ct);
}

#[test]
fn test_custom_clifford_synthesis() {
    let clifford = setup_sample_ct();
    let mut mock = MockCircuit::new();
    CliffordTableauSynthesizer::run_custom(&clifford, &mut mock, vec![0, 1, 2], vec![0, 1, 2]);

    let ref_ct = parse_clifford_commands(3, mock.commands());

    assert_eq!(clifford, ref_ct);
}

#[test]
fn test_custom_clifford_synthesis_large() {
    let clifford = setup_sample_inverse_ct();
    let mut mock = MockCircuit::new();
    CliffordTableauSynthesizer::run_custom(
        &clifford,
        &mut mock,
        vec![0, 1, 2, 3],
        vec![0, 2, 1, 3],
    );

    let mut ref_ct = parse_clifford_commands(4, mock.commands());
    ref_ct.permute(vec![0, 2, 1, 3]);

    assert_eq!(clifford, ref_ct);
}

#[test]
fn test_custom_clifford_synthesis_simple() {
    let mut clifford = CliffordTableau::new(3);
    clifford.cx(0, 1);
    clifford.cx(1, 2);
    let mut mock = MockCircuit::new();
    CliffordTableauSynthesizer::run_custom(&clifford, &mut mock, vec![0, 1, 2], vec![0, 1, 2]);
    let ref_ct = parse_clifford_commands(3, mock.commands());
    assert_eq!(clifford, ref_ct);
}
