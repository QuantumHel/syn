mod common;

use bitvec::bitvec;
use bitvec::prelude::Lsb0;
use common::{parse_clifford_commands, MockCircuit, MockCommand};
use syn::architecture::connectivity::Connectivity;
use syn::architecture::Architecture;
use syn::data_structures::{CliffordTableau, PauliString, PropagateClifford};
use syn::ir::clifford_tableau::{
    CallbackCliffordSynthesizer, NaiveCliffordSynthesizer, PermRowColCliffordSynthesizer,
};
use syn::ir::{AdjointSynthesizer, CliffordGates, Synthesizer};

fn setup_sample_ct() -> CliffordTableau {
    // Stab: ZZZ, -YIY, XIX
    // Destab: -IXI, XXI, IYY
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
    CliffordTableau::from_parts(vec![pauli_1, pauli_2, pauli_3], signs)
}

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

fn setup_2_qubit_clifford() -> CliffordTableau {
    // qubit 1x: ZZXI
    // qubit 1z: YYZX
    let pauli_1 = PauliString::from_text("XIZI");

    // qubit 2x: IZZZ
    // qubit 2z: YYIZ
    let pauli_2 = PauliString::from_text("IXIZ");

    let signs = bitvec![0, 0, 0, 0];
    CliffordTableau::from_parts(vec![pauli_1, pauli_2], signs)
}

#[test]
fn test_id_synthesis() {
    let clifford_tableau = setup_2_qubit_clifford();
    let mut mock = MockCircuit::new();

    let mut synthesizer = NaiveCliffordSynthesizer::default();
    synthesizer.synthesize(clifford_tableau, &mut mock);
    assert_eq!(mock.commands(), &vec![]);
}

#[test]
fn test_s_synthesis() {
    let mut clifford_tableau = setup_2_qubit_clifford();
    clifford_tableau.s(1);
    let mut mock = MockCircuit::new();
    let mut synthesizer = NaiveCliffordSynthesizer::default();
    synthesizer.synthesize(clifford_tableau.clone(), &mut mock);

    assert_eq!(mock.commands(), &vec![MockCommand::S(1)]);
}

#[test]
fn test_s_adjoint_synthesis() {
    let mut clifford_tableau = setup_2_qubit_clifford();
    clifford_tableau.s(1);
    let mut mock = MockCircuit::new();
    let mut synthesizer = NaiveCliffordSynthesizer::default();
    synthesizer.synthesize_adjoint(clifford_tableau.clone(), &mut mock);

    let ref_ct = parse_clifford_commands(2, mock.commands());
    assert_eq!(clifford_tableau * ref_ct, CliffordTableau::new(2));

    assert_eq!(mock.commands(), &vec![MockCommand::S(1), MockCommand::Z(1)]);
}

#[test]
fn test_v_synthesis() {
    let mut clifford_tableau = setup_2_qubit_clifford();
    clifford_tableau.v(1);
    let mut mock = MockCircuit::new();
    let mut synthesizer = NaiveCliffordSynthesizer::default();
    synthesizer.synthesize(clifford_tableau.clone(), &mut mock);

    assert_eq!(
        mock.commands(),
        &vec![
            MockCommand::S(1),
            MockCommand::H(1),
            MockCommand::S(1),
            MockCommand::X(1)
        ]
    );
}

#[test]
fn test_v_adjoint_synthesis() {
    let mut clifford_tableau = setup_2_qubit_clifford();
    clifford_tableau.v(1);
    let mut mock = MockCircuit::new();
    let mut synthesizer = NaiveCliffordSynthesizer::default();
    synthesizer.synthesize_adjoint(clifford_tableau.clone(), &mut mock);

    let ref_ct = parse_clifford_commands(2, mock.commands());
    assert_eq!(clifford_tableau * ref_ct, CliffordTableau::new(2));

    assert_eq!(
        mock.commands(),
        &vec![MockCommand::S(1), MockCommand::H(1), MockCommand::S(1),]
    );
}

#[test]
fn test_cnot_synthesis() {
    let mut clifford_tableau = setup_2_qubit_clifford();
    clifford_tableau.cx(0, 1);
    let mut mock = MockCircuit::new();

    let mut synthesizer = NaiveCliffordSynthesizer::default();
    synthesizer.synthesize(clifford_tableau.clone(), &mut mock);

    assert_eq!(mock.commands(), &vec![MockCommand::CX(0, 1)]);
}

#[test]
fn test_cnot_reverse_synthesis() {
    let mut clifford_tableau = setup_2_qubit_clifford();
    clifford_tableau.cx(1, 0);
    let mut mock = MockCircuit::new();

    let mut synthesizer = NaiveCliffordSynthesizer::default();
    synthesizer.synthesize(clifford_tableau.clone(), &mut mock);

    assert_eq!(mock.commands(), &vec![MockCommand::CX(1, 0)]);
}

#[test]
fn test_clifford_synthesis() {
    let clifford_tableau = setup_sample_ct();
    let mut mock = MockCircuit::new();

    let mut synthesizer = NaiveCliffordSynthesizer::default();
    synthesizer.synthesize(clifford_tableau.clone(), &mut mock);

    let ref_ct = parse_clifford_commands(3, mock.commands());

    assert_eq!(clifford_tableau, ref_ct);
}

#[test]
fn test_clifford_synthesis_large() {
    let clifford_tableau = setup_sample_inverse_ct();
    let mut mock = MockCircuit::new();

    let mut synthesizer = NaiveCliffordSynthesizer::default();
    synthesizer.synthesize(clifford_tableau.clone(), &mut mock);

    let ref_ct = parse_clifford_commands(4, mock.commands());

    assert_eq!(clifford_tableau, ref_ct);
}

#[test]
fn test_clifford_synthesis_simple() {
    let mut clifford_tableau = CliffordTableau::new(3);
    clifford_tableau.cx(0, 1);
    clifford_tableau.cx(1, 2);
    let mut mock = MockCircuit::new();

    let mut synthesizer = NaiveCliffordSynthesizer::default();
    synthesizer.synthesize(clifford_tableau.clone(), &mut mock);

    let ref_ct = parse_clifford_commands(3, mock.commands());
    assert_eq!(clifford_tableau, ref_ct);
}

#[test]
fn test_custom_clifford_synthesis() {
    let clifford_tableau = setup_sample_ct();
    let mut mock = MockCircuit::new();

    let mut synthesizer = CallbackCliffordSynthesizer::custom_pivot(vec![0, 1, 2], vec![0, 1, 2]);
    synthesizer.synthesize(clifford_tableau.clone(), &mut mock);

    let ref_ct = parse_clifford_commands(3, mock.commands());

    assert_eq!(clifford_tableau, ref_ct);
}

#[test]
fn test_custom_clifford_synthesis_large() {
    let clifford_tableau = setup_sample_inverse_ct();
    let mut mock = MockCircuit::new();

    let mut synthesizer =
        CallbackCliffordSynthesizer::custom_pivot(vec![0, 1, 2, 3], vec![0, 2, 1, 3]);

    synthesizer.synthesize(clifford_tableau.clone(), &mut mock);

    let mut ref_ct = parse_clifford_commands(4, mock.commands());
    ref_ct.permute(&[0, 2, 1, 3]);

    assert_eq!(clifford_tableau, ref_ct);
}

#[test]
fn test_custom_clifford_synthesis_simple() {
    let mut clifford_tableau = CliffordTableau::new(3);
    clifford_tableau.cx(0, 1);
    clifford_tableau.cx(1, 2);
    let mut mock = MockCircuit::new();

    let mut synthesizer = CallbackCliffordSynthesizer::custom_pivot(vec![0, 1, 2], vec![0, 1, 2]);
    synthesizer.synthesize(clifford_tableau.clone(), &mut mock);

    let ref_ct = parse_clifford_commands(3, mock.commands());
    assert_eq!(clifford_tableau, ref_ct);
}

#[test]
fn test_prc_clifford_synthesis() {
    let clifford_tableau = setup_sample_ct();
    let num_qubits = clifford_tableau.size();
    let mut mock = MockCircuit::new();
    let connectivity = Connectivity::complete(num_qubits);
    let mut synthesizer = PermRowColCliffordSynthesizer::new(connectivity);
    synthesizer.synthesize(clifford_tableau.clone(), &mut mock);

    let mut ref_ct = parse_clifford_commands(3, mock.commands());
    ref_ct.permute(synthesizer.permutation());
    println!("ref_ct: {}", ref_ct);
    println!("clifford_tableau: {}", clifford_tableau);
    assert_eq!(clifford_tableau, ref_ct);
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
