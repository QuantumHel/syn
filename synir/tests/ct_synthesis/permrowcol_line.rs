use crate::common::mock_circuit::{check_mock_equals_clifford_tableau, MockCircuit, MockCommand};
use crate::common::sample_clifford_tableaus::{
    half_swap_0_1, half_swap_1_0, sample_2cnot_ladder, sample_cnot_gate, sample_cnot_reverse_gate,
    sample_s_dgr_gate, sample_s_gate, sample_swap_ct, sample_v_dgr_gate, sample_v_gate,
    setup_sample_ct, setup_sample_inverse_ct,
};
use synir::architecture::connectivity::Connectivity;
use synir::data_structures::CliffordTableau;
use synir::ir::clifford_tableau::PermRowColCliffordSynthesizer;
use synir::ir::Synthesizer;

fn run_synthesizer(clifford_tableau: &CliffordTableau) -> (MockCircuit, CliffordTableau) {
    let num_qubits = clifford_tableau.size();
    let mut mock = MockCircuit::new();
    let connectivity = Connectivity::line(num_qubits);
    println!("size {}, qpu {}", num_qubits, connectivity.node_count());
    let mut synthesizer = PermRowColCliffordSynthesizer::new(connectivity);
    let new_ct = synthesizer.synthesize(clifford_tableau.clone(), &mut mock);
    return (mock, new_ct);
}

#[test]
fn test_id_synthesis() {
    let clifford_tableau = CliffordTableau::new(2);
    let (mock, new_ct) = run_synthesizer(&clifford_tableau);
    check_mock_equals_clifford_tableau(&clifford_tableau, &mock, new_ct.get_permutation());
    assert_eq!(mock.commands(), &vec![]);
}

#[test]
fn test_s_synthesis() {
    let clifford_tableau = sample_s_gate();
    let (mock, new_ct) = run_synthesizer(&clifford_tableau);
    check_mock_equals_clifford_tableau(&clifford_tableau, &mock, new_ct.get_permutation());
    assert_eq!(mock.commands(), &vec![MockCommand::S(0)]);
}

#[test]
fn test_s_dgr_synthesis() {
    let clifford_tableau = sample_s_dgr_gate();
    let (mock, new_ct) = run_synthesizer(&clifford_tableau);
    check_mock_equals_clifford_tableau(&clifford_tableau, &mock, new_ct.get_permutation());
    assert_eq!(mock.commands(), &vec![MockCommand::S(0), MockCommand::Z(0)]);
}

#[test]
fn test_v_synthesis() {
    let clifford_tableau = sample_v_gate();
    let (mock, new_ct) = run_synthesizer(&clifford_tableau);
    check_mock_equals_clifford_tableau(&clifford_tableau, &mock, new_ct.get_permutation());
    assert_eq!(mock.commands(), &vec![MockCommand::V(0),]);
}

#[test]
fn test_v_dgr_synthesis() {
    let clifford_tableau = sample_v_dgr_gate();
    let (mock, new_ct) = run_synthesizer(&clifford_tableau);
    check_mock_equals_clifford_tableau(&clifford_tableau, &mock, new_ct.get_permutation());
    assert_eq!(mock.commands(), &vec![MockCommand::V(0), MockCommand::X(0)]);
}

#[test]
fn test_cnot_synthesis() {
    let clifford_tableau = sample_cnot_gate();
    let (mock, new_ct) = run_synthesizer(&clifford_tableau);
    check_mock_equals_clifford_tableau(&clifford_tableau, &mock, new_ct.get_permutation());
    assert_eq!(mock.commands(), &vec![MockCommand::CX(0, 1)]);
}

#[test]
fn test_cnot_reverse_synthesis() {
    let clifford_tableau = sample_cnot_reverse_gate();
    let (mock, new_ct) = run_synthesizer(&clifford_tableau);
    check_mock_equals_clifford_tableau(&clifford_tableau, &mock, new_ct.get_permutation());
    assert_eq!(mock.commands(), &vec![MockCommand::CX(1, 0)]);
}

#[test]
fn test_clifford_synthesis() {
    let clifford_tableau = setup_sample_ct();
    let (mock, new_ct) = run_synthesizer(&clifford_tableau);
    check_mock_equals_clifford_tableau(&clifford_tableau, &mock, new_ct.get_permutation());
}

#[test]
fn test_clifford_synthesis_large() {
    let clifford_tableau = setup_sample_inverse_ct();
    let (mock, new_ct) = run_synthesizer(&clifford_tableau);
    check_mock_equals_clifford_tableau(&clifford_tableau, &mock, new_ct.get_permutation());
}

#[test]
fn test_clifford_synthesis_simple() {
    let clifford_tableau = sample_2cnot_ladder();
    let (mock, new_ct) = run_synthesizer(&clifford_tableau);
    check_mock_equals_clifford_tableau(&clifford_tableau, &mock, new_ct.get_permutation());
    assert_eq!(
        mock.commands(),
        &vec![MockCommand::CX(0, 1), MockCommand::CX(1, 2)]
    );
}

#[test]
fn test_prc_swap_to_identity() {
    let clifford_tableau = sample_swap_ct();
    let (mock, new_ct) = run_synthesizer(&clifford_tableau);
    check_mock_equals_clifford_tableau(&clifford_tableau, &mock, new_ct.get_permutation());
    assert_eq!(mock.commands().len(), 0);
}

#[test]
fn test_prc_2_cnot_reduction_v1() {
    let clifford_tableau = half_swap_0_1();
    let (mock, new_ct) = run_synthesizer(&clifford_tableau);
    check_mock_equals_clifford_tableau(&clifford_tableau, &mock, new_ct.get_permutation());
    assert_eq!(mock.commands().len(), 1);
}

#[test]
fn test_prc_2_cnot_reduction_v2() {
    let clifford_tableau = half_swap_1_0();
    let (mock, new_ct) = run_synthesizer(&clifford_tableau);
    check_mock_equals_clifford_tableau(&clifford_tableau, &mock, new_ct.get_permutation());
    assert_eq!(mock.commands().len(), 1);
}
