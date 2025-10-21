use crate::common::mock_circuit::{check_mock_equals_clifford_tableau, MockCircuit, MockCommand};
use crate::common::sample_clifford_tableaus::{
    half_swap_0_1, half_swap_1_0, identity_2qb_ct, sample_2cnot_ladder, sample_cnot_gate,
    sample_cnot_reverse_gate, sample_s_dgr_gate, sample_s_gate, sample_swap_ct, sample_v_dgr_gate,
    sample_v_gate, setup_sample_ct, setup_sample_inverse_ct,
};
use synir::data_structures::CliffordTableau;
use synir::ir::clifford_tableau::NaiveCliffordSynthesizer;
use synir::ir::Synthesizer;

fn run_synthesizer(clifford_tableau: &CliffordTableau) -> (MockCircuit, CliffordTableau) {
    let mut mock = MockCircuit::new();
    let mut synthesizer = NaiveCliffordSynthesizer::default();
    let new_ct = synthesizer.synthesize(clifford_tableau.clone(), &mut mock);
    return (mock, new_ct);
}

macro_rules! test_clifford {
    ($fun:ident, $expected:expr) => {
        paste::item! {
                #[test]
                fn [< synthesize_ $fun>]() {
                    let clifford_tableau = $fun();
                    let (mock, new_ct) = run_synthesizer(&clifford_tableau);
                    if $expected.is_some() {
                        assert_eq!(mock.commands(), $expected.unwrap());
                    }
                    check_mock_equals_clifford_tableau(&clifford_tableau, &mock, new_ct.get_permutation());
                }
            }
    };
}

test_clifford!(identity_2qb_ct, Some::<&Vec<MockCommand>>(&vec![]));
test_clifford!(sample_s_gate, Some(&vec![MockCommand::S(0)]));
test_clifford!(
    sample_s_dgr_gate,
    Some(&vec![MockCommand::S(0), MockCommand::Z(0)])
);
test_clifford!(sample_v_gate, Some(&vec![MockCommand::V(0)]));
test_clifford!(
    sample_v_dgr_gate,
    Some(&vec![MockCommand::V(0), MockCommand::X(0)])
);
test_clifford!(sample_cnot_gate, Some(&vec![MockCommand::CX(0, 1)]));
test_clifford!(sample_cnot_reverse_gate, Some(&vec![MockCommand::CX(1, 0)]));
test_clifford!(setup_sample_ct, None::<&Vec<MockCommand>>);
test_clifford!(setup_sample_inverse_ct, None::<&Vec<MockCommand>>);
test_clifford!(
    sample_2cnot_ladder,
    Some(&vec![MockCommand::CX(0, 1), MockCommand::CX(1, 2)])
);
test_clifford!(
    sample_swap_ct,
    Some(&vec![
        MockCommand::CX(1, 0),
        MockCommand::CX(0, 1),
        MockCommand::CX(1, 0)
    ])
);
test_clifford!(
    half_swap_0_1,
    Some(&vec![
        MockCommand::CX(1, 0),
        MockCommand::CX(0, 1),
        MockCommand::CX(1, 0),
        MockCommand::CX(0, 1)
    ])
);
test_clifford!(
    half_swap_1_0,
    Some(&vec![
        MockCommand::CX(1, 0),
        MockCommand::CX(0, 1),
        MockCommand::CX(1, 0),
        MockCommand::CX(1, 0)
    ])
);
