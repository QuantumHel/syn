extern crate rand;

use rand::seq::SliceRandom;

use crate::common::mock_circuit::{
    check_mock_equals_clifford_tableau, parse_clifford_commands, MockCircuit, MockCommand,
};
use crate::common::sample_clifford_tableaus::{
    half_swap_0_1, half_swap_1_0, identity_2qb_ct, sample_2cnot_ladder, sample_cnot_gate,
    sample_cnot_reverse_gate, sample_s_dgr_gate, sample_s_gate, sample_swap_ct, sample_v_dgr_gate,
    sample_v_gate, setup_sample_ct, setup_sample_inverse_ct,
};
use itertools::Itertools;
use synir::data_structures::{CliffordTableau, PropagateClifford};
use synir::ir::clifford_tableau::CallbackCliffordSynthesizer;
use synir::ir::Synthesizer;

fn run_synthesizer(clifford_tableau: &CliffordTableau) -> (MockCircuit, CliffordTableau) {
    let mut mock = MockCircuit::new();
    let mut rng = rand::rng(); //TODO make this from seed
    let custom_columns = (0..clifford_tableau.size()).collect_vec();
    let mut custom_rows = (0..clifford_tableau.size()).collect_vec();
    custom_rows.shuffle(&mut rng);

    let mut synthesizer = CallbackCliffordSynthesizer::custom_pivot(custom_columns, custom_rows);
    let new_ct = synthesizer.synthesize(clifford_tableau.clone(), &mut mock);
    (mock, new_ct)
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

// Single qubit gates
test_clifford!(identity_2qb_ct, None::<&Vec<MockCommand>>);
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

// Advances Clifford Tableau
test_clifford!(setup_sample_ct, None::<&Vec<MockCommand>>);
test_clifford!(setup_sample_inverse_ct, None::<&Vec<MockCommand>>);

// CNOT synthesis
test_clifford!(sample_cnot_gate, None::<&Vec<MockCommand>>);
test_clifford!(sample_cnot_reverse_gate, None::<&Vec<MockCommand>>);
test_clifford!(sample_2cnot_ladder, None::<&Vec<MockCommand>>);
test_clifford!(sample_swap_ct, None::<&Vec<MockCommand>>);
test_clifford!(half_swap_0_1, None::<&Vec<MockCommand>>);
test_clifford!(half_swap_1_0, None::<&Vec<MockCommand>>);

#[test]
fn test_custom_clifford_synthesis_old() {
    let clifford_tableau = setup_sample_ct();
    let mut mock = MockCircuit::new();

    let mut synthesizer = CallbackCliffordSynthesizer::custom_pivot(vec![0, 1, 2], vec![0, 1, 2]);
    synthesizer.synthesize(clifford_tableau.clone(), &mut mock);

    let ref_ct = parse_clifford_commands(3, mock.commands());

    assert_eq!(clifford_tableau, ref_ct);
}

#[test]
fn test_custom_clifford_synthesis_large_old() {
    let clifford_tableau = setup_sample_inverse_ct();
    let mut mock = MockCircuit::new();

    let mut synthesizer =
        CallbackCliffordSynthesizer::custom_pivot(vec![0, 1, 2, 3], vec![0, 2, 1, 3]);

    synthesizer.synthesize(clifford_tableau.clone(), &mut mock);

    let mut ref_ct = parse_clifford_commands(4, mock.commands());
    ref_ct.permute(vec![0, 2, 1, 3]);

    assert_eq!(clifford_tableau, ref_ct);
}

#[test]
fn test_custom_clifford_synthesis_simple_old() {
    let mut clifford_tableau = CliffordTableau::new(3);
    clifford_tableau.cx(0, 1);
    clifford_tableau.cx(1, 2);
    let mut mock = MockCircuit::new();

    let mut synthesizer = CallbackCliffordSynthesizer::custom_pivot(vec![0, 1, 2], vec![0, 1, 2]);
    synthesizer.synthesize(clifford_tableau.clone(), &mut mock);

    let ref_ct = parse_clifford_commands(3, mock.commands());
    assert_eq!(clifford_tableau, ref_ct);
}
