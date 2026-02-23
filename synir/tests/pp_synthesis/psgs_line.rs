use std::collections::VecDeque;

use crate::common::mock_circuit::{MockCircuit, MockCommand};
use crate::common::sample_pauli_poly::*;
use synir::architecture::connectivity::Connectivity;
use synir::data_structures::{CliffordTableau, PauliPolynomial};
use synir::ir::pauli_polynomial::psgs::PSGSPauliPolynomialSynthesizer;
use synir::ir::Synthesizer;

fn run_synthesizer(pp: VecDeque<PauliPolynomial>) -> (MockCircuit, CliffordTableau) {
    let mut mock = MockCircuit::new();
    let mut synthesizer = PSGSPauliPolynomialSynthesizer::default();
    synthesizer.set_connectivity(Connectivity::line(pp[0].num_qubits()));
    let ct = synthesizer.synthesize(pp, &mut mock);
    assert!(mock.fits_connectivity(synthesizer.get_connectivity()));
    return (mock, ct);
}

macro_rules! test_pp {
    ($fun:ident, $expected:expr) => {
        paste::item! {
            #[test]
            fn [< synthesize_ $fun>]() {
                let pp = $fun();
                let ref_pp_mock = $expected;
                let ref_ct_mock = ref_pp_mock.cliffords_only();
                let (mock, new_ct) = run_synthesizer(pp);
                println!("Synthesized:");
                for c in mock.commands() {
                    println!("{:?}", c);
                }
                assert_eq!(mock, ref_pp_mock);
                assert!(ref_ct_mock.equals_clifford_tableau(&new_ct, new_ct.get_permutation()));
            }
        }
    };
}

test_pp!(
    setup_simple_pp,
    MockCircuit::from_vec(
        [
            MockCommand::CX(3, 2),
            MockCommand::VDgr(2),
            MockCommand::SDgr(1),
            MockCommand::CX(2, 1),
            MockCommand::Ry(1, 0.3),
        ]
        .into_iter()
        .collect()
    )
);

test_pp!(
    setup_complex_pp,
    MockCircuit::from_vec(
        [
            MockCommand::CX(3, 2),
            MockCommand::CX(2, 1),
            MockCommand::Rz(1, 0.3),
            MockCommand::H(1),
            MockCommand::SDgr(0),
            MockCommand::CX(1, 0),
            MockCommand::Ry(0, 0.7),
        ]
        .into_iter()
        .collect()
    )
);

test_pp!(
    setup_parallel_pp,
    MockCircuit::from_vec(vec![
        MockCommand::CX(1, 0), 
        MockCommand::Rx(1, 0.234234),
        MockCommand::Ry(0, 0.53423),
    ])
);
