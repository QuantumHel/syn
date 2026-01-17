extern crate rand;

use crate::common::mock_circuit::{check_mock_equals_clifford_tableau, MockCircuit, MockCommand};
use itertools::Itertools;
use synir::data_structures::{Angle, CliffordTableau, PropagateClifford};
use synir::ir::clifford_tableau::CallbackCliffordSynthesizer;
use synir::ir::Synthesizer;

fn run_synthesizer(clifford_tableau: &CliffordTableau) -> (MockCircuit, CliffordTableau) {
    let mut mock = MockCircuit::new();
    // let mut rng = rand::rng(); //TODO make this from seed
    let custom_columns = (0..clifford_tableau.size()).collect_vec();
    let custom_rows = (0..clifford_tableau.size()).collect_vec();
    // custom_rows.shuffle(&mut rng);

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

fn compose_x_gadget() -> CliffordTableau {
    let mut ct = CliffordTableau::new(2);
    let _ = ct.compose_gadget((
        synir::data_structures::PauliString::from_text("XI"),
        Angle::from_pi4_rotation(2),
    ));
    ct
}

fn compose_z_gadget() -> CliffordTableau {
    let mut ct = CliffordTableau::new(2);
    let _ = ct.compose_gadget((
        synir::data_structures::PauliString::from_text("ZI"),
        Angle::from_pi4_rotation(4),
    ));
    ct
}

fn compose_y_gadget() -> CliffordTableau {
    let mut ct = CliffordTableau::new(2);
    let _ = ct.compose_gadget((
        synir::data_structures::PauliString::from_text("YI"),
        Angle::from_pi4_rotation(6),
    ));
    ct
}

test_clifford!(compose_x_gadget, Some(&vec![MockCommand::V(0),]));
test_clifford!(compose_z_gadget, Some(&vec![MockCommand::Z(0),]));

test_clifford!(
    compose_y_gadget,
    Some(&vec![MockCommand::H(0), MockCommand::Z(0)])
);

fn compose_xx_gadget() -> CliffordTableau {
    let mut ct = CliffordTableau::new(2);
    let _ = ct.compose_gadget((
        synir::data_structures::PauliString::from_text("XX"),
        Angle::from_pi4_rotation(2),
    ));
    ct
}

fn compose_zz_gadget() -> CliffordTableau {
    let mut ct = CliffordTableau::new(2);
    let _ = ct.compose_gadget((
        synir::data_structures::PauliString::from_text("ZZ"),
        Angle::from_pi4_rotation(2),
    ));
    ct
}

fn compose_yy_gadget() -> CliffordTableau {
    let mut ct = CliffordTableau::new(2);
    let _ = ct.compose_gadget((
        synir::data_structures::PauliString::from_text("YY"),
        Angle::from_pi4_rotation(2),
    ));
    ct
}

test_clifford!(
    compose_xx_gadget,
    Some(&vec![
        MockCommand::V(0),
        MockCommand::H(1),
        MockCommand::CX(1, 0),
        MockCommand::H(1),
        MockCommand::V(1)
    ])
);

test_clifford!(
    compose_zz_gadget,
    Some(&vec![
        MockCommand::S(0),
        MockCommand::H(1),
        MockCommand::CX(0, 1),
        MockCommand::S(1),
        MockCommand::V(1),
        MockCommand::Z(1)
    ])
);

test_clifford!(
    compose_yy_gadget,
    Some(&vec![
        MockCommand::H(0),
        MockCommand::S(1),
        MockCommand::CX(0, 1),
        MockCommand::H(1),
        MockCommand::CX(1, 0),
        MockCommand::S(1),
        MockCommand::V(1),
        MockCommand::X(1),
        MockCommand::Z(0),
        MockCommand::Z(1)
    ])
);

fn compose_complex_gadget() -> CliffordTableau {
    let mut ct = CliffordTableau::new(3);
    let _ = ct.compose_gadget((
        synir::data_structures::PauliString::from_text("XYZ"),
        Angle::from_pi4_rotation(2),
    ));
    ct
}

fn manual_compose_complex_gadget() -> CliffordTableau {
    let mut ct = CliffordTableau::new(3);
    ct.h(0);
    ct.v(1);
    ct.cx(0, 1);
    ct.cx(1, 2);
    ct.s(2);
    ct.cx(1, 2);
    ct.cx(0, 1);
    ct.h(0);
    ct.v_dgr(1);
    ct
}

test_clifford!(
    compose_complex_gadget,
    Some(&vec![
        MockCommand::V(0),
        MockCommand::V(1),
        MockCommand::CX(1, 0),
        MockCommand::CX(2, 0),
        MockCommand::S(1),
        MockCommand::H(2),
        MockCommand::CX(1, 2),
        MockCommand::V(1),
        MockCommand::S(2),
        MockCommand::V(2),
        MockCommand::X(1),
        MockCommand::Z(2)
    ])
);

#[test]
fn test_correctness() {
    let ct1 = compose_complex_gadget();
    let ct2 = manual_compose_complex_gadget();
    assert_eq!(ct1, ct2);
}
