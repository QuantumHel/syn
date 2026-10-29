//! End-to-end `PauliExponential` round trips **on a restricted coupling graph**.
//!
//! [`pe_synthesis_tests`](../pe_synthesis_tests.rs) checks that synthesis preserves the unitary;
//! this checks that it does so without ever placing a two-qubit gate off the hardware graph, which
//! `Connectivity::complete` can never catch.
//!
//! The gap this closes: `pp_synthesis/psgs_line.rs` covers the polynomial half alone and
//! `ct_synthesis/permrowcol_line.rs` the Clifford half alone. Nothing ran *both* through
//! [`PauliExponentialSynthesizer`] against the same non-complete graph, which is exactly what
//! `synthesize_to_qiskit` now does when it is handed a connectivity.

mod common;

use common::mock_circuit::{MockCircuit, MockCommand};
use common::unitary::TOL;

use synir::architecture::connectivity::Connectivity;
use synir::import::{import_circuit, ImportedGate};
use synir::ir::{
    clifford_tableau::CliffordTableauSynthStrategy, pauli_exponential::PauliExponentialSynthesizer,
    pauli_polynomial::PauliPolynomialSynthStrategy, Synthesizer,
};

fn to_imported(command: &MockCommand) -> ImportedGate {
    match command {
        MockCommand::H(t) => ImportedGate::H(*t),
        MockCommand::S(t) => ImportedGate::S(*t),
        MockCommand::SDgr(t) => ImportedGate::SDgr(*t),
        MockCommand::V(t) => ImportedGate::V(*t),
        MockCommand::VDgr(t) => ImportedGate::VDgr(*t),
        MockCommand::X(t) => ImportedGate::X(*t),
        MockCommand::Y(t) => ImportedGate::Y(*t),
        MockCommand::Z(t) => ImportedGate::Z(*t),
        MockCommand::CX(c, t) => ImportedGate::CX(*c, *t),
        MockCommand::Rx(t, angle) => ImportedGate::Rx(*t, *angle),
        MockCommand::Rz(t, angle) => ImportedGate::Rz(*t, *angle),
        other => panic!("{other:?} has no ImportedGate equivalent; rewrite the reference circuit"),
    }
}

/// The one strategy pair that makes the *whole* circuit conform.
///
/// Conformance is a property of the emitted circuit, so it needs both halves: `PSGS` for the
/// gadgets and `PermRowCol` for the trailing Clifford. Either `Naive` leaves its half routing
/// wherever it likes, which `naive_half_leaves_gates_off_the_graph` below pins down.
fn conforming_strategies() -> Vec<(
    &'static str,
    PauliPolynomialSynthStrategy,
    CliffordTableauSynthStrategy,
)> {
    vec![(
        "PSGS/PermRowCol",
        PauliPolynomialSynthStrategy::PSGS,
        CliffordTableauSynthStrategy::PermRowCol,
    )]
}

/// Synthesizes `reference` against `connectivity` and asserts both properties that matter:
/// the unitary is preserved, and every two-qubit gate sits on an edge.
///
/// `expect_conforming` is false only for the strategies documented as ignoring the graph, where
/// the point of the test is to pin that behaviour rather than to check it.
fn assert_round_trip_on(
    name: &str,
    num_qubits: usize,
    reference: &[MockCommand],
    connectivity: &Connectivity,
    strategy_name: &str,
    pauli_strategy: PauliPolynomialSynthStrategy,
    tableau_strategy: CliffordTableauSynthStrategy,
    expect_conforming: bool,
) {
    let gates: Vec<ImportedGate> = reference.iter().map(to_imported).collect();
    let pe = import_circuit(num_qubits, &gates);

    let mut mock = MockCircuit::new();
    let mut synthesizer =
        PauliExponentialSynthesizer::from_strategy(pauli_strategy, tableau_strategy);
    synthesizer.set_connectivity(connectivity.clone());
    synthesizer.synthesize(pe, &mut mock);

    if expect_conforming {
        let violations = mock.connectivity_violations(connectivity);
        assert!(
            violations.is_empty(),
            "\n{name} [{strategy_name}]: emitted two-qubit gates off the coupling graph: \
             {violations:?}\ngraph edges: {:?}\nfull output: {:?}\n",
            connectivity.edges(),
            mock.commands()
        );
    }

    let permutation = mock.final_permutation().cloned();
    let expected = MockCircuit::from_vec(reference.to_vec()).to_unitary(num_qubits);
    let got = mock.to_unitary_with_permutation(num_qubits, permutation.clone());
    assert!(
        got.equiv_up_to_global_phase(&expected, TOL),
        "\n{name} [{strategy_name}]: round trip does not preserve the unitary.\n\
         input:  {reference:?}\noutput: {:?}\npermutation: {permutation:?}\n",
        mock.commands()
    );
}

/// Declares one test per conforming strategy pair, reporting every failing pair at once.
macro_rules! connectivity_tests {
    ($name:ident, $num_qubits:expr, $connectivity:expr, $commands:expr) => {
        paste::item! {
            #[test]
            fn [< fits_ $name >]() {
                let commands: Vec<MockCommand> = $commands;
                let connectivity: Connectivity = $connectivity;
                let mut failures = Vec::new();
                for (strategy_name, pauli, tableau) in conforming_strategies() {
                    let result = std::panic::catch_unwind(|| {
                        assert_round_trip_on(
                            stringify!($name),
                            $num_qubits,
                            &commands,
                            &connectivity,
                            strategy_name,
                            pauli,
                            tableau,
                            true,
                        );
                    });
                    if result.is_err() {
                        failures.push(strategy_name);
                    }
                }
                assert!(
                    failures.is_empty(),
                    "{}: failed for strategies {:?} (details above)",
                    stringify!($name),
                    failures
                );
            }
        }
    };
}

// ---------------------------------------------------------------------------
// A line: the harshest common topology, and the one where every routing decision shows up
// ---------------------------------------------------------------------------

connectivity_tests!(
    line_bell,
    2,
    Connectivity::line(2),
    vec![MockCommand::H(0), MockCommand::CX(0, 1)]
);

connectivity_tests!(
    line_pauli_gadget,
    3,
    Connectivity::line(3),
    vec![
        MockCommand::CX(0, 1),
        MockCommand::CX(1, 2),
        MockCommand::Rz(2, std::f64::consts::FRAC_PI_4),
        MockCommand::CX(1, 2),
        MockCommand::CX(0, 1),
    ]
);

connectivity_tests!(
    line_pauli_gadget_half,
    3,
    Connectivity::line(3),
    vec![
        MockCommand::CX(0, 1),
        MockCommand::CX(1, 2),
        MockCommand::Rz(2, std::f64::consts::FRAC_PI_4),
    ]
);

connectivity_tests!(
    line_repeated_s_cx_rz,
    3,
    Connectivity::line(3),
    vec![
        MockCommand::S(0),
        MockCommand::CX(0, 1),
        MockCommand::Rz(1, 1.5),
        MockCommand::S(0),
        MockCommand::CX(0, 1),
        MockCommand::Rz(1, 1.5),
        MockCommand::S(0),
        MockCommand::CX(0, 1),
        MockCommand::Rz(1, 1.5),
    ]
);

// The input names a coupling the hardware does not have, so synthesis has to route across the line
// rather than copy the structure it was given.
connectivity_tests!(
    line_far_cx,
    4,
    Connectivity::line(4),
    vec![
        MockCommand::CX(0, 3),
        MockCommand::Rz(3, 0.4),
        MockCommand::CX(0, 3),
    ]
);

connectivity_tests!(
    line_four_qubit_chain,
    4,
    Connectivity::line(4),
    vec![
        MockCommand::H(0),
        MockCommand::CX(0, 1),
        MockCommand::CX(1, 2),
        MockCommand::CX(2, 3),
        MockCommand::Rz(3, 0.73),
        MockCommand::Rx(1, 0.31),
        MockCommand::CX(0, 2),
    ]
);

// ---------------------------------------------------------------------------
// Other topologies, to prove nothing is accidentally line-specific
// ---------------------------------------------------------------------------

connectivity_tests!(
    grid_2x2,
    4,
    Connectivity::grid(2, 2),
    vec![
        MockCommand::H(0),
        MockCommand::CX(0, 3),
        MockCommand::Rz(3, 0.42),
        MockCommand::CX(1, 2),
        MockCommand::Rx(2, 0.17),
    ]
);

connectivity_tests!(
    ring_4,
    4,
    Connectivity::from_edges_with_size(4, &[(0, 1), (1, 2), (2, 3), (3, 0)]).unwrap(),
    vec![
        MockCommand::CX(0, 2),
        MockCommand::Rz(2, std::f64::consts::FRAC_PI_4),
        MockCommand::CX(1, 3),
        MockCommand::Rz(3, 0.9),
    ]
);

connectivity_tests!(
    star_4,
    4,
    Connectivity::from_edges_with_size(4, &[(0, 1), (0, 2), (0, 3)]).unwrap(),
    vec![
        MockCommand::CX(1, 2),
        MockCommand::Rz(2, 0.55),
        MockCommand::CX(2, 3),
        MockCommand::Rz(3, 0.25),
    ]
);

// ---------------------------------------------------------------------------
// The documented non-property, and the validation failures
// ---------------------------------------------------------------------------

/// A `Naive` half ignores the graph outright, even when the other half honours it.
///
/// Pinned here because it is the reason `synthesize_to_qiskit` warns when a connectivity is paired
/// with a `Naive` strategy: the result is still a correct circuit, just not one the hardware can
/// run. `Naive`/`PermRowCol` is the case that matters most, since it is `PyPauliExponential`'s
/// default pair. If a future change makes either `Naive` synthesizer conform, this test fails and
/// the warning can go.
#[test]
fn naive_half_leaves_gates_off_the_graph() {
    // A CX between the endpoints of a line: no conforming synthesis can emit it directly.
    let commands = vec![
        MockCommand::CX(0, 2),
        MockCommand::Rz(2, std::f64::consts::FRAC_PI_4),
        MockCommand::CX(0, 2),
    ];
    let connectivity = Connectivity::line(3);
    let gates: Vec<ImportedGate> = commands.iter().map(to_imported).collect();

    for (strategy_name, pauli, tableau) in [
        (
            "Naive/Naive",
            PauliPolynomialSynthStrategy::Naive,
            CliffordTableauSynthStrategy::Naive,
        ),
        (
            "Naive/PermRowCol",
            PauliPolynomialSynthStrategy::Naive,
            CliffordTableauSynthStrategy::PermRowCol,
        ),
    ] {
        let mut mock = MockCircuit::new();
        let mut synthesizer = PauliExponentialSynthesizer::from_strategy(pauli, tableau);
        synthesizer.set_connectivity(connectivity.clone());
        synthesizer.synthesize(import_circuit(3, &gates), &mut mock);
        assert!(
            !mock.connectivity_violations(&connectivity).is_empty(),
            "{strategy_name} unexpectedly produced a conforming circuit: {:?}",
            mock.commands()
        );
    }

    // ... and the same input *does* conform under PSGS/PermRowCol, so the difference is the
    // strategy and not the circuit.
    assert_round_trip_on(
        "naive_half_leaves_gates_off_the_graph",
        3,
        &commands,
        &connectivity,
        "PSGS/PermRowCol",
        PauliPolynomialSynthStrategy::PSGS,
        CliffordTableauSynthStrategy::PermRowCol,
        true,
    );
}

#[test]
#[should_panic(expected = "connectivity is disconnected")]
fn disconnected_graph_is_rejected() {
    let commands = vec![MockCommand::CX(0, 1), MockCommand::Rz(1, 0.3)];
    assert_round_trip_on(
        "disconnected",
        4,
        &commands,
        &Connectivity::from_edges_with_size(4, &[(0, 1), (2, 3)]).unwrap(),
        "PSGS/PermRowCol",
        PauliPolynomialSynthStrategy::PSGS,
        CliffordTableauSynthStrategy::PermRowCol,
        true,
    );
}

#[test]
#[should_panic(expected = "the circuit has 3 qubits")]
fn oversized_graph_is_rejected() {
    let commands = vec![MockCommand::CX(0, 1), MockCommand::Rz(1, 0.3)];
    assert_round_trip_on(
        "oversized",
        3,
        &commands,
        &Connectivity::line(5),
        "PSGS/PermRowCol",
        PauliPolynomialSynthStrategy::PSGS,
        CliffordTableauSynthStrategy::PermRowCol,
        true,
    );
}

#[test]
#[should_panic(expected = "the circuit has 4 qubits")]
fn undersized_graph_is_rejected() {
    let commands = vec![MockCommand::CX(0, 1), MockCommand::Rz(1, 0.3)];
    assert_round_trip_on(
        "undersized",
        4,
        &commands,
        &Connectivity::line(2),
        "PSGS/PermRowCol",
        PauliPolynomialSynthStrategy::PSGS,
        CliffordTableauSynthStrategy::PermRowCol,
        true,
    );
}
