//! End-to-end `PauliExponential` round trips.
//!
//! This is the Rust analogue of `round_loop`/`check_equiv` in
//! `synpy/integration_tests/test_qiskit.py`: import a circuit, synthesize it back out, and compare
//! the two dense unitaries up to global phase. The import goes through `synir::import`, which is
//! the same code `synpy`'s `PyPauliExponential` drives, so a failure here is a failure there.

mod common;

use std::collections::VecDeque;

use common::mock_circuit::{MockCircuit, MockCommand};
use common::unitary::TOL;

use synir::data_structures::{CliffordTableau, PauliExponential, PauliPolynomial};
use synir::import::{import_circuit, ImportedGate};
use synir::ir::{
    clifford_tableau::CliffordTableauSynthStrategy, pauli_exponential::PauliExponentialSynthesizer,
    pauli_polynomial::PauliPolynomialSynthStrategy, Synthesizer,
};

fn run_synthesizer(pe: PauliExponential) -> MockCircuit {
    let mut mock: MockCircuit = MockCircuit::new();
    let mut synthesizer = PauliExponentialSynthesizer::from_strategy(
        PauliPolynomialSynthStrategy::Naive,
        CliffordTableauSynthStrategy::Naive,
    );
    synthesizer.synthesize(pe, &mut mock);
    return mock;
}

#[test]
fn test_empty_pe() {
    let pe = PauliExponential::new(VecDeque::new(), CliffordTableau::new(5));
    let mock = run_synthesizer(pe);
    let ref_commands = [];
    assert_eq!(mock.commands(), &ref_commands);
}

#[test]
fn test_empty_pp_in_pe() {
    let pe = PauliExponential::new(
        VecDeque::from([PauliPolynomial::new(5)]),
        CliffordTableau::new(5),
    );
    let mock = run_synthesizer(pe);
    let ref_commands = [];
    assert_eq!(mock.commands(), &ref_commands);
}

// ---------------------------------------------------------------------------
// Round-trip harness
// ---------------------------------------------------------------------------

/// The same strategy matrix pytest parametrises over (`all_strats`).
fn all_strategies() -> Vec<(
    &'static str,
    PauliPolynomialSynthStrategy,
    CliffordTableauSynthStrategy,
)> {
    vec![
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
        (
            "PSGS/Naive",
            PauliPolynomialSynthStrategy::PSGS,
            CliffordTableauSynthStrategy::Naive,
        ),
        (
            "PSGS/PermRowCol",
            PauliPolynomialSynthStrategy::PSGS,
            CliffordTableauSynthStrategy::PermRowCol,
        ),
    ]
}

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

/// Imports `reference`, synthesizes it back out, and asserts the two unitaries agree.
///
/// Mirrors `synpy.qiskit.plugin.circuit_to_circuit`: the reported permutation is appended as a
/// Qiskit `PermutationGate` built from the inverse pattern.
fn assert_round_trip(
    name: &str,
    num_qubits: usize,
    reference: &[MockCommand],
    strategy_name: &str,
    pauli_strategy: PauliPolynomialSynthStrategy,
    tableau_strategy: CliffordTableauSynthStrategy,
) {
    let gates: Vec<ImportedGate> = reference.iter().map(to_imported).collect();
    let pe = import_circuit(num_qubits, &gates);

    let mut mock = MockCircuit::new();
    let mut synthesizer =
        PauliExponentialSynthesizer::from_strategy(pauli_strategy, tableau_strategy.clone());
    synthesizer.synthesize(pe, &mut mock);

    let permutation = mock.final_permutation().cloned();
    let expected = MockCircuit::from_vec(reference.to_vec()).to_unitary(num_qubits);
    let got = mock.to_unitary_with_permutation(num_qubits, permutation.clone());

    if got.equiv_up_to_global_phase(&expected, TOL) {
        // `circuit_to_circuit` drops the permutation entirely for the Naive tableau strategy, so a
        // non-trivial permutation there is a silent wrong answer on the Python side even though
        // the Rust round trip above is fine.
        if matches!(tableau_strategy, CliffordTableauSynthStrategy::Naive) {
            if let Some(pattern) = permutation.as_ref() {
                let identity: Vec<usize> = (0..num_qubits).collect();
                assert_eq!(
                    pattern, &identity,
                    "\n{name} [{strategy_name}]: the Naive tableau strategy reported a \
                     non-identity permutation {pattern:?}, but circuit_to_circuit only appends a \
                     PermutationGate when ct_strat != \"Naive\", so synpy silently drops it."
                );
            }
        }
        return;
    }

    let mut hint = String::new();
    if let Some(pattern) = permutation.as_ref() {
        for (label, candidate) in mock.permutation_variants(num_qubits, pattern) {
            if candidate.equiv_up_to_global_phase(&expected, TOL) {
                hint = format!("\nhint: this WOULD match with the permutation as: {label}");
            }
        }
        if mock
            .to_unitary(num_qubits)
            .equiv_up_to_global_phase(&expected, TOL)
        {
            hint.push_str("\nhint: the circuit matches with NO permutation applied at all");
        }
    }
    if got.adjoint().equiv_up_to_global_phase(&expected, TOL) {
        hint.push_str("\nhint: the ADJOINT of the synthesized circuit matches the input");
    }
    let mut reversed = mock.commands().clone();
    reversed.reverse();
    if MockCircuit::from_vec(reversed)
        .to_unitary(num_qubits)
        .equiv_up_to_global_phase(&expected, TOL)
    {
        hint.push_str("\nhint: the REVERSED synthesized circuit matches the input");
    }

    panic!(
        "\n{name} [{strategy_name}]: round trip does not preserve the unitary.\n\
         input:  {reference:?}\noutput: {:?}\npermutation: {permutation:?}{hint}\n",
        mock.commands()
    );
}

/// Declares one round-trip test per strategy combination, like pytest's `@parametrize`.
macro_rules! round_trip_tests {
    ($name:ident, $num_qubits:expr, $commands:expr) => {
        paste::item! {
            #[test]
            fn [< round_trip_ $name >]() {
                let commands: Vec<MockCommand> = $commands;
                let mut failures = Vec::new();
                for (strategy_name, pauli, tableau) in all_strategies() {
                    let result = std::panic::catch_unwind(|| {
                        assert_round_trip(
                            stringify!($name),
                            $num_qubits,
                            &commands,
                            strategy_name,
                            pauli,
                            tableau,
                        );
                    });
                    if result.is_err() {
                        failures.push(strategy_name);
                    }
                }
                assert!(
                    failures.is_empty(),
                    "{}: round trip failed for strategies {:?} (details above)",
                    stringify!($name),
                    failures
                );
            }
        }
    };
}

// ---------------------------------------------------------------------------
// Cases ported from synpy/integration_tests/test_qiskit.py
// ---------------------------------------------------------------------------

round_trip_tests!(bell, 2, vec![MockCommand::H(0), MockCommand::CX(0, 1)]);

round_trip_tests!(single_cx, 2, vec![MockCommand::CX(0, 1)]);

round_trip_tests!(
    rz_at_start,
    2,
    vec![MockCommand::Rz(0, 0.234), MockCommand::CX(0, 1)]
);

// `test_sqg`: a bare T gate, i.e. rz(pi/4).
round_trip_tests!(
    t_gate,
    1,
    vec![MockCommand::Rz(0, std::f64::consts::FRAC_PI_4)]
);

// `test_pauli_gadget`
round_trip_tests!(
    pauli_gadget,
    3,
    vec![
        MockCommand::CX(0, 1),
        MockCommand::CX(1, 2),
        MockCommand::Rz(2, std::f64::consts::FRAC_PI_4),
        MockCommand::CX(1, 2),
        MockCommand::CX(0, 1),
    ]
);

// `test_pauli_gadget_half`
round_trip_tests!(
    pauli_gadget_half,
    3,
    vec![
        MockCommand::CX(0, 1),
        MockCommand::CX(1, 2),
        MockCommand::Rz(2, std::f64::consts::FRAC_PI_4),
    ]
);

// `test_qiskit_loop`
round_trip_tests!(
    repeated_s_cx_rz,
    3,
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

// `test_2sqg`
round_trip_tests!(
    two_single_qubit_gadgets,
    2,
    vec![
        MockCommand::Rz(0, 0.59156),
        MockCommand::Rx(0, 0.8513),
        MockCommand::CX(0, 1),
        MockCommand::Rx(0, 0.8513 * 2.0),
        MockCommand::CX(0, 1),
    ]
);

// `test_toffoli`: the standard ccx decomposition over {h, cx, rz}.
fn toffoli(a: usize, b: usize, c: usize) -> Vec<MockCommand> {
    let t = std::f64::consts::FRAC_PI_4;
    vec![
        MockCommand::H(c),
        MockCommand::CX(b, c),
        MockCommand::Rz(c, -t),
        MockCommand::CX(a, c),
        MockCommand::Rz(c, t),
        MockCommand::CX(b, c),
        MockCommand::Rz(c, -t),
        MockCommand::CX(a, c),
        MockCommand::Rz(b, t),
        MockCommand::Rz(c, t),
        MockCommand::H(c),
        MockCommand::CX(a, b),
        MockCommand::Rz(a, t),
        MockCommand::Rz(b, -t),
        MockCommand::CX(a, b),
    ]
}

round_trip_tests!(toffoli, 3, toffoli(0, 1, 2));

// `test_merge_repeats`: two Toffolis cancel to the identity.
round_trip_tests!(toffoli_twice, 3, {
    let mut commands = toffoli(0, 1, 2);
    commands.extend(toffoli(0, 1, 2));
    commands
});

// `test_each_pauli_push`: one leading Clifford, then a gadget pair.
macro_rules! pauli_push_test {
    ($name:ident, $lead:expr) => {
        round_trip_tests!(
            $name,
            2,
            vec![
                $lead,
                MockCommand::CX(0, 1),
                MockCommand::Rx(0, 0.234234),
                MockCommand::Rz(1, 0.53423),
            ]
        );
    };
}

pauli_push_test!(push_s, MockCommand::S(0));
pauli_push_test!(push_sdg, MockCommand::SDgr(0));
pauli_push_test!(push_v, MockCommand::V(0));
pauli_push_test!(push_vdg, MockCommand::VDgr(0));
pauli_push_test!(push_y, MockCommand::Y(0));
pauli_push_test!(push_x, MockCommand::X(0));
pauli_push_test!(push_z, MockCommand::Z(0));

// A fixed stand-in for `quantum_volume(3, 1)`: a Clifford+T circuit dense enough to exercise the
// tableau strategies without depending on Qiskit's random generator.
round_trip_tests!(
    dense_3q,
    3,
    vec![
        MockCommand::H(0),
        MockCommand::Rz(0, 0.31),
        MockCommand::CX(0, 1),
        MockCommand::H(1),
        MockCommand::Rz(1, 1.07),
        MockCommand::CX(1, 2),
        MockCommand::Rz(2, -0.62),
        MockCommand::H(2),
        MockCommand::CX(2, 0),
        MockCommand::Rz(0, 2.4),
        MockCommand::CX(0, 1),
    ]
);

round_trip_tests!(
    dense_4q,
    4,
    vec![
        MockCommand::H(0),
        MockCommand::CX(0, 1),
        MockCommand::Rz(1, 0.77),
        MockCommand::H(2),
        MockCommand::CX(2, 3),
        MockCommand::Rz(3, -1.21),
        MockCommand::CX(1, 2),
        MockCommand::Rz(2, 0.44),
        MockCommand::H(1),
        MockCommand::CX(3, 0),
        MockCommand::Rz(0, 1.9),
    ]
);

// ---------------------------------------------------------------------------
// Seeded randomised round trips
// ---------------------------------------------------------------------------

/// A pseudo-random Clifford+T circuit. Seeded so a failure is reproducible -- note that
/// `ct_synthesis/custom_callback.rs` still uses an unseeded `rand::rng()`.
fn random_circuit(rng: &mut impl rand::Rng, num_qubits: usize, depth: usize) -> Vec<MockCommand> {
    let mut commands = Vec::with_capacity(depth);
    for _ in 0..depth {
        let target = rng.random_range(0..num_qubits);
        match rng.random_range(0..5) {
            0 => commands.push(MockCommand::H(target)),
            1 => commands.push(MockCommand::S(target)),
            2 => commands.push(MockCommand::Rz(target, rng.random_range(0.1..3.0))),
            3 => commands.push(MockCommand::Rx(target, rng.random_range(0.1..3.0))),
            _ if num_qubits > 1 => {
                let mut control = rng.random_range(0..num_qubits);
                while control == target {
                    control = rng.random_range(0..num_qubits);
                }
                commands.push(MockCommand::CX(control, target));
            }
            _ => commands.push(MockCommand::H(target)),
        }
    }
    commands
}

#[test]
fn round_trip_random_circuits() {
    use rand::SeedableRng;

    let mut failures = Vec::new();
    for num_qubits in 1..=4usize {
        for seed in 0..3u64 {
            let mut rng = rand::rngs::StdRng::seed_from_u64(seed * 100 + num_qubits as u64);
            let commands = random_circuit(&mut rng, num_qubits, 4 * num_qubits);
            for (strategy_name, pauli, tableau) in all_strategies() {
                let name = format!("random n={num_qubits} seed={seed}");
                let result = std::panic::catch_unwind(|| {
                    assert_round_trip(&name, num_qubits, &commands, strategy_name, pauli, tableau);
                });
                if result.is_err() {
                    failures.push(format!("{name} [{strategy_name}]"));
                }
            }
        }
    }
    assert!(
        failures.is_empty(),
        "random round trips failed for {} of 48 cases: {:?}",
        failures.len(),
        failures
    );
}

/// A gate list implementing Qiskit's `PermutationGate(pattern)`, i.e. output bit `k` is input bit
/// `pattern[k]`. Verified against the dense permutation matrix before use, so the fixture cannot
/// silently encode the wrong direction.
fn swap_network(pattern: &[usize]) -> Vec<MockCommand> {
    let mut current: Vec<usize> = (0..pattern.len()).collect();
    let mut gates = Vec::new();
    for position in 0..pattern.len() {
        let wanted = pattern[position];
        let found = current.iter().position(|&q| q == wanted).unwrap();
        if found != position {
            gates.push(MockCommand::CX(position, found));
            gates.push(MockCommand::CX(found, position));
            gates.push(MockCommand::CX(position, found));
            current.swap(position, found);
        }
    }

    let n = pattern.len();
    assert!(
        MockCircuit::from_vec(gates.clone())
            .to_unitary(n)
            .equiv_up_to_global_phase(&common::unitary::Unitary::permutation(pattern, n), TOL),
        "swap_network({pattern:?}) does not implement PermutationGate({pattern:?})"
    );
    gates
}

// `test_final_permutation`: a bare PermutationGate, which qiskit transpiles into a swap network.
round_trip_tests!(final_permutation, 5, swap_network(&[4, 2, 3, 0, 1]));
