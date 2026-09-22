//! Tests that anchor `synir`'s tableau algebra to an absolute, Qiskit-compatible meaning.
//!
//! Everything else under `synir/tests` checks the synthesizers against the same `PropagateClifford`
//! code the synthesizers use, which stays self-consistent under a global adjoint or transpose flip.
//! These tests compare against the dense reference in `common::unitary` instead, so a flip has
//! nowhere to hide. If these are green, the Python boundary in `synpy` is being fed a tableau that
//! means what `qiskit.quantum_info.Clifford` means.

mod common;

use common::mock_circuit::{MockCircuit, MockCommand};
use common::tableau_from_unitary::{pauli_word_of, tableau_from_unitary};
use common::unitary::{self, Unitary, TOL};

use bitvec::bitvec;
use bitvec::prelude::Lsb0;
use synir::data_structures::{CliffordTableau, HasAdjoint, PauliLetter, PauliString};

// ---------------------------------------------------------------------------
// 0. The reference itself. If these fail, nothing below means anything.
// ---------------------------------------------------------------------------

#[test]
fn reference_identity_tableau_is_identity() {
    for n in 1..=3 {
        assert_eq!(
            tableau_from_unitary(&Unitary::identity(n)),
            CliffordTableau::new(n),
            "identity unitary on {n} qubits must give the identity tableau"
        );
    }
}

#[test]
fn reference_reads_back_every_pauli_word() {
    // Every signed Pauli word on 2 qubits must round-trip through `pauli_word_of`.
    for text in ["II", "IX", "IY", "IZ", "XI", "YZ", "ZY", "XX", "YY", "ZZ"] {
        let pauli = PauliString::from_text(text);
        let matrix = Unitary::from_pauli_string(&pauli);
        let word = pauli_word_of(&matrix);
        assert!(!word.negative, "{text} should read back positive");

        let negated = matrix.scale(-unitary::ONE);
        let negated_word = pauli_word_of(&negated);
        assert!(negated_word.negative, "-{text} should read back negative");
        assert_eq!(word.x_mask, negated_word.x_mask);
        assert_eq!(word.z_mask, negated_word.z_mask);

        // And the masks must agree with the letters, little-endian.
        for q in 0..2 {
            let (x, z) = (word.x_mask >> q & 1 == 1, word.z_mask >> q & 1 == 1);
            assert_eq!(
                PauliLetter::new(x, z),
                pauli.pauli(q),
                "{text}: qubit {q} decoded incorrectly"
            );
        }
    }
}

#[test]
fn reference_cx_matches_hand_computed_tableau() {
    // Independently of any `synir` code: CX(0,1) maps X0 -> X0X1, X1 -> X1, Z0 -> Z0, Z1 -> Z0Z1.
    let mut u = Unitary::identity(2);
    u.apply_cx(0, 1);
    let ct = tableau_from_unitary(&u);

    // Written out in synir's storage layout: one PauliString per *qubit*, indexed by generator
    // row. Reading down a column gives the letters that qubit carries in rows
    // [X0, X1, Z0, Z1] = [X0X1, X1, Z0, Z0Z1].
    let expected = CliffordTableau::from_parts(
        vec![
            PauliString::from_text("XIZZ"), // qubit 0
            PauliString::from_text("XXIZ"), // qubit 1
        ],
        bitvec![0, 0, 0, 0],
    );
    assert_eq!(
        ct, expected,
        "\nCX(0,1) tableau:\n{ct}\nhand-computed:\n{expected}\n"
    );
}

// ---------------------------------------------------------------------------
// 1. Per-gate: does `PropagateClifford::<gate>` mean the same thing as the gate matrix?
// ---------------------------------------------------------------------------

/// Builds the one-gate circuit both ways and compares the resulting tableaus.
fn assert_gate_convention(name: &str, num_qubits: usize, command: MockCommand) {
    let circuit = MockCircuit::from_vec(vec![command.clone()]);
    let from_matrix = tableau_from_unitary(&circuit.to_unitary(num_qubits));
    let from_propagation = circuit.to_clifford_tableau(num_qubits);

    if from_matrix == from_propagation {
        return;
    }
    let adjoint_hint = if tableau_from_unitary(&circuit.to_unitary(num_qubits).adjoint())
        == from_propagation
    {
        "\nhint: PropagateClifford matches the ADJOINT of the gate matrix (U^dag P U instead of U P U^dag)"
    } else {
        ""
    };
    panic!(
        "\n{name}: PropagateClifford disagrees with the gate matrix.\n\
         from matrix:\n{from_matrix}\nfrom propagation:\n{from_propagation}{adjoint_hint}\n"
    );
}

#[test]
fn single_qubit_gate_conventions() {
    assert_gate_convention("s", 1, MockCommand::S(0));
    assert_gate_convention("s_dgr", 1, MockCommand::SDgr(0));
    assert_gate_convention("v", 1, MockCommand::V(0));
    assert_gate_convention("v_dgr", 1, MockCommand::VDgr(0));
    assert_gate_convention("x", 1, MockCommand::X(0));
    assert_gate_convention("y", 1, MockCommand::Y(0));
    assert_gate_convention("z", 1, MockCommand::Z(0));
    assert_gate_convention("h", 1, MockCommand::H(0));
}

#[test]
fn two_qubit_gate_conventions() {
    assert_gate_convention("cx(0,1)", 2, MockCommand::CX(0, 1));
    assert_gate_convention("cx(1,0)", 2, MockCommand::CX(1, 0));
    assert_gate_convention("cz(0,1)", 2, MockCommand::CZ(0, 1));
    assert_gate_convention("cz(1,0)", 2, MockCommand::CZ(1, 0));
}

#[test]
fn single_qubit_gates_act_on_the_named_qubit() {
    // A gate on qubit 1 of a 3-qubit register must leave qubits 0 and 2 alone. This catches a
    // reversed qubit index (big- vs little-endian) that a 1-qubit test cannot see.
    for command in [
        MockCommand::S(1),
        MockCommand::V(1),
        MockCommand::H(1),
        MockCommand::Y(1),
    ] {
        assert_gate_convention(&format!("{command:?} on 3 qubits"), 3, command);
    }
}

// ---------------------------------------------------------------------------
// 2. Ordering: is a gate list read first-gate-first?
// ---------------------------------------------------------------------------

/// Gate lists whose unitary is *not* self-adjoint and *not* order-symmetric.
///
/// This matters: the existing samples are mostly order-blind. `sample_s_dgr_gate` expects
/// `[S, Z]` and `sample_v_dgr_gate` expects `[V, X]`, but S and Z commute and so do V and X, so
/// those assertions hold under either reading. `sample_swap_ct` is a palindrome. Only the lists
/// below can distinguish a reversed circuit.
fn order_sensitive_circuits() -> Vec<(&'static str, usize, Vec<MockCommand>)> {
    vec![
        ("h then s", 1, vec![MockCommand::H(0), MockCommand::S(0)]),
        ("s then h", 1, vec![MockCommand::S(0), MockCommand::H(0)]),
        ("v then s", 1, vec![MockCommand::V(0), MockCommand::S(0)]),
        ("bell", 2, vec![MockCommand::H(0), MockCommand::CX(0, 1)]),
        (
            "bell reversed",
            2,
            vec![MockCommand::CX(0, 1), MockCommand::H(0)],
        ),
        (
            "half swap 0-1",
            2,
            vec![MockCommand::CX(0, 1), MockCommand::CX(1, 0)],
        ),
        (
            "half swap 1-0",
            2,
            vec![MockCommand::CX(1, 0), MockCommand::CX(0, 1)],
        ),
        (
            "mixed 3q",
            3,
            vec![
                MockCommand::H(0),
                MockCommand::S(1),
                MockCommand::CX(0, 1),
                MockCommand::V(2),
                MockCommand::CX(2, 0),
                MockCommand::SDgr(1),
                MockCommand::CZ(1, 2),
                MockCommand::Y(0),
            ],
        ),
    ]
}

#[test]
fn gate_lists_are_read_first_gate_first() {
    for (name, n, commands) in order_sensitive_circuits() {
        let circuit = MockCircuit::from_vec(commands.clone());
        let from_matrix = tableau_from_unitary(&circuit.to_unitary(n));
        let from_propagation = circuit.to_clifford_tableau(n);

        if from_matrix == from_propagation {
            continue;
        }

        let mut reversed = commands.clone();
        reversed.reverse();
        let reversed_circuit = MockCircuit::from_vec(reversed);
        let mut hint = String::new();
        if tableau_from_unitary(&reversed_circuit.to_unitary(n)) == from_propagation {
            hint.push_str(
                "\nhint: the REVERSED gate list matches, so the list is being read last-gate-first",
            );
        }
        if tableau_from_unitary(&circuit.to_unitary(n).adjoint()) == from_propagation {
            hint.push_str("\nhint: the ADJOINT of the circuit matches");
        }
        panic!(
            "\n{name}: to_clifford_tableau disagrees with the dense circuit unitary.\n\
             commands: {commands:?}\nfrom matrix:\n{from_matrix}\n\
             from propagation:\n{from_propagation}{hint}\n"
        );
    }
}

// ---------------------------------------------------------------------------
// 3. compose / adjoint
// ---------------------------------------------------------------------------

#[test]
fn adjoint_matches_matrix_adjoint() {
    for (name, n, commands) in order_sensitive_circuits() {
        let circuit = MockCircuit::from_vec(commands);
        let ct = circuit.to_clifford_tableau(n);
        let expected = tableau_from_unitary(&circuit.to_unitary(n).adjoint());
        assert_eq!(
            ct.adjoint(),
            expected,
            "\n{name}: CliffordTableau::adjoint disagrees with the matrix adjoint.\n\
             got:\n{}\nexpected:\n{expected}\n",
            ct.adjoint()
        );
    }
}

#[test]
fn compose_matches_matrix_product() {
    // `compose` is documented as `-[lhs]-[rhs]-`, i.e. lhs runs first, so the matrix product is
    // U_rhs * U_lhs.
    let cases = order_sensitive_circuits();
    for (lhs_name, n, lhs_commands) in cases.iter() {
        for (rhs_name, rhs_n, rhs_commands) in cases.iter() {
            if n != rhs_n {
                continue;
            }
            let lhs = MockCircuit::from_vec(lhs_commands.clone());
            let rhs = MockCircuit::from_vec(rhs_commands.clone());

            let composed = lhs
                .to_clifford_tableau(*n)
                .compose(&rhs.to_clifford_tableau(*n));
            let expected = tableau_from_unitary(&rhs.to_unitary(*n).mul(&lhs.to_unitary(*n)));

            assert_eq!(
                composed, expected,
                "\ncompose({lhs_name}, {rhs_name}): disagrees with the matrix product \
                 U_rhs * U_lhs.\ngot:\n{composed}\nexpected:\n{expected}\n"
            );
        }
    }
}

#[test]
fn concatenating_gate_lists_equals_composing_tableaus() {
    let cases = order_sensitive_circuits();
    for (lhs_name, n, lhs_commands) in cases.iter() {
        for (rhs_name, rhs_n, rhs_commands) in cases.iter() {
            if n != rhs_n {
                continue;
            }
            let mut joined = lhs_commands.clone();
            joined.extend(rhs_commands.iter().cloned());

            let direct = MockCircuit::from_vec(joined).to_clifford_tableau(*n);
            let composed = MockCircuit::from_vec(lhs_commands.clone())
                .to_clifford_tableau(*n)
                .compose(&MockCircuit::from_vec(rhs_commands.clone()).to_clifford_tableau(*n));

            assert_eq!(
                direct, composed,
                "\n{lhs_name} ++ {rhs_name}: running both gate lists differs from composing their \
                 tableaus.\ndirect:\n{direct}\ncomposed:\n{composed}\n"
            );
        }
    }
}

// ---------------------------------------------------------------------------
// 4. Permutations
// ---------------------------------------------------------------------------

#[test]
fn get_permutation_round_trips_through_the_dense_permutation() {
    // Build a tableau that is a pure qubit permutation using SWAPs, ask for its permutation, and
    // check the reported pattern reproduces the same unitary under one documented reading.
    let patterns: Vec<Vec<usize>> = vec![
        vec![1, 0],
        vec![1, 2, 0],
        vec![2, 0, 1],
        vec![0, 2, 1],
        vec![3, 2, 1, 0],
    ];

    for pattern in patterns {
        let n = pattern.len();
        let gates = swap_network(&pattern);
        let circuit = MockCircuit::from_vec(gates);
        let ct = circuit.to_clifford_tableau(n);

        let reported = ct.get_permutation().unwrap_or_else(|| {
            panic!("tableau for {pattern:?} is not recognised as a permutation")
        });

        let actual = circuit.to_unitary(n);
        let matches: Vec<&str> = [
            ("forward", Unitary::permutation(&reported, n)),
            (
                "inverse",
                Unitary::permutation(&Unitary::invert_pattern(&reported), n),
            ),
        ]
        .into_iter()
        .filter(|(_, candidate)| candidate.equiv_up_to_global_phase(&actual, TOL))
        .map(|(label, _)| label)
        .collect();

        assert!(
            !matches.is_empty(),
            "\nget_permutation returned {reported:?} for the swap network implementing \
             {pattern:?}, but neither PermutationGate({reported:?}) nor its inverse reproduces \
             the circuit unitary"
        );
    }
}

#[test]
fn final_permutation_is_stored_not_expanded() {
    // Guards the alignment between MockCircuit and QiskitSynIR: both must *store* the permutation.
    // If MockCircuit ever falls back to the default CX-swap expansion, the Rust tests stop
    // exercising the path that synpy takes.
    use synir::ir::CliffordGates;

    let mut circuit = MockCircuit::new();
    circuit.add_final_permutation(vec![1, 2, 0]);
    assert!(
        circuit.commands().is_empty(),
        "add_final_permutation must not emit gates, but emitted {:?}",
        circuit.commands()
    );
    assert_eq!(circuit.final_permutation(), Some(&vec![1, 2, 0]));
}

/// A gate list implementing the qubit permutation `pattern` in Qiskit's `PermutationGate` reading.
fn swap_network(pattern: &[usize]) -> Vec<MockCommand> {
    let mut current: Vec<usize> = (0..pattern.len()).collect();
    let mut gates = Vec::new();
    // Target: after the network, position k holds what qubit pattern[k] held.
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
    gates
}
