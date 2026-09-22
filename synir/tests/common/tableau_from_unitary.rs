//! Builds a [`CliffordTableau`] from a dense unitary, **without using `PropagateClifford`**.
//!
//! This is the anchor that gives the Rust test suite an absolute, Qiskit-compatible opinion about
//! what a tableau means. Every other correctness check in `synir/tests` compares the synthesizer
//! against the very propagation code the synthesizer uses, so a globally consistent adjoint or
//! transpose flip stays invisible. Going through the dense matrix removes that freedom.
//!
//! # Conventions pinned here
//!
//! * A tableau stores the **Heisenberg images** `U P U^dagger` of the generators, matching
//!   `qiskit.quantum_info.Clifford`.
//! * Rows `0..n` are the destabilizers (images of `X_r`), rows `n..2n` the stabilizers
//!   (images of `Z_{r-n}`).
//! * `synir` stores this transposed relative to Qiskit: `pauli_columns[q]` is indexed by
//!   *generator row*, and `q` selects the *qubit*. So `pauli_columns[q].pauli(r)` is the Pauli
//!   letter acting on qubit `q` in the image of generator `r`.

#![allow(dead_code)]

use std::collections::VecDeque;

use bitvec::prelude::BitVec;
use num_complex::Complex64 as C;

use synir::data_structures::{CliffordTableau, HasAdjoint, PauliPolynomial, PauliString};

use super::mock_circuit::MockCircuit;
use super::unitary::{Unitary, ONE, TOL};

/// A Pauli word read off a dense matrix: `sign * i^popcount(x & z) * X^x * Z^z`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SignedPauliWord {
    /// Bitmask of qubits carrying an `X` component.
    pub x_mask: usize,
    /// Bitmask of qubits carrying a `Z` component.
    pub z_mask: usize,
    /// `true` when the word carries a `-1`.
    pub negative: bool,
}

fn i_power(exponent: usize) -> C {
    match exponent % 4 {
        0 => C::new(1.0, 0.0),
        1 => C::new(0.0, 1.0),
        2 => C::new(-1.0, 0.0),
        _ => C::new(0.0, -1.0),
    }
}

/// Reads a signed Pauli word off a dense matrix that is known to be one.
///
/// A Pauli word `W = i^popcount(a & b) * X^a * Z^b` satisfies `W|j> = i^k (-1)^(b.j) |j ^ a>`,
/// so it has exactly one nonzero per column and:
///
/// * `a` is the row index of the nonzero in column `0`;
/// * `b_k` is `1` exactly when `M[2^k ^ a, 2^k] / M[a, 0] == -1`;
/// * the sign is `M[a, 0] / i^popcount(a & b)`.
pub fn pauli_word_of(matrix: &Unitary) -> SignedPauliWord {
    try_pauli_word_of(matrix).expect("matrix is not a signed Pauli word")
}

/// Fallible [`pauli_word_of`], for callers that need to ask "is this Clifford at all?".
pub fn try_pauli_word_of(matrix: &Unitary) -> Option<SignedPauliWord> {
    let num_qubits = matrix.num_qubits();
    let dim = matrix.dim();

    let x_mask = (0..dim).find(|&row| matrix.get(row, 0).norm() > TOL)?;
    let anchor = matrix.get(x_mask, 0);

    let mut z_mask = 0usize;
    for k in 0..num_qubits {
        let column = 1usize << k;
        let ratio = matrix.get(column ^ x_mask, column) / anchor;
        if (ratio + ONE).norm() <= TOL {
            z_mask |= column;
        } else if (ratio - ONE).norm() > TOL {
            return None;
        }
    }

    let sign = anchor / i_power((x_mask & z_mask).count_ones() as usize);
    let negative = if (sign + ONE).norm() <= TOL {
        true
    } else if (sign - ONE).norm() <= TOL {
        false
    } else {
        return None;
    };

    // The reads above only look at n + 1 of the dim columns, so confirm the whole matrix really is
    // the word we decoded before reporting success.
    let mut reconstructed = Unitary::identity(num_qubits);
    for q in 0..num_qubits {
        match (x_mask >> q & 1 == 1, z_mask >> q & 1 == 1) {
            (false, false) => {}
            (true, false) => reconstructed.apply_1q(q, super::unitary::x_matrix()),
            (true, true) => reconstructed.apply_1q(q, super::unitary::y_matrix()),
            (false, true) => reconstructed.apply_1q(q, super::unitary::z_matrix()),
        }
    }
    if negative {
        reconstructed = reconstructed.scale(-ONE);
    }
    if !reconstructed.approx_eq(matrix, TOL) {
        return None;
    }

    Some(SignedPauliWord {
        x_mask,
        z_mask,
        negative,
    })
}

/// Computes the tableau of `unitary`, which must be Clifford.
pub fn tableau_from_unitary(unitary: &Unitary) -> CliffordTableau {
    try_tableau_from_unitary(unitary)
        .expect("unitary is not Clifford: it does not map every Pauli generator to a Pauli word")
}

/// Fallible [`tableau_from_unitary`]; returns `None` when `unitary` is not Clifford.
pub fn try_tableau_from_unitary(unitary: &Unitary) -> Option<CliffordTableau> {
    let num_qubits = unitary.num_qubits();
    let adjoint = unitary.adjoint();

    let mut x_bits = vec![BitVec::repeat(false, 2 * num_qubits); num_qubits];
    let mut z_bits = vec![BitVec::repeat(false, 2 * num_qubits); num_qubits];
    let mut signs = BitVec::repeat(false, 2 * num_qubits);

    for row in 0..2 * num_qubits {
        // Rows 0..n are the images of X_row, rows n..2n the images of Z_(row - n).
        let mut generator = Unitary::identity(num_qubits);
        if row < num_qubits {
            generator.apply_1q(row, super::unitary::x_matrix());
        } else {
            generator.apply_1q(row - num_qubits, super::unitary::z_matrix());
        }

        let image = unitary.mul(&generator).mul(&adjoint);
        let word = try_pauli_word_of(&image)?;

        for qubit in 0..num_qubits {
            x_bits[qubit].set(row, word.x_mask >> qubit & 1 == 1);
            z_bits[qubit].set(row, word.z_mask >> qubit & 1 == 1);
        }
        signs.set(row, word.negative);
    }

    let columns = x_bits
        .into_iter()
        .zip(z_bits)
        .map(|(x, z)| PauliString::new(x, z))
        .collect();

    Some(CliffordTableau::from_parts(columns, signs))
}

/// Asserts that the circuit `mock` emitted by a synthesizer really implements the unitary that
/// `expected` denotes, in Qiskit's convention and up to global phase.
///
/// `permutation` is the leftover qubit permutation a synthesizer may report through
/// `CliffordTableau::get_permutation`. It is applied the same way
/// `synpy.qiskit.plugin.circuit_to_circuit` applies it: as a Qiskit `PermutationGate` built from
/// the *inverse* pattern, appended after the synthesized circuit.
///
/// Unlike a reference unitary this needs nothing but the tableau itself, so it also covers the
/// samples that are built with `CliffordTableau::from_parts` and have no originating gate list.
pub fn assert_realises_tableau(
    mock: &MockCircuit,
    expected: &CliffordTableau,
    permutation: Option<Vec<usize>>,
    context: &str,
) {
    let num_qubits = expected.size();
    let realised =
        tableau_from_unitary(&mock.to_unitary_with_permutation(num_qubits, permutation.clone()));
    if realised == *expected {
        return;
    }

    // Report which other reading of the permutation (if any) would have matched: that turns a bare
    // "not equal" into an actionable statement about the convention.
    let mut hint = String::new();
    if let Some(pattern) = permutation.as_ref() {
        for (label, candidate) in mock.permutation_variants(num_qubits, pattern) {
            if tableau_from_unitary(&candidate) == *expected {
                hint = format!("\nhint: this WOULD match with the permutation as: {label}");
            }
        }
    }
    let adjoint_matches = tableau_from_unitary(
        &mock
            .to_unitary_with_permutation(num_qubits, permutation.clone())
            .adjoint(),
    ) == *expected;
    if adjoint_matches {
        hint.push_str("\nhint: the ADJOINT of the synthesized circuit matches the input tableau");
    }

    panic!(
        "\n{context}: the synthesized circuit does not implement the input tableau.\n\
         circuit: {:?}\npermutation: {:?}\nrealised:\n{realised}\nexpected:\n{expected}{hint}\n",
        mock.commands(),
        permutation,
    );
}

/// Dense unitary of a queue of Pauli polynomials, read front-to-back and, within a polynomial,
/// index 0 first.
///
/// Each gadget is `exp(-i * theta / 2 * P)`, so a `Z` gadget with angle `theta` is `rz(theta)` on
/// that qubit -- Qiskit's convention.
pub fn pauli_polynomials_to_unitary(
    polynomials: &VecDeque<PauliPolynomial>,
    num_qubits: usize,
) -> Unitary {
    let mut u = Unitary::identity(num_qubits);
    for polynomial in polynomials.iter() {
        for i in 0..polynomial.length() {
            let gadget = Unitary::pauli_rotation(
                &polynomial.pauli_gadget(i),
                polynomial.angle(i).to_radians(),
            );
            u = gadget.mul(&u);
        }
    }
    u
}

/// Asserts that a Pauli polynomial synthesizer produced a circuit that, together with the residual
/// Clifford it reports, reproduces the input.
///
/// The contract: with `C` the emitted circuit and `T` the returned tableau, the synthesizer skipped
/// the un-computation of its basis changes and accumulated them into `T`, so the input satisfies
/// `U_pp = T^dagger * C`, equivalently `C * U_pp^dagger == T`.
///
/// This is the only assertion in the suite that looks at rotation *angles*: `equals_clifford_tableau`
/// compares `cliffords_only()` and never sees them.
pub fn assert_pauli_polynomials_round_trip(
    context: &str,
    polynomials: &VecDeque<PauliPolynomial>,
    num_qubits: usize,
    mock: &MockCircuit,
    residual: &CliffordTableau,
) {
    let target = pauli_polynomials_to_unitary(polynomials, num_qubits);
    let circuit = mock.to_unitary(num_qubits);
    let leftover = circuit.mul(&target.adjoint());

    if let Some(tableau) = try_tableau_from_unitary(&leftover) {
        if tableau == *residual {
            return;
        }
    }

    // Work out which relation does hold, so the failure names the convention rather than just
    // reporting inequality.
    let negated: VecDeque<PauliPolynomial> = polynomials
        .iter()
        .map(|polynomial| {
            let mut flipped = polynomial.clone();
            for angle in flipped.mut_angles().iter_mut() {
                *angle = -*angle;
            }
            flipped
        })
        .collect();
    let target_negated = pauli_polynomials_to_unitary(&negated, num_qubits);

    let reversed: VecDeque<PauliPolynomial> = polynomials.iter().rev().cloned().collect();
    let target_reversed = pauli_polynomials_to_unitary(&reversed, num_qubits);

    let candidates: Vec<(&str, Unitary)> = vec![
        ("C * U_pp^dag", circuit.mul(&target.adjoint())),
        ("U_pp^dag * C", target.adjoint().mul(&circuit)),
        ("C^dag * U_pp", circuit.adjoint().mul(&target)),
        ("U_pp * C^dag", target.mul(&circuit.adjoint())),
        (
            "C * U_pp^dag with NEGATED angles",
            circuit.mul(&target_negated.adjoint()),
        ),
        (
            "U_pp^dag * C with NEGATED angles",
            target_negated.adjoint().mul(&circuit),
        ),
        (
            "C * U_pp^dag with REVERSED polynomial order",
            circuit.mul(&target_reversed.adjoint()),
        ),
        (
            "U_pp^dag * C with REVERSED polynomial order",
            target_reversed.adjoint().mul(&circuit),
        ),
    ];

    let mut hint = String::new();
    for (label, candidate) in candidates {
        match try_tableau_from_unitary(&candidate) {
            None => hint.push_str(&format!("\n  {label}: not Clifford")),
            Some(tableau) if tableau == *residual => {
                hint.push_str(&format!("\n  {label}: MATCHES the residual tableau"))
            }
            Some(tableau) if tableau == residual.adjoint() => {
                hint.push_str(&format!("\n  {label}: matches the ADJOINT of the residual"))
            }
            Some(_) => hint.push_str(&format!("\n  {label}: Clifford, but a different one")),
        }
    }

    panic!(
        "\n{context}: the synthesized circuit plus its residual Clifford do not reproduce the \
         input Pauli polynomials.\ncircuit: {:?}\nresidual:\n{residual}\ncandidate relations:{hint}\n",
        mock.commands()
    );
}
