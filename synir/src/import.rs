//! Building a [`PauliExponential`] from a gate list.
//!
//! This is the Rust home of the logic that used to live only in `synpy`: the per-gate handling in
//! `PyPauliExponential::add_*`, the Clifford-angle special-casing in `add_rz`, and the
//! preprocessing passes that run before synthesis. Keeping it here means the Rust test suite can
//! exercise the exact same code path that `synpy.qiskit.plugin.qiskit_to_synir` drives, rather
//! than a lookalike that can drift away from it.
//!
//! # Traversal order and the two handednesses
//!
//! A circuit is imported **back to front**: the last gate executed is added first. Each `add_rz`
//! then does a `push_front`, so the resulting deque of Pauli polynomials ends up in circuit order,
//! and each Clifford conjugates the gadgets that come after it. [`import_circuit`] does the
//! reversal for you; the incremental [`CircuitImporter`] API leaves it to the caller, because
//! `qiskit_to_synir` already iterates `reversed(circuit.data)` on the Python side.
//!
//! The target normal form is `U = T * U_pp`: the gadgets run first, then the Clifford `T`. Walking
//! `U = G_m ... G_1` backwards and meeting a Clifford `C` therefore requires
//!
//! * `Q -> C^dag Q C` on every gadget collected so far, and
//! * `T <- T * C`, a **prepend**.
//!
//! `PropagateClifford for PauliExponential` cannot be used for this. It *appends* to the tableau
//! (`T <- C * T`), which is right for its other caller -- `PSGSPauliPolynomialSynthesizer` pushes
//! an already-*emitted* gate through the remaining work -- but for a backward walk it yields
//! `C_1 * C_2 * ... * C_m`, the reverse-order product. So the importer drives the two fields
//! itself, below.
//!
//! Prepending needs no new algorithm: `(T * C)^dag = C^dag * T^dag`, so the importer keeps the
//! tableau as `T^dag` and *appends the adjoint gate*, then inverts once in [`CircuitImporter::finish`].
//! While that inverse window is open the tableau is NOT the value the `clifford_tableau` field
//! promises; `finish` closes it.

use std::collections::VecDeque;

use crate::data_structures::pauli_exponential::{
    merge_commuting_pp, merge_repeats, push_clifford_angles_inverse,
};
use crate::data_structures::{
    Angle, CliffordTableau, HasAdjoint, PauliExponential, PauliPolynomial, PropagateClifford,
};
use crate::IndexType;

/// A gate as it arrives from a foreign circuit representation.
///
/// This mirrors the basis `qiskit_to_synir` transpiles to (`cx`, `h`, `rz`), plus the extra cases
/// that function accepts directly.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ImportedGate {
    H(IndexType),
    S(IndexType),
    SDgr(IndexType),
    V(IndexType),
    VDgr(IndexType),
    X(IndexType),
    Y(IndexType),
    Z(IndexType),
    CX(IndexType, IndexType),
    Rx(IndexType, f64),
    Rz(IndexType, f64),
}

/// Incrementally builds a [`PauliExponential`] from gates supplied in reverse circuit order.
pub struct CircuitImporter {
    pe: PauliExponential,
}

impl CircuitImporter {
    pub fn new(num_qubits: usize) -> Self {
        Self {
            pe: PauliExponential::new(VecDeque::new(), CliffordTableau::new(num_qubits)),
        }
    }

    pub fn num_qubits(&self) -> usize {
        self.pe.size()
    }

    // Each Clifford below conjugates the gadgets collected so far (`Q -> C^dag Q C`, which is what
    // `PropagateClifford for PauliPolynomial` does) and appends the *adjoint* gate to the tableau,
    // which is held in inverse form until `finish`. `h`, `x`, `y`, `z` and `cx` are self-adjoint,
    // so only `s` and `v` differ between the two sides. The pairs are written out one by one
    // rather than leaning on the trait's default expansions, whose call order is append-shaped.

    pub fn add_h(&mut self, target: IndexType) {
        self.pe.mut_chains().h(target);
        self.pe.mut_clifford_tableau().h(target);
    }

    pub fn add_v(&mut self, target: IndexType) {
        self.pe.mut_chains().v(target);
        self.pe.mut_clifford_tableau().v_dgr(target);
    }

    pub fn add_v_dgr(&mut self, target: IndexType) {
        self.pe.mut_chains().v_dgr(target);
        self.pe.mut_clifford_tableau().v(target);
    }

    pub fn add_s(&mut self, target: IndexType) {
        self.pe.mut_chains().s(target);
        self.pe.mut_clifford_tableau().s_dgr(target);
    }

    pub fn add_s_dgr(&mut self, target: IndexType) {
        self.pe.mut_chains().s_dgr(target);
        self.pe.mut_clifford_tableau().s(target);
    }

    pub fn add_x(&mut self, target: IndexType) {
        self.pe.mut_chains().x(target);
        self.pe.mut_clifford_tableau().x(target);
    }

    pub fn add_y(&mut self, target: IndexType) {
        self.pe.mut_chains().y(target);
        self.pe.mut_clifford_tableau().y(target);
    }

    pub fn add_z(&mut self, target: IndexType) {
        self.pe.mut_chains().z(target);
        self.pe.mut_clifford_tableau().z(target);
    }

    pub fn add_cx(&mut self, control: IndexType, target: IndexType) {
        self.pe.mut_chains().cx(control, target);
        self.pe.mut_clifford_tableau().cx(control, target);
    }

    /// An `rx(angle)` as `h . rz(angle) . h`, matching `qiskit_to_synir`.
    ///
    /// The surrounding `h` calls are emitted in the same order the Python loop uses.
    pub fn add_rx(&mut self, target: IndexType, angle: f64) {
        self.add_h(target);
        self.add_rz(target, angle);
        self.add_h(target);
    }

    /// Adds a `Z` rotation, folding Clifford angles into the tableau instead of keeping a gadget.
    pub fn add_rz(&mut self, target: IndexType, angle: f64) {
        let mut angle = Angle::from_angle(angle);
        if let Ok(pi4_rotations) = angle.to_pi4_rotation() {
            match pi4_rotations {
                0 => return,
                2 => return self.add_s(target),
                4 => return self.add_z(target),
                6 => return self.add_s_dgr(target),
                n => angle = Angle::from_pi4_rotation(n), // Non-Clifford
            }
        }
        let size = self.pe.size();
        let gadget =
            PauliPolynomial::from_hamiltonian(vec![(&single_pauli_term(size, target, 'Z'), angle)]);
        self.pe.mut_chains().push_front(gadget);
    }

    pub fn add_gate(&mut self, gate: ImportedGate) {
        match gate {
            ImportedGate::H(t) => self.add_h(t),
            ImportedGate::S(t) => self.add_s(t),
            ImportedGate::SDgr(t) => self.add_s_dgr(t),
            ImportedGate::V(t) => self.add_v(t),
            ImportedGate::VDgr(t) => self.add_v_dgr(t),
            ImportedGate::X(t) => self.add_x(t),
            ImportedGate::Y(t) => self.add_y(t),
            ImportedGate::Z(t) => self.add_z(t),
            ImportedGate::CX(c, t) => self.add_cx(c, t),
            ImportedGate::Rx(t, angle) => self.add_rx(t, angle),
            ImportedGate::Rz(t, angle) => self.add_rz(t, angle),
        }
    }

    /// The Pauli exponential as built so far, without preprocessing.
    ///
    /// Its tableau is still in inverse form -- only [`Self::finish`] closes that window.
    pub fn into_raw(self) -> PauliExponential {
        self.pe
    }

    /// Borrows the Pauli exponential as built so far, for inspection and debug printing.
    ///
    /// Its tableau is still in inverse form -- only [`Self::finish`] closes that window.
    pub fn raw(&self) -> &PauliExponential {
        &self.pe
    }

    /// Runs the preprocessing passes that `synpy::wrapper::synthesize` runs before synthesis, then
    /// closes the inverse window so the returned `clifford_tableau` is the honest `T`.
    pub fn finish(mut self) -> PauliExponential {
        merge_commuting_pp(&mut self.pe); // polynomials only
        merge_repeats(&mut self.pe); // polynomials only
        push_clifford_angles_inverse(&mut self.pe); // needs the tableau still inverted
        let inverse = std::mem::take(self.pe.mut_clifford_tableau());
        *self.pe.mut_clifford_tableau() = inverse.adjoint();
        self.pe
    }
}

/// Imports a circuit given in **execution order** and returns the preprocessed
/// [`PauliExponential`].
///
/// The reversal that `qiskit_to_synir` performs with `reversed(circuit.data)` happens here.
pub fn import_circuit(num_qubits: usize, gates: &[ImportedGate]) -> PauliExponential {
    let mut importer = CircuitImporter::new(num_qubits);
    for gate in gates.iter().rev() {
        importer.add_gate(*gate);
    }
    importer.finish()
}

/// `"II..P..II"` with `pauli` on `target`. Index 0 of the string is qubit 0.
pub fn single_pauli_term(size: usize, target: IndexType, pauli: char) -> String {
    (0..size)
        .map(|i| if i == target { pauli } else { 'I' })
        .collect()
}
