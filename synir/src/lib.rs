//! **synIR** — an intermediate representation for quantum circuits, together with the
//! architecture-aware passes that synthesize it back into gates.
//!
//! # The normal form
//!
//! synIR rewrites a circuit into a [`PauliExponential`](data_structures::PauliExponential): a
//! run of Pauli polynomials followed by a single Clifford operator.
//!
//! ```text
//! U = T * U_pp
//! ```
//!
//! This is *matrix* order, so `U_pp` runs first and the Clifford `T` runs last. `U_pp` is a
//! sequence of [`PauliPolynomial`](data_structures::PauliPolynomial)s — each a bag of commuting
//! Pauli gadgets `exp(-i * theta/2 * P)` — and `T` is a
//! [`CliffordTableau`](data_structures::CliffordTableau).
//!
//! The split is what makes the representation useful: every Clifford in the input circuit is
//! *absorbed* into `T` by conjugating the gadgets that follow it, so the only things left in
//! `U_pp` are the genuinely non-Clifford rotations. Synthesis then gets to choose fresh Clifford
//! structure for both halves, subject to a hardware connectivity graph.
//!
//! Circuit order and matrix order are opposites, and nearly every sign and ordering bug in this
//! crate traces back to confusing them. The convention used throughout is:
//!
//! * **append** a gate `G` — place it at the *end* of the circuit, i.e. `T <- G * T`.
//! * **prepend** a gate `G` — place it at the *start*, i.e. `T <- T * G`.
//! * [`CliffordTableau::compose`](data_structures::CliffordTableau::compose) is written in
//!   circuit order: `a.compose(&b)` is the circuit `-[a]-[b]-`, whose matrix is `U_b * U_a`.
//!
//! # Module map
//!
//! | Module | Contents |
//! |---|---|
//! | [`data_structures`] | The IR itself: `PauliString`, `PauliPolynomial`, `CliffordTableau`, `PauliExponential`, `Angle` |
//! | [`import`] | Building a `PauliExponential` from a foreign gate list |
//! | [`ir`] | Synthesis: turning the IR back into gates through the [`CliffordGates`](ir::CliffordGates) / [`Gates`](ir::Gates) sinks |
//! | [`architecture`] | Hardware connectivity graphs that constrain where two-qubit gates may go |
//!
//! # Example
//!
//! Import a small circuit and synthesize it into a gate sink. A sink is any type implementing
//! [`CliffordGates`](ir::CliffordGates) and [`Gates`](ir::Gates) — here one that just counts.
//!
//! ```
//! use synir::import::{import_circuit, ImportedGate};
//! use synir::ir::pauli_exponential::PauliExponentialSynthesizer;
//! use synir::ir::{CliffordGates, Gates, Synthesizer};
//! use synir::IndexType;
//!
//! #[derive(Default)]
//! struct GateCounter {
//!     two_qubit: usize,
//!     rotations: usize,
//! }
//!
//! impl CliffordGates for GateCounter {
//!     fn s(&mut self, _t: IndexType) {}
//!     fn v(&mut self, _t: IndexType) {}
//!     fn s_dgr(&mut self, _t: IndexType) {}
//!     fn v_dgr(&mut self, _t: IndexType) {}
//!     fn x(&mut self, _t: IndexType) {}
//!     fn y(&mut self, _t: IndexType) {}
//!     fn z(&mut self, _t: IndexType) {}
//!     fn h(&mut self, _t: IndexType) {}
//!     fn cx(&mut self, _c: IndexType, _t: IndexType) { self.two_qubit += 1; }
//!     fn cz(&mut self, _c: IndexType, _t: IndexType) { self.two_qubit += 1; }
//! }
//!
//! impl Gates for GateCounter {
//!     fn rx(&mut self, _t: IndexType, _angle: f64) { self.rotations += 1; }
//!     fn ry(&mut self, _t: IndexType, _angle: f64) { self.rotations += 1; }
//!     fn rz(&mut self, _t: IndexType, _angle: f64) { self.rotations += 1; }
//! }
//!
//! // Gates are given in execution order; `import_circuit` reverses internally.
//! let pe = import_circuit(2, &[
//!     ImportedGate::H(0),
//!     ImportedGate::CX(0, 1),
//!     ImportedGate::Rz(1, 0.3),   // non-Clifford: becomes a gadget
//!     ImportedGate::CX(0, 1),
//!     ImportedGate::H(0),
//! ]);
//!
//! let mut circuit = GateCounter::default();
//! PauliExponentialSynthesizer::default().synthesize(pe, &mut circuit);
//!
//! // The single non-Clifford angle survives import as exactly one rotation.
//! assert_eq!(circuit.rotations, 1);
//! ```

pub mod architecture;
pub mod data_structures;
pub mod import;
pub mod ir;

/// Index of a qubit, a tableau row, or a gadget within a polynomial.
///
/// Every index in this crate is 0-based and dense: an `n`-qubit object uses exactly `0..n`.
pub type IndexType = usize;
