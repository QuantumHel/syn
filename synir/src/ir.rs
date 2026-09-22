//! Synthesis: turning the IR back into a gate sequence.
//!
//! # Sinks
//!
//! Synthesis does not build a circuit type. It *emits* gates into a caller-supplied sink — any
//! type implementing [`CliffordGates`] (and [`Gates`], for non-Clifford rotations). That is how
//! `synpy` hands results straight to a Qiskit `QuantumCircuit` without synIR knowing Qiskit
//! exists, and how the test suite collects gates into a mock for unitary comparison.
//!
//! # The two synthesizer traits
//!
//! [`Synthesizer`] emits a circuit for the IR it is given. [`AdjointSynthesizer`] emits one for
//! the **inverse**.
//!
//! The inverse is the more natural direction for tableau synthesis, which works by *reducing* a
//! tableau to the identity: the gates that reduce `T` compose to `T^dag`. Rather than reverse and
//! dagger that sequence, the algorithms implement `AdjointSynthesizer` directly, and a blanket
//! impl in [`clifford_tableau`] upgrades any of them to a `Synthesizer` by inverting the input
//! first.
//!
//! # Available strategies
//!
//! | IR | Strategy | Notes |
//! |---|---|---|
//! | [`CliffordTableau`](crate::data_structures::CliffordTableau) | [`NaiveCliffordSynthesizer`](clifford_tableau::NaiveCliffordSynthesizer) | Ignores connectivity |
//! | | [`PermRowColCliffordSynthesizer`](clifford_tableau::PermRowColCliffordSynthesizer) | Connectivity-aware, picks its own qubit order |
//! | | [`CallbackCliffordSynthesizer`](clifford_tableau::CallbackCliffordSynthesizer) | Pivot choice delegated to a closure |
//! | [`PauliPolynomial`](crate::data_structures::PauliPolynomial) | [`NaivePauliPolynomialSynthesizer`](pauli_polynomial::NaivePauliPolynomialSynthesizer) | One gadget at a time |
//! | | [`PSGSPauliPolynomialSynthesizer`](pauli_polynomial::psgs::PSGSPauliPolynomialSynthesizer) | Recursive partitioning, shares CNOT structure |
//! | [`PauliExponential`](crate::data_structures::PauliExponential) | [`PauliExponentialSynthesizer`](pauli_exponential::PauliExponentialSynthesizer) | Drives both halves; pick each via a strategy enum |

use crate::{data_structures::HasAdjoint, IndexType};

pub mod clifford_tableau;
pub mod pauli_exponential;
pub mod pauli_polynomial;

/// A sink that Clifford gates can be emitted into, in circuit order.
///
/// Implement this on whatever circuit representation you want synthesis to produce. Unlike
/// [`PropagateClifford`](crate::data_structures::PropagateClifford), every gate is required —
/// a sink that silently expanded `H` into `S . V . S` would emit a needlessly long circuit.
pub trait CliffordGates {
    /// Phase gate `S = diag(1, i)`.
    fn s(&mut self, target: IndexType);
    /// `V = sqrt(X)`.
    fn v(&mut self, target: IndexType);
    /// `S^dag`.
    fn s_dgr(&mut self, target: IndexType);
    /// `V^dag`.
    fn v_dgr(&mut self, target: IndexType);
    /// Pauli `X`.
    fn x(&mut self, target: IndexType);
    /// Pauli `Y`.
    fn y(&mut self, target: IndexType);
    /// Pauli `Z`.
    fn z(&mut self, target: IndexType);
    /// Hadamard.
    fn h(&mut self, target: IndexType);
    /// Controlled-X.
    fn cx(&mut self, control: IndexType, target: IndexType);
    /// Controlled-Z.
    fn cz(&mut self, control: IndexType, target: IndexType);
    /// Exchanges two qubits.
    ///
    /// Defaults to three alternating `CX`s. Override it if the target hardware has a native swap,
    /// or if you would rather track the permutation than pay for the gates.
    fn swap(&mut self, control: IndexType, target: IndexType) {
        self.cx(control, target);
        self.cx(target, control);
        self.cx(control, target);
    }
    /// Emits the qubit relabelling that synthesis has left over.
    ///
    /// Tableau synthesis stops once the residual tableau is a bare permutation (see
    /// [`CliffordTableau::get_permutation`](crate::data_structures::CliffordTableau::get_permutation)),
    /// and this closes the gap. `permutation[physical] = logical`.
    ///
    /// The default expands it into swaps, which costs real two-qubit gates. A sink that can simply
    /// *record* the permutation and relabel downstream — as a transpiler usually can — should
    /// override this and emit nothing.
    fn add_final_permutation(&mut self, permutation: Vec<IndexType>) {
        let mut perm = permutation.clone();
        for i in 0..permutation.len() {
            // TODO: Rewrite so that it doesn't need to create 3 clones
            let tmp_perm = perm.clone();
            let j = tmp_perm.iter().find(|&x| *x == i).unwrap();
            if i != *j {
                // Swap the qubits
                self.cx(i, *j);
                self.cx(*j, i);
                self.cx(i, *j);
                // Update the permutation
                perm.swap(i, *j);
            }
        }
    }
}

/// A sink that arbitrary-angle rotations can be emitted into.
///
/// The non-Clifford companion to [`CliffordGates`]. Angles are in radians.
pub trait Gates {
    /// Rotation about X by `angle` radians.
    fn rx(&mut self, target: IndexType, angle: f64);
    /// Rotation about Y by `angle` radians.
    fn ry(&mut self, target: IndexType, angle: f64);
    /// Rotation about Z by `angle` radians.
    fn rz(&mut self, target: IndexType, angle: f64);
}

/// Emits a gate sequence implementing the **inverse** of the given IR.
///
/// Tableau algorithms reduce a tableau to the identity, and the gates that do the reducing
/// compose to `ir^dag` — so this, not [`Synthesizer`], is what they naturally implement. The
/// blanket impl in [`clifford_tableau`] derives the forward direction by inverting the input.
///
/// `Returns` is the residual IR left after synthesis; for the tableau synthesizers it is the
/// leftover permutation tableau.
pub trait AdjointSynthesizer<From, To, Returns = ()>
where
    From: HasAdjoint,
{
    /// Emits `ir^dag` into `external_repr`, returning whatever residual the algorithm defines.
    fn synthesize_adjoint(&mut self, ir: From, external_repr: &mut To) -> Returns;
}

/// Emits a gate sequence implementing the given IR.
///
/// `To` is the sink — see [`CliffordGates`]. `Returns` is the residual IR, if the algorithm leaves
/// one; it is `()` for a synthesizer that consumes its input entirely.
pub trait Synthesizer<From, To, Returns = ()> {
    /// Emits `ir` into `repr`.
    fn synthesize(&mut self, ir: From, repr: &mut To) -> Returns;
}
