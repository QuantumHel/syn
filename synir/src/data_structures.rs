//! The intermediate representation itself.
//!
//! Four types, layered from the bottom up:
//!
//! * [`PauliString`] — a tensor product of Pauli letters, stored as a pair of bit vectors.
//! * [`PauliPolynomial`] — a bag of Pauli *gadgets* `exp(-i * theta/2 * P)` with their angles.
//! * [`CliffordTableau`] — a Clifford operator, as the images of the Pauli generators.
//! * [`PauliExponential`] — the full normal form `U = T * U_pp`: a run of polynomials plus a
//!   trailing tableau.
//!
//! # Propagating Cliffords
//!
//! A Clifford `C` is absorbed into the IR by *conjugating* what it meets rather than being stored
//! as a gate. That is what [`PropagateClifford`] does, and both `PauliPolynomial` and
//! `CliffordTableau` implement it. Note the two implementations act in opposite directions, by
//! design:
//!
//! * on a tableau it **appends**, `T <- C * T`;
//! * on a polynomial it conjugates the gadgets, `Q -> C^dag Q C`.
//!
//! So a caller that wants both sides to describe the same rewrite hands one of them the adjoint
//! gate. [`crate::import`] and [`PauliExponential`] each document which way they go and why.
//!
//! [`MaskedPropagateClifford`] is the same operation restricted to a subset of a polynomial's
//! gadgets, which is what lets a synthesizer work on part of a polynomial without disturbing the
//! rest.

use std::fmt;

use crate::IndexType;
use bitvec::vec::BitVec;

pub mod angle;
mod clifford_tableau;
pub mod pauli_exponential;
mod pauli_polynomial;
mod pauli_string;

pub use angle::Angle;
pub use clifford_tableau::CliffordTableau;
use itertools::Itertools;
pub use pauli_exponential::PauliExponential;
pub use pauli_polynomial::PauliPolynomial;
pub use pauli_string::PauliString;

/// An IR object that can be inverted.
///
/// Implemented by [`CliffordTableau`]; required by
/// [`AdjointSynthesizer`](crate::ir::AdjointSynthesizer), which synthesizes the *inverse* of what
/// it is given.
pub trait HasAdjoint {
    /// The inverse operation, `U^dag`.
    fn adjoint(&self) -> Self;
}
/// Absorbs a Clifford gate into an IR object by conjugation.
///
/// Only [`cx`](Self::cx), [`s`](Self::s) and [`v`](Self::v) are primitive; every other gate has a
/// default expansion in terms of those three, so an implementor writes three methods. The
/// expansions are exact, not approximate — e.g. `h = s . v . s` and `x = v . v`.
///
/// Each method returns `&mut Self` so calls chain: `tableau.h(0).cx(0, 1).s(1);`.
///
/// **Direction is implementation-defined.** On [`CliffordTableau`] these *append*
/// (`T <- G * T`); on [`PauliPolynomial`] they conjugate the gadgets (`Q -> G^dag Q G`). See the
/// module documentation.
pub trait PropagateClifford
where
    Self: Sized,
{
    /// Controlled-X from `control` onto `target`.
    fn cx(&mut self, control: IndexType, target: IndexType) -> &mut Self;
    /// Phase gate `S = diag(1, i)`, a quarter turn about Z.
    fn s(&mut self, target: IndexType) -> &mut Self;
    /// `V = sqrt(X)`, a quarter turn about X.
    fn v(&mut self, target: IndexType) -> &mut Self;

    /// `S^dag`, as `Z . S`.
    fn s_dgr(&mut self, target: IndexType) -> &mut Self {
        self.z(target).s(target)
    }

    /// `V^dag`, as `X . V`.
    fn v_dgr(&mut self, target: IndexType) -> &mut Self {
        self.x(target).v(target)
    }

    /// Pauli `X`, as `V . V`.
    fn x(&mut self, target: IndexType) -> &mut Self {
        self.v(target).v(target)
    }

    /// Pauli `Y`, as `S . X . S^dag`.
    fn y(&mut self, target: IndexType) -> &mut Self {
        self.s_dgr(target).x(target).s(target)
    }

    /// Pauli `Z`, as `S . S`.
    fn z(&mut self, target: IndexType) -> &mut Self {
        self.s(target).s(target)
    }

    /// Hadamard, as `S . V . S`.
    fn h(&mut self, target: IndexType) -> &mut Self {
        self.s(target).v(target).s(target)
    }

    /// Controlled-Z, as a `CX` conjugated by `H` on the target.
    fn cz(&mut self, control: IndexType, target: IndexType) -> &mut Self {
        self.h(target);
        self.cx(control, target);
        self.h(target)
    }

    /// Exchanges two qubits, as the usual three alternating `CX`s.
    fn swap(&mut self, control: IndexType, target: IndexType) -> &mut Self {
        self.cx(control, target)
            .cx(target, control)
            .cx(control, target)
    }

    /// Absorbs a whole Pauli gadget `exp(-i * theta/2 * P)`, provided its angle is Clifford.
    ///
    /// The gadget is decomposed the textbook way: basis-change each non-identity leg into `Z`,
    /// fold the legs together with a CNOT ladder, apply a single `Z` rotation on the last leg,
    /// then undo the ladder and the basis changes. Because the angle is a multiple of `pi/2`, that
    /// rotation is itself Clifford (`S`, `Z` or `S^dag`) and the whole gadget stays inside the
    /// group.
    ///
    /// # Errors
    ///
    /// Returns `Err` if `angle` is not a multiple of `pi/2` — a non-Clifford gadget cannot be
    /// absorbed and has to stay in a [`PauliPolynomial`].
    fn propate_gadget(&mut self, gadget: PauliString, angle: Angle) -> Result<(), String> {
        let pi4rotations = angle.to_pi4_rotation();
        let pi2rotations = match pi4rotations {
            Err(_) => Err(format!(
                "Cannot compose Clifford tableau with non-Clifford angle: {}",
                angle
            )),
            Ok(n) => {
                if n % 2 == 1 {
                    Err(format!("Cannot compose Clifford tableau with non-Clifford angle: {} pi/4 rotations", n))
                } else {
                    Ok(n / 2)
                }
            }
        }?;
        let mut leg_numbers = Vec::with_capacity(gadget.len());
        for i in 0..gadget.len() {
            match gadget.pauli(i) {
                PauliLetter::I => {}
                PauliLetter::X => {
                    self.h(i);
                    leg_numbers.push(i);
                }
                PauliLetter::Y => {
                    self.v(i);
                    leg_numbers.push(i);
                }
                PauliLetter::Z => {
                    leg_numbers.push(i);
                }
            }
        }

        for (control, target) in leg_numbers.iter().tuple_windows() {
            self.cx(*control, *target);
        }
        match pi2rotations {
            0 => {}
            1 => {
                let target = *leg_numbers.last().unwrap();
                self.s(target);
            }
            2 => {
                let target = *leg_numbers.last().unwrap();
                self.z(target);
            }
            3 => {
                let target = *leg_numbers.last().unwrap();
                self.s_dgr(target);
            }
            _ => unreachable!(),
        }

        for (control, target) in leg_numbers
            .iter()
            .tuple_windows()
            .collect_vec()
            .iter()
            .rev()
        {
            self.cx(**control, **target);
        }
        for i in 0..gadget.len() {
            match gadget.pauli(i) {
                PauliLetter::I => {}
                PauliLetter::X => {
                    self.h(i);
                }
                PauliLetter::Y => {
                    self.v_dgr(i);
                }
                PauliLetter::Z => {}
            }
        }
        Ok(())
    }
}

/// [`PropagateClifford`] restricted to a subset of a polynomial's gadgets.
///
/// `mask` is indexed by *gadget*, not by qubit: bit `i` selects gadget `i` of the polynomial, so
/// the mask has [`PauliPolynomial::length`] bits. Gadgets whose bit is clear are left untouched.
///
/// This is what lets a synthesizer diagonalize one block of a polynomial while the remaining
/// gadgets keep their original form — the partitioning synthesizers in
/// [`crate::ir::pauli_polynomial`] rely on it.
///
/// As with `PropagateClifford`, only the three primitives are required; the rest expand in terms
/// of them.
pub trait MaskedPropagateClifford
where
    Self: Sized,
{
    /// Controlled-X, applied only to the gadgets selected by `mask`.
    fn masked_cx(&mut self, control: IndexType, target: IndexType, mask: &BitVec) -> &mut Self;
    /// Phase gate `S`, applied only to the gadgets selected by `mask`.
    fn masked_s(&mut self, target: IndexType, mask: &BitVec) -> &mut Self;
    /// `V = sqrt(X)`, applied only to the gadgets selected by `mask`.
    fn masked_v(&mut self, target: IndexType, mask: &BitVec) -> &mut Self;

    fn masked_s_dgr(&mut self, target: IndexType, mask: &BitVec) -> &mut Self {
        self.masked_z(target, mask).masked_s(target, mask)
    }

    fn masked_v_dgr(&mut self, target: IndexType, mask: &BitVec) -> &mut Self {
        self.masked_x(target, mask).masked_v(target, mask)
    }

    fn masked_x(&mut self, target: IndexType, mask: &BitVec) -> &mut Self {
        self.masked_v(target, mask).masked_v(target, mask)
    }

    fn masked_y(&mut self, target: IndexType, mask: &BitVec) -> &mut Self {
        self.masked_s_dgr(target, mask)
            .masked_x(target, mask)
            .masked_s(target, mask)
    }

    fn masked_z(&mut self, target: IndexType, mask: &BitVec) -> &mut Self {
        self.masked_s(target, mask).masked_s(target, mask)
    }

    fn masked_h(&mut self, target: IndexType, mask: &BitVec) -> &mut Self {
        self.masked_s(target, mask)
            .masked_v(target, mask)
            .masked_s(target, mask)
    }

    fn masked_cz(&mut self, control: IndexType, target: IndexType, mask: &BitVec) -> &mut Self {
        self.masked_h(target, mask);
        self.masked_cx(control, target, mask);
        self.masked_h(target, mask)
    }
}

/// A single-qubit Pauli matrix, ignoring phase.
///
/// Stored elsewhere as a pair of bits `(x, z)`; this enum is the readable view of that pair.
/// `Y` is `(1, 1)` — i.e. `Y ~ XZ` up to the phase this type deliberately does not track. Signs
/// live separately, on [`CliffordTableau`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PauliLetter {
    /// Identity, `(x, z) = (0, 0)`.
    I,
    /// Pauli X, `(x, z) = (1, 0)`.
    X,
    /// Pauli Y, `(x, z) = (1, 1)`.
    Y,
    /// Pauli Z, `(x, z) = (0, 1)`.
    Z,
}

impl PauliLetter {
    /// The letter with the given `x` and `z` bits.
    pub fn new(x: bool, z: bool) -> Self {
        match (x, z) {
            (false, false) => PauliLetter::I,
            (true, false) => PauliLetter::X,
            (true, true) => PauliLetter::Y,
            (false, true) => PauliLetter::Z,
        }
    }

    /// Parses one of `'I'`, `'X'`, `'Y'`, `'Z'`.
    ///
    /// # Panics
    ///
    /// Panics on any other character. Parsing is uppercase-only.
    pub fn from_char(c: char) -> PauliLetter {
        let (x, z) = match c {
            'I' => (false, false),
            'X' => (true, false),
            'Y' => (true, true),
            'Z' => (false, true),
            _ => panic!("Unknown Pauli letter {c}"),
        };
        PauliLetter::new(x, z)
    }
}

impl fmt::Display for PauliLetter {
    // This trait requires `fmt` with this exact signature.
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let s = match self {
            PauliLetter::I => "I",
            PauliLetter::X => "X",
            PauliLetter::Y => "Y",
            PauliLetter::Z => "Z",
        };
        write!(f, "{}", s)
    }
}
