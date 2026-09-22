//! Dense unitary reference simulator for tests.
//!
//! This module is deliberately **independent ground truth**: every gate matrix below is written out
//! literally, matching the definitions used by `qiskit.circuit.library` / `qiskit.quantum_info`.
//! Nothing here may be derived from `synir`'s own tableau or propagation code — the whole point is
//! to have a second, external opinion about what a circuit means.
//!
//! # Conventions
//!
//! * **Qubit ordering is little-endian**, exactly like `qiskit.quantum_info.Operator`: bit `k` of a
//!   basis index is the state of qubit `k`. So for two qubits the basis order is
//!   `|q1 q0> = |00>, |01>, |10>, |11>` and a gate on qubit 0 acts on the *low* bit.
//! * **Matrices are row-major**, `data[row * dim + col]`.
//! * **Gates are applied by left-multiplication**, `U <- G * U`. Folding a gate list in *time*
//!   order (first gate executed first) therefore yields the circuit unitary `Gk * ... * G1`.

#![allow(dead_code)]

use num_complex::Complex64 as C;
use std::f64::consts::FRAC_1_SQRT_2;

use synir::data_structures::{PauliLetter, PauliString};

/// Tolerance used by every approximate comparison in the test suite.
pub const TOL: f64 = 1e-9;

fn c(re: f64, im: f64) -> C {
    C::new(re, im)
}

pub const ZERO: C = C::new(0.0, 0.0);
pub const ONE: C = C::new(1.0, 0.0);
pub const I: C = C::new(0.0, 1.0);

/// A dense `2^n x 2^n` unitary.
#[derive(Clone, Debug)]
pub struct Unitary {
    num_qubits: usize,
    dim: usize,
    data: Vec<C>,
}

impl Unitary {
    pub fn identity(num_qubits: usize) -> Self {
        let dim = 1usize << num_qubits;
        let mut data = vec![ZERO; dim * dim];
        for i in 0..dim {
            data[i * dim + i] = ONE;
        }
        Unitary {
            num_qubits,
            dim,
            data,
        }
    }

    pub fn num_qubits(&self) -> usize {
        self.num_qubits
    }

    pub fn dim(&self) -> usize {
        self.dim
    }

    pub fn get(&self, row: usize, col: usize) -> C {
        self.data[row * self.dim + col]
    }

    pub fn set(&mut self, row: usize, col: usize, value: C) {
        let dim = self.dim;
        self.data[row * dim + col] = value;
    }

    /// Left-multiplies by the single qubit gate `m` acting on `target`.
    pub fn apply_1q(&mut self, target: usize, m: [[C; 2]; 2]) {
        assert!(
            target < self.num_qubits,
            "qubit {target} out of range for a {}-qubit unitary",
            self.num_qubits
        );
        let bit = 1usize << target;
        let dim = self.dim;
        for row in 0..dim {
            if row & bit != 0 {
                continue;
            }
            let row1 = row | bit;
            for col in 0..dim {
                let a = self.data[row * dim + col];
                let b = self.data[row1 * dim + col];
                self.data[row * dim + col] = m[0][0] * a + m[0][1] * b;
                self.data[row1 * dim + col] = m[1][0] * a + m[1][1] * b;
            }
        }
    }

    /// Left-multiplies by `CX(control, target)`.
    ///
    /// `CX` is a basis permutation, so this is a row permutation of the matrix. The map
    /// `f(i) = i ^ (bit_control(i) << target)` is an involution, hence `(CX * U)[i] = U[f(i)]`.
    pub fn apply_cx(&mut self, control: usize, target: usize) {
        assert_ne!(control, target, "CX control and target must differ");
        assert!(control < self.num_qubits && target < self.num_qubits);
        let control_bit = 1usize << control;
        let target_bit = 1usize << target;
        let dim = self.dim;
        for row in 0..dim {
            // Only touch rows where the control is set and the target is clear, and swap the pair.
            if row & control_bit == 0 || row & target_bit != 0 {
                continue;
            }
            let other = row | target_bit;
            for col in 0..dim {
                self.data.swap(row * dim + col, other * dim + col);
            }
        }
    }

    /// Left-multiplies by `CZ(a, b)` (symmetric in its arguments).
    pub fn apply_cz(&mut self, a: usize, b: usize) {
        assert_ne!(a, b, "CZ qubits must differ");
        assert!(a < self.num_qubits && b < self.num_qubits);
        let mask = (1usize << a) | (1usize << b);
        let dim = self.dim;
        for row in 0..dim {
            if row & mask == mask {
                for col in 0..dim {
                    self.data[row * dim + col] = -self.data[row * dim + col];
                }
            }
        }
    }

    /// Ordinary matrix product `self * rhs`.
    pub fn mul(&self, rhs: &Unitary) -> Unitary {
        assert_eq!(self.num_qubits, rhs.num_qubits);
        let dim = self.dim;
        let mut data = vec![ZERO; dim * dim];
        for row in 0..dim {
            for k in 0..dim {
                let left = self.data[row * dim + k];
                if left == ZERO {
                    continue;
                }
                for col in 0..dim {
                    data[row * dim + col] += left * rhs.data[k * dim + col];
                }
            }
        }
        Unitary {
            num_qubits: self.num_qubits,
            dim,
            data,
        }
    }

    pub fn adjoint(&self) -> Unitary {
        let dim = self.dim;
        let mut data = vec![ZERO; dim * dim];
        for row in 0..dim {
            for col in 0..dim {
                data[col * dim + row] = self.data[row * dim + col].conj();
            }
        }
        Unitary {
            num_qubits: self.num_qubits,
            dim,
            data,
        }
    }

    pub fn scale(&self, factor: C) -> Unitary {
        Unitary {
            num_qubits: self.num_qubits,
            dim: self.dim,
            data: self.data.iter().map(|v| v * factor).collect(),
        }
    }

    pub fn add(&self, rhs: &Unitary) -> Unitary {
        assert_eq!(self.num_qubits, rhs.num_qubits);
        Unitary {
            num_qubits: self.num_qubits,
            dim: self.dim,
            data: self
                .data
                .iter()
                .zip(rhs.data.iter())
                .map(|(a, b)| a + b)
                .collect(),
        }
    }

    /// The analogue of `qiskit.quantum_info.Operator.equiv`: equality up to a global phase.
    pub fn equiv_up_to_global_phase(&self, other: &Unitary, tol: f64) -> bool {
        if self.num_qubits != other.num_qubits {
            return false;
        }
        let anchor = self
            .data
            .iter()
            .zip(other.data.iter())
            .find(|(a, _)| a.norm() > tol);
        let phase = match anchor {
            // Both all-zero is not a thing for unitaries, but be defensive.
            None => return other.data.iter().all(|b| b.norm() <= tol),
            Some((a, b)) => b / a,
        };
        if (phase.norm() - 1.0).abs() > tol {
            return false;
        }
        self.data
            .iter()
            .zip(other.data.iter())
            .all(|(a, b)| (a * phase - b).norm() <= tol)
    }

    /// Exact equality (no phase freedom), used where the phase is part of the contract.
    pub fn approx_eq(&self, other: &Unitary, tol: f64) -> bool {
        self.num_qubits == other.num_qubits
            && self
                .data
                .iter()
                .zip(other.data.iter())
                .all(|(a, b)| (a - b).norm() <= tol)
    }

    /// The global phase relating `self` to `other`, i.e. `phase` with `other == phase * self`.
    pub fn global_phase_to(&self, other: &Unitary, tol: f64) -> Option<C> {
        let (a, b) = self
            .data
            .iter()
            .zip(other.data.iter())
            .find(|(a, _)| a.norm() > tol)?;
        Some(b / a)
    }

    /// Dense matrix of a Pauli word. Qubit `i` of the word is `pauli.pauli(i)`, matching
    /// `PauliString::from_text` where the *leftmost* character is qubit 0. Note this is the
    /// reverse of Qiskit's `Pauli("...")` string order.
    pub fn from_pauli_string(pauli: &PauliString) -> Unitary {
        let num_qubits = pauli.len();
        let mut u = Unitary::identity(num_qubits);
        for q in 0..num_qubits {
            match pauli.pauli(q) {
                PauliLetter::I => {}
                PauliLetter::X => u.apply_1q(q, x_matrix()),
                PauliLetter::Y => u.apply_1q(q, y_matrix()),
                PauliLetter::Z => u.apply_1q(q, z_matrix()),
            }
        }
        u
    }

    /// `exp(-i * theta / 2 * P) = cos(theta/2) * I - i * sin(theta/2) * P`.
    ///
    /// This is the same convention as Qiskit's `PauliEvolutionGate` / `RZGate`: a `Z` gadget with
    /// angle `theta` equals `rz(theta)` on that qubit.
    pub fn pauli_rotation(pauli: &PauliString, theta: f64) -> Unitary {
        let half = theta / 2.0;
        let identity = Unitary::identity(pauli.len()).scale(c(half.cos(), 0.0));
        let rotation = Unitary::from_pauli_string(pauli).scale(c(0.0, -half.sin()));
        identity.add(&rotation)
    }

    /// Dense matrix of Qiskit's `PermutationGate(pattern)`.
    ///
    /// Qiskit defines `pattern[k] = m` to mean "qubit `m` is permuted to position `k`", so the
    /// output basis bit `k` is the input basis bit `pattern[k]`.
    pub fn permutation(pattern: &[usize], num_qubits: usize) -> Unitary {
        assert_eq!(pattern.len(), num_qubits);
        let mut sorted = pattern.to_vec();
        sorted.sort_unstable();
        assert_eq!(
            sorted,
            (0..num_qubits).collect::<Vec<_>>(),
            "permutation pattern {pattern:?} is not a permutation of 0..{num_qubits}"
        );
        let dim = 1usize << num_qubits;
        let mut data = vec![ZERO; dim * dim];
        for input in 0..dim {
            let mut output = 0usize;
            for k in 0..num_qubits {
                if input >> pattern[k] & 1 == 1 {
                    output |= 1 << k;
                }
            }
            data[output * dim + input] = ONE;
        }
        Unitary {
            num_qubits,
            dim,
            data,
        }
    }

    /// The inverse of a permutation pattern, i.e. `inverse[pattern[k]] == k`.
    pub fn invert_pattern(pattern: &[usize]) -> Vec<usize> {
        let mut inverse = vec![0usize; pattern.len()];
        for (position, &qubit) in pattern.iter().enumerate() {
            inverse[qubit] = position;
        }
        inverse
    }
}

// ---------------------------------------------------------------------------
// Gate matrices, transcribed from Qiskit's definitions.
// ---------------------------------------------------------------------------

pub fn x_matrix() -> [[C; 2]; 2] {
    [[ZERO, ONE], [ONE, ZERO]]
}

pub fn y_matrix() -> [[C; 2]; 2] {
    [[ZERO, -I], [I, ZERO]]
}

pub fn z_matrix() -> [[C; 2]; 2] {
    [[ONE, ZERO], [ZERO, -ONE]]
}

pub fn h_matrix() -> [[C; 2]; 2] {
    let s = c(FRAC_1_SQRT_2, 0.0);
    [[s, s], [s, -s]]
}

pub fn s_matrix() -> [[C; 2]; 2] {
    [[ONE, ZERO], [ZERO, I]]
}

pub fn sdg_matrix() -> [[C; 2]; 2] {
    [[ONE, ZERO], [ZERO, -I]]
}

/// Qiskit's `SXGate`, including its global phase: `0.5 * [[1+i, 1-i], [1-i, 1+i]]`.
pub fn sx_matrix() -> [[C; 2]; 2] {
    [[c(0.5, 0.5), c(0.5, -0.5)], [c(0.5, -0.5), c(0.5, 0.5)]]
}

/// Qiskit's `SXdgGate`: `0.5 * [[1-i, 1+i], [1+i, 1-i]]`.
pub fn sxdg_matrix() -> [[C; 2]; 2] {
    [[c(0.5, -0.5), c(0.5, 0.5)], [c(0.5, 0.5), c(0.5, -0.5)]]
}

pub fn rx_matrix(theta: f64) -> [[C; 2]; 2] {
    let (sin, cos) = (theta / 2.0).sin_cos();
    [[c(cos, 0.0), c(0.0, -sin)], [c(0.0, -sin), c(cos, 0.0)]]
}

pub fn ry_matrix(theta: f64) -> [[C; 2]; 2] {
    let (sin, cos) = (theta / 2.0).sin_cos();
    [[c(cos, 0.0), c(-sin, 0.0)], [c(sin, 0.0), c(cos, 0.0)]]
}

pub fn rz_matrix(theta: f64) -> [[C; 2]; 2] {
    let half = theta / 2.0;
    [
        [c(half.cos(), -half.sin()), ZERO],
        [ZERO, c(half.cos(), half.sin())],
    ]
}
