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

pub trait HasAdjoint {
    fn adjoint(&self) -> Self;
}
pub trait PropagateClifford
where
    Self: Sized,
{
    fn cx(&mut self, control: IndexType, target: IndexType) -> &mut Self;
    fn s(&mut self, target: IndexType) -> &mut Self;
    fn v(&mut self, target: IndexType) -> &mut Self;

    fn s_dgr(&mut self, target: IndexType) -> &mut Self {
        self.z(target).s(target)
    }

    fn v_dgr(&mut self, target: IndexType) -> &mut Self {
        self.x(target).v(target)
    }

    fn x(&mut self, target: IndexType) -> &mut Self {
        self.v(target).v(target)
    }

    fn y(&mut self, target: IndexType) -> &mut Self {
        self.s_dgr(target).x(target).s(target)
    }

    fn z(&mut self, target: IndexType) -> &mut Self {
        self.s(target).s(target)
    }

    fn h(&mut self, target: IndexType) -> &mut Self {
        self.s(target).v(target).s(target)
    }

    fn cz(&mut self, control: IndexType, target: IndexType) -> &mut Self {
        self.h(target);
        self.cx(control, target);
        self.h(target)
    }

    fn swap(&mut self, control: IndexType, target: IndexType) -> &mut Self {
        self.cx(control, target)
            .cx(target, control)
            .cx(control, target)
    }

    /// Composes a gadget onto a Clifford tableau if the angle is Clifford
    /// Decomposes the Pauli gadget by performing naive decomposition into mapping to Z legs, CNOT walls and Z-rotations
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

pub trait MaskedPropagateClifford
where
    Self: Sized,
{
    fn masked_cx(&mut self, control: IndexType, target: IndexType, mask: &BitVec) -> &mut Self;
    fn masked_s(&mut self, target: IndexType, mask: &BitVec) -> &mut Self;
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

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PauliLetter {
    I,
    X,
    Y,
    Z,
}

impl PauliLetter {
    pub fn new(x: bool, z: bool) -> Self {
        match (x, z) {
            (false, false) => PauliLetter::I,
            (true, false) => PauliLetter::X,
            (true, true) => PauliLetter::Y,
            (false, true) => PauliLetter::Z,
        }
    }

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
