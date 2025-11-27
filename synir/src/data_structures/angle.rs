use std::ops::{AddAssign, SubAssign, Add, Sub};

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Angle {
    Angle(f64),
    Pi4Rotations(usize),
}

impl Angle {
    pub fn from_angle(rad: f64) -> Self {
        Angle::Angle(rad)
    }

    pub fn from_angles(angles: &[f64]) -> Vec<Self> {
        angles
            .into_iter()
            .map(|rad| Angle::from_angle(*rad))
            .collect()
    }

    pub fn from_pi4_rotations(n: usize) -> Self {
        Angle::Pi4Rotations(n % 8)
    }

    pub fn forpi4_rotations(ns: &[usize]) -> Vec<Self> {
        ns.into_iter()
            .map(|n| Angle::from_pi4_rotations(*n))
            .collect()
    }

    pub fn to_radians(&self) -> f64 {
        match self {
            Angle::Angle(rad) => *rad,
            Angle::Pi4Rotations(n) => (*n as f64) * (std::f64::consts::FRAC_PI_4),
        }
    }

    pub fn flip(&mut self) {
        match self {
            Angle::Angle(rad) => *rad = -*rad,
            Angle::Pi4Rotations(n) => *n = (8 - *n) % 8,
        }
    }
}

impl AddAssign for Angle {
    fn add_assign(&mut self, other: Self) {
        match (self, other) {
            (Angle::Angle(rad1), Angle::Angle(rad2)) => {
                *rad1 += rad2;
            }
            (Angle::Pi4Rotations(n1), Angle::Pi4Rotations(n2)) => {
                *n1 = (*n1 + n2) % 8;
            }
            _ => panic!("Cannot add different types of Angles"),
        }
    }
}

impl SubAssign for Angle {
    fn sub_assign(&mut self, other: Self) {
        match (self, other) {
            (Angle::Angle(rad1), Angle::Angle(rad2)) => {
                *rad1 -= rad2;
            }
            (Angle::Pi4Rotations(n1), Angle::Pi4Rotations(n2)) => {
                *n1 = (*n1 + (8 - n2)) % 8;
            }
            _ => panic!("Cannot subtract different types of Angles"),
        }
    }
}

impl Add for Angle {
    type Output = Angle;

    fn add(self, other: Angle) -> Angle {
        match (self, other) {
            (Angle::Angle(rad1), Angle::Angle(rad2)) => Angle::Angle(rad1 + rad2),
            (Angle::Pi4Rotations(n1), Angle::Pi4Rotations(n2)) => {
                Angle::Pi4Rotations((n1 + n2) % 8)
            }
            _ => panic!("Cannot add different types of Angles"),
        }
    }
}

impl Sub for Angle {
    type Output = Angle;

    fn sub(self, other: Angle) -> Angle {
        match (self, other) {
            (Angle::Angle(rad1), Angle::Angle(rad2)) => Angle::Angle(rad1 - rad2),
            (Angle::Pi4Rotations(n1), Angle::Pi4Rotations(n2)) => {
                Angle::Pi4Rotations((n1 + 8 - n2) % 8)
            }
            _ => panic!("Cannot add different types of Angles"),
        }
    }
}
