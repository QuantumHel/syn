use std::{
    f64::{consts::PI, EPSILON},
    fmt,
    ops::{Add, AddAssign, Neg, Sub, SubAssign},
};

use crate::data_structures::angle;

#[derive(Debug, Clone, Copy, PartialEq)]
enum AngleType {
    // Hiding the Angle type so we can enforce u8 mod 8 and -pi < f64 < pi
    Arbitrary(f64),
    Pi4Rotations(u8),
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Angle(AngleType);

impl Angle {
    pub fn from_angle(rad: f64) -> Self {
        let mut bounded_rad = rad;
        while bounded_rad < -PI {
            bounded_rad += 2. * PI;
        }
        while bounded_rad > PI {
            bounded_rad -= 2. * PI;
        }
        Angle(AngleType::Arbitrary(bounded_rad))
    }

    pub fn from_angles(angles: &[f64]) -> Vec<Self> {
        angles
            .into_iter()
            .map(|rad| Angle::from_angle(*rad))
            .collect()
    }

    pub fn from_pi4_rotation(n: u8) -> Self {
        Angle(AngleType::Pi4Rotations(n % 8))
    }

    pub fn from_pi4_rotations(ns: &[u8]) -> Vec<Self> {
        ns.into_iter()
            .map(|n| Angle::from_pi4_rotation(*n))
            .collect()
    }

    pub fn to_radians(&self) -> f64 {
        match self {
            Angle(AngleType::Arbitrary(rad)) => *rad,
            Angle(AngleType::Pi4Rotations(n)) => (*n as f64) * (std::f64::consts::FRAC_PI_4),
        }
    }

    pub fn to_pi4_rotation(&self) -> Result<u8, String> {
        match self {
            Angle(AngleType::Pi4Rotations(n)) => Ok(*n),
            Angle(AngleType::Arbitrary(rad)) => {
                let pi4_rot = rad * std::f64::consts::FRAC_2_PI * 2.;
                let mut n = pi4_rot.round() as i64;
                let diff = (pi4_rot - n as f64).abs();
                if diff > 0. {
                    Err(format!(
                        "Can only cast Angles that are multiples of pi/4. Fraction part is {}",
                        pi4_rot.fract()
                    ))
                } else {
                    while n < 0 {
                        n += 8;
                    }
                    Ok((n % 8) as u8)
                }
            }
        }
    }

    fn to_pi_rotation(&self) -> f64 {
        self.to_radians() * std::f64::consts::FRAC_1_PI
    }

    pub fn flip(&mut self) {
        match self {
            Angle(AngleType::Arbitrary(rad)) => *rad = -*rad,
            Angle(AngleType::Pi4Rotations(n)) => *n = (8 - *n) % 8,
        }
    }

    pub fn is_clifford(&self) -> bool {
        match self.to_pi4_rotation() {
            Ok(n) => n % 2 == 0,
            Err(_) => false,
        }
    }
}

impl AddAssign for Angle {
    fn add_assign(&mut self, other: Self) {
        *self = *self + other;
    }
}

impl SubAssign for Angle {
    fn sub_assign(&mut self, other: Self) {
        *self = *self - other;
    }
}

impl Add for Angle {
    type Output = Angle;

    fn add(self, other: Angle) -> Angle {
        return match self {
            Angle(AngleType::Arbitrary(rad1)) => Angle::from_angle(rad1 + other.to_radians()),
            Angle(AngleType::Pi4Rotations(n1)) => {
                let maybe_pi4_rot = other.to_pi4_rotation();
                match maybe_pi4_rot {
                    Ok(n2) => Angle::from_pi4_rotation(n1 + n2),
                    Err(_) => Angle::from_angle(self.to_radians() + other.to_radians()),
                }
            }
        };
    }
}

impl Neg for Angle {
    type Output = Angle;

    fn neg(self) -> Angle {
        match self {
            Angle(AngleType::Arbitrary(rad)) => Angle::from_angle(-rad),
            Angle(AngleType::Pi4Rotations(n)) => Angle::from_pi4_rotation(8 - n),
        }
    }
}

impl Sub for Angle {
    type Output = Angle;

    fn sub(self, other: Angle) -> Angle {
        return match self {
            Angle(AngleType::Arbitrary(rad1)) => Angle::from_angle(rad1 - other.to_radians()),
            Angle(AngleType::Pi4Rotations(n1)) => {
                let maybe_pi4_rot = other.to_pi4_rotation();
                match maybe_pi4_rot {
                    Ok(n2) => Angle::from_pi4_rotation((8 + n1) - n2),
                    Err(_) => Angle::from_angle(self.to_radians() - other.to_radians()),
                }
            }
        };
    }
}

impl fmt::Display for Angle {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Angle(AngleType::Arbitrary(_)) => write!(f, "{}*PI", self.to_pi_rotation()),
            Angle(AngleType::Pi4Rotations(n)) => write!(f, "{}*PI/4", n),
        }
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    fn check_angle_approx(angle1: Angle, angle2: Angle) -> bool {
        let diff = angle1 - angle2;
        match diff {
            Angle(AngleType::Arbitrary(rad)) => rad.abs() < 1e-9,
            Angle(AngleType::Pi4Rotations(n)) => n == 0,
        }
    }

    #[test]
    fn test_angle_simple_add() {
        let n1 = 1;
        let n2 = 2;

        let a1 = Angle::from_pi4_rotation(n1);
        let a2 = Angle::from_pi4_rotation(n2);

        assert_eq!(a1 + a2, Angle::from_pi4_rotation(3));

        let mut a3 = Angle::from_pi4_rotation(n1);
        a3 += a2;

        assert_eq!(a3, Angle::from_pi4_rotation(3));
    }

    #[test]
    fn test_angle_overflow_add() {
        let n1 = 5;
        let n2 = 6;

        let a1 = Angle::from_pi4_rotation(n1);
        let a2 = Angle::from_pi4_rotation(n2);

        assert_eq!(a1 + a2, Angle::from_pi4_rotation(3));

        let mut a3 = Angle::from_pi4_rotation(n1);
        a3 += a2;

        assert_eq!(a3, Angle::from_pi4_rotation(3));
    }

    #[test]
    fn test_angle_simple_sub() {
        let n1 = 4;
        let n2 = 2;

        let a1 = Angle::from_pi4_rotation(n1);
        let a2 = Angle::from_pi4_rotation(n2);

        assert_eq!(a1 - a2, Angle::from_pi4_rotation(2));

        let mut a3 = Angle::from_pi4_rotation(n1);
        a3 -= a2;

        assert_eq!(a3, Angle::from_pi4_rotation(2));
    }

    #[test]
    fn test_angle_overflow_sub() {
        let n1 = 2;
        let n2 = 6;

        let a1 = Angle::from_pi4_rotation(n1);
        let a2 = Angle::from_pi4_rotation(n2);

        assert_eq!(a1 - a2, Angle::from_pi4_rotation(4));

        let mut a3 = Angle::from_pi4_rotation(n1);
        a3 -= a2;

        assert_eq!(a3, Angle::from_pi4_rotation(4));
    }

    #[test]
    fn test_angle_float_simple_add() {
        let n1 = 0.32;
        let n2 = 0.64;

        let a1 = Angle::from_angle(n1);
        let a2 = Angle::from_angle(n2);

        let ref_a = Angle::from_angle(0.96);

        assert!(check_angle_approx(a1 + a2, ref_a));

        let mut a3 = Angle::from_angle(n1);
        a3 += a2;

        assert!(check_angle_approx(a3, ref_a));
    }

    #[test]
    fn test_angle_float_simple_sub() {
        let n1 = 0.32;
        let n2 = 0.64;

        let a1 = Angle::from_angle(n1);
        let a2 = Angle::from_angle(n2);

        let ref_a = Angle::from_angle(-0.32);

        assert!(check_angle_approx(a1 - a2, ref_a));

        let mut a3 = Angle::from_angle(n1);
        a3 -= a2;

        assert!(check_angle_approx(a3, ref_a));
    }

    #[test]
    fn test_angle_mixed_simple_add() {
        let n1 = 0.32;
        let n2 = 2;

        let a1 = Angle::from_angle(n1);
        let a2 = Angle::from_pi4_rotation(n2);

        let ref_a = Angle::from_angle(1.8907963268);
        assert!(check_angle_approx(a1 + a2, ref_a));
        assert!(check_angle_approx(a2 + a1, ref_a));

        let mut a3 = Angle::from_angle(n1);
        a3 += a2;

        assert!(check_angle_approx(a3, ref_a));
        let mut a4 = Angle::from_pi4_rotation(n2);
        a4 += a1;
        assert!(check_angle_approx(a4, ref_a));
    }

    #[test]
    fn test_angle_mixed_simple_sub() {
        let n1 = 0.32;
        let n2 = 2;

        let a1 = Angle::from_angle(n1);
        let a2 = Angle::from_pi4_rotation(n2);

        let ref_a1 = Angle::from_angle(-1.2507963268);
        let ref_a2 = Angle::from_angle(1.2507963268);

        assert!(check_angle_approx(a1 - a2, ref_a1));
        assert!(check_angle_approx(a2 - a1, ref_a2));

        let mut a3 = Angle::from_angle(n1);
        a3 -= a2;

        assert!(check_angle_approx(a3, ref_a1));

        let mut a4 = Angle::from_pi4_rotation(n2);
        a4 -= a1;
        assert!(check_angle_approx(a4, ref_a2));
    }

    #[test]
    fn test_angle_to_pi4() {
        for n in 0..16 {
            let angle = Angle::from_pi4_rotation(n);
            let a_angle = Angle::from_angle(angle.to_radians());
            let res = a_angle.to_pi4_rotation();
            match res {
                Ok(val) => assert!((n as u8) % 8 == val),
                // floating point errors
                //Err(msg) => assert!(msg.chars().counts_by(|c| c)[&'0'] > 14)
                Err(msg) => {
                    println!("{}", msg);
                    assert!(false)
                }
            }
        }
        let alt_angle = Angle::from_angle(-PI / 4.);
        let res2 = alt_angle.to_pi4_rotation();
        assert!(res2.is_ok(), "{}", res2.err().unwrap());
        assert!(res2.unwrap() == 7);
    }
}
