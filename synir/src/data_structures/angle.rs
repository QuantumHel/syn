use std::{
    f64::{EPSILON, consts::PI}, fmt, ops::{Add, AddAssign, Neg, Sub, SubAssign}
};

use crate::data_structures::angle;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Angle {
    Arbitrary(f64),
    Pi4Rotations(u8),
}

impl Angle {
    pub fn from_angle(rad: f64) -> Self {
        Angle::Arbitrary(rad)
    }

    pub fn from_angles(angles: &[f64]) -> Vec<Self> {
        angles
            .into_iter()
            .map(|rad| Angle::from_angle(*rad))
            .collect()
    }

    pub fn from_pi4_rotation(n: u8) -> Self {
        Angle::Pi4Rotations(n % 8)
    }

    pub fn from_pi4_rotations(ns: &[u8]) -> Vec<Self> {
        ns.into_iter()
            .map(|n| Angle::from_pi4_rotation(*n))
            .collect()
    }

    pub fn to_radians(&self) -> f64 {
        match self {
            Angle::Arbitrary(rad) => *rad,
            Angle::Pi4Rotations(n) => (*n as f64) * (std::f64::consts::FRAC_PI_4),
        }
    }

    pub fn to_pi4_rotation(&self) -> Result<u8, String>{
        match self {
            Angle::Pi4Rotations(n) => Ok(*n),
            Angle::Arbitrary(rad) => {
                let pi4_rot = rad * std::f64::consts::FRAC_2_PI * 2.;
                if pi4_rot.fract() > 1.0e-14 { // Catch floating point errors
                    Err(format!("Can only cast Angles that are multiples of pi/4. Fraction part is {}", pi4_rot.fract()))
                } else {
                    let mut n = pi4_rot.trunc() as i64;
                    while n < 0 {
                        n += 8;
                    }
                    Ok((n % 8) as u8)
                }
            }
        }
    }

    pub fn flip(&mut self) {
        match self {
            Angle::Arbitrary(rad) => *rad = -*rad,
            Angle::Pi4Rotations(n) => *n = (8 - *n) % 8,
        }
    }

    pub fn is_clifford(&self) -> bool {
        match self.to_pi4_rotation() {
            Ok(n) => n%2 == 0,
            Err(_) => false
        }
    }

}

impl AddAssign for Angle {
    fn add_assign(&mut self, other: Self) {
        match (self, other) {
            (Angle::Arbitrary(rad1), Angle::Arbitrary(rad2)) => {
                *rad1 += rad2;
            }
            (Angle::Pi4Rotations(n1), Angle::Pi4Rotations(n2)) => {
                *n1 = (*n1 + n2) % 8;
            }
            (Angle::Arbitrary(rad1), Angle::Pi4Rotations(n2)) => {
                *rad1 += n2 as f64 * PI / 4.0;
            }
            _ => panic!("Cannot add Arbitrary Angle to Pi4 rotation"),
        }
    }
}

impl SubAssign for Angle {
    fn sub_assign(&mut self, other: Self) {
        match (self, other) {
            (Angle::Arbitrary(rad1), Angle::Arbitrary(rad2)) => {
                *rad1 -= rad2;
            }
            (Angle::Pi4Rotations(n1), Angle::Pi4Rotations(n2)) => {
                *n1 = (*n1 + (8 - n2)) % 8;
            }
            (Angle::Arbitrary(rad1), Angle::Pi4Rotations(n2)) => {
                *rad1 -= n2 as f64 * PI / 4.0;
            }
            _ => panic!("Cannot subtract different types of Angles"),
        }
    }
}

impl Add for Angle {
    type Output = Angle;

    fn add(self, other: Angle) -> Angle {
        match (self, other) {
            (Angle::Arbitrary(rad1), Angle::Arbitrary(rad2)) => Angle::Arbitrary(rad1 + rad2),
            (Angle::Pi4Rotations(n1), Angle::Pi4Rotations(n2)) => {
                Angle::Pi4Rotations((n1 + n2) % 8)
            }
            (Angle::Arbitrary(rad1), Angle::Pi4Rotations(n2)) => {
                Angle::Arbitrary(rad1 + n2 as f64 * PI / 4.0)
            }
            (Angle::Pi4Rotations(n1), Angle::Arbitrary(rad2)) => {
                Angle::Arbitrary(rad2 + n1 as f64 * PI / 4.0)
            }
        }
    }
}

impl Neg for Angle {
    type Output = Angle;

    fn neg(self) -> Angle {
        match self {
            Angle::Arbitrary(rad) => Angle::Arbitrary(-rad),
            Angle::Pi4Rotations(n) => Angle::Pi4Rotations((8 - n) % 8),
        }
    }
}

impl Sub for Angle {
    type Output = Angle;

    fn sub(self, other: Angle) -> Angle {
        match (self, other) {
            (Angle::Arbitrary(rad1), Angle::Arbitrary(rad2)) => Angle::Arbitrary(rad1 - rad2),
            (Angle::Pi4Rotations(n1), Angle::Pi4Rotations(n2)) => {
                Angle::Pi4Rotations((n1 + 8 - n2) % 8)
            }
            (Angle::Arbitrary(rad1), Angle::Pi4Rotations(n2)) => {
                Angle::Arbitrary(rad1 - n2 as f64 * PI / 4.0)
            }
            (Angle::Pi4Rotations(n1), Angle::Arbitrary(rad2)) => {
                Angle::Arbitrary(n1 as f64 * PI / 4.0 - rad2)
            }
        }
    }
}

impl fmt::Display for Angle{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Angle::Arbitrary(rad) => write!(f, "{}", rad),
            Angle::Pi4Rotations(n) => write!(f, "{} PI / 4", n)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn check_angle_approx(angle1: Angle, angle2: Angle) -> bool {
        match (angle1, angle2) {
            (Angle::Arbitrary(a1), Angle::Arbitrary(a2)) => (a1 - a2).abs() < 1e-9,
            (Angle::Pi4Rotations(a1), Angle::Pi4Rotations(a2)) => a1 == a2,
            _ => panic!("Not defined for Arbitrary Angles and Pi4 rotations"),
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
    }

    #[test]
    #[should_panic]
    fn test_angle_bad_mixed_simple_add() {
        let n1 = 0.32;
        let n2 = 2;

        let mut a2 = Angle::from_pi4_rotation(n2);
        let a3 = Angle::from_angle(n1);
        a2 += a3
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
    }

    #[test]
    #[should_panic]
    fn test_angle_bad_mixed_simple_sub() {
        let n1 = 0.32;
        let n2 = 2;

        let mut a2 = Angle::from_pi4_rotation(n2);
        let a3 = Angle::from_angle(n1);
        a2 -= a3
    }

    #[test]
    fn test_angle_to_pi4(){
        for n in 0..16{
            let angle = Angle::Pi4Rotations(n);
            let a_angle = Angle::Arbitrary(angle.to_radians());
            println!("{} {}", angle, a_angle);
            let res = a_angle.to_pi4_rotation();
            assert!(res.is_ok(), "{}", res.err().unwrap());
            println!("{}", res.as_ref().unwrap());
            assert!((n as u8) % 8 == res.unwrap());
        }
        let alt_angle = Angle::Arbitrary(-PI/4.);
        println!("{}", alt_angle);
        let res2 = alt_angle.to_pi4_rotation();
        assert!(res2.is_ok(), "{}", res2.err().unwrap());
        println!("{}", res2.as_ref().unwrap());
        assert!(res2.unwrap() == 7);
    }
}
