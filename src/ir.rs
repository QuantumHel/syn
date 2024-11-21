use crate::IndexType;

pub mod clifford_tableau;
pub mod pauli_exponential;

pub trait CliffordGates {
    fn s(&mut self, target: IndexType);
    fn v(&mut self, target: IndexType);
    fn s_dgr(&mut self, target: IndexType);
    fn v_dgr(&mut self, target: IndexType);
    fn x(&mut self, target: IndexType);
    fn y(&mut self, target: IndexType);
    fn z(&mut self, target: IndexType);
    fn h(&mut self, target: IndexType);
    fn cx(&mut self, control: IndexType, target: IndexType);
    fn cz(&mut self, control: IndexType, target: IndexType);
}

pub trait Gates {
    fn rx(&mut self, target: IndexType, angle: f64);
    fn ry(&mut self, target: IndexType, angle: f64);
    fn rz(&mut self, target: IndexType, angle: f64);
}

impl CliffordGates for &mut CliffordGatesPrinter {
    fn s(&mut self, target: IndexType) {
        self.gates.push(format!("S {}", target));
    }

    fn v(&mut self, target: IndexType) {
        self.gates.push(format!("V {}", target));
    }

    fn s_dgr(&mut self, target: IndexType) {
        self.gates.push(format!("S† {}", target));
    }

    fn v_dgr(&mut self, target: IndexType) {
        self.gates.push(format!("V† {}", target));
    }

    fn x(&mut self, target: IndexType) {
        self.gates.push(format!("X {}", target));
    }

    fn y(&mut self, target: IndexType) {
        self.gates.push(format!("Y {}", target));
    }

    fn z(&mut self, target: IndexType) {
        self.gates.push(format!("Z {}", target));
    }

    fn h(&mut self, target: IndexType) {
        self.gates.push(format!("H {}", target));
    }

    fn cx(&mut self, control: IndexType, target: IndexType) {
        self.gates.push(format!("CX {} {}", control, target));
    }

    fn cz(&mut self, control: IndexType, target: IndexType) {
        self.gates.push(format!("CZ {} {}", control, target));
    }
}