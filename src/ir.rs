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

pub struct CliffordGatesPrinter {
    pub size: usize,
    pub gates: Vec<String>,
}

impl CliffordGatesPrinter {
    pub fn new(size: usize) -> Self {
        CliffordGatesPrinter {
            size,
            gates: vec![
                "OPENQASM 2.0;".to_string(),
                "include \"qelib1.inc\";".to_string(),
                format!("qreg q[{}];", size),
            ],
        }
    }
}

impl CliffordGates for CliffordGatesPrinter {
    fn s(&mut self, target: IndexType) {
        self.gates.push(format!("s q[{}];", target));
    }

    fn v(&mut self, target: IndexType) {
        self.gates.push(format!("v q[{}];", target));
    }

    fn s_dgr(&mut self, target: IndexType) {
        self.gates.push(format!("sdg q[{}];", target));
    }

    fn v_dgr(&mut self, target: IndexType) {
        self.gates.push(format!("vdg q[{}];", target));
    }

    fn x(&mut self, target: IndexType) {
        self.gates.push(format!("x q[{}];", target));
    }

    fn y(&mut self, target: IndexType) {
        self.gates.push(format!("y q[{}];", target));
    }

    fn z(&mut self, target: IndexType) {
        self.gates.push(format!("z q[{}];", target));
    }

    fn h(&mut self, target: IndexType) {
        self.gates.push(format!("h q[{}];", target));
    }

    fn cx(&mut self, control: IndexType, target: IndexType) {
        self.gates.push(format!("cx q[{}], q[{}];", control, target));
    }

    fn cz(&mut self, control: IndexType, target: IndexType) {
        self.gates.push(format!("cz q[{}], q[{}];", control, target));
    }
}