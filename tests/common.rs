use synir::{
    data_structures::{CliffordTableau, PropagateClifford},
    ir::{CliffordGates, Gates},
    IndexType,
};

type Angle = f64;
#[derive(Debug, Default)]
pub struct MockCircuit {
    commands: Vec<MockCommand>,
}

#[derive(Debug, PartialEq)]
pub enum MockCommand {
    CX(usize, usize),
    CZ(usize, usize),
    X(usize),
    Y(usize),
    Z(usize),
    H(usize),
    S(usize),
    V(usize),
    SDgr(usize),
    VDgr(usize),
    Rx(usize, f64),
    Ry(usize, f64),
    Rz(usize, f64),
}

impl MockCircuit {
    pub fn new() -> Self {
        Self {
            commands: Vec::new(),
        }
    }
    pub fn commands(&self) -> &Vec<MockCommand> {
        &self.commands
    }
}

impl CliffordGates for MockCircuit {
    fn s(&mut self, target: IndexType) {
        self.commands.push(MockCommand::S(target));
    }

    fn v(&mut self, target: IndexType) {
        self.commands.push(MockCommand::V(target));
    }

    fn s_dgr(&mut self, target: IndexType) {
        self.commands.push(MockCommand::SDgr(target));
    }

    fn v_dgr(&mut self, target: IndexType) {
        self.commands.push(MockCommand::VDgr(target));
    }

    fn x(&mut self, target: IndexType) {
        self.commands.push(MockCommand::X(target));
    }

    fn y(&mut self, target: IndexType) {
        self.commands.push(MockCommand::Y(target));
    }

    fn z(&mut self, target: IndexType) {
        self.commands.push(MockCommand::Z(target));
    }

    fn h(&mut self, target: IndexType) {
        self.commands.push(MockCommand::H(target));
    }

    fn cx(&mut self, control: IndexType, target: IndexType) {
        self.commands.push(MockCommand::CX(control, target));
    }

    fn cz(&mut self, control: IndexType, target: IndexType) {
        self.commands.push(MockCommand::CZ(control, target));
    }
}

impl Gates for MockCircuit {
    fn rx(&mut self, target: IndexType, angle: Angle) {
        self.commands.push(MockCommand::Rx(target, angle));
    }

    fn ry(&mut self, target: IndexType, angle: Angle) {
        self.commands.push(MockCommand::Ry(target, angle));
    }

    fn rz(&mut self, target: IndexType, angle: Angle) {
        self.commands.push(MockCommand::Rz(target, angle));
    }
}

pub fn parse_clifford_commands(size: usize, commands: &[MockCommand]) -> CliffordTableau {
    let mut tableau = CliffordTableau::new(size);
    for command in commands.iter() {
        match command {
            MockCommand::H(target) => {
                tableau.h(*target);
            }
            MockCommand::S(target) => {
                tableau.s(*target);
            }
            MockCommand::V(target) => {
                tableau.v(*target);
            }
            MockCommand::CX(control, target) => {
                tableau.cx(*control, *target);
            }
            MockCommand::X(target) => {
                tableau.x(*target);
            }
            MockCommand::Z(target) => {
                tableau.z(*target);
            }
            _ => {
                panic!("not found")
            }
        }
    }
    tableau
}
