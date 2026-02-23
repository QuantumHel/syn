use synir::{
    architecture::connectivity::{self, Connectivity},
    data_structures::{CliffordTableau, PropagateClifford},
    ir::{CliffordGates, Gates},
    IndexType,
};

#[derive(Debug, Default, PartialEq)]
pub struct MockCircuit {
    commands: Vec<MockCommand>,
}

#[derive(Debug, PartialEq, Clone)]
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

    pub fn from_vec(commands: Vec<MockCommand>) -> Self {
        Self { commands }
    }

    pub fn commands(&self) -> &Vec<MockCommand> {
        &self.commands
    }

    pub fn to_clifford_tableau(&self, size: usize) -> CliffordTableau {
        let mut tableau = CliffordTableau::new(size);
        for command in self.commands.iter() {
            match command {
                MockCommand::H(target) => {
                    tableau.h(*target);
                }
                MockCommand::S(target) => {
                    tableau.s(*target);
                }
                MockCommand::SDgr(target) => {
                    tableau.s_dgr(*target);
                }
                MockCommand::V(target) => {
                    tableau.v(*target);
                }
                MockCommand::VDgr(target) => {
                    tableau.v_dgr(*target);
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
                    panic!("Circuit contains non-cliffords")
                }
            }
        }
        tableau
    }

    pub fn equals_clifford_tableau(
        &self,
        clifford_tableau: &CliffordTableau,
        permutation: Option<Vec<usize>>,
    ) -> bool {
        let mut ref_ct = self.to_clifford_tableau(clifford_tableau.size());
        if permutation.is_some() {
            ref_ct.permute(permutation.unwrap());
        }
        *clifford_tableau == ref_ct
    }

    pub fn fits_connectivity(&self, connectivity: &Connectivity) -> bool {
        let mut result = true;
        for command in self.commands.iter() {
            result &= match command {
                MockCommand::CX(i, j) => connectivity.has_edge(*i, *j),
                MockCommand::CZ(i, j) => connectivity.has_edge(*i, *j),
                _ => true,
            };
        }
        result
    }

    pub fn cliffords_only(&self) -> MockCircuit {
        MockCircuit::from_vec(
            self.commands
                .iter()
                .filter(|c| match **c {
                    MockCommand::Rx(_, _) => false,
                    MockCommand::Ry(_, _) => false,
                    MockCommand::Rz(_, _) => false,
                    _ => true,
                })
                .map(|c| (*c).clone())
                .collect(),
        )
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
    fn rx(&mut self, target: IndexType, angle: f64) {
        self.commands.push(MockCommand::Rx(target, angle));
    }

    fn ry(&mut self, target: IndexType, angle: f64) {
        self.commands.push(MockCommand::Ry(target, angle));
    }

    fn rz(&mut self, target: IndexType, angle: f64) {
        self.commands.push(MockCommand::Rz(target, angle));
    }
}
