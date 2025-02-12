use syn::ir::CliffordGates;

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
    fn s(&mut self, target: syn::IndexType) {
        self.commands.push(MockCommand::S(target));
    }

    fn v(&mut self, target: syn::IndexType) {
        self.commands.push(MockCommand::V(target));
    }

    fn s_dgr(&mut self, target: syn::IndexType) {
        self.commands.push(MockCommand::SDgr(target));
    }

    fn v_dgr(&mut self, target: syn::IndexType) {
        self.commands.push(MockCommand::VDgr(target));
    }

    fn x(&mut self, target: syn::IndexType) {
        self.commands.push(MockCommand::X(target));
    }

    fn y(&mut self, target: syn::IndexType) {
        self.commands.push(MockCommand::Y(target));
    }

    fn z(&mut self, target: syn::IndexType) {
        self.commands.push(MockCommand::Z(target));
    }

    fn h(&mut self, target: syn::IndexType) {
        self.commands.push(MockCommand::H(target));
    }

    fn cx(&mut self, control: syn::IndexType, target: syn::IndexType) {
        self.commands.push(MockCommand::CX(control, target));
    }

    fn cz(&mut self, control: syn::IndexType, target: syn::IndexType) {
        self.commands.push(MockCommand::CZ(control, target));
    }
}
