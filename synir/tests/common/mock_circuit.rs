use synir::{
    architecture::connectivity::{self, Connectivity},
    data_structures::{CliffordTableau, PropagateClifford},
    ir::{CliffordGates, Gates},
    IndexType,
};

use super::unitary::{self, Unitary};

#[derive(Debug, Default, PartialEq)]
pub struct MockCircuit {
    commands: Vec<MockCommand>,
    /// Mirrors `QiskitSynIR`, which *stores* the final permutation instead of expanding it into
    /// CX-swaps the way the default `CliffordGates::add_final_permutation` does. Keeping the two
    /// mocks aligned is the whole point: with the default expansion the Rust tests would exercise
    /// a code path that the Python wrapper never takes.
    final_permutation: Option<Vec<IndexType>>,
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
            final_permutation: None,
        }
    }

    pub fn from_vec(commands: Vec<MockCommand>) -> Self {
        Self {
            commands,
            final_permutation: None,
        }
    }

    pub fn commands(&self) -> &Vec<MockCommand> {
        &self.commands
    }

    pub fn final_permutation(&self) -> Option<&Vec<IndexType>> {
        self.final_permutation.as_ref()
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
                MockCommand::Y(target) => {
                    tableau.y(*target);
                }
                MockCommand::CZ(control, target) => {
                    tableau.cz(*control, *target);
                }
                MockCommand::Rx(_, _) | MockCommand::Ry(_, _) | MockCommand::Rz(_, _) => {
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

    /// Dense unitary of the emitted gates, in time order (first command executed first), ignoring
    /// any final permutation.
    pub fn to_unitary(&self, size: usize) -> Unitary {
        let mut u = Unitary::identity(size);
        for command in self.commands.iter() {
            match command {
                MockCommand::H(t) => u.apply_1q(*t, unitary::h_matrix()),
                MockCommand::S(t) => u.apply_1q(*t, unitary::s_matrix()),
                MockCommand::SDgr(t) => u.apply_1q(*t, unitary::sdg_matrix()),
                MockCommand::V(t) => u.apply_1q(*t, unitary::sx_matrix()),
                MockCommand::VDgr(t) => u.apply_1q(*t, unitary::sxdg_matrix()),
                MockCommand::X(t) => u.apply_1q(*t, unitary::x_matrix()),
                MockCommand::Y(t) => u.apply_1q(*t, unitary::y_matrix()),
                MockCommand::Z(t) => u.apply_1q(*t, unitary::z_matrix()),
                MockCommand::Rx(t, angle) => u.apply_1q(*t, unitary::rx_matrix(*angle)),
                MockCommand::Ry(t, angle) => u.apply_1q(*t, unitary::ry_matrix(*angle)),
                MockCommand::Rz(t, angle) => u.apply_1q(*t, unitary::rz_matrix(*angle)),
                MockCommand::CX(c, t) => u.apply_cx(*c, *t),
                MockCommand::CZ(c, t) => u.apply_cz(*c, *t),
            }
        }
        u
    }

    /// Dense unitary of the emitted gates followed by a final qubit permutation.
    ///
    /// The permutation is applied exactly the way `synpy.qiskit.plugin.circuit_to_circuit` applies
    /// it: a Qiskit `PermutationGate` built from the *inverse* of the reported pattern, appended
    /// after the synthesized circuit. See [`permutation_variants`] for the alternatives, which the
    /// convention tests compare against.
    pub fn to_unitary_with_permutation(
        &self,
        size: usize,
        permutation: Option<Vec<IndexType>>,
    ) -> Unitary {
        let gates = self.to_unitary(size);
        match permutation {
            None => gates,
            Some(pattern) => {
                let inverse = Unitary::invert_pattern(&pattern);
                Unitary::permutation(&inverse, size).mul(&gates)
            }
        }
    }

    /// The same circuit under each of the four plausible readings of a reported permutation.
    ///
    /// Returned as `(label, unitary)` pairs so a failing assertion can say *which* reading would
    /// have worked instead of only that the chosen one did not.
    pub fn permutation_variants(
        &self,
        size: usize,
        pattern: &[IndexType],
    ) -> Vec<(&'static str, Unitary)> {
        let gates = self.to_unitary(size);
        let inverse = Unitary::invert_pattern(pattern);
        let forward = Unitary::permutation(pattern, size);
        let backward = Unitary::permutation(&inverse, size);
        vec![
            ("inverse pattern, appended after", backward.mul(&gates)),
            ("forward pattern, appended after", forward.mul(&gates)),
            ("inverse pattern, prepended before", gates.mul(&backward)),
            ("forward pattern, prepended before", gates.mul(&forward)),
        ]
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

    /// Stores the permutation rather than expanding it into CX-swaps, mirroring
    /// `QiskitSynIR::add_final_permutation` in `synpy/src/wrapper/qiskit.rs`.
    fn add_final_permutation(&mut self, permutation: Vec<IndexType>) {
        self.final_permutation = Some(permutation);
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
