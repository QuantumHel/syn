use crate::connectivity::connectivity_benchmark;
use bitvec::bitvec;
use bitvec::prelude::Lsb0;
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use syn::data_structures::CliffordTableau;
use syn::data_structures::PauliString;
use syn::ir::clifford_tableau::CliffordTableauSynthesizer;
use syn::ir::CliffordGates;
use syn::synthesis_methods::custom::Custom;
use syn::synthesis_methods::naive::Naive;

mod connectivity;

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

fn setup_sample_ct() -> CliffordTableau {
    // Stab: ZZZ, -YIY, XIX
    // Destab: -IXI, XXI, IYY
    let ct_size = 3;
    // qubit 1x: ZYI
    // qubit 1z: IZZ
    let pauli_1 = PauliString::from_text("ZYIIZZ");

    // qubit 2x: ZIX
    // qubit 2z: XII
    let pauli_2 = PauliString::from_text("ZIXXII");

    // qubit 3x: ZYY
    // qubit 3z: IIZ
    let pauli_3 = PauliString::from_text("ZYYIIZ");

    let signs = bitvec![0, 1, 0, 1, 0, 0];
    CliffordTableau::from_parts(vec![pauli_1, pauli_2, pauli_3], signs, ct_size)
}

fn setup_sample_inverse_ct() -> CliffordTableau {
    // Stab: -ZIYZ, -ZZYZ, -XZXI, IZXX
    // Destab: -YYIZ, -YYXZ, ZIXX, -XZXZ
    let ct_size = 4;
    // qubit 1x: ZZXI
    // qubit 1z: YYZX
    let pauli_1 = PauliString::from_text("ZZXIYYZX");

    // qubit 2x: IZZZ
    // qubit 2z: YYIZ
    let pauli_2 = PauliString::from_text("IZZZYYIZ");

    // qubit 3x: YYXX
    // qubit 3z: IXXX
    let pauli_3 = PauliString::from_text("YYXXIXXX");

    // qubit 3x: ZZIX
    // qubit 3z: ZZXZ
    let pauli_4 = PauliString::from_text("ZZIXZZXZ");

    let signs = bitvec![1, 1, 1, 0, 1, 1, 0, 1];
    CliffordTableau::from_parts(vec![pauli_1, pauli_2, pauli_3, pauli_4], signs, ct_size)
}

fn naive_ct(clifford: CliffordTableau) {
    let mut mock = MockCircuit::new();
    CliffordTableauSynthesizer::run_naive(&clifford, &mut mock);
}

fn custom_ct(clifford: CliffordTableau) {
    let num_qubits = clifford.size();

    let mut mock = MockCircuit::new();
    CliffordTableauSynthesizer::run_custom(
        &clifford,
        &mut mock,
        (0..num_qubits).collect(),
        (0..num_qubits).collect(),
    );
}

pub fn ct_bench(c: &mut Criterion) {
    let mut group = c.benchmark_group("Fibonacci");

    group.bench_function(BenchmarkId::new("Naive", "naive"), |b| {
        b.iter(|| naive_ct(black_box(setup_sample_inverse_ct())))
    });
    group.bench_function(BenchmarkId::new("Custom", "custom"), |b| {
        b.iter(|| custom_ct(black_box(setup_sample_inverse_ct())))
    });
}

criterion_group!(benches, ct_bench);
criterion_main!(benches, connectivity_benchmark);
