use std::collections::VecDeque;

use bitvec::bitvec;
use bitvec::prelude::Lsb0;
use syn::data_structures::{CliffordTableau, PauliPolynomial, PauliString};
use syn::ir::pauli_exponential::PauliExponential;

fn main() {
    // test tableaus
    // Stab: ZZZ, -YIY, XIX
    // Destab: -IXI, XXI, IYY
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
    let my_tableaus = CliffordTableau::from_parts(vec![pauli_1, pauli_2, pauli_3], signs);

    println!("{}", my_tableaus);
    //test pauli polynomial

    let ham = vec![("IXYZ", 0.3), ("XXII", 0.7), ("YYII", 0.12)];
    let pp = PauliPolynomial::from_hamiltonian(ham);
    println!("{}", pp);

    // // visualize_pauli_exponential_simple(&pe);
    let ham = vec![("IZZZ", 0.3)];
    let pp = PauliPolynomial::from_hamiltonian(ham);
    let ct = CliffordTableau::new(4);
    let pe = PauliExponential::new(VecDeque::from([pp]), ct);
    println!("{}", pe);

    // //visualize_pauli_exponential complex
    let ham = vec![("IXYZ", 0.3), ("XXII", 0.7), ("YYII", 0.12)];

    let pauli_polynomial = PauliPolynomial::from_hamiltonian(ham);
    let clifford_tableau = CliffordTableau::new(4);
    let complex_pe = PauliExponential::new(VecDeque::from([pauli_polynomial]), clifford_tableau);
    println!("{}", complex_pe);
}
