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
    // qubit 4x: ZYX
    // qubit 4z: IZI
    let signs = bitvec![0, 1, 0, 1, 0, 0, 1, 1];
    let my_tableaux = CliffordTableau::from_parts(vec![pauli_1, pauli_2, pauli_3], signs);
    println!("Test clifford tableaux small");
    println!("{}", my_tableaux);
    let big_tableaux = CliffordTableau::new(10);
    println!("Test clifford tableaux big");
    println!("{}", big_tableaux);

    let ham = vec![("IXYZ", 0.3), ("XXII", 0.7), ("YYII", 0.12)];
    //pauli chain IXY; XXY; YII; ZII
    let pp = PauliPolynomial::from_hamiltonian(ham);
    println!("Test pauli polynomial");
    println!("{}", pp);

    let ham1 = vec![("IZZZ", 0.3)];
    let pp1 = PauliPolynomial::from_hamiltonian(ham1);
    let ct = CliffordTableau::new(4);
    let ham2 = vec![("XIII", 0.7)];
    let pp2 = PauliPolynomial::from_hamiltonian(ham2);
    let pe = PauliExponential::new(VecDeque::from([pp, pp1, pp2]), ct);

    println!("Test pauli exponential");
    println!("{}", pe);
}
