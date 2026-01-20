use std::collections::VecDeque;

use synir::data_structures::PauliPolynomial;

pub fn setup_complex_pp() -> VecDeque<PauliPolynomial> {
    let ham_1 = vec![("IZZZ", 0.3)];
    let ham_2 = vec![("XXII", 0.7)];

    let pp_1 = PauliPolynomial::from_hamiltonian(ham_1);
    let pp_2 = PauliPolynomial::from_hamiltonian(ham_2);
    VecDeque::from([pp_1, pp_2])
}

pub fn setup_simple_pp() -> VecDeque<PauliPolynomial> {
    let ham = vec![("IXYZ", 0.3)];

    let pauli_polynomial = PauliPolynomial::from_hamiltonian(ham);

    VecDeque::from([pauli_polynomial])
}
