use std::collections::VecDeque;

use synir::{
    data_structures::{Angle, PauliPolynomial, PauliString, PropagateClifford},
    ir::pauli_polynomial,
};

pub fn setup_complex_pp() -> VecDeque<PauliPolynomial> {
    let ham_1 = vec![("IZZZ", Angle::from_angle(0.3))];
    let ham_2 = vec![("XXII", Angle::from_angle(0.7))];

    let pp_1 = PauliPolynomial::from_hamiltonian(ham_1);
    let pp_2 = PauliPolynomial::from_hamiltonian(ham_2);
    VecDeque::from([pp_1, pp_2])
}

pub fn setup_simple_pp() -> VecDeque<PauliPolynomial> {
    let ham = vec![("IXYZ", Angle::from_angle(0.3))];

    let pauli_polynomial = PauliPolynomial::from_hamiltonian(ham);

    VecDeque::from([pauli_polynomial])
}

pub fn setup_parallel_pp() -> VecDeque<PauliPolynomial> {
    let mut pauli_polynomial = PauliPolynomial::new(2);
    pauli_polynomial.append_gadget(PauliString::from_text("IZ"), Angle::from_angle(0.53423));
    pauli_polynomial.append_gadget(PauliString::from_text("XI"), Angle::from_angle(0.234234));
    pauli_polynomial.cx(0, 1);
    pauli_polynomial.v(0);
    VecDeque::from([pauli_polynomial])
}
