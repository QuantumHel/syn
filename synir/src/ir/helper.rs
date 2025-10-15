use crate::data_structures::PauliLetter;

#[allow(dead_code)]
pub(super) fn is_i(pauli_letter: PauliLetter) -> bool {
    pauli_letter == PauliLetter::I
}

pub(super) fn is_not_i(pauli_letter: PauliLetter) -> bool {
    pauli_letter != PauliLetter::I
}

pub(super) fn is_x(pauli_letter: PauliLetter) -> bool {
    pauli_letter == PauliLetter::X
}

pub(super) fn is_not_x(pauli_letter: PauliLetter) -> bool {
    pauli_letter != PauliLetter::X
}

pub(super) fn is_y(pauli_letter: PauliLetter) -> bool {
    pauli_letter == PauliLetter::Y
}

#[allow(dead_code)]
pub(super) fn is_not_y(pauli_letter: PauliLetter) -> bool {
    pauli_letter != PauliLetter::Y
}

pub(super) fn is_z(pauli_letter: PauliLetter) -> bool {
    pauli_letter == PauliLetter::Z
}

pub(super) fn is_not_z(pauli_letter: PauliLetter) -> bool {
    pauli_letter != PauliLetter::Z
}
