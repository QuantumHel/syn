use std::iter::zip;

use crate::{
    data_structures::{CliffordTableau, PauliLetter, PauliVec, PropagateClifford},
    ir::CliffordGates,
};

fn get_pauli(pauli_string: &PauliVec, row: usize) -> PauliLetter {
    PauliLetter::new(pauli_string.x(row), pauli_string.z(row))
}

#[allow(dead_code)]
fn is_i(pauli_letter: PauliLetter) -> bool {
    pauli_letter == PauliLetter::I
}

fn is_not_i(pauli_letter: PauliLetter) -> bool {
    pauli_letter != PauliLetter::I
}

fn is_x(pauli_letter: PauliLetter) -> bool {
    pauli_letter == PauliLetter::X
}

fn is_not_x(pauli_letter: PauliLetter) -> bool {
    pauli_letter != PauliLetter::X
}

fn is_y(pauli_letter: PauliLetter) -> bool {
    pauli_letter == PauliLetter::Y
}

#[allow(dead_code)]
fn is_not_y(pauli_letter: PauliLetter) -> bool {
    pauli_letter != PauliLetter::Y
}

fn is_z(pauli_letter: PauliLetter) -> bool {
    pauli_letter == PauliLetter::Z
}

fn is_not_z(pauli_letter: PauliLetter) -> bool {
    pauli_letter != PauliLetter::Z
}

pub(super) fn clean_pivot<G>(
    repr: &mut G,
    ct: &mut CliffordTableau,
    pivot_column: usize,
    pivot_row: usize,
) where
    G: CliffordGates,
{
    let num_qubits = ct.size();
    if check_pauli(&*ct, pivot_row, pivot_column + num_qubits, is_y) {
        ct.s(pivot_row);
        repr.s(pivot_row);
    }

    if check_pauli(&*ct, pivot_row, pivot_column + num_qubits, is_not_z) {
        ct.h(pivot_row);
        repr.h(pivot_row);
    }

    if check_pauli(&*ct, pivot_row, pivot_column, is_not_x) {
        ct.s(pivot_row);
        repr.s(pivot_row);
    }
}

pub(super) fn clean_x_pivot<G>(
    repr: &mut G,
    clifford_tableau: &mut CliffordTableau,
    pivot_column: usize,
    pivot_row: usize,
) where
    G: CliffordGates,
{
    // These are switched around because of implementation
    if check_pauli(&*clifford_tableau, pivot_row, pivot_column, is_y) {
        clifford_tableau.s(pivot_row);
        repr.s(pivot_row);
    }

    // These are switched around because of implementation
    if check_pauli(&*clifford_tableau, pivot_row, pivot_column, is_z) {
        clifford_tableau.h(pivot_row);
        repr.h(pivot_row);
    }
}

pub(super) fn clean_z_pivot<G>(
    repr: &mut G,
    clifford_tableau: &mut CliffordTableau,
    pivot_column: usize,
    pivot_row: usize,
) where
    G: CliffordGates,
{
    let num_qubits = clifford_tableau.size();
    // These are switched around because of implementation
    if check_pauli(
        &*clifford_tableau,
        pivot_row,
        pivot_column + num_qubits,
        is_y,
    ) {
        clifford_tableau.v(pivot_row);
        repr.v(pivot_row);
    }

    // These are switched around because of implementation
    if check_pauli(
        &*clifford_tableau,
        pivot_row,
        pivot_column + num_qubits,
        is_x,
    ) {
        clifford_tableau.h(pivot_row);
        repr.h(pivot_row);
    }
}

pub(super) fn clean_x_observables<G>(
    repr: &mut G,
    clifford_tableau: &mut CliffordTableau,
    remaining_rows: &[usize],
    pivot_column: usize,
    pivot_row: usize,
) where
    G: CliffordGates,
{
    let affected_cols =
        check_across_columns(&*clifford_tableau, remaining_rows, pivot_column, is_y);

    for col in affected_cols {
        repr.s(col);
        clifford_tableau.s(col);
    }

    let affected_cols =
        check_across_columns(&*clifford_tableau, remaining_rows, pivot_column, is_z);

    for col in affected_cols {
        repr.h(col);
        clifford_tableau.h(col);
    }

    let affected_cols =
        check_across_columns(&*clifford_tableau, remaining_rows, pivot_column, is_not_i);

    for col in affected_cols {
        repr.cx(pivot_row, col);
        clifford_tableau.cx(pivot_row, col);
    }
}

pub(super) fn clean_z_observables<G>(
    repr: &mut G,
    clifford_tableau: &mut CliffordTableau,
    remaining_rows: &[usize],
    pivot_column: usize,
    pivot_row: usize,
) where
    G: CliffordGates,
{
    let num_qubits = clifford_tableau.size();
    let affected_cols = check_across_columns(
        &*clifford_tableau,
        remaining_rows,
        pivot_column + num_qubits,
        is_y,
    );
    for col in affected_cols {
        repr.v(col);
        clifford_tableau.v(col);
    }

    let affected_cols = check_across_columns(
        &*clifford_tableau,
        remaining_rows,
        pivot_column + num_qubits,
        is_x,
    );
    for col in affected_cols {
        repr.h(col);
        clifford_tableau.h(col);
    }

    let affected_cols = check_across_columns(
        &*clifford_tableau,
        remaining_rows,
        pivot_column + num_qubits,
        is_not_i,
    );
    for col in affected_cols {
        repr.cx(col, pivot_row);
        clifford_tableau.cx(col, pivot_row);
    }
}

pub(super) fn clean_signs<G>(
    repr: &mut G,
    clifford_tableau: &mut CliffordTableau,
    row_permutation: &[usize],
) where
    G: CliffordGates,
{
    let z_signs = clifford_tableau.z_signs();

    for (sign, row) in zip(z_signs, row_permutation) {
        if sign {
            repr.x(*row);
            clifford_tableau.x(*row);
        }
    }

    let x_signs = clifford_tableau.x_signs();

    for (sign, row) in zip(x_signs, row_permutation) {
        if sign {
            repr.z(*row);
            clifford_tableau.z(*row);
        }
    }
}

pub(super) fn swap<G>(
    repr: &mut G,
    clifford_tableau: &mut CliffordTableau,
    row: usize,
    pivot_col: usize,
) where
    G: CliffordGates,
{
    repr.cx(pivot_col, row);
    repr.cx(row, pivot_col);
    repr.cx(pivot_col, row);

    clifford_tableau.cx(pivot_col, row);
    clifford_tableau.cx(row, pivot_col);
    clifford_tableau.cx(pivot_col, row);
}

pub(super) fn naive_pivot_search(
    clifford_tableau: &CliffordTableau,
    num_qubits: usize,
    row: usize,
) -> usize {
    let mut pivot_col = 0;

    for col in 0..num_qubits {
        let column = clifford_tableau.column(col);
        let x_pauli = get_pauli(column, row);
        let z_pauli = get_pauli(column, row + num_qubits);
        if x_pauli != PauliLetter::I && z_pauli != PauliLetter::I && x_pauli != z_pauli {
            pivot_col = col;
            break;
        }
    }

    pivot_col
}

pub(super) fn check_pauli(
    clifford_tableau: &CliffordTableau,
    column: usize,
    row: usize,
    pauli_check: fn(PauliLetter) -> bool,
) -> bool {
    let pauli_string = clifford_tableau.column(column);
    pauli_check(get_pauli(pauli_string, row))
}

/// Helper function that returns indices for a particular column `col` of `clifford_tableau` that match the provided closure `pauli_check`.
pub(super) fn check_across_columns(
    clifford_tableau: &CliffordTableau,
    columns: &[usize],
    checked_row: usize,
    pauli_check: fn(PauliLetter) -> bool,
) -> Vec<usize> {
    let mut affected_cols = Vec::new();

    for column in columns {
        let pauli_string = clifford_tableau.column(*column);

        if pauli_check(get_pauli(pauli_string, checked_row)) {
            affected_cols.push(*column);
        }
    }
    affected_cols
}
