use core::num;
use std::iter::{self, zip};

use itertools::{iproduct, Itertools};

use crate::{
    data_structures::{CliffordTableau, PauliLetter, PauliString, PropagateClifford},
    synthesis_methods::{custom::Custom, naive::Naive},
};

use super::CliffordGates;

fn get_pauli(pauli_string: &PauliString, row: usize) -> PauliLetter {
    PauliLetter::new(pauli_string.x(row), pauli_string.z(row))
}

fn is_i(pauli_letter: PauliLetter) -> bool {
    return pauli_letter == PauliLetter::I;
}

fn is_not_i(pauli_letter: PauliLetter) -> bool {
    return pauli_letter != PauliLetter::I;
}

fn is_x(pauli_letter: PauliLetter) -> bool {
    return pauli_letter == PauliLetter::X;
}

fn is_not_x(pauli_letter: PauliLetter) -> bool {
    return pauli_letter != PauliLetter::X;
}

fn is_y(pauli_letter: PauliLetter) -> bool {
    return pauli_letter == PauliLetter::Y;
}

fn is_not_y(pauli_letter: PauliLetter) -> bool {
    return pauli_letter != PauliLetter::Y;
}

fn is_z(pauli_letter: PauliLetter) -> bool {
    return pauli_letter == PauliLetter::Z;
}

fn is_not_z(pauli_letter: PauliLetter) -> bool {
    return pauli_letter != PauliLetter::Z;
}

pub struct CliffordTableauSynthesizer {
    custom_pivots: Option<Vec<usize>>,
}

impl<G> Naive<&CliffordTableau, G> for CliffordTableauSynthesizer
where
    G: CliffordGates,
{
    fn run_naive(clifford_tableau: &CliffordTableau, repr: &mut G) {
        let mut ct = clifford_tableau.adjoint();

        let num_qubits = ct.size();
        for row in 0..num_qubits {
            let pivot_col = naive_pivot_search(&ct, num_qubits, row);

            if pivot_col != row {
                swap(repr, &mut ct, row, pivot_col);
            }

            // Cleanup pivot column
            clean_pivot(repr, &mut ct, row, row);

            let checked_rows = (row + 1..num_qubits).collect::<Vec<_>>();

            // Use the pivot to remove all other terms in the X observable.
            clean_x_observables(repr, &mut ct, &checked_rows, row, row);

            // Use the pivot to remove all other terms in the Z observable.
            clean_z_observables(repr, &mut ct, &checked_rows, row, row);
        }

        clean_signs(repr, &mut ct, &(0..num_qubits).collect::<Vec<_>>());
    }
}

impl<G> Custom<&CliffordTableau, G> for CliffordTableauSynthesizer
where
    G: CliffordGates,
{
    fn run_custom(
        clifford_tableau: &CliffordTableau,
        repr: &mut G,
        custom_columns: Vec<usize>,
        custom_rows: Vec<usize>,
    ) {
        let mut ct = clifford_tableau.adjoint();
        let num_qubits = ct.size();
        let mut remaining_columns = (0..num_qubits).collect::<Vec<_>>();
        let mut remaining_rows = (0..num_qubits).collect::<Vec<_>>();

        assert_eq!(
            custom_columns.iter().copied().sorted().collect::<Vec<_>>(),
            remaining_columns
        );

        assert_eq!(
            custom_rows.iter().copied().sorted().collect::<Vec<_>>(),
            remaining_rows
        );

        for (&pivot_column, &pivot_row) in zip(&custom_columns, &custom_rows) {
            // Cleanup pivot column

            remaining_columns.retain(|&x| x != pivot_column);
            remaining_rows.retain(|&x| x != pivot_row);
            {
                clean_x_pivot(repr, &mut ct, pivot_column, pivot_row);

                // Use the pivot to remove all other terms in the X observable.
                clean_x_observables(repr, &mut ct, &remaining_rows, pivot_column, pivot_row);

                clean_z_pivot(repr, &mut ct, pivot_column, pivot_row);

                // Use the pivot to remove all other terms in the Z observable.
                clean_z_observables(repr, &mut ct, &remaining_rows, pivot_column, pivot_row);
            }
        }
        let final_permutation = zip(custom_columns.clone(), custom_rows.clone())
            .sorted_by_key(|a| a.1)
            .map(|a| a.0)
            .collect::<Vec<_>>();

        clean_signs(repr, &mut ct, &final_permutation);
    }
}

fn clean_pivot<G>(repr: &mut G, ct: &mut CliffordTableau, pivot_column: usize, row: usize)
where
    G: CliffordGates,
{
    let num_qubits = ct.size();
    if check_pauli(&*ct, pivot_column, row + num_qubits, is_y) {
        ct.s(row);
        repr.s(row);
    }

    if check_pauli(&*ct, pivot_column, row + num_qubits, is_not_z) {
        ct.h(row);
        repr.h(row);
    }

    if check_pauli(&*ct, pivot_column, row, is_not_x) {
        ct.s(row);
        repr.s(row);
    }
}

fn clean_x_pivot<G>(repr: &mut G, ct: &mut CliffordTableau, pivot_column: usize, pivot_row: usize)
where
    G: CliffordGates,
{
    if check_pauli(&*ct, pivot_row, pivot_column, is_y) {
        ct.s(pivot_row);
        repr.s(pivot_row);
    }

    if check_pauli(&*ct, pivot_row, pivot_column, is_z) {
        ct.h(pivot_row);
        repr.h(pivot_row);
    }
}

fn clean_z_pivot<G>(repr: &mut G, ct: &mut CliffordTableau, pivot_column: usize, pivot_row: usize)
where
    G: CliffordGates,
{
    let num_qubits = ct.size();
    if check_pauli(&*ct, pivot_row, pivot_column + num_qubits, is_y) {
        ct.v(pivot_row);
        repr.v(pivot_row);
    }

    if check_pauli(&*ct, pivot_row, pivot_column + num_qubits, is_x) {
        ct.h(pivot_row);
        repr.h(pivot_row);
    }
}

fn clean_x_observables<G>(
    repr: &mut G,
    ct: &mut CliffordTableau,
    remaining_rows: &[usize],
    pivot_col: usize,
    pivot_row: usize,
) where
    G: CliffordGates,
{
    let affected_cols = check_across_columns(&*ct, remaining_rows, pivot_col, is_y);

    for col in affected_cols {
        repr.s(col);
        ct.s(col);
    }

    let affected_cols = check_across_columns(&*ct, remaining_rows, pivot_col, is_z);

    for col in affected_cols {
        repr.h(col);
        ct.h(col);
    }

    let affected_cols = check_across_columns(&*ct, remaining_rows, pivot_col, is_not_i);

    for col in affected_cols {
        repr.cx(pivot_row, col);
        ct.cx(pivot_row, col);
    }
}

fn clean_z_observables<G>(
    repr: &mut G,
    ct: &mut CliffordTableau,
    remaining_rows: &[usize],
    pivot_col: usize,
    pivot_row: usize,
) where
    G: CliffordGates,
{
    let num_qubits = ct.size();
    let affected_cols = check_across_columns(&*ct, remaining_rows, pivot_col + num_qubits, is_y);
    for col in affected_cols {
        repr.v(col);
        ct.v(col);
    }

    let affected_cols = check_across_columns(&*ct, remaining_rows, pivot_col + num_qubits, is_x);
    for col in affected_cols {
        repr.h(col);
        ct.h(col);
    }

    let affected_cols =
        check_across_columns(&*ct, remaining_rows, pivot_col + num_qubits, is_not_i);
    for col in affected_cols {
        repr.cx(col, pivot_row);
        ct.cx(col, pivot_row);
    }
}

fn clean_signs<G>(repr: &mut G, ct: &mut CliffordTableau, row_permutation: &[usize])
where
    G: CliffordGates,
{
    let z_signs = ct.z_signs();

    for (sign, row) in zip(z_signs, row_permutation) {
        if sign {
            repr.x(*row);
            ct.x(*row);
        }
    }

    let x_signs = ct.x_signs();

    for (sign, row) in zip(x_signs, row_permutation) {
        if sign {
            repr.z(*row);
            ct.z(*row);
        }
    }
}

fn swap<G>(repr: &mut G, clifford_tableau: &mut CliffordTableau, row: usize, pivot_col: usize)
where
    G: CliffordGates,
{
    repr.cx(pivot_col, row);
    repr.cx(row, pivot_col);
    repr.cx(pivot_col, row);

    clifford_tableau.cx(pivot_col, row);
    clifford_tableau.cx(row, pivot_col);
    clifford_tableau.cx(pivot_col, row);
}

fn naive_pivot_search(clifford_tableau: &CliffordTableau, num_qubits: usize, row: usize) -> usize {
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

fn check_pauli(
    clifford_tableau: &CliffordTableau,
    column: usize,
    row: usize,
    pauli_check: fn(PauliLetter) -> bool,
) -> bool {
    let pauli_string = clifford_tableau.column(column);
    pauli_check(get_pauli(pauli_string, row))
}

/// Helper function that returns indices for a particular column `col` of `clifford_tableau` that match the provided closure `pauli_check`.
fn check_across_columns(
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
