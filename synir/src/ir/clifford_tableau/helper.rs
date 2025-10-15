use std::iter::zip;

use itertools::Itertools;

use crate::{
    architecture::{connectivity::Connectivity, Architecture},
    data_structures::{CliffordTableau, PauliLetter, PauliString, PropagateClifford},
    ir::CliffordGates,
};

fn get_pauli(pauli_string: &PauliString, row: usize) -> PauliLetter {
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

pub(super) fn clean_naive_pivot<G>(
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

pub(super) fn clean_pivot<G>(
    repr: &mut G,
    clifford_tableau: &mut CliffordTableau,
    pivot_column: usize,
    pivot_row: usize,
    letter: PauliLetter,
) where
    G: CliffordGates,
{
    match letter {
        PauliLetter::X => clean_x_pivot(repr, clifford_tableau, pivot_column, pivot_row),
        PauliLetter::Z => clean_z_pivot(repr, clifford_tableau, pivot_column, pivot_row),
        _ => panic!("Invalid Pauli letter for pivot cleaning"),
    }
}

pub(super) fn clean_observables<G>(
    repr: &mut G,
    clifford_tableau: &mut CliffordTableau,
    remaining_rows: &[usize],
    pivot_column: usize,
    pivot_row: usize,
    letter: PauliLetter,
) where
    G: CliffordGates,
{
    match letter {
        PauliLetter::X => clean_x_observables(
            repr,
            clifford_tableau,
            remaining_rows,
            pivot_column,
            pivot_row,
        ),
        PauliLetter::Z => clean_z_observables(
            repr,
            clifford_tableau,
            remaining_rows,
            pivot_column,
            pivot_row,
        ),
        _ => panic!("Invalid Pauli letter for observable cleaning"),
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
    if check_pauli(&*clifford_tableau, pivot_column, pivot_row, is_y) {
        clifford_tableau.s(pivot_column);
        repr.s(pivot_column);
    }

    // These are switched around because of implementation
    if check_pauli(&*clifford_tableau, pivot_column, pivot_row, is_z) {
        clifford_tableau.h(pivot_column);
        repr.h(pivot_column);
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
        pivot_column,
        pivot_row + num_qubits,
        is_y,
    ) {
        clifford_tableau.v(pivot_column);
        repr.v(pivot_column);
    }

    // These are switched around because of implementation
    if check_pauli(
        &*clifford_tableau,
        pivot_column,
        pivot_row + num_qubits,
        is_x,
    ) {
        clifford_tableau.h(pivot_column);
        repr.h(pivot_column);
    }
}

pub(super) fn clean_x_observables<G>(
    repr: &mut G,
    clifford_tableau: &mut CliffordTableau,
    remaining_columns: &[usize],
    pivot_column: usize,
    pivot_row: usize,
) where
    G: CliffordGates,
{
    let affected_cols =
        check_across_columns(&*clifford_tableau, remaining_columns, pivot_row, is_y);

    for col in affected_cols {
        repr.s(col);
        clifford_tableau.s(col);
    }

    let affected_cols =
        check_across_columns(&*clifford_tableau, remaining_columns, pivot_row, is_z);

    for col in affected_cols {
        repr.h(col);
        clifford_tableau.h(col);
    }

    let affected_cols =
        check_across_columns(&*clifford_tableau, remaining_columns, pivot_row, is_not_i);

    for col in affected_cols {
        repr.cx(pivot_column, col);
        clifford_tableau.cx(pivot_column, col);
    }
}

pub(super) fn clean_z_observables<G>(
    repr: &mut G,
    clifford_tableau: &mut CliffordTableau,
    remaining_columns: &[usize],
    pivot_column: usize,
    pivot_row: usize,
) where
    G: CliffordGates,
{
    let num_qubits = clifford_tableau.size();
    let affected_cols = check_across_columns(
        &*clifford_tableau,
        remaining_columns,
        pivot_row + num_qubits,
        is_y,
    );
    for col in affected_cols {
        repr.v(col);
        clifford_tableau.v(col);
    }

    let affected_cols = check_across_columns(
        &*clifford_tableau,
        remaining_columns,
        pivot_row + num_qubits,
        is_x,
    );
    for col in affected_cols {
        repr.h(col);
        clifford_tableau.h(col);
    }

    let affected_cols = check_across_columns(
        &*clifford_tableau,
        remaining_columns,
        pivot_row + num_qubits,
        is_not_i,
    );
    for col in affected_cols {
        repr.cx(col, pivot_column);
        clifford_tableau.cx(col, pivot_column);
    }
}

pub(super) fn clean_signs<G>(repr: &mut G, clifford_tableau: &mut CliffordTableau)
where
    G: CliffordGates,
{
    let z_signs = clifford_tableau.z_signs();
    let inv_perm = match clifford_tableau.get_permutation() {
        None => panic!("Cleaning signs but tableau is not a permutation matrix: \n{}", clifford_tableau),
        Some(perm) => perm
    };
    let row_permutation = (0..clifford_tableau.size())
        .into_iter()
        .map(|i| inv_perm.iter().find_position(|&&x| x == i))
        .map(|x| x.unwrap().0)
        .collect_vec();
    for (sign, row) in zip(z_signs, row_permutation.iter()) {
        if sign {
            repr.x(*row);
            clifford_tableau.x(*row);
        }
    }

    let x_signs = clifford_tableau.x_signs();

    for (sign, row) in zip(x_signs, row_permutation.iter()) {
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

/// function to pick a stabilizer / destabilizer to set to identity in Clifford tableau.
pub(super) fn pick_row(
    clifford_tableau: &CliffordTableau,
    connectivity: &Connectivity,
    remaining_rows: &[usize],
) -> usize {
    let mut row_weights = vec![usize::MAX; clifford_tableau.size()];
    for row in remaining_rows {
        row_weights[*row] = 0;
        for qubit in connectivity.nodes() {
            if is_not_i(clifford_tableau.stabilizer(qubit, *row)) {
                row_weights[*row] += 1;
            }
            if is_not_i(clifford_tableau.destabilizer(qubit, *row)) {
                row_weights[*row] += 1;
            }
        }
    }

    row_weights
        .into_iter()
        .enumerate()
        .min_by_key(|&(_, weight)| weight)
        .map(|(index, _)| index)
        .unwrap()
}

/// function to pick a qubit to disconnect in Clifford tableau.
pub(super) fn pick_column(
    clifford_tableau: &CliffordTableau,
    connectivity: &Connectivity,
    pivot_row: usize,
) -> usize {
    let mut column_weights = vec![usize::MAX; clifford_tableau.size()];

    let non_cutting = connectivity.non_cutting();

    for qubit in non_cutting {
        column_weights[*qubit] = 0;
        for interaction in connectivity.nodes() {
            if interaction != pivot_row {
                let mult_z =
                    (clifford_tableau.stabilizer(*qubit, interaction) != PauliLetter::I) as usize;
                let mult_x =
                    (clifford_tableau.destabilizer(*qubit, interaction) != PauliLetter::I) as usize;
                column_weights[*qubit] +=
                    connectivity.distance(*qubit, interaction) * (mult_x + mult_z);
            }
        }
    }
    column_weights
        .iter()
        .enumerate()
        .min_by_key(|&(_, &weight)| weight)
        .map(|(index, _)| index)
        .unwrap()
}

pub(super) fn clean_prc<G>(
    repr: &mut G,
    clifford_tableau: &mut CliffordTableau,
    connectivity: &Connectivity,
    remaining_rows: &[usize],
    pivot_column: usize,
    pivot_row: usize,
    letter: PauliLetter,
) where
    G: CliffordGates,
{
    match letter {
        PauliLetter::X => clean_x_prc(
            repr,
            clifford_tableau,
            connectivity,
            remaining_rows,
            pivot_column,
            pivot_row,
        ),
        PauliLetter::Z => clean_z_prc(
            repr,
            clifford_tableau,
            connectivity,
            remaining_rows,
            pivot_column,
            pivot_row,
        ),
        _ => panic!("Invalid Pauli letter for observable cleaning"),
    }
}

pub(super) fn clean_x_prc<G>(
    repr: &mut G,
    clifford_tableau: &mut CliffordTableau,
    connectivity: &Connectivity,
    remaining_columns: &[usize],
    pivot_column: usize,
    pivot_row: usize,
) where
    G: CliffordGates,
{
    let mut terminals = remaining_columns
        .iter()
        .filter_map(|qubit| {
            if is_not_i(clifford_tableau.destabilizer(*qubit, pivot_row)) {
                Some(*qubit)
            } else {
                None
            }
        })
        .collect::<Vec<_>>();

    if terminals.is_empty() {
        return;
    }
    terminals.push(pivot_column);

    let traversal = connectivity
        .get_cx_ladder(&terminals, &pivot_column)
        .unwrap();

    let affected_cols = check_across_columns(&*clifford_tableau, &terminals, pivot_row, is_y);
    for col in affected_cols {
        repr.s(col);
        clifford_tableau.s(col);
    }

    let affected_cols = check_across_columns(&*clifford_tableau, &terminals, pivot_row, is_z);
    for col in affected_cols {
        repr.h(col);
        clifford_tableau.h(col);
    }

    for (parent, child) in traversal.iter().rev() {
        if is_i(clifford_tableau.destabilizer(*parent, pivot_row)) {
            repr.cx(*child, *parent);
            clifford_tableau.cx(*child, *parent);
        }
    }

    for (parent, child) in traversal.iter().rev() {
        repr.cx(*parent, *child);
        clifford_tableau.cx(*parent, *child);
    }
}

pub(super) fn clean_z_prc<G>(
    repr: &mut G,
    clifford_tableau: &mut CliffordTableau,
    connectivity: &Connectivity,
    remaining_columns: &[usize],
    pivot_column: usize,
    pivot_row: usize,
) where
    G: CliffordGates,
{
    let num_qubits = clifford_tableau.size();
    let mut terminals = remaining_columns
        .iter()
        .filter_map(|qubit| {
            if is_not_i(clifford_tableau.stabilizer(*qubit, pivot_row)) {
                Some(*qubit)
            } else {
                None
            }
        })
        .collect::<Vec<_>>();
    if terminals.is_empty() {
        return;
    }
    terminals.push(pivot_column);

    let traversal = connectivity
        .get_cx_ladder(&terminals, &pivot_column)
        .unwrap();

    let affected_cols =
        check_across_columns(&*clifford_tableau, &terminals, pivot_row + num_qubits, is_y);

    for col in affected_cols {
        repr.v(col);
        clifford_tableau.v(col);
    }

    let affected_cols =
        check_across_columns(&*clifford_tableau, &terminals, pivot_row + num_qubits, is_x);
    for col in affected_cols {
        repr.h(col);
        clifford_tableau.h(col);
    }

    for (parent, child) in traversal.iter().rev() {
        if is_i(clifford_tableau.stabilizer(*parent, pivot_row)) {
            repr.cx(*parent, *child);
            clifford_tableau.cx(*parent, *child);
        }
    }

    for (parent, child) in traversal.iter().rev() {
        repr.cx(*child, *parent);
        clifford_tableau.cx(*child, *parent);
    }
}
