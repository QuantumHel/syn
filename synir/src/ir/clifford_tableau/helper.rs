use core::num;
use std::{borrow::BorrowMut, iter::zip};

use bitvec::prelude::{BitVec, Msb0};
use bitvec::vec;
use itertools::Itertools;

use crate::{
    architecture::{connectivity::Connectivity, Architecture},
    data_structures::{CliffordTableau, PauliLetter, PauliString, PropagateClifford},
    ir::CliffordGates,
};

pub(super) fn clean_naive_pivot<G>(
    repr: &mut G,
    ct: &mut CliffordTableau,
    pivot_column: usize,
    pivot_row: usize,
) where
    G: CliffordGates,
{
    let num_qubits = ct.size();
    match ct.column(pivot_row).pauli(pivot_column + num_qubits) {
        PauliLetter::Y => {
            ct.s(pivot_row);
            repr.s(pivot_row);
        }
        PauliLetter::X => {
            ct.v(pivot_row);
            repr.v(pivot_row);
        }
        PauliLetter::Z => {
            ct.h(pivot_row);
            repr.h(pivot_row);
        }
        PauliLetter::I => (),
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

/// Sets destabilizer entry at (pivot_column, pivot_row) to X if it is not I, leaves I terms unchanged.
pub(super) fn clean_x_pivot<G>(
    repr: &mut G,
    clifford_tableau: &mut CliffordTableau,
    pivot_column: usize,
    pivot_row: usize,
) where
    G: CliffordGates,
{
    match clifford_tableau.column(pivot_column).pauli(pivot_row) {
        PauliLetter::Y => {
            clifford_tableau.s(pivot_column);
            repr.s(pivot_column);
        }
        PauliLetter::Z => {
            clifford_tableau.h(pivot_column);
            repr.h(pivot_column);
        }
        _ => (),
    }
}

/// Sets destabilizer entry at (pivot_column, pivot_row) to Z if it is not I, leaves I terms unchanged.
pub(super) fn clean_z_pivot<G>(
    repr: &mut G,
    clifford_tableau: &mut CliffordTableau,
    pivot_column: usize,
    pivot_row: usize,
) where
    G: CliffordGates,
{
    let num_qubits = clifford_tableau.size();
    match clifford_tableau
        .column(pivot_column)
        .pauli(pivot_row + num_qubits)
    {
        PauliLetter::Y => {
            clifford_tableau.v(pivot_column);
            repr.v(pivot_column);
        }
        PauliLetter::X => {
            clifford_tableau.h(pivot_column);
            repr.h(pivot_column);
        }
        _ => (),
    }
}

/// Cleans Ithe destabilizer observables for `pivot_row` in the Clifford tableau using (pivot_column, pivot_row) as the entry for elimination..
/// Assumes that (pivot_column, pivot_row) is either an I term or a X term.
/// Only removes entries from columns in `remaining_columns` and assumes `pivot_column` has already been removed from `remaining_columns`.
/// If (pivot_column, pivot_row) is an I term, set it to X first using a non-I term in (pivot_row, remaining_columns).
pub(super) fn clean_x_observables<G>(
    repr: &mut G,
    clifford_tableau: &mut CliffordTableau,
    remaining_columns: &[usize],
    pivot_column: usize,
    pivot_row: usize,
) where
    G: CliffordGates,
{
    let affected_cols = remaining_columns
        .iter()
        .filter(|col| clifford_tableau.column(**col).pauli(pivot_row) == PauliLetter::Y)
        .collect_vec();
    for col in affected_cols {
        repr.s(*col);
        clifford_tableau.s(*col);
    }

    let affected_cols = remaining_columns
        .iter()
        .filter(|col| clifford_tableau.column(**col).pauli(pivot_row) == PauliLetter::Z)
        .collect_vec();

    for col in affected_cols {
        repr.h(*col);
        clifford_tableau.h(*col);
    }

    let affected_cols = remaining_columns
        .iter()
        .filter(|col| clifford_tableau.column(**col).pauli(pivot_row) != PauliLetter::I)
        .collect_vec();
    if clifford_tableau.column(pivot_column).pauli(pivot_row) == PauliLetter::I {
        repr.cx(*affected_cols[0], pivot_column);
        clifford_tableau.cx(*affected_cols[0], pivot_column);
    }

    for col in affected_cols {
        repr.cx(pivot_column, *col);
        clifford_tableau.cx(pivot_column, *col);
    }
}

/// Cleans the destabilizer observables for `pivot_row` in the Clifford tableau using (pivot_column, pivot_row) as the entry for elimination..
/// Assumes that (pivot_column, pivot_row) is either an I term or a X term.
/// Only removes entries from columns in `remaining_columns` and assumes `pivot_column` has already been removed from `remaining_columns`.
/// If (pivot_column, pivot_row) is an I term, set it to X first using a non-I term in (pivot_row, remaining_columns).
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
    let affected_cols = remaining_columns
        .iter()
        .filter(|col| {
            clifford_tableau.column(**col).pauli(pivot_row + num_qubits) == PauliLetter::Y
        })
        .collect_vec();
    for col in affected_cols {
        repr.v(*col);
        clifford_tableau.v(*col);
    }

    let affected_cols = remaining_columns
        .iter()
        .filter(|col| {
            clifford_tableau.column(**col).pauli(pivot_row + num_qubits) == PauliLetter::X
        })
        .collect_vec();
    for col in affected_cols {
        repr.h(*col);
        clifford_tableau.h(*col);
    }

    let affected_cols = remaining_columns
        .iter()
        .filter(|col| {
            clifford_tableau.column(**col).pauli(pivot_row + num_qubits) != PauliLetter::I
        })
        .collect_vec();
    if clifford_tableau
        .column(pivot_column)
        .pauli(pivot_row + num_qubits)
        == PauliLetter::I
    {
        repr.cx(pivot_column, *affected_cols[0]);
        clifford_tableau.cx(pivot_column, *affected_cols[0]);
    }

    for col in affected_cols {
        repr.cx(*col, pivot_column);
        clifford_tableau.cx(*col, pivot_column);
    }
}

pub(super) fn clean_signs<G>(repr: &mut G, clifford_tableau: &mut CliffordTableau)
where
    G: CliffordGates,
{
    let n = clifford_tableau.size();
    let inv_perm = match clifford_tableau.get_permutation() {
        None => panic!(
            "Cleaning signs but tableau is not a permutation matrix: \n{}",
            clifford_tableau
        ),
        Some(perm) => perm,
    };
    let row_permutation = (0..n)
        .map(|i| inv_perm.iter().find_position(|&&x| x == i))
        .map(|x| x.unwrap().0)
        .collect_vec();
    let z_signs = clifford_tableau.z_signs();
    for r in (0..n).filter(|r| z_signs[*r]) {
        let row = row_permutation[r];
        repr.x(row);
        clifford_tableau.x(row);
        println!("X on {}", row);
    }
    assert_eq!(
        clifford_tableau.z_signs(),
        BitVec::<u8, Msb0>::repeat(false, n)
    );
    let x_signs = clifford_tableau.x_signs();
    for r in (0..n).filter(|r| x_signs[*r]) {
        let row = row_permutation[r];
        repr.z(row);
        clifford_tableau.z(row);
    }
    assert_eq!(
        clifford_tableau.x_signs(),
        BitVec::<u8, Msb0>::repeat(false, n)
    );
}

pub(super) fn naive_pivot_search(
    clifford_tableau: &CliffordTableau,
    num_qubits: usize,
    row: usize,
) -> usize {
    let mut pivot_col = 0;

    for col in 0..num_qubits {
        let column = clifford_tableau.column(col);
        let x_pauli = column.pauli(row);
        let z_pauli = column.pauli(row + num_qubits);
        if x_pauli != PauliLetter::I && z_pauli != PauliLetter::I && x_pauli != z_pauli {
            pivot_col = col;
            break;
        }
    }
    pivot_col
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
            if clifford_tableau.stabilizer(qubit, *row) != PauliLetter::I {
                row_weights[*row] += 1;
            }
            if clifford_tableau.destabilizer(qubit, *row) != PauliLetter::I {
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
            } else {
                column_weights[*qubit] +=
                    (clifford_tableau.stabilizer(*qubit, interaction) == PauliLetter::I) as usize;
                column_weights[*qubit] +=
                    (clifford_tableau.destabilizer(*qubit, interaction) == PauliLetter::I) as usize;
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
            if clifford_tableau.destabilizer(*qubit, pivot_row) != PauliLetter::I {
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

    let affected_cols = terminals
        .iter()
        .filter(|col| clifford_tableau.column(**col).pauli(pivot_row) == PauliLetter::Y)
        .collect_vec();
    for col in affected_cols {
        repr.s(*col);
        clifford_tableau.s(*col);
    }

    let affected_cols = terminals
        .iter()
        .filter(|col| clifford_tableau.column(**col).pauli(pivot_row) == PauliLetter::Z)
        .collect_vec();
    for col in affected_cols {
        repr.h(*col);
        clifford_tableau.h(*col);
    }

    for (parent, child) in traversal.iter().rev() {
        if clifford_tableau.destabilizer(*parent, pivot_row) == PauliLetter::I {
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
            if clifford_tableau.stabilizer(*qubit, pivot_row) != PauliLetter::I {
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

    let affected_cols = terminals
        .iter()
        .filter(|col| {
            clifford_tableau.column(**col).pauli(pivot_row + num_qubits) == PauliLetter::Y
        })
        .collect_vec();
    for col in affected_cols {
        repr.v(*col);
        clifford_tableau.v(*col);
    }

    let affected_cols = terminals
        .iter()
        .filter(|col| {
            clifford_tableau.column(**col).pauli(pivot_row + num_qubits) == PauliLetter::X
        })
        .collect_vec();
    for col in affected_cols {
        repr.h(*col);
        clifford_tableau.h(*col);
    }

    for (parent, child) in traversal.iter().rev() {
        if clifford_tableau.stabilizer(*parent, pivot_row) == PauliLetter::I {
            repr.cx(*parent, *child);
            clifford_tableau.cx(*parent, *child);
        }
    }

    for (parent, child) in traversal.iter().rev() {
        repr.cx(*child, *parent);
        clifford_tableau.cx(*child, *parent);
    }
}
