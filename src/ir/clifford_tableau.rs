use std::iter::zip;

use itertools::{iproduct, Itertools};
use petgraph::algo::connected_components;

use crate::{
    architecture::{connectivity::Connectivity, Architecture},
    data_structures::{CliffordTableau, PauliLetter, PauliString, PropagateClifford},
    synthesis_methods::{architectureaware::PermRowCol, custom::Custom, naive::Naive},
};

use super::CliffordGates;

fn get_pauli(pauli_string: &PauliString, row: usize) -> PauliLetter {
    PauliLetter::new(pauli_string.x(row), pauli_string.z(row))
}

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

fn is_not_y(pauli_letter: PauliLetter) -> bool {
    pauli_letter != PauliLetter::Y
}

fn is_z(pauli_letter: PauliLetter) -> bool {
    pauli_letter == PauliLetter::Z
}

fn is_not_z(pauli_letter: PauliLetter) -> bool {
    pauli_letter != PauliLetter::Z
}

pub struct CliffordTableauSynthesizer;
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

impl<G, C> PermRowCol<CliffordTableau, G, C> for CliffordTableauSynthesizer
where
    G: CliffordGates,
    C: Architecture,
{
    fn run_prc(
        clifford_tableau: &CliffordTableau,
        repr: &mut G,
        mut connectivity: C,
        pick_row: Option<fn(&C, &[usize], &CliffordTableau) -> usize>,
        pick_column: Option<fn(&C, usize, &[usize], &CliffordTableau) -> usize>,
    ) {
        let mut ct = clifford_tableau.adjoint();

        let num_qubits = ct.size();

        let mut remaining_columns = (0..num_qubits).collect::<Vec<_>>();
        let mut remaining_rows = (0..num_qubits).collect::<Vec<_>>();

        // let mut qubit_mapping = (0..num_qubits).collect::<Vec<_>>();
        // let mut permutation = (0..num_qubits).collect::<Vec<_>>();

        // Identify row that causes the least amount of fill-in
        let pick_pivot_row: fn(&C, &[usize], &CliffordTableau) -> usize =
            if let Some(pivot_row) = pick_row {
                pivot_row
            } else {
                default_pick_row
            };
        // Identify column that is easiest to clean up
        let pick_pivot_col: fn(&C, usize, &[usize], &CliffordTableau) -> usize =
            if let Some(pivot_col) = pick_column {
                pivot_col
            } else {
                default_pick_column
            };

        let mut column_ordering = Vec::new();
        let mut row_ordering = Vec::new();
        while !connectivity.nodes().is_empty() {
            let pivot_row = pick_pivot_row(&connectivity, &remaining_rows, &ct);
            let pivot_column = pick_pivot_col(&connectivity, pivot_row, &remaining_columns, &ct);
            column_ordering.push(pivot_column);
            row_ordering.push(pivot_row);
            // let non_cutting = connectivity.non_cutting();

            // if !non_cutting.contains(&pivot_column) {
            //     // Pick nearest non-cutting node -> How do we optimize this? Could we randomize this as well?
            //     let non_cutting_node = non_cutting
            //         .iter()
            //         .map(|node| connectivity.distance(*node, pivot_column))
            //         .min()
            //         .unwrap();
            //     (qubit_mapping[pivot_column], qubit_mapping[non_cutting_node]) =
            //         (qubit_mapping[non_cutting_node], qubit_mapping[pivot_column]);
            // }
            connectivity.remove_node(pivot_column);
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

        let column_permutation = zip((0..num_qubits).collect::<Vec<_>>(), column_ordering)
            .sorted_by_key(|a| a.1)
            .map(|a| a.0)
            .collect::<Vec<_>>();

        let row_permutation = zip((0..num_qubits).collect::<Vec<_>>(), row_ordering)
            .sorted_by_key(|a| a.1)
            .map(|a| a.0)
            .collect::<Vec<_>>();
        println!("column_permutation: {:?}", column_permutation);
        println!("row_permutation: {:?}", row_permutation);
        println!("ct after: {}", ct);
        clean_signs(repr, &mut ct, &row_permutation);
    }

    // for (&pivot_column, &pivot_row) in zip(&custom_columns, &custom_rows) {
    //     // Cleanup pivot column
}

fn clean_pivot<G>(repr: &mut G, ct: &mut CliffordTableau, pivot_column: usize, pivot_row: usize)
where
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

fn clean_x_pivot<G>(repr: &mut G, ct: &mut CliffordTableau, pivot_column: usize, pivot_row: usize)
where
    G: CliffordGates,
{
    // These are switched around because of implementation
    if check_pauli(&*ct, pivot_row, pivot_column, is_y) {
        ct.s(pivot_row);
        repr.s(pivot_row);
    }

    // These are switched around because of implementation
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
    // These are switched around because of implementation
    if check_pauli(&*ct, pivot_row, pivot_column + num_qubits, is_y) {
        ct.v(pivot_row);
        repr.v(pivot_row);
    }

    // These are switched around because of implementation
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

/// Pick a non-cutting vertex in connectivity that requires the fewest operations to remove from the graph
/// Ensures that picked row is valid for Gaussian Elimination.
fn default_pick_column<C>(
    connectivity: &C,
    pivot_row: usize,
    remaining_cols: &[usize],
    clifford_tableau: &CliffordTableau,
) -> usize
where
    C: Architecture,
{
    let mut scores = Vec::with_capacity(remaining_cols.len());

    for &col in connectivity.non_cutting() {
        // If both stabilizer and destabilizer is I for this column, skip.
        if check_pauli(clifford_tableau, pivot_row, col, is_i)
            && check_pauli(
                clifford_tableau,
                pivot_row,
                col + clifford_tableau.size(),
                is_i,
            )
        {
            continue;
        }
        let mut score = 0;
        for row in connectivity.nodes() {
            let node_distance = connectivity.distance(col, row);
            if check_pauli(clifford_tableau, row, col, is_not_i) {
                score += node_distance;
            }
            if check_pauli(
                clifford_tableau,
                row,
                col + clifford_tableau.size(),
                is_not_i,
            ) {
                score += node_distance;
            }
        }
        scores.push((col, score));
    }
    scores.iter().min_by_key(|a| a.1).unwrap().0
}

/// Pick row for cancellation that results in the least fill-in
fn default_pick_row<C>(
    connectivity: &C,
    remaining_rows: &[usize],
    clifford_tableau: &CliffordTableau,
) -> usize
where
    C: Architecture,
{
    let mut scores = Vec::with_capacity(remaining_rows.len());
    for &row in remaining_rows {
        let mut score = 0;
        for col in connectivity.nodes() {
            if check_pauli(clifford_tableau, row, col, is_not_i) {
                score += 1
            }

            if check_pauli(
                clifford_tableau,
                row,
                col + clifford_tableau.size(),
                is_not_i,
            ) {
                score += 1
            }
        }
        scores.push((row, score));
    }
    scores.iter().min_by_key(|a| a.1).unwrap().0
}

#[cfg(test)]
mod tests {
    use crate::{
        architecture::connectivity::Connectivity,
        data_structures::{CliffordTableau, PropagateClifford},
        ir::clifford_tableau::default_pick_column,
    };

    use super::default_pick_row;

    #[test]
    fn test_pick_column() {
        let num_qubits = 5;
        let weighted_edges = [(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1), (2, 4, 1)];
        let connectivity = Connectivity::from_weighted_edges(&weighted_edges);

        let mut clifford_tableau = CliffordTableau::new(num_qubits);
        clifford_tableau.cx(0, 1);
        clifford_tableau.cx(2, 3);
        clifford_tableau.cx(3, 4);
        clifford_tableau.cx(1, 2);
        let remaining_columns = [0, 1, 2, 3, 4];
        assert_eq!(
            default_pick_column(&connectivity, 1, &remaining_columns, &clifford_tableau),
            0
        );
    }

    #[test]
    fn test_pick_row() {
        let num_qubits = 5;
        let weighted_edges = [(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1), (2, 4, 1)];
        let connectivity = Connectivity::from_weighted_edges(&weighted_edges);

        let mut clifford_tableau = CliffordTableau::new(num_qubits);
        clifford_tableau.cx(0, 1);
        clifford_tableau.cx(2, 3);
        clifford_tableau.cx(3, 4);
        clifford_tableau.cx(1, 2);
        let remaining_rows = [0, 1, 2, 3, 4];

        assert_eq!(
            default_pick_row(&connectivity, &remaining_rows, &clifford_tableau),
            0
        );
    }
}
