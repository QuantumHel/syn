use crate::{
    data_structures::{CliffordTableau, PauliString, PropagateClifford},
    synthesis_methods::naive::Naive,
};

use super::CliffordGates;

fn get_pauli(pauli_string: &PauliString, i: usize) -> usize {
    pauli_string.x(i) as usize + pauli_string.z(i) as usize
}
pub struct CliffordTableauSynthesizer;
impl<G> Naive<CliffordTableau, G> for CliffordTableauSynthesizer
where
    G: CliffordGates,
{
    fn run(mut clifford_tableau: CliffordTableau, mut repr: G) {
        let n = clifford_tableau.size();
        for col in 0..n {
            let mut pivot_row = col;
            let mut x_pauli;
            let mut z_pauli;
            {
                let column = clifford_tableau.column(col);
                for row in 0..n {
                    x_pauli = get_pauli(column, row);
                    z_pauli = get_pauli(column, row + n);
                    if x_pauli != 0 && z_pauli != 0 && x_pauli != z_pauli {
                        pivot_row = row;
                        break;
                    }
                }
            }

            if pivot_row != col {
                repr.cx(pivot_row, col);
                repr.cx(col, pivot_row);
                repr.cx(pivot_row, col);

                clifford_tableau.cx(pivot_row, col);
                clifford_tableau.cx(col, pivot_row);
                clifford_tableau.cx(pivot_row, col);
            }

            {
                let column = clifford_tableau.column(col);
                x_pauli = get_pauli(column, col);
                z_pauli = get_pauli(column, col + n);
            }

            // Transform the pivot to the XZ plane
            if z_pauli == 3 {
                repr.s(col);
                clifford_tableau.s(col);
            }

            if z_pauli == 2 {
                repr.h(col);
                clifford_tableau.h(col);
            }

            if x_pauli == 2 {
                repr.s(col);
                clifford_tableau.s(col);
            }

            // Use the pivot to remove all other terms in the X observable.
            let affected_rows = check_rows(&clifford_tableau, col, n, |p| p == 3);
            for row in affected_rows {
                clifford_tableau.s(row);
            }

            let affected_rows = check_rows(&clifford_tableau, col, n, |p| p == 2);
            for row in affected_rows {
                repr.h(row);
                clifford_tableau.h(row);
            }

            let affected_rows = check_rows(&clifford_tableau, col, n, |p| p != 0);
            for row in affected_rows {
                repr.cx(col, row);
                clifford_tableau.cx(col, row);
            }

            // Use the pivot to remove all other terms in the Z observable.
            let affected_rows = check_rows(&clifford_tableau, col + n, 2 * n, |p| p == 3);
            for row in affected_rows {
                repr.s(row);
                clifford_tableau.s(row);
            }

            let affected_rows = check_rows(&clifford_tableau, col + n, 2 * n, |p| p == 1);
            for row in affected_rows {
                repr.h(row);
                clifford_tableau.h(row);
            }

            let affected_rows = check_rows(&clifford_tableau, col + n, 2 * n, |p| p != 0);
            for row in affected_rows {
                repr.cx(row - n, col);
                clifford_tableau.cx(row - n, col);
            }
        }

        let z_signs = clifford_tableau.z_signs();

        for (col, sign) in z_signs.iter().enumerate() {
            if *sign {
                repr.h(col);
                repr.s(col);
                repr.s(col);
                repr.h(col);
                clifford_tableau.h(col);
                clifford_tableau.s(col);
                clifford_tableau.s(col);
                clifford_tableau.h(col);
            }
        }

        let x_signs = clifford_tableau.x_signs();

        for (col, sign) in x_signs.iter().enumerate() {
            if *sign {
                repr.s(col);
                repr.s(col);
                clifford_tableau.s(col);
                clifford_tableau.s(col);
            }
        }
    }
}

/// Helper function that returns indices for a particular column `col` of `clifford_tableau` that match the provided closure `pauli_check`.
fn check_rows(
    clifford_tableau: &CliffordTableau,
    col: usize,
    n: usize,
    pauli_check: fn(usize) -> bool,
) -> Vec<usize> {
    let column = clifford_tableau.column(col);
    let mut affected_rows = Vec::new();
    for row in (col + 1)..n {
        if pauli_check(get_pauli(column, row)) {
            affected_rows.push(row);
        }
    }
    affected_rows
}
