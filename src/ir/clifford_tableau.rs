use crate::{
    data_structures::{CliffordTableau, PauliString, PropagateClifford},
    synthesis_methods::naive::Naive,
};

use super::CliffordGates;

fn get_pauli(pauli_string: &PauliString, row: usize) -> usize {
    pauli_string.x(row) as usize + pauli_string.z(row) as usize * 2
}
pub struct CliffordTableauSynthesizer;
impl<G> Naive<&CliffordTableau, G> for CliffordTableauSynthesizer
where
    G: CliffordGates,
{
    fn run(clifford_tableau: &CliffordTableau, repr: &mut G) {
        let mut clifford_tableau = clifford_tableau.adjoint();

        let num_qubits = clifford_tableau.size();
        for row in 0..num_qubits {
            let mut pivot_col = row;
            let mut x_pauli;
            let mut z_pauli;
            {
                for col in 0..num_qubits {
                    let column = clifford_tableau.column(col);
                    x_pauli = get_pauli(column, row);
                    z_pauli = get_pauli(column, row + num_qubits);
                    if x_pauli != 0 && z_pauli != 0 && x_pauli != z_pauli {
                        pivot_col = col;
                        break;
                    }
                }
            }

            if pivot_col != row {
                repr.cx(pivot_col, row);
                repr.cx(row, pivot_col);
                repr.cx(pivot_col, row);

                clifford_tableau.cx(pivot_col, row);
                clifford_tableau.cx(row, pivot_col);
                clifford_tableau.cx(pivot_col, row);
            }

            {
                let column = clifford_tableau.column(row);
                z_pauli = get_pauli(column, row + num_qubits);
            }

            // Transform the pivot to the XZ plane
            if z_pauli == 3 {
                repr.s(row);
                clifford_tableau.s(row);
            }

            {
                let column = clifford_tableau.column(row);
                z_pauli = get_pauli(column, row + num_qubits);
            }

            if z_pauli != 2 {
                repr.h(row);
                clifford_tableau.h(row);
            }

            {
                let column = clifford_tableau.column(row);
                x_pauli = get_pauli(column, row);
            }

            if x_pauli != 1 {
                repr.s(row);
                clifford_tableau.s(row);
            }

            // Use the pivot to remove all other terms in the X observable.
            let affected_cols = check_x_cols(&clifford_tableau, row, |p| p == 3);
            for col in affected_cols {
                repr.s(col);
                clifford_tableau.s(col);
            }

            let affected_cols = check_x_cols(&clifford_tableau, row, |p| p == 2);
            for col in affected_cols {
                repr.h(col);
                clifford_tableau.h(col);
            }

            let affected_cols = check_x_cols(&clifford_tableau, row, |p| p != 0);
            for col in affected_cols {
                repr.cx(row, col);
                clifford_tableau.cx(row, col);
            }

            // Use the pivot to remove all other terms in the Z observable.
            let affected_cols = check_z_cols(&clifford_tableau, row, |p| p == 3);
            for col in affected_cols {
                repr.s(col);
                clifford_tableau.s(col);
            }

            let affected_cols = check_z_cols(&clifford_tableau, row, |p| p == 1);
            for col in affected_cols {
                repr.h(col);
                clifford_tableau.h(col);
            }

            let affected_cols = check_z_cols(&clifford_tableau, row, |p| p != 0);
            for col in affected_cols {
                repr.cx(col, row);
                clifford_tableau.cx(col, row);
            }
        }

        let z_signs = clifford_tableau.z_signs();

        for (row, sign) in z_signs.iter().enumerate() {
            if *sign {
                repr.x(row);
                clifford_tableau.x(row);
            }
        }

        let x_signs = clifford_tableau.x_signs();

        for (row, sign) in x_signs.iter().enumerate() {
            if *sign {
                repr.z(row);
                clifford_tableau.z(row);
            }
        }
    }
}

/// Helper function that returns indices for a particular column `col` of `clifford_tableau` that match the provided closure `pauli_check`.
fn check_x_cols(
    clifford_tableau: &CliffordTableau,
    row: usize,
    pauli_check: fn(usize) -> bool,
) -> Vec<usize> {
    let mut affected_cols = Vec::new();
    let num_qubits = clifford_tableau.size();
    for col in row + 1..num_qubits {
        let pauli_string = clifford_tableau.column(col);
        if pauli_check(get_pauli(pauli_string, row)) {
            affected_cols.push(col);
        }
    }
    affected_cols
}

/// Helper function that returns indices for a particular column `col` of `clifford_tableau` that match the provided closure `pauli_check`.
fn check_z_cols(
    clifford_tableau: &CliffordTableau,
    row: usize,
    pauli_check: fn(usize) -> bool,
) -> Vec<usize> {
    let mut affected_cols = Vec::new();
    let num_qubits = clifford_tableau.size();
    for col in row + 1..num_qubits {
        let pauli_string = clifford_tableau.column(col);

        if pauli_check(get_pauli(pauli_string, row + num_qubits)) {
            affected_cols.push(col);
        }
    }
    affected_cols
}
