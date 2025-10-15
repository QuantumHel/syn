use crate::{
    architecture::connectivity::Connectivity,
    data_structures::{CliffordTableau, PauliLetter},
    ir::{
        clifford_tableau::helper::{clean_pivot, clean_prc, pick_column, pick_row},
        AdjointSynthesizer, CliffordGates,
    },
};

use super::helper::clean_signs;

// #[derive(Default)]
pub struct PermRowColCliffordSynthesizer {
    connectivity: Connectivity,
    row_strategy: fn(&CliffordTableau, &Connectivity, &[usize]) -> usize,
    column_strategy: fn(&CliffordTableau, &Connectivity, usize) -> usize,
}

impl PermRowColCliffordSynthesizer {
    pub fn new(connectivity: Connectivity) -> Self {
        let size = connectivity.node_count();

        Self {
            connectivity,
            row_strategy: pick_row,
            column_strategy: pick_column,
        }
    }

    pub fn set_row_strategy(
        &mut self,
        row_strategy: fn(&CliffordTableau, &Connectivity, &[usize]) -> usize,
    ) {
        (self.row_strategy) = row_strategy;
    }

    pub fn set_column_strategy(
        &mut self,
        column_strategy: fn(&CliffordTableau, &Connectivity, usize) -> usize,
    ) {
        (self.column_strategy) = column_strategy;
    }
}

impl<G> AdjointSynthesizer<CliffordTableau, G, CliffordTableau> for PermRowColCliffordSynthesizer
where
    G: CliffordGates,
{
    fn synthesize_adjoint(
        &mut self,
        mut clifford_tableau: CliffordTableau,
        repr: &mut G,
    ) -> CliffordTableau {
        let num_qubits = clifford_tableau.size();
        let machine_size = self.connectivity.node_count();
        assert!(
            num_qubits <= machine_size,
            "Number of qubits {} exceeds machine size {}",
            num_qubits,
            machine_size
        );
        // logical qubit remaining to be disconnected
        let mut remaining_columns = (0..num_qubits).collect::<Vec<_>>();
        // stabilizers / destabilizers that are not yet identity rows
        let mut remaining_rows = (0..num_qubits).collect::<Vec<_>>();

        while !remaining_columns.is_empty() {
            let pivot_row =
                (self.row_strategy)(&clifford_tableau, &self.connectivity, &remaining_rows);
            let pivot_column =
                (self.column_strategy)(&clifford_tableau, &self.connectivity, pivot_row);
            let column = clifford_tableau.column(pivot_column);
            let x_weight = column.x_weight();
            let z_weight = column.z_weight();

            let (first_letter, second_letter) = if z_weight > x_weight {
                (PauliLetter::Z, PauliLetter::X)
            } else {
                (PauliLetter::X, PauliLetter::Z)
            };

            remaining_columns.retain(|&x| x != pivot_column);
            remaining_rows.retain(|&x| x != pivot_row);

            clean_pivot(
                repr,
                &mut clifford_tableau,
                pivot_column,
                pivot_row,
                first_letter,
            );

            // Use the pivot to remove all other terms in the X observable.
            clean_prc(
                repr,
                &mut clifford_tableau,
                &self.connectivity,
                &remaining_columns,
                pivot_column,
                pivot_row,
                first_letter,
            );

            clean_pivot(
                repr,
                &mut clifford_tableau,
                pivot_column,
                pivot_row,
                second_letter,
            );

            // Use the pivot to remove all other terms in the Z observable.
            clean_prc(
                repr,
                &mut clifford_tableau,
                &self.connectivity,
                &remaining_columns,
                pivot_column,
                pivot_row,
                second_letter,
            );

            // If the pivot row is now an identity row, we can remove it from the tableau.
            self.connectivity.remove_node(pivot_column);
        }

        clean_signs(repr, &mut clifford_tableau);
        return clifford_tableau;
    }
}
