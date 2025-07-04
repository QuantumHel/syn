use std::fmt::Debug;

use crate::{
    architecture::connectivity::Connectivity,
    data_structures::{CliffordTableau, PauliLetter},
    ir::{
        clifford_tableau::{
            self,
            helper::{
                clean_observables, clean_pivot, clean_prc, clean_x_pivot, clean_z_pivot,
                pick_column, pick_row,
            },
        },
        AdjointSynthesizer, CliffordGates,
    },
};

use super::helper::{
    clean_naive_pivot, clean_signs, clean_x_observables, clean_z_observables, naive_pivot_search,
    swap,
};

#[derive(Default)]
pub struct PermRowColCliffordSynthesizer {
    connectivity: Connectivity,
    permutation: Vec<usize>,
}

impl PermRowColCliffordSynthesizer {
    pub fn new(connectivity: Connectivity) -> Self {
        let size = connectivity.node_bound();

        Self {
            connectivity,
            permutation: (0..size).collect(),
        }
    }

    pub fn permutation(&self) -> &[usize] {
        &self.permutation
    }
}

impl<G> AdjointSynthesizer<CliffordTableau, G> for PermRowColCliffordSynthesizer
where
    G: CliffordGates + Debug,
{
    fn synthesize_adjoint(&mut self, mut clifford_tableau: CliffordTableau, repr: &mut G) {
        let num_qubits = clifford_tableau.size();
        let machine_size = self.connectivity.node_bound();
        assert!(
            num_qubits <= machine_size,
            "Number of qubits {} exceeds machine size {}",
            num_qubits,
            machine_size
        );
        // Mapping between logical qubit to physical qubit
        let mut permutation = (0..num_qubits).collect::<Vec<_>>();
        // logical qubit remaining to be disconnected
        let mut remaining_columns = (0..num_qubits).collect::<Vec<_>>();
        // stabilizers / destabilizers that are not yet identity rows
        let mut remaining_rows = (0..num_qubits).collect::<Vec<_>>();

        while !remaining_columns.is_empty() {
            let pivot_row = pick_row(&clifford_tableau, &self.connectivity, &remaining_rows);
            let pivot_column = pick_column(&clifford_tableau, &self.connectivity);
            println!("pivot row: {}", pivot_row);
            println!("pivot column: {}", pivot_column);
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
            {
                clean_pivot(
                    repr,
                    &mut clifford_tableau,
                    pivot_column,
                    pivot_row,
                    first_letter,
                );
                println!("Cleaning: {:?}", first_letter);
                println!("repr: {:?}", repr);
                println!("clifford_tableau: {}", clifford_tableau);

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
                println!("PRC: {:?}", first_letter);
                println!("repr: {:?}", repr);
                println!("clifford_tableau: {}", clifford_tableau);

                clean_pivot(
                    repr,
                    &mut clifford_tableau,
                    pivot_column,
                    pivot_row,
                    second_letter,
                );
                println!("Cleaning: {:?}", second_letter);
                println!("repr: {:?}", repr);
                println!("clifford_tableau: {}", clifford_tableau);

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
                println!("PRC: {:?}", second_letter);
                println!("repr: {:?}", repr);
                println!("clifford_tableau: {}", clifford_tableau);
            }

            // If the pivot row is now an identity row, we can remove it from the tableau.

            permutation[pivot_row] = pivot_column;
            self.connectivity.remove_node(pivot_column);
            println!("nodes: {:?}", self.connectivity.nodes());
        }

        clean_signs(repr, &mut clifford_tableau, &permutation);

        self.permutation = permutation;
    }
}
