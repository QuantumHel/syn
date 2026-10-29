use crate::{
    architecture::connectivity::{self, Connectivity},
    data_structures::{CliffordTableau, PauliLetter},
    ir::{
        clifford_tableau::helper::{clean_pivot, clean_prc, pick_column, pick_row},
        AdjointSynthesizer, CliffordGates,
    },
};

use super::helper::clean_signs;

/// Connectivity-aware tableau synthesis: PermRowCol.
///
/// Where [`NaiveCliffordSynthesizer`](super::NaiveCliffordSynthesizer) walks the qubits in index
/// order and emits `CX` wherever it likes, this picks each pivot from the current connectivity
/// graph, eliminates along architecture-conforming CNOT ladders, and then *removes* the finished
/// qubit from the graph. Removing it is what keeps later pivots honest: the remaining qubits must
/// stay connected among themselves, so the pivot heuristics prefer non-cutting vertices.
///
/// Because it is free to choose which logical qubit lands where, the output realizes the input
/// tableau only up to a final permutation — returned as the residual tableau, for the caller to
/// hand to
/// [`CliffordGates::add_final_permutation`](crate::ir::CliffordGates::add_final_permutation).
///
/// The pivot heuristics are swappable via [`set_row_strategy`](Self::set_row_strategy) and
/// [`set_column_strategy`](Self::set_column_strategy).
pub struct PermRowColCliffordSynthesizer {
    connectivity: Connectivity,
    row_strategy: fn(&CliffordTableau, &Connectivity, &[usize]) -> usize,
    column_strategy: fn(&CliffordTableau, &Connectivity, usize) -> usize,
}

impl PermRowColCliffordSynthesizer {
    /// Builds a synthesizer for the given hardware graph, with the default pivot heuristics.
    ///
    /// The graph must cover exactly the tableau's qubits and be connected; that is checked at
    /// synthesis time, via [`Connectivity::validate_for`].
    pub fn new(connectivity: Connectivity) -> Self {
        Self {
            connectivity,
            row_strategy: pick_row,
            column_strategy: pick_column,
        }
    }

    /// Replaces the heuristic that picks the next pivot row.
    ///
    /// Called with the tableau, the remaining connectivity, and the rows still to be eliminated;
    /// it returns one of those rows.
    pub fn set_row_strategy(
        &mut self,
        row_strategy: fn(&CliffordTableau, &Connectivity, &[usize]) -> usize,
    ) {
        (self.row_strategy) = row_strategy;
    }

    /// Replaces the heuristic that picks the pivot column for a chosen row.
    pub fn set_column_strategy(
        &mut self,
        column_strategy: fn(&CliffordTableau, &Connectivity, usize) -> usize,
    ) {
        (self.column_strategy) = column_strategy;
    }

    /// The hardware graph this synthesizer targets.
    pub fn get_connectivity(&self) -> &Connectivity {
        &self.connectivity
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
        // `num_qubits <= node_count()` used to be the check here, but a *wider* graph is just as
        // broken as a narrower one: `pick_row` reads `clifford_tableau.stabilizer(qubit, ..)` for
        // every graph node, so a node beyond the tableau panics on the tableau's own bounds check.
        // What synthesis actually needs is qubits numbered exactly `0..num_qubits`.
        if let Err(error) = self.connectivity.validate_for(num_qubits) {
            panic!("{error}");
        }
        let mut connectivity = self.connectivity.clone();
        // logical qubit remaining to be disconnected
        let mut remaining_columns = (0..num_qubits).collect::<Vec<_>>();
        // stabilizers / destabilizers that are not yet identity rows
        let mut remaining_rows = (0..num_qubits).collect::<Vec<_>>();

        while !remaining_columns.is_empty() {
            let pivot_row = (self.row_strategy)(&clifford_tableau, &connectivity, &remaining_rows);
            let pivot_column = (self.column_strategy)(&clifford_tableau, &connectivity, pivot_row);
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
                &connectivity,
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
                &connectivity,
                &remaining_columns,
                pivot_column,
                pivot_row,
                second_letter,
            );

            // If the pivot row is now an identity row, we can remove it from the tableau.
            connectivity.remove_node(pivot_column);
        }

        clean_signs(repr, &mut clifford_tableau);
        clifford_tableau
    }
}
