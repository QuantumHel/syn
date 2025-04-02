use crate::{
    data_structures::CliffordTableau,
    ir::{AdjointSynthesizer, CliffordGates, HasAdjoint},
};

use super::helper::{
    clean_pivot, clean_signs, clean_x_observables, clean_z_observables, naive_pivot_search, swap,
};

#[derive(Default)]
pub struct NaiveCliffordSynthesizer {}

impl<G> AdjointSynthesizer<CliffordTableau, G> for NaiveCliffordSynthesizer
where
    G: CliffordGates,
{
    fn synthesize_adjoint(&mut self, clifford_tableau: CliffordTableau, repr: &mut G) {
        let mut clifford_tableau = clifford_tableau.adjoint();
        let num_qubits = clifford_tableau.size();

        for row in 0..num_qubits {
            let pivot_col = naive_pivot_search(&clifford_tableau, num_qubits, row);

            if pivot_col != row {
                swap(repr, &mut clifford_tableau, row, pivot_col);
            }

            // Cleanup pivot column
            clean_pivot(repr, &mut clifford_tableau, row, row);

            let checked_rows = (row + 1..num_qubits).collect::<Vec<_>>();

            // Use the pivot to remove all other terms in the X observable.
            clean_x_observables(repr, &mut clifford_tableau, &checked_rows, row, row);

            // Use the pivot to remove all other terms in the Z observable.
            clean_z_observables(repr, &mut clifford_tableau, &checked_rows, row, row);
        }

        clean_signs(
            repr,
            &mut clifford_tableau,
            &(0..num_qubits).collect::<Vec<_>>(),
        );
    }
}
