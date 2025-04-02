use std::iter::zip;

use itertools::Itertools;

use crate::{
    data_structures::CliffordTableau,
    ir::{
        clifford_tableau::helper::{clean_signs, clean_x_pivot, clean_z_pivot},
        CliffordGates, HasAdjoint,
    },
};

use super::helper::{clean_x_observables, clean_z_observables};

pub struct CustomCallbackCliffordSynthesizer {
    custom_callback: fn(&[usize], &[usize], &CliffordTableau) -> (usize, usize),
}

impl Default for CustomCallbackCliffordSynthesizer {
    fn default() -> Self {
        Self {
            custom_callback: |c: &[usize], r: &[usize], _ct: &CliffordTableau| (c[0], r[0]),
        }
    }
}

impl CustomCallbackCliffordSynthesizer {
    pub fn set_custom_callback(
        &mut self,
        callback: fn(&[usize], &[usize], &CliffordTableau) -> (usize, usize),
    ) -> &mut Self {
        self.custom_callback = callback;
        self
    }
}

impl<G> HasAdjoint<CliffordTableau, G> for CustomCallbackCliffordSynthesizer
where
    G: CliffordGates,
{
    fn synthesize_adjoint(&mut self, ct: CliffordTableau, repr: &mut G) {
        let mut ct = ct.adjoint();
        let num_qubits = ct.size();

        let mut remaining_columns = (0..num_qubits).collect::<Vec<_>>();
        let mut remaining_rows = (0..num_qubits).collect::<Vec<_>>();

        let mut custom_columns = vec![];
        let mut custom_rows = vec![];

        // for (&pivot_column, &pivot_row) in zip(custom_columns, custom_rows) {
        while !remaining_columns.is_empty() {
            let (pivot_column, pivot_row) =
                (self.custom_callback)(&remaining_columns, &remaining_rows, &ct);
            // Cleanup pivot column
            custom_columns.push(pivot_column);
            custom_rows.push(pivot_row);

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
