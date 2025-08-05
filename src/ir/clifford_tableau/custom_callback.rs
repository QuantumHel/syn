use std::iter::zip;

use itertools::Itertools;

use crate::{
    data_structures::CliffordTableau,
    ir::{
        clifford_tableau::helper::{clean_signs, clean_x_pivot, clean_z_pivot},
        AdjointSynthesizer, CliffordGates,
    },
};

use super::helper::{clean_x_observables, clean_z_observables};

pub struct CallbackCliffordSynthesizer {
    custom_callback: Box<dyn FnMut(&[usize], &[usize], &CliffordTableau) -> (usize, usize)>,
}

impl CallbackCliffordSynthesizer {
    pub fn new(
        custom_callback: Box<dyn FnMut(&[usize], &[usize], &CliffordTableau) -> (usize, usize)>,
    ) -> Self {
        Self { custom_callback }
    }

    pub fn custom_pivot(custom_columns: Vec<usize>, custom_rows: Vec<usize>) -> Self {
        let mut loc = 0;
        Self {
            custom_callback: Box::new(
                move |_cc: &[usize], _cr: &[usize], _ct: &CliffordTableau| {
                    let next = (custom_columns[loc], custom_rows[loc]);
                    loc += 1;
                    next
                },
            ),
        }
    }
}

impl Default for CallbackCliffordSynthesizer {
    fn default() -> Self {
        Self::new(Box::new(
            |c: &[usize], r: &[usize], _ct: &CliffordTableau| (c[0], r[0]),
        ))
    }
}

impl CallbackCliffordSynthesizer {
    pub fn set_custom_callback(
        &mut self,
        callback: Box<dyn FnMut(&[usize], &[usize], &CliffordTableau) -> (usize, usize)>,
    ) -> &mut Self {
        self.custom_callback = callback;
        self
    }
}

impl<'a, G> AdjointSynthesizer<CliffordTableau, G> for CallbackCliffordSynthesizer
where
    G: CliffordGates,
{
    fn synthesize_adjoint(&mut self, mut clifford_tableau: CliffordTableau, repr: &mut G) {
        let num_qubits = clifford_tableau.size();

        let mut remaining_columns = (0..num_qubits).collect::<Vec<_>>();
        let mut remaining_rows = (0..num_qubits).collect::<Vec<_>>();

        let mut custom_columns = vec![];
        let mut custom_rows = vec![];

        // for (&pivot_column, &pivot_row) in zip(custom_columns, custom_rows) {
        while !remaining_columns.is_empty() {
            let (pivot_column, pivot_row) =
                (self.custom_callback)(&remaining_columns, &remaining_rows, &clifford_tableau);
            // Cleanup pivot column
            custom_columns.push(pivot_column);
            custom_rows.push(pivot_row);

            remaining_columns.retain(|&x| x != pivot_column);
            remaining_rows.retain(|&x| x != pivot_row);
            {
                clean_x_pivot(repr, &mut clifford_tableau, pivot_column, pivot_row);
                // Use the pivot to remove all other terms in the X observable.
                clean_x_observables(
                    repr,
                    &mut clifford_tableau,
                    &remaining_rows,
                    pivot_column,
                    pivot_row,
                );

                clean_z_pivot(repr, &mut clifford_tableau, pivot_column, pivot_row);
                // Use the pivot to remove all other terms in the Z observable.
                clean_z_observables(
                    repr,
                    &mut clifford_tableau,
                    &remaining_rows,
                    pivot_column,
                    pivot_row,
                );
            }
        }
        let final_permutation = zip(custom_columns, custom_rows)
            .sorted_by_key(|a| a.1)
            .map(|a| a.0)
            .collect::<Vec<_>>();

        clean_signs(repr, &mut clifford_tableau, &final_permutation);
    }
}
