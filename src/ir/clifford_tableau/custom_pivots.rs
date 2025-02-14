use std::iter::zip;

use itertools::Itertools;

use crate::{
    data_structures::CliffordTableau,
    ir::{
        clifford_tableau::helper::{
            clean_signs, clean_x_observables, clean_x_pivot, clean_z_observables, clean_z_pivot,
        },
        CliffordGates, Synthesizer,
    },
};

#[derive(Default)]
pub struct CustomPivotCliffordSynthesizer {
    custom_rows: Vec<usize>,
    custom_columns: Vec<usize>,
}

impl CustomPivotCliffordSynthesizer {
    pub fn set_custom_columns(&mut self, custom_columns: Vec<usize>) -> &mut Self {
        self.custom_columns = custom_columns;
        self
    }

    pub fn set_custom_rows(&mut self, custom_rows: Vec<usize>) -> &mut Self {
        self.custom_rows = custom_rows;
        self
    }
}

impl<G> Synthesizer<CliffordTableau, G, ()> for CustomPivotCliffordSynthesizer
where
    G: CliffordGates,
{
    fn synthesize(&mut self, mut clifford_tableau: CliffordTableau, repr: &mut G) {
        clifford_tableau = clifford_tableau.adjoint();
        self.synthesize_adjoint(clifford_tableau, repr);
    }

    fn synthesize_adjoint(&mut self, mut clifford_tableau: CliffordTableau, repr: &mut G) {
        let num_qubits = clifford_tableau.size();

        let mut remaining_columns = (0..num_qubits).collect::<Vec<_>>();
        let mut remaining_rows = (0..num_qubits).collect::<Vec<_>>();

        let custom_columns = &self.custom_columns;
        let custom_rows = &self.custom_rows;

        assert_eq!(
            custom_columns.iter().copied().sorted().collect::<Vec<_>>(),
            remaining_columns, "custom_columns is not a valid permutation, use `set_custom_columns` to define custom column pivots"
        );

        assert_eq!(
            custom_rows.iter().copied().sorted().collect::<Vec<_>>(),
            remaining_rows, "custom_rows is not a valid permutation, use `set_custom_rows` to define custom row pivots"
        );

        for (&pivot_column, &pivot_row) in zip(custom_columns, custom_rows) {
            // Cleanup pivot column

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
        let final_permutation = zip(custom_columns.clone(), custom_rows.clone())
            .sorted_by_key(|a| a.1)
            .map(|a| a.0)
            .collect::<Vec<_>>();

        clean_signs(repr, &mut clifford_tableau, &final_permutation);
    }
}
