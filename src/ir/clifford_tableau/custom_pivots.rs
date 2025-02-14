use std::iter::zip;

use itertools::Itertools;

use crate::{
    data_structures::CliffordTableau,
    ir::{
        clifford_tableau::helper::{
            clean_signs, clean_x_observables, clean_x_pivot, clean_z_observables, clean_z_pivot,
        },
        CliffordGates,
    },
};

use super::CliffordTableauSynthesizer;

pub struct CustomPivotCliffordSynthesizer {
    clifford_tableau: CliffordTableau,
    custom_rows: Vec<usize>,
    custom_columns: Vec<usize>,
}

impl CustomPivotCliffordSynthesizer {
    pub fn new(
        clifford_tableau: CliffordTableau,
        custom_rows: Vec<usize>,
        custom_columns: Vec<usize>,
    ) -> Self {
        Self {
            clifford_tableau,
            custom_rows,
            custom_columns,
        }
    }
}

impl<G> CliffordTableauSynthesizer<G> for CustomPivotCliffordSynthesizer
where
    G: CliffordGates,
{
    fn synthesize(&mut self, repr: &mut G) {
        self.clifford_tableau = self.clifford_tableau.adjoint();
        CliffordTableauSynthesizer::<G>::synthesize_adjoint(self, repr)
    }

    fn synthesize_adjoint(&mut self, repr: &mut G) {
        let num_qubits = self.clifford_tableau.size();
        let ct = &mut self.clifford_tableau;

        let mut remaining_columns = (0..num_qubits).collect::<Vec<_>>();
        let mut remaining_rows = (0..num_qubits).collect::<Vec<_>>();

        let custom_columns = &self.custom_columns;
        let custom_rows = &self.custom_rows;

        assert_eq!(
            custom_columns.iter().copied().sorted().collect::<Vec<_>>(),
            remaining_columns
        );

        assert_eq!(
            custom_rows.iter().copied().sorted().collect::<Vec<_>>(),
            remaining_rows
        );

        for (&pivot_column, &pivot_row) in zip(custom_columns, custom_rows) {
            // Cleanup pivot column

            remaining_columns.retain(|&x| x != pivot_column);
            remaining_rows.retain(|&x| x != pivot_row);
            {
                clean_x_pivot(repr, ct, pivot_column, pivot_row);

                // Use the pivot to remove all other terms in the X observable.
                clean_x_observables(repr, ct, &remaining_rows, pivot_column, pivot_row);

                clean_z_pivot(repr, ct, pivot_column, pivot_row);

                // Use the pivot to remove all other terms in the Z observable.
                clean_z_observables(repr, ct, &remaining_rows, pivot_column, pivot_row);
            }
        }
        let final_permutation = zip(custom_columns.clone(), custom_rows.clone())
            .sorted_by_key(|a| a.1)
            .map(|a| a.0)
            .collect::<Vec<_>>();

        clean_signs(repr, ct, &final_permutation);
    }
}

// impl<G> Custom<&CliffordTableau, G> for CliffordTableauSynthesizer
// where
//     G: CliffordGates,
// {
//     fn run_custom(
//         clifford_tableau: &CliffordTableau,
//         repr: &mut G,
//         custom_columns: Vec<usize>,
//         custom_rows: Vec<usize>,
//     ) {
//         let mut ct = clifford_tableau.adjoint();
//         let num_qubits = ct.size();
//         let mut remaining_columns = (0..num_qubits).collect::<Vec<_>>();
//         let mut remaining_rows = (0..num_qubits).collect::<Vec<_>>();

//         assert_eq!(
//             custom_columns.iter().copied().sorted().collect::<Vec<_>>(),
//             remaining_columns
//         );

//         assert_eq!(
//             custom_rows.iter().copied().sorted().collect::<Vec<_>>(),
//             remaining_rows
//         );

//         for (&pivot_column, &pivot_row) in zip(&custom_columns, &custom_rows) {
//             // Cleanup pivot column

//             remaining_columns.retain(|&x| x != pivot_column);
//             remaining_rows.retain(|&x| x != pivot_row);
//             {
//                 clean_x_pivot(repr, &mut ct, pivot_column, pivot_row);

//                 // Use the pivot to remove all other terms in the X observable.
//                 clean_x_observables(repr, &mut ct, &remaining_rows, pivot_column, pivot_row);

//                 clean_z_pivot(repr, &mut ct, pivot_column, pivot_row);

//                 // Use the pivot to remove all other terms in the Z observable.
//                 clean_z_observables(repr, &mut ct, &remaining_rows, pivot_column, pivot_row);
//             }
//         }
//         let final_permutation = zip(custom_columns.clone(), custom_rows.clone())
//             .sorted_by_key(|a| a.1)
//             .map(|a| a.0)
//             .collect::<Vec<_>>();

//         clean_signs(repr, &mut ct, &final_permutation);
//     }
// }
