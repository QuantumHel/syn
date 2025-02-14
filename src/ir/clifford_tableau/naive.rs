use crate::{data_structures::CliffordTableau, ir::CliffordGates};

use super::{
    helper::{
        clean_pivot, clean_signs, clean_x_observables, clean_z_observables, naive_pivot_search,
        swap,
    },
    CliffordTableauSynthesizer,
};

pub struct NaiveCliffordSynthesizer {
    clifford_tableau: CliffordTableau,
}

impl NaiveCliffordSynthesizer {
    pub fn new(clifford_tableau: CliffordTableau) -> Self {
        Self { clifford_tableau }
    }
}

impl<G> CliffordTableauSynthesizer<G> for NaiveCliffordSynthesizer
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

        for row in 0..num_qubits {
            let pivot_col = naive_pivot_search(ct, num_qubits, row);

            if pivot_col != row {
                swap(repr, ct, row, pivot_col);
            }

            // Cleanup pivot column
            clean_pivot(repr, ct, row, row);

            let checked_rows = (row + 1..num_qubits).collect::<Vec<_>>();

            // Use the pivot to remove all other terms in the X observable.
            clean_x_observables(repr, ct, &checked_rows, row, row);

            // Use the pivot to remove all other terms in the Z observable.
            clean_z_observables(repr, ct, &checked_rows, row, row);
        }

        clean_signs(repr, ct, &(0..num_qubits).collect::<Vec<_>>());
    }
}

// impl<G> Naive<&CliffordTableau, G> for CliffordTableauSynthesizer
// where
//     G: CliffordGates,
// {
//     fn run_naive(clifford_tableau: &CliffordTableau, repr: &mut G) {
//         let mut ct = clifford_tableau.adjoint();

//         let num_qubits = ct.size();
//         for row in 0..num_qubits {
//             let pivot_col = naive_pivot_search(&ct, num_qubits, row);

//             if pivot_col != row {
//                 swap(repr, &mut ct, row, pivot_col);
//             }

//             // Cleanup pivot column
//             clean_pivot(repr, &mut ct, row, row);

//             let checked_rows = (row + 1..num_qubits).collect::<Vec<_>>();

//             // Use the pivot to remove all other terms in the X observable.
//             clean_x_observables(repr, &mut ct, &checked_rows, row, row);

//             // Use the pivot to remove all other terms in the Z observable.
//             clean_z_observables(repr, &mut ct, &checked_rows, row, row);
//         }

//         clean_signs(repr, &mut ct, &(0..num_qubits).collect::<Vec<_>>());
//     }
// }

// impl<G> NaiveAdjoint<&CliffordTableau, G> for CliffordTableauSynthesizer
// where
//     G: CliffordGates,
// {
//     fn run_naive_adjoint(clifford_tableau: &CliffordTableau, repr: &mut G) {
//         let mut ct = clifford_tableau.clone();
//         let num_qubits = ct.size();
//         for row in 0..num_qubits {
//             let pivot_col = naive_pivot_search(&ct, num_qubits, row);

//             if pivot_col != row {
//                 swap(repr, &mut ct, row, pivot_col);
//             }

//             // Cleanup pivot column
//             clean_pivot(repr, &mut ct, row, row);

//             let checked_rows = (row + 1..num_qubits).collect::<Vec<_>>();

//             // Use the pivot to remove all other terms in the X observable.
//             clean_x_observables(repr, &mut ct, &checked_rows, row, row);

//             // Use the pivot to remove all other terms in the Z observable.
//             clean_z_observables(repr, &mut ct, &checked_rows, row, row);
//         }

//         clean_signs(repr, &mut ct, &(0..num_qubits).collect::<Vec<_>>());
//     }
// }
