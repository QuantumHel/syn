use crate::{
    data_structures::{CliffordTableau, PropagateClifford},
    ir::{AdjointSynthesizer, CliffordGates},
};

use super::helper::{
    clean_pivot, clean_signs, clean_x_observables, clean_z_observables, naive_pivot_search,
};

use crate::data_structures::PauliLetter;

/// Reduces a tableau qubit by qubit, ignoring hardware connectivity.
///
/// For each row in turn it finds a pivot column, swaps it into place, clears the pivot to a bare
/// `X` and then to a bare `Z`, and uses it to eliminate that row's remaining entries in the
/// destabilizer and stabilizer halves. A final pass fixes the signs.
///
/// Because it emits `CX` between arbitrary qubit pairs, the result needs routing before it can run
/// on hardware with limited connectivity — use
/// [`PermRowColCliffordSynthesizer`](super::PermRowColCliffordSynthesizer) to account for
/// connectivity during synthesis instead.
///
/// Returns the residual permutation tableau.
#[derive(Default, Debug)]
pub struct NaiveCliffordSynthesizer {}

impl NaiveCliffordSynthesizer {
    /// The short name of this strategy, `"naive"`.
    pub fn name(&self) -> &str {
        return "naive";
    }
}

impl<G> AdjointSynthesizer<CliffordTableau, G, CliffordTableau> for NaiveCliffordSynthesizer
where
    G: CliffordGates,
{
    fn synthesize_adjoint(
        &mut self,
        mut clifford_tableau: CliffordTableau,
        repr: &mut G,
    ) -> CliffordTableau {
        let num_qubits = clifford_tableau.size();

        for row in 0..num_qubits {
            let pivot_col = naive_pivot_search(&clifford_tableau, num_qubits, row);

            if pivot_col != row {
                repr.swap(pivot_col, row);
                clifford_tableau.swap(row, pivot_col);
            }

            // Cleanup pivot column
            // clean_naive_pivot(repr, &mut clifford_tableau, row, row);
            clean_pivot(repr, &mut clifford_tableau, row, row, PauliLetter::X);

            let checked_rows = (row + 1..num_qubits).collect::<Vec<_>>();

            // Use the pivot to remove all other terms in the X observable.
            clean_x_observables(repr, &mut clifford_tableau, &checked_rows, row, row);

            clean_pivot(repr, &mut clifford_tableau, row, row, PauliLetter::Z);
            // Use the pivot to remove all other terms in the Z observable.
            clean_z_observables(repr, &mut clifford_tableau, &checked_rows, row, row);
        }

        clean_signs(repr, &mut clifford_tableau);
        clifford_tableau
    }
}
