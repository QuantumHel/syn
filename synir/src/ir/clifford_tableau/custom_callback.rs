use crate::{
    data_structures::CliffordTableau,
    ir::{
        clifford_tableau::helper::{clean_signs, clean_x_pivot, clean_z_pivot},
        AdjointSynthesizer, CliffordGates,
    },
};

use super::helper::{clean_x_observables, clean_z_observables};

/// Chooses the next `(column, row)` pivot, given the columns and rows still outstanding and the
/// current tableau.
type CliffordCallBack = Box<dyn FnMut(&[usize], &[usize], &CliffordTableau) -> (usize, usize)>;

/// Tableau synthesis with the pivot choice delegated to a caller-supplied closure.
///
/// The elimination itself is the same as
/// [`NaiveCliffordSynthesizer`](super::NaiveCliffordSynthesizer); only the order is external. That
/// makes this the hook for experimenting with pivot heuristics, and — via
/// [`custom_pivot`](Self::custom_pivot) — for replaying a fixed schedule so a run is reproducible.
///
/// The closure is called once per pivot and must return a `(column, row)` drawn from the slices it
/// is given; returning anything else desynchronizes the elimination.
pub struct CallbackCliffordSynthesizer {
    custom_callback: CliffordCallBack,
}

impl CallbackCliffordSynthesizer {
    /// Builds a synthesizer driven by `custom_callback`.
    pub fn new(custom_callback: CliffordCallBack) -> Self {
        Self { custom_callback }
    }

    /// Replays a fixed pivot schedule, taking the `i`-th pivot from
    /// `(custom_columns[i], custom_rows[i])`.
    ///
    /// # Panics
    ///
    /// Panics during synthesis if either vector is shorter than the qubit count.
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
    /// Swaps in a different pivot callback.
    pub fn set_custom_callback(&mut self, callback: CliffordCallBack) -> &mut Self {
        self.custom_callback = callback;
        self
    }
}

impl<G> AdjointSynthesizer<CliffordTableau, G, CliffordTableau> for CallbackCliffordSynthesizer
where
    G: CliffordGates,
{
    fn synthesize_adjoint(
        &mut self,
        mut clifford_tableau: CliffordTableau,
        repr: &mut G,
    ) -> CliffordTableau {
        let num_qubits = clifford_tableau.size();

        let mut remaining_columns = (0..num_qubits).collect::<Vec<_>>();
        let mut remaining_rows = (0..num_qubits).collect::<Vec<_>>();

        // let mut custom_columns = vec![];
        // let mut custom_rows = vec![];

        while !remaining_columns.is_empty() {
            let (pivot_column, pivot_row) =
                (self.custom_callback)(&remaining_columns, &remaining_rows, &clifford_tableau);

            remaining_columns.retain(|&x| x != pivot_column);
            remaining_rows.retain(|&x| x != pivot_row);
            {
                clean_x_pivot(repr, &mut clifford_tableau, pivot_column, pivot_row);
                // Use the pivot to remove all other terms in the X observable.
                clean_x_observables(
                    repr,
                    &mut clifford_tableau,
                    &remaining_columns,
                    pivot_column,
                    pivot_row,
                );

                clean_z_pivot(repr, &mut clifford_tableau, pivot_column, pivot_row);
                // Use the pivot to remove all other terms in the Z observable.
                clean_z_observables(
                    repr,
                    &mut clifford_tableau,
                    &remaining_columns,
                    pivot_column,
                    pivot_row,
                );
            }
        }
        clean_signs(repr, &mut clifford_tableau);
        clifford_tableau
    }
}
