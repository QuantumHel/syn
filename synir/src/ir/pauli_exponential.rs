//! Whole-circuit synthesis: the polynomial half, then the Clifford half.
//!
//! See [`PauliExponentialSynthesizer`].

use core::num;
use std::collections::VecDeque;

use crate::{
    architecture::connectivity::Connectivity,
    data_structures::{CliffordTableau, HasAdjoint, PauliExponential, PauliPolynomial},
    ir::AdjointSynthesizer,
};

use crate::ir::pauli_polynomial::psgs::PSGSPauliPolynomialSynthesizer;
use crate::ir::{CliffordGates, Gates, Synthesizer};

use crate::ir::clifford_tableau::NaiveCliffordSynthesizer;
use crate::ir::clifford_tableau::{CallbackCliffordSynthesizer, PermRowColCliffordSynthesizer};
use crate::ir::{
    clifford_tableau::CliffordTableauSynthStrategy,
    pauli_polynomial::{naive::NaivePauliPolynomialSynthesizer, PauliPolynomialSynthStrategy},
};

/// Synthesizes a [`PauliExponential`] end to end, into a sink implementing
/// [`CliffordGates`] + [`Gates`].
///
/// The two halves of the normal form are handled by separate algorithms, chosen independently:
///
/// 1. the polynomials, via [`PauliPolynomialSynthStrategy`];
/// 2. the trailing Clifford, via [`CliffordTableauSynthStrategy`].
///
/// # Why the two halves are not independent
///
/// The polynomial synthesizer does not un-compute the basis changes it introduces; it returns them
/// as a tableau `ct`, so what it actually emitted is `C_pp` with `ct = C_pp * U_pp^dag`. The
/// Clifford half therefore has to clean up after it, synthesizing `T * ct^dag` rather than `T`
/// alone. Since the tableau algorithms emit inverses, the value handed to them is
/// `(T * ct^dag)^dag = ct * T^dag`, spelled `T^dag.compose(&ct)` in circuit order.
///
/// Folding the fix-up into a tableau that has to be synthesized anyway makes it nearly free —
/// which is the whole reason the polynomial synthesizers skip the un-computation.
///
/// Whatever permutation is left over at the end goes to
/// [`CliffordGates::add_final_permutation`].
///
/// # Panics
///
/// Panics if the residual tableau after Clifford synthesis is not a bare permutation, which would
/// mean one of the two halves failed to reduce.
#[derive(Default)]
pub struct PauliExponentialSynthesizer {
    pauli_strategy: PauliPolynomialSynthStrategy,
    clifford_strategy: CliffordTableauSynthStrategy,
}

impl PauliExponentialSynthesizer {
    /// Builds a synthesizer with both strategies chosen up front.
    pub fn from_strategy(
        pauli_strategy: PauliPolynomialSynthStrategy,
        clifford_strategy: CliffordTableauSynthStrategy,
    ) -> Self {
        Self {
            pauli_strategy,
            clifford_strategy,
        }
    }

    /// Chooses the algorithm for the polynomial half.
    pub fn set_pauli_strategy(
        &mut self,
        pauli_strategy: PauliPolynomialSynthStrategy,
    ) -> &mut Self {
        self.pauli_strategy = pauli_strategy;
        self
    }

    /// Chooses the algorithm for the Clifford half.
    pub fn set_clifford_strategy(
        &mut self,
        clifford_strategy: CliffordTableauSynthStrategy,
    ) -> &mut Self {
        self.clifford_strategy = clifford_strategy;
        self
    }
}

impl<G> Synthesizer<PauliExponential, G> for PauliExponentialSynthesizer
where
    G: CliffordGates + Gates,
{
    fn synthesize(&mut self, pauli_exponential: PauliExponential, repr: &mut G) {
        let PauliExponential {
            pauli_polynomials,
            clifford_tableau,
        } = pauli_exponential;
        let num_qubits = clifford_tableau.size();
        let ct = match pauli_polynomials.is_empty() {
            true => CliffordTableau::new(num_qubits), // Skip PauliPolynomial synthesis
            false => match self.pauli_strategy {
                PauliPolynomialSynthStrategy::Naive => {
                    let mut pauli_synthesizer = NaivePauliPolynomialSynthesizer::default();
                    pauli_synthesizer.synthesize(pauli_polynomials, repr)
                }
                PauliPolynomialSynthStrategy::PSGS => {
                    let mut pauli_synthesizer = PSGSPauliPolynomialSynthesizer::default();
                    pauli_synthesizer.set_connectivity(Connectivity::complete(num_qubits));
                    pauli_synthesizer.synthesize(pauli_polynomials, repr)
                }
            },
        };
        // The emitted circuit is [gadget circuit][Clifford circuit], so the whole thing is
        // `U_cliff * C_pp`, and the IR's normal form is `U = T * U_pp` with `T = clifford_tableau`.
        //
        // The Pauli-polynomial synthesizer skips the un-computation of its basis changes and
        // returns it as `ct`, so `ct = C_pp * U_pp^dag`, i.e. `U_pp = ct^dag * C_pp`. Substituting:
        //
        //     U = T * ct^dag * C_pp   =>   U_cliff = T * ct^dag
        //
        // `synthesize_adjoint(X)` emits `X^dag`, so it must be handed `(T * ct^dag)^dag = ct * T^dag`.
        // `a.compose(&b)` is the matrix `U_b * U_a` (pinned by
        // `convention_tests::compose_matches_matrix_product`), so that is `T^dag.compose(&ct)`.
        let combined_ct = clifford_tableau.adjoint().compose(&ct);
        let final_ct = match &self.clifford_strategy {
            CliffordTableauSynthStrategy::Naive => {
                let mut clifford_synthesizer = NaiveCliffordSynthesizer::default();
                clifford_synthesizer.synthesize_adjoint(combined_ct, repr)
            }
            CliffordTableauSynthStrategy::Custom(custom_rows, custom_columns) => {
                let mut clifford_synthesizer = CallbackCliffordSynthesizer::custom_pivot(
                    custom_columns.to_owned(),
                    custom_rows.to_owned(),
                );
                clifford_synthesizer.synthesize_adjoint(combined_ct, repr)
            }
            CliffordTableauSynthStrategy::PermRowCol => {
                let mut clifford_synthesizer =
                    PermRowColCliffordSynthesizer::new(Connectivity::complete(num_qubits));
                clifford_synthesizer.synthesize_adjoint(combined_ct, repr)
            }
        };
        let final_perm = final_ct.get_permutation();
        match final_perm {
            Some(perm) => repr.add_final_permutation(perm),
            None => panic!("Final state was not a permutation: {final_ct}"),
        }
    }
}

/// Prints each gadget and angle of every polynomial to stdout, one per line, with `--` between
/// polynomials.
///
/// A debugging aid, used by the test suite to make a failing comparison readable.
pub fn print_pp_help(pauli_polynomials: &VecDeque<PauliPolynomial>) {
    for pp in pauli_polynomials {
        for i in 0..pp.length() {
            println!("{}, {}", pp.pauli_gadget(i), pp.angle(i));
        }
        println!("--");
    }
}
