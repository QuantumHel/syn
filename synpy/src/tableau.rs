//! The Python-facing Clifford tableau.

use crate::synthesis::{CommandCollector, PyCommand, Synthesize};

use bitvec::prelude::BitVec;
use pyo3::basic::CompareOp;
use pyo3::prelude::*;
use pyo3::PyResult;
use synir::data_structures::CliffordTableau;
use synir::data_structures::PauliString;
use synir::data_structures::PropagateClifford;
use synir::ir::clifford_tableau::NaiveCliffordSynthesizer;
use synir::ir::Synthesizer;

/// A Clifford operator on `n` qubits, exposed to Python.
///
/// Wraps [`synir::data_structures::CliffordTableau`]. Build one from Qiskit's `Clifford` with
/// [`from_parts`](Self::from_parts), or start from the identity and propagate gates into it.
///
/// `synthesize` always uses the naive strategy, which ignores connectivity — go through
/// [`PyPauliExponential`](crate::wrapper::PyPauliExponential) if you need a
/// connectivity-aware result.
///
/// Supports `==` and `!=`.
#[pyclass(unsendable)]
pub struct PyCliffordTableau {
    tableau: CliffordTableau,
}

#[pymethods]
impl PyCliffordTableau {
    /// The identity operation on `n` qubits.
    #[new]
    fn new(n: usize) -> Self {
        PyCliffordTableau {
            tableau: CliffordTableau::new(n),
        }
    }

    /// Builds a tableau from its columns and signs.
    ///
    /// `pauli_strings` holds one string per qubit, each `2n` letters long — rows `0..n` are the
    /// images of the `X` generators and rows `n..2n` those of the `Z` generators. `signs` has `2n`
    /// booleans in the same order, `True` meaning negative.
    ///
    /// Note this is the *transpose* of how Qiskit's `Clifford.tableau` is laid out; see
    /// `synpy.qiskit.plugin.SynPyCliffordPlugin` for the conversion.
    #[staticmethod]
    pub fn from_parts(pauli_strings: Vec<String>, signs: Vec<bool>) -> Self {
        let pauli_columns: Vec<PauliString> = pauli_strings
            .iter()
            .map(|pauli_string| PauliString::from_text(pauli_string))
            .collect();
        let signs_bitvec: BitVec = signs.iter().copied().collect();
        let tableau = CliffordTableau::from_parts(pauli_columns, signs_bitvec);

        PyCliffordTableau { tableau }
    }

    /// The number of qubits.
    pub fn size(&self) -> usize {
        self.tableau.size()
    }

    /// Composes in circuit order: `self` runs first, then `rhs`.
    pub(crate) fn compose(&self, rhs: &Self) -> Self {
        PyCliffordTableau {
            tableau: self.tableau.compose(&rhs.tableau),
        }
    }

    fn __richcmp__(&self, other: &PyCliffordTableau, op: CompareOp) -> PyResult<bool> {
        match op {
            CompareOp::Eq => Ok(self.tableau == other.tableau),
            CompareOp::Ne => Ok(self.tableau != other.tableau),
            _ => Ok(false),
        }
    }

    /// Synthesizes this tableau into a list of [`PyCommand`]s, in circuit order.
    ///
    /// Uses the naive strategy, so the result assumes all-to-all connectivity.
    pub fn synthesize(&self) -> Vec<PyCommand> {
        <Self as Synthesize>::synthesize(self)
    }
}

impl PropagateClifford for PyCliffordTableau {
    fn cx(&mut self, control: usize, target: usize) -> &mut Self {
        self.tableau.cx(control, target);
        self
    }

    fn s(&mut self, target: usize) -> &mut Self {
        self.tableau.s(target);
        self
    }

    fn v(&mut self, target: usize) -> &mut Self {
        self.tableau.v(target);
        self
    }
}

impl Synthesize for PyCliffordTableau {
    fn synthesize(&self) -> Vec<PyCommand> {
        let mut tracker = CommandCollector::new();
        let mut synthesizer = NaiveCliffordSynthesizer::default();
        synthesizer.synthesize(self.tableau.clone(), &mut tracker);
        tracker.commands()
    }
}
