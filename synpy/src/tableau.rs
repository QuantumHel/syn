use crate::synthesis::{CommandCollector, PyCommand, Synthesize};

use bitvec::prelude::BitVec;
use pyo3::basic::CompareOp;
use pyo3::prelude::*;
use pyo3::PyResult;
use syn::data_structures::CliffordTableau as CliffordTableau;
use syn::data_structures::PauliString;
use syn::data_structures::PropagateClifford;

use syn::ir::clifford_tableau::NaiveCliffordSynthesizer;
use syn::ir::Synthesizer;

#[pyclass(unsendable)]
pub struct PyCliffordTableau {
    tableau: CliffordTableau,
}

#[pymethods]
impl PyCliffordTableau {
    #[new]
    fn new(n: usize) -> Self {
        PyCliffordTableau {
            tableau: CliffordTableau::new(n),
        }
    }
    fn __str__(&self) -> PyResult<String> {
        Ok(self.tableau.to_string())
    }

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

    pub fn size(&self) -> usize {
        self.tableau.size()
    }

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