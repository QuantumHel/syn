use crate::synthesis::{CommandCollector, PyCommand};

use bitvec::prelude::BitVec;
use pyo3::basic::CompareOp;
use pyo3::prelude::*;
use pyo3::PyResult;
use syn::data_structures::clifford_tableau::CliffordTableau as SynCliffordTableau;
use syn::data_structures::pauli_string::PauliString;
use syn::data_structures::PropagateClifford;

use syn::ir::clifford_tableau::CliffordTableauSynthesizer;
use syn::ir::CliffordGatesPrinter;

use syn::synthesis_methods::naive::Naive;

#[pyclass(unsendable)]
pub struct CliffordTableau {
    tableau: SynCliffordTableau,
}

#[pymethods]
impl CliffordTableau {
    #[new]
    fn new(n: usize) -> Self {
        CliffordTableau {
            tableau: SynCliffordTableau::new(n),
        }
    }

    #[staticmethod]
    pub fn from_parts(pauli_strings: Vec<String>, signs: Vec<bool>) -> Self {
        let pauli_columns: Vec<PauliString> = pauli_strings
            .iter()
            .map(|pauli_string| PauliString::from_text(pauli_string))
            .collect();
        let signs_bitvec: BitVec = signs.iter().copied().collect();
        let tableau = SynCliffordTableau::from_parts(pauli_columns, signs_bitvec);

        CliffordTableau { tableau }
    }

    pub fn size(&self) -> usize {
        self.tableau.size()
    }

    pub(crate) fn compose(&self, rhs: &Self) -> Self {
        CliffordTableau {
            tableau: self.tableau.compose(&rhs.tableau),
        }
    }

    fn __richcmp__(&self, other: &CliffordTableau, op: CompareOp) -> PyResult<bool> {
        match op {
            CompareOp::Eq => Ok(self.tableau == other.tableau),
            CompareOp::Ne => Ok(self.tableau != other.tableau),
            _ => Ok(false),
        }
    }

    // pub fn synthesize(&self) -> Vec<PyCommand> {
    //     let mut tracker = CommandCollector::new();
    //     CliffordTableauSynthesizer::run(&self.tableau, &mut tracker);
    //     tracker.commands()
    // }
    pub fn synthesize(&self) -> Vec<String> {
        let mut tracker = CliffordGatesPrinter::new(self.tableau.size());
        CliffordTableauSynthesizer::run(&self.tableau, &mut tracker);
        tracker.gates
    }
}

impl PropagateClifford for CliffordTableau {
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