use bitvec::prelude::BitVec;
use pyo3::prelude::*;
use pyo3::basic::CompareOp;
use syn::data_structures::clifford_tableau::CliffordTableau as SynCliffordTableau;
use syn::data_structures::pauli_string::PauliString;

use syn::ir::clifford_tableau::CliffordTableauSynthesizer;
use syn::ir::CliffordGatesPrinter;

use syn::synthesis_methods::naive::Naive;

#[pyclass]
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
    pub fn from_parts(pauli_strings: Vec<String>, signs: Vec<bool>, n: usize) -> Self {
        let pauli_columns: Vec<PauliString> = pauli_strings
            .iter()
            .map(|pauli_string| PauliString::from_text(pauli_string))
            .collect();
        let signs_bitvec: BitVec = signs.iter().copied().collect();
        let tableau = SynCliffordTableau::from_parts(pauli_columns, signs_bitvec, n);

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

    fn __richcmp__(&self, other: &PyAny, op: CompareOp) -> PyResult<bool> {
        if let Ok(other) = other.extract::<PyRef<CliffordTableau>>() {
            match op {
                CompareOp::Eq => Ok(self.tableau == other.tableau),
                CompareOp::Ne => Ok(self.tableau != other.tableau),
                _ => Ok(false),
            }
        } else {
            Ok(false)
        }
    }

    pub fn synthesize(&self) -> Vec<String> {
        let mut printer = CliffordGatesPrinter::new(self.size());
        CliffordTableauSynthesizer::run(&self.tableau, &mut printer);
        printer.gates
    }
}

#[pyfunction]
fn add(x: i32, y: i32) -> PyResult<i32> {
    Ok(x + y)
}

#[pymodule]
fn synpy(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(add, m)?)?;
    m.add_class::<CliffordTableau>()?;
    Ok(())
}