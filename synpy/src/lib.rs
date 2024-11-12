use bitvec::prelude::BitVec;
use pyo3::prelude::*;
use pyo3::basic::CompareOp;
use syn::datastructures::clifford_tableau::CliffordTableau as SynCliffordTableau;
use syn::datastructures::pauli_string::PauliString;


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
    pub fn with_pauli_columns(n: usize, pauli_strings: Vec<String>, signs: Vec<bool>) -> Self {
        let pauli_columns: Vec<PauliString> = pauli_strings
            .iter()
            .map(|pauli_string| PauliString::from_text(pauli_string))
            .collect();
        let signs_bitvec: BitVec = signs.iter().copied().collect();
        let tableau = SynCliffordTableau::with_pauli_columns(n, pauli_columns, signs_bitvec);

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