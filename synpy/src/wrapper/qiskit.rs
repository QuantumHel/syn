extern crate pyo3;
extern crate pyo3_ffi;

use pyo3::{prelude::*, types::PyList};
use synir::ir::{CliffordGates, Gates};

#[pyclass]
pub struct QiskitSynIR {
    circuit: Py<PyAny>,
    final_permutation: Option<Py<PyList>>,
}

#[pymethods]
impl QiskitSynIR {
    #[new]
    pub fn new(qiskit_circuit: Py<PyAny>) -> Self {
        QiskitSynIR {
            circuit: qiskit_circuit,
            final_permutation: None,
        }
    }

    fn get_circuit(&self, py: Python) -> Py<PyAny> {
        self.circuit.clone_ref(py)
    }

    pub fn get_permutation(&self) -> Option<&Py<PyList>> {
        match &self.final_permutation {
            Some(perm) => Some(perm),
            None => None,
        }
    }
}

impl CliffordGates for QiskitSynIR {
    fn s(&mut self, target: synir::IndexType) {
        Python::attach(|py| -> PyResult<()> {
            self.circuit.call_method1(py, "s", (target,))?;
            Ok(())
        })
        .unwrap();
    }

    fn v(&mut self, target: synir::IndexType) {
        Python::attach(|py| -> PyResult<()> {
            self.circuit.call_method1(py, "sx", (target,))?;
            Ok(())
        })
        .unwrap();
    }

    fn s_dgr(&mut self, target: synir::IndexType) {
        Python::attach(|py| -> PyResult<()> {
            self.circuit.call_method1(py, "sdg", (target,))?;
            Ok(())
        })
        .unwrap();
    }

    fn v_dgr(&mut self, target: synir::IndexType) {
        Python::attach(|py| -> PyResult<()> {
            self.circuit.call_method1(py, "sxdg", (target,))?;
            Ok(())
        })
        .unwrap();
    }

    fn x(&mut self, target: synir::IndexType) {
        Python::attach(|py| -> PyResult<()> {
            self.circuit.call_method1(py, "x", (target,))?;
            Ok(())
        })
        .unwrap();
    }

    fn y(&mut self, target: synir::IndexType) {
        Python::attach(|py| -> PyResult<()> {
            self.circuit.call_method1(py, "y", (target,))?;
            Ok(())
        })
        .unwrap();
    }

    fn z(&mut self, target: synir::IndexType) {
        Python::attach(|py| -> PyResult<()> {
            self.circuit.call_method1(py, "z", (target,))?;
            Ok(())
        })
        .unwrap();
    }

    fn h(&mut self, target: synir::IndexType) {
        Python::attach(|py| -> PyResult<()> {
            self.circuit.call_method1(py, "h", (target,))?;
            Ok(())
        })
        .unwrap();
    }

    fn cx(&mut self, control: synir::IndexType, target: synir::IndexType) {
        Python::attach(|py| -> PyResult<()> {
            self.circuit.call_method1(py, "cx", (control, target))?;
            Ok(())
        })
        .unwrap();
    }

    fn cz(&mut self, control: synir::IndexType, target: synir::IndexType) {
        Python::attach(|py| -> PyResult<()> {
            self.circuit.call_method1(py, "cz", (control, target))?;
            Ok(())
        })
        .unwrap();
    }

    fn add_final_permutation(&mut self, permutation: Vec<synir::IndexType>) {
        Python::attach(|py| -> () {
            match PyList::new(py, permutation) {
                Ok(list) => self.final_permutation = Some(list.unbind()),
                _ => (),
            }
        })
    }
}

impl Gates for QiskitSynIR {
    fn rx(&mut self, target: synir::IndexType, angle: f64) {
        Python::attach(|py| -> PyResult<()> {
            self.circuit.call_method1(py, "rx", (angle, target))?;
            Ok(())
        })
        .unwrap();
    }

    fn ry(&mut self, target: synir::IndexType, angle: f64) {
        Python::attach(|py| -> PyResult<()> {
            self.circuit.call_method1(py, "ry", (angle, target))?;
            Ok(())
        })
        .unwrap();
    }

    fn rz(&mut self, target: synir::IndexType, angle: f64) {
        println!("RZ on {}", target);
        Python::attach(|py| -> PyResult<()> {
            self.circuit.call_method1(py, "rz", (angle, target))?;
            Ok(())
        })
        .unwrap();
    }
}
