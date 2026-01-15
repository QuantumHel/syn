extern crate pyo3;
extern crate pyo3_ffi;

use pyo3::prelude::*;
use synir::ir::{CliffordGates, Gates};

#[pyclass]
pub struct QiskitSynIR {
    circuit: Py<PyAny>,
}

#[pymethods]
impl QiskitSynIR {
    #[new]
    pub fn new(qiskit_circuit: Py<PyAny>) -> Self {
        QiskitSynIR {
            circuit: qiskit_circuit,
        }
    }

    pub fn s(&mut self, target: synir::IndexType) {
        Python::attach(|py| -> PyResult<()> {
            self.circuit.call_method1(py, "s", (target,))?;
            Ok(())
        })
        .unwrap();
    }

    pub fn v(&mut self, target: synir::IndexType) {
        Python::attach(|py| -> PyResult<()> {
            self.circuit.call_method1(py, "sx", (target,))?;
            Ok(())
        })
        .unwrap();
    }

    pub fn s_dgr(&mut self, target: synir::IndexType) {
        Python::attach(|py| -> PyResult<()> {
            self.circuit.call_method1(py, "sdg", (target,))?;
            Ok(())
        })
        .unwrap();
    }

    pub fn v_dgr(&mut self, target: synir::IndexType) {
        Python::attach(|py| -> PyResult<()> {
            self.circuit.call_method1(py, "sxdg", (target,))?;
            Ok(())
        })
        .unwrap();
    }

    pub fn x(&mut self, target: synir::IndexType) {
        Python::attach(|py| -> PyResult<()> {
            self.circuit.call_method1(py, "x", (target,))?;
            Ok(())
        })
        .unwrap();
    }

    pub fn y(&mut self, target: synir::IndexType) {
        Python::attach(|py| -> PyResult<()> {
            self.circuit.call_method1(py, "y", (target,))?;
            Ok(())
        })
        .unwrap();
    }

    pub fn z(&mut self, target: synir::IndexType) {
        Python::attach(|py| -> PyResult<()> {
            self.circuit.call_method1(py, "z", (target,))?;
            Ok(())
        })
        .unwrap();
    }

    pub fn h(&mut self, target: synir::IndexType) {
        Python::attach(|py| -> PyResult<()> {
            self.circuit.call_method1(py, "h", (target,))?;
            Ok(())
        })
        .unwrap();
    }

    pub fn cx(&mut self, control: synir::IndexType, target: synir::IndexType) {
        Python::attach(|py| -> PyResult<()> {
            self.circuit.call_method1(py, "cx", (control, target))?;
            Ok(())
        })
        .unwrap();
    }

    pub fn cz(&mut self, control: synir::IndexType, target: synir::IndexType) {
        Python::attach(|py| -> PyResult<()> {
            self.circuit.call_method1(py, "cz", (control, target))?;
            Ok(())
        })
        .unwrap();
    }

    pub fn rx(&mut self, target: synir::IndexType, angle: f64) {
        Python::attach(|py| -> PyResult<()> {
            self.circuit.call_method1(py, "rx", (angle, target))?;
            Ok(())
        })
        .unwrap();
    }

    pub fn ry(&mut self, target: synir::IndexType, angle: f64) {
        Python::attach(|py| -> PyResult<()> {
            self.circuit.call_method1(py, "ry", (angle, target))?;
            Ok(())
        })
        .unwrap();
    }

    pub fn rz(&mut self, target: synir::IndexType, angle: f64) {
        Python::attach(|py| -> PyResult<()> {
            self.circuit.call_method1(py, "rz", (angle, target))?;
            Ok(())
        })
        .unwrap();
    }
}

impl CliffordGates for QiskitSynIR {
    fn s(&mut self, target: synir::IndexType) {
        self.s(target);
    }

    fn v(&mut self, target: synir::IndexType) {
        self.v(target);
    }

    fn s_dgr(&mut self, target: synir::IndexType) {
        self.s_dgr(target);
    }

    fn v_dgr(&mut self, target: synir::IndexType) {
        self.v_dgr(target);
    }

    fn x(&mut self, target: synir::IndexType) {
        self.x(target);
    }

    fn y(&mut self, target: synir::IndexType) {
        self.y(target);
    }

    fn z(&mut self, target: synir::IndexType) {
        self.z(target);
    }

    fn h(&mut self, target: synir::IndexType) {
        self.h(target);
    }

    fn cx(&mut self, control: synir::IndexType, target: synir::IndexType) {
        self.cx(control, target);
    }

    fn cz(&mut self, control: synir::IndexType, target: synir::IndexType) {
        self.cz(control, target);
    }
}

impl Gates for QiskitSynIR {
    fn rx(&mut self, target: synir::IndexType, angle: f64) {
        self.rx(target, angle);
    }

    fn ry(&mut self, target: synir::IndexType, angle: f64) {
        self.ry(target, angle);
    }

    fn rz(&mut self, target: synir::IndexType, angle: f64) {
        self.rz(target, angle);
    }
}
