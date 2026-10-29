//! A gate sink that writes straight into a Qiskit `QuantumCircuit`.

extern crate pyo3;
extern crate pyo3_ffi;

use pyo3::{prelude::*, types::PyList};
use synir::ir::{CliffordGates, Gates};

/// A [`CliffordGates`] + [`Gates`] sink backed by a Qiskit `QuantumCircuit`.
///
/// Each emitted gate is forwarded as a method call on the wrapped circuit object — `s`, `sx`,
/// `sdg`, `sxdg`, `x`, `y`, `z`, `h`, `cx`, `cz`, `rx`, `ry`, `rz` — so nothing is buffered and
/// the circuit is complete as soon as synthesis returns.
///
/// # The final permutation is stored, not emitted
///
/// Synthesis may finish with the qubits relabelled. The default
/// `CliffordGates::add_final_permutation` would expand that into CX-triple swaps; this type
/// instead **records** it, retrievable via [`get_permutation`](Self::get_permutation), and emits
/// no gates. A caller that ignores the permutation therefore gets a circuit that is correct only
/// up to that relabelling.
///
/// ```python
/// from qiskit import QuantumCircuit
/// from synpy.synpy_rust import QiskitSynIR
///
/// sink = QiskitSynIR(QuantumCircuit(n))
/// pe.synthesize_to_qiskit(sink)
/// circuit = sink.get_circuit()
/// perm = sink.get_permutation()      # None, or a list with perm[physical] = logical
/// ```
///
/// # Panics
///
/// Every gate method unwraps the Python call, so a circuit object missing one of those methods
/// aborts rather than raising.
#[pyclass]
pub struct QiskitSynIR {
    circuit: Py<PyAny>,
    final_permutation: Option<Py<PyList>>,
}

#[pymethods]
impl QiskitSynIR {
    /// Wraps a Qiskit `QuantumCircuit`, which must already have enough qubits.
    #[new]
    pub fn new(qiskit_circuit: Py<PyAny>) -> Self {
        QiskitSynIR {
            circuit: qiskit_circuit,
            final_permutation: None,
        }
    }

    /// The wrapped circuit, with everything emitted so far.
    fn get_circuit(&self, py: Python) -> Py<PyAny> {
        self.circuit.clone_ref(py)
    }

    /// The final qubit relabelling, if synthesis produced one.
    ///
    /// `perm[physical] = logical`. `None` means no relabelling was needed. This is *not* applied
    /// to the circuit — see the note on the type.
    pub fn get_permutation(&self) -> Option<&Py<PyList>> {
        match &self.final_permutation {
            Some(perm) => Some(perm),
            None => None,
        }
    }

    /// Appends a phase gate `S`.
    pub fn s(&mut self, target: synir::IndexType) {
        CliffordGates::s(self, target)
    }

    /// Appends `V = sqrt(X)`, as Qiskit's `sx`.
    pub fn v(&mut self, target: synir::IndexType) {
        CliffordGates::v(self, target)
    }

    /// Appends `S^dag`, as Qiskit's `sdg`.
    pub fn s_dgr(&mut self, target: synir::IndexType) {
        CliffordGates::s_dgr(self, target)
    }

    /// Appends `V^dag`, as Qiskit's `sxdg`.
    pub fn v_dgr(&mut self, target: synir::IndexType) {
        CliffordGates::v_dgr(self, target)
    }

    /// Appends Pauli `X`.
    pub fn x(&mut self, target: synir::IndexType) {
        CliffordGates::x(self, target)
    }

    /// Appends Pauli `Y`.
    pub fn y(&mut self, target: synir::IndexType) {
        CliffordGates::y(self, target)
    }

    /// Appends Pauli `Z`.
    pub fn z(&mut self, target: synir::IndexType) {
        CliffordGates::z(self, target)
    }

    /// Appends a Hadamard.
    pub fn h(&mut self, target: synir::IndexType) {
        CliffordGates::h(self, target)
    }

    /// Appends a `CX`.
    fn cx(&mut self, control: synir::IndexType, target: synir::IndexType) {
        CliffordGates::cx(self, control, target)
    }

    /// Appends a `CZ`.
    fn cz(&mut self, control: synir::IndexType, target: synir::IndexType) {
        CliffordGates::cz(self, control, target)
    }

    /// Records the final qubit relabelling without emitting any gates.
    fn add_final_permutation(&mut self, permutation: Vec<synir::IndexType>) {
        CliffordGates::add_final_permutation(self, permutation);
    }

    /// Appends an X rotation by `angle` radians.
    fn rx(&mut self, target: synir::IndexType, angle: f64) {
        Gates::rx(self, target, angle)
    }

    /// Appends a Y rotation by `angle` radians.
    fn ry(&mut self, target: synir::IndexType, angle: f64) {
        Gates::ry(self, target, angle)
    }

    /// Appends a Z rotation by `angle` radians.
    fn rz(&mut self, target: synir::IndexType, angle: f64) {
        Gates::rz(self, target, angle)
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
        Python::attach(|py| -> PyResult<()> {
            self.circuit.call_method1(py, "rz", (angle, target))?;
            Ok(())
        })
        .unwrap();
    }
}
