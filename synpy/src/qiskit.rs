extern crate pyo3;
extern crate pyo3_ffi;

use std::collections::VecDeque;

use pyo3::{prelude::*, types::{PyInt, PyList, PyString, PyTuple}};
use synir::{data_structures::CliffordTableau, ir::{CliffordGates, Gates, Synthesizer, clifford_tableau::CliffordTableauSynthStrategy, pauli_exponential::{PauliExponential, PauliExponentialSynthesizer}, pauli_polynomial::PauliPolynomialSynthStrategy}};

#[pyclass]
pub struct QiskitSynIR{
    circuit: Py<PyAny>
}

#[pymethods]
impl QiskitSynIR {

    #[new]
    pub fn new(qiskit_circuit: Py<PyAny>) -> Self {
        return QiskitSynIR{
            circuit: qiskit_circuit
        }
     }

    pub fn add_h(&self, tgt: usize){
        Python::attach(
            |py| -> PyResult<()>{
                self.circuit.call_method1(py, "h", (tgt,))?;
                Ok(())
            }
        ).unwrap();
     }
}

impl CliffordGates for QiskitSynIR{
    fn s(&mut self, target: synir::IndexType) {
        todo!()
    }

    fn v(&mut self, target: synir::IndexType) {
        todo!()
    }

    fn s_dgr(&mut self, target: synir::IndexType) {
        todo!()
    }

    fn v_dgr(&mut self, target: synir::IndexType) {
        todo!()
    }

    fn x(&mut self, target: synir::IndexType) {
        todo!()
    }

    fn y(&mut self, target: synir::IndexType) {
        todo!()
    }

    fn z(&mut self, target: synir::IndexType) {
        todo!()
    }

    fn h(&mut self, target: synir::IndexType) {
        todo!()
    }

    fn cx(&mut self, control: synir::IndexType, target: synir::IndexType) {
        todo!()
    }

    fn cz(&mut self, control: synir::IndexType, target: synir::IndexType) {
        todo!()
    }
}

impl Gates for QiskitSynIR {
    fn rx(&mut self, target: synir::IndexType, angle: f64) {
        todo!()
    }

    fn ry(&mut self, target: synir::IndexType, angle: f64) {
        todo!()
    }

    fn rz(&mut self, target: synir::IndexType, angle: f64) {
        todo!()
    }
}

#[pyfunction]
#[pyo3(signature = (circuit, num_qubits), text_signature = "(circuit: QuantumCircuit, num_qubits: int)")]
pub fn qiskit_to_synir(circuit: Py<PyAny>, num_qubits: usize) -> PyResult<PauliExponentialWrap> {
    let mut pe = PauliExponentialWrap::new(num_qubits);
    Python::attach(
    |py| -> PyResult<()> {
        let fun: Py<PyAny> = PyModule::from_code(
            py, 
c"
from qiskit import QuantumCircuit
from qiskit import transpile as qiskit_transpile

def transpile(circuit: QuantumCircuit) -> list[tuple[str, list[int]]]:
    circ = qiskit_transpile(circuit, basis_gates=['cx', 'h', 'rz'])
    return [(i.name, [circuit.find_bit(q).index for q in i.qubits]) for i in circ.data]
", 
            c"transpile.py", c"")?
            .getattr("transpile")?
            .into();
        let gates = fun.call1(py, (circuit,))?;
        let gates_list = gates.cast_bound::<PyList>(py)?;
        for gate in gates_list.iter() {
            let gate_tuple = gate.cast_into::<PyTuple>()?;
            let name = gate_tuple.get_item(0)?.cast_into::<PyString>()?;
            let qubits = gate_tuple.get_item(1)?.cast_into::<PyList>()?;
            println!("{:?}, {:?}", name, qubits);
            break;
            match format!("{:?}", name).as_str() {
                "cx" => {
                    let ctrl = qubits.get_item(0)?.cast::<PyInt>()?;
                    let trgt = qubits.get_item(1)?.cast::<PyInt>()?;
                    todo!("Implement the necessary functions");
                    //pe.cx(ctrl, trgt);
                },
                _ => todo!("Throw error"),
            };
        }
        Ok(())
    })?;
    Ok(pe)
}

// TODO Move below class to synpy generic stuff
#[pyclass]
pub struct PauliExponentialWrap {
    pe: PauliExponential
}

#[pymethods]
impl PauliExponentialWrap {

    #[new]
    pub fn new(num_qubits: usize) -> PauliExponentialWrap{
        let pe = PauliExponential::new(VecDeque::from(vec![]), CliffordTableau::new(num_qubits));
        PauliExponentialWrap { pe }
    }

}

impl FromPyObject<'_,'_> for PauliExponentialWrap{
    type Error = PyErr;

    fn extract(obj: Borrowed<'_, '_, PyAny>) -> Result<Self, Self::Error> {
        todo!("Impl FromPyObject for PauliExponentialWrap")
    }
}

// Keep this function here - Qiskit specific
#[pyfunction]
pub fn synthesize_to_qiskit(mut pe: PauliExponentialWrap, circuit: &mut QiskitSynIR){
    synthesize(pe, circuit);
}

// Move this function with PauliWrap - Can be used in by others.
pub fn synthesize<G>(mut pe: PauliExponentialWrap, circuit: &mut G) where G: CliffordGates + Gates{
    let mut synth = PauliExponentialSynthesizer::from_strategy(PauliPolynomialSynthStrategy::Naive, CliffordTableauSynthStrategy::PermRowCol);
    synth.synthesize(pe.pe, circuit)
}
