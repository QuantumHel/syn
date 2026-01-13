from qiskit.quantum_info import Clifford
from qiskit import QuantumCircuit

from synpy.qiskit.plugin import SynPyCliffordPlugin


def test_qiskit_bell() -> None:
    qc = QuantumCircuit(2)
    qc.h(0)
    qc.cx(0, 1)
    cliff = Clifford(qc)

    plugin = SynPyCliffordPlugin()
    circ = plugin.run(cliff, None, None, [])

    assert circ == qc

def test_qiskit_circuit() -> None:
    qc = QuantumCircuit(3)
    qc.h(0)
    qc.cx(0, 1)
    qc.t(1)
    qc.cx(1, 2)
    qc.h(2)
    qc.cz(0, 2)

    cliff = Clifford(qc)

    plugin = SynPyCliffordPlugin()
    circ = plugin.run(cliff, None, None, [])

    assert circ == qc
