from qiskit.quantum_info import Clifford
from qiskit import QuantumCircuit

from synpy.qiskit.plugin import SynPyPlugin


def test_qiskit_bell() -> None:
    qc = QuantumCircuit(2)
    qc.h(0)
    qc.cx(0, 1)
    cliff = Clifford(qc)

    plugin = SynPyPlugin()
    circ = plugin.run(cliff, None, None, [])
    # circ.draw()

    assert circ == qc
