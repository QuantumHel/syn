from qiskit.quantum_info import Clifford
from qiskit import QuantumCircuit

from synpy.qiskit.plugin import SynPyCliffordPlugin
from synpy.synpy_rust import QiskitSynIR


def test_qiskit_synir() -> None:
    qc = QuantumCircuit(2)
    synir = QiskitSynIR(qc)

    synir.s(0)
    synir.v(0)
    synir.s_dgr(0)
    synir.v_dgr(0)
    synir.x(0)
    synir.y(0)
    synir.z(0)
    synir.h(0)
    synir.cx(0, 1)
    synir.cz(0, 1)
    synir.rx(0, 1.23)
    synir.ry(0, 1.23)
    synir.rz(0, 1.23)


def test_qiskit_bell() -> None:
    qc = QuantumCircuit(2)
    qc.h(0)
    qc.cx(0, 1)

    cliff = Clifford(qc)

    plugin = SynPyCliffordPlugin()
    circ = plugin.run(cliff, None, None, [])
    # circ.draw()

    assert circ == qc
