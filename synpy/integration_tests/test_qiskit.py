from qiskit.quantum_info import Clifford, Operator
from qiskit import QuantumCircuit


from synpy.qiskit.plugin import SynPyCliffordPlugin, qiskit_to_synir
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

    reference_circuit = ["s", "sx", "sdg", "sxdg", "x", "y", "z", "h", "cx", "cz", "rx", "ry", "rz"]
    reference_angles = [None, None, None, None, None, None, None, None, None, None, 1.23, 1.23, 1.23]

    for inst in qc.data:
        assert inst.name == reference_circuit.pop(0)
        reference_param = reference_angles.pop(0)
        if inst.params:
            assert inst.params[0] == reference_param


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

    assert circ == qc

def test_qiskit_loop() -> None:
    circuit = QuantumCircuit(3)
    circuit.h(0)
    circuit.cx(0, 1)
    circuit.rz(1.5, 1)

    import synpy
    print(synpy.__file__)
    print(dir(synpy.qiskit.plugin))

    pe_wrap = qiskit_to_synir(circuit)
    
    synir_result = QiskitSynIR(QuantumCircuit(3))
    pe_wrap.synthesize_to_qiskit(synir_result)
    circuit = synir_result.get_circuit()
    
    op = Operator.from_circuit(circuit)

    sample_circuit = QuantumCircuit(3)
    sample_circuit.h(0)
    sample_circuit.cx(0, 1)
    sample_circuit.rz(1.5, 1)

    op2 = Operator.from_circuit(sample_circuit)
    assert op.equiv(op2)