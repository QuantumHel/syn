from qiskit.quantum_info import Clifford, Operator
from qiskit import QuantumCircuit, QuantumRegister
from qiskit.circuit.library import PermutationGate, QuantumVolume


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
    round_loop(circuit)


def check_equiv(circuit: QuantumCircuit, circuit2: QuantumCircuit):
    op1 = Operator.from_circuit(circuit)
    op2 = Operator.from_circuit(circuit2)
    assert op1.equiv(op2)

def circuit_to_circuit(circuit: QuantumCircuit) -> QuantumCircuit:
    pe_wrap = qiskit_to_synir(circuit)
    synir_result = QiskitSynIR(circuit.copy_empty_like())
    pe_wrap.synthesize_to_qiskit(synir_result)
    new_circuit = synir_result.get_circuit()
    new_circuit.append(PermutationGate(synir_result.get_permutation()), new_circuit.qubits, new_circuit.clbits)
    return synir_result.get_circuit()


def test_qiskit_multiple_registers():
    reg1 = QuantumRegister(1)
    reg2 = QuantumRegister(1)
    circuit = QuantumCircuit(reg1, reg2)
    circuit.cx(reg1, reg2)
    round_loop(circuit)

def round_loop(circuit):
    new_circuit = circuit_to_circuit(circuit)
    check_equiv(circuit, new_circuit)

def test_rz_at_start_of_circuit():
    circuit = QuantumCircuit(2)
    circuit.rz(0.234, 0)
    circuit.cx(0,1)
    round_loop(circuit)

def test_quantum_volume():
    circuit = QuantumVolume(3)
    round_loop(circuit)