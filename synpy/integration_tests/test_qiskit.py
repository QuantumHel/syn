from qiskit.quantum_info import Clifford, Operator
from qiskit import QuantumCircuit, QuantumRegister
from qiskit.circuit.library import PermutationGate, QuantumVolume
import pytest

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

all_strats = (
    ("Naive", "Naive"),
    ("Naive", "PermRowCol"),
    ("PSGS", "Naive"),
    ("PSGS", "PermRowCol")
)

@pytest.mark.parametrize(("pauli_strat", "ct_strat"), all_strats)
def test_qiskit_loop(pauli_strat, ct_strat) -> None:
    circuit = QuantumCircuit(3)
    circuit.s(0)
    circuit.cx(0, 1)
    circuit.rz(1.5, 1)
    circuit.s(0)
    circuit.cx(0, 1)
    circuit.rz(1.5, 1)
    circuit.s(0)
    circuit.cx(0, 1)
    circuit.rz(1.5, 1)
    round_loop(circuit, pauli_strat, ct_strat)

def check_equiv(circuit: QuantumCircuit, circuit2: QuantumCircuit):
    op1 = Operator.from_circuit(circuit)
    op2 = Operator.from_circuit(circuit2)
    print(circuit)
    # if circuit.num_qubits < 4:
    #     print("Is not the same as")
    #     print(circuit2)
    assert op1.equiv(op2)

def circuit_to_circuit(circuit: QuantumCircuit, pauli_strat = None, ct_strat = None) -> QuantumCircuit:
    pe_wrap = qiskit_to_synir(circuit)
    if pauli_strat:
        pe_wrap.set_pauli_strategy(pauli_strat)
    if ct_strat:
        pe_wrap.set_tableau_strategy(ct_strat)
    synir_result = QiskitSynIR(circuit.copy_empty_like())
    pe_wrap.synthesize_to_qiskit(synir_result)
    new_circuit = synir_result.get_circuit()
    if ct_strat != "Naive":
        perm1 = synir_result.get_permutation()
        perm2 = [perm1.index(i) for i in range(len(perm1))]
        new_circuit.append(PermutationGate(perm2), new_circuit.qubits, [])
    return synir_result.get_circuit()


@pytest.mark.parametrize(("pauli_strat", "ct_strat"), all_strats)
def test_qiskit_multiple_registers(pauli_strat, ct_strat):
    reg1 = QuantumRegister(1)
    reg2 = QuantumRegister(1)
    circuit = QuantumCircuit(reg1, reg2)
    circuit.cx(reg1, reg2)
    round_loop(circuit, pauli_strat, ct_strat)

def round_loop(circuit, pauli_strat = None, ct_strat = None):
    new_circuit = circuit_to_circuit(circuit, pauli_strat, ct_strat)
    check_equiv(circuit, new_circuit)
    return new_circuit

@pytest.mark.parametrize(("pauli_strat", "ct_strat"), all_strats)
def test_rz_at_start_of_circuit(pauli_strat, ct_strat):
    circuit = QuantumCircuit(2)
    circuit.rz(0.234, 0)
    circuit.cx(0,1)
    round_loop(circuit, pauli_strat, ct_strat)

@pytest.mark.parametrize(("pauli_strat", "ct_strat"), all_strats)
def test_quantum_volume2(pauli_strat, ct_strat):
    circuit = QuantumVolume(2,1)

    round_loop(circuit, pauli_strat, ct_strat)

@pytest.mark.parametrize(("pauli_strat", "ct_strat"), all_strats)
def test_quantum_volume3(pauli_strat, ct_strat):
    circuit = QuantumVolume(3,1)
    round_loop(circuit, pauli_strat, ct_strat)

@pytest.mark.parametrize(("pauli_strat", "ct_strat"), all_strats)
def test_quantum_volume4(pauli_strat, ct_strat):
    circuit = QuantumVolume(4)
    round_loop(circuit, pauli_strat, ct_strat)

@pytest.mark.parametrize(("pauli_strat", "ct_strat"), all_strats)
def test_final_permutation(pauli_strat, ct_strat):
    circuit = QuantumCircuit(5)
    circuit.append(PermutationGate([4,2,3,0,1]), circuit.qubits)
    round_loop(circuit, pauli_strat, ct_strat)

@pytest.mark.parametrize(("pauli_strat", "ct_strat"), all_strats)
def test_toffoli(pauli_strat, ct_strat):
    circuit = QuantumCircuit(3)
    circuit.ccx(0,1,2)
    round_loop(circuit, pauli_strat, ct_strat)

@pytest.mark.parametrize(("pauli_strat", "ct_strat"), all_strats)
def test_sqg(pauli_strat, ct_strat):
    circuit = QuantumCircuit(1)
    circuit.t(0)
    round_loop(circuit, pauli_strat, ct_strat)

@pytest.mark.parametrize(("pauli_strat", "ct_strat"), all_strats)
def test_pauli_gadget(pauli_strat, ct_strat):
    circuit = QuantumCircuit(3)
    circuit.cx(0,1)
    circuit.cx(1,2)
    circuit.t(2)
    circuit.cx(1,2)
    circuit.cx(0,1)
    round_loop(circuit, pauli_strat, ct_strat)


@pytest.mark.parametrize(("pauli_strat", "ct_strat"), all_strats)
def test_pauli_gadget_half(pauli_strat, ct_strat):
    circuit = QuantumCircuit(3)
    circuit.cx(0,1)
    circuit.cx(1,2)
    circuit.t(2)
    round_loop(circuit, pauli_strat, ct_strat)