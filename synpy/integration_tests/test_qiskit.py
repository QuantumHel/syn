from typing import Optional
from qiskit.quantum_info import Clifford, Operator
from qiskit import QuantumCircuit, QuantumRegister
from qiskit.circuit.library import PermutationGate, quantum_volume
import pytest

from math import pi

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
def test_qiskit_loop(pauli_strat: str, ct_strat: str) -> None:
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
    # circuit.y(2)
    round_loop(circuit, pauli_strat, ct_strat)
    print("fail for check")
    # assert False


def check_equiv(circuit: QuantumCircuit, circuit2: QuantumCircuit) -> None:
    op1 = Operator.from_circuit(circuit)
    op2 = Operator.from_circuit(circuit2)
    check = op1.equiv(op2)
    if not check:
        print(circuit2)
        qiskit_to_synir(circuit2).print()
        print("Equivalent with smaller tolerance?", op1.equiv(op2, 1e-4, 1e-4))
    assert check


def circuit_to_circuit(circuit: QuantumCircuit, pauli_strat: Optional[str] = None, ct_strat: Optional[str] = None) -> QuantumCircuit:
    pe_wrap = qiskit_to_synir(circuit)
    print(pe_wrap.print())
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
def test_qiskit_multiple_registers(pauli_strat: str, ct_strat: str) -> None:
    reg1 = QuantumRegister(1)
    reg2 = QuantumRegister(1)
    circuit = QuantumCircuit(reg1, reg2)
    circuit.cx(reg1, reg2)
    round_loop(circuit, pauli_strat, ct_strat)


def round_loop(circuit: QuantumCircuit, pauli_strat: Optional[str] = None, ct_strat: Optional[str] = None) -> QuantumCircuit:
    new_circuit = circuit_to_circuit(circuit, pauli_strat, ct_strat)
    check_equiv(circuit, new_circuit)
    return new_circuit


@pytest.mark.parametrize(("pauli_strat", "ct_strat"), all_strats)
def test_rz_at_start_of_circuit(pauli_strat: str, ct_strat: str) -> None:
    circuit = QuantumCircuit(2)
    circuit.rz(0.234, 0)
    circuit.cx(0, 1)
    round_loop(circuit, pauli_strat, ct_strat)


@pytest.mark.parametrize(("pauli_strat", "ct_strat"), all_strats)
def test_quantum_volume2(pauli_strat: str, ct_strat: str) -> None:
    circuit = quantum_volume(2, 1)
    round_loop(circuit, pauli_strat, ct_strat)


@pytest.mark.parametrize(("pauli_strat", "ct_strat"), all_strats)
def test_quantum_volume3(pauli_strat: str, ct_strat: str) -> None:
    circuit = quantum_volume(3, 1)
    round_loop(circuit, pauli_strat, ct_strat)


@pytest.mark.parametrize(("pauli_strat", "ct_strat"), all_strats)
def test_quantum_volume4(pauli_strat: str, ct_strat: str) -> None:
    circuit = quantum_volume(4, 1)
    round_loop(circuit, pauli_strat, ct_strat)

qv5_circuit = quantum_volume(5)

@pytest.mark.parametrize(("pauli_strat", "ct_strat"), all_strats)
def test_quantum_volume5(pauli_strat: str, ct_strat: str) -> None:
    circuit = qv5_circuit
    round_loop(circuit, pauli_strat, ct_strat)

@pytest.mark.parametrize(("pauli_strat", "ct_strat"), all_strats)
def test_final_permutation(pauli_strat: str, ct_strat: str) -> None:
    circuit = QuantumCircuit(5)
    circuit.append(PermutationGate([4, 2, 3, 0, 1]), circuit.qubits)
    round_loop(circuit, pauli_strat, ct_strat)


@pytest.mark.parametrize(("pauli_strat", "ct_strat"), all_strats)
def test_toffoli(pauli_strat: str, ct_strat: str) -> None:
    circuit = QuantumCircuit(3)
    circuit.ccx(0, 1, 2)
    round_loop(circuit, pauli_strat, ct_strat)


@pytest.mark.parametrize(("pauli_strat", "ct_strat"), all_strats)
def test_sqg(pauli_strat: str, ct_strat: str) -> None:
    circuit = QuantumCircuit(1)
    circuit.t(0)
    round_loop(circuit, pauli_strat, ct_strat)


@pytest.mark.parametrize(("pauli_strat", "ct_strat"), all_strats)
def test_pauli_gadget(pauli_strat: str, ct_strat: str) -> None:
    circuit = QuantumCircuit(3)
    circuit.cx(0, 1)
    circuit.cx(1, 2)
    circuit.t(2)
    circuit.cx(1, 2)
    circuit.cx(0, 1)
    round_loop(circuit, pauli_strat, ct_strat)


@pytest.mark.parametrize(("pauli_strat", "ct_strat"), all_strats)
def test_pauli_gadget_half(pauli_strat: str, ct_strat: str) -> None:
    circuit = QuantumCircuit(3)
    circuit.cx(0, 1)
    circuit.cx(1, 2)
    circuit.t(2)
    round_loop(circuit, pauli_strat, ct_strat)

@pytest.mark.parametrize(("pauli_strat", "ct_strat"), all_strats)
def test_merge_repeats(pauli_strat: str, ct_strat: str):
    circuit = QuantumCircuit(3)
    circuit.ccx(0,1,2)
    circuit.ccx(0,1,2)
    new_circuit = round_loop(circuit, pauli_strat, ct_strat)
    all_ops = new_circuit.count_ops() 
    if ct_strat == "Naive":
        assert len(all_ops) == 0
    else:
        assert len(all_ops) == 1
        assert "permutation" in all_ops
 
@pytest.mark.parametrize(("pauli_strat", "ct_strat"), all_strats)
def test_2sqg(pauli_strat: str, ct_strat: str):
    circuit = QuantumCircuit(2)
    circuit.rz(0.59156, 0)
    circuit.rx(0.8513, 0)
    circuit.cx(0, 1)
    circuit.rx(0.8513*2, 0)
    circuit.cx(0, 1)
    round_loop(circuit, pauli_strat, ct_strat)

@pytest.mark.parametrize("i", range(7))
def test_each_pauli_push(i):
    gate = [
        lambda c, i: c.s(i),
        lambda c, i: c.sdg(i),
        lambda c, i: c.rx(pi/2, i),
        lambda c, i: c.rx(-pi/2, i),
        lambda c, i: c.y(i),
        lambda c, i: c.x(i),
        lambda c, i: c.z(i),
    ]
    circuit = QuantumCircuit(2)
    gate[i](circuit, 0)
    circuit.cx(0,1)
    circuit.rx(0.234234, 0)
    circuit.rz(0.53423, 1)
    round_loop(circuit, "PSGS", "Naive")