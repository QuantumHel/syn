import warnings
from typing import Optional
from qiskit.quantum_info import Clifford, Operator
from qiskit import QuantumCircuit, QuantumRegister
from qiskit.circuit.library import CXGate, PermutationGate, SXGate, quantum_volume
from qiskit.transpiler import CouplingMap, Target
import pytest

from math import pi

from synpy.qiskit.connectivity import coupling_map_to_connectivity, target_to_connectivity
from synpy.qiskit.plugin import SynPyCliffordPlugin, qiskit_to_synir
from synpy.synpy_rust import PyConnectivity, QiskitSynIR


def test_qiskit_synir() -> None:
    qc = QuantumCircuit(2)
    synir_obj = QiskitSynIR(qc)

    synir_obj.s(0)
    synir_obj.v(0)
    synir_obj.s_dgr(0)
    synir_obj.v_dgr(0)
    synir_obj.x(0)
    synir_obj.y(0)
    synir_obj.z(0)
    synir_obj.h(0)
    synir_obj.cx(0, 1)
    synir_obj.cz(0, 1)
    synir_obj.rx(0, 1.23)
    synir_obj.ry(0, 1.23)
    synir_obj.rz(0, 1.23)

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


all_strats = (("Naive", "Naive"), ("Naive", "PermRowCol"), ("PSGS", "Naive"), ("PSGS", "PermRowCol"))


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


def circuit_to_circuit(
    circuit: QuantumCircuit,
    pauli_strat: Optional[str] = None,
    ct_strat: Optional[str] = None,
    connectivity: Optional[PyConnectivity] = None,
) -> QuantumCircuit:
    pe_wrap = qiskit_to_synir(circuit)
    print(pe_wrap.print())
    if pauli_strat:
        pe_wrap.set_pauli_strategy(pauli_strat)
    if ct_strat:
        pe_wrap.set_tableau_strategy(ct_strat)
    synir_result = QiskitSynIR(circuit.copy_empty_like())
    pe_wrap.synthesize_to_qiskit(synir_result, connectivity)
    new_circuit = synir_result.get_circuit()
    if connectivity is not None:
        # Checked before the PermutationGate is appended below, since that instruction spans every
        # qubit and is a relabelling rather than a gate the hardware has to run.
        assert_fits_connectivity(new_circuit, connectivity)
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


def round_loop(
    circuit: QuantumCircuit,
    pauli_strat: Optional[str] = None,
    ct_strat: Optional[str] = None,
    connectivity: Optional[PyConnectivity] = None,
) -> QuantumCircuit:
    new_circuit = circuit_to_circuit(circuit, pauli_strat, ct_strat, connectivity)
    check_equiv(circuit, new_circuit)
    return new_circuit


def assert_fits_connectivity(circuit: QuantumCircuit, connectivity: PyConnectivity) -> None:
    """Every two-qubit gate in ``circuit`` lies on an edge of ``connectivity``."""
    for inst in circuit.data:
        if len(inst.qubits) != 2:
            continue
        i, j = (circuit.find_bit(qubit).index for qubit in inst.qubits)
        assert connectivity.has_edge(i, j), f"{inst.operation.name} on ({i}, {j}) is not a coupling: {connectivity}"


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
    circuit.ccx(0, 1, 2)
    circuit.ccx(0, 1, 2)
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
    circuit.rx(0.8513 * 2, 0)
    circuit.cx(0, 1)
    round_loop(circuit, pauli_strat, ct_strat)


@pytest.mark.parametrize("i", range(7))
def test_each_pauli_push(i):
    gate = [
        lambda c, i: c.s(i),
        lambda c, i: c.sdg(i),
        lambda c, i: c.rx(pi / 2, i),
        lambda c, i: c.rx(-pi / 2, i),
        lambda c, i: c.y(i),
        lambda c, i: c.x(i),
        lambda c, i: c.z(i),
    ]
    circuit = QuantumCircuit(2)
    gate[i](circuit, 0)
    circuit.cx(0, 1)
    circuit.rx(0.234234, 0)
    circuit.rz(0.53423, 1)
    round_loop(circuit, "PSGS", "Naive")


# ---------------------------------------------------------------------------
# Connectivity-constrained synthesis
# ---------------------------------------------------------------------------

# Only this pair makes the *whole* circuit conform: PSGS routes the gadgets, PermRowCol the
# trailing Clifford. A "Naive" half ignores the graph, which is covered separately below.
conforming_strats = (("PSGS", "PermRowCol"),)


def test_py_connectivity_constructors() -> None:
    assert PyConnectivity.line(4).edges() == [(0, 1), (1, 2), (2, 3)]
    assert PyConnectivity.grid(2, 3).num_qubits() == 6
    assert len(PyConnectivity.complete(4).edges()) == 6
    assert PyConnectivity.line(4).is_connected()
    assert PyConnectivity.grid(2, 2).is_connected()


def test_py_connectivity_keeps_isolated_qubits() -> None:
    # The register width is explicit, so qubits 2 and 3 survive even though no edge names them.
    connectivity = PyConnectivity(4, [(0, 1)])
    assert connectivity.num_qubits() == 4
    assert not connectivity.is_connected()


def test_py_connectivity_edges_are_undirected() -> None:
    connectivity = PyConnectivity(3, [(0, 1), (1, 0), (2, 2), (1, 2)])
    assert connectivity.edges() == [(0, 1), (1, 2)]
    assert connectivity.has_edge(1, 0) and connectivity.has_edge(0, 1)


def test_py_connectivity_rejects_out_of_range_edge() -> None:
    with pytest.raises(ValueError):
        PyConnectivity(2, [(0, 5)])


def test_target_to_connectivity_symmetrizes() -> None:
    # (3, 2) only: a Target records each direction separately, synIR's graph does not.
    target = Target(num_qubits=4)
    target.add_instruction(CXGate(), {(0, 1): None, (1, 2): None, (3, 2): None})
    connectivity = target_to_connectivity(target)
    assert connectivity is not None
    assert connectivity.edges() == [(0, 1), (1, 2), (2, 3)]
    assert connectivity.num_qubits() == 4
    assert connectivity.is_connected()


def test_target_to_connectivity_from_configuration() -> None:
    target = Target.from_configuration(
        basis_gates=["cx", "h", "rz", "x"],
        num_qubits=5,
        coupling_map=CouplingMap.from_line(5, bidirectional=False),
    )
    connectivity = target_to_connectivity(target)
    assert connectivity is not None
    assert connectivity.edges() == PyConnectivity.line(5).edges()


def test_target_to_connectivity_ideal_is_unconstrained() -> None:
    # A globally available cx constrains nothing, and None is what synthesize_to_qiskit reads as
    # all-to-all, so it can be passed straight through.
    target = Target(num_qubits=3)
    target.add_instruction(CXGate())
    assert target_to_connectivity(target) is None


def test_target_to_connectivity_without_two_qubit_gates() -> None:
    target = Target(num_qubits=3)
    target.add_instruction(SXGate(), {(0,): None, (1,): None, (2,): None})
    connectivity = target_to_connectivity(target)
    assert connectivity is not None
    assert connectivity.edges() == []
    assert not connectivity.is_connected()


def test_coupling_map_to_connectivity_keeps_idle_tail() -> None:
    connectivity = coupling_map_to_connectivity(CouplingMap([[0, 1], [1, 2]]), num_qubits=5)
    assert connectivity.num_qubits() == 5
    assert connectivity.edges() == [(0, 1), (1, 2)]


@pytest.mark.parametrize(("pauli_strat", "ct_strat"), conforming_strats)
@pytest.mark.parametrize("num_qubits", [2, 3, 4])
def test_line_connectivity_respected(pauli_strat: str, ct_strat: str, num_qubits: int) -> None:
    circuit = QuantumCircuit(num_qubits)
    circuit.h(0)
    circuit.cx(0, num_qubits - 1)
    circuit.rz(0.3, num_qubits - 1)
    circuit.rx(0.7, 0)
    round_loop(circuit, pauli_strat, ct_strat, PyConnectivity.line(num_qubits))


@pytest.mark.parametrize(("pauli_strat", "ct_strat"), conforming_strats)
def test_grid_connectivity_respected(pauli_strat: str, ct_strat: str) -> None:
    circuit = QuantumCircuit(4)
    circuit.h(0)
    circuit.cx(0, 3)
    circuit.rz(0.42, 3)
    circuit.cx(1, 2)
    circuit.rx(0.17, 2)
    round_loop(circuit, pauli_strat, ct_strat, PyConnectivity.grid(2, 2))


@pytest.mark.parametrize(("pauli_strat", "ct_strat"), conforming_strats)
def test_ring_connectivity_respected(pauli_strat: str, ct_strat: str) -> None:
    circuit = QuantumCircuit(4)
    circuit.cx(0, 2)
    circuit.rz(pi / 4, 2)
    circuit.cx(1, 3)
    circuit.rz(0.9, 3)
    round_loop(circuit, pauli_strat, ct_strat, PyConnectivity(4, [(0, 1), (1, 2), (2, 3), (3, 0)]))


@pytest.mark.parametrize(("pauli_strat", "ct_strat"), conforming_strats)
def test_quantum_volume_on_a_line(pauli_strat: str, ct_strat: str) -> None:
    round_loop(quantum_volume(3, 1), pauli_strat, ct_strat, PyConnectivity.line(3))


def test_end_to_end_from_target() -> None:
    target = Target(num_qubits=3)
    target.add_instruction(CXGate(), {(0, 1): None, (1, 2): None})
    circuit = QuantumCircuit(3)
    circuit.h(0)
    circuit.cx(0, 2)
    circuit.rz(0.55, 2)
    round_loop(circuit, "PSGS", "PermRowCol", target_to_connectivity(target))


def test_connectivity_wrong_size_raises() -> None:
    pe_wrap = qiskit_to_synir(QuantumCircuit(3))
    pe_wrap.set_pauli_strategy("PSGS")
    with pytest.raises(ValueError, match="circuit has 3 qubits"):
        pe_wrap.synthesize_to_qiskit(QiskitSynIR(QuantumCircuit(3)), PyConnectivity.line(5))


def test_disconnected_connectivity_raises() -> None:
    circuit = QuantumCircuit(4)
    circuit.cx(0, 1)
    pe_wrap = qiskit_to_synir(circuit)
    pe_wrap.set_pauli_strategy("PSGS")
    with pytest.raises(ValueError, match="disconnected"):
        pe_wrap.synthesize_to_qiskit(QiskitSynIR(circuit.copy_empty_like()), PyConnectivity(4, [(0, 1), (2, 3)]))


def test_rejected_connectivity_leaves_the_builder_intact() -> None:
    # The rejected call must not consume the accumulated gates, or the user has lost their circuit.
    circuit = QuantumCircuit(3)
    circuit.h(0)
    circuit.cx(0, 2)
    circuit.rz(0.4, 2)
    pe_wrap = qiskit_to_synir(circuit)
    pe_wrap.set_pauli_strategy("PSGS")
    with pytest.raises(ValueError):
        pe_wrap.synthesize_to_qiskit(QiskitSynIR(circuit.copy_empty_like()), PyConnectivity.line(5))

    synir_result = QiskitSynIR(circuit.copy_empty_like())
    pe_wrap.synthesize_to_qiskit(synir_result, PyConnectivity.line(3))
    new_circuit = synir_result.get_circuit()
    # The default tableau strategy is PermRowCol, so the result is only correct up to the
    # permutation the sink records instead of emitting -- the same fix-up circuit_to_circuit does.
    perm = synir_result.get_permutation()
    assert perm is not None
    new_circuit.append(PermutationGate([perm.index(i) for i in range(len(perm))]), new_circuit.qubits, [])
    check_equiv(circuit, new_circuit)


def test_connectivity_with_two_naive_strategies_raises() -> None:
    pe_wrap = qiskit_to_synir(QuantumCircuit(3))
    pe_wrap.set_pauli_strategy("Naive")
    pe_wrap.set_tableau_strategy("Naive")
    with pytest.raises(ValueError, match="both strategies"):
        pe_wrap.synthesize_to_qiskit(QiskitSynIR(QuantumCircuit(3)), PyConnectivity.line(3))


@pytest.mark.parametrize(("pauli_strat", "ct_strat"), (("Naive", "PermRowCol"), ("PSGS", "Naive")))
def test_connectivity_with_one_naive_strategy_warns(pauli_strat: str, ct_strat: str) -> None:
    # A Naive half ignores the graph, so the circuit comes out correct but not runnable. Legitimate
    # for measuring one half in isolation, hence a warning rather than an error.
    pe_wrap = qiskit_to_synir(QuantumCircuit(3))
    pe_wrap.set_pauli_strategy(pauli_strat)
    pe_wrap.set_tableau_strategy(ct_strat)
    with pytest.warns(UserWarning, match="partly respected"):
        pe_wrap.synthesize_to_qiskit(QiskitSynIR(QuantumCircuit(3)), PyConnectivity.line(3))


def test_no_connectivity_does_not_warn() -> None:
    # Every existing call site passes no connectivity and must stay silent.
    circuit = QuantumCircuit(3)
    circuit.h(0)
    circuit.cx(0, 2)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        round_loop(circuit, "Naive", "PermRowCol")
