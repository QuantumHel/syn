"""Synthesize a Shor factoring circuit with synIR, all-to-all and onto real hardware.

The second half targets IQM's 20-qubit Apollo/Garnet coupling graph, so the circuit that comes
out can be run without a routing pass afterwards: every two-qubit gate already sits on a coupling
the device has.

Run from anywhere::

    python "examples/shor demo.py"
"""

from pathlib import Path

from qiskit import QuantumCircuit, qasm2, transpile
from qiskit.circuit.library import QFT, PermutationGate
from qiskit.quantum_info import Operator

# Apollo carries the same 20-qubit coupling map as IQM's Garnet, and unlike `IQMFakeGarnet` it
# is exported from the package root.
from iqm.qiskit_iqm import IQMFakeApollo

from synpy.qiskit.connectivity import target_to_connectivity
from synpy.qiskit.plugin import qiskit_to_synir
from synpy.synpy_rust import PyConnectivity, QiskitSynIR

HERE = Path(__file__).parent


def check_equivalent(circuit1: QuantumCircuit, circuit2: QuantumCircuit) -> bool:
    op1 = Operator.from_circuit(circuit1)
    op2 = Operator.from_circuit(circuit2)
    b = op1.equiv(op2)
    assert b
    return b


def our_transpile(
    circuit: QuantumCircuit,
    method: str = "Naive",
    connectivity: PyConnectivity | None = None,
) -> QuantumCircuit:
    """Synthesize `circuit` with synIR, optionally constrained to a coupling graph.

    Pass `connectivity` to keep every two-qubit gate on a hardware coupling. That needs the
    architecture-aware strategies on both halves -- `method="PSGS"` here and `"PermRowCol"`
    below -- since a `"Naive"` half routes wherever it likes.
    """
    synir_result = QiskitSynIR(circuit.copy_empty_like())
    pe_wrap = qiskit_to_synir(circuit)
    pe_wrap.set_pauli_strategy(method)
    pe_wrap.set_tableau_strategy("PermRowCol")
    pe_wrap.synthesize_to_qiskit(synir_result, connectivity)
    perm1 = synir_result.get_permutation()
    perm2 = [perm1.index(i) for i in range(len(perm1))]
    new_circuit = synir_result.get_circuit()
    new_circuit.append(PermutationGate(perm2), new_circuit.qubits, [])
    return synir_result.get_circuit()


def count_violations(circuit: QuantumCircuit, connectivity: PyConnectivity) -> int:
    """How many two-qubit gates sit on a pair the hardware does not couple."""
    violations = 0
    for inst in circuit.data:
        if len(inst.qubits) != 2:
            continue
        i, j = (circuit.find_bit(qubit).index for qubit in inst.qubits)
        if not connectivity.has_edge(i, j):
            violations += 1
    return violations


def report(
    label: str, circuit: QuantumCircuit, connectivity: PyConnectivity | None = None
) -> None:
    ops = circuit.count_ops()
    line = f"{label:<34} total {sum(ops.values()):>6}   cx {ops.get('cx', 0):>6}"
    if connectivity is not None:
        line += f"   off-graph 2q gates {count_violations(circuit, connectivity):>5}"
    print(line)


original_shor = qasm2.load(str(HERE / "factor21.qasm"))
print("Original gate count", dict(original_shor.count_ops()))
print("Qubits", original_shor.num_qubits)

# Remove the measure-based QFT code
no_qft = original_shor.copy_empty_like()
for gate in original_shor.data:
    if gate.name not in ["measure", "u"]:
        no_qft.compose(gate, gate.qubits, gate.clbits, inplace=True)

shor_circuit = no_qft.compose(QFT(num_qubits=10), no_qft.qubits[-10:])
basic_shor = transpile(
    shor_circuit, basis_gates=["cx", "h", "x", "rz"], optimization_level=0
)
report("Decomposed", basic_shor)
# check_equivalent(shor_circuit, basic_shor)

qiskit_shor = transpile(
    basic_shor, optimization_level=3, basis_gates=["cx", "h", "x", "rz"]
)
report("Qiskit optimized (all-to-all)", qiskit_shor)
# check_equivalent(shor_circuit, qiskit_shor)

syn_shor = our_transpile(basic_shor)
report("SynPy (all-to-all)", syn_shor)
# check_equivalent(shor_circuit, syn_shor)

# ---------------------------------------------------------------------------
# On hardware: IQM's 20-qubit Apollo/Garnet coupling graph
# ---------------------------------------------------------------------------

device = IQMFakeApollo()
connectivity = target_to_connectivity(device.target)
assert (
    connectivity is not None
), "the device constrains connectivity, so this cannot be None"
print(
    f"\n{device.name}: {connectivity.num_qubits()} qubits, "
    f"{connectivity.edge_count()} couplings"
)

# The coupling graph spans the whole device, so the circuit has to as well: synthesis addresses
# qubits by their physical number. Choosing *which* physical qubits a narrower circuit should
# occupy is layout, a separate concern -- here the identity layout will do.
wide_shor = QuantumCircuit(connectivity.num_qubits())
wide_shor.compose(basic_shor, qubits=list(range(basic_shor.num_qubits)), inplace=True)

coupling_map = device.target.build_coupling_map()
for level in (1, 3):
    report(
        f"Qiskit routed to device (opt {level})",
        transpile(
            wide_shor,
            coupling_map=coupling_map,
            basis_gates=["cx", "h", "x", "rz"],
            optimization_level=level,
            seed_transpiler=7,
        ),
        connectivity,
    )

# No routing pass afterwards: synthesis picks its CNOT structure from the coupling graph itself.
device_shor = our_transpile(wide_shor, "PSGS", connectivity)
report("SynPy synthesized for device", device_shor, connectivity)

violations = count_violations(device_shor, connectivity)
assert violations == 0, f"{violations} two-qubit gates are not device couplings"
print("\nEvery two-qubit gate lies on a device coupling.")

opt_device_shor = transpile(device_shor, basis_gates=["cx", "h", "x", "rz"], routing_method="none", optimization_level=3)
report("SynPy synthesized, Qiskit optimized", opt_device_shor, connectivity)