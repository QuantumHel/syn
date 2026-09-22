from qiskit.quantum_info import Operator
from qiskit import QuantumCircuit, transpile, qasm2
from qiskit.circuit.library import PermutationGate, QFT

from synpy.qiskit.plugin import SynPyCliffordPlugin, qiskit_to_synir
from synpy.synpy_rust import QiskitSynIR

def check_equivalent(circuit1: QuantumCircuit, circuit2: QuantumCircuit) -> bool:
    op1 = Operator.from_circuit(circuit1)
    op2 = Operator.from_circuit(circuit2)
    b = op1.equiv(op2)
    assert b
    return b

def our_transpile(circuit: QuantumCircuit, method = "Naive") -> QuantumCircuit:
    synir_result = QiskitSynIR(circuit.copy_empty_like())
    pe_wrap = qiskit_to_synir(circuit)
    pe_wrap.set_pauli_strategy(method)
    pe_wrap.set_tableau_strategy("PermRowCol")
    pe_wrap.synthesize_to_qiskit(synir_result)
    perm1 = synir_result.get_permutation()
    perm2 = [perm1.index(i) for i in range(len(perm1))]
    new_circuit = synir_result.get_circuit()
    new_circuit.append(PermutationGate(perm2), new_circuit.qubits, [])
    return synir_result.get_circuit()

original_shor = qasm2.load("factor21.qasm")
print("Original gates count", original_shor.count_ops())
print("Qubits", original_shor.qubits)
# Remove the measure-based QFT code
no_qft = original_shor.copy_empty_like()
for gate in original_shor.data: 
    if gate.name not in ["measure", "u"]:
        no_qft.compose(gate, gate.qubits, gate.clbits, inplace=True)
no_qft.count_ops()

shor_circuit = no_qft.compose(QFT(num_qubits=10), no_qft.qubits[-10:])
basic_shor = transpile(shor_circuit, basis_gates=["cx", "h", "x", "rz"], optimization_level=0)
print("Decomposed gate count",basic_shor.count_ops())
#check_equivalent(shor_circuit, basic_shor)
qiskit_shor = transpile(basic_shor, optimization_level=3, basis_gates=["cx", "h", "x", "rz"])
print("Qiskit optimized gatecount", qiskit_shor.count_ops())
#check_equivalent(shor_circuit, qiskit_shor)
syn_shor = our_transpile(basic_shor)
print("SynPy gate count", syn_shor.count_ops())
#check_equivalent(shor_circuit, syn_shor)