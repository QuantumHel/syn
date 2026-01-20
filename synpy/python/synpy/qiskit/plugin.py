from qiskit.transpiler import CouplingMap, Target
from qiskit.transpiler.passes.synthesis.plugin import HighLevelSynthesisPlugin
from qiskit.quantum_info import Clifford
from qiskit import QuantumCircuit, transpile

from synpy.synpy_rust import PyCliffordTableau, PyPauliExponential
from synpy.utils import pycommand_to_qasm


class SynPyCliffordPlugin(HighLevelSynthesisPlugin):
    def __init__(self) -> None:
        super().__init__()

    def run(self, clifford: Clifford, coupling_map: CouplingMap, target: Target, qubits: list) -> QuantumCircuit:
        n = clifford.num_qubits
        tableau_x = clifford.tableau[:, :n]
        tableau_z = clifford.tableau[:, n : 2 * n]

        pauli_columns = ["".join(["IZXY"[2 * x + z] for x, z in zip(col_x, col_z)]) for col_x, col_z in zip(tableau_x.T, tableau_z.T)]
        signs = clifford.tableau[:, -1].tolist()

        # Convert Clifford to CliffordTableau
        synpy_tableau = PyCliffordTableau.from_parts(pauli_columns, signs)

        # Synthesize CliffordTableau
        commands = synpy_tableau.synthesize()
        qasm = pycommand_to_qasm(n, commands)
        return QuantumCircuit.from_qasm_str(qasm)


def qiskit_to_synir(circuit: QuantumCircuit) -> PyPauliExponential:
    new_circuit = transpile(circuit, basis_gates=["cx", "h", "rz"])
    pe = PyPauliExponential(new_circuit.num_qubits)

    for gate in reversed(new_circuit.data):
        qubit_indices = [new_circuit.find_bit(q).index for q in gate.qubits]
        if gate.name == "cx":
            pe.add_cx(*qubit_indices)
        elif gate.name == "h":
            pe.add_h(*qubit_indices)
        elif gate.name == "rz":
            pe.add_rz(*qubit_indices, gate.params[0])
        else:
            raise Exception("Gate is not supported")
    return pe
