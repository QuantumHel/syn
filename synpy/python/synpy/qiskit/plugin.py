from qiskit.transpiler import CouplingMap, Target
from qiskit.transpiler.passes.synthesis.plugin import HighLevelSynthesisPlugin
from qiskit.quantum_info import Clifford
from qiskit import QuantumCircuit, transpile

from synpy.synpy_rust import PyCliffordTableau, QiskitSynIR, PauliExponentialWrap, synthesize_to_qiskit
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
    
def alt_qiskit_to_synir(circuit:QuantumCircuit) -> PauliExponentialWrap:
    new_circuit = transpile(circuit, basis_gates=['cx', 'h', 'rz'])
    pe = PauliExponentialWrap(new_circuit.num_qubits)
    for gate in new_circuit.data:
        # TODO add gate to pe
        if gate.name == 'cx':
            pass
        # etc
        pass
    print(pe)
    return pe

def main():
    circuit = QuantumCircuit(3)
    synir = QiskitSynIR(circuit)
    synir.add_h(0)
    print(circuit)
    # Use python for circuit parsing
    pe_wrap2 = alt_qiskit_to_synir(circuit)
    # Do use rust for synthesis.
    # synthesize(pe_wrap2, synir)
    synir_result = QiskitSynIR(QuantumCircuit(3))
    synthesize_to_qiskit(pe_wrap2, synir_result)
    print(synir_result)
