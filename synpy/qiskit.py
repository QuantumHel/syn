from qiskit.transpiler.passes.synthesis.plugin import HighLevelSynthesisPlugin
from qiskit.quantum_info import Clifford
from qiskit import QuantumCircuit

from qiskit.transpiler import CouplingMap, Target


from synpy import CliffordTableau

class SynPyPlugin(HighLevelSynthesisPlugin):
    def __init__(self):
        super().__init__()

    def run(self, clifford: Clifford, coupling_map: CouplingMap, target: Target, qubits: list) -> QuantumCircuit:
        n = clifford.num_qubits
        tableau_x = clifford.tableau[:, :n]
        tableau_z = clifford.tableau[:, n:2*n]

        pauli_columns = [
            "".join(
                ["IZXY"[2 * x + z] for x, z in zip(col_x, col_z)]
            )
            for col_x, col_z in zip(tableau_x.T, tableau_z.T)
        ]
        signs = clifford.tableau[:, -1].tolist()

        # Convert Clifford to CliffordTableau
        synpy_tableau = CliffordTableau.from_parts(pauli_columns, signs, n)
        
        # Synthesize CliffordTableau
        qasm = synpy_tableau.synthesize()
        return QuantumCircuit.from_qasm_str('\n'.join(qasm))


qc = QuantumCircuit(2)
qc.h(0)
qc.cx(0, 1)
cliff = Clifford(qc)

plugin = SynPyPlugin()
circ = plugin.run(cliff, None, None, None)
circ.draw()

assert circ == qc