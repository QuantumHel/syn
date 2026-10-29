"""Qiskit integration: a Clifford synthesis plugin and a circuit importer."""

from qiskit.transpiler import CouplingMap, Target
from qiskit.transpiler.passes.synthesis.plugin import HighLevelSynthesisPlugin
from qiskit.quantum_info import Clifford
from qiskit import QuantumCircuit, transpile

from math import pi

from synpy.synpy_rust import PyCliffordTableau, PyPauliExponential
from synpy.utils import pycommand_to_qasm


class SynPyCliffordPlugin(HighLevelSynthesisPlugin):
    """A Qiskit high-level-synthesis plugin that synthesizes ``Clifford`` objects with synIR.

    Register it as a ``clifford`` synthesis plugin to have Qiskit's ``HighLevelSynthesis`` pass
    route Clifford blocks through :class:`PyCliffordTableau`.

    Synthesis goes through ``PyCliffordTableau.synthesize()``, which uses the naive strategy and so
    assumes all-to-all connectivity; the ``coupling_map`` and ``target`` arguments are accepted for
    the plugin interface but not consulted. Routing them through
    :func:`synpy.qiskit.connectivity.coupling_map_to_connectivity` would need
    ``PyCliffordTableau`` to gain a PermRowCol path first -- the connectivity-aware synthesis
    reachable today goes through :func:`qiskit_to_synir` and
    ``PyPauliExponential.synthesize_to_qiskit``.
    """

    def __init__(self) -> None:
        super().__init__()

    def run(self, clifford: Clifford, coupling_map: CouplingMap, target: Target, qubits: list) -> QuantumCircuit:
        """Synthesize a Qiskit ``Clifford`` into a ``QuantumCircuit``.

        Qiskit stores its tableau as rows of ``[X | Z | sign]``; synIR stores columns per qubit.
        The transpose happens here, packing each ``(x, z)`` bit pair into a letter of ``IZXY``.

        Parameters
        ----------
        clifford
            The Clifford operator to synthesize.
        coupling_map, target
            Accepted for the plugin interface; not currently used.
        qubits
            Accepted for the plugin interface; not currently used.

        Returns
        -------
        QuantumCircuit
            An equivalent circuit over ``clifford.num_qubits`` qubits.
        """
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
    """Import a Qiskit circuit into a :class:`PyPauliExponential`, ready for synthesis.

    The circuit is first transpiled to the ``cx``/``h``/``rz``/``x`` basis if it is not already
    in it, then walked **backwards** -- ``PyPauliExponential`` expects gates last-executed-first,
    so that each Clifford conjugates the gadgets that follow it.

    ``x`` and ``rx`` are re-expressed as ``h . rz . h`` on the way in; ``s`` arrives as an
    ``rz(pi/2)``, which the importer recognizes as Clifford and folds into the tableau.

    Parameters
    ----------
    circuit
        Any circuit whose gates survive transpilation to the basis above.

    Returns
    -------
    PyPauliExponential
        A builder holding the whole circuit. Synthesize it with ``synthesize_to_qiskit``, which
        consumes it. Pass that method a connectivity -- see
        :func:`synpy.qiskit.connectivity.target_to_connectivity` -- to constrain where two-qubit
        gates may be placed; note the graph spans the whole device, so ``circuit`` must be
        device-width for a device's ``Target`` to fit it.

    Raises
    ------
    Exception
        If a gate outside the supported set survives transpilation.
    """
    # basis_gates=["cx", "h", "rz", "rx", "s"]
    basis_gates = ["cx", "h", "rz", "x"]

    if not all(gate in basis_gates for gate in circuit.count_ops()):
        new_circuit = transpile(circuit, basis_gates=basis_gates, routing_method="none", optimization_level=0)
    else:
        new_circuit = circuit
    pe = PyPauliExponential(new_circuit.num_qubits)
    for gate in reversed(new_circuit.data):
        qubit_indices = [new_circuit.find_bit(q).index for q in gate.qubits]
        if gate.name == "cx":
            pe.add_cx(*qubit_indices)
        elif gate.name == "h":
            pe.add_h(*qubit_indices)
        elif gate.name == "rz":
            pe.add_rz(*qubit_indices, gate.params[0])
        elif gate.name == "s":
            pe.add_rz(*qubit_indices, pi / 2)
        elif gate.name == "x":
            pe.add_h(*qubit_indices)
            pe.add_rz(*qubit_indices, pi)
            pe.add_h(*qubit_indices)
        elif gate.name == "rx":
            pe.add_h(*qubit_indices)
            pe.add_rz(*qubit_indices, gate.params[0])
            pe.add_h(*qubit_indices)
        else:
            raise Exception("Gate is not supported")
    return pe
