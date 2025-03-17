from pauliopt.circuits import Circuit
from pauliopt.pauli.pauli_gadget import PauliGadget
from pauliopt.pauli.pauli_polynomial import PauliPolynomial
from synpy.synpy_rust import PyPauliString, PyCommand, synthesize_pauli_exponential

from synpy.utils import pycommand_to_tuple


class PauliOptSynthesizer:
    def __init__(self, pauli_polynomial_strategy: str = "NAIVE", clifford_strategy: str = "NAIVE") -> None:
        self.pauli_polynomial_strategy = pauli_polynomial_strategy
        self.clifford_strategy = clifford_strategy

    @staticmethod
    def _to_py_pauli_string(pauli_gadget: PauliGadget) -> PyPauliString:
        pauli_string = "".join([pauli.value for pauli in pauli_gadget.paulis])
        return PyPauliString(pauli_string=pauli_string, phase=pauli_gadget.angle)

    @staticmethod
    def _py_commands_to_circuit(py_commends: list[PyCommand], num_qubits: int) -> Circuit:
        circuit = Circuit(num_qubits)

        for command in py_commends:
            t = pycommand_to_tuple(command)
            gate = t[0]
            if gate == "H":
                circuit.h(t[1])
            elif gate == "X":
                circuit.x(t[1])
            elif gate == "Y":
                circuit.y(t[1])
            elif gate == "Z":
                circuit.z(t[1])
            elif gate == "S":
                circuit.s(t[1])
            elif gate == "SDgr":
                circuit.sdg(t[1])
            elif gate == "V":
                circuit.v(t[1])
            elif gate == "VDgr":
                circuit.vdg(t[1])
            elif gate == "Rx":
                circuit.rx(t[2], t[1])
            elif gate == "Ry":
                circuit.ry(t[2], t[1])
            elif gate == "Rz":
                circuit.rz(t[2], t[1])
            elif gate == "CX":
                circuit.cx(t[2], t[1])
            elif gate == "CZ":
                circuit.cz(t[2], t[1])
            else:
                raise ValueError("Unhandled gate type: " + gate)

        return circuit

    def synthesize(self, pauli_polynomial: PauliPolynomial) -> Circuit:
        hamiltonian_repr = [self._to_py_pauli_string(pauli_gadget) for pauli_gadget in pauli_polynomial]

        py_commands = synthesize_pauli_exponential(hamiltonian=[hamiltonian_repr], clifford_gates=[], nr_qubits=pauli_polynomial.num_qubits)

        return self._py_commands_to_circuit(py_commands, pauli_polynomial.num_qubits)
