from typing import Optional

import numpy as np
import pytest
from pauliopt.pauli.pauli_gadget import PauliGadget, PPhase
from pauliopt.pauli.pauli_polynomial import PauliPolynomial
from pauliopt.pauli_strings import I, X, Y, Z
from qiskit import QuantumCircuit
from qiskit.quantum_info import Operator

from synpy.pauliopt.synthesizer import PauliOptSynthesizer


def create_random_phase_gadget(num_qubits: int, min_legs: int, max_legs: int, allowed_angels: Optional[list]) -> PauliGadget:
    """
    Generate a random phase gadget.
    :param num_qubits:
    :param min_legs:
    :param max_legs:
    :param allowed_angels:
    :return:
    """
    angle = np.random.choice(allowed_angels)
    nr_legs = np.random.randint(min_legs, max_legs)
    legs = np.random.choice([i for i in range(num_qubits)], size=nr_legs, replace=False)
    phase_gadget = [I for _ in range(num_qubits)]
    for leg in legs:
        phase_gadget[leg] = np.random.choice([X, Y, Z])
    return PPhase(angle) @ phase_gadget


def generate_random_pauli_polynomial(
    num_qubits: int, num_gadgets: int, min_legs: Optional[int] = None, max_legs: Optional[int] = None, allowed_angles: Optional[list] = None
) -> PauliPolynomial:
    """
    Generate a random pauli polynomial.
    :param num_qubits:
    :param num_gadgets:
    :param min_legs:
    :param max_legs:
    :param allowed_angles:
    :return:
    """
    if min_legs is None:
        min_legs = 1
    if max_legs is None:
        max_legs = num_qubits
    if allowed_angles is None:
        allowed_angles = [2 * np.pi, np.pi, 0.5 * np.pi, 0.25 * np.pi, 0.125 * np.pi]

    pp = PauliPolynomial(num_qubits)
    for _ in range(num_gadgets):
        pp >>= create_random_phase_gadget(num_qubits, min_legs, max_legs, allowed_angles)

    return pp


def verify_equality(qc_in: QuantumCircuit, qc_out: QuantumCircuit) -> bool:
    """
    Verify the equality up to a global phase
    :param qc_in:
    :param qc_out:
    :return:
    """
    return Operator.from_circuit(qc_in).equiv(Operator.from_circuit(qc_out))


class TestBasicSyn:
    @pytest.mark.parametrize("nr_qubits", [3])
    @pytest.mark.parametrize("nr_gadgets", [3])
    def test_naive_synthesis(self, nr_qubits: int, nr_gadgets: int) -> None:
        pauli_polynomial = generate_random_pauli_polynomial(nr_qubits, nr_gadgets)

        qc = pauli_polynomial.copy().to_qiskit()

        qc_syn = PauliOptSynthesizer().synthesize(pauli_polynomial).to_qiskit()
        print()
        print(qc)
        print(qc_syn)
        assert verify_equality(qc, qc_syn)
