"""Python bindings for synIR: Pauli-based quantum circuit synthesis.

synIR rewrites a circuit into the normal form ``U = T * U_pp`` -- a run of Pauli gadgets
followed by a single Clifford operator -- and then synthesizes both halves back into gates,
optionally respecting a hardware coupling map.

The Rust extension module is :mod:`synpy.synpy_rust`; the subpackages here adapt it to
existing frameworks.

Subpackages
-----------
``synpy.qiskit``
    A Qiskit high-level-synthesis plugin and a ``QuantumCircuit`` importer.
``synpy.pauliopt``
    A synthesizer taking a ``pauliopt`` ``PauliPolynomial``.
``synpy.utils``
    Conversions between :class:`PyCommand` and tuples or OpenQASM 2.

Examples
--------
Synthesize a Qiskit circuit::

    from qiskit import QuantumCircuit
    from synpy.synpy_rust import QiskitSynIR
    from synpy.qiskit.plugin import qiskit_to_synir

    pe = qiskit_to_synir(circuit)
    pe.set_pauli_strategy("PSGS")

    sink = QiskitSynIR(QuantumCircuit(circuit.num_qubits))
    pe.synthesize_to_qiskit(sink)
    result = sink.get_circuit()

The result is correct up to ``sink.get_permutation()``, a final qubit relabelling that the
sink records rather than emitting as swap gates.
"""
