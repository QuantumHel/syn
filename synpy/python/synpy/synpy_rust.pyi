# encoding: utf-8
# module synpy.synpy_rust
"""Type stubs for the ``synpy_rust`` Rust extension module.

Docstrings here mirror the ``///`` comments on the Rust side, which pyo3 also installs as the
runtime ``__doc__`` of each object.
"""

from typing import Any

def synthesize_pauli_exponential(
    hamiltonian: list[list["PyPauliString"]],
    clifford_gates: list["PyCommand"],
    nr_qubits: int,
) -> list["PyCommand"]:
    """Synthesize a Hamiltonian plus a trailing Clifford into a flat command list.

    The one-shot alternative to building a :class:`PyPauliExponential` gate by gate.
    ``hamiltonian`` is a list of polynomials, each a list of gadgets; ``clifford_gates`` is
    replayed into the trailing tableau.

    Uses the naive strategy for both halves, so the result assumes all-to-all connectivity.

    Raises:
        Exception: If the Pauli strings differ in length, disagree with ``nr_qubits``, the
            Hamiltonian is empty, a Clifford gate names an out-of-range qubit, or
            ``clifford_gates`` contains a rotation.
    """
    ...

class PyCommand(object):
    """A single emitted gate.

    A complex enum: each variant is a callable class carrying its qubit indices, and rotations
    carry an angle in radians as the last field. Build one with e.g. ``PyCommand.CX(0, 1)`` or
    ``PyCommand.Rz(0, 0.5)``; read the fields back by index, ``cmd[0]``.

    Use :func:`synpy.utils.pycommand_to_tuple` to flatten one into ``(name, *args)``.
    """

    def __init__(self, *args: Any, **kwargs: Any) -> None: ...
    def __new__(*args: Any, **kwargs: Any) -> Any: ...

    CX: Any
    """Controlled-X, ``(control, target)``."""
    CZ: Any
    """Controlled-Z, ``(control, target)``."""
    H: Any
    """Hadamard, ``(target,)``."""
    Rx: Any
    """X rotation, ``(target, angle)``."""
    Ry: Any
    """Y rotation, ``(target, angle)``."""
    Rz: Any
    """Z rotation, ``(target, angle)``."""
    S: Any
    """Phase gate ``S``, ``(target,)``."""
    SDgr: Any
    """``S^dag``, ``(target,)``."""
    V: Any
    """``V = sqrt(X)``, ``(target,)``."""
    VDgr: Any
    """``V^dag``, ``(target,)``."""
    X: Any
    """Pauli X, ``(target,)``."""
    Y: Any
    """Pauli Y, ``(target,)``."""
    Z: Any
    """Pauli Z, ``(target,)``."""

class PyPauliString(object):
    """A Pauli string with a rotation angle -- one gadget of a Hamiltonian."""

    pauli_string: str
    """One letter from ``IXYZ`` per qubit, e.g. ``"ZZII"``."""
    phase: float
    """The rotation angle, in radians."""

    def __init__(self, pauli_string: str, phase: float) -> None: ...
    def __new__(*args: Any, **kwargs: Any) -> Any: ...
    def as_tuple(self) -> tuple[str, float]:
        """Return ``(pauli_string, phase)``."""
        ...
    def get_qubits(self) -> int:
        """Return the number of qubits, i.e. the length of the string."""
        ...

class PyCliffordTableau(object):
    """A Clifford operator on ``n`` qubits.

    Supports ``==`` and ``!=``.
    """

    def __init__(self, n: int) -> None:
        """Create the identity operation on ``n`` qubits."""
        ...
    def __new__(*args: Any, **kwargs: Any) -> Any: ...
    @staticmethod
    def from_parts(pauli_strings: list[str], signs: list[bool]) -> "PyCliffordTableau":
        """Build a tableau from its columns and signs.

        ``pauli_strings`` holds one string per qubit, each ``2n`` letters long: rows ``0..n``
        are the images of the ``X`` generators and rows ``n..2n`` those of the ``Z``
        generators. ``signs`` has ``2n`` booleans in the same order, ``True`` meaning negative.

        This is the transpose of Qiskit's ``Clifford.tableau`` layout; see
        :class:`synpy.qiskit.plugin.SynPyCliffordPlugin` for the conversion.
        """
        ...
    def synthesize(self) -> list[PyCommand]:
        """Synthesize this tableau into commands, in circuit order.

        Uses the naive strategy, so the result assumes all-to-all connectivity.
        """
        ...
    def size(self) -> int:
        """Return the number of qubits."""
        ...
    def compose(self, rhs: "PyCliffordTableau") -> "PyCliffordTableau":
        """Compose in circuit order: ``self`` runs first, then ``rhs``."""
        ...

class PyPauliExponential(object):
    """Build a circuit gate by gate, then synthesize it.

    Gates must be added in **reverse circuit order** -- last executed first. That is what
    :func:`synpy.qiskit.plugin.qiskit_to_synir` does with ``reversed(circuit.data)``, and it is
    what lets each Clifford conjugate the gadgets that follow it.

    Clifford gates never become gadgets; they are folded into the tableau on arrival.
    ``add_rz`` does the same whenever its angle turns out to be a multiple of ``pi/2``, so only
    genuinely non-Clifford rotations survive.

    ``synthesize_to_qiskit`` consumes the accumulated gates, leaving the object empty and ready
    to build up again.
    """

    def __init__(self, num_qubits: int) -> None:
        """Create an empty builder over ``num_qubits`` qubits.

        Defaults to the ``"Naive"`` polynomial strategy and ``"PermRowCol"`` tableau strategy.
        """
        ...
    def __new__(*args: Any, **kwargs: Any) -> Any: ...
    def set_pauli_strategy(self, strategy: str) -> None:
        """Choose the Pauli polynomial algorithm: ``"Naive"`` or ``"PSGS"``."""
        ...
    def set_tableau_strategy(self, strategy: str) -> None:
        """Choose the Clifford tableau algorithm: ``"Naive"`` or ``"PermRowCol"``."""
        ...
    def add_h(self, target: int) -> None:
        """Add a Hadamard on ``target``."""
        ...
    def add_v(self, target: int) -> None:
        """Add ``V = sqrt(X)`` on ``target``."""
        ...
    def add_s(self, target: int) -> None:
        """Add the phase gate ``S`` on ``target``."""
        ...
    def add_s_dgr(self, target: int) -> None:
        """Add ``S^dag`` on ``target``."""
        ...
    def add_x(self, target: int) -> None:
        """Add Pauli ``X`` on ``target``."""
        ...
    def add_y(self, target: int) -> None:
        """Add Pauli ``Y`` on ``target``."""
        ...
    def add_z(self, target: int) -> None:
        """Add Pauli ``Z`` on ``target``."""
        ...
    def add_cx(self, control: int, target: int) -> None:
        """Add a ``CX`` from ``control`` onto ``target``."""
        ...
    def add_rz(self, target: int, angle: float) -> None:
        """Add a ``Z`` rotation of ``angle`` radians on ``target``.

        A Clifford angle is folded into the tableau rather than kept as a gadget; an angle of 0
        adds nothing.
        """
        ...
    def synthesize_to_qiskit(self, circuit: "QiskitSynIR", connectivity: "PyConnectivity | None" = None) -> None:
        """Synthesize everything accumulated so far into ``circuit``, and reset this builder.

        ``connectivity`` constrains where two-qubit gates may go; ``None`` means all-to-all. Its
        qubit count must match this builder's.

        Only ``"PSGS"`` (polynomial half) and ``"PermRowCol"`` (Clifford half) consult it, and
        conformance takes both. One ``"Naive"`` half raises ``UserWarning``; two raise
        ``ValueError``, since then nothing would consult it.

        Raises
        ------
        ValueError
            If ``connectivity`` does not cover exactly this builder's qubits, is disconnected, or
            both strategies are ``"Naive"``. The builder is left intact in that case.
        """
        ...
    def print(self) -> None:
        """Print the accumulated gadgets and tableau, for debugging.

        The tableau shown is still in the inverse form the importer keeps it in.
        """
        ...

class PyConnectivity(object):
    """A hardware coupling graph: which qubit pairs may host a two-qubit gate.

    Qubits are ``0`` through ``num_qubits - 1`` and couplings are undirected, since a ``cx`` can
    be reversed with four Hadamards. Pass one to
    :meth:`PyPauliExponential.synthesize_to_qiskit` and every two-qubit gate synthesis emits will
    lie on an edge of it.

    Synthesis requires the graph to cover exactly ``0..num_qubits`` and to be connected: it routes
    two-qubit gates along paths, and there is no path between components.

    Build one from a Qiskit backend with
    :func:`synpy.qiskit.connectivity.target_to_connectivity`.
    """

    def __init__(self, num_qubits: int, edges: list[tuple[int, int]] | None = None) -> None:
        """Create a graph over ``num_qubits`` qubits, coupled by ``edges``.

        Couplings are undirected, so ``(i, j)`` and ``(j, i)`` are one edge and duplicates are
        dropped; self-loops are ignored. Qubits named by no edge are kept as isolated nodes, so the
        register width is always exactly ``num_qubits``.

        Raises
        ------
        ValueError
            If a coupling names a qubit outside ``0..num_qubits``.
        """
        ...
    def __new__(*args: Any, **kwargs: Any) -> Any: ...
    @staticmethod
    def line(num_qubits: int) -> "PyConnectivity":
        """A path graph: ``0 - 1 - ... - (num_qubits - 1)``."""
        ...
    @staticmethod
    def grid(num_rows: int, num_cols: int) -> "PyConnectivity":
        """A rectangular lattice, with qubit ``r * num_cols + c`` at row ``r``, column ``c``."""
        ...
    @staticmethod
    def complete(num_qubits: int) -> "PyConnectivity":
        """All-to-all: every qubit coupled to every other. The "no constraints" case."""
        ...
    def num_qubits(self) -> int:
        """How many qubits the graph covers."""
        ...
    def edges(self) -> list[tuple[int, int]]:
        """The couplings. Undirected, so the order within a pair means nothing."""
        ...
    def has_edge(self, i: int, j: int) -> bool:
        """Whether ``i`` and ``j`` are directly coupled."""
        ...
    def is_connected(self) -> bool:
        """Whether every qubit can reach every other one. Synthesis requires this."""
        ...
    def edge_count(self) -> int:
        """The number of couplings.

        Deliberately not ``__len__``: a graph with no couplings is still a graph, and making it
        falsy would turn ``if connectivity:`` into a silent bug.
        """
        ...
    def __repr__(self) -> str: ...

class QiskitSynIR(object):
    """A gate sink backed by a Qiskit ``QuantumCircuit``.

    Each emitted gate is forwarded as a method call on the wrapped circuit, so nothing is
    buffered and the circuit is complete as soon as synthesis returns.

    The final qubit permutation is **stored, not emitted**: rather than expanding it into swap
    gates, this sink records it for :meth:`get_permutation`. A caller that ignores it gets a
    circuit correct only up to that relabelling.
    """

    def __init__(self, qiskit_circuit: Any) -> None:
        """Wrap a Qiskit ``QuantumCircuit``, which must already have enough qubits."""
        ...
    def __new__(*args: Any, **kwargs: Any) -> Any: ...
    def s(self, target: int) -> None:
        """Append a phase gate ``S``."""
        ...
    def v(self, target: int) -> None:
        """Append ``V = sqrt(X)``, as Qiskit's ``sx``."""
        ...
    def s_dgr(self, target: int) -> None:
        """Append ``S^dag``, as Qiskit's ``sdg``."""
        ...
    def v_dgr(self, target: int) -> None:
        """Append ``V^dag``, as Qiskit's ``sxdg``."""
        ...
    def x(self, target: int) -> None:
        """Append Pauli ``X``."""
        ...
    def y(self, target: int) -> None:
        """Append Pauli ``Y``."""
        ...
    def z(self, target: int) -> None:
        """Append Pauli ``Z``."""
        ...
    def h(self, target: int) -> None:
        """Append a Hadamard."""
        ...
    def cx(self, control: int, target: int) -> None:
        """Append a ``CX``."""
        ...
    def cz(self, control: int, target: int) -> None:
        """Append a ``CZ``."""
        ...
    def add_final_permutation(self, permutation: list[int]) -> None:
        """Record the final qubit relabelling without emitting any gates."""
        ...
    def rx(self, target: int, angle: float) -> None:
        """Append an X rotation by ``angle`` radians."""
        ...
    def ry(self, target: int, angle: float) -> None:
        """Append a Y rotation by ``angle`` radians."""
        ...
    def rz(self, target: int, angle: float) -> None:
        """Append a Z rotation by ``angle`` radians."""
        ...
    def get_circuit(self) -> Any:
        """Return the wrapped ``QuantumCircuit``, with everything emitted so far."""
        ...
    def get_permutation(self) -> list[int] | None:
        """Return the final qubit relabelling, or ``None`` if none was needed.

        ``perm[physical] = logical``. This is not applied to the circuit.
        """
        ...

__all__ = [
    "PyCliffordTableau",
    "PyPauliString",
    "PyCommand",
    "synthesize_pauli_exponential",
    "PyPauliExponential",
    "PyConnectivity",
    "QiskitSynIR",
]

__loader__: Any = ...
__spec__: Any = ...
