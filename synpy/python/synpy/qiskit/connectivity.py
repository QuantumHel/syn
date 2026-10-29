"""Turning Qiskit hardware descriptions into synIR coupling graphs."""

from typing import Iterable, Optional

from qiskit.transpiler import CouplingMap, Target

from synpy.synpy_rust import PyConnectivity


def _undirected_edges(edges: Iterable[tuple[int, int]]) -> list[tuple[int, int]]:
    """Symmetrize and deduplicate directed couplings.

    A Qiskit coupling map is directed -- ``(0, 1)`` means a ``cx`` with control 0 -- while synIR
    routes over an undirected graph, since a ``cx`` can be reversed by conjugating it with
    Hadamards. A pair coupled in *either* direction therefore becomes one edge.

    Self-loops are dropped: they carry no routing information.
    """
    return sorted({(a, b) if a < b else (b, a) for a, b in edges if a != b})


def coupling_map_to_connectivity(coupling_map: CouplingMap, num_qubits: Optional[int] = None) -> PyConnectivity:
    """Convert a Qiskit ``CouplingMap`` into a :class:`PyConnectivity`.

    Parameters
    ----------
    coupling_map
        The device couplings. Treated as undirected; see :func:`_undirected_edges`.
    num_qubits
        The register width, defaulting to ``coupling_map.size()``. Pass it explicitly when the map
        omits idle qubits at the end of the register, so that those qubits keep their indices --
        synthesis addresses qubits by number, so a shifted index is a wrong circuit.

    Returns
    -------
    PyConnectivity
        A graph over ``num_qubits`` qubits.
    """
    if num_qubits is None:
        num_qubits = coupling_map.size()
    return PyConnectivity(num_qubits, _undirected_edges(coupling_map.get_edges()))


def target_to_connectivity(target: Target) -> Optional[PyConnectivity]:
    """Convert a Qiskit ``Target`` into a :class:`PyConnectivity`, or ``None`` if unconstrained.

    ``None`` means the target places no connectivity constraint at all: it has no instructions
    with qargs, or its two-qubit gates are declared globally available rather than per pair, or it
    has operations on more than two qubits. ``None`` is also what
    ``PyPauliExponential.synthesize_to_qiskit`` takes to mean all-to-all, so the result can be
    passed straight through either way.

    A target that *does* constrain connectivity but has no two-qubit gates yields an edgeless
    graph. Synthesis rejects that with ``ValueError`` for being disconnected, which is the honest
    answer: such a device cannot host a ``cx``.

    Notes
    -----
    The graph spans ``target.num_qubits`` -- the whole device. The circuit handed to
    :func:`~synpy.qiskit.plugin.qiskit_to_synir` must therefore be device-width too, or
    ``synthesize_to_qiskit`` raises ``ValueError`` over the mismatch. Choosing *which* physical
    qubits a smaller circuit should occupy is layout, a separate concern: run Qiskit's layout pass
    first and synthesize the device-width result, or build a :class:`PyConnectivity` over just the
    qubits you intend to use.

    Parameters
    ----------
    target
        The backend description to read couplings from.

    Returns
    -------
    PyConnectivity or None
        The device's coupling graph, or ``None`` if it constrains nothing.
    """
    coupling_map = target.build_coupling_map()
    if coupling_map is None:
        return None
    return coupling_map_to_connectivity(coupling_map, target.num_qubits)
