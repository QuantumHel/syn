"""Qiskit integration: a synthesis plugin, a circuit importer, and a coupling-map converter.

See :mod:`synpy.qiskit.plugin` for the plugin and importer, and
:mod:`synpy.qiskit.connectivity` for turning a ``Target`` or ``CouplingMap`` into the graph that
``PyPauliExponential.synthesize_to_qiskit`` takes.
"""

from synpy.qiskit.connectivity import coupling_map_to_connectivity, target_to_connectivity

__all__ = ["coupling_map_to_connectivity", "target_to_connectivity"]
