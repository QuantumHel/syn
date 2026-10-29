//! The Python-facing hardware coupling graph.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use synir::architecture::connectivity::Connectivity;

/// A hardware coupling graph: which qubit pairs may host a two-qubit gate.
///
/// Qubits are `0` through `num_qubits - 1` and couplings are undirected, since a `CX` can be
/// reversed with four Hadamards. Hand one to
/// [`PyPauliExponential.synthesize_to_qiskit`](crate::wrapper::PyPauliExponential) and every
/// two-qubit gate synthesis emits will lie on an edge of it.
///
/// Synthesis requires the graph to cover exactly `0..num_qubits` and to be connected: it routes
/// two-qubit gates along paths, and there is no path between components.
///
/// ```python
/// from synpy.synpy_rust import PyConnectivity
///
/// PyConnectivity(5, [(0, 1), (1, 2), (2, 3), (3, 4)])  # a line, spelled out
/// PyConnectivity.line(5)                               # the same graph
/// PyConnectivity(5, [(0, 1)])                          # qubits 2, 3, 4 stay isolated
/// ```
///
/// To build one from a Qiskit `Target`, use `synpy.qiskit.connectivity.target_to_connectivity`.
#[pyclass]
#[derive(Clone)]
pub struct PyConnectivity {
    pub(crate) inner: Connectivity,
}

impl PyConnectivity {
    /// Wraps a `synir` graph, for the other conversions in this crate.
    pub fn from_inner(inner: Connectivity) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl PyConnectivity {
    /// `num_qubits` qubits, coupled by `edges`.
    ///
    /// Couplings are undirected, so `(i, j)` and `(j, i)` are the same edge and duplicates are
    /// dropped — which is what a symmetrized Qiskit coupling map hands over. Self-loops carry no
    /// routing information and are ignored. Qubits named by no edge are kept as isolated nodes,
    /// so the register width is always exactly `num_qubits`.
    ///
    /// Raises `ValueError` if a coupling names a qubit outside `0..num_qubits`.
    #[new]
    #[pyo3(
        signature = (num_qubits, edges=None),
        text_signature = "(num_qubits: int, edges: list[tuple[int, int]] | None = None)"
    )]
    pub fn new(num_qubits: usize, edges: Option<Vec<(usize, usize)>>) -> PyResult<Self> {
        let edges = edges.unwrap_or_default();
        Connectivity::from_edges_with_size(num_qubits, &edges)
            .map(Self::from_inner)
            .map_err(|error| PyValueError::new_err(error.to_string()))
    }

    /// A path graph: `0 - 1 - ... - (num_qubits - 1)`.
    #[staticmethod]
    pub fn line(num_qubits: usize) -> Self {
        Self::from_inner(Connectivity::line(num_qubits))
    }

    /// A rectangular lattice, with qubit `r * num_cols + c` at row `r`, column `c`.
    #[staticmethod]
    pub fn grid(num_rows: usize, num_cols: usize) -> Self {
        Self::from_inner(Connectivity::grid(num_rows, num_cols))
    }

    /// All-to-all: every qubit coupled to every other. The "no constraints" case.
    #[staticmethod]
    pub fn complete(num_qubits: usize) -> Self {
        Self::from_inner(Connectivity::complete(num_qubits))
    }

    /// How many qubits the graph covers.
    pub fn num_qubits(&self) -> usize {
        self.inner.node_count()
    }

    /// The couplings, as `(i, j)` pairs. Undirected, so the order within a pair means nothing.
    pub fn edges(&self) -> Vec<(usize, usize)> {
        self.inner.edges()
    }

    /// Whether `i` and `j` are directly coupled.
    pub fn has_edge(&self, i: usize, j: usize) -> bool {
        self.inner.has_edge(i, j)
    }

    /// Whether every qubit can reach every other one.
    ///
    /// Synthesis requires this, and rejects a graph where it does not hold.
    pub fn is_connected(&self) -> bool {
        self.inner.is_connected()
    }

    /// The number of couplings.
    ///
    /// Deliberately not `__len__`: a graph with no couplings is still a graph, and making it falsy
    /// would turn `if connectivity:` into a silent bug.
    pub fn edge_count(&self) -> usize {
        self.inner.edge_count()
    }

    pub fn __repr__(&self) -> String {
        format!(
            "PyConnectivity({}, {:?})",
            self.inner.node_count(),
            self.inner.edges()
        )
    }
}
