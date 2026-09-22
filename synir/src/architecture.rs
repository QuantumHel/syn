//! Hardware connectivity: which qubit pairs may host a two-qubit gate.
//!
//! Real devices do not let every qubit talk to every other one. A `CX` between distant qubits has
//! to be routed through intervening ones, which costs gates. The synthesizers in
//! [`crate::ir`] take an [`Architecture`] so they can account for that *while* choosing their
//! structure, rather than emitting an all-to-all circuit and paying a routing pass afterwards.
//!
//! [`Connectivity`](connectivity::Connectivity) is the implementation; it precomputes
//! all-pairs shortest paths and caches which vertices may be removed without disconnecting the
//! graph.

pub mod connectivity;

/// A qubit, as a node index into the connectivity graph.
type GraphIndex = usize;
/// The cost of routing across an edge. Unweighted graphs use 1.
type EdgeWeight = usize;
/// Nodes carry no data; the index is the qubit.
type NodeWeight = ();

/// Why [`Architecture::get_cx_ladder`] could not build a ladder.
#[derive(Debug, PartialEq)]
pub enum LadderError {
    /// The requested root is not in the tree spanning the given nodes.
    RootNotFound,
    /// These nodes are not in the graph at all.
    NodesNotFound(Vec<GraphIndex>),
}

/// A hardware coupling graph, as synthesis needs to query it.
///
/// Implemented by [`Connectivity`](connectivity::Connectivity).
pub trait Architecture {
    /// A shortest path from `i` to `j`, inclusive of both endpoints.
    ///
    /// Returns an empty vector if the two are in different components.
    fn best_path(&self, i: GraphIndex, j: GraphIndex) -> Vec<GraphIndex>;

    /// The total edge weight along a shortest path from `i` to `j`.
    fn distance(&self, i: GraphIndex, j: GraphIndex) -> usize;

    /// The qubits directly coupled to `i`.
    fn neighbors(&self, i: GraphIndex) -> Vec<GraphIndex>;

    /// The qubits that can be removed without disconnecting the graph.
    ///
    /// These are the non-articulation vertices. Synthesis eliminates qubits one at a time and
    /// needs the survivors to stay mutually reachable, so a pivot heuristic picks from this set.
    fn non_cutting(&self) -> &Vec<GraphIndex>;

    /// A CNOT ladder over `nodes`, rooted at `root`.
    ///
    /// Returns the edges of a breadth-first traversal of a Steiner tree spanning `nodes` —
    /// architecture-conforming pairs that fold the given qubits together onto the root.
    ///
    /// # Errors
    ///
    /// [`LadderError::NodesNotFound`] if any of `nodes` is absent from the graph, or
    /// [`LadderError::RootNotFound`] if `root` is not in the resulting tree.
    fn get_cx_ladder(
        &self,
        nodes: &[GraphIndex],
        root: &GraphIndex,
    ) -> Result<Vec<(GraphIndex, GraphIndex)>, LadderError>;

    /// A copy of the graph with qubit `i` removed.
    fn disconnect(&self, i: GraphIndex) -> Self;
}
