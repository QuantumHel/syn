//! A concrete coupling graph.
//!
//! See [`Connectivity`].

use super::{Architecture, EdgeWeight, GraphIndex, LadderError, NodeWeight};
use itertools::Itertools;
use petgraph::algo::floyd_warshall::floyd_warshall_path;
use petgraph::algo::steiner_tree::stable_steiner_tree;
use petgraph::prelude::{EdgeRef, StableUnGraph};
use petgraph::visit::{Bfs, IntoEdgeReferences, VisitMap, Visitable};
use petgraph::{
    algo::articulation_points::articulation_points,
    graph::NodeIndex,
    visit::{IntoNodeReferences, NodeIndexable, NodeRef},
};
use std::collections::HashMap;

/// Get all the vertices in a graph that are non-cutting (won't make the graph disconnected)
fn get_non_cutting_vertices(
    graph: &StableUnGraph<NodeWeight, EdgeWeight, GraphIndex>,
) -> Vec<GraphIndex> {
    let art_points = articulation_points(&graph);
    graph
        .node_indices()
        .filter_map(|node| {
            if !art_points.contains(&node) {
                Some(node.index())
            } else {
                None
            }
        })
        .collect()
}

/// Why a [`Connectivity`] cannot be used as given.
#[derive(Debug, PartialEq)]
pub enum ConnectivityError {
    /// A coupling names a qubit outside `0..num_qubits`.
    EdgeOutOfRange {
        edge: (GraphIndex, GraphIndex),
        num_qubits: usize,
    },
    /// The graph's qubits are not exactly `0..num_qubits`.
    ///
    /// Synthesis indexes tableau columns and polynomial chains *by qubit number*, so the right
    /// count is not enough: [`remove_node`](Connectivity::remove_node) leaves index gaps, and a
    /// gapped graph would have synthesis index its arrays with numbers that are out of range.
    SizeMismatch {
        nodes: Vec<GraphIndex>,
        num_qubits: usize,
    },
    /// Some qubits are unreachable from others, so a two-qubit gate between them cannot be routed.
    Disconnected {
        root: GraphIndex,
        unreachable: Vec<GraphIndex>,
    },
}

impl std::fmt::Display for ConnectivityError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ConnectivityError::EdgeOutOfRange { edge, num_qubits } => {
                write!(f, "coupling {edge:?} names a qubit outside 0..{num_qubits}")
            }
            ConnectivityError::SizeMismatch { nodes, num_qubits } => write!(
                f,
                "connectivity has qubits {nodes:?} but the circuit has {num_qubits} qubits; \
                 synthesis indexes by qubit number, so the graph must cover exactly 0..{num_qubits}"
            ),
            ConnectivityError::Disconnected { root, unreachable } => write!(
                f,
                "connectivity is disconnected: qubits {unreachable:?} cannot be reached from \
                 qubit {root}; a two-qubit gate between components cannot be routed"
            ),
        }
    }
}

impl std::error::Error for ConnectivityError {}

/// A hardware coupling graph, with shortest paths and cut vertices precomputed.
///
/// Qubits are node indices, edges are couplings, and an edge weight is the cost of routing across
/// it (1 for an unweighted graph).
///
/// # Precomputation
///
/// Construction runs Floyd–Warshall over the whole graph and finds the articulation points, so
/// [`distance`](Architecture::distance), [`best_path`](Architecture::best_path) and
/// [`non_cutting`](Architecture::non_cutting) are lookups rather than searches. Synthesis queries
/// those in its inner loops, which is what makes the trade worthwhile.
///
/// The flip side: every mutation ([`remove_node`](Self::remove_node),
/// [`add_edge`](Self::add_edge), [`add_weighted_edge`](Self::add_weighted_edge)) redoes that work
/// for the whole graph. Build the graph in one go with [`from_edges`](Self::from_edges) or
/// [`from_edges_with_size`](Self::from_edges_with_size) rather than edge by edge.
///
/// ```
/// use synir::architecture::connectivity::Connectivity;
/// use synir::architecture::Architecture;
///
/// let line = Connectivity::line(4);              // 0 - 1 - 2 - 3
/// assert_eq!(line.distance(0, 3), 3);
/// assert_eq!(line.best_path(0, 2), vec![0, 1, 2]);
///
/// // The interior qubits of a line are cut vertices; only the endpoints are removable.
/// assert_eq!(line.non_cutting(), &vec![0, 3]);
/// ```
#[derive(Debug, Default, Clone)]
pub struct Connectivity {
    graph: StableUnGraph<NodeWeight, EdgeWeight, GraphIndex>,
    non_cutting: Vec<GraphIndex>,
    prev: Vec<Vec<Option<GraphIndex>>>,
    distance: HashMap<(NodeIndex<GraphIndex>, NodeIndex<GraphIndex>), EdgeWeight>,
}

impl Connectivity {
    /// `num_qubits` qubits and no couplings at all.
    pub fn new(num_qubits: usize) -> Self {
        let mut graph = StableUnGraph::with_capacity(num_qubits, 0);
        for _ in 0..num_qubits {
            graph.add_node(());
        }
        Connectivity::from_graph(graph)
    }

    /// A path graph: `0 - 1 - ... - (num_qubits - 1)`.
    pub fn line(num_qubits: usize) -> Self {
        let edges: Vec<(usize, usize)> = (1..num_qubits).map(|i| (i - 1, i)).collect();
        Connectivity::from_edges(&edges)
    }

    /// A rectangular lattice, with qubit `r * num_cols + c` at row `r`, column `c`.
    pub fn grid(num_rows: usize, num_cols: usize) -> Self {
        let mut edges = Vec::with_capacity(2 * num_rows * num_cols);

        for r in 0..num_rows {
            for c in 0..num_cols {
                if r < num_rows - 1 {
                    edges.push((num_cols * r + c, num_cols * (r + 1) + c));
                }
                if c < num_cols - 1 {
                    edges.push((num_cols * r + c, num_cols * r + (c + 1)));
                }
            }
        }
        Connectivity::from_edges(&edges)
    }

    /// All-to-all connectivity — every qubit coupled to every other.
    ///
    /// The "no constraints" case, and what the synthesis strategies default to when no
    /// architecture is supplied.
    pub fn complete(num_qubits: usize) -> Self {
        let mut edges = Vec::with_capacity(num_qubits * num_qubits / 2);
        for i in 0..num_qubits {
            for j in (i + 1)..num_qubits {
                edges.push((i, j));
            }
        }
        Connectivity::from_edges(&edges)
    }

    /// Builds a graph from unweighted couplings, each given weight 1.
    ///
    /// Node count is inferred from the largest index mentioned, so an isolated qubit at the end of
    /// the register will not appear. An empty slice yields a single-qubit graph.
    pub fn from_edges(edges: &[(GraphIndex, GraphIndex)]) -> Self {
        if edges.len() > 0 {
            let mut graph = StableUnGraph::from_edges(edges);
            graph.edge_weights_mut().for_each(|weight| *weight = 1); // Default weight of 1 for unweighted edges
            Connectivity::from_graph(graph)
        } else {
            Connectivity::new(1)
        }
    }

    /// Builds a graph with exactly `num_qubits` qubits, coupling each listed pair with weight 1.
    ///
    /// Unlike [`from_edges`](Self::from_edges), the register width is given rather than inferred,
    /// so an isolated qubit at the end of the register survives. That is what a Qiskit `Target`
    /// hands over: a qubit count and an edge list, separately.
    ///
    /// Couplings are undirected, so `(i, j)` and `(j, i)` are the same edge and only one is kept;
    /// self-loops carry no routing information and are ignored. The precomputation runs once,
    /// however many edges there are — unlike a loop over [`add_edge`](Self::add_edge), which
    /// redoes it per edge.
    ///
    /// # Errors
    ///
    /// [`ConnectivityError::EdgeOutOfRange`] if a coupling names a qubit `>= num_qubits`.
    pub fn from_edges_with_size(
        num_qubits: usize,
        edges: &[(GraphIndex, GraphIndex)],
    ) -> Result<Self, ConnectivityError> {
        let mut graph: StableUnGraph<NodeWeight, EdgeWeight, GraphIndex> =
            StableUnGraph::with_capacity(num_qubits, edges.len());
        for _ in 0..num_qubits {
            graph.add_node(());
        }
        let mut seen = std::collections::HashSet::with_capacity(edges.len());
        for &(i, j) in edges {
            if i >= num_qubits || j >= num_qubits {
                return Err(ConnectivityError::EdgeOutOfRange {
                    edge: (i, j),
                    num_qubits,
                });
            }
            if i == j {
                continue;
            }
            // `StableUnGraph` happily stores parallel edges; a symmetrized coupling map would
            // otherwise leave `edges()` and `edge_count()` reporting each coupling twice.
            if !seen.insert((i.min(j), i.max(j))) {
                continue;
            }
            let (a, b) = (graph.from_index(i), graph.from_index(j));
            graph.add_edge(a, b, 1);
        }
        Ok(Connectivity::from_graph(graph))
    }

    /// Builds a graph from `(i, j, weight)` couplings, for hardware where routing costs differ.
    ///
    /// An empty slice yields a single-qubit graph.
    pub fn from_weighted_edges(edges: &[(GraphIndex, GraphIndex, EdgeWeight)]) -> Self {
        if edges.len() > 0 {
            let graph = StableUnGraph::from_edges(edges);
            Connectivity::from_graph(graph)
        } else {
            Connectivity::new(1)
        }
    }

    /// Wraps an existing petgraph graph, running the shortest-path and cut-vertex precomputation.
    pub fn from_graph(graph: StableUnGraph<NodeWeight, EdgeWeight, GraphIndex>) -> Self {
        let non_cutting = get_non_cutting_vertices(&graph);
        let (distance, prev) = floyd_warshall_path(&graph, |e| *e.weight()).unwrap();
        let distance = distance.into_iter().collect();

        Connectivity {
            graph,
            non_cutting,
            prev,
            distance,
        }
    }

    /// The qubits present in the graph.
    ///
    /// Not necessarily `0..node_count()`: [`remove_node`](Self::remove_node) leaves gaps, since
    /// the remaining qubits keep their original indices.
    pub fn nodes(&self) -> Vec<GraphIndex> {
        self.graph
            .node_references()
            .map(|(node, _)| node.id().index())
            .collect()
    }

    /// How many qubits are in the graph.
    pub fn node_count(&self) -> usize {
        self.graph.node_count()
    }

    /// How many couplings are in the graph.
    pub fn edge_count(&self) -> usize {
        self.graph.edge_count()
    }

    /// The couplings, as `(source, target)` pairs. The graph is undirected, so the order within
    /// each pair carries no meaning.
    pub fn edges(&self) -> Vec<(GraphIndex, GraphIndex)> {
        self.graph
            .edge_references()
            .map(|node| (node.source().index(), node.target().index()))
            .collect()
    }

    fn update(&mut self) {
        let graph = std::mem::take(&mut self.graph);
        let updated_self = Self::from_graph(graph);
        *self = updated_self;
    }

    /// Removes qubit `i` and every coupling touching it, in place.
    ///
    /// Surviving qubits keep their indices. Recomputes all shortest paths — see the type-level
    /// note on mutation cost.
    pub fn remove_node(&mut self, i: GraphIndex) {
        self.graph.remove_node(self.graph.from_index(i));
        self.update();
    }

    /// Couples `i` and `j` with weight 1. Recomputes all shortest paths.
    pub fn add_edge(&mut self, i: GraphIndex, j: GraphIndex) {
        self.add_weighted_edge(i, j, 1);
    }

    /// Couples `i` and `j` with the given weight. Recomputes all shortest paths.
    pub fn add_weighted_edge(&mut self, i: GraphIndex, j: GraphIndex, weight: EdgeWeight) {
        self.graph
            .add_edge(self.graph.from_index(i), self.graph.from_index(j), weight);
        self.update();
    }

    fn path_from_shortest_path_tree(&self, u: GraphIndex, mut v: GraphIndex) -> Vec<GraphIndex> {
        if self.prev[u][v].is_none() {
            return Vec::new();
        }

        let mut path = vec![v];
        while u != v {
            let Some(new_v) = self.prev[u][v] else {
                panic!("broken path from {u} to {v}");
            };
            v = new_v;
            path.push(v);
        }

        path.reverse();
        path
    }

    /// Whether every qubit can reach every other one.
    ///
    /// A lookup, not a search: unreachability is already recorded in the shortest-path tree that
    /// construction computes. The empty graph is vacuously connected.
    pub fn is_connected(&self) -> bool {
        let mut nodes = self.graph.node_indices();
        let Some(root) = nodes.next() else {
            return true;
        };
        let root = root.index();
        nodes.all(|node| self.prev[root][node.index()].is_some())
    }

    /// The qubits that cannot be reached from `i`.
    pub fn unreachable_from(&self, i: GraphIndex) -> Vec<GraphIndex> {
        self.graph
            .node_indices()
            .map(|node| node.index())
            .filter(|&node| node != i && self.prev[i][node].is_none())
            .collect()
    }

    /// Checks this graph can host synthesis of a `num_qubits`-qubit circuit.
    ///
    /// Both requirements come from how synthesis uses the graph: it indexes chains and tableau
    /// columns by qubit number, and it folds qubits together along paths.
    ///
    /// # Errors
    ///
    /// [`ConnectivityError::SizeMismatch`] if the qubits are not exactly `0..num_qubits`, or
    /// [`ConnectivityError::Disconnected`] if some qubit cannot reach some other.
    pub fn validate_for(&self, num_qubits: usize) -> Result<(), ConnectivityError> {
        let nodes = self.nodes();
        if nodes.len() != num_qubits || nodes.iter().enumerate().any(|(i, node)| i != *node) {
            return Err(ConnectivityError::SizeMismatch { nodes, num_qubits });
        }
        if !self.is_connected() {
            let root = nodes.first().copied().unwrap_or(0);
            return Err(ConnectivityError::Disconnected {
                root,
                unreachable: self.unreachable_from(root),
            });
        }
        Ok(())
    }

    /// Whether `v1` and `v2` are directly coupled.
    pub fn has_edge(&self, v1: GraphIndex, v2: GraphIndex) -> bool {
        self.edges().contains(&(v1, v2)) || self.edges().contains(&(v2, v1))
    }
}

impl Architecture for Connectivity {
    fn best_path(&self, i: GraphIndex, j: GraphIndex) -> Vec<GraphIndex> {
        self.path_from_shortest_path_tree(i, j)
    }

    fn distance(&self, i: GraphIndex, j: GraphIndex) -> usize {
        self.distance[&(self.graph.from_index(i), self.graph.from_index(j))]
    }

    fn neighbors(&self, i: GraphIndex) -> Vec<GraphIndex> {
        self.graph
            .neighbors(self.graph.from_index(i))
            .map(|neighbor| neighbor.index())
            .collect()
    }

    fn non_cutting(&self) -> &Vec<GraphIndex> {
        &self.non_cutting
    }

    /// Builds the ladder from a Steiner tree spanning `nodes`, traversed breadth-first from
    /// `root`.
    fn get_cx_ladder(
        &self,
        nodes: &[GraphIndex],
        root: &GraphIndex,
    ) -> Result<Vec<(GraphIndex, GraphIndex)>, LadderError> {
        let mut nodes = nodes.to_vec();
        let terminals = self
            .graph
            .node_references()
            .filter_map(|(node_index, _)| {
                // Try to remove node from `nodes`
                nodes
                    .iter()
                    .position(|x| *x == node_index.index())
                    .map(|pos| {
                        nodes.swap_remove(pos);
                        node_index
                    })
            })
            .collect_vec();

        if !nodes.is_empty() {
            return Err(LadderError::NodesNotFound(nodes));
        }

        let tree = stable_steiner_tree(&self.graph, &terminals);

        let root_node = tree
            .node_references()
            .find_map(|(item, _)| {
                if item.index() == *root {
                    Some(item)
                } else {
                    None
                }
            })
            .ok_or(LadderError::RootNotFound)?;

        let mut bfs = Bfs::new(&tree, root_node);
        let mut edge_list = Vec::new();
        let mut visited = tree.visit_map();
        visited.visit(root_node);

        while let Some(node) = bfs.next(&tree) {
            for neighbor in tree.neighbors(node) {
                if !visited.is_visited(&neighbor) {
                    visited.visit(neighbor);
                    edge_list.push((node.index(), neighbor.index()));
                }
            }
        }
        Ok(edge_list)
    }

    fn disconnect(&self, i: GraphIndex) -> Connectivity {
        let mut graph = self.graph.clone();
        graph.remove_node(graph.from_index(i));
        Connectivity::from_graph(graph)
    }
}

#[cfg(test)]
mod tests {
    use crate::architecture::{Architecture, EdgeWeight, GraphIndex, LadderError};

    use super::{Connectivity, ConnectivityError};
    fn setup_weighted() -> Vec<(GraphIndex, GraphIndex, EdgeWeight)> {
        vec![
            (0, 1, 7),
            (0, 5, 6),
            (1, 2, 1),
            (1, 5, 7),
            (2, 3, 1),
            (2, 4, 3),
            (3, 4, 1),
            (3, 5, 4),
            (4, 5, 10),
        ]
    }

    fn setup_simple() -> Vec<(GraphIndex, GraphIndex)> {
        vec![
            (0, 1),
            (0, 5),
            (1, 2),
            (1, 5),
            (2, 3),
            (2, 4),
            (3, 4),
            (3, 5),
            (4, 5),
        ]
    }

    #[test]
    fn test_line_creation() {
        let line_architecture = Connectivity::line(5);

        assert_eq!(
            line_architecture.edges(),
            vec![(0, 1), (1, 2), (2, 3), (3, 4)]
        );
    }

    #[test]
    fn test_grid_creation() {
        let grid_architecture = Connectivity::grid(3, 3);
        let mut edges = grid_architecture.edges();
        edges.sort();
        assert_eq!(
            edges,
            vec![
                (0, 1),
                (0, 3),
                (1, 2),
                (1, 4),
                (2, 5),
                (3, 4),
                (3, 6),
                (4, 5),
                (4, 7),
                (5, 8),
                (6, 7),
                (7, 8)
            ]
        );
    }

    #[test]
    fn test_root_is_not_present() {
        let new_architecture = Connectivity::from_edges(&setup_simple());
        assert_eq!(
            new_architecture
                .get_cx_ladder(&[1, 2, 3, 4], &42)
                .expect_err("Should return a Error that the root was not found"),
            LadderError::RootNotFound
        );
    }

    #[test]
    fn test_cx_ladder_error() {
        let new_architecture = Connectivity::from_edges(&setup_simple());
        assert_eq!(
            new_architecture
                .get_cx_ladder(&[6, 1, 2, 3], &0)
                .expect_err(
                    "Should return an error that the nodes are not part of the architecture"
                ),
            LadderError::NodesNotFound(vec![6])
        );
    }

    #[test]
    fn test_simple_constructor() {
        let new_architecture = Connectivity::from_edges(&setup_simple());
        assert_eq!(new_architecture.nodes(), vec![0, 1, 2, 3, 4, 5]);
    }

    #[test]
    fn test_cx_ladder_line_setup() {
        let new_architecture = Connectivity::from_edges(&setup_simple());
        assert_eq!(
            new_architecture.get_cx_ladder(&[0, 1, 2, 3], &0).unwrap(),
            vec![(0, 1), (1, 2), (2, 3)]
        );
    }

    #[test]
    fn test_cx_ladder_extended_triangle() {
        let new_architecture = Connectivity::from_edges(&setup_simple());
        assert_eq!(
            new_architecture
                .get_cx_ladder(&[0, 1, 2, 4, 5], &1)
                .unwrap()
                .len(),
            4
        );
    }

    #[test]
    fn test_cx_ladder_weighted_extended_triangle() {
        let new_architecture = Connectivity::from_weighted_edges(&setup_weighted());
        assert_eq!(
            new_architecture
                .get_cx_ladder(&[0, 1, 2, 4, 5], &1)
                .unwrap(),
            vec![(1, 2), (2, 3), (3, 5), (3, 4), (5, 0)]
        );
    }

    #[test]
    fn test_cx_ladder_small_triangle() {
        let new_architecture = Connectivity::from_edges(&setup_simple());
        assert_eq!(
            new_architecture
                .get_cx_ladder(&[2, 3, 4], &2)
                .unwrap()
                .len(),
            2
        );
        assert_eq!(
            new_architecture
                .get_cx_ladder(&[2, 3, 4], &4)
                .unwrap()
                .len(),
            2
        );
    }

    #[test]
    fn test_cx_ladder_weighted_small_triangle() {
        let new_architecture = Connectivity::from_weighted_edges(&setup_weighted());
        assert_eq!(
            new_architecture.get_cx_ladder(&[2, 3, 4], &2).unwrap(),
            vec![(2, 3), (3, 4)]
        );
        assert_eq!(
            new_architecture.get_cx_ladder(&[2, 3, 4], &4).unwrap(),
            vec![(4, 3), (3, 2)]
        );
    }

    #[test]
    fn test_weighted_constructor() {
        let new_architecture = Connectivity::from_weighted_edges(&setup_weighted());
        assert_eq!(new_architecture.nodes(), vec![0, 1, 2, 3, 4, 5]);
    }

    #[test]
    fn test_best_simple_path() {
        let new_architecture = Connectivity::from_edges(&setup_simple());

        assert_eq!(new_architecture.best_path(0, 4), vec![0, 5, 4]);
    }

    #[test]
    fn test_best_weighted_path() {
        let new_architecture = Connectivity::from_weighted_edges(&setup_weighted());

        assert_eq!(new_architecture.best_path(0, 4), vec![0, 1, 2, 3, 4]);
    }

    #[test]
    #[should_panic = "index out of bounds: the len is 6 but the index is 6"]
    fn test_best_path_missing() {
        let new_architecture = Connectivity::from_edges(&setup_simple());
        new_architecture.best_path(5, 6);
    }

    #[test]
    fn test_distance() {
        let new_architecture = Connectivity::from_weighted_edges(&setup_weighted());
        assert_eq!(new_architecture.distance(2, 4), 2);
        assert_eq!(new_architecture.distance(4, 2), 2);
        assert_eq!(new_architecture.distance(0, 4), 10);
    }

    #[test]
    fn test_simple_distance() {
        let new_architecture = Connectivity::from_edges(&setup_simple());
        assert_eq!(1, new_architecture.distance(0, 1));
        assert_eq!(2, new_architecture.distance(1, 4));
    }

    #[test]
    #[should_panic = "no entry found for key"]
    fn test_distance_missing() {
        let new_architecture = Connectivity::from_edges(&setup_simple());
        new_architecture.distance(5, 6);
    }

    #[test]
    fn test_neighbors() {
        let new_architecture = Connectivity::from_edges(&setup_simple());
        assert_eq!(new_architecture.neighbors(2), vec![4, 3, 1]);
    }

    #[test]
    #[should_panic = "no entry found for key"]
    fn test_neighbor_missing() {
        let new_architecture = Connectivity::from_edges(&setup_simple());
        new_architecture.distance(2, 7);
    }

    #[test]
    fn test_non_cutting() {
        let new_architecture = Connectivity::from_edges(&setup_simple());
        assert_eq!(&new_architecture.nodes(), new_architecture.non_cutting());
    }

    #[test]
    fn test_non_cutting_line() {
        let line_architecture = Connectivity::line(5);
        assert_eq!(*line_architecture.non_cutting(), vec![0, 4]);
    }

    #[test]
    fn test_non_cutting_grid() {
        let line_architecture = Connectivity::grid(3, 3);
        assert_eq!(
            *line_architecture.non_cutting(),
            vec![0, 1, 2, 3, 4, 5, 6, 7, 8]
        );
    }

    #[test]
    fn test_non_cutting_complete() {
        let line_architecture = Connectivity::complete(3);
        assert_eq!(*line_architecture.non_cutting(), vec![0, 1, 2]);
    }

    #[test]
    fn test_remove_node() {
        let mut architecture = Connectivity::from_edges(&setup_simple());
        assert_eq!(architecture.nodes(), vec![0, 1, 2, 3, 4, 5]);

        assert_eq!(
            architecture.edges(),
            vec![
                (0, 1),
                (0, 5),
                (1, 2),
                (1, 5),
                (2, 3),
                (2, 4),
                (3, 4),
                (3, 5),
                (4, 5),
            ]
        );

        architecture.remove_node(1);

        assert_eq!(architecture.nodes(), vec![0, 2, 3, 4, 5]);

        assert_eq!(
            architecture.edges(),
            vec![(0, 5), (2, 3), (2, 4), (3, 4), (3, 5), (4, 5)]
        );
    }

    #[test]
    fn test_remove_node_line() {
        let mut architecture = Connectivity::line(5);
        assert_eq!(architecture.nodes(), vec![0, 1, 2, 3, 4]);

        assert_eq!(architecture.edges(), vec![(0, 1), (1, 2), (2, 3), (3, 4)]);

        architecture.remove_node(1);

        assert_eq!(architecture.nodes(), vec![0, 2, 3, 4]);

        assert_eq!(architecture.edges(), vec![(2, 3), (3, 4)]);
    }

    #[test]
    fn test_remove_node_grid() {
        let mut architecture = Connectivity::grid(3, 3);
        assert_eq!(architecture.nodes(), vec![0, 1, 2, 3, 4, 5, 6, 7, 8]);

        assert_eq!(
            architecture.edges(),
            vec![
                (0, 3),
                (0, 1),
                (1, 4),
                (1, 2),
                (2, 5),
                (3, 6),
                (3, 4),
                (4, 7),
                (4, 5),
                (5, 8),
                (6, 7),
                (7, 8)
            ]
        );

        architecture.remove_node(1);

        assert_eq!(architecture.nodes(), vec![0, 2, 3, 4, 5, 6, 7, 8]);

        assert_eq!(
            architecture.edges(),
            vec![
                (0, 3),
                (2, 5),
                (3, 6),
                (3, 4),
                (4, 7),
                (4, 5),
                (5, 8),
                (6, 7),
                (7, 8),
            ]
        );
    }

    #[test]
    fn test_disconnect() {
        let architecture = Connectivity::from_edges(&setup_simple());
        let new_architecture = architecture.disconnect(1);

        assert_eq!(architecture.nodes(), vec![0, 1, 2, 3, 4, 5]);

        assert_eq!(
            architecture.edges(),
            vec![
                (0, 1),
                (0, 5),
                (1, 2),
                (1, 5),
                (2, 3),
                (2, 4),
                (3, 4),
                (3, 5),
                (4, 5),
            ]
        );

        assert_eq!(new_architecture.nodes(), vec![0, 2, 3, 4, 5]);

        assert_eq!(
            new_architecture.edges(),
            vec![(0, 5), (2, 3), (2, 4), (3, 4), (3, 5), (4, 5)]
        );
    }

    #[test]
    fn test_disconnect_line() {
        let architecture = Connectivity::line(5);
        let new_architecture = architecture.disconnect(1);

        assert_eq!(architecture.nodes(), vec![0, 1, 2, 3, 4]);

        assert_eq!(architecture.edges(), vec![(0, 1), (1, 2), (2, 3), (3, 4)]);

        assert_eq!(new_architecture.nodes(), vec![0, 2, 3, 4]);

        assert_eq!(new_architecture.edges(), vec![(2, 3), (3, 4)]);
    }

    #[test]
    fn test_disconnect_grid() {
        let architecture = Connectivity::grid(3, 3);
        let new_architecture = architecture.disconnect(1);

        assert_eq!(architecture.nodes(), vec![0, 1, 2, 3, 4, 5, 6, 7, 8]);

        assert_eq!(
            architecture.edges(),
            vec![
                (0, 3),
                (0, 1),
                (1, 4),
                (1, 2),
                (2, 5),
                (3, 6),
                (3, 4),
                (4, 7),
                (4, 5),
                (5, 8),
                (6, 7),
                (7, 8)
            ]
        );

        assert_eq!(new_architecture.nodes(), vec![0, 2, 3, 4, 5, 6, 7, 8]);

        assert_eq!(
            new_architecture.edges(),
            vec![
                (0, 3),
                (2, 5),
                (3, 6),
                (3, 4),
                (4, 7),
                (4, 5),
                (5, 8),
                (6, 7),
                (7, 8)
            ]
        );
    }

    #[test]
    fn test_disconnected_cx_ladder() {
        let architecture = Connectivity::from_weighted_edges(&setup_weighted());
        let new_architecture = architecture.disconnect(3);

        assert_eq!(new_architecture.nodes(), vec![0, 1, 2, 4, 5]);

        assert_eq!(
            new_architecture.edges(),
            vec![(0, 1), (0, 5), (1, 2), (1, 5), (2, 4), (4, 5)]
        );

        assert_eq!(
            new_architecture.get_cx_ladder(&[1, 2, 4, 5], &2).unwrap(),
            vec![(2, 4), (2, 1), (1, 5)]
        );

        let new_architecture = new_architecture.disconnect(2);

        assert_eq!(new_architecture.nodes(), vec![0, 1, 4, 5]);

        assert_eq!(
            new_architecture.edges(),
            vec![(0, 1), (0, 5), (1, 5), (4, 5)]
        );

        assert_eq!(
            new_architecture.get_cx_ladder(&[1, 4, 5], &1).unwrap(),
            vec![(1, 5), (5, 4)]
        );
    }

    #[test]
    fn test_from_edges_with_size_keeps_isolated_tail() {
        // `from_edges` would infer three qubits here and lose qubit 3.
        let architecture = Connectivity::from_edges_with_size(4, &[(0, 1), (1, 2)]).unwrap();
        assert_eq!(architecture.nodes(), vec![0, 1, 2, 3]);
        assert_eq!(architecture.edge_count(), 2);
        assert_eq!(Connectivity::from_edges(&[(0, 1), (1, 2)]).node_count(), 3);
    }

    #[test]
    fn test_from_edges_with_size_matches_line() {
        let architecture =
            Connectivity::from_edges_with_size(4, &[(0, 1), (1, 2), (2, 3)]).unwrap();
        assert_eq!(architecture.edges(), Connectivity::line(4).edges());
    }

    #[test]
    fn test_from_edges_with_size_rejects_out_of_range() {
        assert_eq!(
            Connectivity::from_edges_with_size(2, &[(0, 5)]).unwrap_err(),
            ConnectivityError::EdgeOutOfRange {
                edge: (0, 5),
                num_qubits: 2
            }
        );
    }

    #[test]
    fn test_from_edges_with_size_ignores_self_loops_and_duplicates() {
        // A symmetrized Qiskit coupling map hands over both directions of every pair.
        let architecture =
            Connectivity::from_edges_with_size(3, &[(1, 1), (0, 1), (1, 0), (1, 2)]).unwrap();
        assert_eq!(architecture.edge_count(), 2);
        assert_eq!(architecture.nodes(), vec![0, 1, 2]);
    }

    #[test]
    fn test_is_connected() {
        assert!(Connectivity::line(5).is_connected());
        assert!(Connectivity::grid(2, 3).is_connected());
        assert!(Connectivity::complete(3).is_connected());
        assert!(Connectivity::new(1).is_connected());
        assert!(Connectivity::default().is_connected()); // vacuous: no qubits
        assert!(!Connectivity::new(2).is_connected()); // two isolated qubits
        assert!(!Connectivity::from_edges_with_size(4, &[(0, 1), (2, 3)])
            .unwrap()
            .is_connected());
        assert!(!Connectivity::from_edges_with_size(3, &[(0, 1)])
            .unwrap()
            .is_connected());
    }

    #[test]
    fn test_unreachable_from() {
        let architecture = Connectivity::from_edges_with_size(4, &[(0, 1), (2, 3)]).unwrap();
        assert_eq!(architecture.unreachable_from(0), vec![2, 3]);
        assert_eq!(architecture.unreachable_from(2), vec![0, 1]);
        assert!(Connectivity::line(4).unreachable_from(0).is_empty());
    }

    #[test]
    fn test_is_connected_after_remove_node() {
        let mut architecture = Connectivity::line(3);
        architecture.remove_node(1);
        assert!(!architecture.is_connected());
        assert_eq!(architecture.unreachable_from(0), vec![2]);
    }

    #[test]
    fn test_validate_for() {
        assert!(Connectivity::line(4).validate_for(4).is_ok());
        assert!(Connectivity::complete(1).validate_for(1).is_ok());
        assert!(matches!(
            Connectivity::line(4).validate_for(5),
            Err(ConnectivityError::SizeMismatch { .. })
        ));
        assert!(matches!(
            Connectivity::line(4).validate_for(3),
            Err(ConnectivityError::SizeMismatch { .. })
        ));
        assert!(matches!(
            Connectivity::from_edges_with_size(4, &[(0, 1), (2, 3)])
                .unwrap()
                .validate_for(4),
            Err(ConnectivityError::Disconnected { .. })
        ));
    }

    #[test]
    fn test_validate_for_rejects_gapped_indices() {
        // The right *count* is not enough: `remove_node` leaves the surviving qubits with their
        // original numbers, and synthesis would index its arrays with the out-of-range ones.
        let mut architecture = Connectivity::complete(4);
        architecture.remove_node(1);
        assert_eq!(architecture.node_count(), 3);
        assert!(matches!(
            architecture.validate_for(3),
            Err(ConnectivityError::SizeMismatch { .. })
        ));
    }

    #[test]
    fn test_connectivity_error_messages_name_the_problem() {
        let message = Connectivity::line(4)
            .validate_for(3)
            .unwrap_err()
            .to_string();
        assert!(message.contains("the circuit has 3 qubits"), "{message}");

        let message = Connectivity::from_edges_with_size(4, &[(0, 1), (2, 3)])
            .unwrap()
            .validate_for(4)
            .unwrap_err()
            .to_string();
        assert!(message.contains("disconnected"), "{message}");
    }
}
