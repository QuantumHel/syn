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

#[derive(Debug, Default)]
pub struct Connectivity {
    graph: StableUnGraph<NodeWeight, EdgeWeight, GraphIndex>,
    non_cutting: Vec<GraphIndex>,
    prev: Vec<Vec<Option<GraphIndex>>>,
    distance: HashMap<(NodeIndex<GraphIndex>, NodeIndex<GraphIndex>), EdgeWeight>,
}

impl Connectivity {
    pub fn new(num_qubits: usize) -> Self {
        let mut graph = StableUnGraph::with_capacity(num_qubits, 0);
        for _ in 0..num_qubits {
            graph.add_node(());
        }
        Connectivity::from_graph(graph)
    }

    pub fn line(num_qubits: usize) -> Self {
        let edges: Vec<(usize, usize)> = (1..num_qubits).map(|i| (i - 1, i)).collect();
        Connectivity::from_edges(&edges)
    }

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

    pub fn complete(num_qubits: usize) -> Self {
        let mut edges = Vec::with_capacity(num_qubits * num_qubits / 2);
        for i in 0..num_qubits {
            for j in (i + 1)..num_qubits {
                edges.push((i, j));
            }
        }
        Connectivity::from_edges(&edges)
    }

    pub fn from_edges(edges: &[(GraphIndex, GraphIndex)]) -> Self {
        if edges.len() > 0 {
            let mut graph = StableUnGraph::from_edges(edges);
            graph.edge_weights_mut().for_each(|weight| *weight = 1); // Default weight of 1 for unweighted edges
            Connectivity::from_graph(graph)
        } else {
            Connectivity::new(1)
        }
    }

    pub fn from_weighted_edges(edges: &[(GraphIndex, GraphIndex, EdgeWeight)]) -> Self {
        if edges.len() > 0 {
            let graph = StableUnGraph::from_edges(edges);
            Connectivity::from_graph(graph)
        } else {
            Connectivity::new(1)
        }
    }

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

    pub fn nodes(&self) -> Vec<GraphIndex> {
        self.graph
            .node_references()
            .map(|(node, _)| node.id().index())
            .collect()
    }

    pub fn node_count(&self) -> usize {
        self.graph.node_count()
    }

    pub fn edge_count(&self) -> usize {
        self.graph.edge_count()
    }

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

    pub fn remove_node(&mut self, i: GraphIndex) {
        self.graph.remove_node(self.graph.from_index(i));
        self.update();
    }

    pub fn add_edge(&mut self, i: GraphIndex, j: GraphIndex) {
        self.add_weighted_edge(i, j, 1);
    }

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

    /// Obtain cx ladder that is architecture conforming that is rooted at `root`
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

    use super::Connectivity;
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
}
