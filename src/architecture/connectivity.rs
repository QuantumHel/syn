use super::{Architecture, EdgeWeight, GraphIndex, LadderError, NodeWeight};
use petgraph::algo::floyd_warshall::floyd_warshall_path;
use petgraph::algo::steiner_tree::steiner_tree;
use petgraph::prelude::EdgeRef;
use petgraph::visit::{Bfs, Dfs, GraphBase, IntoNeighbors, VisitMap, Visitable, Walker};
use petgraph::{
    algo::articulation_points::articulation_points,
    graph::{NodeIndex, UnGraph},
    visit::{IntoNodeReferences, NodeIndexable, NodeRef},
};
use std::borrow::Borrow;
use std::collections::HashMap;
use std::ops::Index;

/// Get all the vertices in a graph that are non-cutting (won't make the graph disconnected)
fn get_non_cutting_vertices(
    graph: &UnGraph<NodeWeight, EdgeWeight, GraphIndex>,
) -> Vec<GraphIndex> {
    let art_points = articulation_points(&graph);
    (0..graph.node_count())
        .filter(|node| !art_points.contains(&graph.from_index(*node)))
        .collect()
}

#[derive(Debug)]
pub struct Connectivity {
    graph: UnGraph<NodeWeight, EdgeWeight, GraphIndex>,
    non_cutting: Vec<GraphIndex>,
    prev: Vec<Vec<Option<GraphIndex>>>,
    distance: HashMap<(NodeIndex<GraphIndex>, NodeIndex<GraphIndex>), EdgeWeight>,
}

impl Connectivity {
    pub fn new(num_qubits: usize) -> Self {
        Connectivity {
            graph: UnGraph::with_capacity(num_qubits, 0),
            non_cutting: Default::default(),
            prev: Default::default(),
            distance: HashMap::new(),
        }
    }

    pub fn line(nr_qubits: usize) -> Self {
        let edges: Vec<(usize, usize)> = (0..nr_qubits - 1).map(|i| (i, i + 1)).collect();
        println!("{:?}", edges);
        Connectivity::from_edges(&edges)
    }

    pub fn grid(num_rows: usize, num_cols: usize) -> Self {
        let mut edges = Vec::new();

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
        let mut edges = Vec::new();
        for i in 0..num_qubits {
            for j in (i + 1)..num_qubits {
                edges.push((i, j));
            }
        }
        Connectivity::from_edges(&edges)
    }

    pub fn from_edges(edges: &[(GraphIndex, GraphIndex)]) -> Self {
        let graph = UnGraph::from_edges(edges);
        Connectivity::from_graph(graph)
    }

    pub fn from_weighted_edges(edges: &[(GraphIndex, GraphIndex, EdgeWeight)]) -> Self {
        let graph = UnGraph::from_edges(edges);
        Connectivity::from_graph(graph)
    }

    pub fn from_graph(graph: UnGraph<NodeWeight, EdgeWeight, GraphIndex>) -> Self {
        let non_cutting = get_non_cutting_vertices(&graph);
        let (distance, prev) = floyd_warshall_path(&graph, |e| *e.weight()).unwrap();
        let distance = distance.iter().map(|(k, v)| (*k, *v as usize)).collect();
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
            .map(|node| self.graph.to_index(node.id()))
            .collect()
    }

    pub fn edges(&self) -> Vec<(GraphIndex, GraphIndex)> {
        let graph_edges: Vec<(GraphIndex, GraphIndex)> = self
            .graph
            .edge_references()
            .map(|node| {
                (
                    self.graph.to_index(node.source()),
                    self.graph.to_index(node.target()),
                )
            })
            .collect();
        graph_edges
    }

    fn update(&mut self) {
        let non_cutting = get_non_cutting_vertices(&self.graph);

        let (distance, prev) = floyd_warshall_path(&self.graph, |e| *e.weight()).unwrap();

        let distance = distance.iter().map(|(k, v)| (*k, *v as usize)).collect();

        self.non_cutting = non_cutting;
        self.distance = distance;
        self.prev = prev;
    }

    pub fn remove_node(&mut self, i: GraphIndex) {
        self.graph.remove_node(self.graph.from_index(i));
        self.update();
    }

    pub fn add_node(&mut self) {
        self.graph.add_node(());
    }

    pub fn add_edge(&mut self, i: GraphIndex, j: GraphIndex) {
        self.graph
            .add_edge(self.graph.from_index(i), self.graph.from_index(j), 1);
        self.update();
    }

    pub fn add_weighted_edge(&mut self, i: GraphIndex, j: GraphIndex, weight: EdgeWeight) {
        self.graph
            .add_edge(self.graph.from_index(i), self.graph.from_index(j), weight);
        self.update();
    }

    fn path_from_shortest_path_tree(&self, u: GraphIndex, mut v: GraphIndex) -> Vec<GraphIndex> {
        let mut path = vec![v];

        if self.prev[u][v].is_none() {
            return Vec::new();
        }

        while u != v {
            if let Some(new_v) = self.prev[u][v] {
                v = new_v;
                path.push(v);
            }
        }

        path.reverse();
        path
    }
}

impl Architecture for Connectivity {
    fn best_path(&self, i: GraphIndex, j: GraphIndex) -> Vec<GraphIndex> {
        assert!(
            i < self.graph.node_count(),
            "architecture does not contain node {i}"
        );
        assert!(
            j < self.graph.node_count(),
            "architecture does not contain node {j}"
        );
        self.path_from_shortest_path_tree(i, j)
    }

    fn distance(&self, i: GraphIndex, j: GraphIndex) -> usize {
        assert!(
            i < self.graph.node_count(),
            "architecture does not contain node {i}"
        );
        assert!(
            j < self.graph.node_count(),
            "architecture does not contain node {j}"
        );
        self.distance[&(self.graph.from_index(i), self.graph.from_index(j))] as usize
    }

    fn neighbors(&self, i: GraphIndex) -> Vec<usize> {
        assert!(
            i < self.graph.node_count(),
            "architecture does not contain node {i}"
        );
        self.graph
            .neighbors(NodeIndex::new(i))
            .map(|neighbor| neighbor.index())
            .collect()
    }

    fn non_cutting(&mut self) -> &Vec<GraphIndex> {
        &self.non_cutting
    }

    fn get_cx_ladder(
        &self,
        nodes: &[GraphIndex],
        root: &GraphIndex,
    ) -> Result<Vec<(usize, usize)>, LadderError> {
        // TODO fix me we currently have to convert the graph
        let terminals: Result<Vec<_>, LadderError> = self
            .graph
            .node_indices()
            .filter(|node_index| nodes.contains(&self.graph.to_index(*node_index)))
            .map(|node_index| {
                let idx = node_index.index() as u32;
                idx.try_into()
                    .map(NodeIndex::new)
                    .map_err(|_| LadderError::ConversionError)
            })
            .collect();
        let terminals = terminals?;

        let tree = steiner_tree(&self.graph, &terminals);

        let root_node = tree
            .node_indices()
            .find(|item| item.index() == *root as usize)
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
            (1, 5, 5),
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
        let mut line_architecture = Connectivity::line(5);
        assert_eq!(
            line_architecture.edges(),
            vec![(0, 1), (1, 2), (2, 3), (3, 4)]
        );
    }

    #[test]
    fn test_grid_creation() {
        let mut line_architecture = Connectivity::grid(3, 3);
        let mut edges = line_architecture.edges();
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
    fn test_weight_is_considered() {
        let new_architecture = Connectivity::from_weighted_edges(&setup_weighted());
        assert_eq!(
            new_architecture
                .get_cx_ladder(&vec![1, 2, 3, 4], &2)
                .unwrap()
                .sort(),
            vec![(2, 1), (2, 3), (3, 4)].sort()
        );
    }

    #[test]
    fn test_root_is_not_present() {
        let new_architecture = Connectivity::from_edges(&setup_simple());
        assert_eq!(
            new_architecture
                .get_cx_ladder(&vec![1, 2, 3, 4], &42)
                .expect_err("Should return a Error that the root was not found"),
            LadderError::RootNotFound
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
            new_architecture
                .get_cx_ladder(&vec![0, 1, 2, 3], &0)
                .unwrap(),
            vec![(0, 1), (1, 2), (2, 3)]
        );
    }

    #[test]
    fn test_cx_ladder_extended_triangle() {
        let new_architecture = Connectivity::from_edges(&setup_simple());
        assert_eq!(
            new_architecture
                .get_cx_ladder(&vec![0, 1, 2, 4, 5], &1)
                .unwrap()
                .sort(),
            vec![(0, 1), (1, 5), (1, 2), (2, 4)].sort()
        );
    }

    #[test]
    fn test_cx_ladder_small_triangle() {
        let new_architecture = Connectivity::from_edges(&setup_simple());
        assert_eq!(
            new_architecture
                .get_cx_ladder(&vec![2, 3, 4], &2)
                .unwrap()
                .sort(),
            vec![(2, 4), (2, 3)].sort()
        );
        assert_eq!(
            new_architecture
                .get_cx_ladder(&vec![2, 3, 4], &4)
                .unwrap()
                .sort(),
            vec![(4, 2), (4, 3)].sort()
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

        assert_eq!(vec![0, 1, 2, 4], new_architecture.best_path(0, 4));
    }

    #[test]
    fn test_best_weighted_path() {
        let new_architecture = Connectivity::from_weighted_edges(&setup_weighted());

        assert_eq!(vec![0, 1, 2, 3, 4], new_architecture.best_path(0, 4));
    }

    #[test]
    #[should_panic = "architecture does not contain node 6"]
    fn test_best_path_missing() {
        let new_architecture = Connectivity::from_edges(&setup_simple());
        new_architecture.best_path(5, 6);
    }

    #[test]
    fn test_distance() {
        let new_architecture = Connectivity::from_weighted_edges(&setup_weighted());
        assert_eq!(2, new_architecture.distance(2, 4));
        assert_eq!(2, new_architecture.distance(4, 2));
        assert_eq!(10, new_architecture.distance(0, 4));
    }

    #[test]
    #[should_panic = "architecture does not contain node 6"]
    fn test_distance_missing() {
        let new_architecture = Connectivity::from_edges(&setup_simple());
        new_architecture.distance(5, 6);
    }

    #[test]
    fn test_neighbors() {
        let new_architecture = Connectivity::from_edges(&setup_simple());
        assert_eq!(vec![4, 3, 1], new_architecture.neighbors(2));
    }

    #[test]
    #[should_panic = "architecture does not contain node 7"]
    fn test_neighbor_missing() {
        let new_architecture = Connectivity::from_edges(&setup_simple());
        new_architecture.distance(2, 7);
    }

    #[test]
    fn test_non_cutting() {
        let mut new_architecture = Connectivity::from_edges(&setup_simple());
        assert_eq!(&new_architecture.nodes(), new_architecture.non_cutting());
    }

    #[test]
    fn test_non_cutting_line() {
        let mut line_architecture = Connectivity::line(5);
        assert_eq!(*line_architecture.non_cutting(), vec![0, 4]);
    }
}
