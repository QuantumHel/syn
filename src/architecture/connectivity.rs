use super::{Architecture, EdgeWeight, GraphIndex, NodeWeight};
use petgraph::algo::floyd_warshall::floyd_warshall_path;
use petgraph::{
    algo::articulation_points::articulation_points,
    graph::{NodeIndex, UnGraph},
    visit::{IntoNodeReferences, NodeIndexable, NodeRef},
};
use std::collections::HashMap;

#[derive(Debug)]
pub struct Connectivity {
    graph: UnGraph<NodeWeight, EdgeWeight, GraphIndex>,
    non_cutting: Vec<GraphIndex>,
    prev: Vec<Vec<Option<GraphIndex>>>,
    distance: HashMap<(NodeIndex<GraphIndex>, NodeIndex<GraphIndex>), usize>,
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

    pub fn from_edges(edges: &[(GraphIndex, GraphIndex)]) -> Self {
        let graph = UnGraph::from_edges(edges);
        let art_points = articulation_points(&graph);

        let non_cutting = (0..graph.node_count())
            .filter(|node| art_points.contains(&graph.from_index(*node)))
            .collect();

        let (distance, prev) = floyd_warshall_path(&graph, |e| *e.weight()).unwrap();

        Connectivity {
            graph,
            non_cutting,
            prev,
            distance,
        }
    }

    pub fn from_weighted_edges(edges: &[(GraphIndex, GraphIndex, EdgeWeight)]) -> Self {
        let graph = UnGraph::from_edges(edges);
        let art_points = articulation_points(&graph);

        let non_cutting = (0..graph.node_count())
            .filter(|node| art_points.contains(&graph.from_index(*node)))
            .collect();

        let (distance, prev) = floyd_warshall_path(&graph, |e| *e.weight()).unwrap();

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

    fn update(&mut self) {
        let art_points = articulation_points(&self.graph);

        let non_cutting = (0..self.graph.node_count())
            .filter(|node| art_points.contains(&self.graph.from_index(*node)))
            .collect();

        let (distance, prev) = floyd_warshall_path(&self.graph, |e| *e.weight()).unwrap();
        println!("prev: {:?}", prev);
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
        self.distance[&(self.graph.from_index(i), self.graph.from_index(j))]
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
}

#[cfg(test)]
mod tests {
    use crate::architecture::{Architecture, EdgeWeight, GraphIndex};

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
    fn test_simple_constuctor() {
        let new_architecture = Connectivity::from_edges(&setup_simple());
        assert_eq!(new_architecture.nodes(), vec![0, 1, 2, 3, 4, 5]);
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
        assert_eq!(&Vec::<GraphIndex>::new(), new_architecture.non_cutting());
    }
}
