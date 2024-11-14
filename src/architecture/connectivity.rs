use petgraph::{
    algo::{articulation_points::articulation_points, floyd_warshall_path},
    graph::{NodeIndex, UnGraph},
    visit::NodeIndexable,
};
use std::collections::{HashMap, HashSet};

use super::{Architecture, EdgeWeight, GraphIndex, NodeWeight};

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
        Connectivity {
            graph: UnGraph::from_edges(edges),
            non_cutting: Default::default(),
            prev: Default::default(),
            distance: HashMap::new(),
        }
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
            i <= self.graph.node_count(),
            "architecture does not contain node {i}"
        );
        assert!(
            j <= self.graph.node_count(),
            "architecture does not contain node {j}"
        );
        self.path_from_shortest_path_tree(i, j)
    }

    fn distance(&self, i: GraphIndex, j: GraphIndex) -> usize {
        assert!(
            i <= self.graph.node_count(),
            "architecture does not contain node {i}"
        );
        assert!(
            j <= self.graph.node_count(),
            "architecture does not contain node {j}"
        );
        self.distance[&(self.graph.from_index(i), self.graph.from_index(j))]
    }

    fn neighbors(&self, i: GraphIndex) -> Vec<usize> {
        assert!(
            i <= self.graph.node_count(),
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
