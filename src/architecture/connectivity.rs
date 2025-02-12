use petgraph::{
    algo::{articulation_points::articulation_points, floyd_warshall::floyd_warshall_path}, graph::{Graph, NodeIndex}, visit::{IntoNodeReferences, NodeIndexable, NodeRef}, Directed, Direction, EdgeType, Undirected
};
use std::collections::HashMap;

use super::{Architecture, EdgeWeight, GraphIndex, NodeWeight};

#[derive(Debug)]
pub struct Connectivity<T: EdgeType> {
    graph: Graph<NodeWeight, EdgeWeight, T, GraphIndex>,
    non_cutting: Vec<GraphIndex>,
    prev: Vec<Vec<Option<GraphIndex>>>,
    distance: HashMap<(NodeIndex<GraphIndex>, NodeIndex<GraphIndex>), usize>,

}

type DirectedConnectivity = Connectivity<Directed>;
type UndirectedConnectivity = Connectivity<Undirected>;

impl DirectedConnectivity {

    pub fn add_node(&mut self, connections: Vec<(GraphIndex, Direction)>) -> GraphIndex {
        self.add_node_weighted_edges(connections.iter().map(|(i, d)| (*i, *d, 1 as EdgeWeight)).collect())
    }

    pub fn add_node_weighted_edges(&mut self, connections: Vec<(GraphIndex, Direction, EdgeWeight)>) -> GraphIndex {
        let new_node = self.graph.add_node(());
        for (other, direction, weight) in connections{
            match direction{
                Direction::Outgoing => self.add_weighted_edge(new_node.index(), other, weight),
                Direction::Incoming => self.add_weighted_edge(other, new_node.index(), weight),
            }
        }
        self.update();
        new_node.index()
    }
}

impl UndirectedConnectivity {

    pub fn add_node(&mut self, connections: Vec<GraphIndex>) -> GraphIndex {
        self.add_node_weighted_edges(connections.iter().map(|i| (*i, 1 as EdgeWeight)).collect())
    }

    pub fn add_node_weighted_edges(&mut self, connections: Vec<(GraphIndex, EdgeWeight)>) -> GraphIndex {
        let new_node = self.graph.add_node(());
        for (other, weight) in connections{
            self.add_weighted_edge(new_node.index(), other, weight);
        }
        self.update();
        new_node.index()
    }
}

impl<T: EdgeType> Connectivity<T>{

    pub fn new() -> Self {
        Connectivity {
            graph: Graph::with_capacity(1, 0),
            non_cutting: Default::default(),
            prev: Default::default(),
            distance: HashMap::new(),

        }
    }

    pub fn from_edges(edges: &[(GraphIndex, GraphIndex)]) -> Self {
        let mut connectivity = Connectivity {
            graph: Graph::from_edges(edges),
            non_cutting: Default::default(),
            prev: Default::default(),
            distance: HashMap::new(),
        };
        connectivity.update();
        connectivity
    }

    pub fn from_weighted_edges(edges: &[(GraphIndex, GraphIndex, EdgeWeight)]) -> Self {
        let mut connectivity = Connectivity {
            graph: Graph::from_edges(edges),
            non_cutting: Default::default(),
            prev: Default::default(),
            distance: HashMap::new(),
        };
        connectivity.update();
        connectivity
    }

    pub fn nodes(&self) -> Vec<GraphIndex> {
        self.graph
            .node_references()
            .map(|node| self.graph.to_index(node.id()))
            .collect()
    }

    fn update(&mut self) {
        let art_points = articulation_points(&self.graph);
        println!("Articulation points: {:?}", art_points.iter().map(|i| self.graph.to_index(*i)).collect::<Vec<usize>>());

        let non_cutting = self.nodes().iter()
            .filter_map(|node| match art_points.contains(&self.graph.from_index(*node)) {
                false => Some(*node),
                true => None,
            })
            .collect(); // For some reason, using filter here makes |node: &&usize| which collect does not like.

        let (distance, prev) = floyd_warshall_path(&self.graph, |e| *e.weight()).unwrap();
        println!("prev: {:?}", prev);
        println!("non cutting: {:?}", non_cutting);
        self.non_cutting = non_cutting;
        self.distance = distance;
        self.prev = prev;
    }

    pub fn remove_node(&mut self, i: GraphIndex) {
        if !self.non_cutting.contains(&i){
            panic!("Removing node {} will disconnect the graph!", i);
        }
        self.graph.remove_node(self.graph.from_index(i));
        self.update();
    }

    pub fn remove_edge(&mut self, i: GraphIndex, j:GraphIndex){
        self.graph.remove_edge(self.graph.find_edge(self.graph.from_index(i), self.graph.from_index(j)).unwrap());
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

impl<T: EdgeType>  Architecture for Connectivity<T> {
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

    use crate::architecture::{connectivity::{DirectedConnectivity, UndirectedConnectivity}, Architecture, EdgeWeight, GraphIndex};

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
        let new_diarchitecture: DirectedConnectivity = Connectivity::from_edges(&setup_simple());
        assert_eq!(new_diarchitecture.nodes(), vec![0, 1, 2, 3, 4, 5]);
        let new_unarchitecture: UndirectedConnectivity = Connectivity::from_edges(&setup_simple());
        assert_eq!(new_unarchitecture.nodes(), vec![0, 1, 2, 3, 4, 5]);
    }

    #[test]
    fn test_weighted_constructor() {
        let new_diarchitecture: DirectedConnectivity = Connectivity::from_weighted_edges(&setup_weighted());
        assert_eq!(new_diarchitecture.nodes(), vec![0, 1, 2, 3, 4, 5]);
        let new_unarchitecture: UndirectedConnectivity = Connectivity::from_weighted_edges(&setup_weighted());
        assert_eq!(new_unarchitecture.nodes(), vec![0, 1, 2, 3, 4, 5]);
    }

    #[test]
    fn test_best_simple_path() {
        let new_unarchitecture = UndirectedConnectivity::from_edges(&setup_simple());
        assert_eq!(vec![0, 1, 2, 4], new_unarchitecture.best_path(0, 4));

        let new_diarchitecture = DirectedConnectivity::from_edges(&setup_simple());
        assert_eq!(vec![0, 1, 2, 4], new_diarchitecture.best_path(0, 4));
    }

    #[test]
    fn test_best_weighted_path() {
        let new_diarchitecture = DirectedConnectivity::from_weighted_edges(&setup_weighted());
        assert_eq!(vec![0, 1, 2, 3, 4], new_diarchitecture.best_path(0, 4));
        let new_unarchitecture = UndirectedConnectivity::from_weighted_edges(&setup_weighted());
        assert_eq!(vec![0, 1, 2, 3, 4], new_unarchitecture.best_path(0, 4));
    }

    #[test]
    #[should_panic = "architecture does not contain node 6"]
    fn test_best_path_missing() {
        let new_architecture = UndirectedConnectivity::from_edges(&setup_simple());
        new_architecture.best_path(5, 6);
    }

    #[test]
    fn test_distance() {
        let directed_architecture = DirectedConnectivity::from_weighted_edges(&setup_weighted());
        assert_eq!(2, directed_architecture.distance(2, 4));
        assert_eq!(10, directed_architecture.distance(4, 5));
        assert_eq!(10, directed_architecture.distance(0, 4));
        let undirected_architecture = UndirectedConnectivity::from_weighted_edges(&setup_weighted());
        assert_eq!(2, undirected_architecture.distance(2, 4));
        assert_eq!(2, undirected_architecture.distance(4, 2));
        assert_eq!(5, undirected_architecture.distance(4, 5));
        assert_eq!(10, undirected_architecture.distance(0, 4));
    }

    #[test]
    #[should_panic = "architecture does not contain node 6"]
    fn test_distance_missing() {
        let directed_architecture = DirectedConnectivity::from_edges(&setup_simple());
        directed_architecture.distance(5, 6);
    }

    #[test]
    fn test_neighbors() {
        let directed_architecture = DirectedConnectivity::from_edges(&setup_simple());
        assert_eq!(vec![4, 3], directed_architecture.neighbors(2));
        let undirected_architecture = UndirectedConnectivity::from_edges(&setup_simple());
        assert_eq!(vec![4, 3, 1], undirected_architecture.neighbors(2));
    }

    #[test]
    #[should_panic = "architecture does not contain node 7"]
    fn test_neighbor_missing() {
        let directed_architecture = DirectedConnectivity::from_edges(&setup_simple());
        directed_architecture.distance(2, 7);
    }

    #[test]
    fn test_non_cutting() {
        let mut undirected_architecture = UndirectedConnectivity::from_edges(&setup_simple());
        assert_eq!(&vec![0,1,2,3,4,5], undirected_architecture.non_cutting());
        // Below failure is caused by a bug in petgraph articulation points for directed graphs: 0 is not an articulation point.
        let mut directed_architecture = DirectedConnectivity::from_edges(&setup_simple());
        assert_eq!(&vec![0,1,2,3,4,5], directed_architecture.non_cutting());
    }
}
