use crate::edge;
use edge::Edge;
use std::collections::{HashMap, VecDeque};

mod edge {
    use std::hash::Hash;

    #[derive(Debug, Clone, Eq)]
    pub struct Edge {
        pub(crate) edge: [usize; 2],
    }

    impl Edge {
        pub fn contains(&self, vertex: usize) -> bool {
            self.edge.contains(&vertex)
        }

        #[must_use]
        pub fn new(v1: usize, v2: usize) -> Self {
            Edge { edge: [v1, v2] }
        }
    }

    impl Hash for Edge {
        fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
            let mut edge = self.edge;
            edge.sort();
            edge.hash(state);
        }
    }

    impl PartialEq for Edge {
        fn eq(&self, other: &Self) -> bool {
            let mut self_edge = self.edge;
            let mut other_edge = other.edge;
            self_edge.sort();
            other_edge.sort();

            self_edge == other_edge
        }
    }

    #[macro_export]
    macro_rules! edges {
            ( $( $x:expr ),* $(,)? ) => {
                {
                    &[$(Edge{edge: ($x.into())}),*]
                }
            };
        }

    #[macro_export]
    macro_rules! edge {
        ( $x:expr, $y:expr ) => {{
            Edge { edge: [$x, $y] }
        }};
    }

    #[cfg(test)]
    mod tests {
        use std::{
            collections::HashMap,
            hash::{DefaultHasher, Hasher},
        };

        use super::*;

        #[test]
        fn test_constructor() {
            let e1 = Edge { edge: [1, 2] };
            let e2 = Edge::new(2, 1);

            assert_eq!(e1, e2);
        }

        #[test]
        fn test_hash() {
            let e1 = Edge { edge: [2, 1] };
            let e2 = Edge::new(2, 1);

            // Manually contructing Edge means vertices can be out of order and still pas equality.
            assert_eq!(e1, e2);

            let mut h1 = DefaultHasher::new();
            let mut h2 = DefaultHasher::new();
            e1.hash(&mut h1);
            e2.hash(&mut h2);
            // Hash should ensure that Edges containing the same vertices are treated as the same hashmap key.
            assert_eq!(h1.finish(), h2.finish());
        }

        #[test]
        fn test_hashmap() {
            let e1 = Edge { edge: [2, 1] };
            let e2 = Edge::new(1, 2);

            let mut edge_weights = HashMap::new();

            edge_weights.entry(e1.clone()).or_insert(1);
            edge_weights.entry(e2).and_modify(|v| *v = 5);

            assert_eq!(5, *edge_weights.get(&e1).unwrap());
        }

        #[test]
        fn test_edge_macro() {
            let edges = edges!((1, 2), (2, 3));

            assert_eq!(vec![Edge::new(2, 1), Edge::new(3, 2)], edges);
        }
    }
}

#[derive(Debug)]
pub struct Connectivity {
    size: usize,
    edges: Vec<Edge>,
    adjacency: HashMap<usize, Vec<usize>>,
    distance: HashMap<Edge, usize>,
}

impl Connectivity {
    pub fn line(size: usize, edge_weights: Option<HashMap<Edge, usize>>) -> Self {
        let edges = (0..(size - 1))
            .map(|q| Edge { edge: [q, (q + 1)] })
            .collect::<Vec<Edge>>();

        let adjacency = setup_adjacency(&edges);
        let distance = setup_distance(size, &adjacency, edge_weights);

        Connectivity {
            size,
            edges,
            adjacency,
            distance,
        }
    }

    pub fn cycle(size: usize, edge_weights: Option<HashMap<Edge, usize>>) -> Self {
        let edges = (0..size)
            .map(|q| Edge {
                edge: [q, (q + 1) % size],
            })
            .collect::<Vec<Edge>>();

        let adjacency = setup_adjacency(&edges);
        let distance = setup_distance(size, &adjacency, edge_weights);

        Connectivity {
            size,
            edges,
            adjacency,
            distance,
        }
    }

    pub fn complete(size: usize, edge_weights: Option<HashMap<Edge, usize>>) -> Self {
        let mut edges = Vec::new();

        for i in 0..size {
            edges.extend(((i + 1)..size).map(|q| Edge { edge: [i, q] }));
        }

        let adjacency = setup_adjacency(&edges);
        let distance = setup_distance(size, &adjacency, edge_weights);

        Connectivity {
            size,
            edges,
            adjacency,
            distance,
        }
    }

    pub fn grid(length: usize, height: usize, edge_weights: Option<HashMap<Edge, usize>>) -> Self {
        let mut edges = Vec::new();
        for h in 0..height {
            let vertical_shift = h * length;
            edges.extend(
                (vertical_shift..(length + vertical_shift - 1)).map(|q| Edge { edge: [q, q + 1] }),
            );
        }
        for h in 0..(height - 1) {
            let vertical_shift = h * length;
            edges.extend((vertical_shift..(length + vertical_shift)).map(|q| Edge {
                edge: [q, q + length],
            }));
        }

        let adjacency = setup_adjacency(&edges);
        let distance = setup_distance(length * height, &adjacency, edge_weights);

        Connectivity {
            size: length * height,
            edges,
            adjacency,
            distance,
        }
    }
}

fn setup_adjacency(edges: &[Edge]) -> HashMap<usize, Vec<usize>> {
    let mut adjacency = HashMap::new();

    for &Edge { edge: [i, j] } in edges.iter() {
        adjacency
            .entry(i)
            .and_modify(|nodes: &mut Vec<usize>| nodes.push(j))
            .or_insert(vec![j]);
        adjacency
            .entry(j)
            .and_modify(|nodes: &mut Vec<usize>| nodes.push(i))
            .or_insert(vec![i]);
    }

    adjacency
}

/// Breadth-first search to identify shortest distances in unweighted graph.
fn bfs(size: usize, adjacency: &HashMap<usize, Vec<usize>>) -> HashMap<Edge, usize> {
    let mut distance = HashMap::new();
    for i in 0..size {
        let mut queue = VecDeque::new();
        let mut visited = vec![false; size];
        visited[i] = true;
        distance.entry(Edge::new(i, i)).or_insert(0);
        queue.push_back(i);
        while !queue.is_empty() {
            let current = queue.pop_front().unwrap();
            for neighbor in adjacency[&current].iter() {
                if !visited[*neighbor] {
                    let dist = distance[&Edge::new(i, current)] + 1;
                    visited[*neighbor] = true;
                    distance.entry(Edge::new(i, *neighbor)).or_insert(dist);

                    queue.push_back(*neighbor);
                }
            }
        }
    }
    distance
}

/// Floyd-Warshall to identify shortest distances in a weighted, undirected graph.
fn floyd_warshall(
    size: usize,
    _adjacency: &HashMap<usize, Vec<usize>>,
    mut distance: HashMap<Edge, usize>,
) -> HashMap<Edge, usize> {
    // Set distance to self as 0.
    for i in 0..size {
        distance.entry(edge!(i, i)).or_insert(0);
    }

    // Set all edges not defined in `distance` as max/2 to prevent overflow.
    for i in 0..size {
        for j in i..size {
            distance.entry(edge!(i, j)).or_insert(usize::MAX / 2);
        }
    }

    for k in 0..size {
        for i in 0..size {
            for j in i..size {
                if distance[&edge!(i, j)] > distance[&edge!(i, k)] + distance[&edge!(k, j)] {
                    distance.insert(edge!(i, j), distance[&edge!(i, k)] + distance[&edge!(k, j)]);
                }
            }
        }
    }

    distance
}

fn setup_distance(
    size: usize,
    adjacency: &HashMap<usize, Vec<usize>>,
    edge_weights: Option<HashMap<Edge, usize>>,
) -> HashMap<Edge, usize> {
    if let Some(edge_weights) = edge_weights {
        floyd_warshall(size, adjacency, edge_weights)
    } else {
        bfs(size, adjacency)
    }
}

#[cfg(test)]
mod tests {
    use std::iter::zip;

    use edge::Edge;

    use crate::{edge, edges};

    use super::*;

    #[test]
    fn test_adjacency_line() {
        let edge_vec = edges![(0, 1), (1, 2), (2, 3), (3, 4)].to_vec();
        let adjacency = setup_adjacency(&edge_vec);

        let ref_adjacency = HashMap::from([
            (0, vec![1]),
            (1, vec![0, 2]),
            (2, vec![1, 3]),
            (3, vec![2, 4]),
            (4, vec![3]),
        ]);
        assert_eq!(ref_adjacency, adjacency);
    }

    #[test]
    fn test_adjacency_cycle() {
        let edges = edges![(0, 1), (1, 2), (2, 3), (3, 0)].to_vec();
        let adjacency = setup_adjacency(&edges);

        let ref_adjacency = HashMap::from([
            (0, vec![1, 3]),
            (1, vec![0, 2]),
            (2, vec![1, 3]),
            (3, vec![2, 0]),
        ]);
        assert_eq!(ref_adjacency, adjacency);
    }

    #[test]
    fn test_adjacency_complete() {
        let edges = edges![(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)].to_vec();
        let adjacency = setup_adjacency(&edges);

        let ref_adjacency = HashMap::from([
            (0, vec![1, 2, 3]),
            (1, vec![0, 2, 3]),
            (2, vec![0, 1, 3]),
            (3, vec![0, 1, 2]),
        ]);
        assert_eq!(ref_adjacency, adjacency);
    }

    #[test]
    fn test_adjacency_grid() {
        let edges = edges![
            (0, 1),
            (1, 2),
            (3, 4),
            (4, 5),
            (6, 7),
            (7, 8),
            (0, 3),
            (1, 4),
            (2, 5),
            (3, 6),
            (4, 7),
            (5, 8)
        ]
        .to_vec();
        let mut adjacency = setup_adjacency(&edges);

        // Entries are the same but permuted due to algorithm ordering
        let ref_adjacency = HashMap::from([
            (0, vec![1, 3]),
            (1, vec![0, 2, 4]),
            (2, vec![1, 5]),
            (3, vec![0, 4, 6]),
            (4, vec![1, 3, 5, 7]),
            (5, vec![2, 4, 8]),
            (6, vec![3, 7]),
            (7, vec![4, 6, 8]),
            (8, vec![5, 7]),
        ]);

        for (key, entry) in adjacency.iter_mut() {
            entry.sort();
            assert_eq!(ref_adjacency[key], *entry);
        }
        assert_eq!(ref_adjacency, adjacency);
    }

    #[test]
    fn test_distance_line() {
        let size = 5;
        let edges = edges![(0, 1), (1, 2), (2, 3), (3, 4)].to_vec();
        let adjacency = setup_adjacency(&edges);
        let distance = setup_distance(size, &adjacency, None);
        let ref_distance = HashMap::from_iter(zip(
            edges![
                (0, 0),
                (0, 1),
                (0, 2),
                (0, 3),
                (0, 4),
                (1, 1),
                (1, 2),
                (1, 3),
                (1, 4),
                (2, 2),
                (2, 3),
                (2, 4),
                (3, 3),
                (3, 4),
                (4, 4),
            ]
            .to_owned(),
            [0, 1, 2, 3, 4, 0, 1, 2, 3, 0, 1, 2, 0, 1, 0],
        ));

        assert_eq!(ref_distance, distance);
    }

    #[test]
    fn test_distance_cycle() {
        let size = 4;
        let edges = edges![(0, 1), (1, 2), (2, 3), (3, 0)].to_vec();
        let adjacency = setup_adjacency(&edges);
        let distance = setup_distance(size, &adjacency, None);
        let ref_distance = HashMap::from_iter(zip(
            edges![
                (0, 0),
                (0, 1),
                (0, 2),
                (0, 3),
                (1, 1),
                (1, 2),
                (1, 3),
                (2, 2),
                (2, 3),
                (3, 3),
            ]
            .to_owned(),
            [0, 1, 2, 1, 0, 1, 2, 0, 1, 0],
        ));

        assert_eq!(ref_distance, distance);
    }

    #[test]
    fn test_distance_complete() {
        let size = 4;
        let edges = edges![(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)].to_vec();
        let adjacency = setup_adjacency(&edges);
        let distance = setup_distance(size, &adjacency, None);
        let ref_distance = HashMap::from_iter(zip(
            edges![
                (0, 0),
                (0, 1),
                (0, 2),
                (0, 3),
                (1, 1),
                (1, 2),
                (1, 3),
                (2, 2),
                (2, 3),
                (3, 3),
            ]
            .to_owned(),
            [0, 1, 1, 1, 0, 1, 1, 0, 1, 0],
        ));

        assert_eq!(ref_distance, distance);
    }

    #[test]
    fn test_distance_grid() {
        let size = 9;
        let edges = edges![
            (0, 1),
            (1, 2),
            (3, 4),
            (4, 5),
            (6, 7),
            (7, 8),
            (0, 3),
            (1, 4),
            (2, 5),
            (3, 6),
            (4, 7),
            (5, 8),
        ]
        .to_owned();
        let adjacency = setup_adjacency(&edges);
        let distance = setup_distance(size, &adjacency, None);

        let ref_iter = zip(
            edges![
                (0, 0),
                (0, 1),
                (0, 2),
                (0, 3),
                (0, 4),
                (0, 5),
                (0, 6),
                (0, 7),
                (0, 8),
            ]
            .to_owned(),
            [0, 1, 2, 1, 2, 3, 2, 3, 4],
        );

        let ref_iter = ref_iter.chain(zip(
            edges![
                (1, 1),
                (1, 2),
                (1, 3),
                (1, 4),
                (1, 5),
                (1, 6),
                (1, 7),
                (1, 8),
            ]
            .to_owned(),
            [0, 1, 2, 1, 2, 3, 2, 3],
        ));

        let ref_iter = ref_iter.chain(zip(
            edges![(2, 2), (2, 3), (2, 4), (2, 5), (2, 6), (2, 7), (2, 8),].to_owned(),
            [0, 3, 2, 1, 4, 3, 2],
        ));

        let ref_iter = ref_iter.chain(zip(
            edges![(3, 3), (3, 4), (3, 5), (3, 6), (3, 7), (3, 8),].to_owned(),
            [0, 1, 2, 1, 2, 3],
        ));

        let ref_iter = ref_iter.chain(zip(
            edges![(4, 4), (4, 5), (4, 6), (4, 7), (4, 8),].to_owned(),
            [0, 1, 2, 1, 2],
        ));

        let ref_iter = ref_iter.chain(zip(
            edges![(5, 5), (5, 6), (5, 7), (5, 8),].to_owned(),
            [0, 3, 2, 1],
        ));

        let ref_iter = ref_iter.chain(zip(edges![(6, 6), (6, 7), (6, 8),].to_owned(), [0, 1, 2]));

        let ref_iter = ref_iter.chain(zip(edges![(7, 7), (7, 8),].to_owned(), [0, 1]));

        let ref_iter = ref_iter.chain(zip(edges![(8, 8),].to_owned(), [0]));

        let ref_distance = HashMap::from_iter(ref_iter);

        assert_eq!(ref_distance, distance);
    }

    #[test]
    fn test_distance_weighted_line() {
        let size = 5;
        let edges = edges![(0, 1), (1, 2), (2, 3), (3, 4)].to_vec();
        let adjacency = setup_adjacency(&edges);

        let edge_weights = HashMap::from_iter(zip(edges, [2, 3, 4, 5]));

        let distance = setup_distance(size, &adjacency, Some(edge_weights));
        let ref_distance = HashMap::from_iter(zip(
            edges![
                (0, 0),
                (0, 1),
                (0, 2),
                (0, 3),
                (0, 4),
                (1, 1),
                (1, 2),
                (1, 3),
                (1, 4),
                (2, 2),
                (2, 3),
                (2, 4),
                (3, 3),
                (3, 4),
                (4, 4),
            ]
            .to_owned(),
            [0, 2, 5, 9, 14, 0, 3, 7, 12, 0, 4, 9, 0, 5, 0],
        ));

        assert_eq!(ref_distance, distance);
    }

    #[test]
    fn test_distance_weighted_grid() {
        let size = 9;
        let edges = edges![
            (0, 1),
            (1, 2),
            (3, 4),
            (4, 5),
            (6, 7),
            (7, 8),
            (0, 3),
            (1, 4),
            (2, 5),
            (3, 6),
            (4, 7),
            (5, 8),
        ]
        .to_owned();
        let adjacency = setup_adjacency(&edges);
        let edge_weights = HashMap::from_iter(zip(edges, [2, 2, 1, 1, 2, 2, 3, 1, 3, 3, 1, 3]));
        let distance = setup_distance(size, &adjacency, Some(edge_weights));

        let ref_iter = zip(
            edges![
                (0, 0),
                (0, 1),
                (0, 2),
                (0, 3),
                (0, 4),
                (0, 5),
                (0, 6),
                (0, 7),
                (0, 8),
            ]
            .to_owned(),
            [0, 2, 4, 3, 3, 4, 6, 4, 6],
        );

        let ref_iter = ref_iter.chain(zip(
            edges![
                (1, 1),
                (1, 2),
                (1, 3),
                (1, 4),
                (1, 5),
                (1, 6),
                (1, 7),
                (1, 8),
            ]
            .to_owned(),
            [0, 2, 2, 1, 2, 4, 2, 4],
        ));

        let ref_iter = ref_iter.chain(zip(
            edges![(2, 2), (2, 3), (2, 4), (2, 5), (2, 6), (2, 7), (2, 8),].to_owned(),
            [0, 4, 3, 3, 6, 4, 6],
        ));

        let ref_iter = ref_iter.chain(zip(
            edges![(3, 3), (3, 4), (3, 5), (3, 6), (3, 7), (3, 8),].to_owned(),
            [0, 1, 2, 3, 2, 4],
        ));

        let ref_iter = ref_iter.chain(zip(
            edges![(4, 4), (4, 5), (4, 6), (4, 7), (4, 8),].to_owned(),
            [0, 1, 3, 1, 3],
        ));

        let ref_iter = ref_iter.chain(zip(
            edges![(5, 5), (5, 6), (5, 7), (5, 8),].to_owned(),
            [0, 4, 2, 3],
        ));

        let ref_iter = ref_iter.chain(zip(edges![(6, 6), (6, 7), (6, 8),].to_owned(), [0, 2, 4]));

        let ref_iter = ref_iter.chain(zip(edges![(7, 7), (7, 8),].to_owned(), [0, 2]));

        let ref_iter = ref_iter.chain(zip(edges![(8, 8),].to_owned(), [0]));

        let ref_distance = HashMap::from_iter(ref_iter);

        assert_eq!(ref_distance, distance);
    }

    #[test]
    fn test_distance_weighted_cycle() {
        let size = 4;
        let edges = edges![(0, 1), (1, 2), (2, 3), (3, 0)].to_vec();
        let adjacency = setup_adjacency(&edges);

        let mut edge_weights = HashMap::from_iter(zip(edges, [2, 3, 4, 4]));
        let entry = edge_weights.entry(edge![0, 3]);
        println!("entry: {:?}", entry);
        let distance = setup_distance(size, &adjacency, Some(edge_weights));
        let ref_distance = HashMap::from_iter(zip(
            edges![
                (0, 0),
                (0, 1),
                (0, 2),
                (0, 3),
                (1, 1),
                (1, 2),
                (1, 3),
                (2, 2),
                (2, 3),
                (3, 3),
            ]
            .to_owned(),
            [0, 2, 5, 4, 0, 3, 6, 0, 4, 0],
        ));

        assert_eq!(ref_distance, distance);
    }

    #[test]
    fn test_connectivity_line() {
        let num_qubits = 5;
        let connectivity = Connectivity::line(num_qubits, None);

        let ref_connectivity = Connectivity {
            size: num_qubits,
            edges: edges![(0, 1), (1, 2), (2, 3), (3, 4)].to_vec(),
            adjacency: Default::default(),
            distance: Default::default(),
        };
        assert_eq!(ref_connectivity.edges, connectivity.edges);
    }

    #[test]
    fn test_connectivity_cycle() {
        let num_qubits = 4;
        let connectivity = Connectivity::cycle(num_qubits, None);

        let ref_connectivity = Connectivity {
            size: num_qubits,
            edges: edges![(0, 1), (1, 2), (2, 3), (3, 0)].to_vec(),
            adjacency: Default::default(),
            distance: Default::default(),
        };
        assert_eq!(ref_connectivity.edges, connectivity.edges);
    }

    #[test]
    fn test_connectivity_complete() {
        let num_qubits = 4;
        let connectivity = Connectivity::complete(num_qubits, None);

        let ref_connectivity = Connectivity {
            size: num_qubits,
            edges: edges![(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)].to_vec(),
            adjacency: Default::default(),
            distance: Default::default(),
        };
        assert_eq!(ref_connectivity.edges, connectivity.edges);
    }

    #[test]
    fn test_connectivity_grid() {
        let length = 3;
        let height = 3;
        let connectivity = Connectivity::grid(length, height, None);

        let ref_connectivity = Connectivity {
            size: length * height,
            edges: edges![
                (0, 1),
                (1, 2),
                (3, 4),
                (4, 5),
                (6, 7),
                (7, 8),
                (0, 3),
                (1, 4),
                (2, 5),
                (3, 6),
                (4, 7),
                (5, 8)
            ]
            .to_vec(),
            adjacency: Default::default(),
            distance: Default::default(),
        };
        assert_eq!(ref_connectivity.edges, connectivity.edges);
    }
}
