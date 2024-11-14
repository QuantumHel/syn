use super::{Architecture, GraphIndex};

#[derive(Debug)]
pub struct Complete {
    nodes: Vec<GraphIndex>,
}

impl Complete {
    pub fn new(num_qubits: usize) -> Self {
        Complete {
            nodes: (0..num_qubits).collect(),
        }
    }

    pub fn remove(&mut self, i: usize) {
        assert!(
            self.nodes.contains(&i),
            "architecture does not contain node {i}"
        );
        self.nodes.retain(|x| *x != i);
    }

    pub fn insert(&mut self, i: usize) {
        assert!(
            !self.nodes.contains(&i),
            "architecture already contains node {i}"
        );
        self.nodes.push(i);
    }
}

impl Architecture for Complete {
    fn best_path(&self, i: GraphIndex, j: GraphIndex) -> Vec<GraphIndex> {
        assert!(
            self.nodes.contains(&i),
            "architecture does not contain node {i}"
        );
        assert!(
            self.nodes.contains(&j),
            "architecture does not contain node {j}"
        );
        vec![i, j]
    }

    fn distance(&self, i: GraphIndex, j: GraphIndex) -> GraphIndex {
        assert!(
            self.nodes.contains(&i),
            "architecture does not contain node {i}"
        );
        assert!(
            self.nodes.contains(&j),
            "architecture does not contain node {j}"
        );
        1
    }

    fn neighbors(&self, i: GraphIndex) -> Vec<GraphIndex> {
        assert!(
            self.nodes.contains(&i),
            "architecture does not contain node {i}"
        );
        self.nodes.iter().filter(|x| **x != i).copied().collect()
    }

    fn non_cutting(&mut self) -> &Vec<GraphIndex> {
        &self.nodes
    }
}

#[cfg(test)]
mod tests {
    use crate::architecture::Architecture;

    use super::Complete;

    #[test]
    fn test_complete() {
        let new_architecture = Complete::new(5);

        assert_eq!(new_architecture.nodes, vec![0, 1, 2, 3, 4]);
    }

    #[test]
    fn test_insert() {
        let mut new_architecture = Complete::new(5);
        new_architecture.insert(5);
        assert_eq!(new_architecture.nodes, vec![0, 1, 2, 3, 4, 5]);
    }

    #[test]
    #[should_panic = "architecture already contains node 1"]
    fn test_bad_insert() {
        let mut new_architecture = Complete::new(5);
        new_architecture.insert(1);
    }

    #[test]
    fn test_remove() {
        let mut new_architecture = Complete::new(5);
        new_architecture.remove(3);
        assert_eq!(new_architecture.nodes, vec![0, 1, 2, 4]);
    }

    #[test]
    #[should_panic = "architecture does not contain node 4"]
    fn test_bad_remove() {
        let mut new_architecture = Complete::new(3);
        new_architecture.remove(4);
    }

    #[test]
    fn test_best_path() {
        let new_architecture = Complete::new(3);
        assert_eq!(vec![1, 2], new_architecture.best_path(1, 2));
    }

    #[test]
    #[should_panic = "architecture does not contain node 4"]
    fn test_best_path_missing() {
        let new_architecture = Complete::new(3);
        new_architecture.best_path(4, 5);
    }

    #[test]
    fn test_distance() {
        let new_architecture = Complete::new(3);
        assert_eq!(1, new_architecture.distance(1, 2));
    }

    #[test]
    #[should_panic = "architecture does not contain node 4"]
    fn test_distance_missing() {
        let new_architecture = Complete::new(3);
        new_architecture.distance(4, 5);
    }

    #[test]
    fn test_neighbors() {
        let new_architecture = Complete::new(4);
        assert_eq!(vec![0, 1, 3], new_architecture.neighbors(2));
    }

    #[test]
    #[should_panic = "architecture does not contain node 3"]
    fn test_neighbor_missing() {
        let new_architecture = Complete::new(3);
        new_architecture.distance(2, 3);
    }

    #[test]
    fn test_non_cutting() {
        let mut new_architecture = Complete::new(4);
        assert_eq!(&vec![0, 1, 2, 3], new_architecture.non_cutting());
    }
}
