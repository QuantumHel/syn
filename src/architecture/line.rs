use super::Architecture;

/// Describes a line architecture. This ensures that consecutive numbers are adjacent. Allows for disjoint breaking of architecture in case this is required.
#[derive(Debug)]
pub struct Line {
    nodes: Vec<usize>,
    non_cutting: Vec<usize>,
    updated: bool,
}

impl Line {
    pub fn new(num_qubits: usize) -> Self {
        Line {
            nodes: (0..num_qubits).collect(),
            non_cutting: Vec::new(),
            updated: false,
        }
    }

    pub fn remove(&mut self, i: usize) {
        assert!(
            self.nodes.contains(&i),
            "architecture does not contain node {i}"
        );
        self.updated = false;
        self.nodes.retain(|x| *x != i);
    }

    /// Ensures that strict ordering is always enforced
    pub fn insert(&mut self, i: usize) {
        match self.nodes.binary_search(&i) {
            Ok(_) => panic!("architecture already contains node {i}"),
            Err(pos) => self.nodes.insert(pos, i),
        }
        self.updated = false;
    }
}

impl Architecture for Line {
    fn best_path(&self, i: usize, j: usize) -> Vec<usize> {
        let (smaller, larger) = (i.min(j), i.max(j));
        assert!(
            self.nodes.contains(&j),
            "architecture does not contain node {i}"
        );
        assert!(
            self.nodes.contains(&j),
            "architecture does not contain node {j}"
        );
        assert!(
            ((smaller + 1)..larger).all(|node| self.nodes.contains(&node)),
            "no path exists between {i} and {j}"
        );
        if i == smaller {
            (smaller..=larger).collect()
        } else {
            (smaller..=larger).rev().collect()
        }
    }

    fn distance(&self, i: usize, j: usize) -> usize {
        let (smaller, larger) = (i.min(j), i.max(j));
        assert!(
            self.nodes.contains(&i),
            "architecture does not contain node {i}"
        );
        assert!(
            self.nodes.contains(&j),
            "architecture does not contain node {j}"
        );
        assert!(
            ((smaller + 1)..larger).all(|node| self.nodes.contains(&node)),
            "no path exists between {i} and {j}"
        );
        larger - smaller
    }

    fn neighbors(&self, i: usize) -> Vec<usize> {
        assert!(
            self.nodes.contains(&i),
            "architecture does not contain node {i}"
        );
        let mut neighbors = Vec::new();
        if let Some(left_neighbor) = i.checked_sub(1) {
            if self.nodes.contains(&left_neighbor) {
                neighbors.push(left_neighbor);
            }
        }

        if let Some(right_neighbor) = i.checked_add(1) {
            if self.nodes.contains(&right_neighbor) {
                neighbors.push(right_neighbor);
            }
        }

        neighbors
    }

    fn non_cutting(&mut self) -> &Vec<usize> {
        if !self.updated {
            let mut non_cutting = Vec::new();
            non_cutting.push(self.nodes[0]);
            for nodes in self.nodes.windows(3) {
                if let &[node1, node2, node3] = nodes {
                    // Since strict ordering is enforced during insertion and creation, adding 1 to smaller node should not cause overflow.
                    if node1 + 1 != node2 || node2 + 1 != node3 {
                        non_cutting.push(node2);
                    }
                }
            }

            if self.nodes.len() > 1 {
                non_cutting.push(*self.nodes.last().unwrap());
            }
            self.non_cutting = non_cutting;
            self.updated = true;
        }
        &self.non_cutting
    }
}

#[cfg(test)]
mod tests {
    use crate::architecture::Architecture;

    use super::Line;

    #[test]
    fn test_complete() {
        let new_architecture = Line::new(5);

        assert_eq!(new_architecture.nodes, vec![0, 1, 2, 3, 4]);
    }

    #[test]
    fn test_insert() {
        let mut new_architecture = Line::new(5);
        new_architecture.insert(5);
        assert_eq!(new_architecture.nodes, vec![0, 1, 2, 3, 4, 5]);
    }

    #[test]
    #[should_panic = "architecture already contains node 1"]
    fn test_bad_insert() {
        let mut new_architecture = Line::new(5);
        new_architecture.insert(1);
    }

    #[test]
    fn test_remove() {
        let mut new_architecture = Line::new(5);
        new_architecture.remove(3);
        assert_eq!(new_architecture.nodes, vec![0, 1, 2, 4]);
    }

    #[test]
    fn test_remove_remove_insert() {
        let mut new_architecture = Line::new(7);
        new_architecture.remove(3);
        new_architecture.remove(4);
        new_architecture.remove(5);
        new_architecture.insert(4);
        assert_eq!(new_architecture.nodes, vec![0, 1, 2, 4, 6]);
    }

    #[test]
    #[should_panic = "architecture does not contain node 4"]
    fn test_bad_remove() {
        let mut new_architecture = Line::new(3);
        new_architecture.remove(4);
    }

    #[test]
    fn test_best_path() {
        let new_architecture = Line::new(6);
        assert_eq!(vec![2, 3, 4], new_architecture.best_path(2, 4));
    }

    #[test]
    #[should_panic = "architecture does not contain node 4"]
    fn test_best_path_missing() {
        let new_architecture = Line::new(3);
        new_architecture.best_path(4, 5);
    }

    #[test]
    fn test_distance() {
        let new_architecture = Line::new(6);
        assert_eq!(2, new_architecture.distance(2, 4));
    }

    #[test]
    #[should_panic = "architecture does not contain node 4"]
    fn test_distance_missing() {
        let new_architecture = Line::new(3);
        new_architecture.distance(4, 5);
    }

    #[test]
    fn test_neighbors() {
        let new_architecture = Line::new(4);
        assert_eq!(vec![1, 3], new_architecture.neighbors(2));
    }

    #[test]
    fn test_neighbors_0() {
        let new_architecture = Line::new(4);
        assert_eq!(vec![1], new_architecture.neighbors(0));
    }

    #[test]
    #[should_panic = "architecture does not contain node 3"]
    fn test_neighbor_missing() {
        let new_architecture = Line::new(3);
        new_architecture.distance(2, 3);
    }

    #[test]
    fn test_non_cutting() {
        let mut new_architecture = Line::new(5);
        assert_eq!(&vec![0, 4], new_architecture.non_cutting());
    }

    #[test]
    fn test_non_cutting_complex() {
        let mut new_architecture = Line::new(7);
        new_architecture.remove(1);
        new_architecture.remove(5);
        assert_eq!(&vec![0, 2, 4, 6], new_architecture.non_cutting());
    }

    #[test]
    fn test_non_cutting_cache() {
        let mut new_architecture = Line::new(7);
        assert_eq!(&vec![0, 6], new_architecture.non_cutting());
        new_architecture.remove(1);
        assert_eq!(&vec![0, 2, 6], new_architecture.non_cutting());
        new_architecture.remove(5);
        assert_eq!(&vec![0, 2, 4, 6], new_architecture.non_cutting());
    }
}
