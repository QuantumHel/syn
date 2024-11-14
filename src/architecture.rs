pub mod complete;
pub mod connectivity;
pub mod line;

type GraphIndex = usize;
type EdgeWeight = usize;
type NodeWeight = ();

pub trait Architecture {
    fn best_path(&self, i: GraphIndex, j: GraphIndex) -> Vec<GraphIndex>;
    fn distance(&self, i: GraphIndex, j: GraphIndex) -> GraphIndex;
    fn neighbors(&self, i: GraphIndex) -> Vec<GraphIndex>;
    fn non_cutting(&mut self) -> &Vec<GraphIndex>;
}
