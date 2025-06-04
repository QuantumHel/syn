pub mod connectivity;

type GraphIndex = usize;
type EdgeWeight = usize;
type NodeWeight = usize;

#[derive(Debug, PartialEq)]
pub enum LadderError {
    RootNotFound,
}

pub trait Architecture {
    fn best_path(&self, i: GraphIndex, j: GraphIndex) -> Vec<GraphIndex>;
    fn distance(&self, i: GraphIndex, j: GraphIndex) -> GraphIndex;
    fn neighbors(&self, i: GraphIndex) -> Vec<GraphIndex>;
    fn non_cutting(&mut self) -> &Vec<GraphIndex>;
    fn get_cx_ladder(
        &self,
        nodes: &[GraphIndex],
        root: &GraphIndex,
    ) -> Result<Vec<(usize, usize)>, LadderError>;
}
