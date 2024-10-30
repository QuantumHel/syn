pub mod complete;
pub mod line;

pub trait Architecture {
    fn best_path(&self, i: usize, j: usize) -> Vec<usize>;
    fn distance(&self, i: usize, j: usize) -> usize;
    fn neighbors(&self, i: usize) -> Vec<usize>;
    fn non_cutting(&mut self) -> &Vec<usize>;
}
