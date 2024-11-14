use crate::{datastructures::PropagateClifford, ir::CliffordGates};

pub trait Naive<T, G> {
    fn run(program: T, external_repr: G);
}
