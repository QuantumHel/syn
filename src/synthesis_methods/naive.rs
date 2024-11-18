use crate::{data_structures::PropagateClifford, ir::CliffordGates};

pub trait Naive<T, G> {
    fn run(program: T, external_repr: G);
}
