pub trait Naive<T, G> {
    fn run(program: T, external_repr: &mut G);
}
