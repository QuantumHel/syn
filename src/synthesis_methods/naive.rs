pub trait Naive<T, G> {
    fn run_naive(program: T, external_repr: &mut G);
}
