pub trait Naive<T, G> {
    fn run_naive(program: T, external_repr: &mut G);
}

pub trait NaiveAdjoint<T, G> {
    fn run_naive_adjoint(program: T, external_repr: &mut G);
}
