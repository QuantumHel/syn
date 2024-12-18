pub trait Custom<T, G> {
    fn run_custom(
        program: T,
        external_repr: &mut G,
        custom_columns: Vec<usize>,
        custom_rows: Vec<usize>,
    );
}
