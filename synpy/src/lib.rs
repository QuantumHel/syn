use pyo3::prelude::*;

#[pyfunction]
fn add(x: i32, y: i32) -> PyResult<i32> {
    Ok(x + y)
}

#[pymodule]
fn synpy(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(add, m)?)?;
    Ok(())
}