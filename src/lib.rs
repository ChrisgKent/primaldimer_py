use primaldimer_rs::{do_pools_interact, do_seqs_interact};
use pyo3::prelude::*;

/// Formats the sum of two numbers as string.
#[pyfunction]
fn do_seqs_interact_py(seq1: &str, seq2: &str, t: i64) -> bool {
    return do_seqs_interact(seq1, seq2, t);
}
#[pyfunction]
fn do_pools_interact_py(pool1: Vec<&str>, pool2: Vec<&str>, t: i64) -> bool {
    return do_pools_interact(pool1, pool2, t);
}

/// A Python module implemented in Rust.
#[pymodule]
fn primaldimer_py(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(do_pools_interact_py, m)?)?;
    m.add_function(wrap_pyfunction!(do_seqs_interact_py, m)?)?;
    Ok(())
}
