use primaldimer_rs::{calc_at_offset, do_pools_interact, do_seqs_interact, encode_base};
use pyo3::prelude::*;

#[pyfunction]
fn calc_at_offset_py(seq1: &str, seq2: &str, offset: i32) -> f64 {
    //Provide strings in 5'-3'
    // This will return the score for this offset
    let seq1 = encode_base(seq1);
    let mut seq2 = encode_base(seq2);
    seq2.reverse();

    match calc_at_offset(&seq1, &seq2, offset) {
        Some(score) => return score,
        None => return 100.,
    };
}
#[pyfunction]
fn do_seqs_interact_py(seq1: &str, seq2: &str, t: f64) -> bool {
    return do_seqs_interact(seq1, seq2, t);
}
#[pyfunction]
fn do_pools_interact_py(pool1: Vec<&str>, pool2: Vec<&str>, t: f64) -> bool {
    return do_pools_interact(pool1, pool2, t);
}

/// A Python module implemented in Rust.
#[pymodule]
fn primaldimer_py(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(do_pools_interact_py, m)?)?;
    m.add_function(wrap_pyfunction!(do_seqs_interact_py, m)?)?;
    m.add_function(wrap_pyfunction!(calc_at_offset_py, m)?)?;
    Ok(())
}
