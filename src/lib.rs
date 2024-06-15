use pyo3::prelude::*;

mod primaldimer;
use primaldimer::*;

fn get_r_window<'a>(
    rkmers: &'a Vec<&'a DebruijnRKmer>,
    start: usize,
    end: usize,
) -> Vec<&'a DebruijnRKmer> {
    // Binary search to return Rkmers in the window
    let mut pos_rkmers: Vec<&DebruijnRKmer> = Vec::new();
    let n_kmers = rkmers.len();
    let mut high = n_kmers - 1;
    let mut low = 0;
    let mut mid = 0;

    while low <= high {
        mid = (low + high) / 2;
        // If the midpoint is inside the window, then we need to find the first
        if start <= rkmers[mid].start {
            loop {
                if mid == 0 {
                    break;
                } else if rkmers[mid - 1].start >= start {
                    mid -= 1;
                } else {
                    break;
                }
            }
            // Mid is now the first value to walk forwards
            loop {
                if mid < n_kmers && rkmers[mid].start <= end {
                    pos_rkmers.push(rkmers[mid]);
                    mid += 1;
                } else {
                    return pos_rkmers;
                }
            }
        } else if rkmers[mid].start < start {
            low = mid + 1;
        } else if rkmers[mid].start > end {
            high = mid - 1;
        }
    }
    return pos_rkmers;
}

#[pyfunction]
pub fn generate_primerpairs<'a>(
    py: Python<'_>,
    fkmers: Vec<Py<DebruijnFKmer>>,
    rkmers: Vec<Py<DebruijnRKmer>>,
    amp_size_min: usize,
    amp_size_max: usize,
    dimerscore: f64,
) -> Vec<(DebruijnFKmer, DebruijnRKmer)> {
    let mut primerpairs: Vec<(&DebruijnFKmer, &DebruijnRKmer)> = Vec::new();

    for fkmer in fkmers.iter() {
        // Find all possible Rkmers in the window
        let pos_rkmers = get_r_window(
            &rkmers,
            fkmer.borrow(py).start + amp_size_min,
            fkmer.borrow(py).start + amp_size_max,
        );
        // Write a single threaded version of interchecker
        for rkmer in pos_rkmers.iter() {
            if !do_kmers_interact(fkmer, rkmer, dimerscore) {
                primerpairs.push((fkmer, rkmer));
            }
        }
    }

    return primerpairs;
}

#[pyclass]
struct FKmer {
    #[pyo3(get)]
    end: usize,
    #[pyo3(get)]
    seqs: Vec<String>,
    encoded_seqs: Vec<Vec<usize>>,
}
#[pymethods]
impl FKmer {
    #[new]
    fn new(end: usize, seqs: Vec<String>) -> Self {
        let encoded_seqs = seqs.iter().map(|x| encode_base(x)).collect();

        FKmer {
            end,
            seqs,
            encoded_seqs,
        }
    }
    pub fn to_bed(
        &self,
        chromname: &str,
        amplicon_prefix: &str,
        amplicon_number: usize,
        pool: usize,
    ) -> String {
        let mut bed: String = String::new();
        let mut counter: u32 = 0;

        for seq in self.seqs.iter() {
            bed.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                chromname,
                self.end - seq.len(),
                self.end,
                format!("{}_{}_LEFT_{}", amplicon_prefix, amplicon_number, counter),
                pool,
                "+",
                seq
            ));
            counter += 1;
        }
        bed
    }
    pub fn encoded_seqs(&self) -> Vec<Vec<usize>> {
        self.encoded_seqs.clone()
    }
}

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
    m.add_class::<FKmer>()?;
    m.add_class::<DebruijnFKmer>()?;
    Ok(())
}
