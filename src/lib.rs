mod primaldimer;

use pyo3::prelude::*;

use debruijn::dna_string::*;

#[pyclass(subclass)]
pub struct Kmer {
    #[pyo3(get)]
    pub encodedseqs: Vec<Vec<usize>>,
}
#[pymethods]
impl Kmer {
    #[new]
    pub fn new(_idx: usize, seqs: Vec<String>) -> Self {
        // Check that the sequences are valid
        for seq in &seqs {
            if !seq.chars().all(|c| "ATCG".contains(c)) {
                panic!("Sequence contains not ACGT bases: {}", seq);
            }
        }
        // Encode the sequences
        let mut encoded_seqs: Vec<Vec<usize>> = seqs
            .iter()
            .map(|s| {
                DnaString::from_dna_string(s)
                    .to_bytes()
                    .iter()
                    .map(|&x| x as usize)
                    .collect()
            })
            .collect();
        // Sort and dedup the sequences
        encoded_seqs.sort_unstable();
        encoded_seqs.dedup();

        let encodedseqs = encoded_seqs;

        Kmer { encodedseqs }
    }

    fn into_bytes(&self) -> Vec<Vec<u8>> {
        // Return the sequences in bytes
        self.encodedseqs
            .iter()
            .map(|s| s.iter().map(|&x| x as u8).collect())
            .collect()
    }

    #[getter]
    pub fn seqs(&self) -> Vec<String> {
        // Return the sequences in ATCG format
        self.into_bytes()
            .iter()
            .map(|s| DnaString::from_bytes(s).to_string())
            .collect()
    }

    pub fn lens(&self) -> Vec<usize> {
        // Return the lengths of the sequences
        self.encodedseqs.iter().map(|s| s.len()).collect()
    }
}

fn do_kmers_interact(kmer1: &Kmer, kmer2: &Kmer, t: f64) -> bool {
    // Check if two kmers interact
    for seq1 in &kmer1.encodedseqs {
        for seq2 in &kmer2.encodedseqs {
            if primaldimer::does_seq1_extend(&seq1, &seq2, t)
                | primaldimer::does_seq1_extend(&seq2, &seq1, t)
            {
                return true;
            }
        }
    }
    false
}

#[pyfunction]
fn which_kmers_pools_interact(
    py: Python<'_>,
    kmers1: Vec<Py<Kmer>>,
    kmers2: Vec<Py<Kmer>>,
    t: f64,
    calc_all: bool,
) -> PyResult<Vec<(Py<Kmer>, Py<Kmer>)>> {
    // Interaction tuples
    let mut interacting_kmers: Vec<(Py<Kmer>, Py<Kmer>)> = Vec::new();

    // Check if two pools of kmers interact
    for kmer1 in &kmers1 {
        for kmer2 in &kmers2 {
            if do_kmers_interact(&kmer1.as_ref(py).borrow(), &kmer2.as_ref(py).borrow(), t) {
                interacting_kmers.push((kmer1.clone(), kmer2.clone()));
                // Early return if we only want to know if any interact
                if !calc_all {
                    return Ok(interacting_kmers);
                }
            }
        }
    }
    return Ok(interacting_kmers);
}

#[pyfunction]
fn calc_at_offset_py(seq1: &str, seq2: &str, offset: i32) -> f64 {
    //Provide strings in 5'-3'
    // This will return the score for this offset
    let seq1 = primaldimer::encode_base(seq1);
    let mut seq2 = primaldimer::encode_base(seq2);
    seq2.reverse();

    match primaldimer::calc_at_offset(&seq1, &seq2, offset) {
        Some(score) => return score,
        None => return 100.,
    };
}
#[pyfunction]
fn do_seqs_interact_py(seq1: &str, seq2: &str, t: f64) -> bool {
    return primaldimer::do_seqs_interact(seq1, seq2, t);
}
#[pyfunction]
fn do_pools_interact_py(pool1: Vec<&str>, pool2: Vec<&str>, t: f64) -> bool {
    return primaldimer::do_pools_interact(pool1, pool2, t);
}

/// A Python module implemented in Rust.
#[pymodule]
fn primaldimer_py(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(do_pools_interact_py, m)?)?;
    m.add_function(wrap_pyfunction!(do_seqs_interact_py, m)?)?;
    m.add_function(wrap_pyfunction!(calc_at_offset_py, m)?)?;
    m.add_function(wrap_pyfunction!(which_kmers_pools_interact, m)?)?;
    m.add_class::<Kmer>()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_kmer_new() {
        // Test creating a new Kmer instance with valid sequences
        let kmer = Kmer::new(0, vec!["ATCG".to_string(), "GCTA".to_string()]);
        assert_eq!(kmer.seqs(), vec!["ATCG", "GCTA"]);
        assert_eq!(kmer.lens(), vec![4, 4]);
    }

    #[test]
    fn test_kmer_dedupe_order() {
        // Test creating a new Kmer instance with sequences in different order
        let kmer = Kmer::new(
            0,
            vec![
                "G".to_string(),
                "C".to_string(),
                "A".to_string(),
                "T".to_string(),
                "T".to_string(),
            ],
        );
        assert_eq!(kmer.seqs(), vec!["A", "C", "G", "T"]);
    }

    #[test]
    #[should_panic(expected = "Sequence contains not ACGT bases: ATCGX")]
    fn test_kmer_new_invalid_seq() {
        // Test creating a new Kmer instance with an invalid sequence
        Kmer::new(0, vec!["ATCG".to_string(), "ATCGX".to_string()]);
    }

    #[test]
    fn test_kmer_seqs() {
        // Test getting the sequences in ATCG format
        let kmer = Kmer::new(0, vec!["ATCG".to_string(), "GCTA".to_string()]);
        assert_eq!(kmer.seqs(), vec!["ATCG", "GCTA"]);
    }

    #[test]
    fn test_kmer_lens() {
        // Test getting the lengths of the sequences
        let kmer = Kmer::new(0, vec!["ATCG".to_string(), "GCTAA".to_string()]);
        assert_eq!(kmer.lens(), vec![4, 5]);
    }
}
