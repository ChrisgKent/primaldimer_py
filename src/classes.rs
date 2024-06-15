use debruijn::dna_string::*;
use pyo3::prelude::*;

#[pyclass]
pub struct DebruijnFKmer {
    #[pyo3(get)]
    pub end: usize,
    pub seqs: Vec<DnaString>,
    #[pyo3(get)]
    pub starts: Vec<usize>,
    #[pyo3(get)]
    pub start: usize,
    #[pyo3(get)]
    pub encodedseqs: Vec<Vec<u8>>,
}
#[pymethods]
impl DebruijnFKmer {
    #[new]
    pub fn new(end: usize, seqs: Vec<String>) -> Self {
        let starts = seqs.iter().map(|s| end - s.len()).collect();
        let dnaseqs: Vec<DnaString> = seqs.iter().map(|s| DnaString::from_dna_string(s)).collect();
        let encoded_seqs: Vec<Vec<u8>> = dnaseqs.iter().map(|s| s.to_bytes()).collect();
        DebruijnFKmer {
            end,
            seqs: dnaseqs,
            starts,
            encodedseqs: encoded_seqs,
            start: *starts.iter().min().unwrap(),
        }
    }
    pub fn seqs(&self) -> Vec<String> {
        self.seqs.iter().map(|s| s.to_string()).collect()
    }
}
#[pyclass]
#[derive(Debug)]
pub struct DebruijnRKmer {
    #[pyo3(get)]
    pub start: usize,
    seqs: Vec<DnaString>,
    #[pyo3(get)]
    pub _ends: Vec<usize>,
    #[pyo3(get)]
    pub _end: usize,
    #[pyo3(get)]
    pub _encodedseqs: Vec<Vec<u8>>,
}
#[pymethods]
impl DebruijnRKmer {
    #[new]
    pub fn new(start: usize, seqs: Vec<String>) -> Self {
        let _ends = seqs.iter().map(|s| start + s.len()).collect();
        let dnaseqs: Vec<DnaString> = seqs.iter().map(|s| DnaString::from_dna_string(s)).collect();
        let encoded_seqs: Vec<Vec<u8>> = dnaseqs.iter().map(|s| s.to_bytes()).collect();
        DebruijnRKmer {
            start,
            seqs: dnaseqs,
            _ends,
            _encodedseqs: encoded_seqs,
            _end: *_ends.iter().max().unwrap(),
        }
    }
    pub fn seqs(&self) -> Vec<String> {
        self.seqs.iter().map(|s| s.to_string()).collect()
    }
}
