use std::hash::{Hash, Hasher};

#[derive(Debug, Clone, Eq)]
pub struct Kmer {
    sequence: Vec<u8>,
}

impl Kmer {
    pub fn new(sequence: &[u8]) -> Self {
        Kmer {
            sequence: sequence.to_vec(),
        }
    }

    pub fn sequence(&self) -> String {
        String::from_utf8_lossy(&self.sequence).into_owned()
    }
}

impl Hash for Kmer {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.sequence.hash(state);
    }
}

impl PartialEq for Kmer {
    fn eq(&self, other: &Self) -> bool {
        self.sequence == other.sequence
    }
}