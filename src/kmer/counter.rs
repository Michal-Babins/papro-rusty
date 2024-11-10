use std::collections::HashMap;
use anyhow::Result;
use dashmap::DashMap;
use rayon::prelude::*;

use super::types::Kmer;

pub struct KmerCounter {
    k: usize,
    counts: DashMap<Kmer, usize>,
}

impl KmerCounter {
    /// Create a new KmerCounter with specified k-mer size
    pub fn new(k: usize) -> Self {
        KmerCounter {
            k,
            counts: DashMap::new(),
        }
    }

    /// Count k-mers in a sequence
    pub fn count_sequence(&self, sequence: &[u8]) -> Result<()> {
        if sequence.len() < self.k {
            return Ok(());
        }

        // Create windows of size k and count them
        sequence.windows(self.k).for_each(|window| {
            let kmer = Kmer::new(window);
            self.counts.entry(kmer).and_modify(|count| *count += 1).or_insert(1);
        });

        Ok(())
    }

    /// Process sequences in parallel using rayon
    pub fn count_sequences<I>(&self, sequences: I) -> Result<()>
    where
        I: ParallelIterator<Item = Vec<u8>>,
    {
        sequences.try_for_each(|seq| self.count_sequence(&seq))?;
        Ok(())
    }

    /// Get k-mer counts as a regular HashMap
    pub fn get_counts(&self) -> HashMap<String, usize> {
        self.counts
            .iter()
            .map(|entry| (entry.key().sequence(), *entry.value()))
            .collect()
    }

    /// Get the k-mer size
    pub fn kmer_size(&self) -> usize {
        self.k
    }

    /// Get the number of unique k-mers
    pub fn unique_kmers(&self) -> usize {
        self.counts.len()
    }

    /// Get the total number of k-mers (including duplicates)
    pub fn total_kmers(&self) -> usize {
        self.counts.iter().map(|entry| *entry.value()).sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basic_counting() {
        let counter = KmerCounter::new(3);
        counter.count_sequence(b"ATGATG").unwrap();
        
        let counts = counter.get_counts();
        assert_eq!(counts.get("ATG").unwrap(), &2);
        assert_eq!(counts.get("TGA").unwrap(), &1);
        assert_eq!(counts.get("GAT").unwrap(), &1);
    }

    #[test]
    fn test_short_sequence() {
        let counter = KmerCounter::new(3);
        counter.count_sequence(b"AT").unwrap(); // shorter than k
        assert_eq!(counter.unique_kmers(), 0);
        assert_eq!(counter.total_kmers(), 0);
    }

    #[test]
    fn test_multiple_sequences() {
        let counter = KmerCounter::new(2);
        counter.count_sequence(b"ATCG").unwrap();
        counter.count_sequence(b"CGAT").unwrap();
        
        let counts = counter.get_counts();
        assert_eq!(counts.get("AT").unwrap(), &2);
        assert_eq!(counts.get("TC").unwrap(), &1);
        assert_eq!(counts.get("CG").unwrap(), &2);
        assert_eq!(counts.get("GA").unwrap(), &1);
    }

    #[test]
    fn test_empty_sequence() {
        let counter = KmerCounter::new(3);
        counter.count_sequence(b"").unwrap();
        assert_eq!(counter.unique_kmers(), 0);
        assert_eq!(counter.total_kmers(), 0);
    }

    #[test]
    fn test_parallel_counting() {
        let counter = KmerCounter::new(2);
        let sequences = vec![
            b"ATCG".to_vec(),
            b"CGAT".to_vec(),
        ];
        
        counter.count_sequences(sequences.into_par_iter()).unwrap();
        
        let counts = counter.get_counts();
        assert_eq!(counts.get("AT").unwrap(), &2);
        assert_eq!(counts.get("TC").unwrap(), &1);
        assert_eq!(counts.get("CG").unwrap(), &2);
        assert_eq!(counts.get("GA").unwrap(), &1);
    }
}