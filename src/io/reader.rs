use std::path::{Path, PathBuf};
use anyhow::{Result, Context};
use needletail::{parse_fastx_file, Sequence};
use log::{info, warn};

/// Represents a FASTA/FASTQ sequence reader that can handle multiple files
pub struct FastxReader {
    files: Vec<PathBuf>
}

impl FastxReader {
    /// Create a new FastxReader from a list of input files
    pub fn new<P: AsRef<Path>>(files: Vec<P>) -> Self {
        let files = files.into_iter()
            .map(|p| p.as_ref().to_owned())
            .collect();
        
        FastxReader {
            files
        }
    }

    /// Process each sequence in all input files
    pub fn process_all<F>(&self, mut callback: F) -> Result<()>
    where
        F: FnMut(&[u8], &str) -> Result<()>
    {
        for file in &self.files {
            self.process_file(file, &mut callback)
                .with_context(|| format!("Failed to process file: {}", file.display()))?;
        }
        Ok(())
    }

    /// Process a single FASTA/FASTQ file
    fn process_file<F>(&self, path: &Path, callback: &mut F) -> Result<()>
    where
        F: FnMut(&[u8], &str) -> Result<()>
    {
        info!("Processing file: {}", path.display());
        
        let mut reader = parse_fastx_file(path)
            .with_context(|| format!("Failed to open file: {}", path.display()))?;
        
        let mut num_sequences = 0;
        let mut num_invalid = 0;

        while let Some(record) = reader.next() {
            let record = record.with_context(|| "Failed to parse sequence record")?;
            
            // Normalize sequence to uppercase and process
            let sequence = record.normalize(false);
            let id = String::from_utf8_lossy(record.id());
            
            // Check for invalid characters (non-ACGT)
            if sequence.iter().any(|&b| !matches!(b, b'A' | b'C' | b'G' | b'T')) {
                num_invalid += 1;
                continue;
            }

            callback(&sequence, &id)?;
            num_sequences += 1;
        }

        info!("Processed {} sequences from {}", num_sequences, path.display());
        if num_invalid > 0 {
            warn!("Skipped {} sequences containing invalid characters", num_invalid);
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use std::io::Write;
    use tempfile::tempdir;

    #[test]
    fn test_process_valid_fasta() -> Result<()> {
        // Create a temporary directory and fasta file
        let dir = tempdir()?;
        let file_path = dir.path().join("test.fasta");
        let mut file = File::create(&file_path)?;

        // Write test data
        writeln!(file, ">seq1\nACGT\n>seq2\nGTCA")?;

        let reader = FastxReader::new(vec![file_path]);
        let mut sequences = Vec::new();
        let mut ids = Vec::new();

        reader.process_all(|seq, id| {
            sequences.push(seq.to_vec());
            ids.push(id.to_string());
            Ok(())
        })?;

        assert_eq!(sequences.len(), 2);
        assert_eq!(ids.len(), 2);
        assert_eq!(sequences[0], b"ACGT");
        assert_eq!(sequences[1], b"GTCA");
        assert_eq!(ids[0], "seq1");
        assert_eq!(ids[1], "seq2");

        Ok(())
    }

    #[test]
    fn test_process_invalid_sequences() -> Result<()> {
        // Create a temporary directory and fasta file
        let dir = tempdir()?;
        let file_path = dir.path().join("test.fasta");
        let mut file = File::create(&file_path)?;

        // Write test data with invalid sequences
        writeln!(file, ">seq1\nACGT\n>seq2\nNNNN\n>seq3\nGTCA")?;

        let reader = FastxReader::new(vec![file_path]);
        let mut sequences = Vec::new();

        reader.process_all(|seq, _id| {
            sequences.push(seq.to_vec());
            Ok(())
        })?;

        assert_eq!(sequences.len(), 2);
        assert_eq!(sequences[0], b"ACGT");
        assert_eq!(sequences[1], b"GTCA");

        Ok(())
    }
}