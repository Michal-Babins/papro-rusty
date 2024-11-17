use std::path::PathBuf;
use std::fs::File;
use std::io::{Seek, Write};
use anyhow::Result;
use crate::profile::types::ProfileMatch;
use crate::profile::analyzer::ProfileAnalyzer;
use crate::kmer::KmerCounter;

pub fn output_analysis(
    sample_name: &str,
    counter: &KmerCounter,
    matches: &[ProfileMatch],
    detailed: bool,
    analyzer: &ProfileAnalyzer,
    sample_writer: &mut impl Write,
    matches_writer: &mut (impl Write + Seek),
) -> Result<()> {
    // Check if we need to write headers (if file is empty)
    if matches_writer.stream_position()? == 0 {
        writeln!(sample_writer, "{:<30}\t{}", "Metric", "Value")?;
        writeln!(sample_writer, "{}", "-".repeat(50))?;
    
        writeln!(matches_writer, "{:<40}\t{:<40}\t{:>10}\t{:>10}\t{:>10}\t{:>10}\t{:>10}",
            "Name", "Sample", "Sample%", "Shared", "Unique%", "Size", "Confidence")?;
        writeln!(matches_writer, "{}", "-".repeat(140))?;
    }

    // Write sample information
    writeln!(sample_writer, "\nSample: {}", sample_name)?;
    writeln!(sample_writer, "{:<30}\t{}", "Total k-mers", counter.total_kmers())?;
    writeln!(sample_writer, "{:<30}\t{}", "Unique k-mers", counter.unique_kmers())?;
    writeln!(sample_writer, "{:<30}\t{}", "K-mer size", counter.kmer_size())?;

    // Write matches for this sample
    for m in matches {
        writeln!(matches_writer, "{:<40}\t{:<40}\t{:>10.2}\t{:>10}\t{:>10.2}\t{:>10.3}\t{:>10.3}",
            m.name,
            sample_name,
            m.sample_coverage * 100.0,
            m.shared_kmers,
            m.uniqueness_score * 100.0,
            m.size_ratio,
            m.confidence_score,
        )?;

        // Write detailed analysis if requested
        if detailed {
            if let Some(analysis) = analyzer.get_detailed_analysis(counter, &m.name)? {
                let detailed_path = PathBuf::from(format!("{}_{}_detailed.tsv", sample_name, m.name));
                let mut detailed_writer = File::create(detailed_path)?;
                
                writeln!(detailed_writer, "Profile: {}", m.name)?;
                writeln!(detailed_writer, "{}", "-".repeat(75))?;
                
                // Profile statistics
                writeln!(detailed_writer, "Statistic\tValue")?;
                writeln!(detailed_writer, "Total shared k-mers\t{}", analysis.statistics.total_shared)?;
                writeln!(detailed_writer, "Total unique to reference\t{}", analysis.statistics.total_unique_reference)?;
                writeln!(detailed_writer, "Total unique to sample\t{}", analysis.statistics.total_unique_sample)?;
                writeln!(detailed_writer, "Sample coverage\t{:.6}", analysis.statistics.sample_coverage)?;
                writeln!(detailed_writer, "Uniqueness score\t{:.6}", analysis.statistics.uniqueness_score)?;
                writeln!(detailed_writer, "Size ratio\t{:.6}", analysis.statistics.size_ratio)?;
                writeln!(detailed_writer, "Confidence score\t{:.6}", analysis.statistics.confidence_score)?;
                writeln!(detailed_writer, "Profile unique k-mers\t{}", analysis.statistics.profile_unique_kmers)?;
                writeln!(detailed_writer, "Shared unique k-mers\t{}", analysis.statistics.shared_unique_kmers)?;
                writeln!(detailed_writer)?;

                // Top shared k-mers
                writeln!(detailed_writer, "\nTop Shared K-mers")?;
                writeln!(detailed_writer, "K-mer\tSample%\tUnique\tFrequency")?;
                let mut shared_kmers: Vec<_> = analysis.shared_kmers.iter().collect();
                shared_kmers.sort_by(|a, b| b.sample_frequency.partial_cmp(&a.sample_frequency).unwrap());
                for kmer in shared_kmers.iter().take(10) {
                    writeln!(detailed_writer, "{}\t{:.6}\t{}\t{:.6}",
                        kmer.sequence,
                        kmer.sample_frequency * 100.0,
                        if kmer.is_unique { "Yes" } else { "No" },
                        kmer.sample_frequency
                    )?;
                }
            }
        }
    }

    Ok(())
}