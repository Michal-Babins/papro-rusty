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
    
        writeln!(matches_writer, "{:<40}\t{:<40}\t{:>10}\t{:>10}\t{:>10}\t{:>10}\t{:>10}\t{:>10}",
            "Name", "Sample", "Sample%", "Shared", "Core", "Rare", "Avg_Freq", "Size_Ratio")?;
        writeln!(matches_writer, "{}", "-".repeat(140))?;
    }

    // Write sample information
    writeln!(sample_writer, "\nSample: {}", sample_name)?;
    writeln!(sample_writer, "{:<30}\t{}", "Total k-mers", counter.total_kmers())?;
    writeln!(sample_writer, "{:<30}\t{}", "Unique k-mers", counter.unique_kmers())?;
    writeln!(sample_writer, "{:<30}\t{}", "K-mer size", counter.kmer_size())?;

    // Write matches for this sample
    for m in matches {
        let profile_size = analyzer.get_profile_kmer_count(m.name.clone())?;
        if profile_size < 0 {
            return Err(anyhow::anyhow!("Invalid negative kmer count in database"));
        }
        
        if let Some(analysis) = analyzer.get_detailed_analysis(counter, &m.name)? {
            let sample_coverage = (m.shared_kmers as f64 / counter.unique_kmers() as f64) * 100.0;
            
            // Calculate new metrics
            let core_matches = analysis.shared_kmers.iter()
                .filter(|kmer| kmer.reference_frequency >= 0.95)
                .count();
                
            let rare_matches = analysis.shared_kmers.iter()
                .filter(|kmer| kmer.reference_frequency <= 0.10)
                .count();
                
            let avg_frequency = if !analysis.shared_kmers.is_empty() {
                analysis.shared_kmers.iter()
                    .map(|k| k.reference_frequency)
                    .sum::<f64>() / analysis.shared_kmers.len() as f64
            } else {
                0.0
            };
            
            let size_ratio = counter.unique_kmers() as f64 / profile_size as f64;

            writeln!(matches_writer, "{:<40}\t{:<40}\t{:>10.2}\t{:>10}\t{:>10}\t{:>10}\t{:>10.3}\t{:>10.3}",
                m.name,
                sample_name,
                sample_coverage,
                m.shared_kmers,
                core_matches,
                rare_matches,
                avg_frequency,
                size_ratio,
            )?;

            // Write detailed analysis if requested
            if detailed {
                let detailed_path = PathBuf::from(format!("{}_{}_detailed.tsv", sample_name, m.name));
                let mut detailed_writer = File::create(detailed_path)?;
                
                writeln!(detailed_writer, "Profile: {}", m.name)?;
                writeln!(detailed_writer, "{}", "-".repeat(75))?;
                
                // Profile statistics
                writeln!(detailed_writer, "Statistic\tValue")?;
                writeln!(detailed_writer, "Total shared k-mers\t{}", analysis.statistics.total_shared)?;
                writeln!(detailed_writer, "Unique to reference\t{}", analysis.statistics.total_unique_reference)?;
                writeln!(detailed_writer, "Unique to sample\t{}", analysis.statistics.total_unique_sample)?;
                writeln!(detailed_writer, "Avg frequency diff\t{:.6}", analysis.statistics.average_frequency_difference)?;
                writeln!(detailed_writer, "Core matches\t{}", core_matches)?;
                writeln!(detailed_writer, "Rare matches\t{}", rare_matches)?;
                writeln!(detailed_writer, "Average frequency\t{:.6}", avg_frequency)?;
                writeln!(detailed_writer, "Size ratio\t{:.6}", size_ratio)?;
                writeln!(detailed_writer)?;

                // K-mer frequency distribution
                let high_freq = analysis.shared_kmers.iter()
                    .filter(|k| k.reference_frequency >= 0.75)
                    .count();
                let mid_freq = analysis.shared_kmers.iter()
                    .filter(|k| k.reference_frequency >= 0.25 && k.reference_frequency < 0.75)
                    .count();
                let low_freq = analysis.shared_kmers.iter()
                    .filter(|k| k.reference_frequency < 0.25)
                    .count();

                writeln!(detailed_writer, "\nK-mer Frequency Distribution")?;
                writeln!(detailed_writer, "High frequency (>=75%)\t{}", high_freq)?;
                writeln!(detailed_writer, "Medium frequency (25-75%)\t{}", mid_freq)?;
                writeln!(detailed_writer, "Low frequency (<25%)\t{}", low_freq)?;
                writeln!(detailed_writer)?;

                // Top shared k-mers
                writeln!(detailed_writer, "\nTop Shared K-mers")?;
                writeln!(detailed_writer, "K-mer\tSample%\tRef%\tRatio")?;
                let mut shared_kmers: Vec<_> = analysis.shared_kmers.iter().collect();
                shared_kmers.sort_by(|a, b| b.sample_frequency.partial_cmp(&a.sample_frequency).unwrap());
                for kmer in shared_kmers.iter().take(10) {
                    let ratio = kmer.sample_frequency / kmer.reference_frequency.max(f64::EPSILON);
                    writeln!(detailed_writer, "{}\t{:.6}\t{:.6}\t{:.6}",
                        kmer.sequence,
                        kmer.sample_frequency * 100.0,
                        kmer.reference_frequency * 100.0,
                        ratio
                    )?;
                }
            }
        }
    }

    Ok(())
}