use std::path::PathBuf;
use std::fs::File;
use std::io::Write;
use anyhow::Result;
use crate::profile::types::ProfileMatch;
use crate::profile::analyzer::ProfileAnalyzer;
use crate::kmer::KmerCounter;

pub fn output_analysis(
    matches: &[ProfileMatch],
    detailed: bool,
    analyzer: &ProfileAnalyzer,
    counter: &KmerCounter,
    sample_info_path: &PathBuf,
    matches_path: &PathBuf
) -> Result<()> {
    // Write sample information
    let mut sample_writer = File::create(sample_info_path)?;
    writeln!(sample_writer, "{:<30}\t{}", "Metric", "Value")?;
    writeln!(sample_writer, "{}", "-".repeat(50))?;
    writeln!(sample_writer, "{:<30}\t{}", "Total k-mers", counter.total_kmers())?;
    writeln!(sample_writer, "{:<30}\t{}", "Unique k-mers", counter.unique_kmers())?;
    writeln!(sample_writer, "{:<30}\t{}", "K-mer size", counter.kmer_size())?;

    // Write matches summary
    let mut matches_writer = File::create(matches_path)?;
    writeln!(matches_writer, "{:<40}\t{:>10}\t{:>10}\t{:>10}\t{:>10}\t{:>10}",
        "Name", "Sample%", "Profile%", "Shared", "Unique", "Jaccard%")?;
    writeln!(matches_writer, "{}", "-".repeat(100))?;

    for m in matches {
        let profile_size = analyzer.get_profile_kmer_count(m.name.clone())?;
        if profile_size < 0 {
            return Err(anyhow::anyhow!("Invalid negative kmer count in database"));
        }
        
        let sample_coverage = (m.shared_kmers as f64 / counter.unique_kmers() as f64) * 100.0;
        let profile_coverage = (m.shared_kmers as f64 / profile_size as f64) * 100.0;
        let profile_size = profile_size as usize;
        let union_size = counter.unique_kmers() + profile_size - m.shared_kmers;
        let jaccard = (m.shared_kmers as f64 / union_size as f64) * 100.0;
    
        writeln!(matches_writer, "{:<40}\t{:>10.2}\t{:>10.2}\t{:>10}\t{:>10}\t{:>10.2}",
            m.name,
            sample_coverage,
            profile_coverage, 
            m.shared_kmers,
            m.unique_matches,
            jaccard,
        )?;
    }

    // Write detailed analysis if requested
    if detailed {
        let stem = matches_path.file_stem().unwrap().to_string_lossy();
        let detailed_path = PathBuf::from(format!("{}_detailed.tsv", stem));
        let mut detailed_writer = File::create(&detailed_path)?;
        
        for m in matches {
            if let Some(analysis) = analyzer.get_detailed_analysis(counter, &m.name)? {
                writeln!(detailed_writer, "Profile: {}", m.name)?;
                writeln!(detailed_writer, "{}", "-".repeat(75))?;
                
                // Profile statistics
                writeln!(detailed_writer, "Statistic\tValue")?;
                writeln!(detailed_writer, "Total shared k-mers\t{}", analysis.statistics.total_shared)?;
                writeln!(detailed_writer, "Unique to reference\t{}", analysis.statistics.total_unique_reference)?;
                writeln!(detailed_writer, "Unique to sample\t{}", analysis.statistics.total_unique_sample)?;
                writeln!(detailed_writer, "Avg frequency diff\t{:.6}", analysis.statistics.average_frequency_difference)?;
                writeln!(detailed_writer, "Marker matches\t{}", analysis.statistics.marker_kmer_matches)?;
                writeln!(detailed_writer)?;

                // Shared k-mers
                writeln!(detailed_writer, "Shared K-mers\tSample%\tRef%\tRatio")?;
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
                writeln!(detailed_writer)?;
            }
        }
    }

    Ok(())
}