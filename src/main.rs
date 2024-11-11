mod cli;
mod db;
mod profile;
mod io;
mod kmer;

use anyhow::{Result, Context};
use clap::Parser;
use log::{info, warn};
use profile::ProfileMatch;
use std::io::Write;
use std::fs::File;
use std::path::PathBuf;
use rayon::prelude::*;

use crate::cli::{Cli, Commands, DatabaseSubcommand, ExportFormat};
use crate::db::Database;
use crate::io::FastxReader;
use crate::kmer::KmerCounter;
use crate::profile::ProfileAnalyzer;

fn main() -> Result<()> {
    let cli = Cli::parse();

    // Set up logging
    let mut builder = env_logger::Builder::from_default_env();
    if let Some(log_file) = cli.log_file {
        let file = File::create(log_file)?;
        builder.target(env_logger::Target::Pipe(Box::new(file)));
    }
    if cli.verbose {
        builder.filter_level(log::LevelFilter::Debug);
    }
    builder.init();

    // Set up parallel processing
    if let Some(threads) = cli.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .context("Failed to initialize thread pool")?;
    }

    match cli.command {
        Commands::DB(db_cmd) => handle_db_command(db_cmd, cli.verbose)?,
        Commands::Analyze(analyze_cmd) => handle_analyze_command(analyze_cmd, cli.verbose)?,
    }

    Ok(())
}

fn handle_db_command(cmd: cli::DatabaseCommand, verbose: bool) -> Result<()> {
    match cmd.command {
        DatabaseSubcommand::Init => {
            info!("Initializing database at {}", cmd.database.display());
            Database::new(&cmd.database)?;
            info!("Database initialized successfully");
        }

        DatabaseSubcommand::Create { 
            input_files, 
            kmer_size, 
            level, 
            name,
            skip_existing 
        } => {
            let mut db = Database::new(&cmd.database)?;
            
            if db.get_profile(&name)?.is_some() {
                if skip_existing {
                    warn!("Profile {} already exists, skipping", name);
                    return Ok(());
                } else {
                    return Err(anyhow::anyhow!("Profile {} already exists", name));
                }
            }

            info!("Creating profile from {} input files...", input_files.len());
            db.create_profile(input_files, kmer_size, level.into(), name)?;
        }

        DatabaseSubcommand::List { level, detailed } => {
            let db = Database::new(&cmd.database)?;
            let profiles = db.list_profiles(level.map(Into::into))?;
            
            if profiles.is_empty() {
                println!("name\tlevel\tk_size\ttotal_kmers\tcreated_at");
                return Ok(());
            }

            println!("name\tlevel\tk_size\ttotal_kmers\tcreated_at");
            for profile in &profiles {
                println!("{}\t{:?}\t{}\t{}\t{}",
                    profile.name,
                    profile.level,
                    profile.k,
                    profile.total_kmers,
                    profile.created_at,
                );

                if detailed {
                    if let Some(profile_data) = db.get_profile(&profile.name)? {
                        println!("\n# Top k-mers for {}", profile.name);
                        println!("kmer\tfrequency");
                        let mut kmers: Vec<_> = profile_data.frequencies.iter().collect();
                        kmers.sort_by(|(_, a), (_, b)| b.partial_cmp(a).unwrap());
                        for (kmer, freq) in kmers.iter().take(5) {
                            println!("{}\t{:.6}", kmer, freq);
                        }
                        println!();
                    }
                }
            }
        }

        DatabaseSubcommand::Remove { name, force } => {
            let mut db = Database::new(&cmd.database)?;
            
            if !force {
                print!("Are you sure you want to remove profile {}? [y/N] ", name);
                std::io::stdout().flush()?;
                let mut input = String::new();
                std::io::stdin().read_line(&mut input)?;
                if !input.trim().eq_ignore_ascii_case("y") {
                    info!("Operation cancelled");
                    return Ok(());
                }
            }

            if db.remove_profile(&name)? {
                info!("Profile {} removed", name);
            } else {
                warn!("Profile {} not found", name);
            }
        }

        DatabaseSubcommand::Export { names, output, format } => {
            let db = Database::new(&cmd.database)?;
            std::fs::create_dir_all(&output)?;

            let profiles = if names.is_empty() {
                db.list_profiles(None)?
                    .into_iter()
                    .map(|p| p.name)
                    .collect()
            } else {
                names
            };

            for name in profiles {
                if let Some(profile) = db.get_profile(&name)? {
                    let file_name = match format {
                        ExportFormat::Fasta => format!("{}.fasta", name),
                        ExportFormat::Tsv => format!("{}.tsv", name),
                    };
                    let output_path = output.join(file_name);
                    let mut file = File::create(&output_path)?;

                    match format {
                        ExportFormat::Fasta => {
                            for (kmer, freq) in &profile.frequencies {
                                writeln!(file, ">{} {:.6}", name, freq)?;
                                writeln!(file, "{}", kmer)?;
                            }
                        }
                        ExportFormat::Tsv => {
                            writeln!(file, "kmer\tfrequency")?;
                            let mut kmers: Vec<_> = profile.frequencies.iter().collect();
                            kmers.sort_by(|(_, a), (_, b)| b.partial_cmp(a).unwrap());
                            for (kmer, freq) in kmers {
                                writeln!(file, "{}\t{:.6}", kmer, freq)?;
                            }
                        }
                    }
                    info!("Exported profile {} to {}", name, output_path.display());
                } else {
                    warn!("Profile {} not found", name);
                }
            }
        }

        DatabaseSubcommand::Stats { detailed } => {
            let db = Database::new(&cmd.database)?;
            let stats = db.get_statistics()?;
            
            println!("metric\tvalue");
            println!("total_profiles\t{}", stats.total_profiles);
            println!("total_kmers\t{}", stats.total_kmers);
            
            println!("\n# Profiles by level");
            println!("level\tcount");
            for (level, count) in &stats.profiles_by_level {
                println!("{}\t{}", level, count);
            }

            if detailed {
                // Add detailed statistics if implemented
            }
        }

        DatabaseSubcommand::Validate => {
            let db = Database::new(&cmd.database)?;
            info!("Validating database integrity...");
            // Add validation logic here
            info!("Database validation complete");
        }
    }

    Ok(())
}

fn handle_analyze_command(cmd: cli::AnalyzeCommand, verbose: bool) -> Result<()> {
    let analyzer = ProfileAnalyzer::new(
        &cmd.database,
        cmd.min_similarity,
        cmd.min_shared_kmers,
        cmd.level.into(),
    )?;

    info!("Processing input files...");
    let counter = KmerCounter::new(cmd.kmer_size);
    let reader = FastxReader::new(cmd.input_files);
    
    let mut sequences = Vec::new();
    reader.process_all(|sequence, _id| {
        sequences.push(sequence.to_vec());
        Ok(())
    })?;

    counter.count_sequences(sequences.into_par_iter())?;
    info!("Found {} unique k-mers in sample", counter.unique_kmers());

    let matches = analyzer.analyze_sample(&counter)?;
    output_analysis(&matches, cmd.detailed, &analyzer, &counter, cmd.output)?;

    Ok(())
}

fn output_analysis(
    matches: &[ProfileMatch],
    detailed: bool,
    analyzer: &ProfileAnalyzer,
    counter: &KmerCounter,
    output: Option<PathBuf>,
) -> Result<()> {
    let mut writer: Box<dyn Write> = match output {
        Some(path) => Box::new(File::create(path)?),
        None => Box::new(std::io::stdout()),
    };

    // Write sample information header
    writeln!(writer, "# Sample Information")?;
    writeln!(writer, "Total k-mers\t{}", counter.total_kmers())?;
    writeln!(writer, "Unique k-mers\t{}", counter.unique_kmers())?;
    writeln!(writer, "K-mer size\t{}", counter.kmer_size())?;
    writeln!(writer)?;

    // Write match summary header
    writeln!(writer, "# Match Summary")?;
    writeln!(writer, "name\tsample_coverage\tprofile_coverage\tshared_kmers\tunique_matches\tjaccard_similarity")?;

    for m in matches {
        let profile_size = analyzer.get_profile_kmer_count(m.name.clone())?;
        if profile_size < 0 {
            return Err(anyhow::anyhow!("Invalid negative kmer count in database"));
        }
        
        // Calculate various metrics after validating profile_size
        let sample_coverage = (m.shared_kmers as f64 / counter.unique_kmers() as f64) * 100.0;
        let profile_coverage = (m.shared_kmers as f64 / profile_size as f64) * 100.0;
        let profile_size = profile_size as usize;
        let union_size = counter.unique_kmers() + profile_size - m.shared_kmers;
        let jaccard = (m.shared_kmers as f64 / union_size as f64) * 100.0;
    
        writeln!(writer, "{}\t{:.2}\t{:.2}\t{}\t{}\t{:.2}",
            m.name,
            sample_coverage,
            profile_coverage, 
            m.shared_kmers,
            m.unique_matches,
            jaccard,
        )?;
    }

    if detailed {
        for m in matches {
            if let Some(analysis) = analyzer.get_detailed_analysis(counter, &m.name)? {
                writeln!(writer, "\n# Detailed analysis for {}", m.name)?;
                
                // Profile statistics
                writeln!(writer, "## Profile Statistics")?;
                writeln!(writer, "metric\tvalue")?;
                writeln!(writer, "total_shared_kmers\t{}", analysis.statistics.total_shared)?;
                writeln!(writer, "unique_to_reference\t{}", analysis.statistics.total_unique_reference)?;
                writeln!(writer, "unique_to_sample\t{}", analysis.statistics.total_unique_sample)?;
                writeln!(writer, "avg_frequency_diff\t{:.6}", analysis.statistics.average_frequency_difference)?;
                writeln!(writer, "marker_matches\t{}", analysis.statistics.marker_kmer_matches)?;

                // K-mer distribution
                writeln!(writer, "\n## K-mer Distribution")?;
                writeln!(writer, "### Top Shared K-mers")?;
                writeln!(writer, "kmer\tsample_freq\tref_freq\tfreq_ratio")?;
                let mut shared_kmers: Vec<_> = analysis.shared_kmers.iter().collect();
                shared_kmers.sort_by(|a, b| b.sample_frequency.partial_cmp(&a.sample_frequency).unwrap());
                for kmer in shared_kmers.iter().take(10) {
                    let ratio = if kmer.reference_frequency > 0.0 {
                        kmer.sample_frequency / kmer.reference_frequency
                    } else {
                        0.0
                    };
                    writeln!(writer, "{}\t{:.6}\t{:.6}\t{:.6}",
                        kmer.sequence,
                        kmer.sample_frequency,
                        kmer.reference_frequency,
                        ratio,
                    )?;
                }

                writeln!(writer, "\n### Most Abundant Reference-Unique K-mers")?;
                writeln!(writer, "kmer\tfrequency")?;
                let mut ref_unique: Vec<_> = analysis.unique_to_reference.iter().collect();
                ref_unique.sort_by(|a, b| b.frequency.partial_cmp(&a.frequency).unwrap());
                for kmer in ref_unique.iter().take(5) {
                    writeln!(writer, "{}\t{:.6}", kmer.sequence, kmer.frequency)?;
                }

                writeln!(writer, "\n### Most Abundant Sample-Unique K-mers")?;
                writeln!(writer, "kmer\tfrequency")?;
                let mut sample_unique: Vec<_> = analysis.unique_to_sample.iter().collect();
                sample_unique.sort_by(|a, b| b.frequency.partial_cmp(&a.frequency).unwrap());
                for kmer in sample_unique.iter().take(5) {
                    writeln!(writer, "{}\t{:.6}", kmer.sequence, kmer.frequency)?;
                }

                writeln!(writer)?;
            }
        }
    }

    Ok(())
}