use clap::{Parser, Subcommand, ValueEnum};
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(author, version, about = "K-mer based pathogen profiling tool")]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,

    /// Enable verbose output
    #[arg(short, long, global = true)]
    pub verbose: bool,

    /// Number of threads to use
    #[arg(short, long, global = true)]
    pub threads: Option<usize>,

    /// Path to log file
    #[arg(long, global = true)]
    pub log_file: Option<PathBuf>,
}

#[derive(Subcommand, Debug)]
pub enum Commands {
    /// Create and manage reference profiles
    DB(DatabaseCommand),

    /// Analyze samples against reference profiles
    Analyze(AnalyzeCommand),
}

#[derive(Parser, Debug)]
pub struct DatabaseCommand {

    /// Path to the SQLite database file
    #[arg(short, long, default_value = "profiles.db")]
    pub database: PathBuf,

    #[command(subcommand)]
    pub command: DatabaseSubcommand,
}

#[derive(Subcommand, Debug)]
pub enum DatabaseSubcommand {

    /// Initialize a new database, maybe reconstruct his approach later
    Init,

    /// Create a new profile
    Create {
        /// Input FASTA/FASTQ files
        #[arg(required = true)]
        input_files: Vec<PathBuf>,

        /// K-mer size to use
        #[arg(short, long, default_value = "31")]
        kmer_size: usize,

        /// Taxonomic level
        #[arg(short, long, value_enum)]
        level: TaxonomyLevel,

        /// Name of the organism (e.g., "e_coli")
        #[arg(short, long)]
        name: String,

        /// Skip existing files instead of erroring
        #[arg(long)]
        skip_existing: bool,
    },

    /// List profiles in database
    List {
        /// Filter by taxonomic level
        #[arg(short, long, value_enum)]
        level: Option<TaxonomyLevel>,

        /// Show detailed k-mer information
        #[arg(long)]
        detailed: bool,
    },

    /// Remove a profile
    Remove {
        /// Name of profile to remove
        name: String,

        /// Force removal without confirmation
        #[arg(short, long)]
        force: bool,
    },

    /// Export profiles
    Export {
        /// Names of profiles to export (exports all if none specified)
        names: Vec<String>,

        /// Output directory
        #[arg(short, long)]
        output: PathBuf,

        /// Export format (fasta or json)
        #[arg(short, long, value_enum, default_value = "fasta")]
        format: ExportFormat,
    },

    /// Show database statistics
    Stats {
        /// Include k-mer distribution statistics
        #[arg(long)]
        detailed: bool,
    },

    /// Validate database integrity
    Validate,
}

#[derive(Parser, Debug)]
pub struct AnalyzeCommand {
    /// Input FASTA/FASTQ files to analyze
    #[arg(required = true)]
    pub input_files: Vec<PathBuf>,

    /// Path to reference profile database
    #[arg(short, long)]
    pub database: PathBuf,

    /// K-mer size to use
    #[arg(short, long, default_value = "31")]
    pub kmer_size: usize,

    /// Taxonomic level to analyze
    #[arg(short, long, value_enum, default_value = "species")]
    pub level: TaxonomyLevel,

    /// Minimum similarity score (0.0-1.0)
    #[arg(long, default_value = "0.80")]
    pub min_similarity: f64,

    /// Minimum number of shared k-mers
    #[arg(long, default_value = "100")]
    pub min_shared_kmers: usize,

    /// Generate detailed report
    #[arg(long)]
    pub detailed: bool,

    /// Output format
    #[arg(short = 'f', long, value_enum, default_value = "text")]
    pub format: OutputFormat,

    /// Output file (defaults to stdout)
    #[arg(short, long)]
    pub output: Option<PathBuf>,
}

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum, Debug)]
pub enum TaxonomyLevel {
    Genus,
    Species,
    Strain,
}

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum, Debug)]
pub enum ExportFormat {
    Fasta,
    Tsv,
}

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum, Debug)]
pub enum OutputFormat {
    Text,
    Json,
    Tsv,
}

impl From<TaxonomyLevel> for crate::profile::TaxonomyLevel {
    fn from(level: TaxonomyLevel) -> Self {
        match level {
            TaxonomyLevel::Genus => Self::Genus,
            TaxonomyLevel::Species => Self::Species,
            TaxonomyLevel::Strain => Self::Strain,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn verify_cli() {
        Cli::command().debug_assert();
    }

    #[test]
    fn test_taxonomy_level_conversion() {
        assert!(matches!(
            TaxonomyLevel::Species.into(),
            crate::profile::TaxonomyLevel::Species
        ));
    }
}