# 🧬 PAPRO-RUSTY
## A Pathothogen based K-mer Profiler written in Rust (beta)
### Package Under Counstruction 🚧
[![Rust](https://github.com/Michal-Babins/papro-rusty/workflows/Rust/badge.svg)](https://github.com/Michal-Babins/papro-rusty/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A lightweight k-mer based tool for pathogen profiling and taxonomic classification. This tool allows you to create, manage, and analyze k-mer profiles from genomic sequences for pathogen identification and characterization.

## ✨ Features

### 🔬 Profile Management
- Create k-mer profiles from FASTA/FASTQ files
- Support for multiple taxonomic levels (Genus, Species, Strain)
- Efficient SQLite-based profile storage
- Profile import/export capabilities
- Built-in profile statistics and analysis

### 🧪 Sample Analysis
- Fast k-mer counting and profile comparison
- Parallel processing support
- Multiple similarity metrics
- Detailed analysis reports
- Flexible output options

### 🚀 Performance
- Optimized k-mer counting algorithms
- Memory-efficient data structures
- Multi-threaded processing
- Efficient database queries

## 📦 Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/papro-rusty.git
cd papro-rusty

# Build the project
cargo build --release

# Run tests
cargo test
```

## 🚀 Quick Start

### Database Management

```bash
# Initialize a new profile database
papro-rusty db -d profiles.db init

# Add a species-level profile from a single file
papro-rusty db -d profiles.db add \
    -l species \
    -n "e_coli" \
    -k 21 \
    ecoli_genome.fasta

# List all profiles
papro-rusty db -d profiles.db list

# Show database statistics
papro-rusty db -d profiles.db stats
```

### Sample Analysis

```bash
# Analyze a sample against the database
papro-rusty analyze \
    -d profiles.db \
    -l species \
    --min-similarity 0.7 \
    --min-shared-kmers 100 \
    sample.fasta

# Detailed analysis with verbose output
papro-rusty -v analyze \
    -d profiles.db \
    -o results.txt \
    sample.fasta
```

## 📋 Command Reference

### Global Options
```
-v, --verbose    Enable verbose output
-t, --threads    Specify number of threads (default: all available)
```

### Database Commands
```bash
# Initialize database
db init

# Add profile
db add [options] <files>...
  -l, --level <LEVEL>     Taxonomic level (genus|species|strain)
  -n, --name <NAME>       Profile name
  -k, --kmer-size <SIZE>  K-mer size (default: 21)

# List profiles
db list [options]
  -l, --level <LEVEL>     Filter by taxonomic level

# Remove profile
db remove <name>

# Export profile
db export -o <file> <name>

# Show statistics
db stats
```

### Analysis Commands
```bash
analyze [options] <files>...
  -d, --database <FILE>         Reference database
  -l, --level <LEVEL>           Taxonomic level
  --min-similarity <FLOAT>      Minimum similarity score (0.0-1.0)
  --min-shared-kmers <INT>      Minimum shared k-mers
  -o, --output <FILE>           Output file for detailed results
```

## 📊 Output Format

### Profile List
```
Profiles in database:
Name                 Level      K-mer size     Total k-mers    Created
--------------------------------------------------------------------------------
Escherichia_coli             Species    31             1234567         2024-03-09 10:30:45
Salmonella_enterica          Species    31             987654          2024-03-09 11:15:30
```

### Analysis Results
```
Results at Species level:
Name                                     Similarity      Shared     Unique
---------------------------------------------------------------------------
Escherichia_coli                        0.945          45678      789
Salmonella_enterica                     0.432          12345      234
```

## 🧪 Testing

Run the test suite:
```bash
# Run all tests
cargo test

# Run tests with output
cargo test -- --nocapture

# Run specific test
cargo test test_name
```

## 🛠 Development

### Project Structure
```
src/
├── main.rs           # Entry point
├── cli.rs           # CLI implementation
├── db/              # Database management
│   ├── mod.rs
│   ├── database.rs
│   ├── schema.rs
│   └── types.rs
├── profile/         # Profile handling
│   ├── mod.rs
│   ├── analyzer.rs
│   └── types.rs
├── io/             # Input/Output handling
│   ├── mod.rs
│   ├── reader.rs
│   └── writer.rs
└── kmer/           # K-mer processing
    ├── mod.rs
    └── counter.rs
```

### Adding New Features
1. Create feature branch
2. Implement changes
3. Add tests
4. Update documentation
5. Submit pull request

## 📚 Technical Details

### K-mer Processing
- Canonical k-mer representation
- Rolling hash implementation
- Memory-efficient storage
- Parallel processing support

### Database Schema
```sql
CREATE TABLE profiles (
    id INTEGER PRIMARY KEY,
    name TEXT NOT NULL UNIQUE,
    taxonomy_level TEXT NOT NULL,
    k INTEGER NOT NULL,
    total_kmers INTEGER NOT NULL,
    created_at DATETIME DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE kmers (
    profile_id INTEGER,
    kmer TEXT NOT NULL,
    frequency REAL NOT NULL,
    FOREIGN KEY(profile_id) REFERENCES profiles(id),
    PRIMARY KEY(profile_id, kmer)
);
```

## 📝 Contributing

Contributions are welcome! No formal guide yet.

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- [Rust Bioinformatics Working Group](https://github.com/rust-bio)
- [Rusqlite](https://github.com/rusqlite/rusqlite)
- [Clap](https://github.com/clap-rs/clap)
- [Rayon](https://github.com/rayon-rs/rayon)

## 📧 Contact

Your Name - Michal Babinski - mbabinski17@gmail.com

Project Link: [https://github.com/Michal-Babins/papro-rusty](https://github.com/Michal-Babins/papro-rusty)
