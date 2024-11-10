paptro-rusty/
├── Cargo.toml
├── src/
│   ├── main.rs         # Entry point
│   ├── cli.rs          # CLI argument handling
│   ├── config.rs       # Configuration management
│   ├── kmer/
│   │   ├── mod.rs      # K-mer module exports
│   │   ├── counter.rs  # K-mer counting logic
│   │   ├── profile.rs  # Profile generation
│   │   └── types.rs    # K-mer related types
│   ├── io/
│   │   ├── mod.rs      # I/O module exports
│   │   ├── reader.rs   # FASTA/FASTQ file reading
│   │   └── writer.rs   # Output writing
│   └── error.rs        # Error types
├── tests/              # Integration tests
└── benches/           # Benchmarks
Key Features and Components:

Input Processing:

Support for single or multiple FASTA/FASTQ files
Parallel processing support
Streaming for memory efficiency


K-mer Analysis:

Configurable k-mer size
Memory-efficient k-mer storage
Canonical k-mer representation


Profile Generation:

K-mer frequency counting
Position tracking (optional)
Profile comparison capabilities


Output Options:

Standard k-mer counts
Profile visualization
Machine-readable formats (JSON/CSV)


Performance Features:

Parallel processing
Memory-efficient algorithms
Progress reporting