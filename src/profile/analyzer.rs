use std::collections::{HashMap, HashSet};
use anyhow::{Result, Context};
use log::{debug, info, warn};
use rusqlite::{Connection, params, OptionalExtension};
use super::types::{ProfileMatch, TaxonomyLevel};
use crate::kmer::KmerCounter;

pub struct ProfileAnalyzer {
    conn: Connection,
    min_similarity: f64,
    min_shared_kmers: usize,
    taxonomy_level: TaxonomyLevel,
}

impl ProfileAnalyzer {
    pub fn new<P: AsRef<std::path::Path>>(
        database_path: P,
        min_similarity: f64,
        min_shared_kmers: usize,
        taxonomy_level: TaxonomyLevel,
    ) -> Result<Self> {
        let conn = Connection::open(database_path)
            .context("Failed to open database connection")?;
        
        Ok(ProfileAnalyzer {
            conn,
            min_similarity,
            min_shared_kmers,
            taxonomy_level,
        })
    }

    /// Analyze a sample against the database at the current taxonomy level
    pub fn analyze_sample(&self, counter: &KmerCounter) -> Result<Vec<ProfileMatch>> {
        info!(
            "Analyzing sample against reference profiles at {} level",
            self.taxonomy_level
        );
    
        let profile_count: i64 = self.conn.query_row(
            "SELECT COUNT(*) FROM profiles WHERE taxonomy_level = ?",
            params![self.taxonomy_level.to_string()],
            |row| row.get(0)
        )?;
    
        info!("Found {} profiles at {} level", profile_count, self.taxonomy_level);
    
        if profile_count == 0 {
            warn!("No profiles found at {} level in the database", self.taxonomy_level);
            return Ok(Vec::new());
        }
    
        let mut profile_stmt = self.conn.prepare(
            "SELECT id, name, k, total_kmers 
             FROM profiles 
             WHERE taxonomy_level = ?"
        )?;
    
        let sample_kmers = counter.get_counts();
        info!("Sample has {} unique k-mers of size {}", 
            sample_kmers.len(), counter.kmer_size());
    
        let mut matches = Vec::new();
        let profiles = profile_stmt.query_map(
            params![self.taxonomy_level.to_string()],
            |row| {
                Ok((
                    row.get::<_, i64>(0)?,
                    row.get::<_, String>(1)?,
                    row.get::<_, i64>(2)?,
                    row.get::<_, i64>(3)?,
                ))
            }
        )?;
    
        for profile_result in profiles {
            let (profile_id, name, k, total_kmers) = profile_result?;
            info!("Checking profile '{}' (id={}, k={}, total_kmers={})", 
                name, profile_id, k, total_kmers);
    
            if k as usize != counter.kmer_size() {
                warn!("K-mer size mismatch: profile {} has k={}, sample has k={}", 
                    name, k, counter.kmer_size());
                continue;
            }
    
            match self.compare_with_profile(
                profile_id,
                &name,
                &sample_kmers,
                counter.total_kmers(),
                total_kmers as usize,
            )? {
                Some(profile_match) => {
                    info!("Found match: {} (coverage={:.4}%, shared={}, core={}, rare={})",
                        name, 
                        profile_match.sample_coverage * 100.0,
                        profile_match.shared_kmers,
                        profile_match.core_matches,
                        profile_match.rare_matches
                    );
                    matches.push(profile_match);
                }
                None => {
                    info!("Profile {} did not meet thresholds (min_similarity={}, min_shared_kmers={})",
                        name, self.min_similarity, self.min_shared_kmers);
                }
            }
        }
    
        // Sort by sample coverage instead of similarity score
        matches.sort_by(|a, b| b.sample_coverage.partial_cmp(&a.sample_coverage).unwrap());
        info!("Found {} potential matches", matches.len());
        Ok(matches)
    }
    
    fn compare_with_profile(
        &self,
        profile_id: i64,
        profile_name: &str,
        sample_kmers: &HashMap<String, usize>,
        total_sample_kmers: usize,
        total_profile_kmers: usize,
    ) -> Result<Option<ProfileMatch>> {
        info!("Comparing profile {} (id={})", profile_name, profile_id);
    
        let mut kmer_stmt = self.conn.prepare(
            "SELECT kmer, frequency FROM kmers WHERE profile_id = ?"
        )?;
    
        let mut shared_kmers = 0;
        let mut core_matches = 0;
        let mut rare_matches = 0;
        let mut freq_sum = 0.0;
    
        let mut profile_unique_kmers = HashSet::new();
    
        for kmer_result in kmer_stmt.query_map(params![profile_id], |row| {
            Ok((
                row.get::<_, String>(0)?,
                row.get::<_, f64>(1)?,
            ))
        })? {
            let (kmer, ref_freq) = kmer_result?;
            profile_unique_kmers.insert(kmer.clone());
            
            if let Some(&_sample_count) = sample_kmers.get(&kmer) {
                shared_kmers += 1;
                freq_sum += ref_freq;
                
                if ref_freq >= 0.95 {
                    core_matches += 1;
                }
                if ref_freq <= 0.10 {
                    rare_matches += 1;
                }
            }
        }
    
        let sample_size = sample_kmers.len();
        let profile_size = profile_unique_kmers.len();
        let avg_frequency = if shared_kmers > 0 {
            freq_sum / shared_kmers as f64
        } else {
            0.0
        };
        let size_ratio = sample_size as f64 / profile_size as f64;
        let sample_coverage = shared_kmers as f64 / sample_size as f64;
    
        info!(
            "Comparison summary for {}:
            Shared k-mers: {}
            Core matches: {}
            Rare matches: {}
            Average frequency: {:.6}
            Size ratio: {:.6}
            Sample coverage: {:.6}",
            profile_name, 
            shared_kmers,
            core_matches,
            rare_matches,
            avg_frequency,
            size_ratio,
            sample_coverage
        );
    
        if sample_coverage >= self.min_similarity && shared_kmers >= self.min_shared_kmers {
            Ok(Some(ProfileMatch::new(
                profile_name.to_string(),
                sample_coverage,
                shared_kmers,
                core_matches,
                rare_matches,
                avg_frequency,
                size_ratio,
            )))
        } else {
            info!(
                "Profile {} did not meet thresholds:
                Sample coverage: {:.6} (minimum: {})
                Shared k-mers: {} (minimum: {})",
                profile_name, 
                sample_coverage, 
                self.min_similarity,
                shared_kmers, 
                self.min_shared_kmers
            );
            Ok(None)
        }
    }

    pub fn get_detailed_analysis(
        &self,
        counter: &KmerCounter,
        profile_name: &str,
    ) -> Result<Option<DetailedAnalysis>> {
        let profile_id: Option<i64> = self.conn.query_row(
            "SELECT id FROM profiles WHERE name = ?",
            params![profile_name],
            |row| row.get(0)
        ).optional()?;
    
        let Some(profile_id) = profile_id else {
            return Ok(None);
        };
    
        let mut kmer_stmt = self.conn.prepare(
            "SELECT kmer, frequency FROM kmers WHERE profile_id = ?"
        )?;
    
        let sample_kmers = counter.get_counts();
        let total_sample_kmers = counter.total_kmers() as f64;
    
        let mut analysis = DetailedAnalysis::new();
    
        // Get total profile k-mers for size ratio calculation
        let total_profile_kmers = self.get_profile_kmer_count(profile_name.to_string())?;
    
        for kmer_result in kmer_stmt.query_map(params![profile_id], |row| {
            Ok((
                row.get::<_, String>(0)?,
                row.get::<_, f64>(1)?,
            ))
        })? {
            let (kmer, ref_freq) = kmer_result?;
            
            if let Some(&sample_count) = sample_kmers.get(&kmer) {
                let sample_freq = sample_count as f64 / total_sample_kmers;
                analysis.add_shared_kmer(kmer, sample_freq, ref_freq);
            } else {
                analysis.add_reference_unique_kmer(kmer, ref_freq);
            }
        }
    
        // Add sample-unique k-mers
        for (kmer, count) in sample_kmers {
            if !analysis.has_kmer(&kmer) {
                let freq = count as f64 / total_sample_kmers;
                analysis.add_sample_unique_kmer(kmer.clone(), freq);
            }
        }
    
        analysis.calculate_statistics();
        Ok(Some(analysis))
    }

    pub fn get_profile_kmer_count(&self, name: String) -> Result<i64> {
        // Query total_kmers directly from profiles table and return error if not found
        let total_kmers: i64 = self.conn.query_row(
            "SELECT total_kmers FROM profiles WHERE name = ?",
            params![name],
            |row| row.get(0)
        )?;
    
        debug!("Profile '{}' has {} total k-mers from database", name, total_kmers);
    
        Ok(total_kmers)
    }
}


#[derive(Debug, Clone)]
pub struct SharedKmer {
    pub sequence: String,
    pub sample_frequency: f64,
    pub reference_frequency: f64,
}

#[derive(Debug, Clone)]
pub struct UniqueKmer {
    pub sequence: String,
    pub frequency: f64,
}

#[derive(Debug, Clone)]
pub struct DetailedAnalysis {
    pub shared_kmers: Vec<SharedKmer>,
    pub unique_to_reference: Vec<UniqueKmer>,
    pub unique_to_sample: Vec<UniqueKmer>,
    pub statistics: AnalysisStatistics,
}

#[derive(Debug, Clone)]
pub struct AnalysisStatistics {
    pub total_shared: usize,
    pub total_unique_reference: usize,
    pub total_unique_sample: usize,
    pub core_matches: usize,
    pub rare_matches: usize,
    pub avg_frequency: f64,
    pub size_ratio: f64,
    pub frequency_distribution: FrequencyDistribution,
    pub average_frequency_difference: f64,  // Added this back to maintain compatibility
    pub marker_kmer_matches: usize,        // Added this back to maintain compatibility
}

#[derive(Debug, Clone)]
pub struct FrequencyDistribution {
    pub high_freq: usize,   // >= 0.75
    pub mid_freq: usize,    // 0.25-0.75
    pub low_freq: usize,    // < 0.25
}

impl DetailedAnalysis {
    fn new() -> Self {
        DetailedAnalysis {
            shared_kmers: Vec::new(),
            unique_to_reference: Vec::new(),
            unique_to_sample: Vec::new(),
            statistics: AnalysisStatistics {
                total_shared: 0,
                total_unique_reference: 0,
                total_unique_sample: 0,
                core_matches: 0,
                rare_matches: 0,
                avg_frequency: 0.0,
                size_ratio: 0.0,
                frequency_distribution: FrequencyDistribution {
                    high_freq: 0,
                    mid_freq: 0,
                    low_freq: 0,
                },
                average_frequency_difference: 0.0,
                marker_kmer_matches: 0,
            },
        }
    }

    fn has_kmer(&self, kmer: &str) -> bool {
        self.shared_kmers.iter().any(|sk| sk.sequence == kmer)
    }

    fn add_shared_kmer(&mut self, sequence: String, sample_freq: f64, ref_freq: f64) {
        info!(
            "Adding shared k-mer: {} (sample_freq={:.6}, ref_freq={:.6})",
            sequence, sample_freq, ref_freq
        );
        self.shared_kmers.push(SharedKmer {
            sequence,
            sample_frequency: sample_freq,
            reference_frequency: ref_freq,
        });
    }

    fn add_reference_unique_kmer(&mut self, sequence: String, frequency: f64) {
        self.unique_to_reference.push(UniqueKmer {
            sequence,
            frequency,
        });
    }

    fn add_sample_unique_kmer(&mut self, sequence: String, frequency: f64) {
        self.unique_to_sample.push(UniqueKmer {
            sequence,
            frequency,
        });
    }

    fn calculate_statistics(&mut self) {
        self.statistics.total_shared = self.shared_kmers.len();
        self.statistics.total_unique_reference = self.unique_to_reference.len();
        self.statistics.total_unique_sample = self.unique_to_sample.len();

        // Calculate core and rare matches
        self.statistics.core_matches = self.shared_kmers.iter()
            .filter(|sk| sk.reference_frequency >= 0.95)
            .count();
        
        self.statistics.rare_matches = self.shared_kmers.iter()
            .filter(|sk| sk.reference_frequency <= 0.10)
            .count();

        // Calculate average frequency
        self.statistics.avg_frequency = if !self.shared_kmers.is_empty() {
            self.shared_kmers.iter()
                .map(|sk| sk.reference_frequency)
                .sum::<f64>() / self.shared_kmers.len() as f64
        } else {
            0.0
        };

        // Calculate size ratio
        if self.statistics.total_unique_reference > 0 {
            self.statistics.size_ratio = self.statistics.total_unique_sample as f64 / 
                self.statistics.total_unique_reference as f64;
        }

        // Calculate frequency distribution
        let dist = &mut self.statistics.frequency_distribution;
        for kmer in &self.shared_kmers {
            match kmer.reference_frequency {
                f if f >= 0.75 => dist.high_freq += 1,
                f if f >= 0.25 => dist.mid_freq += 1,
                _ => dist.low_freq += 1,
            }
        }

        // Calculate average frequency
        self.statistics.avg_frequency = if !self.shared_kmers.is_empty() {
            self.shared_kmers.iter()
                .map(|sk| sk.reference_frequency)
                .sum::<f64>() / self.shared_kmers.len() as f64
        } else {
            0.0
        };
    }
}