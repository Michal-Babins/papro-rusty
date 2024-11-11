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
                    info!("Found match: {} (similarity={:.4}, shared={}, unique={})",
                        name, 
                        profile_match.similarity_score,
                        profile_match.shared_kmers,
                        profile_match.unique_matches
                    );
                    matches.push(profile_match);
                }
                None => {
                    info!("Profile {} did not meet thresholds (min_similarity={}, min_shared_kmers={})",
                        name, self.min_similarity, self.min_shared_kmers);
                }
            }
        }
    
        matches.sort_by(|a, b| b.similarity_score.partial_cmp(&a.similarity_score).unwrap());
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
        let mut unique_matches = 0;
        let mut total_kmers_checked = 0;
    
        // Count profile kmers for accurate union calculation
        let mut profile_unique_kmers = HashSet::new();
    
        for kmer_result in kmer_stmt.query_map(params![profile_id], |row| {
            Ok((
                row.get::<_, String>(0)?,
                row.get::<_, f64>(1)?,
            ))
        })? {
            let (kmer, ref_freq) = kmer_result?;
            total_kmers_checked += 1;
            profile_unique_kmers.insert(kmer.clone());
            
            if let Some(&sample_count) = sample_kmers.get(&kmer) {
                shared_kmers += 1;
                if ref_freq >= 0.001 {
                    unique_matches += 1;
                }
            }
        }
    
        let sample_size = sample_kmers.len();
        let profile_size = profile_unique_kmers.len();
        let union_size = sample_size + profile_size - shared_kmers;
    
        // Calculate various similarity metrics
        let jaccard_similarity = shared_kmers as f64 / union_size as f64;
        let sample_coverage = shared_kmers as f64 / sample_size as f64;
        let profile_coverage = shared_kmers as f64 / profile_size as f64;
    
        info!(
            "Comparison summary for {}:
            Total k-mers checked: {}
            Shared k-mers: {}
            Sample unique k-mers: {}
            Profile unique k-mers: {}
            Union size: {}
            
            Similarity metrics:
            Jaccard similarity: {:.6}
            Sample coverage: {:.6} (shared/sample)
            Profile coverage: {:.6} (shared/profile)
            Unique marker matches: {}",
            profile_name, 
            total_kmers_checked,
            shared_kmers,
            sample_size,
            profile_size,
            union_size,
            jaccard_similarity,
            sample_coverage,
            profile_coverage,
            unique_matches
        );
    
        
        if sample_coverage >= self.min_similarity && 
        profile_coverage >= 0.05 && // At least 5% of the profile should match
        shared_kmers >= self.min_shared_kmers 
        {
            Ok(Some(ProfileMatch::new(
                profile_name.to_string(),
                sample_coverage,  // Keep using sample coverage as main score
                shared_kmers,
                unique_matches,
                sample_size,
            )))
        } else {
            info!(
                "Profile {} did not meet thresholds:
                Sample coverage: {:.6} (minimum: {})
                Profile coverage: {:.6} (minimum: 0.05)
                Shared k-mers: {} (minimum: {})",
                profile_name, 
                sample_coverage, 
                self.min_similarity,
                profile_coverage,
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
        info!("Adding shared k-mer: {:?}", sample_kmers);
        let mut analysis = DetailedAnalysis::new();

        let profile_kmers = kmer_stmt.query_map(params![profile_id], |row| {
            Ok((
                row.get::<_, String>(0)?,
                row.get::<_, f64>(1)?,
            ))
        })?;

        // Here we compare the kmers
        let total_sample_kmers = counter.total_kmers() as f64;
        for kmer_result in profile_kmers {
            let (kmer, ref_freq) = kmer_result?;
            
            if let Some(&sample_count) = sample_kmers.get(&kmer) {
                let sample_freq = sample_count as f64 / total_sample_kmers;
                analysis.add_shared_kmer(kmer, sample_freq, ref_freq);
            } else {
                analysis.add_reference_unique_kmer(kmer, ref_freq);
            }
        }

        // For a sample get unique kmers
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


#[derive(Debug)]
pub struct DetailedAnalysis {
    pub shared_kmers: Vec<SharedKmer>,
    pub unique_to_reference: Vec<UniqueKmer>,
    pub unique_to_sample: Vec<UniqueKmer>,
    pub statistics: AnalysisStatistics,
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
pub struct AnalysisStatistics {
    pub total_shared: usize,
    pub total_unique_reference: usize,
    pub total_unique_sample: usize,
    pub average_frequency_difference: f64,
    pub marker_kmer_matches: usize,
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

        let total_diff: f64 = self.shared_kmers.iter()
            .map(|sk| (sk.sample_frequency - sk.reference_frequency).abs())
            .sum();

        self.statistics.average_frequency_difference = if self.shared_kmers.is_empty() {
            0.0
        } else {
            total_diff / self.shared_kmers.len() as f64
        };

        self.statistics.marker_kmer_matches = self.shared_kmers.iter()
            .filter(|sk| sk.reference_frequency >= 0.001)
            .count();
    }

}
