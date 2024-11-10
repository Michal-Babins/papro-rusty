use std::collections::HashMap;
use anyhow::{Result, Context};
use log::{info, warn};
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

        // Get all profiles at the specified taxonomy level
        let mut profile_stmt = self.conn.prepare(
            "SELECT id, name, k, total_kmers 
             FROM profiles 
             WHERE taxonomy_level = ?"
        )?;

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

        // Parellel process
        let sample_kmers = counter.get_counts();
        let total_sample_kmers = counter.total_kmers();

        let mut matches = Vec::new();
        for profile_result in profiles {
            let (profile_id, name, k, total_kmers) = profile_result?;

            // K-mer size mismatch doesn't skip
            if k as usize != counter.kmer_size() {
                warn!("Skipping profile {} due to k-mer size mismatch ({} vs {})",
                    name, k, counter.kmer_size());
                continue;
            }

            if let Some(profile_match) = self.compare_with_profile(
                profile_id,
                &name,
                &sample_kmers,
                total_sample_kmers,
                total_kmers as usize,
            )? {
                matches.push(profile_match);
            }
        }

        // Sort records by similirity score
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
        let mut kmer_stmt = self.conn.prepare(
            "SELECT kmer, frequency FROM kmers WHERE profile_id = ?"
        )?;

        let mut shared_kmers = 0;
        let mut similarity_score = 0.0;
        let mut unique_matches = 0;

        let profile_kmers = kmer_stmt.query_map(params![profile_id], |row| {
            Ok((
                row.get::<_, String>(0)?,
                row.get::<_, f64>(1)?,
            ))
        })?;

        // Calculate similarity metrics
        for kmer_result in profile_kmers {
            let (kmer, ref_freq) = kmer_result?;
            
            if let Some(&sample_count) = sample_kmers.get(&kmer) {
                shared_kmers += 1;
                
                // Let's use frequencies 
                let sample_freq = sample_count as f64 / total_sample_kmers as f64;
                
                // Use a wieghted Jaccard Similarty approach
                similarity_score += (sample_freq * ref_freq).sqrt();

                // Consider as marker k-mer if reference frequency is significant
                if ref_freq >= 0.001 {
                    unique_matches += 1;
                }
            }
        }

        // Normalize similarity score
        if shared_kmers > 0 {
            similarity_score /= (total_sample_kmers as f64 * total_profile_kmers as f64).sqrt();
        }

        
        if similarity_score >= self.min_similarity && shared_kmers >= self.min_shared_kmers {
            Ok(Some(ProfileMatch::new(
                profile_name.to_string(),
                similarity_score,
                shared_kmers,
                unique_matches,
                total_sample_kmers,
            )))
        } else {
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
