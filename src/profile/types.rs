use serde::{Serialize, Deserialize};
use std::collections::HashMap;

/// Represents the taxonomic level for a profile
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub enum TaxonomyLevel {
    Genus,
    Species,
    Strain,
}

impl std::fmt::Display for TaxonomyLevel {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            TaxonomyLevel::Genus => write!(f, "Genus"),
            TaxonomyLevel::Species => write!(f, "Species"),
            TaxonomyLevel::Strain => write!(f, "Strain"),
        }
    }
}

/// Represents a profile match with its similarity metrics
#[derive(Debug, Clone)]
pub struct ProfileMatch {
    /// Name of the matched profile
    pub name: String,
    /// Sample coverage score (0.0 to 1.0)
    pub sample_coverage: f64,
    /// Number of k-mers shared between sample and reference
    pub shared_kmers: usize,
    /// Number of core k-mers (frequency >= 0.95)
    pub core_matches: usize,
    /// Number of rare k-mers (frequency <= 0.10)
    pub rare_matches: usize,
    /// Average frequency of matched k-mers
    pub avg_frequency: f64,
    /// Size ratio of sample to profile
    pub size_ratio: f64,
}

/// Represents a k-mer profile
#[derive(Debug, Clone)]
pub struct Profile {
    /// Profile identifier
    pub name: String,
    /// Taxonomic level of this profile
    pub level: TaxonomyLevel,
    /// K-mer size used
    pub k: usize,
    /// K-mer frequencies (k-mer sequence -> frequency)
    pub frequencies: HashMap<String, f64>,
    /// Total number of k-mers
    pub total_kmers: usize,
}

impl Profile {
    /// Create a new profile
    pub fn new(name: String, level: TaxonomyLevel, k: usize) -> Self {
        Profile {
            name,
            level,
            k,
            frequencies: HashMap::new(),
            total_kmers: 0,
        }
    }

}

impl ProfileMatch {
    pub fn new(
        name: String,
        sample_coverage: f64,
        shared_kmers: usize,
        core_matches: usize,
        rare_matches: usize,
        avg_frequency: f64,
        size_ratio: f64,
    ) -> Self {
        ProfileMatch {
            name,
            sample_coverage,
            shared_kmers,
            core_matches,
            rare_matches,
            avg_frequency,
            size_ratio,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_profile_creation() {
        let profile = Profile::new(
            "Test Profile".to_string(),
            TaxonomyLevel::Species,
            31,
        );
        assert_eq!(profile.name, "Test Profile");
        assert_eq!(profile.level, TaxonomyLevel::Species);
        assert_eq!(profile.k, 31);
        assert_eq!(profile.total_kmers, 0);
        assert!(profile.frequencies.is_empty());
    }

    #[test]
    fn test_profile_match_creation() {
        let match_result = ProfileMatch::new(
            "Test Match".to_string(),
            0.95,  // sample_coverage
            1000,  // shared_kmers
            500,   // core_matches (k-mers with freq >= 0.95)
            100,   // rare_matches (k-mers with freq <= 0.10)
            0.45,  // avg_frequency
            0.85,  // size_ratio
        );

        assert_eq!(match_result.name, "Test Match");
        assert_eq!(match_result.sample_coverage, 0.95);
        assert_eq!(match_result.shared_kmers, 1000);
        assert_eq!(match_result.core_matches, 500);
        assert_eq!(match_result.rare_matches, 100);
        assert!((match_result.avg_frequency - 0.45).abs() < f64::EPSILON);
        assert!((match_result.size_ratio - 0.85).abs() < f64::EPSILON);
    }

    #[test]
    fn test_taxonomy_level_display() {
        assert_eq!(TaxonomyLevel::Genus.to_string(), "Genus");
        assert_eq!(TaxonomyLevel::Species.to_string(), "Species");
        assert_eq!(TaxonomyLevel::Strain.to_string(), "Strain");
    }
}