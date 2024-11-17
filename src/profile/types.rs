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
    /// Percentage of sample k-mers found in profile
    pub sample_coverage: f64,
    /// Number of k-mers shared between sample and reference
    pub shared_kmers: usize,
    /// Ratio of sample size to profile size
    pub size_ratio: f64,
    /// Score for how unique these matches are to this profile
    pub uniqueness_score: f64,
    /// Confidence score for this match
    pub confidence_score: f64,
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
        size_ratio: f64,
        uniqueness_score: f64,
        confidence_score: f64,
    ) -> Self {
        ProfileMatch {
            name,
            sample_coverage,
            shared_kmers,
            size_ratio,
            uniqueness_score,
            confidence_score,
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
        0.85,  // size_ratio
        0.75,  // uniqueness_score
        0.82,  // confidence_score
    );

    assert_eq!(match_result.name, "Test Match");
    assert_eq!(match_result.sample_coverage, 0.95);
    assert_eq!(match_result.shared_kmers, 1000);
    assert!((match_result.size_ratio - 0.85).abs() < f64::EPSILON);
    assert!((match_result.uniqueness_score - 0.75).abs() < f64::EPSILON);
    assert!((match_result.confidence_score - 0.82).abs() < f64::EPSILON);
    }

    #[test]
    fn test_taxonomy_level_display() {
        assert_eq!(TaxonomyLevel::Genus.to_string(), "Genus");
        assert_eq!(TaxonomyLevel::Species.to_string(), "Species");
        assert_eq!(TaxonomyLevel::Strain.to_string(), "Strain");
    }
}