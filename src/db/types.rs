use crate::profile::TaxonomyLevel;

/// Summary of a profile for listing
#[derive(Debug)]
pub struct ProfileSummary {
    pub name: String,
    pub level: TaxonomyLevel,
    pub k: usize,
    pub total_kmers: usize,
    pub created_at: String,
}

/// Database statistics
#[derive(Debug)]
pub struct DatabaseStats {
    pub total_profiles: usize,
    pub total_kmers: usize,
    pub profiles_by_level: Vec<(String, usize)>,
}