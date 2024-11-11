use rayon::iter::IntoParallelIterator;
use rusqlite::{params, Connection, OptionalExtension};
use anyhow::Result;
use std::path::{Path, PathBuf};
use log::{info, warn};

use super::schemas::initialize_schema;
use super::types::{DatabaseStats, ProfileSummary};
use crate::io::FastxReader;
use crate::kmer::KmerCounter;
use crate::profile::{Profile, TaxonomyLevel};

pub struct Database {
    conn: Connection,
}

impl Database {
    /// Create a new database or open existing one
    pub fn new<P: AsRef<Path>>(path: P) -> Result<Self> {
        let conn = Connection::open(path)?;
        initialize_schema(&conn)?;
        Ok(Database { conn })
    }

    /// Create a profile from multiple FASTA/FASTQ files
    pub fn create_profile(
        &mut self,
        input_files: Vec<PathBuf>,
        kmer_size: usize,
        level: TaxonomyLevel,
        name: String,
    ) -> Result<Profile> {
        // Initialize k-mer counter
        let counter = KmerCounter::new(kmer_size);
        
        // Process all input files
        info!("Processing {} input files...", input_files.len());
        for (idx, file) in input_files.iter().enumerate() {
            info!("Processing file {}/{}: {}", 
                idx + 1, 
                input_files.len(), 
                file.display()
            );
            
            let reader = FastxReader::new(vec![file.clone()]);
            let mut sequences = Vec::new();
            reader.process_all(|sequence, _id| {
                sequences.push(sequence.to_vec());
                Ok(())
            })?;

            counter.count_sequences(sequences.into_par_iter())?;
        }

        info!("Found {} unique k-mers across all files", counter.unique_kmers());

        // Check if profile already exists
        let exists: bool = self.conn.query_row(
            "SELECT 1 FROM profiles WHERE name = ?",
            params![&name],
            |_| Ok(true)
        ).unwrap_or(false);

        if exists {
            return Err(anyhow::anyhow!("Profile {} already exists in database", name));
        }

        // Create profile
        let mut profile = Profile::new(name, level, kmer_size);

        // Calculate frequencies from total counts
        let total_kmers = counter.total_kmers() as f64;
        for (kmer, count) in counter.get_counts() {
            let frequency = count as f64 / total_kmers;
            profile.frequencies.insert(kmer, frequency);
        }
        profile.total_kmers = counter.total_kmers();

        info!(
            "Created profile with {} k-mers from {} files", 
            profile.frequencies.len(),
            input_files.len()
        );

        // Add profile to database
        self.add_profile(&profile)?;
        
        Ok(profile)
    }

    /// Add a new profile to the database
    pub fn add_profile(&mut self, profile: &Profile) -> Result<()> {
        // Check if profile already exists
        let exists: bool = self.conn.query_row(
            "SELECT 1 FROM profiles WHERE name = ?",
            params![profile.name],
            |_| Ok(true)
        ).unwrap_or(false);

        if exists {
            warn!("Profile {} already exists in database", profile.name);
            return Ok(());
        }

        let tx = self.conn.transaction()?;
        
        // Insert profile
        tx.execute(
            "INSERT INTO profiles (name, taxonomy_level, k, total_kmers)
             VALUES (?1, ?2, ?3, ?4)",
            params![
                profile.name,
                profile.level.to_string(),
                profile.k,
                profile.total_kmers,
            ],
        )?;

        let profile_id = tx.last_insert_rowid();

        // Insert k-mers
        {
            let mut stmt = tx.prepare(
                "INSERT INTO kmers (profile_id, kmer, frequency)
                 VALUES (?1, ?2, ?3)"
            )?;

            for (kmer, frequency) in &profile.frequencies {
                stmt.execute(params![profile_id, kmer, frequency])?;
            }
        }

        tx.commit()?;
        info!("Added profile {} to database", profile.name);
        Ok(())
    }

    /// Remove a profile from the database
    pub fn remove_profile(&mut self, name: &str) -> Result<bool> {
        let tx = self.conn.transaction()?;
        
        let profile_id: Option<i64> = tx.query_row(
            "SELECT id FROM profiles WHERE name = ?",
            params![name],
            |row| row.get(0)
        ).optional()?;

        if let Some(id) = profile_id {
            // Delete k-mers first (foreign key constraint)
            tx.execute(
                "DELETE FROM kmers WHERE profile_id = ?",
                params![id]
            )?;
            
            // Delete profile
            tx.execute(
                "DELETE FROM profiles WHERE id = ?",
                params![id]
            )?;
            
            tx.commit()?;
            info!("Removed profile {} from database", name);
            Ok(true)
        } else {
            Ok(false)
        }
    }

    /// Get a profile by name
    pub fn get_profile(&self, name: &str) -> Result<Option<Profile>> {
        let profile_result = self.conn.query_row(
            "SELECT taxonomy_level, k, total_kmers 
             FROM profiles WHERE name = ?",
            params![name],
            |row| {
                let level_str: String = row.get(0)?;
                let level = match level_str.as_str() {
                    "Genus" => TaxonomyLevel::Genus,
                    "Species" => TaxonomyLevel::Species,
                    "Strain" => TaxonomyLevel::Strain,
                    _ => return Err(rusqlite::Error::InvalidParameterName(level_str)),
                };

                Ok((level, row.get::<_, i64>(1)?, row.get::<_, i64>(2)?))
            }
        ).optional()?;

        if let Some((level, k, total_kmers)) = profile_result {
            let mut profile = Profile::new(
                name.to_string(),
                level,
                k as usize,
            );
            profile.total_kmers = total_kmers as usize;

            // Get k-mers
            let mut stmt = self.conn.prepare(
                "SELECT kmer, frequency 
                 FROM kmers 
                 WHERE profile_id = (SELECT id FROM profiles WHERE name = ?)"
            )?;

            let kmers = stmt.query_map(params![name], |row| {
                Ok((row.get::<_, String>(0)?, row.get::<_, f64>(1)?))
            })?;

            for kmer in kmers {
                let (kmer, freq) = kmer?;
                profile.frequencies.insert(kmer, freq);
            }

            Ok(Some(profile))
        } else {
            Ok(None)
        }
    }

    /// List all profiles, optionally filtered by taxonomy level
    pub fn list_profiles(&self, level: Option<TaxonomyLevel>) -> Result<Vec<ProfileSummary>> {
        let query = match level {
            Some(_) => 
                "SELECT name, taxonomy_level, k, total_kmers, created_at 
                 FROM profiles 
                 WHERE taxonomy_level = ?
                 ORDER BY name",
            None => 
                "SELECT name, taxonomy_level, k, total_kmers, created_at 
                 FROM profiles 
                 ORDER BY name",
        };

        let mut stmt = self.conn.prepare(query)?;
        let mut rows = match level {
            Some(l) => stmt.query(params![l.to_string()])?,
            None => stmt.query([])?,
        };

        let mut profiles = Vec::new();
        while let Some(row) = rows.next()? {
            let level = match row.get::<_, String>(1)?.as_str() {
                "Genus" => TaxonomyLevel::Genus,
                "Species" => TaxonomyLevel::Species,
                "Strain" => TaxonomyLevel::Strain,
                l => return Err(anyhow::anyhow!("Invalid taxonomy level in database: {}", l)),
            };

            profiles.push(ProfileSummary {
                name: row.get(0)?,
                level,
                k: row.get::<_, i64>(2)? as usize,
                total_kmers: row.get::<_, i64>(3)? as usize,
                created_at: row.get(4)?,
            });
        }

        Ok(profiles)
    }

    /// Get database statistics
    pub fn get_statistics(&self) -> Result<DatabaseStats> {
        let total_profiles: i64 = self.conn.query_row(
            "SELECT COUNT(*) FROM profiles",
            [],
            |row| row.get(0)
        )?;

        let total_kmers: i64 = self.conn.query_row(
            "SELECT COUNT(*) FROM kmers",
            [],
            |row| row.get(0)
        )?;

        let level_counts: Vec<(String, i64)> = self.conn.prepare(
            "SELECT taxonomy_level, COUNT(*) 
             FROM profiles 
             GROUP BY taxonomy_level"
        )?.query_map([], |row| {
            Ok((row.get(0)?, row.get(1)?))
        })?.collect::<rusqlite::Result<_>>()?;

        Ok(DatabaseStats {
            total_profiles: total_profiles as usize,
            total_kmers: total_kmers as usize,
            profiles_by_level: level_counts.into_iter()
                .map(|(level, count)| (level, count as usize))
                .collect(),
        })
    }

    pub fn validate(&self) -> Result<ValidationReport> {
        let mut report = ValidationReport::default();

        // 1. Check table existence and schema
        self.validate_schema(&mut report)?;

        // 2. Check data integrity
        self.validate_data(&mut report)?;

        // 3. Check referential integrity
        self.validate_references(&mut report)?;

        Ok(report)
    }

    fn validate_schema(&self, report: &mut ValidationReport) -> Result<()> {
        // Check if tables exist
        let tables = self.conn.query_row(
            "SELECT COUNT(*) FROM sqlite_master WHERE type='table' AND name IN ('profiles', 'kmers')",
            [],
            |row| row.get::<_, i64>(0)
        )?;

        if tables != 2 {
            report.add_error("Missing required tables (profiles and/or kmers)");
            return Ok(());
        }

        // Verify profiles table schema
        let profile_cols = self.conn.query_row(
            "SELECT COUNT(*) FROM pragma_table_info('profiles') WHERE 
             name IN ('id', 'name', 'taxonomy_level', 'k', 'total_kmers', 'created_at')",
            [],
            |row| row.get::<_, i64>(0)
        )?;

        if profile_cols != 6 {
            report.add_error("Profiles table is missing required columns");
        }

        // Verify kmers table schema
        let kmer_cols = self.conn.query_row(
            "SELECT COUNT(*) FROM pragma_table_info('kmers') WHERE 
             name IN ('profile_id', 'kmer', 'frequency')",
            [],
            |row| row.get::<_, i64>(0)
        )?;

        if kmer_cols != 3 {
            report.add_error("Kmers table is missing required columns");
        }

        Ok(())
    }

    fn validate_data(&self, report: &mut ValidationReport) -> Result<()> {
        // Check taxonomy levels are valid
        let invalid_levels: Vec<String> = self.conn.prepare(
            "SELECT DISTINCT taxonomy_level FROM profiles 
             WHERE taxonomy_level NOT IN ('Species', 'Genus', 'Strain')"
        )?.query_map([], |row| row.get(0))?
        .collect::<rusqlite::Result<_>>()?;

        if !invalid_levels.is_empty() {
            report.add_error(format!(
                "Invalid taxonomy levels found: {}", 
                invalid_levels.join(", ")
            ));
        }

        // Check for negative k-mer sizes or total counts
        let invalid_counts = self.conn.query_row(
            "SELECT COUNT(*) FROM profiles WHERE k <= 0 OR total_kmers <= 0",
            [],
            |row| row.get::<_, i64>(0)
        )?;

        if invalid_counts > 0 {
            report.add_error("Found profiles with invalid k-mer size or total count");
        }

        // Check k-mer frequencies are valid (between 0 and 1)
        let invalid_freqs = self.conn.query_row(
            "SELECT COUNT(*) FROM kmers WHERE frequency <= 0 OR frequency > 1",
            [],
            |row| row.get::<_, i64>(0)
        )?;

        if invalid_freqs > 0 {
            report.add_error("Found k-mers with invalid frequencies");
        }

        // Check frequency sums per profile approximately equal 1
        let mut stmt = self.conn.prepare(
            "SELECT profile_id, SUM(frequency) FROM kmers GROUP BY profile_id"
        )?;

        let sums = stmt.query_map([], |row| {
            Ok((row.get::<_, i64>(0)?, row.get::<_, f64>(1)?))
        })?;

        for sum in sums {
            let (profile_id, freq_sum) = sum?;
            if (freq_sum - 1.0).abs() > 0.01 {
                report.add_error(format!(
                    "Profile {} has total frequency sum of {:.6} (expected ≈1.0)",
                    profile_id, freq_sum
                ));
            }
        }

        Ok(())
    }

    fn validate_references(&self, report: &mut ValidationReport) -> Result<()> {
        // Check for orphaned k-mers (no matching profile)
        let orphaned = self.conn.query_row(
            "SELECT COUNT(*) FROM kmers k 
             LEFT JOIN profiles p ON k.profile_id = p.id 
             WHERE p.id IS NULL",
            [],
            |row| row.get::<_, i64>(0)
        )?;

        if orphaned > 0 {
            report.add_error(format!("Found {} orphaned k-mer entries", orphaned));
        }

        // Check each profile has k-mers
        let empty_profiles = self.conn.prepare(
            "SELECT name FROM profiles p 
             LEFT JOIN kmers k ON p.id = k.profile_id 
             GROUP BY p.id HAVING COUNT(k.kmer) = 0"
        )?.query_map([], |row| row.get::<_, String>(0))?
        .collect::<rusqlite::Result<Vec<_>>>()?;

        if !empty_profiles.is_empty() {
            report.add_warning(format!(
                "Found profiles with no k-mers: {}", 
                empty_profiles.join(", ")
            ));
        }

        Ok(())
    }
}

#[derive(Default, Debug)]
pub struct ValidationReport {
    errors: Vec<String>,
    warnings: Vec<String>,
}

impl ValidationReport {
    fn add_error<S: Into<String>>(&mut self, msg: S) {
        self.errors.push(msg.into());
    }

    fn add_warning<S: Into<String>>(&mut self, msg: S) {
        self.warnings.push(msg.into());
    }

    pub fn has_errors(&self) -> bool {
        !self.errors.is_empty()
    }

    pub fn has_warnings(&self) -> bool {
        !self.warnings.is_empty()
    }

    pub fn errors(&self) -> &[String] {
        &self.errors
    }

    pub fn warnings(&self) -> &[String] {
        &self.warnings
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    #[test]
    fn test_database_creation() -> Result<()> {
        let dir = tempdir()?;
        let db_path = dir.path().join("test.db");
        
        let mut db = Database::new(&db_path)?;
        
        let mut profile = Profile::new(
            "Test_Species".to_string(),
            TaxonomyLevel::Species,
            4,
        );
        profile.frequencies.insert("AAAA".to_string(), 0.5);
        profile.frequencies.insert("TTTT".to_string(), 0.5);
        profile.total_kmers = 2;

        db.add_profile(&profile)?;

        let retrieved = db.get_profile("Test_Species")?.unwrap();
        assert_eq!(retrieved.name, profile.name);
        assert_eq!(retrieved.level, profile.level);
        assert_eq!(retrieved.frequencies.len(), profile.frequencies.len());

        Ok(())
    }
}