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