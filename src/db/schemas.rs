use rusqlite::{Connection, Result};

pub(crate) fn initialize_schema(conn: &Connection) -> Result<()> {
    // Create profiles table
    conn.execute(
        "CREATE TABLE IF NOT EXISTS profiles (
            id INTEGER PRIMARY KEY,
            name TEXT NOT NULL UNIQUE,
            taxonomy_level TEXT NOT NULL,
            k INTEGER NOT NULL,
            total_kmers INTEGER NOT NULL,
            created_at DATETIME DEFAULT CURRENT_TIMESTAMP
        )",
        [],
    )?;

    // Create kmers table
    conn.execute(
        "CREATE TABLE IF NOT EXISTS kmers (
            profile_id INTEGER,
            kmer TEXT NOT NULL,
            frequency REAL NOT NULL,
            FOREIGN KEY(profile_id) REFERENCES profiles(id),
            PRIMARY KEY(profile_id, kmer)
        )",
        [],
    )?;

    // Create indices
    conn.execute(
        "CREATE INDEX IF NOT EXISTS idx_kmers_profile 
         ON kmers(profile_id)",
        [],
    )?;

    conn.execute(
        "CREATE INDEX IF NOT EXISTS idx_profiles_taxonomy 
         ON profiles(taxonomy_level)",
        [],
    )?;

    Ok(())
}