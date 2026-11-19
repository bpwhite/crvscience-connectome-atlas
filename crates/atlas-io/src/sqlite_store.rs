// Copyright (C) 2026 CRV Science Network (https://crvscience.com)
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://gnu.org>.

use atlas_core::{AtlasEnvelope, ConnectomeAtlas};
use rusqlite::{params, Connection};
use serde::{Deserialize, Serialize};
use std::path::Path;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum VaultError {
    #[error("SQLite storage error: {0}")]
    Sqlite(#[from] rusqlite::Error),
    #[error("JSON serialization error: {0}")]
    Json(#[from] serde_json::Error),
    #[error("Snapshot envelope failed SHA-256 integrity check")]
    IntegrityMismatch,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SnapshotMeta {
    pub id: i64,
    pub slot_name: String,
    pub atlas_id: String,
    pub organism: String,
    pub sha256_digest: String,
    pub ablated_count: usize,
    pub severed_count: usize,
    pub saved_at_utc: String,
}

pub struct AtlasVault {
    conn: Connection,
}

impl AtlasVault {
    pub fn open(path: impl AsRef<Path>) -> Result<Self, VaultError> {
        let conn = Connection::open(path)?;
        Self::init(conn)
    }

    pub fn open_in_memory() -> Result<Self, VaultError> {
        let conn = Connection::open_in_memory()?;
        Self::init(conn)
    }

    fn init(conn: Connection) -> Result<Self, VaultError> {
        conn.execute_batch(
            "PRAGMA journal_mode = WAL;
             PRAGMA synchronous = NORMAL;
             PRAGMA foreign_keys = ON;
             CREATE TABLE IF NOT EXISTS atlas_snapshots (
                 id INTEGER PRIMARY KEY AUTOINCREMENT,
                 slot_name TEXT NOT NULL UNIQUE,
                 atlas_id TEXT NOT NULL,
                 organism TEXT NOT NULL,
                 sha256_digest TEXT NOT NULL,
                 ablated_count INTEGER NOT NULL,
                 severed_count INTEGER NOT NULL,
                 envelope_json TEXT NOT NULL,
                 saved_at_utc TEXT NOT NULL DEFAULT (strftime('%Y-%m-%dT%H:%M:%SZ', 'now'))
             );
             CREATE TABLE IF NOT EXISTS audit_log (
                 id INTEGER PRIMARY KEY AUTOINCREMENT,
                 event_kind TEXT NOT NULL,
                 detail TEXT NOT NULL,
                 logged_at_utc TEXT NOT NULL DEFAULT (strftime('%Y-%m-%dT%H:%M:%SZ', 'now'))
             );",
        )?;
        Ok(Self { conn })
    }

    pub fn save_snapshot(&self, slot_name: &str, atlas: &ConnectomeAtlas) -> Result<String, VaultError> {
        let envelope = AtlasEnvelope::wrap(atlas.clone())?;
        let json = serde_json::to_string_pretty(&envelope)?;
        self.conn.execute(
            "INSERT INTO atlas_snapshots (
                 slot_name, atlas_id, organism, sha256_digest, ablated_count, severed_count, envelope_json
             ) VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7)
             ON CONFLICT(slot_name) DO UPDATE SET
                 atlas_id = excluded.atlas_id,
                 organism = excluded.organism,
                 sha256_digest = excluded.sha256_digest,
                 ablated_count = excluded.ablated_count,
                 severed_count = excluded.severed_count,
                 envelope_json = excluded.envelope_json,
                 saved_at_utc = strftime('%Y-%m-%dT%H:%M:%SZ', 'now')",
            params![
                slot_name,
                atlas.id,
                atlas.organism,
                envelope.sha256_digest,
                atlas.ablated_regions.len() as i64,
                atlas.severed_tracts.len() as i64,
                json,
            ],
        )?;
        self.conn.execute(
            "INSERT INTO audit_log (event_kind, detail) VALUES ('SAVE_SNAPSHOT', ?1)",
            params![format!("slot={slot_name} sha256={}", envelope.sha256_digest)],
        )?;
        Ok(envelope.sha256_digest)
    }

    pub fn load_snapshot(&self, slot_name: &str) -> Result<ConnectomeAtlas, VaultError> {
        let json: String = self.conn.query_row(
            "SELECT envelope_json FROM atlas_snapshots WHERE slot_name = ?1",
            params![slot_name],
            |row| row.get(0),
        )?;
        let env: AtlasEnvelope = serde_json::from_str(&json)?;
        if !env.verify_digest() {
            return Err(VaultError::IntegrityMismatch);
        }
        Ok(env.atlas)
    }

    pub fn list_snapshots(&self) -> Result<Vec<SnapshotMeta>, VaultError> {
        let mut stmt = self.conn.prepare(
            "SELECT id, slot_name, atlas_id, organism, sha256_digest, ablated_count, severed_count, saved_at_utc
             FROM atlas_snapshots ORDER BY id DESC",
        )?;
        let rows = stmt.query_map([], |r| {
            Ok(SnapshotMeta {
                id: r.get(0)?,
                slot_name: r.get(1)?,
                atlas_id: r.get(2)?,
                organism: r.get(3)?,
                sha256_digest: r.get(4)?,
                ablated_count: r.get::<_, i64>(5)? as usize,
                severed_count: r.get::<_, i64>(6)? as usize,
                saved_at_utc: r.get(7)?,
            })
        })?;
        let mut out = Vec::new();
        for item in rows {
            out.push(item?);
        }
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use atlas_core::presets::human_mni152_atlas;

    #[test]
    fn roundtrips_sqlite_wal_snapshot_with_sha256_verification() {
        let vault = AtlasVault::open_in_memory().expect("in-memory sqlite");
        let mut atlas = human_mni152_atlas();
        atlas.toggle_tract_sever("Corpus Callosum");
        let digest = vault.save_snapshot("callosotomy_study", &atlas).expect("save");
        assert_eq!(digest.len(), 64);

        let loaded = vault.load_snapshot("callosotomy_study").expect("load");
        assert!(loaded.severed_tracts.contains("Corpus Callosum"));
        assert_eq!(vault.list_snapshots().expect("list").len(), 1);
    }
}
