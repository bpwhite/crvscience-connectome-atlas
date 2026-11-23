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

use glam::Vec3;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::collections::BTreeSet;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum AtlasScale {
    MacroParcellation,
    MesoscaleProjection,
    MicroNeuropilSynaptic,
    CellularCompleteConnectome,
}

impl AtlasScale {
    pub fn label(self) -> &'static str {
        match self {
            Self::MacroParcellation => "Macro-Connectome (Stereotaxic MNI mm)",
            Self::MesoscaleProjection => "Mesoscale Atlas (Allen CCFv3 10 µm ref)",
            Self::MicroNeuropilSynaptic => "Whole-Brain Neuropil Connectome (µm)",
            Self::CellularCompleteConnectome => "Complete Cellular Connectome (Soma/Synapse)",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Hemisphere {
    Left,
    Right,
    Midline,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum PathwayModality {
    Association,
    Commissural,
    Projection,
    ChemicalSynapse,
    ElectricalGapJunction,
}

impl PathwayModality {
    pub fn label(self) -> &'static str {
        match self {
            Self::Association => "Association Tract",
            Self::Commissural => "Commissural Tract",
            Self::Projection => "Projection / Subcortical Tract",
            Self::ChemicalSynapse => "Chemical Synapse (Directed)",
            Self::ElectricalGapJunction => "Electrical Gap Junction",
        }
    }
}

/// Directionally encoded color (DEC) following diffusion MRI conventions:
/// Red = Left-Right |dx|, Green = Anterior-Posterior |dy|, Blue = Superior-Inferior |dz|.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct DecColor {
    pub r: u8,
    pub g: u8,
    pub b: u8,
}

impl DecColor {
    pub fn from_endpoints(a: [f32; 3], b: [f32; 3]) -> Self {
        let dir = Vec3::from_array(b) - Vec3::from_array(a);
        let v = dir.abs().normalize_or_zero();
        let boost = |c: f32| -> u8 { ((c.powf(0.75) * 235.0 + 20.0).clamp(20.0, 255.0)) as u8 };
        Self {
            r: boost(v.x),
            g: boost(v.y),
            b: boost(v.z),
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct ReceptorProfile {
    pub da: f32,
    pub serotonin: f32,
    pub ach: f32,
    pub glutamate: f32,
    pub gaba: f32,
}

impl ReceptorProfile {
    pub fn channel_value(&self, channel: &str) -> f64 {
        match channel {
            "da" | "dopamine" => self.da as f64,
            "ht" | "5ht" | "serotonin" => self.serotonin as f64,
            "ach" | "acetylcholine" => self.ach as f64,
            "glu" | "glutamate" => self.glutamate as f64,
            "gaba" => self.gaba as f64,
            _ => self.serotonin as f64,
        }
    }
}

impl Default for ReceptorProfile {
    fn default() -> Self {
        Self {
            da: 0.45,
            serotonin: 0.50,
            ach: 0.50,
            glutamate: 0.65,
            gaba: 0.60,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BrainRegion {
    pub id: String,
    pub name: String,
    pub network: String,
    pub hemisphere: Hemisphere,
    /// Stereotaxic [x, y, z] coordinates (e.g. MNI152 mm: +X Right, +Y Anterior, +Z Superior).
    pub coord: [f32; 3],
    pub radius_mm: f32,
    pub color_hex: String,
    pub citation_url: String,
    pub notes: String,
    #[serde(default)]
    pub what_it_does: String,
    #[serde(default)]
    pub inputs_prose: String,
    #[serde(default)]
    pub outputs_prose: String,
    #[serde(default)]
    pub receptors: ReceptorProfile,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PathwayEdge {
    pub source: String,
    pub target: String,
    pub weight: f64,
    pub tract_name: String,
    pub modality: PathwayModality,
    pub directed: bool,
    pub synapse_count_est: Option<u64>,
}

fn default_conduction_velocity_m_s() -> f64 {
    9.6
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConnectomeAtlas {
    pub id: String,
    pub title: String,
    pub organism: String,
    pub coordinate_space: String,
    pub scale: AtlasScale,
    pub reference_summary: String,
    pub primary_citation_url: String,
    pub benchmark_neurons_or_parcels: String,
    pub benchmark_synapses_or_tracts: String,
    #[serde(default = "default_conduction_velocity_m_s")]
    pub conduction_velocity_m_s: f64,
    pub regions: Vec<BrainRegion>,
    pub edges: Vec<PathwayEdge>,
    pub ablated_regions: BTreeSet<String>,
    pub severed_tracts: BTreeSet<String>,
}

impl ConnectomeAtlas {
    pub fn is_edge_intact(&self, edge: &PathwayEdge, min_weight: f64) -> bool {
        edge.weight >= min_weight
            && !self.ablated_regions.contains(&edge.source)
            && !self.ablated_regions.contains(&edge.target)
            && !self.severed_tracts.contains(&edge.tract_name)
    }

    pub fn unique_tracts(&self) -> Vec<String> {
        let mut set = BTreeSet::new();
        for e in &self.edges {
            set.insert(e.tract_name.clone());
        }
        set.into_iter().collect()
    }

    pub fn unique_networks(&self) -> Vec<String> {
        let mut set = BTreeSet::new();
        for r in &self.regions {
            set.insert(r.network.clone());
        }
        set.into_iter().collect()
    }

    pub fn toggle_region_ablation(&mut self, region_id: &str) {
        if !self.ablated_regions.remove(region_id) {
            self.ablated_regions.insert(region_id.to_string());
        }
    }

    pub fn toggle_tract_sever(&mut self, tract_name: &str) {
        if !self.severed_tracts.remove(tract_name) {
            self.severed_tracts.insert(tract_name.to_string());
        }
    }

    pub fn heal_all(&mut self) {
        self.ablated_regions.clear();
        self.severed_tracts.clear();
    }

    /// Returns `(incoming_afferents, outgoing_efferents, coupled_bidirectional)` for `region_id`.
    pub fn incident_io_breakdown(
        &self,
        region_id: &str,
        min_weight: f64,
    ) -> (Vec<&PathwayEdge>, Vec<&PathwayEdge>, Vec<&PathwayEdge>) {
        let mut incoming = Vec::new();
        let mut outgoing = Vec::new();
        let mut coupled = Vec::new();
        for edge in &self.edges {
            if !self.is_edge_intact(edge, min_weight) {
                continue;
            }
            if edge.source != region_id && edge.target != region_id {
                continue;
            }
            if edge.directed {
                if edge.target == region_id {
                    incoming.push(edge);
                } else {
                    outgoing.push(edge);
                }
            } else {
                coupled.push(edge);
            }
        }
        (incoming, outgoing, coupled)
    }
}

/// Synchronized `.connectome.json` envelope shared between the Rust workspace and 3D Web Preview.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AtlasEnvelope {
    pub schema_version: String,
    pub app_version: String,
    pub sha256_digest: String,
    pub atlas: ConnectomeAtlas,
}

impl AtlasEnvelope {
    pub fn wrap(atlas: ConnectomeAtlas) -> Result<Self, serde_json::Error> {
        let payload = serde_json::to_vec(&atlas)?;
        let mut hasher = Sha256::new();
        hasher.update(&payload);
        let digest = format!("{:x}", hasher.finalize());
        Ok(Self {
            schema_version: "connectome-atlas-v1".to_string(),
            app_version: "0.8.0".to_string(),
            sha256_digest: digest,
            atlas,
        })
    }

    pub fn verify_digest(&self) -> bool {
        if let Ok(payload) = serde_json::to_vec(&self.atlas) {
            let mut hasher = Sha256::new();
            hasher.update(&payload);
            let actual = format!("{:x}", hasher.finalize());
            actual == self.sha256_digest
        } else {
            false
        }
    }
}
