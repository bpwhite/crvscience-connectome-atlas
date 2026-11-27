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

use crate::model::ConnectomeAtlas;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum CommandAction {
    JumpToRegion(String),
    SwitchPreset(String),
    ToggleZenMode,
    SetControlScope(String),
    SetCutaway(String),
    ToggleFlyThrough,
    HealAllLesions,
    ExportEnvelopeJson,
    ExportSelectedSwc,
    TuneCriticalityKappa(i32),
    InspectPersistentCavities,
    InspectTopologyDigest,
    ExportTopologyDigestJson,
    ExportTopologyDigestCsv,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CommandEntry {
    pub category: String,
    pub title: String,
    pub subtitle: String,
    pub action: CommandAction,
    search_blob: String,
}

/// Builds the quick-jump `Ctrl/Cmd+K` command index across all regions in the active connectome
/// plus workspace layout, preset, cutaway, and export actions.
pub fn build_command_index(atlas: &ConnectomeAtlas) -> Vec<CommandEntry> {
    let mut entries = Vec::with_capacity(atlas.regions.len() + 16);

    for r in &atlas.regions {
        let title = format!("{} — {}", r.id, r.name);
        let subtitle = format!(
            "{} · ({:.1}, {:.1}, {:.1})",
            r.network, r.coord[0], r.coord[1], r.coord[2]
        );
        let search_blob = format!(
            "{} {} {} {}",
            r.id.to_lowercase(),
            r.name.to_lowercase(),
            r.network.to_lowercase(),
            r.what_it_does.to_lowercase()
        );
        entries.push(CommandEntry {
            category: "Region".to_string(),
            title,
            subtitle,
            action: CommandAction::JumpToRegion(r.id.clone()),
            search_blob,
        });
    }

    let builtins = [
        (
            "Criticality",
            "Criticality: Poised Homeostatic Regime (kappa = 1.00)",
            "Maximize Kinouchi-Copelli dynamic range Delta (dB) near sigma = 1.0",
            CommandAction::TuneCriticalityKappa(100),
        ),
        (
            "Topology",
            "Algebraic Topology: Inspect Persistent beta-1 Cavities & Cliques",
            "Enumerate 3/4-cliques and Z2 persistent homology loop cavities",
            CommandAction::InspectPersistentCavities,
        ),
        (
            "Layout",
            "Toggle Zen Full-Canvas 3D Mode (Z)",
            "Collapse sidebars and 2D slice dock for 100% 3D viewport",
            CommandAction::ToggleZenMode,
        ),
        (
            "Scope",
            "Scope: Explore (Camera & Palette)",
            "Progressive-disclosure mode showing primary exploration controls",
            CommandAction::SetControlScope("explore".to_string()),
        ),
        (
            "Scope",
            "Scope: 3D Shaders, Cutaway & Tractotomy",
            "Focus sidebar on surgical cutaways, fascicles, and tract bundles",
            CommandAction::SetControlScope("surgical".to_string()),
        ),
        (
            "Scope",
            "Scope: Dynamics, Receptors & Harmonics",
            "Focus sidebar on receptor gain alpha and Laplacian harmonics",
            CommandAction::SetControlScope("dynamics".to_string()),
        ),
        (
            "3D Cutaway",
            "Cutaway: Anterior-Superior Surgical Octant",
            "Carve front-left cortical octant window",
            CommandAction::SetCutaway("octant".to_string()),
        ),
        (
            "3D Cutaway",
            "Cutaway: Left Hemisphere Window",
            "Expose deep left subcortical nuclei and commissures",
            CommandAction::SetCutaway("left_hemi".to_string()),
        ),
        (
            "3D Cutaway",
            "Cutaway: Superior Axial Dome",
            "Slice dorsal cortical dome to view ventricles and thalamus",
            CommandAction::SetCutaway("axial_dome".to_string()),
        ),
        (
            "3D Cutaway",
            "Cutaway: Intact Outer Shell",
            "Restore complete bilateral outer shell",
            CommandAction::SetCutaway("intact".to_string()),
        ),
        (
            "Camera",
            "Fly Active Signaling Pathway",
            "Spline-guided 3D chase camera along traced route",
            CommandAction::ToggleFlyThrough,
        ),
        (
            "Preset",
            "Organism: Human MNI152 Macro-Connectome",
            "40 ROIs across 9 functional systems",
            CommandAction::SwitchPreset("human_mni152".to_string()),
        ),
        (
            "Preset",
            "Organism: Mouse Allen CCFv3 Mesoscale",
            "28 structures with directed viral-tracing pathways",
            CommandAction::SwitchPreset("mouse_ccfv3".to_string()),
        ),
        (
            "Preset",
            "Organism: Drosophila FlyWire FAFB v783",
            "22 neuropils with branched 3D SWC arbors",
            CommandAction::SwitchPreset("drosophila_flywire".to_string()),
        ),
        (
            "Preset",
            "Organism: C. elegans WormAtlas Wiring",
            "26 hermaphrodite neurons with chemical and gap junctions",
            CommandAction::SwitchPreset("celegans_wormatlas".to_string()),
        ),
        (
            "Surgery",
            "Heal All Ablated Regions & Severed Tracts",
            "Restore full intact structural connectome",
            CommandAction::HealAllLesions,
        ),
        (
            "Export",
            "Export Selected .SWC 3D Neuronal Arbor",
            "Write standardized 7-column NeuroMorpho/INCF .swc file",
            CommandAction::ExportSelectedSwc,
        ),
        (
            "Export",
            "Export Synchronized .connectome.json Envelope",
            "Write SHA-256 verified connectome-atlas-v1 JSON",
            CommandAction::ExportEnvelopeJson,
        ),
        (
            "Topology",
            "Comparative Connectome Graph-Topology Digest",
            "Inspect sigma, clustering C, path length L, rich-club Phi, and modularity Q across presets",
            CommandAction::InspectTopologyDigest,
        ),
        (
            "Export",
            "Export Comparative Topology Digest (.json)",
            "Write comparative small-world, rich-club, and modularity Q digest JSON",
            CommandAction::ExportTopologyDigestJson,
        ),
        (
            "Export",
            "Export Comparative Topology Digest (.csv)",
            "Write comparative small-world, rich-club, and modularity Q digest CSV",
            CommandAction::ExportTopologyDigestCsv,
        ),
    ];

    for (cat, title, sub, act) in builtins {
        let blob = format!("{} {} {}", cat.to_lowercase(), title.to_lowercase(), sub.to_lowercase());
        entries.push(CommandEntry {
            category: cat.to_string(),
            title: title.to_string(),
            subtitle: sub.to_string(),
            action: act,
            search_blob: blob,
        });
    }

    entries
}

/// Filters and ranks command palette entries by multi-token substring relevance.
pub fn search_command_index<'a>(
    entries: &'a [CommandEntry],
    query: &str,
    limit: usize,
) -> Vec<&'a CommandEntry> {
    let trimmed = query.trim().to_lowercase();
    if trimmed.is_empty() {
        return entries.iter().take(limit.max(1)).collect();
    }
    let tokens: Vec<&str> = trimmed.split_whitespace().collect();
    let mut scored: Vec<(i32, &CommandEntry)> = Vec::new();
    for entry in entries {
        let mut matched = true;
        let mut score = 0_i32;
        let title_lc = entry.title.to_lowercase();
        for &tk in &tokens {
            if let Some(pos) = entry.search_blob.find(tk) {
                if title_lc.starts_with(tk) {
                    score += 35;
                } else if pos < 16 {
                    score += 18;
                } else {
                    score += 6;
                }
            } else {
                matched = false;
                break;
            }
        }
        if matched {
            scored.push((score, entry));
        }
    }
    scored.sort_by(|a, b| b.0.cmp(&a.0));
    scored
        .into_iter()
        .take(limit.max(1))
        .map(|(_, e)| e)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::presets::load_preset_by_id;

    #[test]
    fn searches_regions_and_commands_cleanly() {
        let human = load_preset_by_id("human_mni152").expect("human preset");
        let idx = build_command_index(&human);
        let hits = search_command_index(&idx, "broca", 5);
        assert!(!hits.is_empty());
        assert!(hits[0].title.to_lowercase().contains("broca"));

        let zen_hits = search_command_index(&idx, "zen", 5);
        assert!(!zen_hits.is_empty());
        assert_eq!(zen_hits[0].action, CommandAction::ToggleZenMode);
    }

    #[test]
    fn command_index_includes_topology_digest_actions() {
        let human = load_preset_by_id("human_mni152").expect("human preset");
        let idx = build_command_index(&human);
        let hits = search_command_index(&idx, "modularity", 5);
        assert!(hits.iter().any(|h| matches!(h.action, CommandAction::InspectTopologyDigest)));
        let csv_hits = search_command_index(&idx, "digest csv", 5);
        assert!(csv_hits.iter().any(|h| matches!(h.action, CommandAction::ExportTopologyDigestCsv)));
        let json_hits = search_command_index(&idx, "digest json", 5);
        assert!(json_hits.iter().any(|h| matches!(h.action, CommandAction::ExportTopologyDigestJson)));
    }
}
