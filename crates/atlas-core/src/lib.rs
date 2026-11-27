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

#![forbid(unsafe_code)]

pub mod command_index;
pub mod control_harmonics;
pub mod criticality_tda;
pub mod metrics;
pub mod model;
pub mod presets;
pub mod systems_neuro;

pub use control_harmonics::{
    jacobi_eigen_symmetric, node_controllability_from_spectrum, normalize_adjacency_schur,
    ControlHarmonicsReport, NodeControlHarmonicMetrics,
};
pub use metrics::{
    compute_all_presets_topology_digest, compute_comparative_topology_digest, AtlasSummaryMetrics,
    ComparativeTopologyDigest, LesionDeltaReport, NodeTopologyMetrics, PathResult,
    TopologyDigestEntry,
};
pub use model::{
    AtlasEnvelope, AtlasScale, BrainRegion, ConnectomeAtlas, DecColor, Hemisphere, PathwayEdge,
    PathwayModality, ReceptorProfile,
};
pub use presets::{all_presets, load_preset_by_id};

pub use command_index::{build_command_index, search_command_index, CommandAction, CommandEntry};
pub use criticality_tda::{
    CriticalityRegime, CriticalityReport, KappaDynamicRangePoint, NodeCriticalityProfile,
    NodeTopologyTdaProfile, PersistentCavityBar, ResponseCurvePoint, TopologyTdaReport,
};
pub use systems_neuro::{
    region_transcriptomic_profile, DomainDecodeEntry, NodeSystemsProfile, NullGenerativeBenchmark,
    SystemsNeuroReport, TranscriptomicProfile,
};

#[cfg(test)]
mod workspace_preset_tests {
    use super::*;

    #[test]
    fn all_four_presets_compute_metrics_and_envelope_roundtrip() {
        let presets = all_presets();
        assert_eq!(presets.len(), 4);
        for atlas in presets {
            let metrics = atlas.compute_summary_metrics(0.0);
            assert!(metrics.global_efficiency.is_finite());
            let env = AtlasEnvelope::wrap(atlas.clone()).expect("AtlasEnvelope::wrap succeeds");
            let json = serde_json::to_string(&env).expect("Serialize AtlasEnvelope");
            let decoded: AtlasEnvelope =
                serde_json::from_str(&json).expect("Deserialize AtlasEnvelope");
            assert_eq!(decoded.sha256_digest, env.sha256_digest);
            assert_eq!(decoded.atlas.id, atlas.id);
        }
    }
}
