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

use atlas_core::PathwayModality;
use glam::Vec3;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ColorSchemeMode {
    FunctionalNetwork,
    DtiDirectionalRgb,
    BetweennessHeatmap,
    HarmonicModePsi,
    CriticalitySusceptibility,
    TdaBeta1Loops,
}

impl ColorSchemeMode {
    pub fn label(self) -> &'static str {
        match self {
            Self::FunctionalNetwork => "Functional Network",
            Self::DtiDirectionalRgb => "DTI Directional RGB",
            Self::BetweennessHeatmap => "Betweenness Heatmap",
            Self::HarmonicModePsi => "Harmonic Mode psi_k",
            Self::CriticalitySusceptibility => "Criticality sigma_i",
            Self::TdaBeta1Loops => "TDA beta-1 Loops",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum TractRenderStyle {
    BraidedFascicles,
    SingleStreamline,
}

/// Generates a smooth, anatomically arched 3D cubic Bezier streamline between two stereotaxic ROIs.
pub fn generate_arched_tract_curve(
    a: [f32; 3],
    b: [f32; 3],
    modality: PathwayModality,
    segments: usize,
) -> Vec<[f32; 3]> {
    let p0 = Vec3::from_array(a);
    let p3 = Vec3::from_array(b);
    let mid = (p0 + p3) * 0.5;
    let span = p0.distance(p3);

    let (p1, p2) = match modality {
        PathwayModality::Commissural => {
            let callosal_arch = Vec3::new(0.0, mid.y * 0.88, mid.z.max(12.0) + span * 0.14);
            (p0.lerp(callosal_arch, 0.65), p3.lerp(callosal_arch, 0.65))
        }
        PathwayModality::Association => {
            let inward = Vec3::new(mid.x * 0.78, mid.y, mid.z + span * 0.12);
            (p0.lerp(inward, 0.55), p3.lerp(inward, 0.55))
        }
        PathwayModality::Projection => {
            let capsule = Vec3::new(mid.x * 0.68, mid.y * 0.90, mid.z);
            (p0.lerp(capsule, 0.55), p3.lerp(capsule, 0.55))
        }
        PathwayModality::ChemicalSynapse | PathwayModality::ElectricalGapJunction => {
            let ring_pull = Vec3::new(mid.x * 0.82, mid.y, mid.z + span * 0.08);
            (p0.lerp(ring_pull, 0.50), p3.lerp(ring_pull, 0.50))
        }
    };

    let steps = segments.max(6);
    let mut pts = Vec::with_capacity(steps + 1);
    for i in 0..=steps {
        let t = (i as f32) / (steps as f32);
        let u = 1.0 - t;
        let pos = p0 * (u * u * u)
            + p1 * (3.0 * u * u * t)
            + p2 * (3.0 * u * t * t)
            + p3 * (t * t * t);
        pts.push(pos.to_array());
    }
    pts
}

/// Generates a multi-strand braided white-matter fascicle bundle with callosal waist compaction,
/// cortical gray-matter insertion fanning, and gentle helical strand braiding.
pub fn generate_braided_fascicle_bundle(
    a: [f32; 3],
    b: [f32; 3],
    modality: PathwayModality,
    weight: f64,
    strand_count: usize,
    segments: usize,
) -> Vec<Vec<[f32; 3]>> {
    let core = generate_arched_tract_curve(a, b, modality, segments);
    let n_strands = strand_count.clamp(1, 9);
    let mut bundle = Vec::with_capacity(n_strands);
    bundle.push(core.clone());
    if n_strands == 1 || core.len() < 3 {
        return bundle;
    }

    let spread_mm = 1.45 + (weight as f32).clamp(0.0, 1.0) * 1.35;
    let n_pts = core.len();

    for s in 1..n_strands {
        let base_angle =
            ((s - 1) as f32 / ((n_strands - 1) as f32)) * std::f32::consts::TAU;
        let mut strand = Vec::with_capacity(n_pts);
        for (k, &pt) in core.iter().enumerate() {
            let t = (k as f32) / ((n_pts - 1) as f32);
            let p_prev = Vec3::from_array(core[k.saturating_sub(1)]);
            let p_next = Vec3::from_array(core[(k + 1).min(n_pts - 1)]);
            let tang = (p_next - p_prev).normalize_or_zero();
            let ref_up = if tang.z.abs() < 0.88 {
                Vec3::Z
            } else {
                Vec3::X
            };
            let binorm = tang.cross(ref_up).normalize_or_zero();
            let norm = binorm.cross(tang).normalize_or_zero();

            let envelope = match modality {
                PathwayModality::Commissural => {
                    0.42 + 0.68 * (2.0 * t - 1.0).abs().powf(1.25)
                }
                _ => 0.32 + 0.78 * (std::f32::consts::PI * t).sin(),
            };
            let end_taper =
                (std::f32::consts::PI * (t.min(1.0 - t) * 4.5).min(1.0)).sin();
            let r_off = spread_mm * envelope * (0.28 + 0.72 * end_taper);
            let twist = base_angle + 1.05 * t;
            let pos = Vec3::from_array(pt)
                + norm * (twist.cos() * r_off)
                + binorm * (twist.sin() * r_off);
            strand.push(pos.to_array());
        }
        bundle.push(strand);
    }
    bundle
}
