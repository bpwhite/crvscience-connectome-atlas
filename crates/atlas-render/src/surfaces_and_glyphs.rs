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

use atlas_core::{jacobi_eigen_symmetric, BrainRegion, ConnectomeAtlas};
use glam::Vec3;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum SurgicalCutawayMode {
    Intact,
    AnteriorSuperiorOctant,
    LeftHemisphere,
    SuperiorAxialDome,
}

impl SurgicalCutawayMode {
    pub fn shader_code(self) -> u32 {
        match self {
            Self::Intact => 0,
            Self::AnteriorSuperiorOctant => 1,
            Self::LeftHemisphere => 2,
            Self::SuperiorAxialDome => 3,
        }
    }
}

/// Evaluates constructive solid geometry (CSG) surgical shell cutaway discard test in stereotaxic mm.
pub fn is_point_clipped_by_cutaway(
    pt: [f32; 3],
    mode: SurgicalCutawayMode,
    offset_mm: f32,
) -> bool {
    match mode {
        SurgicalCutawayMode::Intact => false,
        SurgicalCutawayMode::AnteriorSuperiorOctant => {
            let cx = 2.0 + offset_mm * 0.45;
            let cy = -6.0 - offset_mm * 0.55;
            let cz = 0.0 - offset_mm * 0.45;
            pt[0] < cx && pt[1] > cy && pt[2] > cz
        }
        SurgicalCutawayMode::LeftHemisphere => pt[0] < offset_mm,
        SurgicalCutawayMode::SuperiorAxialDome => pt[2] > (16.0 - offset_mm * 0.65),
    }
}

/// Projects parcellated ROI scalar values onto 3D cortical surface mesh vertices via Gaussian kernel
/// registration fusion (`neuromaps` continuous spatial projection).
pub fn project_parcel_scalars_to_surface(
    vertices: &[[f32; 3]],
    parcel_coords: &[[f32; 3]],
    parcel_values: &[f64],
    sigma_mm: f32,
) -> Vec<f32> {
    let inv_two_sigma_sq = 1.0 / (2.0 * sigma_mm.max(1.0).powi(2));
    vertices
        .iter()
        .map(|&v| {
            let mut w_sum = 0.0_f64;
            let mut acc = 0.0_f64;
            for (pc, &val) in parcel_coords.iter().zip(parcel_values.iter()) {
                let dx = v[0] - pc[0];
                let dy = v[1] - pc[1];
                let dz = v[2] - pc[2];
                let d2 = (dx * dx + dy * dy + dz * dz) as f64;
                let w = (-d2 * (inv_two_sigma_sq as f64)).exp();
                w_sum += w;
                acc += w * val;
            }
            if w_sum > 1e-8 {
                (acc / w_sum) as f32
            } else {
                0.0
            }
        })
        .collect()
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum OrganoidClass {
    CaudateArch,
    LentiformPutamen,
    CoiledHippocampus,
    AlmondAmygdala,
    PulvinarThalamus,
    FoliatedCerebellum,
    PontineColumn,
    CorticalGyralCrown,
    BranchedSwcArbor,
}

pub fn classify_region_organoid(region: &BrainRegion, organism_id: &str) -> OrganoidClass {
    if organism_id == "drosophila_flywire" || organism_id == "celegans_wormatlas" {
        return OrganoidClass::BranchedSwcArbor;
    }
    let id = region.id.as_str();
    if id.contains("Caud") || id.contains("STR.CP") {
        OrganoidClass::CaudateArch
    } else if id.contains("Put_") || id.contains("ACB") {
        OrganoidClass::LentiformPutamen
    } else if id.contains("Hipp") || id.contains("CA1") || id.contains("DG") {
        OrganoidClass::CoiledHippocampus
    } else if id.contains("Amyg") || id.contains("BLA") {
        OrganoidClass::AlmondAmygdala
    } else if id.contains("Thal") || id.contains("TH.") || id == "HY" {
        OrganoidClass::PulvinarThalamus
    } else if id.contains("CB.") || region.network.contains("Cerebell") {
        OrganoidClass::FoliatedCerebellum
    } else if id.contains("BS.") || id.contains("MB.") || id == "MY" {
        OrganoidClass::PontineColumn
    } else {
        OrganoidClass::CorticalGyralCrown
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LocalFiberTensor {
    pub region_id: String,
    pub eigenvalues: [f64; 3],
    pub principal_axis: [f32; 3],
    pub fractional_anisotropy: f64,
    pub westin_linear_cl: f64,
    pub westin_planar_cp: f64,
    pub westin_spherical_cs: f64,
    pub superquadric_alpha_beta: [f64; 2],
}

/// Computes local 3x3 symmetric fiber orientation tensors, Fractional Anisotropy (`FA`),
/// Westin shape metrics (`c_l`, `c_p`, `c_s`), and Kindlmann superquadric glyph exponents.
pub fn compute_atlas_fiber_tensors(
    atlas: &ConnectomeAtlas,
    min_weight: f64,
) -> Vec<LocalFiberTensor> {
    let mut out = Vec::with_capacity(atlas.regions.len());
    for reg in &atlas.regions {
        let mut t_mat = vec![
            vec![0.04_f64, 0.0, 0.0],
            vec![0.0, 0.04, 0.0],
            vec![0.0, 0.0, 0.04],
        ];
        for e in &atlas.edges {
            if !atlas.is_edge_intact(e, min_weight) {
                continue;
            }
            let peer_opt = if e.source == reg.id {
                atlas.regions.iter().find(|r| r.id == e.target)
            } else if e.target == reg.id {
                atlas.regions.iter().find(|r| r.id == e.source)
            } else {
                None
            };
            if let Some(peer) = peer_opt {
                let dir = (Vec3::from_array(peer.coord) - Vec3::from_array(reg.coord))
                    .normalize_or_zero();
                let u = [dir.x as f64, dir.y as f64, dir.z as f64];
                for r_i in 0..3 {
                    for c_i in 0..3 {
                        t_mat[r_i][c_i] += e.weight * u[r_i] * u[c_i];
                    }
                }
            }
        }

        let (evals_asc, evecs) = jacobi_eigen_symmetric(&t_mat, 24);
        let l1 = evals_asc[2].max(1e-5);
        let l2 = evals_asc[1].max(1e-5);
        let l3 = evals_asc[0].max(1e-5);
        let sum_l = l1 + l2 + l3;
        let mean_l = sum_l / 3.0;
        let denom = (l1 * l1 + l2 * l2 + l3 * l3).sqrt().max(1e-9);
        let fa = ((1.5_f64).sqrt()
            * ((l1 - mean_l).powi(2) + (l2 - mean_l).powi(2) + (l3 - mean_l).powi(2)).sqrt()
            / denom)
            .clamp(0.0, 1.0);
        let cl = (l1 - l2) / sum_l;
        let cp = 2.0 * (l2 - l3) / sum_l;
        let cs = 3.0 * l3 / sum_l;
        let gamma = 3.0_f64;
        let (alpha_exp, beta_exp) = if cl >= cp {
            ((1.0 - cp).max(0.06).powf(gamma), (1.0 - cl).max(0.06).powf(gamma))
        } else {
            ((1.0 - cl).max(0.06).powf(gamma), (1.0 - cp).max(0.06).powf(gamma))
        };

        out.push(LocalFiberTensor {
            region_id: reg.id.clone(),
            eigenvalues: [l1, l2, l3],
            principal_axis: [
                evecs[0][2] as f32,
                evecs[1][2] as f32,
                evecs[2][2] as f32,
            ],
            fractional_anisotropy: fa,
            westin_linear_cl: cl,
            westin_planar_cp: cp,
            westin_spherical_cs: cs,
            superquadric_alpha_beta: [alpha_exp, beta_exp],
        });
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use atlas_core::load_preset_by_id;

    #[test]
    fn computes_valid_westin_tensors_and_cutaways() {
        let human = load_preset_by_id("human_mni152").expect("human preset");
        let tensors = compute_atlas_fiber_tensors(&human, 0.0);
        assert_eq!(tensors.len(), human.regions.len());
        for t in &tensors {
            assert!((0.0..=1.0).contains(&t.fractional_anisotropy));
            let sum_westin = t.westin_linear_cl + t.westin_planar_cp + t.westin_spherical_cs;
            assert!((sum_westin - 1.0).abs() < 1e-6);
        }
        assert!(is_point_clipped_by_cutaway(
            [-24.0, 30.0, 28.0],
            SurgicalCutawayMode::AnteriorSuperiorOctant,
            0.0
        ));
        assert!(!is_point_clipped_by_cutaway(
            [28.0, -30.0, -12.0],
            SurgicalCutawayMode::AnteriorSuperiorOctant,
            0.0
        ));
    }
}
