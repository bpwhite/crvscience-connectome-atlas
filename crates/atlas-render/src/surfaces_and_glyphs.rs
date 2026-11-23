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
    #[serde(default)]
    pub secondary_axis: [f32; 3],
    #[serde(default)]
    pub tertiary_axis: [f32; 3],
    pub fractional_anisotropy: f64,
    pub westin_linear_cl: f64,
    pub westin_planar_cp: f64,
    pub westin_spherical_cs: f64,
    pub superquadric_alpha_beta: [f64; 2],
}

/// Shaded 3D ellipsoidal diffusion-tensor glyph mesh plus orthogonal principal-axis segments.
#[derive(Debug, Clone)]
pub struct TensorEllipsoidMesh {
    pub vertices: Vec<[f32; 3]>,
    pub normals: Vec<[f32; 3]>,
    pub triangles: Vec<[u32; 3]>,
    pub axis_segments: [([f32; 3], [f32; 3]); 3],
}

/// Generates a 3-axis ellipsoidal diffusion-tensor glyph aligned to `(e1, e2, e3)` and scaled by
/// normalized semi-axes `sqrt(lambda_k / lambda_1)` (with a minimum floor for visibility).
pub fn generate_tensor_ellipsoid_mesh(
    center: [f32; 3],
    tensor: &LocalFiberTensor,
    max_radius_mm: f32,
    lat_rings: usize,
    lon_sectors: usize,
) -> TensorEllipsoidMesh {
    let l1 = tensor.eigenvalues[0].max(1e-5);
    let l2 = tensor.eigenvalues[1].max(1e-5);
    let l3 = tensor.eigenvalues[2].max(1e-5);
    let r1 = max_radius_mm.max(1.0);
    let r2 = (r1 * ((l2 / l1).sqrt() as f32)).max(r1 * 0.24);
    let r3 = (r1 * ((l3 / l1).sqrt() as f32)).max(r1 * 0.20);

    let e1 = Vec3::from_array(tensor.principal_axis).normalize_or_zero();
    let e1 = if e1.length_squared() < 1e-6 { Vec3::X } else { e1 };
    let e2_raw = Vec3::from_array(tensor.secondary_axis);
    let e2 = (e2_raw - e1 * e2_raw.dot(e1)).normalize_or_zero();
    let e2 = if e2.length_squared() < 1e-6 {
        let ref_v = if e1.z.abs() < 0.85 { Vec3::Z } else { Vec3::Y };
        e1.cross(ref_v).normalize()
    } else {
        e2
    };
    let e3 = e1.cross(e2).normalize_or_zero();
    let c = Vec3::from_array(center);

    let rings = lat_rings.max(3);
    let sectors = lon_sectors.max(6);
    let mut vertices = Vec::with_capacity((rings + 1) * sectors);
    let mut normals = Vec::with_capacity((rings + 1) * sectors);
    let mut triangles = Vec::with_capacity(rings * sectors * 2);

    for i in 0..=rings {
        let v = (i as f32) / (rings as f32);
        let phi = (v - 0.5) * std::f32::consts::PI;
        let cp = phi.cos();
        let sp = phi.sin();
        for j in 0..sectors {
            let theta = (j as f32 / sectors as f32) * std::f32::consts::TAU;
            let ct = theta.cos();
            let st = theta.sin();
            let lx = r1 * cp * ct;
            let ly = r2 * cp * st;
            let lz = r3 * sp;
            let pos = c + e1 * lx + e2 * ly + e3 * lz;
            let nrm = (e1 * (lx / (r1 * r1)) + e2 * (ly / (r2 * r2)) + e3 * (lz / (r3 * r3)))
                .normalize_or_zero();
            vertices.push(pos.to_array());
            normals.push(nrm.to_array());
        }
    }

    for i in 0..rings {
        let row0 = (i * sectors) as u32;
        let row1 = ((i + 1) * sectors) as u32;
        for j in 0..sectors {
            let j0 = j as u32;
            let j1 = ((j + 1) % sectors) as u32;
            triangles.push([row0 + j0, row0 + j1, row1 + j1]);
            triangles.push([row0 + j0, row1 + j1, row1 + j0]);
        }
    }

    let axis_segments = [
        ((c - e1 * r1).to_array(), (c + e1 * r1).to_array()),
        ((c - e2 * r2).to_array(), (c + e2 * r2).to_array()),
        ((c - e3 * r3).to_array(), (c + e3 * r3).to_array()),
    ];

    TensorEllipsoidMesh {
        vertices,
        normals,
        triangles,
        axis_segments,
    }
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
            secondary_axis: [
                evecs[0][1] as f32,
                evecs[1][1] as f32,
                evecs[2][1] as f32,
            ],
            tertiary_axis: [
                evecs[0][0] as f32,
                evecs[1][0] as f32,
                evecs[2][0] as f32,
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
            let e1 = Vec3::from_array(t.principal_axis);
            let e2 = Vec3::from_array(t.secondary_axis);
            let e3 = Vec3::from_array(t.tertiary_axis);
            assert!(e1.dot(e2).abs() < 1e-4 && e1.dot(e3).abs() < 1e-4 && e2.dot(e3).abs() < 1e-4);
        }
        let glyph = generate_tensor_ellipsoid_mesh([10.0, 20.0, 30.0], &tensors[0], 6.0, 4, 8);
        assert!(!glyph.vertices.is_empty());
        assert!(!glyph.triangles.is_empty());
        assert_eq!(glyph.axis_segments.len(), 3);
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
