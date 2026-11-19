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

use glam::Vec3;
use serde::{Deserialize, Serialize};
use std::f32::consts::{PI, TAU};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum AnatomicalLayer {
    OuterShell,
    InternalOrgan,
}

#[derive(Debug, Clone)]
pub struct AnatomicalMesh {
    pub name: &'static str,
    pub layer: AnatomicalLayer,
    pub base_rgb: [u8; 3],
    pub rim_rgb: [u8; 3],
    pub opacity_scale: f32,
    pub clippable: bool,
    pub vertices: Vec<[f32; 3]>,
    pub normals: Vec<[f32; 3]>,
    /// Normalized relief scalar in [-1.0, 1.0]: negative in sulci/folia/annuli for ambient occlusion.
    pub crease: Vec<f32>,
    pub triangles: Vec<[u32; 3]>,
}

#[derive(Debug, Clone)]
pub struct AnatomicalLandmark {
    pub label: &'static str,
    pub coord: [f32; 3],
    pub color: [u8; 3],
}

#[derive(Debug, Clone)]
pub struct SliceContour {
    pub points: Vec<[f32; 2]>,
    pub fill_rgba: [u8; 4],
    pub stroke_rgba: [u8; 4],
    pub stroke_width: f32,
    pub closed: bool,
}

#[derive(Debug, Clone)]
pub struct SpeciesAnatomyScene {
    pub atlas_id: String,
    pub meshes: Vec<AnatomicalMesh>,
    pub landmarks: Vec<AnatomicalLandmark>,
    pub default_target: [f32; 3],
    pub default_distance: f32,
    pub default_yaw: f32,
    pub default_pitch: f32,
    pub default_slice: [f32; 3],
    pub slice_extent_mm: f32,
}

fn recompute_smooth_normals(vertices: &[[f32; 3]], triangles: &[[u32; 3]]) -> Vec<[f32; 3]> {
    let mut acc = vec![Vec3::ZERO; vertices.len()];
    for tri in triangles {
        let i0 = tri[0] as usize;
        let i1 = tri[1] as usize;
        let i2 = tri[2] as usize;
        let p0 = Vec3::from_array(vertices[i0]);
        let p1 = Vec3::from_array(vertices[i1]);
        let p2 = Vec3::from_array(vertices[i2]);
        let fnrm = (p1 - p0).cross(p2 - p0);
        acc[i0] += fnrm;
        acc[i1] += fnrm;
        acc[i2] += fnrm;
    }
    acc.into_iter()
        .map(|n| {
            let len = n.length();
            if len > 1e-6 {
                (n / len).to_array()
            } else {
                [0.0, 0.0, 1.0]
            }
        })
        .collect()
}

fn verify_outward_winding(
    vertices: &[[f32; 3]],
    triangles: &mut [[u32; 3]],
    centroid: Vec3,
) {
    let mut signed_sum = 0.0_f64;
    for tri in triangles.iter() {
        let p0 = Vec3::from_array(vertices[tri[0] as usize]);
        let p1 = Vec3::from_array(vertices[tri[1] as usize]);
        let p2 = Vec3::from_array(vertices[tri[2] as usize]);
        let fnrm = (p1 - p0).cross(p2 - p0);
        let tri_c = (p0 + p1 + p2) / 3.0;
        signed_sum += fnrm.dot(tri_c - centroid) as f64;
    }
    if signed_sum < 0.0 {
        for tri in triangles.iter_mut() {
            tri.swap(1, 2);
        }
    }
}

fn build_sculpted_sphere<F>(
    name: &'static str,
    layer: AnatomicalLayer,
    base_rgb: [u8; 3],
    rim_rgb: [u8; 3],
    opacity_scale: f32,
    clippable: bool,
    rings: usize,
    sectors: usize,
    centroid: Vec3,
    mut deform: F,
) -> AnatomicalMesh
where
    F: FnMut(Vec3) -> ([f32; 3], f32),
{
    let mut vertices = Vec::with_capacity((rings + 1) * (sectors + 1));
    let mut crease = Vec::with_capacity((rings + 1) * (sectors + 1));
    let mut triangles = Vec::with_capacity(rings * sectors * 2);

    for r in 0..=rings {
        let v = r as f32 / rings as f32;
        let theta = v * PI;
        let sin_t = theta.sin();
        let cos_t = theta.cos();
        for s in 0..=sectors {
            let u = s as f32 / sectors as f32;
            let phi = u * TAU;
            let dir = Vec3::new(sin_t * phi.cos(), sin_t * phi.sin(), cos_t);
            let (pt, cr) = deform(dir);
            vertices.push(pt);
            crease.push(cr);
        }
    }

    let stride = sectors + 1;
    for r in 0..rings {
        for s in 0..sectors {
            let i00 = (r * stride + s) as u32;
            let i10 = ((r + 1) * stride + s) as u32;
            let i01 = (r * stride + s + 1) as u32;
            let i11 = ((r + 1) * stride + s + 1) as u32;
            if r > 0 {
                triangles.push([i00, i10, i01]);
            }
            if r + 1 < rings {
                triangles.push([i01, i10, i11]);
            }
        }
    }

    verify_outward_winding(&vertices, &mut triangles, centroid);
    let normals = recompute_smooth_normals(&vertices, &triangles);
    AnatomicalMesh {
        name,
        layer,
        base_rgb,
        rim_rgb,
        opacity_scale,
        clippable,
        vertices,
        normals,
        crease,
        triangles,
    }
}

/// Sweeps a closed elliptical tube along `spine` `(pos, rx_local, rz_local, crease)` with smooth end caps.
fn build_swept_tube(
    name: &'static str,
    layer: AnatomicalLayer,
    base_rgb: [u8; 3],
    rim_rgb: [u8; 3],
    opacity_scale: f32,
    clippable: bool,
    spine: &[([f32; 3], f32, f32, f32)],
    radial_segs: usize,
) -> AnatomicalMesh {
    let mut vertices = Vec::new();
    let mut crease = Vec::new();
    let mut triangles = Vec::new();
    if spine.len() < 2 {
        return AnatomicalMesh {
            name,
            layer,
            base_rgb,
            rim_rgb,
            opacity_scale,
            clippable,
            vertices,
            normals: Vec::new(),
            crease,
            triangles,
        };
    }

    let mut prev_right = Vec3::X;
    for (idx, &(pos_arr, rx, rz, cr)) in spine.iter().enumerate() {
        let pos = Vec3::from_array(pos_arr);
        let tang = if idx == 0 {
            Vec3::from_array(spine[1].0) - pos
        } else if idx + 1 == spine.len() {
            pos - Vec3::from_array(spine[idx - 1].0)
        } else {
            Vec3::from_array(spine[idx + 1].0) - Vec3::from_array(spine[idx - 1].0)
        }
        .normalize_or_zero();

        let ref_up = if tang.z.abs() < 0.88 {
            Vec3::Z
        } else {
            Vec3::Y
        };
        let mut right = tang.cross(ref_up).normalize_or_zero();
        if right.length_squared() < 1e-4 {
            right = prev_right;
        } else if right.dot(prev_right) < 0.0 && idx > 0 {
            right = -right;
        }
        prev_right = right;
        let up = right.cross(tang).normalize_or_zero();

        for s in 0..=radial_segs {
            let phi = (s as f32 / radial_segs as f32) * TAU;
            let p = pos + right * (phi.cos() * rx) + up * (phi.sin() * rz);
            vertices.push(p.to_array());
            crease.push(cr);
        }
    }

    let stride = radial_segs + 1;
    for i in 0..(spine.len() - 1) {
        for s in 0..radial_segs {
            let i00 = (i * stride + s) as u32;
            let i01 = (i * stride + s + 1) as u32;
            let i10 = ((i + 1) * stride + s) as u32;
            let i11 = ((i + 1) * stride + s + 1) as u32;
            triangles.push([i00, i01, i10]);
            triangles.push([i01, i11, i10]);
        }
    }

    // End caps
    let start_center_idx = vertices.len() as u32;
    vertices.push(spine[0].0);
    crease.push(spine[0].3);
    for s in 0..radial_segs {
        triangles.push([start_center_idx, (s + 1) as u32, s as u32]);
    }
    let last_ring_base = ((spine.len() - 1) * stride) as u32;
    let end_center_idx = vertices.len() as u32;
    let last = spine.len() - 1;
    vertices.push(spine[last].0);
    crease.push(spine[last].3);
    for s in 0..radial_segs {
        let su = s as u32;
        triangles.push([end_center_idx, last_ring_base + su, last_ring_base + su + 1]);
    }

    let centroid = spine
        .iter()
        .fold(Vec3::ZERO, |acc, s| acc + Vec3::from_array(s.0))
        / (spine.len() as f32);
    verify_outward_winding(&vertices, &mut triangles, centroid);
    let normals = recompute_smooth_normals(&vertices, &triangles);

    AnatomicalMesh {
        name,
        layer,
        base_rgb,
        rim_rgb,
        opacity_scale,
        clippable,
        vertices,
        normals,
        crease,
        triangles,
    }
}

fn build_torus(
    name: &'static str,
    layer: AnatomicalLayer,
    base_rgb: [u8; 3],
    rim_rgb: [u8; 3],
    opacity_scale: f32,
    center: Vec3,
    axis_u: Vec3,
    axis_v: Vec3,
    major_r: f32,
    minor_r: f32,
    major_segs: usize,
    minor_segs: usize,
) -> AnatomicalMesh {
    let mut vertices = Vec::with_capacity((major_segs + 1) * (minor_segs + 1));
    let mut normals = Vec::with_capacity((major_segs + 1) * (minor_segs + 1));
    let mut crease = Vec::with_capacity((major_segs + 1) * (minor_segs + 1));
    let mut triangles = Vec::with_capacity(major_segs * minor_segs * 2);
    let normal_plane = axis_u.cross(axis_v).normalize_or_zero();

    for i in 0..=major_segs {
        let u = (i as f32 / major_segs as f32) * TAU;
        let radial = axis_u * u.cos() + axis_v * u.sin();
        let ring_center = center + radial * major_r;
        for j in 0..=minor_segs {
            let v = (j as f32 / minor_segs as f32) * TAU;
            let n = (radial * v.cos() + normal_plane * v.sin()).normalize_or_zero();
            vertices.push((ring_center + n * minor_r).to_array());
            normals.push(n.to_array());
            crease.push(0.15 + 0.15 * v.cos());
        }
    }

    let stride = minor_segs + 1;
    for i in 0..major_segs {
        for j in 0..minor_segs {
            let i00 = (i * stride + j) as u32;
            let i10 = ((i + 1) * stride + j) as u32;
            let i01 = (i * stride + j + 1) as u32;
            let i11 = ((i + 1) * stride + j + 1) as u32;
            triangles.push([i00, i01, i10]);
            triangles.push([i01, i11, i10]);
        }
    }

    AnatomicalMesh {
        name,
        layer,
        base_rgb,
        rim_rgb,
        opacity_scale,
        clippable: false,
        vertices,
        normals,
        crease,
        triangles,
    }
}

fn gyral_fold_relief(p: Vec3) -> f32 {
    let g1 = (p.x * 0.24 + p.y * 0.19).sin() * (p.z * 0.22 + p.y * 0.15).cos();
    let g2 = (p.y * 0.31 - p.z * 0.27 + p.x * 0.14).sin();
    let g3 = (p.z * 0.34 + p.x * 0.21).cos() * 0.45;
    (g1 * 0.50 + g2 * 0.35 + g3 * 0.30).clamp(-1.0, 1.0)
}

fn build_human_mni152_scene() -> SpeciesAnatomyScene {
    let mut meshes = Vec::new();

    for &sign in &[-1.0_f32, 1.0_f32] {
        let name = if sign < 0.0 {
            "Left Cerebral Hemisphere"
        } else {
            "Right Cerebral Hemisphere"
        };
        let cx = sign * 34.5;
        let mesh = build_sculpted_sphere(
            name,
            AnatomicalLayer::OuterShell,
            [144, 178, 218],
            [186, 230, 253],
            1.0,
            true,
            26,
            34,
            Vec3::new(cx, -16.0, 15.0),
            |dir| {
                // Medial hemisphere wall sits near the interhemispheric fissure (|x| >= 2.2 mm)
                let medial_pull = dir.x * sign;
                let rx_lat = if medial_pull < 0.0 { 32.0 } else { 37.5 };
                let x_local = dir.x * rx_lat;
                let mut x = cx + x_local;
                if sign < 0.0 && x > -2.2 {
                    x = -2.2 + (x + 2.2) * 0.12;
                } else if sign > 0.0 && x < 2.2 {
                    x = 2.2 + (x - 2.2) * 0.12;
                }

                // Anterior-Posterior asymmetry: longer occipital pole, tapered frontal pole
                let ry = if dir.y > 0.0 { 76.0 } else { 87.0 };
                let mut y = -14.0 + dir.y * ry;
                let mut z = 14.0 + dir.z * if dir.z > 0.0 { 59.0 } else { 44.0 };

                // Temporal lobe downward/lateral convexity + orbital frontal upward scoop
                let temp_bump = (-((dir.y + 0.10) / 0.38).powi(2)
                    - ((dir.z + 0.45) / 0.38).powi(2))
                .exp()
                    * 16.0;
                z -= temp_bump * 0.85;
                x += sign * temp_bump * 0.50;

                // Orbital frontal hollow & occipital cerebellar notch (tentorium)
                if dir.y > 0.35 && dir.z < -0.1 {
                    z += (dir.y - 0.35) * (-dir.z) * 26.0;
                }
                if dir.y < -0.22 && dir.z < -0.15 {
                    let tent = ((-dir.y - 0.22) * (-dir.z - 0.15)).min(0.6);
                    z += tent * 34.0;
                }

                // Sylvian fissure indentation along lateral frontotemporal seam
                let sylvian = (-((dir.z + 0.08 - dir.y * 0.18) / 0.12).powi(2)).exp()
                    * (medial_pull.max(0.0))
                    * ((-((dir.y - 0.15) / 0.55).powi(2)).exp());
                x -= sign * sylvian * 6.5;

                // Cortical gyri & sulci folding (muted on medial fissure wall)
                let gyri = gyral_fold_relief(Vec3::new(x, y, z));
                let lat_factor = ((x.abs() - 3.0) / 14.0).clamp(0.25, 1.0);
                let disp = gyri * 2.5 * lat_factor;
                x += dir.x * disp;
                y += dir.y * disp;
                z += dir.z * disp;

                ([x, y, z], gyri * lat_factor)
            },
        );
        meshes.push(mesh);
    }

    // Corpus Callosum midline arch (genu -> body -> splenium)
    let cc_spine: Vec<([f32; 3], f32, f32, f32)> = (0..=16)
        .map(|i| {
            let t = i as f32 / 16.0;
            let y = 37.0 * (1.0 - t) - 44.0 * t;
            let arch = (t * PI).sin();
            let z = 3.0 + arch * 23.0 + t * 11.0;
            let thick_w = 11.5 - arch * 2.0;
            let thick_h = 4.2 + (1.0 - arch) * 2.4;
            ([0.0, y, z], thick_w, thick_h, 0.25)
        })
        .collect();
    meshes.push(build_swept_tube(
        "Corpus Callosum",
        AnatomicalLayer::InternalOrgan,
        [232, 204, 134],
        [254, 243, 199],
        0.92,
        false,
        &cc_spine,
        14,
    ));

    // Foliated Cerebellum (bilateral hemispheres + midline vermis)
    for &sign in &[-1.0_f32, 0.0_f32, 1.0_f32] {
        let is_vermis = sign == 0.0;
        let cx = sign * 25.0;
        let cy = -57.0 - if is_vermis { 3.0 } else { 0.0 };
        let cz = -30.5 + if is_vermis { 1.5 } else { 0.0 };
        let mesh = build_sculpted_sphere(
            "Foliated Cerebellum",
            AnatomicalLayer::InternalOrgan,
            [174, 140, 218],
            [221, 214, 254],
            0.90,
            false,
            18,
            22,
            Vec3::new(cx, cy, cz),
            |dir| {
                let rx = if is_vermis { 9.0 } else { 22.0 };
                let ry = if is_vermis { 15.0 } else { 18.5 };
                let rz = if is_vermis { 13.0 } else { 15.5 };
                let folia = (dir.z * 26.0).sin();
                let mut p = Vec3::new(cx + dir.x * rx, cy + dir.y * ry, cz + dir.z * rz);
                p += dir * (folia * 0.65);
                (p.to_array(), folia * 0.9)
            },
        );
        meshes.push(mesh);
    }

    // Brainstem (Midbrain + Ventral Pons Bulge + Medulla)
    let bs_spine = [
        ([0.0, -14.0, -3.0], 11.0, 9.5, 0.2),
        ([0.0, -18.0, -14.0], 10.5, 9.5, 0.1),
        ([0.0, -20.5, -24.0], 13.8, 12.8, 0.35), // Ventral Pons bulge
        ([0.0, -25.0, -34.0], 10.2, 9.5, 0.1),
        ([0.0, -29.0, -45.0], 8.2, 8.0, 0.0),
        ([0.0, -32.0, -56.0], 6.8, 6.8, -0.1),
    ];
    meshes.push(build_swept_tube(
        "Brainstem & Pons",
        AnatomicalLayer::InternalOrgan,
        [152, 174, 202],
        [203, 213, 225],
        0.94,
        false,
        &bs_spine,
        14,
    ));

    // Deep Subcortical Structures: Thalamus, Caudate Arch, Lentiform Putamen, Hippocampus, Amygdala
    for &sign in &[-1.0_f32, 1.0_f32] {
        // Thalamus ovoid
        let tc = Vec3::new(sign * 13.5, -21.5, -2.0);
        meshes.push(build_sculpted_sphere(
            "Thalamus",
            AnatomicalLayer::InternalOrgan,
            [236, 118, 134],
            [254, 205, 211],
            0.88,
            false,
            12,
            16,
            tc,
            |d| ((tc + Vec3::new(d.x * 9.5, d.y * 14.5, d.z * 9.0)).to_array(), 0.25),
        ));

        // C-shaped Hippocampus curving down into temporal lobe
        let hipp_spine: Vec<([f32; 3], f32, f32, f32)> = (0..=10)
            .map(|i| {
                let t = i as f32 / 10.0;
                let x = sign * (24.0 + (t * PI).sin() * 12.0);
                let y = -35.0 + t * 28.0;
                let z = 4.0 - (t * PI * 0.85).sin() * 26.0;
                let r = 3.2 + (t * PI).sin() * 1.4;
                ([x, y, z], r, r * 0.85, 0.2)
            })
            .collect();
        meshes.push(build_swept_tube(
            "Hippocampus",
            AnatomicalLayer::InternalOrgan,
            [244, 136, 178],
            [251, 207, 232],
            0.90,
            false,
            &hipp_spine,
            10,
        ));

        // Caudate arch
        let caud_spine: Vec<([f32; 3], f32, f32, f32)> = (0..=10)
            .map(|i| {
                let t = i as f32 / 10.0;
                let x = sign * (13.0 + t * 4.5);
                let y = 18.0 - t * 44.0;
                let z = 6.0 + (t * PI).sin() * 11.5;
                let r = 4.4 * (1.0 - t * 0.55);
                ([x, y, z], r, r, 0.2)
            })
            .collect();
        meshes.push(build_swept_tube(
            "Caudate Nucleus",
            AnatomicalLayer::InternalOrgan,
            [74, 212, 198],
            [167, 243, 208],
            0.86,
            false,
            &caud_spine,
            10,
        ));

        // Lentiform Putamen + Amygdala tip
        let pc = Vec3::new(sign * 25.5, -2.5, -4.0);
        meshes.push(build_sculpted_sphere(
            "Putamen",
            AnatomicalLayer::InternalOrgan,
            [56, 189, 178],
            [153, 246, 228],
            0.84,
            false,
            10,
            14,
            pc,
            |d| ((pc + Vec3::new(d.x * 6.5, d.y * 12.5, d.z * 7.5)).to_array(), 0.15),
        ));
    }

    SpeciesAnatomyScene {
        atlas_id: "human_mni152".to_string(),
        meshes,
        landmarks: vec![
            AnatomicalLandmark { label: "Frontal Pole", coord: [26.0, 58.0, 8.0], color: [186, 230, 253] },
            AnatomicalLandmark { label: "Parietal Crown", coord: [24.0, -28.0, 68.0], color: [186, 230, 253] },
            AnatomicalLandmark { label: "Temporal Lobe", coord: [-54.0, 8.0, -26.0], color: [186, 230, 253] },
            AnatomicalLandmark { label: "Occipital Pole", coord: [24.0, -92.0, 10.0], color: [186, 230, 253] },
            AnatomicalLandmark { label: "Corpus Callosum", coord: [0.0, 6.0, 28.0], color: [253, 230, 138] },
            AnatomicalLandmark { label: "Thalamus", coord: [-15.0, -21.0, 4.0], color: [254, 205, 211] },
            AnatomicalLandmark { label: "Cerebellum", coord: [-32.0, -64.0, -34.0], color: [221, 214, 254] },
            AnatomicalLandmark { label: "Pons / Brainstem", coord: [0.0, -24.0, -52.0], color: [203, 213, 225] },
        ],
        default_target: [0.0, -14.0, 8.0],
        default_distance: 238.0,
        default_yaw: -0.68,
        default_pitch: 0.30,
        default_slice: [0.0, -14.0, 10.0],
        slice_extent_mm: 105.0,
    }
}

fn build_mouse_ccfv3_scene() -> SpeciesAnatomyScene {
    let mut meshes = Vec::new();

    // Smooth lissencephalic rodent hemispheres (tapering rostrally)
    for &sign in &[-1.0_f32, 1.0_f32] {
        let cx = sign * 21.5;
        let mesh = build_sculpted_sphere(
            "Rodent Isocortex",
            AnatomicalLayer::OuterShell,
            [138, 186, 218],
            [186, 230, 253],
            1.0,
            true,
            20,
            26,
            Vec3::new(cx, -2.0, 2.0),
            |dir| {
                // Narrower rostrally (anterior taper) & wider caudally
                let taper = (1.0 - 0.32 * dir.y.max(-0.4)).clamp(0.62, 1.22);
                let mut x = cx + dir.x * 20.5 * taper;
                if sign < 0.0 && x > -1.2 {
                    x = -1.2;
                } else if sign > 0.0 && x < 1.2 {
                    x = 1.2;
                }
                let y = -1.0 + dir.y * 39.0;
                let z = 1.5 + dir.z * 21.5;
                // Subtle lateral rhinal fissure crease (lissencephalic smooth brain!)
                let rhinal = (-((dir.z + 0.30) / 0.14).powi(2)).exp() * (dir.x * sign).max(0.0);
                x -= sign * rhinal * 2.2;
                ([x, y, z], -rhinal * 0.8 + 0.15)
            },
        );
        meshes.push(mesh);

        // Prominent paired Olfactory Bulbs (MOB) nestled continuously against rostral frontal poles
        let mc = Vec3::new(sign * 10.5, 41.5, -3.5);
        meshes.push(build_sculpted_sphere(
            "Olfactory Bulb",
            AnatomicalLayer::OuterShell,
            [220, 156, 124],
            [254, 215, 170],
            0.96,
            false,
            14,
            18,
            mc,
            |d| {
                let taper = 1.0 - d.y.max(0.0) * 0.35;
                (
                    (mc + Vec3::new(d.x * 8.6 * taper, d.y * 13.5, d.z * 7.8 * taper)).to_array(),
                    0.22,
                )
            },
        ));

        // Dorsal C-curved Hippocampus (CA1/DG)
        let hc_spine: Vec<([f32; 3], f32, f32, f32)> = (0..=8)
            .map(|i| {
                let t = i as f32 / 8.0;
                let x = sign * (14.0 + t * 10.0);
                let y = -14.0 - t * 14.0;
                let z = 9.0 - t * 12.0;
                ([x, y, z], 4.2, 3.6, 0.2)
            })
            .collect();
        meshes.push(build_swept_tube(
            "Mouse Hippocampus",
            AnatomicalLayer::InternalOrgan,
            [244, 114, 182],
            [251, 207, 232],
            0.88,
            false,
            &hc_spine,
            10,
        ));
    }

    // Large Posterior Foliated Cerebellum + Paraflocculi
    let cb_center = Vec3::new(0.0, -49.5, 0.5);
    meshes.push(build_sculpted_sphere(
        "Mouse Cerebellum",
        AnatomicalLayer::InternalOrgan,
        [168, 138, 224],
        [221, 214, 254],
        0.92,
        false,
        18,
        24,
        cb_center,
        |d| {
            let folia = (d.y * 22.0 + d.z * 14.0).sin();
            let p = cb_center
                + Vec3::new(d.x * 32.0, d.y * 13.5, d.z * 15.0)
                + d * (folia * 0.55);
            (p.to_array(), folia * 0.9)
        },
    ));

    // Caudal Brainstem / Medulla
    let m_bs = [
        ([0.0, -28.0, -10.0], 11.0, 8.5, 0.15),
        ([0.0, -46.0, -13.0], 9.5, 8.0, 0.1),
        ([0.0, -64.0, -16.0], 7.0, 6.2, -0.05),
    ];
    meshes.push(build_swept_tube(
        "Mouse Brainstem",
        AnatomicalLayer::InternalOrgan,
        [148, 172, 200],
        [203, 213, 225],
        0.92,
        false,
        &m_bs,
        12,
    ));

    SpeciesAnatomyScene {
        atlas_id: "mouse_ccfv3".to_string(),
        meshes,
        landmarks: vec![
            AnatomicalLandmark { label: "Olfactory Bulbs (MOB)", coord: [-12.0, 52.0, 4.0], color: [254, 215, 170] },
            AnatomicalLandmark { label: "Smooth Isocortex", coord: [22.0, 10.0, 22.0], color: [186, 230, 253] },
            AnatomicalLandmark { label: "Hippocampus (CA1)", coord: [-19.0, -20.0, 9.0], color: [251, 207, 232] },
            AnatomicalLandmark { label: "Caudoputamen (CP)", coord: [-22.0, 6.0, 0.0], color: [153, 246, 228] },
            AnatomicalLandmark { label: "Foliated Cerebellum", coord: [0.0, -50.0, 15.0], color: [221, 214, 254] },
            AnatomicalLandmark { label: "Medulla / Brainstem", coord: [0.0, -62.0, -16.0], color: [203, 213, 225] },
        ],
        default_target: [-2.0, -2.0, 1.0],
        default_distance: 162.0,
        default_yaw: -0.78,
        default_pitch: 0.34,
        default_slice: [-12.0, 0.0, 4.0],
        slice_extent_mm: 68.0,
    }
}

fn build_drosophila_flywire_scene() -> SpeciesAnatomyScene {
    let mut meshes = Vec::new();

    // Central Supraesophageal Brain Shell (bilobed lateral protocerebrum)
    meshes.push(build_sculpted_sphere(
        "Central Brain Capsule",
        AnatomicalLayer::OuterShell,
        [124, 178, 228],
        [186, 230, 253],
        0.95,
        true,
        22,
        28,
        Vec3::new(0.0, 2.0, 4.0),
        |d| {
            let lobe_pinch = 1.0 - 0.16 * (-((d.x) / 0.25).powi(2)).exp();
            let p = Vec3::new(
                d.x * 38.0,
                2.0 + d.y * 23.0 * lobe_pinch,
                4.0 + d.z * 24.5 * lobe_pinch,
            );
            (p.to_array(), 0.2)
        },
    ));

    // Ventral Gnathal Ganglia (GNG) subesophageal wedge
    meshes.push(build_sculpted_sphere(
        "Gnathal Ganglia (GNG)",
        AnatomicalLayer::OuterShell,
        [96, 196, 174],
        [167, 243, 208],
        0.92,
        false,
        14,
        18,
        Vec3::new(0.0, 11.0, -28.0),
        |d| {
            let p = Vec3::new(d.x * 18.0, 11.0 + d.y * 13.5, -27.5 + d.z * 11.5);
            (p.to_array(), 0.15)
        },
    ));

    // Prominent Bilateral Optic Lobes (tapered Lobula neck + rounded kidney-bean Medulla cup)
    for &sign in &[-1.0_f32, 1.0_f32] {
        // Lobula complex bridge narrowing smoothly outward
        let lo_spine = [
            ([sign * 34.0, -5.0, 2.5], 12.5, 14.5, 0.15),
            ([sign * 45.0, -7.0, 3.5], 11.5, 14.5, 0.20),
            ([sign * 56.0, -6.0, 4.0], 13.0, 17.5, 0.25),
        ];
        meshes.push(build_swept_tube(
            "Optic Lobe Neck",
            AnatomicalLayer::OuterShell,
            [154, 138, 232],
            [221, 214, 254],
            0.92,
            false,
            &lo_spine,
            14,
        ));

        // Rounded crescent Medulla/Lamina lobe dome curving laterally
        let center = Vec3::new(sign * 59.0, -5.5, 4.0);
        meshes.push(build_sculpted_sphere(
            "Optic Lobe Dome (Medulla/Lamina)",
            AnatomicalLayer::OuterShell,
            [164, 146, 240],
            [224, 215, 255],
            0.98,
            false,
            18,
            22,
            center,
            |d| {
                let bend = (d.x * sign).max(-0.2) * 3.2;
                let px = center.x + d.x * 11.5;
                let py = center.y + d.y * 17.5 - bend;
                let pz = center.z + d.z * 21.5;
                let retinotopy = (d.y * 14.0 + d.z * 16.0).sin() * 0.22;
                ([px, py, pz], retinotopy + 0.22)
            },
        ));

        // Paired Antennal Lobes (AL)
        let alc = Vec3::new(sign * 16.0, 25.0, -10.0);
        meshes.push(build_sculpted_sphere(
            "Antennal Lobe",
            AnatomicalLayer::InternalOrgan,
            [244, 152, 76],
            [254, 215, 170],
            0.94,
            false,
            12,
            14,
            alc,
            |d| ((alc + d * 8.2).to_array(), 0.25),
        ));

        // Mushroom Body Pedunculus + L-Lobes (Calyx -> Anterior Vertical/Medial Lobes)
        let mb_spine = [
            ([sign * 22.0, -18.0, 22.0], 6.4, 6.4, 0.3), // Calyx
            ([sign * 20.0, -2.0, 17.0], 3.8, 3.8, 0.2),  // Pedunculus
            ([sign * 18.0, 14.0, 12.0], 5.0, 4.4, 0.3),  // Alpha/Beta junction
            ([sign * 7.0, 15.0, 11.0], 3.4, 3.4, 0.2),   // Medial beta lobe toward midline
        ];
        meshes.push(build_swept_tube(
            "Mushroom Body",
            AnatomicalLayer::InternalOrgan,
            [246, 202, 68],
            [254, 240, 138],
            0.96,
            false,
            &mb_spine,
            10,
        ));
    }

    // Central Complex: Ellipsoid Body Torus + Fan-Shaped Body + Protocerebral Bridge
    meshes.push(build_torus(
        "Ellipsoid Body",
        AnatomicalLayer::InternalOrgan,
        [56, 212, 246],
        [186, 230, 253],
        0.96,
        Vec3::new(0.0, 6.0, 4.0),
        Vec3::X,
        Vec3::Z,
        6.5,
        2.5,
        20,
        10,
    ));
    let fb_center = Vec3::new(0.0, -4.0, 12.0);
    meshes.push(build_sculpted_sphere(
        "Fan-Shaped Body",
        AnatomicalLayer::InternalOrgan,
        [96, 178, 250],
        [191, 219, 254],
        0.92,
        false,
        12,
        16,
        fb_center,
        |d| ((fb_center + Vec3::new(d.x * 10.5, d.y * 6.0, d.z * 7.2)).to_array(), 0.2),
    ));
    let pb_spine = [
        ([-14.0, -20.0, 19.0], 2.2, 2.2, 0.2),
        ([-6.0, -18.5, 21.5], 2.6, 2.5, 0.25),
        ([0.0, -18.0, 22.0], 2.2, 2.2, 0.2),
        ([6.0, -18.5, 21.5], 2.6, 2.5, 0.25),
        ([14.0, -20.0, 19.0], 2.2, 2.2, 0.2),
    ];
    meshes.push(build_swept_tube(
        "Protocerebral Bridge",
        AnatomicalLayer::InternalOrgan,
        [56, 189, 248],
        [186, 230, 253],
        0.94,
        false,
        &pb_spine,
        8,
    ));

    SpeciesAnatomyScene {
        atlas_id: "drosophila_flywire".to_string(),
        meshes,
        landmarks: vec![
            AnatomicalLandmark { label: "Optic Lobe (L)", coord: [-59.0, -5.5, 25.0], color: [224, 215, 255] },
            AnatomicalLandmark { label: "Optic Lobe (R)", coord: [59.0, -5.5, 25.0], color: [224, 215, 255] },
            AnatomicalLandmark { label: "Central Complex", coord: [0.0, 2.0, 21.0], color: [186, 230, 253] },
            AnatomicalLandmark { label: "Mushroom Body", coord: [-22.0, -18.0, 28.0], color: [254, 240, 138] },
            AnatomicalLandmark { label: "Antennal Lobes", coord: [16.0, 28.0, -14.0], color: [254, 215, 170] },
            AnatomicalLandmark { label: "Gnathal Ganglia (GNG)", coord: [0.0, 12.0, -39.0], color: [167, 243, 208] },
        ],
        default_target: [0.0, 1.0, -2.0],
        default_distance: 194.0,
        default_yaw: -0.34,
        default_pitch: 0.32,
        default_slice: [0.0, 2.0, 4.0],
        slice_extent_mm: 78.0,
    }
}

fn celegans_body_rx(y: f32) -> f32 {
    if y > 44.0 {
        let u = ((70.5 - y) / 26.5).clamp(0.0, 1.0);
        3.2 + 14.8 * u.sqrt()
    } else if y < -30.0 {
        let u = ((y + 72.0) / 42.0).clamp(0.0, 1.0);
        (18.0 * u.powf(1.25)).max(1.1)
    } else {
        let mid_t = (y + 30.0) / 74.0;
        18.0 + (mid_t * PI).sin() * 1.4
    }
}

fn build_celegans_wormatlas_scene() -> SpeciesAnatomyScene {
    let mut meshes = Vec::new();

    // 1. Translucent Tapered Nematode Cuticle Body (Anterior Mouth Y=+70.5 -> Slender Tail Whip Y=-72)
    let cuticle_spine: Vec<([f32; 3], f32, f32, f32)> = (0..=42)
        .map(|i| {
            let t = i as f32 / 42.0;
            let y = 70.5 * (1.0 - t) - 72.0 * t;
            let rx = celegans_body_rx(y);
            let annuli = (y * 1.05).sin() * 0.28;
            ([0.0, y, 0.0], rx, rx * 0.86, annuli)
        })
        .collect();
    meshes.push(build_swept_tube(
        "C. elegans Cuticle Body",
        AnatomicalLayer::OuterShell,
        [118, 204, 218],
        [186, 242, 250],
        0.96,
        true,
        &cuticle_spine,
        24,
    ));

    // 2. Two-Bulb Muscular Pharynx seated completely inside anterior head
    let pharynx_spine: Vec<([f32; 3], f32, f32, f32)> = [
        (66.5, 2.0),
        (62.5, 2.8),
        (58.0, 3.8),
        (54.0, 6.0), // Anterior Metacorpus bulb
        (49.0, 3.6),
        (43.5, 2.6), // Narrow pharyngeal isthmus encircled by nerve ring
        (37.5, 6.4),
        (33.5, 7.2), // Terminal grinder bulb
        (30.0, 2.6),
    ]
    .iter()
    .map(|&(y, r)| ([0.0, y, 0.0], r, r * 0.92, 0.25))
    .collect();
    meshes.push(build_swept_tube(
        "Two-Bulb Pharynx",
        AnatomicalLayer::InternalOrgan,
        [236, 120, 144],
        [254, 205, 211],
        0.95,
        false,
        &pharynx_spine,
        14,
    ));

    // 3. Circumpharyngeal Nerve Ring Torus (surrounding isthmus around Y=+41.5)
    meshes.push(build_torus(
        "Nerve Ring",
        AnatomicalLayer::InternalOrgan,
        [250, 204, 78],
        [254, 240, 138],
        0.96,
        Vec3::new(0.0, 41.5, 0.0),
        Vec3::X,
        Vec3::Z,
        9.8,
        2.4,
        24,
        10,
    ));

    // 4. Ventral & Dorsal Nerve Cords curving naturally inside tapered body cuticle
    let vnc_spine: Vec<([f32; 3], f32, f32, f32)> = (0..=16)
        .map(|i| {
            let t = i as f32 / 16.0;
            let y = 41.0 * (1.0 - t) - 52.0 * t;
            let max_z = (celegans_body_rx(y) * 0.86 - 2.6).max(2.0);
            ([0.0, y, -max_z.min(12.0)], 1.9, 1.6, 0.2)
        })
        .collect();
    meshes.push(build_swept_tube(
        "Ventral Nerve Cord",
        AnatomicalLayer::InternalOrgan,
        [56, 189, 248],
        [186, 230, 253],
        0.92,
        false,
        &vnc_spine,
        8,
    ));

    let dnc_spine: Vec<([f32; 3], f32, f32, f32)> = (0..=16)
        .map(|i| {
            let t = i as f32 / 16.0;
            let y = 41.0 * (1.0 - t) - 49.0 * t;
            let max_z = (celegans_body_rx(y) * 0.86 - 2.6).max(2.0);
            ([0.0, y, max_z.min(12.0)], 1.5, 1.3, 0.2)
        })
        .collect();
    meshes.push(build_swept_tube(
        "Dorsal Nerve Cord",
        AnatomicalLayer::InternalOrgan,
        [52, 211, 153],
        [167, 243, 208],
        0.88,
        false,
        &dnc_spine,
        8,
    ));

    // 5. Midbody Intestinal Tube + Posterior Pre-Anal / Lumbar Ganglia
    let gut_spine: Vec<([f32; 3], f32, f32, f32)> = (0..=12)
        .map(|i| {
            let t = i as f32 / 12.0;
            let y = 29.0 * (1.0 - t) - 44.0 * t;
            let r = 5.8 * (1.0 - t * 0.28);
            ([0.0, y, 0.8], r, r * 0.85, 0.05)
        })
        .collect();
    meshes.push(build_swept_tube(
        "Intestinal Lumen",
        AnatomicalLayer::InternalOrgan,
        [134, 158, 178],
        [186, 210, 228],
        0.44,
        false,
        &gut_spine,
        12,
    ));

    let lumbar = Vec3::new(0.0, -49.0, -1.5);
    meshes.push(build_sculpted_sphere(
        "Lumbar / Tail Ganglia",
        AnatomicalLayer::InternalOrgan,
        [246, 146, 86],
        [254, 215, 170],
        0.88,
        false,
        10,
        14,
        lumbar,
        |d| ((lumbar + Vec3::new(d.x * 7.2, d.y * 5.8, d.z * 4.2)).to_array(), 0.2),
    ));

    SpeciesAnatomyScene {
        atlas_id: "celegans_wormatlas".to_string(),
        meshes,
        landmarks: vec![
            AnatomicalLandmark { label: "Mouth / Amphids", coord: [0.0, 69.5, 6.0], color: [254, 215, 170] },
            AnatomicalLandmark { label: "2-Bulb Pharynx", coord: [0.0, 54.0, 16.0], color: [254, 205, 211] },
            AnatomicalLandmark { label: "Nerve Ring", coord: [0.0, 41.5, -16.5], color: [254, 240, 138] },
            AnatomicalLandmark { label: "Dorsal Nerve Cord", coord: [0.0, 0.0, 15.0], color: [167, 243, 208] },
            AnatomicalLandmark { label: "Ventral Nerve Cord", coord: [0.0, -6.0, -16.0], color: [186, 230, 253] },
            AnatomicalLandmark { label: "Lumbar Ganglia", coord: [0.0, -49.0, -11.0], color: [254, 215, 170] },
            AnatomicalLandmark { label: "Tail Whip", coord: [0.0, -70.0, 3.0], color: [186, 242, 250] },
        ],
        default_target: [0.0, 1.0, 0.0],
        default_distance: 162.0,
        default_yaw: -1.12,
        default_pitch: 0.30,
        default_slice: [0.0, 28.0, 0.0],
        slice_extent_mm: 78.0,
    }
}

pub fn build_species_anatomy(atlas_id: &str) -> SpeciesAnatomyScene {
    match atlas_id {
        "mouse_ccfv3" => build_mouse_ccfv3_scene(),
        "drosophila_flywire" => build_drosophila_flywire_scene(),
        "celegans_wormatlas" => build_celegans_wormatlas_scene(),
        _ => build_human_mni152_scene(),
    }
}

fn ellipse_poly(cx: f32, cy: f32, rx: f32, ry: f32, n: usize, gyri_amp: f32) -> Vec<[f32; 2]> {
    (0..n)
        .map(|i| {
            let a = (i as f32 / n as f32) * TAU;
            let mod_r = 1.0 + gyri_amp * ((a * 7.0).sin() * 0.6 + (a * 11.0).cos() * 0.4);
            [cx + a.cos() * rx * mod_r, cy + a.sin() * ry * mod_r]
        })
        .collect()
}

/// Generates species-accurate 2D cross-sectional anatomical contours in slice-plane coordinates `(u, v)`.
pub fn generate_2d_slice_contours(
    atlas_id: &str,
    plane_mode: u8,
    slice_val: f32,
) -> Vec<SliceContour> {
    let mut contours = Vec::new();
    match atlas_id {
        "mouse_ccfv3" => match plane_mode {
            // 0 = Axial (u=X, v=Y, ortho=Z)
            0 => {
                let s = (1.0 - ((slice_val - 2.0) / 24.0).powi(2)).max(0.25).sqrt();
                contours.push(SliceContour {
                    points: ellipse_poly(0.0, -2.0, 42.0 * s, 40.0 * s, 40, 0.015),
                    fill_rgba: [30, 45, 68, 175],
                    stroke_rgba: [125, 211, 252, 165],
                    stroke_width: 1.5,
                    closed: true,
                });
                // Paired Olfactory bulbs anteriorly
                for &sign in &[-1.0_f32, 1.0_f32] {
                    contours.push(SliceContour {
                        points: ellipse_poly(sign * 11.5, 45.0, 8.5 * s, 12.0 * s, 20, 0.0),
                        fill_rgba: [78, 46, 34, 150],
                        stroke_rgba: [251, 146, 60, 155],
                        stroke_width: 1.2,
                        closed: true,
                    });
                }
                // Posterior cerebellum
                contours.push(SliceContour {
                    points: ellipse_poly(0.0, -49.0, 30.0 * s, 12.0 * s, 28, 0.03),
                    fill_rgba: [52, 42, 82, 165],
                    stroke_rgba: [192, 132, 252, 155],
                    stroke_width: 1.3,
                    closed: true,
                });
            }
            // 1 = Coronal (u=X, v=Z, ortho=Y)
            1 => {
                let is_ob = slice_val > 34.0;
                let is_cb = slice_val < -37.0;
                if is_ob {
                    for &sign in &[-1.0_f32, 1.0_f32] {
                        contours.push(SliceContour {
                            points: ellipse_poly(sign * 11.5, -3.5, 8.5, 7.5, 20, 0.0),
                            fill_rgba: [78, 46, 34, 165],
                            stroke_rgba: [251, 146, 60, 170],
                            stroke_width: 1.4,
                            closed: true,
                        });
                    }
                } else if is_cb {
                    contours.push(SliceContour {
                        points: ellipse_poly(0.0, 0.0, 31.0, 15.0, 28, 0.04),
                        fill_rgba: [52, 42, 82, 175],
                        stroke_rgba: [192, 132, 252, 165],
                        stroke_width: 1.4,
                        closed: true,
                    });
                } else {
                    for &sign in &[-1.0_f32, 1.0_f32] {
                        contours.push(SliceContour {
                            points: ellipse_poly(sign * 19.5, 2.0, 19.0, 20.5, 28, 0.01),
                            fill_rgba: [30, 45, 68, 175],
                            stroke_rgba: [125, 211, 252, 165],
                            stroke_width: 1.4,
                            closed: true,
                        });
                    }
                }
            }
            // 2 = Sagittal (u=Y, v=Z, ortho=X)
            _ => {
                let s = (1.0 - (slice_val / 44.0).powi(2)).max(0.25).sqrt();
                contours.push(SliceContour {
                    points: ellipse_poly(-2.0, 2.0, 38.0 * s, 21.0 * s, 32, 0.01),
                    fill_rgba: [30, 45, 68, 175],
                    stroke_rgba: [125, 211, 252, 165],
                    stroke_width: 1.5,
                    closed: true,
                });
                contours.push(SliceContour {
                    points: ellipse_poly(45.5, -3.5, 12.5 * s, 7.8 * s, 20, 0.0),
                    fill_rgba: [78, 46, 34, 160],
                    stroke_rgba: [251, 146, 60, 160],
                    stroke_width: 1.2,
                    closed: true,
                });
                contours.push(SliceContour {
                    points: ellipse_poly(-49.0, 0.5, 13.0 * s, 14.5 * s, 24, 0.04),
                    fill_rgba: [52, 42, 82, 170],
                    stroke_rgba: [192, 132, 252, 160],
                    stroke_width: 1.3,
                    closed: true,
                });
            }
        },
        "drosophila_flywire" => match plane_mode {
            // Axial (X vs Y)
            0 => {
                contours.push(SliceContour {
                    points: ellipse_poly(0.0, 2.0, 37.0, 23.0, 32, 0.02),
                    fill_rgba: [30, 45, 70, 175],
                    stroke_rgba: [125, 211, 252, 165],
                    stroke_width: 1.5,
                    closed: true,
                });
                for &sign in &[-1.0_f32, 1.0_f32] {
                    contours.push(SliceContour {
                        points: ellipse_poly(sign * 56.0, -6.0, 16.0, 19.0, 24, 0.02),
                        fill_rgba: [48, 38, 86, 170],
                        stroke_rgba: [167, 139, 250, 170],
                        stroke_width: 1.4,
                        closed: true,
                    });
                }
                contours.push(SliceContour {
                    points: ellipse_poly(0.0, 6.0, 7.0, 5.5, 18, 0.0),
                    fill_rgba: [14, 116, 144, 150],
                    stroke_rgba: [56, 189, 248, 175],
                    stroke_width: 1.2,
                    closed: true,
                });
            }
            // Coronal (X vs Z, front view showing iconic fly brain + wing optic lobes + GNG!)
            1 => {
                contours.push(SliceContour {
                    points: ellipse_poly(0.0, 4.0, 37.0, 24.0, 32, 0.02),
                    fill_rgba: [30, 45, 70, 175],
                    stroke_rgba: [125, 211, 252, 165],
                    stroke_width: 1.5,
                    closed: true,
                });
                contours.push(SliceContour {
                    points: ellipse_poly(0.0, -27.0, 18.0, 11.5, 20, 0.0),
                    fill_rgba: [24, 62, 58, 165],
                    stroke_rgba: [52, 211, 153, 160],
                    stroke_width: 1.3,
                    closed: true,
                });
                for &sign in &[-1.0_f32, 1.0_f32] {
                    contours.push(SliceContour {
                        points: ellipse_poly(sign * 56.0, 4.0, 15.5, 22.0, 24, 0.02),
                        fill_rgba: [48, 38, 86, 170],
                        stroke_rgba: [167, 139, 250, 170],
                        stroke_width: 1.4,
                        closed: true,
                    });
                }
                contours.push(SliceContour {
                    points: ellipse_poly(0.0, 4.0, 7.0, 7.0, 18, 0.0),
                    fill_rgba: [14, 116, 144, 145],
                    stroke_rgba: [56, 189, 248, 180],
                    stroke_width: 1.3,
                    closed: true,
                });
            }
            // Sagittal (Y vs Z)
            _ => {
                contours.push(SliceContour {
                    points: ellipse_poly(2.0, 4.0, 23.0, 24.0, 28, 0.02),
                    fill_rgba: [30, 45, 70, 175],
                    stroke_rgba: [125, 211, 252, 165],
                    stroke_width: 1.5,
                    closed: true,
                });
                contours.push(SliceContour {
                    points: ellipse_poly(11.0, -27.0, 13.0, 11.0, 20, 0.0),
                    fill_rgba: [24, 62, 58, 165],
                    stroke_rgba: [52, 211, 153, 160],
                    stroke_width: 1.3,
                    closed: true,
                });
            }
        },
        "celegans_wormatlas" => match plane_mode {
            // Axial (X vs Y — top-down elongated worm!)
            0 => {
                contours.push(SliceContour {
                    points: (0..=44)
                        .map(|i| {
                            let t = i as f32 / 44.0;
                            let a = t * TAU;
                            let y = a.sin() * 70.0;
                            let r = if y > 46.0 {
                                18.0 * ((70.0 - y) / 24.0).max(0.08).sqrt()
                            } else if y < -32.0 {
                                18.0 * ((y + 72.0) / 40.0).max(0.05).powf(1.25)
                            } else {
                                18.2
                            };
                            [a.cos().signum() * r * (a.cos().abs()), y]
                        })
                        .collect(),
                    fill_rgba: [24, 52, 64, 175],
                    stroke_rgba: [103, 232, 249, 175],
                    stroke_width: 1.5,
                    closed: true,
                });
                // Pharynx inner contour
                contours.push(SliceContour {
                    points: ellipse_poly(0.0, 49.0, 6.5, 16.0, 20, 0.0),
                    fill_rgba: [82, 34, 48, 160],
                    stroke_rgba: [251, 113, 133, 175],
                    stroke_width: 1.2,
                    closed: true,
                });
                // Nerve ring bar across isthmus
                contours.push(SliceContour {
                    points: vec![[-11.0, 41.5], [11.0, 41.5]],
                    fill_rgba: [0, 0, 0, 0],
                    stroke_rgba: [250, 204, 21, 195],
                    stroke_width: 2.2,
                    closed: false,
                });
            }
            // Coronal cross-section at AP coordinate Y (X vs Z circular cross section!)
            1 => {
                let y = slice_val;
                let r = if y > 46.0 {
                    18.0 * ((70.0 - y) / 24.0).max(0.08).sqrt()
                } else if y < -32.0 {
                    18.0 * ((y + 72.0) / 40.0).max(0.05).powf(1.25)
                } else {
                    18.2
                };
                contours.push(SliceContour {
                    points: ellipse_poly(0.0, 0.0, r, r * 0.86, 28, 0.0),
                    fill_rgba: [24, 52, 64, 175],
                    stroke_rgba: [103, 232, 249, 175],
                    stroke_width: 1.5,
                    closed: true,
                });
                if (32.0..=67.0).contains(&y) {
                    let pr = if (y - 54.5).abs() < 6.0 || (y - 34.5).abs() < 6.0 { 6.5 } else { 3.2 };
                    contours.push(SliceContour {
                        points: ellipse_poly(0.0, 0.0, pr, pr, 18, 0.0),
                        fill_rgba: [82, 34, 48, 170],
                        stroke_rgba: [251, 113, 133, 185],
                        stroke_width: 1.2,
                        closed: true,
                    });
                }
                if (38.0..=45.0).contains(&y) {
                    contours.push(SliceContour {
                        points: ellipse_poly(0.0, 0.0, 11.0, 10.0, 22, 0.0),
                        fill_rgba: [0, 0, 0, 0],
                        stroke_rgba: [250, 204, 21, 195],
                        stroke_width: 1.8,
                        closed: true,
                    });
                }
            }
            // Sagittal (Y vs Z — full worm side profile showing pharynx, nerve ring & cords!)
            _ => {
                contours.push(SliceContour {
                    points: (0..=44)
                        .map(|i| {
                            let t = i as f32 / 44.0;
                            let a = t * TAU;
                            let y = a.sin() * 70.0;
                            let rz = if y > 46.0 {
                                15.5 * ((70.0 - y) / 24.0).max(0.08).sqrt()
                            } else if y < -32.0 {
                                15.5 * ((y + 72.0) / 40.0).max(0.05).powf(1.25)
                            } else {
                                15.6
                            };
                            [y, a.cos().signum() * rz * (a.cos().abs())]
                        })
                        .collect(),
                    fill_rgba: [24, 52, 64, 175],
                    stroke_rgba: [103, 232, 249, 175],
                    stroke_width: 1.5,
                    closed: true,
                });
                contours.push(SliceContour {
                    points: ellipse_poly(49.0, 0.0, 16.0, 6.2, 20, 0.0),
                    fill_rgba: [82, 34, 48, 160],
                    stroke_rgba: [251, 113, 133, 175],
                    stroke_width: 1.2,
                    closed: true,
                });
                // Ventral & Dorsal cord traces
                contours.push(SliceContour {
                    points: vec![[-54.0, -12.0], [41.0, -12.0]],
                    fill_rgba: [0, 0, 0, 0],
                    stroke_rgba: [56, 189, 248, 180],
                    stroke_width: 1.6,
                    closed: false,
                });
                contours.push(SliceContour {
                    points: vec![[-50.0, 12.0], [41.0, 12.0]],
                    fill_rgba: [0, 0, 0, 0],
                    stroke_rgba: [52, 211, 153, 170],
                    stroke_width: 1.4,
                    closed: false,
                });
            }
        },
        // human_mni152
        _ => match plane_mode {
            // 0 = Axial (u=X, v=Y, ortho=Z)
            0 => {
                let s = (1.0 - ((slice_val - 12.0) / 62.0).powi(2)).max(0.22).sqrt();
                for &sign in &[-1.0_f32, 1.0_f32] {
                    contours.push(SliceContour {
                        points: ellipse_poly(sign * 33.5 * s, -15.0, 32.5 * s, 77.0 * s, 36, 0.035),
                        fill_rgba: [30, 44, 66, 180],
                        stroke_rgba: [125, 211, 252, 165],
                        stroke_width: 1.5,
                        closed: true,
                    });
                    // Internal subcortical oval / ventricles when Z is mid-level
                    if (-15.0..=32.0).contains(&slice_val) {
                        contours.push(SliceContour {
                            points: ellipse_poly(sign * 13.5, -14.0, 8.5, 19.0, 18, 0.0),
                            fill_rgba: [18, 28, 44, 175],
                            stroke_rgba: [148, 163, 184, 120],
                            stroke_width: 1.0,
                            closed: true,
                        });
                    }
                }
                if slice_val < -8.0 {
                    contours.push(SliceContour {
                        points: ellipse_poly(0.0, -56.0, 42.0, 19.0, 26, 0.04),
                        fill_rgba: [52, 42, 82, 170],
                        stroke_rgba: [192, 132, 252, 160],
                        stroke_width: 1.3,
                        closed: true,
                    });
                }
            }
            // 1 = Coronal (u=X, v=Z, ortho=Y)
            1 => {
                let s = (1.0 - ((slice_val + 14.0) / 88.0).powi(2)).max(0.22).sqrt();
                for &sign in &[-1.0_f32, 1.0_f32] {
                    contours.push(SliceContour {
                        points: ellipse_poly(sign * 33.5 * s, 14.0, 32.5 * s, 52.0 * s, 32, 0.035),
                        fill_rgba: [30, 44, 66, 180],
                        stroke_rgba: [125, 211, 252, 165],
                        stroke_width: 1.5,
                        closed: true,
                    });
                }
                if (-45.0..=36.0).contains(&slice_val) {
                    contours.push(SliceContour {
                        points: ellipse_poly(0.0, 22.0, 15.0 * s, 5.5, 16, 0.0),
                        fill_rgba: [78, 66, 36, 165],
                        stroke_rgba: [253, 230, 138, 175],
                        stroke_width: 1.2,
                        closed: true,
                    });
                }
                if slice_val < -35.0 {
                    contours.push(SliceContour {
                        points: ellipse_poly(0.0, -31.0, 42.0, 16.0, 24, 0.04),
                        fill_rgba: [52, 42, 82, 170],
                        stroke_rgba: [192, 132, 252, 160],
                        stroke_width: 1.3,
                        closed: true,
                    });
                }
            }
            // 2 = Sagittal (u=Y, v=Z, ortho=X)
            _ => {
                let s = (1.0 - (slice_val / 72.0).powi(2)).max(0.25).sqrt();
                contours.push(SliceContour {
                    points: ellipse_poly(-14.0, 14.0, 81.0 * s, 53.0 * s, 40, 0.035),
                    fill_rgba: [30, 44, 66, 180],
                    stroke_rgba: [125, 211, 252, 165],
                    stroke_width: 1.5,
                    closed: true,
                });
                if slice_val.abs() < 22.0 {
                    // Corpus callosum arch
                    let cc_pts: Vec<[f32; 2]> = (0..=16)
                        .map(|i| {
                            let t = i as f32 / 16.0;
                            let y = 37.0 * (1.0 - t) - 44.0 * t;
                            let z = 3.0 + (t * PI).sin() * 23.0 + t * 11.0;
                            [y, z]
                        })
                        .collect();
                    contours.push(SliceContour {
                        points: cc_pts,
                        fill_rgba: [0, 0, 0, 0],
                        stroke_rgba: [253, 230, 138, 195],
                        stroke_width: 2.5,
                        closed: false,
                    });
                }
                // Cerebellum + Brainstem
                contours.push(SliceContour {
                    points: ellipse_poly(-57.0, -30.5, 19.0, 15.5, 24, 0.04),
                    fill_rgba: [52, 42, 82, 175],
                    stroke_rgba: [192, 132, 252, 165],
                    stroke_width: 1.3,
                    closed: true,
                });
                contours.push(SliceContour {
                    points: vec![[-15.0, -4.0], [-22.0, -24.0], [-28.0, -42.0], [-32.0, -56.0]],
                    fill_rgba: [0, 0, 0, 0],
                    stroke_rgba: [148, 163, 184, 175],
                    stroke_width: 3.0,
                    closed: false,
                });
            }
        },
    }
    contours
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn all_four_species_anatomy_scenes_build_valid_geometry_and_slices() {
        for pid in [
            "human_mni152",
            "mouse_ccfv3",
            "drosophila_flywire",
            "celegans_wormatlas",
        ] {
            let scene = build_species_anatomy(pid);
            assert_eq!(scene.atlas_id, pid);
            assert!(
                scene.meshes.len() >= 5,
                "{pid} expected rich multi-part anatomy"
            );
            assert!(!scene.landmarks.is_empty(), "{pid} landmarks missing");
            for m in &scene.meshes {
                assert!(!m.vertices.is_empty(), "{} empty vertices", m.name);
                assert!(!m.triangles.is_empty(), "{} empty triangles", m.name);
                assert_eq!(m.vertices.len(), m.normals.len());
                assert_eq!(m.vertices.len(), m.crease.len());
                for tri in &m.triangles {
                    for &idx in tri {
                        assert!((idx as usize) < m.vertices.len());
                    }
                }
                for &n in &m.normals {
                    let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
                    assert!((len - 1.0).abs() < 0.05, "{} non-unit normal", m.name);
                }
            }
            for mode in 0..=2_u8 {
                let slice_v = match mode {
                    0 => scene.default_slice[2],
                    1 => scene.default_slice[1],
                    _ => scene.default_slice[0],
                };
                let contours = generate_2d_slice_contours(pid, mode, slice_v);
                assert!(!contours.is_empty(), "{pid} plane {mode} contours empty");
            }
        }
    }
}
