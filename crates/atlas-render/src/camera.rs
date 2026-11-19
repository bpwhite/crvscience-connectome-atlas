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

use glam::{Mat4, Vec2, Vec3};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ViewPreset {
    Perspective3D,
    LeftSagittal,
    SuperiorAxial,
    AnteriorCoronal,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OrbitCamera {
    pub yaw_rad: f32,
    pub pitch_rad: f32,
    pub distance: f32,
    pub target: Vec3,
    pub fov_y_rad: f32,
}

impl Default for OrbitCamera {
    fn default() -> Self {
        Self {
            yaw_rad: -0.68,
            pitch_rad: 0.34,
            distance: 215.0,
            target: Vec3::new(0.0, -8.0, 10.0),
            fov_y_rad: 42.0_f32.to_radians(),
        }
    }
}

impl OrbitCamera {
    pub fn apply_preset(&mut self, preset: ViewPreset) {
        match preset {
            ViewPreset::Perspective3D => {
                self.yaw_rad = -0.68;
                self.pitch_rad = 0.34;
            }
            ViewPreset::LeftSagittal => {
                self.yaw_rad = -std::f32::consts::FRAC_PI_2;
                self.pitch_rad = 0.0;
            }
            ViewPreset::SuperiorAxial => {
                self.yaw_rad = 0.0;
                self.pitch_rad = 1.48;
            }
            ViewPreset::AnteriorCoronal => {
                self.yaw_rad = 0.0;
                self.pitch_rad = 0.0;
            }
        }
    }

    pub fn orbit(&mut self, delta_px: Vec2) {
        self.yaw_rad -= delta_px.x * 0.0075;
        self.pitch_rad = (self.pitch_rad + delta_px.y * 0.0075).clamp(-1.50, 1.50);
    }

    pub fn zoom(&mut self, scroll_delta_y: f32) {
        let factor = (1.0 - scroll_delta_y * 0.0015).clamp(0.82, 1.22);
        self.distance = (self.distance * factor).clamp(65.0, 520.0);
    }

    pub fn eye_position(&self) -> Vec3 {
        let cos_p = self.pitch_rad.cos();
        // Stereotaxic convention: +X Right, +Y Anterior, +Z Superior (Z-up)
        let offset = Vec3::new(
            self.distance * cos_p * self.yaw_rad.sin(),
            self.distance * cos_p * self.yaw_rad.cos(),
            self.distance * self.pitch_rad.sin(),
        );
        self.target + offset
    }

    pub fn view_proj_matrix(&self, aspect: f32) -> Mat4 {
        let view = Mat4::look_at_rh(self.eye_position(), self.target, Vec3::Z);
        let proj = Mat4::perspective_rh(self.fov_y_rad, aspect.max(0.2), 2.0, 1500.0);
        proj * view
    }

    /// Projects a 3D stereotaxic coordinate `[x, y, z]` to viewport pixel coordinates `(px, py, ndc_depth)`.
    pub fn project_point(&self, pt: [f32; 3], viewport_wh: [f32; 2]) -> Option<([f32; 2], f32)> {
        let vp = self.view_proj_matrix(viewport_wh[0] / viewport_wh[1].max(1.0));
        let clip = vp * Vec3::from_array(pt).extend(1.0);
        if clip.w <= 0.1 {
            return None;
        }
        let ndc = clip.truncate() / clip.w;
        if ndc.z < 0.0 || ndc.z > 1.0 {
            return None;
        }
        let sx = (ndc.x * 0.5 + 0.5) * viewport_wh[0];
        let sy = (0.5 - ndc.y * 0.5) * viewport_wh[1];
        Some(([sx, sy], ndc.z))
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FlyThroughController {
    pub active: bool,
    pub progress_u: f32,
    pub speed: f32,
}

impl Default for FlyThroughController {
    fn default() -> Self {
        Self {
            active: false,
            progress_u: 0.0,
            speed: 0.14,
        }
    }
}

impl FlyThroughController {
    /// Samples a smooth chase-camera `(eye_position, look_target)` along a polyline signaling pathway.
    pub fn sample_pathway_pose(&self, polyline: &[[f32; 3]]) -> Option<(Vec3, Vec3)> {
        if polyline.len() < 2 {
            return None;
        }
        let u = self.progress_u.clamp(0.0, 0.9999);
        let scaled = u * ((polyline.len() - 1) as f32);
        let idx = (scaled.floor() as usize).min(polyline.len() - 2);
        let frac = scaled - (idx as f32);
        let p0 = Vec3::from_array(polyline[idx]);
        let p1 = Vec3::from_array(polyline[idx + 1]);
        let pt = p0.lerp(p1, frac);

        let look_u = (u + 0.08).min(1.0) * ((polyline.len() - 1) as f32);
        let l_idx = (look_u.floor() as usize).min(polyline.len() - 2);
        let l_frac = look_u - (l_idx as f32);
        let look_pt = Vec3::from_array(polyline[l_idx]).lerp(Vec3::from_array(polyline[l_idx + 1]), l_frac);

        let tang = (p1 - p0).normalize_or_zero();
        let radial_up = Vec3::new(pt.x * 0.25, pt.y * 0.25, 18.0).normalize_or_zero();
        let eye = pt - tang * 22.0 + radial_up * 10.0;
        Some((eye, look_pt))
    }
}
