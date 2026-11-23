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

#[derive(Debug, Clone, Copy, PartialEq)]
struct CameraTransition {
    target_yaw_rad: f32,
    target_pitch_rad: f32,
    target_distance: Option<f32>,
    target_focus: Option<Vec3>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OrbitCamera {
    pub yaw_rad: f32,
    pub pitch_rad: f32,
    pub distance: f32,
    pub target: Vec3,
    pub fov_y_rad: f32,
    #[serde(skip)]
    transition: Option<CameraTransition>,
}

impl Default for OrbitCamera {
    fn default() -> Self {
        Self {
            yaw_rad: -0.68,
            pitch_rad: 0.34,
            distance: 215.0,
            target: Vec3::new(0.0, -8.0, 10.0),
            fov_y_rad: 42.0_f32.to_radians(),
            transition: None,
        }
    }
}

impl OrbitCamera {
    pub fn preset_angles(preset: ViewPreset) -> (f32, f32) {
        match preset {
            ViewPreset::Perspective3D => (-0.68, 0.34),
            ViewPreset::LeftSagittal => (-std::f32::consts::FRAC_PI_2, 0.0),
            ViewPreset::SuperiorAxial => (0.0, 1.48),
            ViewPreset::AnteriorCoronal => (0.0, 0.0),
        }
    }

    pub fn apply_preset(&mut self, preset: ViewPreset) {
        let (yaw, pitch) = Self::preset_angles(preset);
        self.yaw_rad = yaw;
        self.pitch_rad = pitch;
        self.transition = None;
    }

    pub fn start_preset_transition(&mut self, preset: ViewPreset) {
        let (yaw, pitch) = Self::preset_angles(preset);
        self.start_transition(yaw, pitch, None, None);
    }

    pub fn start_transition(
        &mut self,
        target_yaw_rad: f32,
        target_pitch_rad: f32,
        target_distance: Option<f32>,
        target_focus: Option<Vec3>,
    ) {
        self.transition = Some(CameraTransition {
            target_yaw_rad,
            target_pitch_rad: target_pitch_rad.clamp(-1.50, 1.50),
            target_distance: target_distance.map(|d| d.clamp(65.0, 520.0)),
            target_focus,
        });
    }

    pub fn cancel_transition(&mut self) {
        self.transition = None;
    }

    pub fn finish_transition(&mut self) {
        if let Some(tr) = self.transition.take() {
            self.yaw_rad = tr.target_yaw_rad;
            self.pitch_rad = tr.target_pitch_rad;
            if let Some(d) = tr.target_distance {
                self.distance = d;
            }
            if let Some(f) = tr.target_focus {
                self.target = f;
            }
        }
    }

    /// Advances smooth exponential camera preset interpolation; returns `true` while active.
    pub fn step_transition(&mut self, dt: f32) -> bool {
        let Some(tr) = self.transition else {
            return false;
        };
        let alpha = (1.0 - (-11.0 * dt.clamp(0.0, 0.1)).exp()).clamp(0.05, 1.0);
        let mut d_yaw = (tr.target_yaw_rad - self.yaw_rad).rem_euclid(std::f32::consts::TAU);
        if d_yaw > std::f32::consts::PI {
            d_yaw -= std::f32::consts::TAU;
        }
        let d_pitch = tr.target_pitch_rad - self.pitch_rad;
        self.yaw_rad += d_yaw * alpha;
        self.pitch_rad += d_pitch * alpha;

        let mut max_err = d_yaw.abs().max(d_pitch.abs());
        if let Some(td) = tr.target_distance {
            let dd = td - self.distance;
            self.distance += dd * alpha;
            max_err = max_err.max(dd.abs() / 100.0);
        }
        if let Some(tf) = tr.target_focus {
            let df = tf - self.target;
            self.target += df * alpha;
            max_err = max_err.max(df.length() / 50.0);
        }

        if max_err < 1e-3 {
            self.finish_transition();
            false
        } else {
            true
        }
    }

    pub fn orbit(&mut self, delta_px: Vec2) {
        self.transition = None;
        self.yaw_rad -= delta_px.x * 0.0075;
        self.pitch_rad = (self.pitch_rad + delta_px.y * 0.0075).clamp(-1.50, 1.50);
    }

    /// Pans the camera focus point in the view plane (right-drag interaction).
    pub fn pan(&mut self, delta_px: Vec2, viewport_height: f32) {
        self.transition = None;
        let forward = (self.target - self.eye_position()).normalize_or_zero();
        let mut right = forward.cross(Vec3::Z);
        if right.length_squared() < 1e-6 {
            right = Vec3::new(self.yaw_rad.cos(), -self.yaw_rad.sin(), 0.0);
        }
        right = right.normalize_or_zero();
        let up = right.cross(forward).normalize_or_zero();
        let world_per_px =
            (2.0 * self.distance * (self.fov_y_rad * 0.5).tan()) / viewport_height.max(120.0);
        self.target += (-right * delta_px.x + up * delta_px.y) * world_per_px;
        self.target = self.target.clamp(Vec3::splat(-260.0), Vec3::splat(260.0));
    }

    pub fn zoom(&mut self, scroll_delta_y: f32) {
        self.transition = None;
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

    /// Projects a 3D stereotaxic coordinate `[x, y, z]` using a precomputed `view_proj` matrix.
    pub fn project_point_with_vp(
        vp: Mat4,
        pt: [f32; 3],
        viewport_wh: [f32; 2],
    ) -> Option<([f32; 2], f32)> {
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

    /// Projects a 3D stereotaxic coordinate `[x, y, z]` to viewport pixel coordinates `(px, py, ndc_depth)`.
    pub fn project_point(&self, pt: [f32; 3], viewport_wh: [f32; 2]) -> Option<([f32; 2], f32)> {
        let vp = self.view_proj_matrix(viewport_wh[0] / viewport_wh[1].max(1.0));
        Self::project_point_with_vp(vp, pt, viewport_wh)
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn camera_pan_precomputed_vp_and_smooth_preset_interpolation_work() {
        let mut cam = OrbitCamera::default();
        let vp = cam.view_proj_matrix(1600.0 / 1000.0);
        let p_direct = cam.project_point([0.0, 0.0, 0.0], [1600.0, 1000.0]).unwrap();
        let p_vp = OrbitCamera::project_point_with_vp(vp, [0.0, 0.0, 0.0], [1600.0, 1000.0]).unwrap();
        assert!((p_direct.0[0] - p_vp.0[0]).abs() < 1e-4);
        assert!((p_direct.0[1] - p_vp.0[1]).abs() < 1e-4);

        let initial_target = cam.target;
        cam.pan(Vec2::new(40.0, -25.0), 800.0);
        assert!(cam.target.distance(initial_target) > 1.0);

        cam.start_preset_transition(ViewPreset::SuperiorAxial);
        let mut active = true;
        for _ in 0..80 {
            active = cam.step_transition(1.0 / 60.0);
            if !active {
                break;
            }
        }
        assert!(!active);
        let (ty, tp) = OrbitCamera::preset_angles(ViewPreset::SuperiorAxial);
        assert!((cam.yaw_rad - ty).abs() < 1e-3);
        assert!((cam.pitch_rad - tp).abs() < 1e-3);
    }
}
