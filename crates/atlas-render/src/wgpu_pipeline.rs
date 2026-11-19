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

use crate::surfaces_and_glyphs::SurgicalCutawayMode;
use glam::{Mat4, Vec3};

pub const WGSL_CONNECTOME_SHADER: &str = include_str!("../shaders/connectome.wgsl");

/// Packs the 96-byte `ViewUniforms` buffer for `connectome.wgsl` using 100% safe byte encoding (`#![forbid(unsafe_code)]`).
pub fn encode_view_uniform_bytes(
    view_proj: Mat4,
    camera_eye: Vec3,
    tonic_phase: f32,
    shell_opacity: f32,
    cutaway_offset: f32,
    cutaway_mode: SurgicalCutawayMode,
) -> [u8; 96] {
    let mut out = [0_u8; 96];
    for (i, val) in view_proj.to_cols_array().iter().enumerate() {
        out[i * 4..i * 4 + 4].copy_from_slice(&val.to_le_bytes());
    }
    out[64..68].copy_from_slice(&camera_eye.x.to_le_bytes());
    out[68..72].copy_from_slice(&camera_eye.y.to_le_bytes());
    out[72..76].copy_from_slice(&camera_eye.z.to_le_bytes());
    out[76..80].copy_from_slice(&1.0_f32.to_le_bytes());
    out[80..84].copy_from_slice(&tonic_phase.to_le_bytes());
    out[84..88].copy_from_slice(&shell_opacity.to_le_bytes());
    out[88..92].copy_from_slice(&cutaway_offset.to_le_bytes());
    out[92..96].copy_from_slice(&cutaway_mode.shader_code().to_le_bytes());
    out
}

/// Validates shader module compilation descriptor readiness for `wgpu 24.0`.
pub fn shader_descriptor() -> wgpu::ShaderModuleDescriptor<'static> {
    wgpu::ShaderModuleDescriptor {
        label: Some("connectome_atlas_v050_wgsl"),
        source: wgpu::ShaderSource::Wgsl(std::borrow::Cow::Borrowed(WGSL_CONNECTOME_SHADER)),
    }
}
