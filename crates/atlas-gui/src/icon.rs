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

//! Embedded procedural RGBA window icon (`ViewportBuilder::with_icon`) matching
//! `crates/atlas-gui/assets/connectome_atlas.ico`.

use eframe::egui;

fn blend_px(buf: &mut [u8], idx: usize, r: f32, g: f32, b: f32, a: f32) {
    let alpha = a.clamp(0.0, 1.0);
    if alpha <= 0.001 {
        return;
    }
    let da = f32::from(buf[idx + 3]) / 255.0;
    let out_a = alpha + da * (1.0 - alpha);
    if out_a <= 1e-5 {
        return;
    }
    let br = f32::from(buf[idx]);
    let bg = f32::from(buf[idx + 1]);
    let bb = f32::from(buf[idx + 2]);
    buf[idx] = ((r * alpha + br * da * (1.0 - alpha)) / out_a).round().clamp(0.0, 255.0) as u8;
    buf[idx + 1] = ((g * alpha + bg * da * (1.0 - alpha)) / out_a).round().clamp(0.0, 255.0) as u8;
    buf[idx + 2] = ((b * alpha + bb * da * (1.0 - alpha)) / out_a).round().clamp(0.0, 255.0) as u8;
    buf[idx + 3] = (out_a * 255.0).round().clamp(0.0, 255.0) as u8;
}

fn sd_rounded_box(px: f32, py: f32, bx: f32, by: f32, r: f32) -> f32 {
    let qx = px.abs() - bx + r;
    let qy = py.abs() - by + r;
    qx.max(0.0).hypot(qy.max(0.0)) + qx.max(qy).min(0.0) - r
}

fn dist_seg(px: f32, py: f32, ax: f32, ay: f32, bx: f32, by: f32) -> f32 {
    let vx = bx - ax;
    let vy = by - ay;
    let wx = px - ax;
    let wy = py - ay;
    let len_sq = vx * vx + vy * vy;
    let t = if len_sq > 1e-6 {
        ((wx * vx + wy * vy) / len_sq).clamp(0.0, 1.0)
    } else {
        0.0
    };
    (px - (ax + vx * t)).hypot(py - (ay + vy * t))
}

/// Renders an anti-aliased `size x size` RGBA ConnectomeAtlas badge matching `assets/connectome_atlas.ico`.
pub fn render_connectome_atlas_rgba(size: usize) -> Vec<u8> {
    let mut buf = vec![0u8; size * size * 4];
    let aa = 2.2 / (size as f32).max(1.0);
    let tracts: [([f32; 2], [f32; 2], [f32; 3]); 10] = [
        ([-0.48, -0.24], [0.48, -0.24], [56.0, 189.0, 248.0]),
        ([-0.48, -0.24], [-0.52, 0.26], [52.0, 211.0, 153.0]),
        ([0.48, -0.24], [0.52, 0.26], [52.0, 211.0, 153.0]),
        ([-0.48, -0.24], [0.0, 0.14], [167.0, 139.0, 250.0]),
        ([0.48, -0.24], [0.0, 0.14], [167.0, 139.0, 250.0]),
        ([-0.52, 0.26], [0.0, 0.14], [56.0, 189.0, 248.0]),
        ([0.52, 0.26], [0.0, 0.14], [56.0, 189.0, 248.0]),
        ([0.0, 0.14], [0.0, 0.56], [251.0, 191.0, 36.0]),
        ([-0.52, 0.26], [0.0, 0.56], [129.0, 140.0, 248.0]),
        ([0.52, 0.26], [0.0, 0.56], [129.0, 140.0, 248.0]),
    ];
    let hubs: [(f32, f32, f32, [f32; 3]); 6] = [
        (-0.48, -0.24, 0.125, [56.0, 189.0, 248.0]),
        (0.48, -0.24, 0.125, [56.0, 189.0, 248.0]),
        (-0.52, 0.26, 0.11, [52.0, 211.0, 153.0]),
        (0.52, 0.26, 0.11, [52.0, 211.0, 153.0]),
        (0.0, 0.14, 0.14, [251.0, 191.0, 36.0]),
        (0.0, 0.56, 0.105, [167.0, 139.0, 250.0]),
    ];

    for py in 0..size {
        let y = ((py as f32 + 0.5) / size as f32) * 2.0 - 1.0;
        for px in 0..size {
            let x = ((px as f32 + 0.5) / size as f32) * 2.0 - 1.0;
            let idx = (py * size + px) * 4;
            let box_d = sd_rounded_box(x, y, 0.88, 0.88, 0.24);
            if box_d < aa {
                let bg_a = (0.5 - box_d / aa).clamp(0.0, 1.0);
                blend_px(&mut buf, idx, 11.0, 16.0, 34.0, bg_a);
                if box_d.abs() < 0.035 {
                    blend_px(
                        &mut buf,
                        idx,
                        129.0,
                        140.0,
                        248.0,
                        (1.0 - box_d.abs() / 0.035) * 0.9 * bg_a,
                    );
                }
                let d_hemi = ((x * 1.05).hypot((y - 0.04) * 1.15) - 0.64).abs();
                if d_hemi < 0.032 {
                    blend_px(
                        &mut buf,
                        idx,
                        71.0,
                        85.0,
                        135.0,
                        (1.0 - d_hemi / 0.032) * 0.65 * bg_a,
                    );
                }
            }
            for &(a, b, col) in &tracts {
                let d = dist_seg(x, y, a[0], a[1], b[0], b[1]) - 0.032;
                if d < aa {
                    blend_px(
                        &mut buf,
                        idx,
                        col[0],
                        col[1],
                        col[2],
                        (0.5 - d / aa).clamp(0.0, 1.0) * 0.92,
                    );
                }
            }
            for &(hx, hy, hr, col) in &hubs {
                let d = (x - hx).hypot(y - hy) - hr;
                if d < aa {
                    blend_px(&mut buf, idx, col[0], col[1], col[2], (0.5 - d / aa).clamp(0.0, 1.0));
                }
            }
        }
    }
    buf
}

/// Returns the `64x64` window/taskbar icon for eframe `ViewportBuilder::with_icon`.
pub fn app_egui_icon() -> egui::IconData {
    let size = 64usize;
    egui::IconData {
        rgba: render_connectome_atlas_rgba(size),
        width: size as u32,
        height: size as u32,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn icon_rgba_is_64x64_and_non_empty() {
        let icon = app_egui_icon();
        assert_eq!(icon.width, 64);
        assert_eq!(icon.height, 64);
        assert_eq!(icon.rgba.len(), 64 * 64 * 4);
        assert!(icon.rgba.chunks_exact(4).any(|px| px[3] > 200));
    }
}
