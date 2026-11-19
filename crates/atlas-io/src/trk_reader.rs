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

use atlas_core::{BrainRegion, ConnectomeAtlas};
use serde::{Deserialize, Serialize};
use thiserror::Error;

#[derive(Debug, Error)]
pub enum TractIoError {
    #[error("TrackVis .trk buffer shorter than 1000-byte binary header")]
    BufferTooShort,
    #[error("Invalid TrackVis magic string (expected 'TRACK\\0')")]
    InvalidTrackVisMagic,
    #[error("Truncated streamline polyline record in .trk payload")]
    TruncatedStreamline,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TrkHeaderSummary {
    pub dim: [i16; 3],
    pub voxel_size: [f32; 3],
    pub n_scalars: i16,
    pub n_properties: i16,
    pub n_count: i32,
    pub version: i32,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SwcCompartment {
    pub id: i64,
    pub compartment_type: u8,
    pub coord: [f32; 3],
    pub radius: f32,
    pub parent_id: i64,
}

/// Safe pure-Rust parser for TrackVis `.trk` diffusion MRI tractography streamlines.
pub fn parse_trk_streamlines(
    bytes: &[u8],
    max_streamlines: usize,
) -> Result<(TrkHeaderSummary, Vec<Vec<[f32; 3]>>), TractIoError> {
    if bytes.len() < 1000 {
        return Err(TractIoError::BufferTooShort);
    }
    if &bytes[0..6] != b"TRACK\0" {
        return Err(TractIoError::InvalidTrackVisMagic);
    }

    let rd_i16 = |off: usize| i16::from_le_bytes([bytes[off], bytes[off + 1]]);
    let rd_i32 = |off: usize| {
        i32::from_le_bytes([bytes[off], bytes[off + 1], bytes[off + 2], bytes[off + 3]])
    };
    let rd_f32 = |off: usize| {
        f32::from_le_bytes([bytes[off], bytes[off + 1], bytes[off + 2], bytes[off + 3]])
    };

    let summary = TrkHeaderSummary {
        dim: [rd_i16(6), rd_i16(8), rd_i16(10)],
        voxel_size: [rd_f32(12), rd_f32(16), rd_f32(20)],
        n_scalars: rd_i16(36),
        n_properties: rd_i16(238),
        n_count: rd_i32(988),
        version: rd_i32(992),
    };

    let hdr_size = rd_i32(996).max(1000) as usize;
    if bytes.len() < hdr_size {
        return Err(TractIoError::BufferTooShort);
    }

    let mut cursor = hdr_size;
    let mut streamlines = Vec::new();
    let floats_per_point = 3 + summary.n_scalars.max(0) as usize;
    let prop_bytes = (summary.n_properties.max(0) as usize) * 4;

    while cursor + 4 <= bytes.len() && streamlines.len() < max_streamlines {
        let m_pts = rd_i32(cursor);
        cursor += 4;
        if m_pts <= 0 {
            break;
        }
        let m = m_pts as usize;
        let needed = m * floats_per_point * 4 + prop_bytes;
        if cursor + needed > bytes.len() {
            return Err(TractIoError::TruncatedStreamline);
        }
        let mut pts = Vec::with_capacity(m);
        for p in 0..m {
            let base = cursor + p * floats_per_point * 4;
            pts.push([rd_f32(base), rd_f32(base + 4), rd_f32(base + 8)]);
        }
        streamlines.push(pts);
        cursor += needed;
    }

    Ok((summary, streamlines))
}

/// Parses standard NeuroMorpho / INCF `.swc` single-neuron reconstructions (`id type x y z r parent`).
pub fn parse_swc_morphology(text: &str) -> Vec<SwcCompartment> {
    let mut nodes = Vec::new();
    for line in text.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = trimmed.split_whitespace().collect();
        if cols.len() < 7 {
            continue;
        }
        let (Ok(id), Ok(ctype), Ok(x), Ok(y), Ok(z), Ok(r), Ok(parent)) = (
            cols[0].parse::<i64>(),
            cols[1].parse::<u8>(),
            cols[2].parse::<f32>(),
            cols[3].parse::<f32>(),
            cols[4].parse::<f32>(),
            cols[5].parse::<f32>(),
            cols[6].parse::<i64>(),
        ) else {
            continue;
        };
        nodes.push(SwcCompartment {
            id,
            compartment_type: ctype,
            coord: [x, y, z],
            radius: r,
            parent_id: parent,
        });
    }
    nodes
}

fn norm3(v: [f32; 3]) -> [f32; 3] {
    let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if len > 1e-5 {
        [v[0] / len, v[1] / len, v[2] / len]
    } else {
        [1.0, 0.0, 0.0]
    }
}

fn cross3(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Synthesizes a standardized 3D branched `.swc` neuronal arbor (`1=soma, 2=axon, 3=basal, 4=apical`)
/// oriented toward afferent presynaptic inputs and efferent postsynaptic targets.
pub fn generate_region_swc_arbor(region: &BrainRegion, atlas: &ConnectomeAtlas) -> Vec<SwcCompartment> {
    let mut in_vec = [0.0_f32; 3];
    let mut out_vec = [0.0_f32; 3];
    let mut in_cnt = 0_usize;
    let mut out_cnt = 0_usize;

    for e in &atlas.edges {
        if e.target == region.id {
            if let Some(src) = atlas.regions.iter().find(|r| r.id == e.source) {
                let d = norm3([
                    src.coord[0] - region.coord[0],
                    src.coord[1] - region.coord[1],
                    src.coord[2] - region.coord[2],
                ]);
                in_vec[0] += d[0];
                in_vec[1] += d[1];
                in_vec[2] += d[2];
                in_cnt += 1;
            }
        }
        if e.source == region.id {
            if let Some(dst) = atlas.regions.iter().find(|r| r.id == e.target) {
                let d = norm3([
                    dst.coord[0] - region.coord[0],
                    dst.coord[1] - region.coord[1],
                    dst.coord[2] - region.coord[2],
                ]);
                out_vec[0] += d[0];
                out_vec[1] += d[1];
                out_vec[2] += d[2];
                out_cnt += 1;
            }
        }
    }

    let in_dir = if in_cnt > 0 {
        norm3(in_vec)
    } else {
        norm3([0.45, 0.72, 0.52])
    };
    let out_dir = if out_cnt > 0 {
        norm3(out_vec)
    } else {
        norm3([-in_dir[0], -in_dir[1], -0.35])
    };

    let ortho_raw = cross3(in_dir, [0.0, 0.0, 1.0]);
    let ortho1 = if ortho_raw[0].abs() + ortho_raw[1].abs() + ortho_raw[2].abs() > 1e-3 {
        norm3(ortho_raw)
    } else {
        [1.0, 0.0, 0.0]
    };
    let ortho2 = norm3(cross3(in_dir, ortho1));

    let c = region.coord;
    let r_soma = (region.radius_mm * 0.44).max(1.35);
    let reach = (region.radius_mm * 1.38).max(3.6);

    let mut comps = Vec::with_capacity(13);
    comps.push(SwcCompartment {
        id: 1,
        compartment_type: 1,
        coord: c,
        radius: r_soma,
        parent_id: -1,
    });

    // Primary dendritic trunk (type = 4)
    let trunk = [
        c[0] + in_dir[0] * reach * 0.52,
        c[1] + in_dir[1] * reach * 0.52,
        c[2] + in_dir[2] * reach * 0.52,
    ];
    comps.push(SwcCompartment {
        id: 2,
        compartment_type: 4,
        coord: trunk,
        radius: r_soma * 0.42,
        parent_id: 1,
    });

    for sgn in [-1.0_f32, 1.0] {
        let b1 = [
            trunk[0] + in_dir[0] * reach * 0.46 + ortho1[0] * sgn * reach * 0.38,
            trunk[1] + in_dir[1] * reach * 0.46 + ortho1[1] * sgn * reach * 0.38,
            trunk[2] + in_dir[2] * reach * 0.46 + ortho1[2] * sgn * reach * 0.38,
        ];
        let b1_id = (comps.len() + 1) as i64;
        comps.push(SwcCompartment {
            id: b1_id,
            compartment_type: 4,
            coord: b1,
            radius: r_soma * 0.26,
            parent_id: 2,
        });
        for sub in [-1.0_f32, 1.0] {
            let twig = [
                b1[0] + in_dir[0] * reach * 0.34 + ortho2[0] * sub * reach * 0.28,
                b1[1] + in_dir[1] * reach * 0.34 + ortho2[1] * sub * reach * 0.28,
                b1[2] + in_dir[2] * reach * 0.34 + ortho2[2] * sub * reach * 0.28,
            ];
            let twig_id = (comps.len() + 1) as i64;
            comps.push(SwcCompartment {
                id: twig_id,
                compartment_type: 3,
                coord: twig,
                radius: r_soma * 0.15,
                parent_id: b1_id,
            });
        }
    }

    // Efferent axonal neurite (type = 2) + bifurcated terminal fork
    let ax_mid = [
        c[0] + out_dir[0] * reach * 0.68,
        c[1] + out_dir[1] * reach * 0.68,
        c[2] + out_dir[2] * reach * 0.68,
    ];
    let ax_id = (comps.len() + 1) as i64;
    comps.push(SwcCompartment {
        id: ax_id,
        compartment_type: 2,
        coord: ax_mid,
        radius: r_soma * 0.24,
        parent_id: 1,
    });
    for sgn in [-1.0_f32, 1.0] {
        let term = [
            ax_mid[0] + out_dir[0] * reach * 0.50 + ortho2[0] * sgn * reach * 0.30,
            ax_mid[1] + out_dir[1] * reach * 0.50 + ortho2[1] * sgn * reach * 0.30,
            ax_mid[2] + out_dir[2] * reach * 0.50 + ortho2[2] * sgn * reach * 0.30,
        ];
        let term_id = (comps.len() + 1) as i64;
        comps.push(SwcCompartment {
            id: term_id,
            compartment_type: 2,
            coord: term,
            radius: r_soma * 0.16,
            parent_id: ax_id,
        });
    }

    comps
}

/// Formats an `.swc` compartment tree into a standard 7-column NeuroMorpho / INCF `.swc` text file.
pub fn format_swc_document(
    region: &BrainRegion,
    organism: &str,
    compartments: &[SwcCompartment],
) -> String {
    let mut out = format!(
        "# ConnectomeAtlas v0.5.0 Standardized .SWC Neuronal Arbor\n# Organism: {organism} | Module: {} ({})\n# id type x y z radius parent_id\n",
        region.id, region.name
    );
    for c in compartments {
        out.push_str(&format!(
            "{} {} {:.3} {:.3} {:.3} {:.3} {}\n",
            c.id, c.compartment_type, c.coord[0], c.coord[1], c.coord[2], c.radius, c.parent_id
        ));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use atlas_core::load_preset_by_id;

    #[test]
    fn parses_synthetic_trackvis_trk_and_roundtrips_swc_arbor() {
        let mut buf = vec![0_u8; 1000];
        buf[0..6].copy_from_slice(b"TRACK\0");
        buf[6..8].copy_from_slice(&182_i16.to_le_bytes());
        buf[8..10].copy_from_slice(&218_i16.to_le_bytes());
        buf[10..12].copy_from_slice(&182_i16.to_le_bytes());
        buf[988..992].copy_from_slice(&1_i32.to_le_bytes());
        buf[992..996].copy_from_slice(&2_i32.to_le_bytes());
        buf[996..1000].copy_from_slice(&1000_i32.to_le_bytes());
        buf.extend_from_slice(&2_i32.to_le_bytes());
        for val in [-20.0_f32, -4.0, -15.0, 1.0, 55.0, -3.0] {
            buf.extend_from_slice(&val.to_le_bytes());
        }

        let (hdr, tracks) = parse_trk_streamlines(&buf, 10).expect("valid trk");
        assert_eq!(hdr.dim, [182, 218, 182]);
        assert_eq!(tracks.len(), 1);
        assert_eq!(tracks[0].len(), 2);

        let worm = load_preset_by_id("celegans_wormatlas").expect("worm preset");
        let reg = &worm.regions[0];
        let arbor = generate_region_swc_arbor(reg, &worm);
        assert!(arbor.len() >= 10);
        let doc = format_swc_document(reg, &worm.organism, &arbor);
        let parsed = parse_swc_morphology(&doc);
        assert_eq!(parsed.len(), arbor.len());
        assert_eq!(parsed[0].compartment_type, 1);
        assert_eq!(parsed[0].parent_id, -1);
    }
}
