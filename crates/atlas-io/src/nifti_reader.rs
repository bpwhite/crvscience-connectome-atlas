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

use flate2::read::GzDecoder;
use serde::{Deserialize, Serialize};
use std::io::Read;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum NiftiError {
    #[error("I/O error reading NIfTI stream: {0}")]
    Io(#[from] std::io::Error),
    #[error("Invalid NIfTI-1 header size (expected 348 bytes, got {actual})")]
    InvalidHeaderSize { actual: i32 },
    #[error("Missing NIfTI-1 magic tag ('n+1' or 'ni1')")]
    InvalidMagic,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Nifti1HeaderSummary {
    pub dimensions: [u16; 4],
    pub voxel_spacing_mm: [f32; 3],
    pub datatype_code: i16,
    pub bitpix: i16,
    pub sform_code: i16,
    pub qform_code: i16,
    pub srow_x: [f32; 4],
    pub srow_y: [f32; 4],
    pub srow_z: [f32; 4],
    pub description: String,
}

impl Nifti1HeaderSummary {
    /// Transforms voxel index `(i, j, k)` into RAS stereotaxic millimeters using `sform`.
    pub fn voxel_to_ras_mm(&self, i: f32, j: f32, k: f32) -> [f32; 3] {
        if self.sform_code > 0 {
            [
                self.srow_x[0] * i + self.srow_x[1] * j + self.srow_x[2] * k + self.srow_x[3],
                self.srow_y[0] * i + self.srow_y[1] * j + self.srow_y[2] * k + self.srow_y[3],
                self.srow_z[0] * i + self.srow_z[1] * j + self.srow_z[2] * k + self.srow_z[3],
            ]
        } else {
            [
                i * self.voxel_spacing_mm[0],
                j * self.voxel_spacing_mm[1],
                k * self.voxel_spacing_mm[2],
            ]
        }
    }
}

/// Safe pure-Rust parser for NIfTI-1 `.nii` and gzipped `.nii.gz` 348-byte headers.
pub fn parse_nifti1_header(raw_bytes: &[u8]) -> Result<Nifti1HeaderSummary, NiftiError> {
    let mut header = [0_u8; 352];
    if raw_bytes.len() >= 2 && raw_bytes[0] == 0x1f && raw_bytes[1] == 0x8b {
        let mut decoder = GzDecoder::new(raw_bytes);
        decoder.read_exact(&mut header)?;
    } else {
        if raw_bytes.len() < 348 {
            return Err(NiftiError::InvalidHeaderSize {
                actual: raw_bytes.len() as i32,
            });
        }
        header[..348].copy_from_slice(&raw_bytes[..348]);
    }

    let sizeof_hdr_le = i32::from_le_bytes([header[0], header[1], header[2], header[3]]);
    let little_endian = if sizeof_hdr_le == 348 {
        true
    } else if i32::from_be_bytes([header[0], header[1], header[2], header[3]]) == 348 {
        false
    } else {
        return Err(NiftiError::InvalidHeaderSize {
            actual: sizeof_hdr_le,
        });
    };

    let magic = &header[344..348];
    if magic != b"n+1\0" && magic != b"ni1\0" {
        return Err(NiftiError::InvalidMagic);
    }

    let rd_i16 = |off: usize| -> i16 {
        let b = [header[off], header[off + 1]];
        if little_endian {
            i16::from_le_bytes(b)
        } else {
            i16::from_be_bytes(b)
        }
    };
    let rd_f32 = |off: usize| -> f32 {
        let b = [header[off], header[off + 1], header[off + 2], header[off + 3]];
        if little_endian {
            f32::from_le_bytes(b)
        } else {
            f32::from_be_bytes(b)
        }
    };

    let dimensions = [
        rd_i16(42).max(1) as u16,
        rd_i16(44).max(1) as u16,
        rd_i16(46).max(1) as u16,
        rd_i16(48).max(1) as u16,
    ];
    let datatype_code = rd_i16(70);
    let bitpix = rd_i16(72);
    let voxel_spacing_mm = [rd_f32(80), rd_f32(84), rd_f32(88)];
    let qform_code = rd_i16(252);
    let sform_code = rd_i16(254);

    let srow_x = [rd_f32(280), rd_f32(284), rd_f32(288), rd_f32(292)];
    let srow_y = [rd_f32(296), rd_f32(300), rd_f32(304), rd_f32(308)];
    let srow_z = [rd_f32(312), rd_f32(316), rd_f32(320), rd_f32(324)];

    let desc_slice = &header[148..228];
    let end = desc_slice.iter().position(|&c| c == 0).unwrap_or(desc_slice.len());
    let description = String::from_utf8_lossy(&desc_slice[..end]).trim().to_string();

    Ok(Nifti1HeaderSummary {
        dimensions,
        voxel_spacing_mm,
        datatype_code,
        bitpix,
        sform_code,
        qform_code,
        srow_x,
        srow_y,
        srow_z,
        description,
    })
}
