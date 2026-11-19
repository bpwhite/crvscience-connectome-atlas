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

pub mod nifti_reader;
pub mod sqlite_store;
pub mod trk_reader;

pub use nifti_reader::{
    parse_nifti1_header, Nifti1HeaderSummary, Nifti1HeaderSummary as NiftiHeaderSummary,
    NiftiError, NiftiError as NiftiIoError,
};
pub use sqlite_store::{AtlasVault, SnapshotMeta, SnapshotMeta as VaultSnapshotMeta, VaultError};
pub use trk_reader::{
    format_swc_document, generate_region_swc_arbor, parse_swc_morphology, parse_trk_streamlines,
    SwcCompartment, TractIoError, TrkHeaderSummary,
};
