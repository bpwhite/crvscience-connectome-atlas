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

pub mod anatomy;
pub mod camera;
pub mod streamlines;
pub mod surfaces_and_glyphs;
pub mod wgpu_pipeline;

pub use anatomy::{
    build_species_anatomy, generate_2d_slice_contours, AnatomicalLandmark, AnatomicalLayer,
    AnatomicalMesh, SliceContour, SpeciesAnatomyScene,
};
pub use camera::{FlyThroughController, OrbitCamera, ViewPreset};
pub use streamlines::{
    generate_arched_tract_curve, generate_braided_fascicle_bundle, ColorSchemeMode,
    TractRenderStyle,
};
pub use surfaces_and_glyphs::{
    classify_region_organoid, compute_atlas_fiber_tensors, is_point_clipped_by_cutaway,
    project_parcel_scalars_to_surface, LocalFiberTensor, OrganoidClass, SurgicalCutawayMode,
};
pub use wgpu_pipeline::{encode_view_uniform_bytes, shader_descriptor, WGSL_CONNECTOME_SHADER};
