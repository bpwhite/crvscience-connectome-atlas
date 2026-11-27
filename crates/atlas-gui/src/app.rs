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

use atlas_core::{
    all_presets, build_command_index, compute_comparative_topology_digest, search_command_index,
    AtlasEnvelope, AtlasSummaryMetrics, CommandAction, CommandEntry, ComparativeTopologyDigest,
    ConnectomeAtlas, ControlHarmonicsReport, CriticalityReport, DecColor, DomainDecodeEntry,
    LesionDeltaReport, PathResult, SystemsNeuroReport, TopologyTdaReport,
};
use atlas_io::{parse_nifti1_header, parse_swc_morphology, parse_trk_streamlines, AtlasVault};
use atlas_io::{format_swc_document, generate_region_swc_arbor};
use atlas_render::{
    build_species_anatomy, compute_atlas_fiber_tensors, generate_2d_slice_contours,
    generate_arched_tract_curve, generate_braided_fascicle_bundle, generate_tensor_ellipsoid_mesh,
    is_point_clipped_by_cutaway, AnatomicalLayer, ColorSchemeMode, FlyThroughController,
    LocalFiberTensor, OrbitCamera, SpeciesAnatomyScene, SurgicalCutawayMode, TractRenderStyle,
    ViewPreset,
};
use eframe::egui::epaint::{Mesh as EguiMesh, Vertex as EguiVertex, WHITE_UV};
use eframe::egui::{self, Color32, FontId, Pos2, Rect, Sense, Stroke, Vec2};
use glam::Vec3;
use std::collections::{HashMap, HashSet};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ProbeMode {
    Inspect,
    ShortestPath,
    LesionAblate,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum RightTab {
    Inspector,
    AdjacencyMatrix,
    RichClubHubs,
    CriticalityTda,
    VaultIo,
}

pub struct ConnectomeAtlasApp {
    presets: Vec<ConnectomeAtlas>,
    active_preset_idx: usize,
    atlas: ConnectomeAtlas,
    metrics: AtlasSummaryMetrics,
    lesion_delta: LesionDeltaReport,
    control_harmonics: ControlHarmonicsReport,
    systems_neuro: SystemsNeuroReport,
    topology_digest: ComparativeTopologyDigest,
    criticality_kappa: f64,
    criticality: CriticalityReport,
    topology_tda: TopologyTdaReport,
    domain_decode: Vec<DomainDecodeEntry>,
    fiber_tensors: Vec<LocalFiberTensor>,
    anatomy_scene: SpeciesAnatomyScene,
    tract_style: TractRenderStyle,
    cutaway_mode: SurgicalCutawayMode,
    cutaway_offset_mm: f32,
    show_tensor_glyphs: bool,
    show_internal_organs: bool,
    show_anatomical_labels: bool,
    fly_ctrl: FlyThroughController,
    receptor_channel: String,
    agonist_gain: f64,
    harmonic_mode_idx: usize,
    camera: OrbitCamera,
    color_mode: ColorSchemeMode,
    probe_mode: ProbeMode,
    right_tab: RightTab,
    min_weight: f64,
    shell_opacity: f32,
    flow_speed: f32,
    tonic_phase: f32,
    conduction_velocity_m_s: f32,
    show_slice_planes: bool,
    slice_x_mm: f32,
    slice_y_mm: f32,
    slice_z_mm: f32,
    selected_region: Option<String>,
    path_source: Option<String>,
    path_target: Option<String>,
    active_path: Option<PathResult>,
    vault: Option<AtlasVault>,
    show_left_panel: bool,
    show_right_panel: bool,
    show_bottom_slices: bool,
    zen_mode: bool,
    cmd_palette_open: bool,
    show_about_window: bool,
    cmd_query: String,
    cmd_selected_idx: usize,
    cmd_entries: Vec<CommandEntry>,
    pending_scope_open: Option<String>,
    slot_name_buf: String,
    status_msg: String,
}

impl ConnectomeAtlasApp {
    #[allow(dead_code)]
    pub fn new(cc: &eframe::CreationContext<'_>) -> Self {
        Self::new_with_db(cc, std::path::PathBuf::from("connectome_vault.sqlite"), None)
    }

    pub fn new_with_db(
        cc: &eframe::CreationContext<'_>,
        db_path: std::path::PathBuf,
        initial_preset: Option<&str>,
    ) -> Self {
        let mut visuals = egui::Visuals::dark();
        visuals.panel_fill = Color32::from_rgb(13, 17, 26);
        visuals.window_fill = Color32::from_rgb(16, 21, 32);
        cc.egui_ctx.set_visuals(visuals);

        let presets = all_presets();
        let active_preset_idx = initial_preset
            .and_then(|pid| presets.iter().position(|p| p.id == pid))
            .unwrap_or(0);
        let vault = AtlasVault::open(db_path).ok();
        Self::init_from_preset(
            presets,
            active_preset_idx,
            vault,
            0.28,
            0.0,
            "Ready - Press Ctrl+K to search, F1 for shortcuts, or drop files onto viewport.",
        )
    }

    pub fn new_headless(preset_idx: usize) -> Self {
        let presets = all_presets();
        let idx = preset_idx.min(presets.len().saturating_sub(1));
        let vault = AtlasVault::open_in_memory().ok();
        Self::init_from_preset(
            presets,
            idx,
            vault,
            0.0,
            0.15,
            "Headless ConnectomeAtlas snapshot ready.",
        )
    }

    fn init_from_preset(
        presets: Vec<ConnectomeAtlas>,
        active_preset_idx: usize,
        vault: Option<AtlasVault>,
        flow_speed: f32,
        tonic_phase: f32,
        status_msg: &str,
    ) -> Self {
        let atlas = presets[active_preset_idx].clone();
        let metrics = atlas.compute_summary_metrics(0.0);
        let lesion_delta = LesionDeltaReport {
            intact_efficiency: metrics.global_efficiency,
            lesioned_efficiency: metrics.global_efficiency,
            efficiency_retention_pct: 100.0,
            intact_components: metrics.connected_components,
            lesioned_components: metrics.connected_components,
        };
        let control_harmonics = atlas.compute_control_and_harmonics(0.0, "ht", 0.0);
        let systems_neuro = atlas.compute_systems_neuro_report(0.0, 1337);
        let topology_digest = compute_comparative_topology_digest(&presets, 0.0, 1337);
        let criticality = atlas.compute_criticality_report(0.0, 1.0);
        let topology_tda = atlas.compute_topology_tda_report(0.0);
        let domain_decode = atlas.decode_functional_domains(None, 0.0);
        let fiber_tensors = compute_atlas_fiber_tensors(&atlas, 0.0);
        let anatomy_scene = build_species_anatomy(&atlas.id);
        let selected_region = atlas.regions.first().map(|r| r.id.clone());
        let (path_source, path_target) = if active_preset_idx == 0 {
            (
                Some("LANG.Wernicke_L".to_string()),
                Some("LANG.Broca_L".to_string()),
            )
        } else {
            (
                atlas.regions.first().map(|r| r.id.clone()),
                atlas.regions.last().map(|r| r.id.clone()),
            )
        };
        let active_path = match (&path_source, &path_target) {
            (Some(s), Some(t)) => atlas.shortest_signaling_path(s, t, 0.0),
            _ => None,
        };
        let cmd_entries = build_command_index(&atlas);
        let mut camera = OrbitCamera::default();
        camera.yaw_rad = anatomy_scene.default_yaw;
        camera.pitch_rad = anatomy_scene.default_pitch;
        camera.distance = anatomy_scene.default_distance;
        camera.target = Vec3::from_array(anatomy_scene.default_target);
        let default_slice = anatomy_scene.default_slice;
        let vel = atlas.conduction_velocity_m_s as f32;

        Self {
            presets,
            active_preset_idx,
            atlas,
            metrics,
            lesion_delta,
            control_harmonics,
            systems_neuro,
            topology_digest,
            criticality_kappa: 1.0,
            criticality,
            topology_tda,
            domain_decode,
            fiber_tensors,
            anatomy_scene,
            tract_style: TractRenderStyle::BraidedFascicles,
            cutaway_mode: SurgicalCutawayMode::AnteriorSuperiorOctant,
            cutaway_offset_mm: 0.0,
            show_tensor_glyphs: true,
            show_internal_organs: true,
            show_anatomical_labels: true,
            fly_ctrl: FlyThroughController::default(),
            receptor_channel: "ht".to_string(),
            agonist_gain: 0.0,
            harmonic_mode_idx: 1,
            camera,
            color_mode: ColorSchemeMode::FunctionalNetwork,
            probe_mode: ProbeMode::Inspect,
            right_tab: RightTab::Inspector,
            min_weight: 0.0,
            shell_opacity: 0.45,
            flow_speed,
            tonic_phase,
            conduction_velocity_m_s: vel,
            show_slice_planes: true,
            slice_x_mm: default_slice[0],
            slice_y_mm: default_slice[1],
            slice_z_mm: default_slice[2],
            selected_region,
            path_source,
            path_target,
            active_path,
            vault,
            show_left_panel: true,
            show_right_panel: true,
            show_bottom_slices: true,
            zen_mode: false,
            cmd_palette_open: false,
            show_about_window: false,
            cmd_query: String::new(),
            cmd_selected_idx: 0,
            cmd_entries,
            pending_scope_open: None,
            slot_name_buf: "session_snapshot_01".to_string(),
            status_msg: status_msg.to_string(),
        }
    }

    fn recompute(&mut self) {
        self.atlas.conduction_velocity_m_s = self.conduction_velocity_m_s as f64;
        self.metrics = self.atlas.compute_summary_metrics(self.min_weight);
        self.lesion_delta =
            if self.atlas.ablated_regions.is_empty() && self.atlas.severed_tracts.is_empty() {
                LesionDeltaReport {
                    intact_efficiency: self.metrics.global_efficiency,
                    lesioned_efficiency: self.metrics.global_efficiency,
                    efficiency_retention_pct: 100.0,
                    intact_components: self.metrics.connected_components,
                    lesioned_components: self.metrics.connected_components,
                }
            } else {
                self.atlas.evaluate_lesion_impact(self.min_weight)
            };
        self.control_harmonics = self.atlas.compute_control_and_harmonics(
            self.min_weight,
            &self.receptor_channel,
            self.agonist_gain,
        );
        let max_h = self.control_harmonics.harmonic_modes.len().max(1);
        self.harmonic_mode_idx = self.harmonic_mode_idx.clamp(1, max_h);
        self.fiber_tensors = compute_atlas_fiber_tensors(&self.atlas, self.min_weight);
        self.systems_neuro = self.atlas.compute_systems_neuro_report(self.min_weight, 1337);
        let mut digest_atlases = self.presets.clone();
        if let Some(slot) = digest_atlases.get_mut(self.active_preset_idx) {
            *slot = self.atlas.clone();
        }
        self.topology_digest = compute_comparative_topology_digest(&digest_atlases, self.min_weight, 1337);
        self.criticality = self
            .atlas
            .compute_criticality_report(self.min_weight, self.criticality_kappa);
        self.topology_tda = self.atlas.compute_topology_tda_report(self.min_weight);
        self.domain_decode = self
            .atlas
            .decode_functional_domains(self.selected_region.as_deref(), self.min_weight);
        self.cmd_entries = build_command_index(&self.atlas);
        if let (Some(s), Some(t)) = (&self.path_source, &self.path_target) {
            self.active_path = self.atlas.shortest_signaling_path(s, t, self.min_weight);
        } else {
            self.active_path = None;
        }
    }

    fn set_view_preset(&mut self, preset: ViewPreset, smooth: bool) {
        self.fly_ctrl.active = false;
        if smooth {
            self.camera.start_preset_transition(preset);
        } else {
            self.camera.apply_preset(preset);
        }
        let label = match preset {
            ViewPreset::Perspective3D => "3D Orbit (1)",
            ViewPreset::LeftSagittal => "Left Sagittal (2)",
            ViewPreset::SuperiorAxial => "Superior Axial (3)",
            ViewPreset::AnteriorCoronal => "Anterior Coronal (4)",
        };
        self.status_msg = format!("Camera preset: {label}");
    }

    fn reset_camera_to_species(&mut self, smooth: bool) {
        self.fly_ctrl.active = false;
        let target = Vec3::from_array(self.anatomy_scene.default_target);
        if smooth {
            self.camera.start_transition(
                self.anatomy_scene.default_yaw,
                self.anatomy_scene.default_pitch,
                Some(self.anatomy_scene.default_distance),
                Some(target),
            );
        } else {
            self.camera.cancel_transition();
            self.camera.yaw_rad = self.anatomy_scene.default_yaw;
            self.camera.pitch_rad = self.anatomy_scene.default_pitch;
            self.camera.distance = self.anatomy_scene.default_distance;
            self.camera.target = target;
        }
        self.status_msg = format!("Reset camera to {} default view (R).", self.atlas.organism);
    }

    fn cycle_cutaway(&mut self) {
        self.cutaway_mode = match self.cutaway_mode {
            SurgicalCutawayMode::Intact => SurgicalCutawayMode::AnteriorSuperiorOctant,
            SurgicalCutawayMode::AnteriorSuperiorOctant => SurgicalCutawayMode::LeftHemisphere,
            SurgicalCutawayMode::LeftHemisphere => SurgicalCutawayMode::SuperiorAxialDome,
            SurgicalCutawayMode::SuperiorAxialDome => SurgicalCutawayMode::Intact,
        };
        let label = match self.cutaway_mode {
            SurgicalCutawayMode::Intact => "Intact Shell",
            SurgicalCutawayMode::AnteriorSuperiorOctant => "Anterior-Superior Octant",
            SurgicalCutawayMode::LeftHemisphere => "Left Hemisphere Window",
            SurgicalCutawayMode::SuperiorAxialDome => "Superior Axial Dome",
        };
        self.status_msg = format!("Surgical cutaway: {label} (C).");
    }

    fn toggle_zen_mode(&mut self) {
        self.zen_mode = !self.zen_mode;
        let vis = !self.zen_mode;
        self.show_left_panel = vis;
        self.show_right_panel = vis;
        self.show_bottom_slices = vis;
    }

    fn apply_species_framing(&mut self) {
        self.fly_ctrl.active = false;
        self.camera.cancel_transition();
        self.anatomy_scene = build_species_anatomy(&self.atlas.id);
        self.camera.target = Vec3::from_array(self.anatomy_scene.default_target);
        self.camera.distance = self.anatomy_scene.default_distance;
        self.camera.yaw_rad = self.anatomy_scene.default_yaw;
        self.camera.pitch_rad = self.anatomy_scene.default_pitch;
        self.slice_x_mm = self.anatomy_scene.default_slice[0];
        self.slice_y_mm = self.anatomy_scene.default_slice[1];
        self.slice_z_mm = self.anatomy_scene.default_slice[2];
    }

    fn switch_preset(&mut self, idx: usize) {
        self.active_preset_idx = idx;
        self.atlas = self.presets[idx].clone();
        self.conduction_velocity_m_s = self.atlas.conduction_velocity_m_s as f32;
        self.selected_region = self.atlas.regions.first().map(|r| r.id.clone());
        if idx == 0 {
            self.path_source = Some("LANG.Wernicke_L".to_string());
            self.path_target = Some("LANG.Broca_L".to_string());
        } else {
            self.path_source = self.atlas.regions.first().map(|r| r.id.clone());
            self.path_target = self.atlas.regions.last().map(|r| r.id.clone());
        }
        self.apply_species_framing();
        self.recompute();
    }

    fn load_file_from_path(&mut self, path: &std::path::Path) {
        let name = path
            .file_name()
            .and_then(|s| s.to_str())
            .unwrap_or("file")
            .to_lowercase();
        let bytes = match std::fs::read(path) {
            Ok(b) => b,
            Err(e) => {
                self.status_msg = format!("Failed reading {}: {e}", path.display());
                return;
            }
        };
        if name.ends_with(".json") {
            match String::from_utf8(bytes)
                .map_err(|e| format!("Invalid UTF-8 JSON: {e}"))
                .and_then(|txt| {
                    serde_json::from_str::<AtlasEnvelope>(&txt)
                        .map_err(|e| format!("Invalid connectome envelope: {e}"))
                }) {
                Ok(env) => {
                    self.atlas = env.atlas;
                    self.conduction_velocity_m_s = self.atlas.conduction_velocity_m_s as f32;
                    self.selected_region = self.atlas.regions.first().map(|r| r.id.clone());
                    self.path_source = self.atlas.regions.first().map(|r| r.id.clone());
                    self.path_target = self.atlas.regions.last().map(|r| r.id.clone());
                    self.apply_species_framing();
                    self.recompute();
                    self.status_msg = format!("Loaded envelope {}", path.display());
                }
                Err(e) => self.status_msg = e,
            }
        } else if name.ends_with(".swc") {
            match String::from_utf8(bytes) {
                Ok(txt) => {
                    let comps = parse_swc_morphology(&txt);
                    self.status_msg =
                        format!("Parsed .SWC arbor ({} compartments).", comps.len());
                }
                Err(e) => self.status_msg = format!("Invalid UTF-8 SWC: {e}"),
            }
        } else if name.ends_with(".trk") {
            match parse_trk_streamlines(&bytes, 250) {
                Ok((hdr, tracks)) => {
                    self.status_msg =
                        format!("Parsed .trk ({:?}, {} streamlines).", hdr.dim, tracks.len());
                }
                Err(e) => self.status_msg = format!("Invalid .trk stream: {e}"),
            }
        } else if name.ends_with(".nii") || name.ends_with(".nii.gz") {
            match parse_nifti1_header(&bytes) {
                Ok(hdr) => {
                    self.status_msg =
                        format!("Inspected NIfTI-1 volume {:?}.", hdr.dimensions);
                }
                Err(e) => self.status_msg = format!("Invalid NIfTI-1 file: {e}"),
            }
        } else {
            self.status_msg = format!("Unsupported file type: {}", path.display());
        }
    }

    pub fn select_region(&mut self, rid: String, sync_slices: bool) {
        if sync_slices {
            if let Some(reg) = self.atlas.regions.iter().find(|r| r.id == rid) {
                self.slice_x_mm = reg.coord[0];
                self.slice_y_mm = reg.coord[1];
                self.slice_z_mm = reg.coord[2];
            }
        }
        self.selected_region = Some(rid);
        self.domain_decode = self
            .atlas
            .decode_functional_domains(self.selected_region.as_deref(), self.min_weight);
    }

    fn open_dataset_dialog(&mut self) {
        if let Some(path) = rfd::FileDialog::new()
            .add_filter("Connectome & Neuroimaging", &["json", "swc", "trk", "nii", "gz"])
            .pick_file()
        {
            self.load_file_from_path(&path);
        }
    }

    fn export_envelope_dialog(&mut self) {
        let default_name = format!("{}.connectome.json", self.atlas.id);
        let Some(path) = rfd::FileDialog::new()
            .add_filter("Connectome JSON Envelope", &["json"])
            .set_file_name(&default_name)
            .save_file()
        else {
            return;
        };
        match AtlasEnvelope::wrap(self.atlas.clone())
            .map_err(|e| format!("Envelope wrap error: {e}"))
            .and_then(|env| {
                serde_json::to_string_pretty(&env)
                    .map_err(|e| format!("JSON serialize error: {e}"))
            })
            .and_then(|txt| {
                std::fs::write(&path, txt).map_err(|e| format!("Write failed: {e}"))
            }) {
            Ok(()) => {
                self.status_msg = format!("Exported envelope to {}", path.display());
            }
            Err(e) => self.status_msg = e,
        }
    }

    fn export_swc_dialog(&mut self) {
        let Some(sel_id) = &self.selected_region else {
            self.status_msg = "Select a brain region before exporting an .SWC arbor.".to_string();
            return;
        };
        let Some(reg) = self.atlas.regions.iter().find(|r| &r.id == sel_id) else {
            self.status_msg = format!("Region '{sel_id}' not found in active atlas.");
            return;
        };
        let arbor = generate_region_swc_arbor(reg, &self.atlas);
        let doc = format_swc_document(reg, &self.atlas.organism, &arbor);
        let default_name = format!("{}_arbor.swc", reg.id.replace('.', "_"));
        let Some(path) = rfd::FileDialog::new()
            .add_filter("SWC Morphology", &["swc"])
            .set_file_name(&default_name)
            .save_file()
        else {
            return;
        };
        match std::fs::write(&path, doc) {
            Ok(()) => {
                self.status_msg = format!("Exported 3D SWC arbor to {}", path.display());
            }
            Err(e) => self.status_msg = format!("SWC export failed: {e}"),
        }
    }

    fn export_topology_digest_dialog(&mut self, is_csv: bool) {
        let (default_name, filter_lbl, ext, payload_res) = if is_csv {
            (
                "connectome_topology_digest.csv",
                "CSV Table",
                "csv",
                Ok(self.topology_digest.to_csv()),
            )
        } else {
            (
                "connectome_topology_digest.json",
                "JSON Digest",
                "json",
                self.topology_digest
                    .to_json_pretty()
                    .map_err(|e| format!("JSON serialize error: {e}")),
            )
        };
        let payload = match payload_res {
            Ok(txt) => txt,
            Err(e) => {
                self.status_msg = e;
                return;
            }
        };
        let Some(path) = rfd::FileDialog::new()
            .add_filter(filter_lbl, &[ext])
            .set_file_name(default_name)
            .save_file()
        else {
            return;
        };
        match std::fs::write(&path, payload) {
            Ok(()) => {
                self.status_msg = format!(
                    "Exported comparative topology digest ({}) to {}",
                    ext.to_uppercase(),
                    path.display()
                );
            }
            Err(e) => self.status_msg = format!("Digest export failed: {e}"),
        }
    }

    fn execute_command_action(&mut self, act: CommandAction, animate: bool) {
        match act {
            CommandAction::JumpToRegion(rid) => {
                self.select_region(rid.clone(), true);
                if let Some(reg) = self.atlas.regions.iter().find(|r| r.id == rid) {
                    let focus = Vec3::from_array(reg.coord);
                    self.fly_ctrl.active = false;
                    if animate {
                        self.camera.start_transition(
                            self.camera.yaw_rad,
                            self.camera.pitch_rad,
                            None,
                            Some(focus),
                        );
                    } else {
                        self.camera.target = focus;
                    }
                }
                self.status_msg = format!("Jumped to {rid}.");
            }
            CommandAction::SwitchPreset(pid) => {
                if let Some(idx) = self.presets.iter().position(|p| p.id == pid) {
                    self.switch_preset(idx);
                    self.status_msg = format!("Switched preset to {}.", self.atlas.organism);
                }
            }
            CommandAction::ToggleZenMode => {
                self.toggle_zen_mode();
            }
            CommandAction::SetControlScope(scope) => {
                self.zen_mode = false;
                self.show_left_panel = true;
                self.pending_scope_open = Some(scope.clone());
                self.status_msg = format!("Focused control scope: {scope}.");
            }
            CommandAction::SetCutaway(m) => {
                self.cutaway_mode = match m.as_str() {
                    "octant" => SurgicalCutawayMode::AnteriorSuperiorOctant,
                    "left_hemi" => SurgicalCutawayMode::LeftHemisphere,
                    "axial_dome" => SurgicalCutawayMode::SuperiorAxialDome,
                    _ => SurgicalCutawayMode::Intact,
                };
                self.status_msg = format!("Set surgical cutaway mode: {m}.");
            }
            CommandAction::ToggleFlyThrough => {
                self.fly_ctrl.active = !self.fly_ctrl.active;
                self.fly_ctrl.progress_u = 0.0;
                if self.fly_ctrl.active {
                    self.camera.cancel_transition();
                }
            }
            CommandAction::HealAllLesions => {
                self.atlas.heal_all();
                self.recompute();
                self.status_msg = "Restored all ablated regions and severed tracts.".to_string();
            }
            CommandAction::ExportEnvelopeJson => {
                self.export_envelope_dialog();
            }
            CommandAction::ExportSelectedSwc => {
                self.export_swc_dialog();
            }
            CommandAction::TuneCriticalityKappa(k_pct) => {
                self.criticality_kappa = (k_pct as f64 / 100.0).clamp(0.35, 2.0);
                self.recompute();
                self.zen_mode = false;
                self.show_right_panel = true;
                self.right_tab = RightTab::CriticalityTda;
                self.color_mode = ColorSchemeMode::CriticalitySusceptibility;
                self.status_msg = format!(
                    "Tuned criticality kappa = {:.2} (sigma = {:.3}).",
                    self.criticality_kappa, self.criticality.branching_ratio_sigma
                );
            }
            CommandAction::InspectPersistentCavities => {
                self.zen_mode = false;
                self.show_right_panel = true;
                self.right_tab = RightTab::CriticalityTda;
                self.color_mode = ColorSchemeMode::TdaBeta1Loops;
                self.status_msg = format!(
                    "Inspecting {} persistent beta-1 topological cavities.",
                    self.topology_tda.persistent_beta1_count
                );
            }
            CommandAction::InspectTopologyDigest => {
                self.zen_mode = false;
                self.show_right_panel = true;
                self.right_tab = RightTab::RichClubHubs;
                self.status_msg = format!(
                    "Viewing comparative topology digest ({} presets).",
                    self.topology_digest.entries.len()
                );
            }
            CommandAction::ExportTopologyDigestJson => {
                self.export_topology_digest_dialog(false);
            }
            CommandAction::ExportTopologyDigestCsv => {
                self.export_topology_digest_dialog(true);
            }
        }
    }

    pub fn ui_step(&mut self, ctx: &egui::Context, animate: bool) {
        let dt = ctx.input(|i| i.stable_dt).min(0.05);
        if animate {
            let cam_anim = self.camera.step_transition(dt);
            if self.flow_speed > 0.001 || self.fly_ctrl.active || cam_anim {
                if self.flow_speed > 0.001 || self.fly_ctrl.active {
                    self.tonic_phase = (self.tonic_phase + dt * self.flow_speed) % 1.0;
                }
                if self.fly_ctrl.active {
                    self.fly_ctrl.progress_u =
                        (self.fly_ctrl.progress_u + dt * self.fly_ctrl.speed) % 1.0;
                }
                ctx.request_repaint();
            }
        }

        // Global Keyboard Shortcuts (Escape, F1, Ctrl/Cmd+K/O/S)
        if ctx.input(|i| i.key_pressed(egui::Key::Escape)) {
            if self.cmd_palette_open {
                self.cmd_palette_open = false;
            } else if self.show_about_window {
                self.show_about_window = false;
            } else if self.fly_ctrl.active {
                self.fly_ctrl.active = false;
            } else if self.zen_mode {
                self.zen_mode = false;
                self.show_left_panel = true;
                self.show_right_panel = true;
                self.show_bottom_slices = true;
            }
        }
        if ctx.input(|i| i.key_pressed(egui::Key::F1)) {
            self.show_about_window = !self.show_about_window;
        }
        if ctx.input(|i| (i.modifiers.command || i.modifiers.ctrl) && i.key_pressed(egui::Key::K)) {
            self.cmd_palette_open = !self.cmd_palette_open;
            self.cmd_query.clear();
            self.cmd_selected_idx = 0;
        }
        if ctx.input(|i| (i.modifiers.command || i.modifiers.ctrl) && i.key_pressed(egui::Key::O)) {
            self.open_dataset_dialog();
        }
        if ctx.input(|i| (i.modifiers.command || i.modifiers.ctrl) && i.key_pressed(egui::Key::S)) {
            self.export_envelope_dialog();
        }

        // Viewport Hotkeys guarded against active text inputs & modifier chords
        let plain_key_allowed = !ctx.wants_keyboard_input()
            && ctx.input(|i| !i.modifiers.ctrl && !i.modifiers.command && !i.modifiers.alt);
        if plain_key_allowed {
            if ctx.input(|i| i.key_pressed(egui::Key::Z)) {
                self.toggle_zen_mode();
            } else if ctx.input(|i| i.key_pressed(egui::Key::Num1)) {
                self.set_view_preset(ViewPreset::Perspective3D, animate);
            } else if ctx.input(|i| i.key_pressed(egui::Key::Num2)) {
                self.set_view_preset(ViewPreset::LeftSagittal, animate);
            } else if ctx.input(|i| i.key_pressed(egui::Key::Num3)) {
                self.set_view_preset(ViewPreset::SuperiorAxial, animate);
            } else if ctx.input(|i| i.key_pressed(egui::Key::Num4)) {
                self.set_view_preset(ViewPreset::AnteriorCoronal, animate);
            } else if ctx.input(|i| i.key_pressed(egui::Key::R)) {
                self.reset_camera_to_species(animate);
            } else if ctx.input(|i| i.key_pressed(egui::Key::I)) {
                self.probe_mode = ProbeMode::Inspect;
                self.status_msg = "Probe mode: Inspect Region (I).".to_string();
            } else if ctx.input(|i| i.key_pressed(egui::Key::P)) {
                self.probe_mode = ProbeMode::ShortestPath;
                self.status_msg = "Probe mode: Shortest Signaling Path (P).".to_string();
            } else if ctx.input(|i| i.key_pressed(egui::Key::L)) {
                self.probe_mode = ProbeMode::LesionAblate;
                self.status_msg = "Probe mode: Virtual Lesion / Ablate (L).".to_string();
            } else if ctx.input(|i| i.key_pressed(egui::Key::C)) {
                self.cycle_cutaway();
            }
        }

        // Drag-and-Drop File Ingestion onto Window (.connectome.json, .swc, .trk, .nii)
        let dropped = ctx.input(|i| i.raw.dropped_files.clone());
        for df in dropped {
            if let Some(path) = df.path {
                self.load_file_from_path(&path);
            }
        }

        egui::TopBottomPanel::top("top_menu_bar").show(ctx, |ui| {
            ui.horizontal_wrapped(|ui| {
                ui.heading("ConnectomeAtlas v0.8.0");
                ui.separator();
                let mut switch_to = None;
                for (i, p) in self.presets.iter().enumerate() {
                    let selected = i == self.active_preset_idx;
                    if ui.selectable_label(selected, &p.organism).clicked() && !selected {
                        switch_to = Some(i);
                    }
                }
                if let Some(target_idx) = switch_to {
                    self.switch_preset(target_idx);
                }
                ui.separator();
                ui.selectable_value(&mut self.probe_mode, ProbeMode::Inspect, "Inspect")
                    .on_hover_text("Inspect node properties, receptors, and functional role (Hotkey: I)");
                ui.selectable_value(&mut self.probe_mode, ProbeMode::ShortestPath, "Signaling Path")
                    .on_hover_text("Click source then target node to compute Dijkstra path (Hotkey: P)");
                ui.selectable_value(&mut self.probe_mode, ProbeMode::LesionAblate, "Lesion / Ablate")
                    .on_hover_text("Click region to simulate focal ablation and recalculate efficiency (Hotkey: L)");
                if (!self.atlas.ablated_regions.is_empty() || !self.atlas.severed_tracts.is_empty())
                    && ui.button("Heal All").clicked()
                {
                    self.atlas.heal_all();
                    self.recompute();
                    self.status_msg = "Restored all ablated regions and severed tracts.".to_string();
                }
                ui.separator();
                if ui.button("Jump / Search (Ctrl+K)").clicked() {
                    self.cmd_palette_open = !self.cmd_palette_open;
                    self.cmd_query.clear();
                    self.cmd_selected_idx = 0;
                }
                ui.toggle_value(&mut self.show_left_panel, "Controls");
                ui.toggle_value(&mut self.show_bottom_slices, "2D Slices");
                ui.toggle_value(&mut self.show_right_panel, "Analytics");
                if ui.selectable_label(self.zen_mode, "Zen 3D").clicked() {
                    self.toggle_zen_mode();
                }
            });
        });

        egui::TopBottomPanel::bottom("status_footer").show(ctx, |ui| {
            let delta = &self.lesion_delta;
            ui.horizontal_wrapped(|ui| {
                ui.label(&self.status_msg);
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    ui.label(format!(
                        "Global Eff: {:.3} ({:.1}% intact) | Components: {} | Rich-Club Phi: {:.3}",
                        self.metrics.global_efficiency,
                        delta.efficiency_retention_pct,
                        self.metrics.connected_components,
                        self.metrics.rich_club_coeff_top_quartile
                    ));
                });
            });
            ui.horizontal_wrapped(|ui| {
                ui.small("Copyright CRV Science Network  |  GPL-3.0-or-later  |  ConnectomeAtlas v0.8.0");
                ui.separator();
                if ui.small_button("Shortcuts & About (F1)").clicked() {
                    self.show_about_window = true;
                }
                ui.hyperlink_to("crvscience.com", "https://crvscience.com");
            });
        });

        if self.show_left_panel {
            let scope_override = self.pending_scope_open.take();
            egui::SidePanel::left("left_controls")
                .default_width(295.0)
                .show(ctx, |ui| {
                    egui::ScrollArea::vertical().show(ui, |ui| {
                        ui.label(egui::RichText::new(&self.atlas.coordinate_space).small().weak());
                        ui.add_space(4.0);

                        egui::CollapsingHeader::new("1. Camera & Palette")
                            .default_open(true)
                            .open(scope_override.as_deref().map(|s| s == "explore"))
                            .show(ui, |ui| {
                                ui.horizontal_wrapped(|ui| {
                                    if ui.button("3D Orbit").on_hover_text("Smooth 3D perspective camera (Hotkey: 1)").clicked() {
                                        self.set_view_preset(ViewPreset::Perspective3D, animate);
                                    }
                                    if ui.button("Sagittal").on_hover_text("Left sagittal view (Hotkey: 2)").clicked() {
                                        self.set_view_preset(ViewPreset::LeftSagittal, animate);
                                    }
                                    if ui.button("Axial").on_hover_text("Superior axial top-down view (Hotkey: 3)").clicked() {
                                        self.set_view_preset(ViewPreset::SuperiorAxial, animate);
                                    }
                                    if ui.button("Coronal").on_hover_text("Anterior coronal view (Hotkey: 4)").clicked() {
                                        self.set_view_preset(ViewPreset::AnteriorCoronal, animate);
                                    }
                                    if ui.button("Reset").on_hover_text("Reset camera to species default framing (Hotkey: R)").clicked() {
                                        self.reset_camera_to_species(animate);
                                    }
                                    let fly_lbl = if self.fly_ctrl.active { "Stop Flight" } else { "Fly Path" };
                                    if ui.button(fly_lbl).on_hover_text("Chase camera flying along the active signaling path").clicked() {
                                        self.fly_ctrl.active = !self.fly_ctrl.active;
                                        self.fly_ctrl.progress_u = 0.0;
                                        if self.fly_ctrl.active {
                                            self.camera.cancel_transition();
                                        }
                                    }
                                });
                                ui.radio_value(&mut self.color_mode, ColorSchemeMode::FunctionalNetwork, "Functional Network")
                                    .on_hover_text("Color ROIs and tracts by intrinsic functional module");
                                ui.radio_value(&mut self.color_mode, ColorSchemeMode::DtiDirectionalRgb, "DTI Directional RGB")
                                    .on_hover_text("Color pathways by stereotaxic fiber orientation (R=L-R, G=A-P, B=S-I)");
                                ui.radio_value(&mut self.color_mode, ColorSchemeMode::BetweennessHeatmap, "Betweenness Heatmap")
                                    .on_hover_text("Highlight Brandes shortest-path hub bottlenecks");
                                ui.radio_value(&mut self.color_mode, ColorSchemeMode::HarmonicModePsi, "Harmonic Mode psi_k")
                                    .on_hover_text("Modulate 3D nodes & tracts by active Graph Laplacian eigenmode psi_k");
                                ui.radio_value(&mut self.color_mode, ColorSchemeMode::CriticalitySusceptibility, "Criticality sigma_i")
                                    .on_hover_text("Color nodes & pathways by Kinouchi-Copelli avalanche susceptibility");
                                ui.radio_value(&mut self.color_mode, ColorSchemeMode::TdaBeta1Loops, "TDA beta-1 Loops")
                                    .on_hover_text("Highlight persistent Z2 H1 topological loop cavities and participating hubs");

                                if ui.add(egui::Slider::new(&mut self.min_weight, 0.0..=0.92).text("Min Weight"))
                                    .on_hover_text("Prune pathways below minimum connection weight w_min and recompute topology/TDA")
                                    .changed()
                                {
                                    self.recompute();
                                }
                                ui.add(egui::Slider::new(&mut self.shell_opacity, 0.0..=0.85).text("Shell Alpha"))
                                    .on_hover_text("Adjust translucent Fresnel glass brain/organism outer shell opacity");
                                ui.add(egui::Slider::new(&mut self.flow_speed, 0.0..=1.0).text("River Flow"))
                                    .on_hover_text("Adjust tonic laminar wave speed along pathways (0 pauses animation)");
                            });

                        egui::CollapsingHeader::new("2. 3D Shaders, Cutaway & Slices")
                            .default_open(false)
                            .open(scope_override.as_deref().map(|s| s == "surgical"))
                            .show(ui, |ui| {
                                ui.horizontal(|ui| {
                                    ui.selectable_value(&mut self.tract_style, TractRenderStyle::BraidedFascicles, "Braided Fascicles");
                                    ui.selectable_value(&mut self.tract_style, TractRenderStyle::SingleStreamline, "Single Tube");
                                });
                                ui.horizontal_wrapped(|ui| {
                                    ui.selectable_value(&mut self.cutaway_mode, SurgicalCutawayMode::Intact, "Intact");
                                    ui.selectable_value(&mut self.cutaway_mode, SurgicalCutawayMode::AnteriorSuperiorOctant, "Octant");
                                    ui.selectable_value(&mut self.cutaway_mode, SurgicalCutawayMode::LeftHemisphere, "L-Hemi");
                                    ui.selectable_value(&mut self.cutaway_mode, SurgicalCutawayMode::SuperiorAxialDome, "Dome");
                                });
                                ui.add(egui::Slider::new(&mut self.cutaway_offset_mm, -28.0..=28.0).text("Cutaway Offset"))
                                    .on_hover_text("Shift the CSG surgical shell cutaway plane in stereotaxic mm (Hotkey C cycles modes)");
                                ui.checkbox(&mut self.show_internal_organs, "3D Subcortical & Internal Anatomy");
                                ui.checkbox(&mut self.show_anatomical_labels, "Anatomical Landmark Callouts");
                                ui.checkbox(&mut self.show_tensor_glyphs, "3D Crossing-Fiber Tensor Ellipsoids");
                                ui.checkbox(&mut self.show_slice_planes, "3D Stereotaxic Slice Grid");
                                if self.show_slice_planes {
                                    let ext = self.anatomy_scene.slice_extent_mm.max(40.0);
                                    ui.add(egui::Slider::new(&mut self.slice_x_mm, -ext..=ext).text("Sagittal X"))
                                        .on_hover_text("Move Sagittal YZ plane along Left-Right X axis (mm)");
                                    ui.add(egui::Slider::new(&mut self.slice_y_mm, -ext..=ext).text("Coronal Y"))
                                        .on_hover_text("Move Coronal XZ plane along Posterior-Anterior Y axis (mm)");
                                    ui.add(egui::Slider::new(&mut self.slice_z_mm, -ext..=ext).text("Axial Z"))
                                        .on_hover_text("Move Axial XY plane along Inferior-Superior Z axis (mm)");
                                }
                            });

                        egui::CollapsingHeader::new("3. Receptors, Control & Harmonics")
                            .default_open(false)
                            .open(scope_override.as_deref().map(|s| s == "dynamics"))
                            .show(ui, |ui| {
                                ui.horizontal_wrapped(|ui| {
                                    for (code, label) in [
                                        ("ht", "5-HT"),
                                        ("da", "DA"),
                                        ("ach", "ACh"),
                                        ("glu", "Glu"),
                                        ("gaba", "GABA"),
                                    ] {
                                        if ui
                                            .selectable_value(
                                                &mut self.receptor_channel,
                                                code.to_string(),
                                                label,
                                            )
                                            .changed()
                                        {
                                            self.recompute();
                                        }
                                    }
                                });
                                if ui.add(egui::Slider::new(&mut self.agonist_gain, -0.8..=2.0).text("Receptor Gain alpha"))
                                    .on_hover_text("Modulate regional control input gain B_i = max(0.15, 1 + alpha * rho_i)")
                                    .changed()
                                {
                                    self.recompute();
                                }
                                if ui.add(egui::Slider::new(&mut self.criticality_kappa, 0.40..=1.80).text("Criticality kappa"))
                                    .on_hover_text("Scale structural branching coupling kappa (sigma ~= 1.0 maximizes dynamic range)")
                                    .changed()
                                {
                                    self.recompute();
                                }
                                ui.label(format!(
                                    "Branching sigma: {:.3} ({}) | Delta: {:.1} dB",
                                    self.criticality.branching_ratio_sigma,
                                    self.criticality.regime.label(),
                                    self.criticality.dynamic_range_db
                                ));
                                let max_h = self.control_harmonics.harmonic_modes.len().max(1);
                                ui.add(egui::Slider::new(&mut self.harmonic_mode_idx, 1..=max_h).text("Harmonic psi_k"))
                                    .on_hover_text("Select normalized Graph Laplacian eigenmode psi_k for 3D Harmonic overlay");
                                if ui.add(egui::Slider::new(&mut self.conduction_velocity_m_s, 3.0..=25.0).text("Velocity (m/s)"))
                                    .on_hover_text("Set axonal conduction velocity v (m/s) used for signaling latency and nodal delays")
                                    .changed()
                                {
                                    self.recompute();
                                }
                                ui.label(format!(
                                    "Ctrl Energy dE: {:+.1}% | omega_k: {:.3}",
                                    self.control_harmonics.energy_landscape_delta_pct,
                                    self.control_harmonics.laplacian_wavenumbers.get(self.harmonic_mode_idx.saturating_sub(1)).copied().unwrap_or(0.0)
                                ));
                            });

                        egui::CollapsingHeader::new("4. Virtual Tractotomy Bundles")
                            .default_open(false)
                            .open(scope_override.as_deref().map(|s| s == "surgical"))
                            .show(ui, |ui| {
                                let tracts = self.atlas.unique_tracts();
                                for tract in tracts {
                                    let mut intact = !self.atlas.severed_tracts.contains(&tract);
                                    if ui.checkbox(&mut intact, &tract).changed() {
                                        self.atlas.toggle_tract_sever(&tract);
                                        self.recompute();
                                    }
                                }
                            });
                    });
                });
        }

        if self.show_right_panel {
            egui::SidePanel::right("right_analytics")
                .default_width(375.0)
                .show(ctx, |ui| {
                    ui.horizontal_wrapped(|ui| {
                        ui.selectable_value(&mut self.right_tab, RightTab::Inspector, "Inspector");
                        ui.selectable_value(&mut self.right_tab, RightTab::AdjacencyMatrix, "Matrix");
                        ui.selectable_value(&mut self.right_tab, RightTab::RichClubHubs, "Hubs");
                        ui.selectable_value(&mut self.right_tab, RightTab::CriticalityTda, "Crit / TDA");
                        ui.selectable_value(&mut self.right_tab, RightTab::VaultIo, "Vault / JSON");
                    });
                    ui.separator();

                    egui::ScrollArea::vertical().show(ui, |ui| match self.right_tab {
                        RightTab::Inspector => self.ui_inspector(ui),
                        RightTab::AdjacencyMatrix => self.ui_matrix(ui),
                        RightTab::RichClubHubs => self.ui_hubs(ui),
                        RightTab::CriticalityTda => self.ui_criticality_tda(ui),
                        RightTab::VaultIo => self.ui_vault(ui),
                    });
                });
        }

        if self.show_bottom_slices {
            egui::TopBottomPanel::bottom("triplanar_2d_dock")
                .exact_height(175.0)
                .show(ctx, |ui| {
                    self.draw_triplanar_2d_strip(ui);
                });
        }

        if self.cmd_palette_open {
            let mut open = self.cmd_palette_open;
            let mut chosen_action: Option<CommandAction> = None;
            egui::Window::new("Quick-Jump Command Palette (Ctrl+K)")
                .open(&mut open)
                .collapsible(false)
                .resizable(false)
                .default_width(480.0)
                .anchor(egui::Align2::CENTER_TOP, [0.0, 48.0])
                .show(ctx, |ui| {
                    let resp = ui.text_edit_singleline(&mut self.cmd_query);
                    if resp.changed() {
                        self.cmd_selected_idx = 0;
                    }
                    resp.request_focus();
                    let hits = search_command_index(&self.cmd_entries, &self.cmd_query, 18);
                    if !hits.is_empty() {
                        self.cmd_selected_idx = self.cmd_selected_idx.min(hits.len() - 1);
                        if ctx.input(|i| i.key_pressed(egui::Key::ArrowDown)) {
                            self.cmd_selected_idx = (self.cmd_selected_idx + 1).min(hits.len() - 1);
                        }
                        if ctx.input(|i| i.key_pressed(egui::Key::ArrowUp)) {
                            self.cmd_selected_idx = self.cmd_selected_idx.saturating_sub(1);
                        }
                        if ctx.input(|i| i.key_pressed(egui::Key::Enter)) {
                            chosen_action = Some(hits[self.cmd_selected_idx].action.clone());
                        }
                    } else {
                        self.cmd_selected_idx = 0;
                    }
                    ui.small("Use Up/Down to select, Enter to execute, Esc to close.");
                    ui.separator();
                    egui::ScrollArea::vertical().max_height(260.0).show(ui, |ui| {
                        for (idx, entry) in hits.iter().enumerate() {
                            let is_sel = idx == self.cmd_selected_idx;
                            let lbl = format!("[{}] {} - {}", entry.category, entry.title, entry.subtitle);
                            if ui.selectable_label(is_sel, lbl).clicked() {
                                chosen_action = Some(entry.action.clone());
                            }
                        }
                    });
                });
            self.cmd_palette_open = open;
            if let Some(act) = chosen_action {
                self.cmd_palette_open = false;
                self.execute_command_action(act, animate);
            }
        }

        if self.show_about_window {
            let mut open = self.show_about_window;
            egui::Window::new("Shortcuts & About ConnectomeAtlas (F1)")
                .open(&mut open)
                .collapsible(false)
                .resizable(false)
                .default_width(500.0)
                .anchor(egui::Align2::CENTER_CENTER, Vec2::ZERO)
                .show(ctx, |ui| {
                    ui.heading("ConnectomeAtlas v0.8.0");
                    ui.label("Copyright CRV Science Network");
                    ui.label("License: GNU General Public License v3.0 or later (GPL-3.0-or-later)");
                    ui.separator();
                    ui.label(egui::RichText::new("Interactive Controls & Keyboard Shortcuts").strong());
                    egui::Grid::new("f1_shortcuts_grid")
                        .striped(true)
                        .num_columns(2)
                        .show(ui, |ui| {
                            let bindings = [
                                ("Left-Drag", "Orbit 3D stereotaxic camera"),
                                ("Right / Middle-Drag", "Pan 3D camera target"),
                                ("Mouse Scroll", "Zoom camera distance (65-520 mm)"),
                                ("Click Node (3D/2D/Matrix)", "Run active Probe Mode on region"),
                                ("1 / 2 / 3 / 4", "3D Orbit / Sagittal / Axial / Coronal preset"),
                                ("R", "Reset camera to species default framing"),
                                ("I / P / L", "Probe Mode: Inspect / Path / Lesion"),
                                ("C", "Cycle 3D surgical cutaway mode"),
                                ("Z", "Toggle Zen Full-Canvas 3D mode"),
                                ("Ctrl/Cmd + K", "Command palette (Up/Down, Enter to run)"),
                                ("Ctrl/Cmd + O / S", "Open dataset / Export .connectome.json"),
                                ("F1 / Esc", "Toggle help / Close overlays, flight, or Zen"),
                            ];
                            for (key, desc) in bindings {
                                ui.strong(key);
                                ui.label(desc);
                                ui.end_row();
                            }
                        });
                    ui.separator();
                    ui.small("Standalone safe-Rust workstation featuring Human MNI152, Mouse CCFv3, Drosophila FlyWire FAFB, and C. elegans WormAtlas.");
                    ui.add_space(6.0);
                    ui.horizontal(|ui| {
                        ui.hyperlink_to("CRV Science Network", "https://crvscience.com");
                        ui.separator();
                        ui.hyperlink_to("CRV Science Public Apps", "https://crvscience.com/apps");
                        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                            if ui.button("Close").clicked() {
                                self.show_about_window = false;
                            }
                        });
                    });
                });
            if !open {
                self.show_about_window = false;
            }
        }

        if !animate {
            self.camera.finish_transition();
        }

        egui::CentralPanel::default().show(ctx, |ui| {
            self.draw_3d_viewport(ui);
        });
    }
}

impl eframe::App for ConnectomeAtlasApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        self.ui_step(ctx, true);
    }
}

impl ConnectomeAtlasApp {
    fn ui_inspector(&mut self, ui: &mut egui::Ui) {
        ui.label(egui::RichText::new(&self.atlas.title).strong());
        ui.label(egui::RichText::new(&self.atlas.reference_summary).small());
        ui.label(
            egui::RichText::new(format!(
                "Scale Benchmark: {} | {}",
                self.atlas.benchmark_neurons_or_parcels, self.atlas.benchmark_synapses_or_tracts
            ))
            .small()
            .color(Color32::from_rgb(125, 211, 252)),
        );
        ui.separator();

        if let Some(path) = &self.active_path {
            ui.group(|ui| {
                ui.label(egui::RichText::new("Shortest Signaling Pathway").strong());
                ui.label(path.hops.join(" -> "));
                ui.label(format!(
                    "Hops: {} | Latency: {:.2} ms ({:.1} mm @ {:.1} m/s) | Cost: {:.3} | Bottleneck W: {:.2}",
                    path.hops.len().saturating_sub(1),
                    path.conduction_latency_ms,
                    path.arc_length_mm,
                    self.conduction_velocity_m_s,
                    path.total_distance_cost,
                    path.bottleneck_weight
                ));
            });
            ui.add_space(6.0);
        }

        let mut next_region: Option<String> = None;
        if let Some(sel_id) = &self.selected_region.clone() {
            if let Some(reg) = self.atlas.regions.iter().find(|r| &r.id == sel_id) {
                ui.heading(format!("{} - {}", reg.id, reg.name));
                ui.label(format!("System / Network: {}", reg.network));
                ui.label(format!(
                    "Stereotaxic Coord: ({:.1}, {:.1}, {:.1})",
                    reg.coord[0], reg.coord[1], reg.coord[2]
                ));
                ui.group(|ui| {
                    ui.label(egui::RichText::new("What It Does (Functional Role)").strong().color(Color32::from_rgb(125, 211, 252)));
                    ui.label(&reg.what_it_does);
                });
                let (inc, out, coup) = self.atlas.incident_io_breakdown(sel_id, self.min_weight);
                ui.add_space(4.0);
                ui.label(egui::RichText::new(format!("<- Afferent Inputs ({})", inc.len())).strong().color(Color32::from_rgb(52, 211, 153)));
                for e in inc {
                    if ui.link(format!("<- {} (w={:.2}, {})", e.source, e.weight, e.tract_name)).clicked() {
                        next_region = Some(e.source.clone());
                    }
                }
                ui.label(egui::RichText::new(format!("-> Efferent Outputs ({})", out.len())).strong().color(Color32::from_rgb(251, 146, 60)));
                for e in out {
                    if ui.link(format!("-> {} (w={:.2}, {})", e.target, e.weight, e.tract_name)).clicked() {
                        next_region = Some(e.target.clone());
                    }
                }
                ui.label(egui::RichText::new(format!("<-> Reciprocal / Coupled ({})", coup.len())).strong().color(Color32::from_rgb(56, 189, 248)));
                for e in coup {
                    let peer = if &e.source == sel_id { &e.target } else { &e.source };
                    if ui.link(format!("<-> {} (w={:.2}, {})", peer, e.weight, e.tract_name)).clicked() {
                        next_region = Some(peer.clone());
                    }
                }

                if let Some(nm) = self
                    .metrics
                    .node_metrics
                    .iter()
                    .find(|m| &m.region_id == sel_id)
                {
                    ui.separator();
                    egui::Grid::new("node_telemetry_grid").striped(true).show(ui, |ui| {
                        ui.label("Degree / Strength:");
                        ui.label(format!("{} / {:.2}", nm.degree_total, nm.strength_total));
                        ui.end_row();
                        ui.label("Participation Coeff (P):");
                        ui.label(format!("{:.3}", nm.participation_coeff));
                        ui.end_row();
                        ui.label("Clustering / Betweenness:");
                        ui.label(format!(
                            "{:.3} / {:.3}",
                            nm.clustering_coeff, nm.betweenness_centrality
                        ));
                        ui.end_row();
                        ui.label("Communicability / Delay:");
                        ui.label(format!("{:.3} / {:.1} ms", nm.communicability_centrality, nm.mean_conduction_delay_ms));
                        ui.end_row();
                        ui.label("Receptors (DA/5HT/ACh/Glu/GABA):");
                        ui.label(format!(
                            "{:.2}/{:.2}/{:.2}/{:.2}/{:.2}",
                            reg.receptors.da, reg.receptors.serotonin, reg.receptors.ach, reg.receptors.glutamate, reg.receptors.gaba
                        ));
                        ui.end_row();
                        if let Some(cp) = self.control_harmonics.node_profiles.iter().find(|p| &p.region_id == sel_id) {
                            ui.label("Avg / Modal Control:");
                            ui.label(format!("{:.3} / {:.3}", cp.average_controllability, cp.modal_controllability));
                            ui.end_row();
                            ui.label("Gradients G1 / G2:");
                            ui.label(format!("{:+.3} / {:+.3}", cp.principal_gradient_g1, cp.secondary_gradient_g2));
                            ui.end_row();
                        }
                        if let Some(ft) = self.fiber_tensors.iter().find(|t| &t.region_id == sel_id) {
                            ui.label("Fiber Tensor (FA / cl / cp):");
                            ui.label(format!(
                                "{:.2} / {:.2} / {:.2}",
                                ft.fractional_anisotropy, ft.westin_linear_cl, ft.westin_planar_cp
                            ));
                            ui.end_row();
                        }
                        if let Some(sp) = self.systems_neuro.node_profiles.iter().find(|p| &p.region_id == sel_id) {
                            ui.label("Synergy Rank / AG:");
                            ui.label(format!("{:+.2} / {:.2}", sp.synergy_rank_gradient, sp.aerobic_glycolysis_index));
                            ui.end_row();
                            ui.label("mRNA (SST / PVALB):");
                            ui.label(format!("{:.2} / {:.2}", sp.transcriptomics.sst, sp.transcriptomics.pvalb));
                            ui.end_row();
                        }
                        if let Some(cp) = self.criticality.node_profiles.iter().find(|p| &p.region_id == sel_id) {
                            ui.label("Criticality sigma_i / Delta Range:");
                            ui.label(format!("{:.2} / {:.1} dB", cp.avalanche_susceptibility, self.criticality.dynamic_range_db));
                            ui.end_row();
                        }
                        if let Some(tp) = self.topology_tda.node_profiles.iter().find(|p| &p.region_id == sel_id) {
                            ui.label("Cliques (3/4) / beta-1 Loops:");
                            ui.label(format!("{} / {} / {}", tp.clique_3_count, tp.clique_4_count, tp.beta1_loop_participation));
                            ui.end_row();
                        }
                        ui.label("Rich-Club Hub:");
                        ui.label(if nm.is_rich_club_hub { "Yes (Top Quartile)" } else { "No" });
                        ui.end_row();
                    });
                    ui.separator();
                    ui.label(egui::RichText::new("Functional Domain / Lesion Decoder").strong());
                    for d in &self.domain_decode {
                        if d.is_lesion_disruption {
                            ui.label(format!("{}: -{:.0}% disruption", d.domain_label, d.score * 100.0));
                        } else {
                            ui.label(format!("{}: {:.2}", d.domain_label, d.score));
                        }
                    }
                }
            }
        }
        if let Some(rid) = next_region {
            self.select_region(rid, true);
        }
    }

    fn ui_matrix(&mut self, ui: &mut egui::Ui) {
        let n = self.atlas.regions.len();
        ui.label(format!("Community-Ordered Connectivity Matrix ({n} x {n})"));
        ui.label(
            egui::RichText::new(format!(
                "TDA Cliques: {} 3-cliq, {} 4-cliq | Persistent beta-1 Cavities: {} (max pi={:.2})",
                self.topology_tda.clique_3_total,
                self.topology_tda.clique_4_total,
                self.topology_tda.persistent_beta1_count,
                self.topology_tda.max_persistence
            ))
            .small()
            .color(Color32::from_rgb(34, 211, 238)),
        );
        if n == 0 {
            return;
        }

        // Group region indices contiguously by functional community network
        let mut net_order: Vec<String> = Vec::new();
        for r in &self.atlas.regions {
            if !net_order.iter().any(|net| net == &r.network) {
                net_order.push(r.network.clone());
            }
        }
        let mut sorted_order: Vec<usize> = Vec::with_capacity(n);
        struct CommBlock {
            network: String,
            start: usize,
            end: usize,
            rgb: [u8; 3],
        }
        let mut blocks: Vec<CommBlock> = Vec::with_capacity(net_order.len());
        for net in &net_order {
            let start = sorted_order.len();
            let mut rgb = [56_u8, 189, 248];
            for (idx, r) in self.atlas.regions.iter().enumerate() {
                if &r.network == net {
                    if sorted_order.len() == start {
                        rgb = parse_hex_color(&r.color_hex);
                    }
                    sorted_order.push(idx);
                }
            }
            let end = sorted_order.len();
            if end > start {
                blocks.push(CommBlock {
                    network: net.clone(),
                    start,
                    end,
                    rgb,
                });
            }
        }

        let side = ui.available_width().min(330.0);
        let strip_px = 6.0_f32;
        let mat_side = (side - strip_px).max(120.0);
        let (resp, painter) = ui.allocate_painter(Vec2::splat(side), Sense::click());
        let outer = resp.rect;
        let rect = Rect::from_min_size(
            Pos2::new(outer.left() + strip_px, outer.top() + strip_px),
            Vec2::splat(mat_side),
        );
        painter.rect_filled(outer, 4.0, Color32::from_rgb(10, 14, 22));

        let cell = mat_side / (n as f32);
        let mut idx_map = HashMap::with_capacity(n);
        for (i, r) in self.atlas.regions.iter().enumerate() {
            idx_map.insert(r.id.as_str(), i);
        }
        let mut weights = vec![vec![0.0_f64; n]; n];
        for e in &self.atlas.edges {
            if !self.atlas.is_edge_intact(e, self.min_weight) {
                continue;
            }
            if let (Some(&i), Some(&j)) = (idx_map.get(e.source.as_str()), idx_map.get(e.target.as_str())) {
                weights[i][j] = weights[i][j].max(e.weight);
                if !e.directed {
                    weights[j][i] = weights[j][i].max(e.weight);
                }
            }
        }

        // Draw sorted matrix cells with perceptual thermal-cyan ramp
        for (row_pos, &src_idx) in sorted_order.iter().enumerate() {
            for (col_pos, &dst_idx) in sorted_order.iter().enumerate() {
                let w = weights[src_idx][dst_idx];
                if w > 0.0 {
                    let rgb = matrix_weight_rgb(w as f32);
                    let cell_rect = Rect::from_min_size(
                        Pos2::new(
                            rect.left() + (col_pos as f32) * cell,
                            rect.top() + (row_pos as f32) * cell,
                        ),
                        Vec2::splat((cell - 0.45).max(1.5)),
                    );
                    painter.rect_filled(
                        cell_rect,
                        0.0,
                        Color32::from_rgb(rgb[0], rgb[1], rgb[2]),
                    );
                }
            }
        }

        // Delineate community side color strips, grid dividers, and intra-network diagonal blocks
        for blk in &blocks {
            let bcol = Color32::from_rgb(blk.rgb[0], blk.rgb[1], blk.rgb[2]);
            let x0 = rect.left() + (blk.start as f32) * cell;
            let x1 = rect.left() + (blk.end as f32) * cell;
            let y0 = rect.top() + (blk.start as f32) * cell;
            let y1 = rect.top() + (blk.end as f32) * cell;

            painter.rect_filled(
                Rect::from_min_max(Pos2::new(x0, outer.top()), Pos2::new(x1, rect.top() - 1.0)),
                0.0,
                bcol,
            );
            painter.rect_filled(
                Rect::from_min_max(Pos2::new(outer.left(), y0), Pos2::new(rect.left() - 1.0, y1)),
                0.0,
                bcol,
            );
            painter.rect_stroke(
                Rect::from_min_max(Pos2::new(x0, y0), Pos2::new(x1, y1)),
                0.0,
                Stroke::new(1.25_f32, Color32::from_rgba_unmultiplied(blk.rgb[0], blk.rgb[1], blk.rgb[2], 210)),
                egui::StrokeKind::Inside,
            );
            if blk.start > 0 {
                let div_col = Color32::from_rgba_unmultiplied(148, 163, 184, 52);
                painter.line_segment(
                    [Pos2::new(x0, rect.top()), Pos2::new(x0, rect.bottom())],
                    Stroke::new(0.8_f32, div_col),
                );
                painter.line_segment(
                    [Pos2::new(rect.left(), y0), Pos2::new(rect.right(), y0)],
                    Stroke::new(0.8_f32, div_col),
                );
            }
        }

        let cell_at = |pos: Pos2| -> Option<(usize, usize)> {
            if !rect.contains(pos) {
                return None;
            }
            let col = (((pos.x - rect.left()) / cell).floor() as isize).clamp(0, n as isize - 1) as usize;
            let row = (((pos.y - rect.top()) / cell).floor() as isize).clamp(0, n as isize - 1) as usize;
            Some((row, col))
        };
        let mut clicked_rid: Option<String> = None;
        if resp.clicked() {
            if let Some(pos) = resp.interact_pointer_pos() {
                if let Some((row_pos, _col_pos)) = cell_at(pos) {
                    clicked_rid = Some(self.atlas.regions[sorted_order[row_pos]].id.clone());
                }
            }
        }
        if let Some(pos) = resp.hover_pos() {
            if let Some((row_pos, col_pos)) = cell_at(pos) {
                let hover_rect = Rect::from_min_size(
                    Pos2::new(rect.left() + (col_pos as f32) * cell, rect.top() + (row_pos as f32) * cell),
                    Vec2::splat(cell.max(2.0)),
                );
                painter.rect_stroke(
                    hover_rect,
                    0.0,
                    Stroke::new(1.0_f32, Color32::WHITE),
                    egui::StrokeKind::Inside,
                );
                let src_idx = sorted_order[row_pos];
                let dst_idx = sorted_order[col_pos];
                let src_reg = &self.atlas.regions[src_idx];
                let dst_reg = &self.atlas.regions[dst_idx];
                let tip = if src_idx == dst_idx {
                    format!("{} — {} [{}]", src_reg.id, src_reg.name, src_reg.network)
                } else {
                    format!(
                        "{} [{}] -> {} [{}] | w={:.2}",
                        src_reg.id, src_reg.network, dst_reg.id, dst_reg.network, weights[src_idx][dst_idx]
                    )
                };
                resp.on_hover_text(tip);
            }
        }
        if let Some(rid) = clicked_rid {
            self.select_region(rid.clone(), true);
            self.status_msg = format!("Selected region {rid} from matrix.");
        }

        // Calibrated colorbar legend + community legend chips
        ui.add_space(6.0);
        let bar_w = side;
        let (bar_resp, bar_painter) = ui.allocate_painter(Vec2::new(bar_w, 26.0), Sense::hover());
        let br = bar_resp.rect;
        let grad_rect = Rect::from_min_size(br.min, Vec2::new(bar_w, 10.0));
        let steps = 40_usize;
        let seg_w = grad_rect.width() / (steps as f32);
        for s in 0..steps {
            let t = (s as f32 + 0.5) / (steps as f32);
            let rgb = matrix_weight_rgb(t);
            let seg = Rect::from_min_size(
                Pos2::new(grad_rect.left() + (s as f32) * seg_w, grad_rect.top()),
                Vec2::new(seg_w + 0.5, grad_rect.height()),
            );
            bar_painter.rect_filled(seg, 0.0, Color32::from_rgb(rgb[0], rgb[1], rgb[2]));
        }
        bar_painter.rect_stroke(
            grad_rect,
            0.0,
            Stroke::new(0.8_f32, Color32::from_rgb(71, 85, 105)),
            egui::StrokeKind::Inside,
        );
        let mw_x = grad_rect.left() + (self.min_weight as f32).clamp(0.0, 1.0) * grad_rect.width();
        bar_painter.line_segment(
            [Pos2::new(mw_x, grad_rect.top() - 1.0), Pos2::new(mw_x, grad_rect.bottom() + 2.0)],
            Stroke::new(1.5_f32, Color32::WHITE),
        );
        let lbl_col = Color32::from_rgb(148, 163, 184);
        bar_painter.text(
            Pos2::new(grad_rect.left(), grad_rect.bottom() + 2.0),
            egui::Align2::LEFT_TOP,
            "0.0",
            FontId::proportional(9.5),
            lbl_col,
        );
        bar_painter.text(
            Pos2::new(grad_rect.center().x, grad_rect.bottom() + 2.0),
            egui::Align2::CENTER_TOP,
            format!("Weight w (min={:.2})", self.min_weight),
            FontId::proportional(9.5),
            lbl_col,
        );
        bar_painter.text(
            Pos2::new(grad_rect.right(), grad_rect.bottom() + 2.0),
            egui::Align2::RIGHT_TOP,
            "1.0",
            FontId::proportional(9.5),
            lbl_col,
        );

        ui.horizontal_wrapped(|ui| {
            for blk in &blocks {
                let col = Color32::from_rgb(blk.rgb[0], blk.rgb[1], blk.rgb[2]);
                ui.colored_label(col, format!("■ {}", blk.network));
            }
        });
    }

    fn ui_criticality_tda(&mut self, ui: &mut egui::Ui) {
        let crit = &self.criticality;
        ui.label(
            egui::RichText::new(format!(
                "Kinouchi-Copelli Criticality | kappa={:.2} | sigma={:.3} ({}) | Delta={:.2} dB",
                crit.coupling_kappa,
                crit.branching_ratio_sigma,
                crit.regime.label(),
                crit.dynamic_range_db
            ))
            .small()
            .strong()
            .color(Color32::from_rgb(125, 211, 252)),
        );

        // 1) Response Curve F(r)
        let w = ui.available_width().min(345.0);
        let (rc_resp, rc_p) = ui.allocate_painter(Vec2::new(w, 130.0), Sense::hover());
        let rcr = rc_resp.rect;
        rc_p.rect_filled(rcr, 4.0, Color32::from_rgb(10, 14, 22));
        let pad_l = 32.0_f32;
        let pad_r = 10.0_f32;
        let pad_t = 16.0_f32;
        let pad_b = 20.0_f32;
        let plot_rc = Rect::from_min_max(
            Pos2::new(rcr.left() + pad_l, rcr.top() + pad_t),
            Pos2::new(rcr.right() - pad_r, rcr.bottom() - pad_b),
        );
        rc_p.rect_stroke(
            plot_rc,
            0.0,
            Stroke::new(0.8_f32, Color32::from_rgb(51, 65, 85)),
            egui::StrokeKind::Inside,
        );
        rc_p.text(
            Pos2::new(rcr.left() + 6.0, rcr.top() + 3.0),
            egui::Align2::LEFT_TOP,
            format!("Response Curve F(r)  |  F0={:.3}, Fmax={:.3}", crit.spontaneous_floor_f0, crit.saturation_ceiling_fmax),
            FontId::proportional(10.0),
            Color32::from_rgb(203, 213, 225),
        );
        if crit.response_curve.len() >= 2 {
            let log_min = crit.response_curve.first().map(|p| p.log10_stimulus_r as f32).unwrap_or(-5.2);
            let log_max = crit.response_curve.last().map(|p| p.log10_stimulus_r as f32).unwrap_or(0.6);
            let log_span = (log_max - log_min).max(1e-3);
            let y_max = (crit.saturation_ceiling_fmax as f32 * 1.08).max(0.05);
            let map_rc = |lr: f32, f: f32| -> Pos2 {
                let tx = ((lr - log_min) / log_span).clamp(0.0, 1.0);
                let ty = (f / y_max).clamp(0.0, 1.0);
                Pos2::new(
                    plot_rc.left() + tx * plot_rc.width(),
                    plot_rc.bottom() - ty * plot_rc.height(),
                )
            };

            // Shade [r10, r90] dynamic range window
            let lr10 = (crit.r10_threshold.max(1e-9).log10() as f32).clamp(log_min, log_max);
            let lr90 = (crit.r90_threshold.max(1e-9).log10() as f32).clamp(log_min, log_max);
            let x10 = map_rc(lr10, 0.0).x;
            let x90 = map_rc(lr90, 0.0).x;
            rc_p.rect_filled(
                Rect::from_min_max(Pos2::new(x10, plot_rc.top()), Pos2::new(x90, plot_rc.bottom())),
                0.0,
                Color32::from_rgba_unmultiplied(56, 189, 248, 28),
            );
            let f0_y = map_rc(log_min, crit.spontaneous_floor_f0 as f32).y;
            rc_p.line_segment(
                [Pos2::new(plot_rc.left(), f0_y), Pos2::new(plot_rc.right(), f0_y)],
                Stroke::new(1.0_f32, Color32::from_rgba_unmultiplied(251, 146, 60, 140)),
            );

            let pts: Vec<Pos2> = crit
                .response_curve
                .iter()
                .map(|p| map_rc(p.log10_stimulus_r as f32, p.mean_firing_response as f32))
                .collect();
            rc_p.add(egui::Shape::line(pts.clone(), Stroke::new(2.0_f32, Color32::from_rgb(56, 189, 248))));
            for pt in &pts {
                rc_p.circle_filled(*pt, 2.2, Color32::from_rgb(125, 211, 252));
            }
            if let Some(hp) = rc_resp.hover_pos() {
                if plot_rc.contains(hp) {
                    if let Some(best) = crit
                        .response_curve
                        .iter()
                        .min_by(|a, b| {
                            let ax = map_rc(a.log10_stimulus_r as f32, a.mean_firing_response as f32).x;
                            let bx = map_rc(b.log10_stimulus_r as f32, b.mean_firing_response as f32).x;
                            (ax - hp.x).abs().total_cmp(&(bx - hp.x).abs())
                        })
                    {
                        let sp = map_rc(best.log10_stimulus_r as f32, best.mean_firing_response as f32);
                        rc_p.circle_filled(sp, 4.0, Color32::WHITE);
                        rc_resp.on_hover_text(format!(
                            "log10(r)={:.2} (r={:.2e}) | F(r)={:.4}",
                            best.log10_stimulus_r, best.stimulus_r, best.mean_firing_response
                        ));
                    }
                }
            }
        }
        rc_p.text(
            Pos2::new(plot_rc.center().x, rcr.bottom() - 2.0),
            egui::Align2::CENTER_BOTTOM,
            "log10(stimulus r)  [shaded: r10..r90 Delta]",
            FontId::proportional(9.0),
            Color32::from_rgb(148, 163, 184),
        );

        // 2) Coupling Sweep kappa -> Dynamic Range Delta (dB)
        ui.add_space(6.0);
        let (ks_resp, ks_p) = ui.allocate_painter(Vec2::new(w, 125.0), Sense::click());
        let ksr = ks_resp.rect;
        ks_p.rect_filled(ksr, 4.0, Color32::from_rgb(10, 14, 22));
        let plot_ks = Rect::from_min_max(
            Pos2::new(ksr.left() + pad_l, ksr.top() + pad_t),
            Pos2::new(ksr.right() - pad_r, ksr.bottom() - pad_b),
        );
        ks_p.rect_stroke(
            plot_ks,
            0.0,
            Stroke::new(0.8_f32, Color32::from_rgb(51, 65, 85)),
            egui::StrokeKind::Inside,
        );
        ks_p.text(
            Pos2::new(ksr.left() + 6.0, ksr.top() + 3.0),
            egui::Align2::LEFT_TOP,
            "Coupling Sweep Delta(kappa) [click to tune kappa]",
            FontId::proportional(10.0),
            Color32::from_rgb(203, 213, 225),
        );
        let mut clicked_kappa: Option<f64> = None;
        if crit.kappa_sweep.len() >= 2 {
            let k_min = crit.kappa_sweep.first().map(|p| p.coupling_kappa as f32).unwrap_or(0.45);
            let k_max = crit.kappa_sweep.last().map(|p| p.coupling_kappa as f32).unwrap_or(1.72);
            let k_span = (k_max - k_min).max(1e-3);
            let db_max = crit
                .kappa_sweep
                .iter()
                .map(|p| p.dynamic_range_db as f32)
                .fold(10.0_f32, f32::max)
                * 1.12;
            let map_ks = |kv: f32, db: f32| -> Pos2 {
                let tx = ((kv - k_min) / k_span).clamp(0.0, 1.0);
                let ty = (db / db_max).clamp(0.0, 1.0);
                Pos2::new(
                    plot_ks.left() + tx * plot_ks.width(),
                    plot_ks.bottom() - ty * plot_ks.height(),
                )
            };
            let x_crit = map_ks(1.0, 0.0).x;
            ks_p.line_segment(
                [Pos2::new(x_crit, plot_ks.top()), Pos2::new(x_crit, plot_ks.bottom())],
                Stroke::new(1.0_f32, Color32::from_rgba_unmultiplied(52, 211, 153, 140)),
            );
            let pts: Vec<Pos2> = crit
                .kappa_sweep
                .iter()
                .map(|p| map_ks(p.coupling_kappa as f32, p.dynamic_range_db as f32))
                .collect();
            ks_p.add(egui::Shape::line(pts.clone(), Stroke::new(2.0_f32, Color32::from_rgb(250, 204, 21))));
            for pt in &pts {
                ks_p.circle_filled(*pt, 2.4, Color32::from_rgb(253, 224, 71));
            }
            let cur_pt = map_ks(crit.coupling_kappa as f32, crit.dynamic_range_db as f32);
            ks_p.circle_filled(cur_pt, 4.5, Color32::from_rgb(52, 211, 153));
            ks_p.circle_stroke(cur_pt, 5.5, Stroke::new(1.2_f32, Color32::WHITE));

            if ks_resp.clicked() {
                if let Some(cp) = ks_resp.interact_pointer_pos() {
                    if plot_ks.contains(cp) {
                        let frac = ((cp.x - plot_ks.left()) / plot_ks.width()).clamp(0.0, 1.0);
                        clicked_kappa = Some((k_min + frac * k_span) as f64);
                    }
                }
            }
            if let Some(hp) = ks_resp.hover_pos() {
                if plot_ks.contains(hp) {
                    if let Some(best) = crit
                        .kappa_sweep
                        .iter()
                        .min_by(|a, b| {
                            let ax = map_ks(a.coupling_kappa as f32, a.dynamic_range_db as f32).x;
                            let bx = map_ks(b.coupling_kappa as f32, b.dynamic_range_db as f32).x;
                            (ax - hp.x).abs().total_cmp(&(bx - hp.x).abs())
                        })
                    {
                        ks_resp.on_hover_text(format!(
                            "kappa={:.2} | sigma={:.3} | Delta={:.2} dB",
                            best.coupling_kappa, best.branching_ratio_sigma, best.dynamic_range_db
                        ));
                    }
                }
            }
        }
        ks_p.text(
            Pos2::new(plot_ks.center().x, ksr.bottom() - 2.0),
            egui::Align2::CENTER_BOTTOM,
            "Coupling kappa (green line: critical kappa = 1.00)",
            FontId::proportional(9.0),
            Color32::from_rgb(148, 163, 184),
        );
        if let Some(new_k) = clicked_kappa {
            self.criticality_kappa = new_k.clamp(0.40, 1.80);
            self.recompute();
            self.status_msg = format!(
                "Tuned kappa = {:.2} (sigma = {:.3}).",
                self.criticality_kappa, self.criticality.branching_ratio_sigma
            );
        }

        // 3) TDA Persistent Homology beta-1 Persistence Barcode
        ui.add_space(8.0);
        let tda = &self.topology_tda;
        ui.label(
            egui::RichText::new(format!(
                "TDA Z2 Persistent Homology | 3-Cliq: {} | 4-Cliq: {} | beta-1: {} ({} essential)",
                tda.clique_3_total,
                tda.clique_4_total,
                tda.persistent_beta1_count,
                tda.essential_beta1_count
            ))
            .small()
            .strong()
            .color(Color32::from_rgb(45, 212, 191)),
        );
        let n_bars = tda.cavities.len().min(14);
        let bar_h = (n_bars.max(3) as f32 * 14.0 + 34.0).clamp(76.0, 236.0);
        let (bc_resp, bc_p) = ui.allocate_painter(Vec2::new(w, bar_h), Sense::click());
        let bcr = bc_resp.rect;
        bc_p.rect_filled(bcr, 4.0, Color32::from_rgb(10, 14, 22));
        let plot_bc = Rect::from_min_max(
            Pos2::new(bcr.left() + 36.0, bcr.top() + 16.0),
            Pos2::new(bcr.right() - 12.0, bcr.bottom() - 18.0),
        );
        bc_p.rect_stroke(
            plot_bc,
            0.0,
            Stroke::new(0.8_f32, Color32::from_rgb(51, 65, 85)),
            egui::StrokeKind::Inside,
        );
        bc_p.text(
            Pos2::new(bcr.left() + 6.0, bcr.top() + 3.0),
            egui::Align2::LEFT_TOP,
            "H1 (beta-1) Persistence Barcode [w_birth -> w_death]",
            FontId::proportional(10.0),
            Color32::from_rgb(203, 213, 225),
        );
        let mut clicked_cavity_node: Option<String> = None;
        if n_bars == 0 {
            bc_p.text(
                plot_bc.center(),
                egui::Align2::CENTER_CENTER,
                "No persistent beta-1 cavities at current min_weight.",
                FontId::proportional(10.5),
                Color32::from_rgb(148, 163, 184),
            );
        } else {
            let row_h = plot_bc.height() / (n_bars as f32);
            let map_w = |weight: f64| -> f32 {
                // Descending weight filtration: birth at high weight (left) -> death at low weight (right)
                let t = (1.0 - weight.clamp(0.0, 1.0)) as f32;
                plot_bc.left() + t * plot_bc.width()
            };
            let mut hovered_tip: Option<String> = None;
            let hover_pos = bc_resp.hover_pos();
            let click_pos = if bc_resp.clicked() { bc_resp.interact_pointer_pos() } else { None };
            for (idx, cav) in tda.cavities.iter().take(n_bars).enumerate() {
                let y_mid = plot_bc.top() + (idx as f32 + 0.5) * row_h;
                let x_start = map_w(cav.birth_weight);
                let x_end = map_w(cav.death_weight).max(x_start + 3.0);
                let col = if cav.is_essential {
                    Color32::from_rgb(250, 204, 21)
                } else {
                    Color32::from_rgb(34, 211, 238)
                };
                bc_p.text(
                    Pos2::new(plot_bc.left() - 4.0, y_mid),
                    egui::Align2::RIGHT_CENTER,
                    format!("g{}", cav.cavity_index),
                    FontId::proportional(9.5),
                    col,
                );
                let bar_rect = Rect::from_min_max(
                    Pos2::new(x_start, y_mid - (row_h * 0.28).min(4.0)),
                    Pos2::new(x_end, y_mid + (row_h * 0.28).min(4.0)),
                );
                bc_p.rect_filled(bar_rect, 1.5, col);
                if cav.is_essential {
                    bc_p.circle_filled(Pos2::new(x_end, y_mid), 2.5, Color32::from_rgb(52, 211, 153));
                }
                let row_band = Rect::from_min_max(
                    Pos2::new(bcr.left(), plot_bc.top() + (idx as f32) * row_h),
                    Pos2::new(bcr.right(), plot_bc.top() + ((idx + 1) as f32) * row_h),
                );
                if let Some(hp) = hover_pos {
                    if row_band.contains(hp) {
                        bc_p.rect_stroke(
                            bar_rect.expand(1.5),
                            1.5,
                            Stroke::new(1.0_f32, Color32::WHITE),
                            egui::StrokeKind::Outside,
                        );
                        hovered_tip = Some(format!(
                            "gamma_{} [{:.2} -> {:.2}{}] | pi={:.3} ({} nodes)\n{}",
                            cav.cavity_index,
                            cav.birth_weight,
                            cav.death_weight,
                            if cav.is_essential { " inf" } else { "" },
                            cav.persistence,
                            cav.cycle_nodes.len(),
                            cav.cycle_nodes.join(" -> ")
                        ));
                    }
                }
                if let Some(cp) = click_pos {
                    if row_band.contains(cp) {
                        clicked_cavity_node = cav.cycle_nodes.first().cloned();
                    }
                }
            }
            if let Some(tip) = hovered_tip {
                bc_resp.on_hover_text(tip);
            }
        }
        bc_p.text(
            Pos2::new(plot_bc.left(), bcr.bottom() - 2.0),
            egui::Align2::LEFT_BOTTOM,
            "w=1.0",
            FontId::proportional(9.0),
            Color32::from_rgb(148, 163, 184),
        );
        bc_p.text(
            Pos2::new(plot_bc.center().x, bcr.bottom() - 2.0),
            egui::Align2::CENTER_BOTTOM,
            "Descending Filtration Weight (gold=essential)",
            FontId::proportional(9.0),
            Color32::from_rgb(148, 163, 184),
        );
        bc_p.text(
            Pos2::new(plot_bc.right(), bcr.bottom() - 2.0),
            egui::Align2::RIGHT_BOTTOM,
            "w=0.0",
            FontId::proportional(9.0),
            Color32::from_rgb(148, 163, 184),
        );
        if let Some(rid) = clicked_cavity_node {
            self.select_region(rid.clone(), true);
            self.color_mode = ColorSchemeMode::TdaBeta1Loops;
            self.status_msg = format!("Selected TDA loop participant {rid}.");
        }
    }

    fn ui_hubs(&mut self, ui: &mut egui::Ui) {
        let nb = &self.systems_neuro.null_benchmark;
        let (mod_q, mod_comms) = self.atlas.compute_modularity_q(self.min_weight);
        ui.group(|ui| {
            ui.label(egui::RichText::new(format!(
                "Economy E/C: {:.2} | Small-World sigma: {:.2}x | Rich-Club Phi_norm: {:.2}x | Modularity Q: {:.3} ({} net)",
                self.systems_neuro.economy_efficiency_ratio, nb.small_world_sigma, nb.rich_club_normalized, mod_q, mod_comms
            )).small());
            ui.label(egui::RichText::new(format!(
                "Wiring Cost vs Null: {:+.1}% | Generative (eta*={:.1}, gamma*={:.1}) | |r(mRNA,G1)|={:.2}",
                nb.wiring_cost_saving_pct, nb.generative_eta_distance, nb.generative_gamma_homophily, self.systems_neuro.transcriptomic_g1_pearson_r
            )).small().color(Color32::from_rgb(45, 212, 191)));
        });
        ui.add_space(4.0);
        egui::CollapsingHeader::new("Comparative Connectome Graph-Topology Digest")
            .default_open(true)
            .show(ui, |ui| {
                ui.horizontal_wrapped(|ui| {
                    if ui.small_button("Export Digest JSON...").clicked() {
                        self.export_topology_digest_dialog(false);
                    }
                    if ui.small_button("Export Digest CSV...").clicked() {
                        self.export_topology_digest_dialog(true);
                    }
                });
                egui::Grid::new("comparative_topology_digest_grid").striped(true).show(ui, |ui| {
                    ui.strong("Preset");
                    ui.strong("Sigma");
                    ui.strong("Clust C");
                    ui.strong("Path L");
                    ui.strong("RC Phi");
                    ui.strong("Mod Q");
                    ui.end_row();
                    for e in &self.topology_digest.entries {
                        let is_cur = e.preset_id == self.atlas.id;
                        let lbl = if is_cur {
                            egui::RichText::new(&e.preset_id).strong().color(Color32::from_rgb(56, 189, 248))
                        } else {
                            egui::RichText::new(&e.preset_id)
                        };
                        ui.label(lbl);
                        ui.label(format!("{:.2}", e.small_world_sigma));
                        ui.label(format!("{:.3}", e.mean_clustering_coeff));
                        ui.label(format!("{:.2}", e.characteristic_path_length));
                        ui.label(format!("{:.3}", e.rich_club_coeff));
                        ui.label(format!("{:.3}", e.modularity_q));
                        ui.end_row();
                    }
                });
            });
        ui.add_space(4.0);
        let mut sorted = self.metrics.node_metrics.clone();
        sorted.sort_by(|a, b| b.strength_total.total_cmp(&a.strength_total));
        egui::Grid::new("hub_ranking_grid").striped(true).show(ui, |ui| {
            ui.strong("Region");
            ui.strong("Deg");
            ui.strong("Str");
            ui.strong("Part P");
            ui.strong("Hub");
            ui.end_row();
            for item in sorted {
                if ui.link(&item.region_id).clicked() {
                    self.select_region(item.region_id.clone(), true);
                }
                ui.label(format!("{}", item.degree_total));
                ui.label(format!("{:.2}", item.strength_total));
                ui.label(format!("{:.2}", item.participation_coeff));
                ui.label(if item.is_rich_club_hub { "RC" } else { "-" });
                ui.end_row();
            }
        });
    }

    fn ui_vault(&mut self, ui: &mut egui::Ui) {
        ui.heading("SQLite WAL Vault & Envelope I/O");
        ui.horizontal(|ui| {
            ui.label("Slot:");
            ui.text_edit_singleline(&mut self.slot_name_buf);
        });
        ui.horizontal_wrapped(|ui| {
            if ui.button("Save Snapshot to SQLite WAL").clicked() {
                if let Some(v) = &self.vault {
                    match v.save_snapshot(&self.slot_name_buf, &self.atlas) {
                        Ok(sha) => {
                            self.status_msg = format!(
                                "Saved slot '{}' (SHA-256 {}..)",
                                self.slot_name_buf,
                                &sha[..sha.len().min(12)]
                            );
                        }
                        Err(e) => self.status_msg = format!("Vault error: {e}"),
                    }
                } else {
                    self.status_msg = "No SQLite vault open; click 'Open SQLite Vault...' first.".to_string();
                }
            }
            if ui.button("Open SQLite Vault...").clicked() {
                if let Some(path) = rfd::FileDialog::new()
                    .add_filter("SQLite Vault", &["sqlite", "sqlite3", "db"])
                    .pick_file()
                {
                    match AtlasVault::open(&path) {
                        Ok(v) => {
                            self.vault = Some(v);
                            self.status_msg = format!("Opened SQLite vault {}", path.display());
                        }
                        Err(e) => self.status_msg = format!("Failed opening vault: {e}"),
                    }
                }
            }
            if ui.button("Open File (.json/.swc/.trk/.nii)...").clicked() {
                self.open_dataset_dialog();
            }
            if ui.button("Export Selected .SWC 3D Arbor...").clicked() {
                self.export_swc_dialog();
            }
            if ui.button("Export .connectome.json Envelope...").clicked() {
                self.export_envelope_dialog();
            }
            if ui.button("Export Topology Digest JSON...").clicked() {
                self.export_topology_digest_dialog(false);
            }
            if ui.button("Export Topology Digest CSV...").clicked() {
                self.export_topology_digest_dialog(true);
            }
        });
        ui.separator();
        let mut load_slot_name: Option<String> = None;
        if let Some(v) = &self.vault {
            match v.list_snapshots() {
                Ok(list) => {
                    for snap in list {
                        ui.horizontal(|ui| {
                            if ui.button("Load").clicked() {
                                load_slot_name = Some(snap.slot_name.clone());
                            }
                            ui.label(format!("{} ({})", snap.slot_name, snap.atlas_id));
                        });
                    }
                }
                Err(e) => {
                    ui.small(format!("Vault read error: {e}"));
                }
            }
        }
        if let Some(slot_name) = load_slot_name {
            let loaded_res = self.vault.as_ref().map(|v| v.load_snapshot(&slot_name));
            match loaded_res {
                Some(Ok(loaded)) => {
                    self.atlas = loaded;
                    self.conduction_velocity_m_s = self.atlas.conduction_velocity_m_s as f32;
                    self.selected_region = self.atlas.regions.first().map(|r| r.id.clone());
                    self.path_source = self.atlas.regions.first().map(|r| r.id.clone());
                    self.path_target = self.atlas.regions.last().map(|r| r.id.clone());
                    self.apply_species_framing();
                    self.recompute();
                    self.status_msg = format!("Loaded snapshot '{slot_name}'.");
                }
                Some(Err(e)) => {
                    self.status_msg = format!("Snapshot load error for '{slot_name}': {e}");
                }
                None => {}
            }
        }
    }

    fn draw_3d_viewport(&mut self, ui: &mut egui::Ui) {
        let size = ui.available_size();
        let (response, painter) = ui.allocate_painter(size, Sense::click_and_drag());
        let rect = response.rect;
        painter.rect_filled(rect, 0.0, Color32::from_rgb(9, 13, 21));

        // Optional chase-camera fly-through along shortest signaling pathway
        if self.fly_ctrl.active {
            if let Some(p) = &self.active_path {
                let mut coord_lookup = HashMap::new();
                for r in &self.atlas.regions {
                    coord_lookup.insert(r.id.as_str(), r.coord);
                }
                let poly: Vec<[f32; 3]> = p
                    .hops
                    .iter()
                    .filter_map(|h| coord_lookup.get(h.as_str()).copied())
                    .collect();
                if let Some((eye, look)) = self.fly_ctrl.sample_pathway_pose(&poly) {
                    let diff = eye - look;
                    let dist = diff.length().max(24.0);
                    self.camera.target = look;
                    self.camera.distance = dist.clamp(65.0, 420.0);
                    self.camera.pitch_rad = (diff.z / dist).clamp(-0.98, 0.98).asin();
                    self.camera.yaw_rad = diff.x.atan2(diff.y);
                }
            }
        }

        if response.dragged_by(egui::PointerButton::Secondary)
            || response.dragged_by(egui::PointerButton::Middle)
        {
            let d = response.drag_delta();
            self.fly_ctrl.active = false;
            self.camera.pan(glam::Vec2::new(d.x, d.y), rect.height());
        } else if response.dragged_by(egui::PointerButton::Primary) {
            let d = response.drag_delta();
            self.fly_ctrl.active = false;
            self.camera.orbit(glam::Vec2::new(d.x, d.y));
        }
        if response.hovered() {
            let scroll = ui.input(|i| i.raw_scroll_delta.y);
            if scroll.abs() > 0.01 {
                self.camera.zoom(scroll);
            }
        }

        let wh = [rect.width(), rect.height()];
        let vp = self.camera.view_proj_matrix(wh[0] / wh[1].max(1.0));
        let to_screen = |pt: [f32; 3]| -> Option<(Pos2, f32)> {
            OrbitCamera::project_point_with_vp(vp, pt, wh)
                .map(|([x, y], z)| (Pos2::new(rect.left() + x, rect.top() + y), z))
        };

        let eye = self.camera.eye_position();
        let light_key = Vec3::new(0.42, 0.58, 0.70).normalize();
        let light_rim = Vec3::new(-0.55, -0.35, 0.35).normalize();

        struct ShadedTri {
            depth: f32,
            p: [Pos2; 3],
            c: [Color32; 3],
        }

        let mut shell_back_tris: Vec<ShadedTri> = Vec::new();
        let mut organ_tris: Vec<ShadedTri> = Vec::new();
        let mut shell_front_tris: Vec<ShadedTri> = Vec::new();
        let mut cutaway_edges: Vec<(Pos2, Pos2)> = Vec::new();

        let shade_vertex = |wp: Vec3,
                            wn: Vec3,
                            cr: f32,
                            base_rgb: [u8; 3],
                            rim_rgb: [u8; 3],
                            layer: AnatomicalLayer,
                            op_scale: f32,
                            is_front: bool|
         -> Color32 {
            let n = if is_front { wn } else { -wn }.normalize_or_zero();
            let v = (eye - wp).normalize_or_zero();
            let n_dot_l = n.dot(light_key).max(0.0);
            let n_dot_fill = n.dot(light_rim).max(0.0);
            let half_vec = (light_key + v).normalize_or_zero();
            let spec = n.dot(half_vec).max(0.0).powf(24.0);
            let fresnel = (1.0 - n.dot(v).abs().clamp(0.0, 1.0)).powf(1.85);
            let ao = (0.72 + 0.28 * cr).clamp(0.42, 1.15);

            match layer {
                AnatomicalLayer::InternalOrgan => {
                    let diff = (0.28 + 0.62 * n_dot_l + 0.22 * n_dot_fill) * ao;
                    let r = (base_rgb[0] as f32 * diff
                        + rim_rgb[0] as f32 * (fresnel * 0.42 + spec * 0.36))
                        .clamp(0.0, 255.0) as u8;
                    let g = (base_rgb[1] as f32 * diff
                        + rim_rgb[1] as f32 * (fresnel * 0.42 + spec * 0.36))
                        .clamp(0.0, 255.0) as u8;
                    let b = (base_rgb[2] as f32 * diff
                        + rim_rgb[2] as f32 * (fresnel * 0.42 + spec * 0.36))
                        .clamp(0.0, 255.0) as u8;
                    let alpha = (232.0 * op_scale).clamp(48.0, 244.0) as u8;
                    Color32::from_rgba_unmultiplied(r, g, b, alpha)
                }
                AnatomicalLayer::OuterShell => {
                    let diff = (0.34 + 0.56 * n_dot_l + 0.26 * n_dot_fill) * ao;
                    let rim_boost = fresnel * if is_front { 0.85 } else { 0.42 };
                    let r = (base_rgb[0] as f32 * diff + rim_rgb[0] as f32 * (rim_boost + spec * 0.32))
                        .clamp(0.0, 255.0) as u8;
                    let g = (base_rgb[1] as f32 * diff + rim_rgb[1] as f32 * (rim_boost + spec * 0.32))
                        .clamp(0.0, 255.0) as u8;
                    let b = (base_rgb[2] as f32 * diff + rim_rgb[2] as f32 * (rim_boost + spec * 0.32))
                        .clamp(0.0, 255.0) as u8;
                    let raw_alpha = if is_front {
                        self.shell_opacity * op_scale * (0.11 + 1.05 * fresnel) * 235.0
                    } else {
                        self.shell_opacity * op_scale * (0.24 + 0.35 * fresnel) * 165.0
                    };
                    let alpha = raw_alpha.clamp(8.0, 210.0) as u8;
                    Color32::from_rgba_unmultiplied(r, g, b, alpha)
                }
            }
        };

        for mesh in &self.anatomy_scene.meshes {
            let is_shell = mesh.layer == AnatomicalLayer::OuterShell;
            if is_shell && self.shell_opacity <= 0.01 {
                continue;
            }
            if !is_shell && !self.show_internal_organs {
                continue;
            }
            let n_verts = mesh.vertices.len();
            let mut proj: Vec<Option<(Pos2, f32)>> = Vec::with_capacity(n_verts);
            let mut clipped: Vec<bool> = Vec::with_capacity(n_verts);
            for v in &mesh.vertices {
                let c = mesh.clippable
                    && is_point_clipped_by_cutaway(*v, self.cutaway_mode, self.cutaway_offset_mm);
                clipped.push(c);
                proj.push(to_screen(*v));
            }

            for tri in &mesh.triangles {
                let (i0, i1, i2) = (tri[0] as usize, tri[1] as usize, tri[2] as usize);
                let c0 = clipped[i0];
                let c1 = clipped[i1];
                let c2 = clipped[i2];
                if c0 || c1 || c2 {
                    // Collect glowing surgical incision contour along partially clipped triangles
                    let count_in = (!c0 as u8) + (!c1 as u8) + (!c2 as u8);
                    if count_in >= 1 && count_in <= 2 {
                        let pairs = [(i0, i1, c0, c1), (i1, i2, c1, c2), (i2, i0, c2, c0)];
                        let mut cut_pts: Vec<Pos2> = Vec::with_capacity(2);
                        for (a, b, ca, cb) in pairs {
                            if ca != cb {
                                let mid = Vec3::from_array(mesh.vertices[a])
                                    .lerp(Vec3::from_array(mesh.vertices[b]), 0.5);
                                if let Some((sp, _)) = to_screen(mid.to_array()) {
                                    cut_pts.push(sp);
                                }
                            }
                        }
                        if cut_pts.len() >= 2 {
                            cutaway_edges.push((cut_pts[0], cut_pts[1]));
                        }
                    }
                    continue;
                }

                let (Some((p0, d0)), Some((p1, d1)), Some((p2, d2))) =
                    (proj[i0], proj[i1], proj[i2])
                else {
                    continue;
                };

                let w0 = Vec3::from_array(mesh.vertices[i0]);
                let w1 = Vec3::from_array(mesh.vertices[i1]);
                let w2 = Vec3::from_array(mesh.vertices[i2]);
                let face_n = (w1 - w0).cross(w2 - w0);
                let face_center = (w0 + w1 + w2) / 3.0;
                let is_front = face_n.dot(eye - face_center) > 0.0;
                if !is_shell && !is_front {
                    continue;
                }

                let tri_item = ShadedTri {
                    depth: (d0 + d1 + d2) / 3.0,
                    p: [p0, p1, p2],
                    c: [
                        shade_vertex(
                            w0,
                            Vec3::from_array(mesh.normals[i0]),
                            mesh.crease[i0],
                            mesh.base_rgb,
                            mesh.rim_rgb,
                            mesh.layer,
                            mesh.opacity_scale,
                            is_front,
                        ),
                        shade_vertex(
                            w1,
                            Vec3::from_array(mesh.normals[i1]),
                            mesh.crease[i1],
                            mesh.base_rgb,
                            mesh.rim_rgb,
                            mesh.layer,
                            mesh.opacity_scale,
                            is_front,
                        ),
                        shade_vertex(
                            w2,
                            Vec3::from_array(mesh.normals[i2]),
                            mesh.crease[i2],
                            mesh.base_rgb,
                            mesh.rim_rgb,
                            mesh.layer,
                            mesh.opacity_scale,
                            is_front,
                        ),
                    ],
                };
                if is_shell {
                    if is_front {
                        shell_front_tris.push(tri_item);
                    } else {
                        shell_back_tris.push(tri_item);
                    }
                } else {
                    organ_tris.push(tri_item);
                }
            }
        }

        let flush_tris = |tris: &mut Vec<ShadedTri>| {
            tris.sort_by(|a, b| b.depth.total_cmp(&a.depth));
            let mut emesh = EguiMesh::default();
            emesh.vertices.reserve(tris.len() * 3);
            emesh.indices.reserve(tris.len() * 3);
            for t in tris.iter() {
                let base = emesh.vertices.len() as u32;
                for k in 0..3 {
                    emesh.vertices.push(EguiVertex {
                        pos: t.p[k],
                        uv: WHITE_UV,
                        color: t.c[k],
                    });
                    emesh.indices.push(base + k as u32);
                }
            }
            if !emesh.indices.is_empty() {
                painter.add(egui::Shape::mesh(emesh));
            }
        };

        // Pass 1: Back-facing translucent outer cortical / body shell
        flush_tris(&mut shell_back_tris);

        // Pass 2: Shaded 3D internal anatomical organs (cerebellum, brainstem, callosum, pharynx, central complex)
        flush_tris(&mut organ_tris);

        // Pass 3: Optional Axial Stereotaxic Slicing Grid scaled to active species
        if self.show_slice_planes {
            let ext = self.anatomy_scene.slice_extent_mm;
            let step = if ext > 80.0 { 25 } else { 20 };
            let ext_i = ((ext / step as f32).floor() as i32) * step;
            let grid_col = Color32::from_rgba_unmultiplied(56, 189, 248, 30);
            for g in (-ext_i..=ext_i).step_by(step as usize) {
                let gf = g as f32;
                if let (Some((a, _)), Some((b, _))) = (
                    to_screen([-ext, gf, self.slice_z_mm]),
                    to_screen([ext, gf, self.slice_z_mm]),
                ) {
                    painter.line_segment([a, b], Stroke::new(1.0_f32, grid_col));
                }
                if let (Some((a, _)), Some((b, _))) = (
                    to_screen([gf, -ext, self.slice_z_mm]),
                    to_screen([gf, ext, self.slice_z_mm]),
                ) {
                    painter.line_segment([a, b], Stroke::new(1.0_f32, grid_col));
                }
            }
        }

        let mut coord_map = HashMap::new();
        let mut color_map = HashMap::new();
        for r in &self.atlas.regions {
            coord_map.insert(r.id.as_str(), r.coord);
            color_map.insert(r.id.as_str(), parse_hex_color(&r.color_hex));
        }
        let mut betw_map = HashMap::new();
        for nm in &self.metrics.node_metrics {
            betw_map.insert(nm.region_id.as_str(), nm.betweenness_centrality);
        }

        // Field lookup maps for Harmonic psi_k, Criticality sigma_i, and TDA beta-1 loops
        let mut harmonic_map: HashMap<&str, f32> = HashMap::new();
        if let Some(mode_vec) = self
            .control_harmonics
            .harmonic_modes
            .get(self.harmonic_mode_idx.saturating_sub(1))
        {
            let max_abs = mode_vec
                .iter()
                .map(|v| v.abs())
                .fold(1e-6_f64, f64::max);
            for (r, &val) in self.atlas.regions.iter().zip(mode_vec.iter()) {
                harmonic_map.insert(r.id.as_str(), (val / max_abs) as f32);
            }
        }
        let mut crit_map: HashMap<&str, f32> = HashMap::new();
        for cp in &self.criticality.node_profiles {
            crit_map.insert(cp.region_id.as_str(), cp.avalanche_susceptibility as f32);
        }
        let mut tda_node_map: HashMap<&str, (usize, f32)> = HashMap::new();
        for tp in &self.topology_tda.node_profiles {
            tda_node_map.insert(
                tp.region_id.as_str(),
                (tp.beta1_loop_participation, tp.topological_cavity_score as f32),
            );
        }
        let mut tda_loop_edges: HashSet<(&str, &str)> = HashSet::new();
        for cav in &self.topology_tda.cavities {
            let cyc = &cav.cycle_nodes;
            if cyc.len() >= 2 {
                for w in cyc.windows(2) {
                    tda_loop_edges.insert((w[0].as_str(), w[1].as_str()));
                    tda_loop_edges.insert((w[1].as_str(), w[0].as_str()));
                }
                let first = cyc[0].as_str();
                let last = cyc[cyc.len() - 1].as_str();
                tda_loop_edges.insert((first, last));
                tda_loop_edges.insert((last, first));
            }
        }

        let mut path_pairs = Vec::new();
        if let Some(p) = &self.active_path {
            for w in p.hops.windows(2) {
                path_pairs.push((w[0].as_str(), w[1].as_str()));
            }
        }

        // Pass 4: Curved 3D White-Matter / Synaptic Streamlines with Calm Tonic Laminar Flow
        for edge in &self.atlas.edges {
            if !self.atlas.is_edge_intact(edge, self.min_weight) {
                continue;
            }
            let (Some(&c_a), Some(&c_b)) = (
                coord_map.get(edge.source.as_str()),
                coord_map.get(edge.target.as_str()),
            ) else {
                continue;
            };
            let on_shortest_path = path_pairs.iter().any(|&(u, v)| {
                (u == edge.source && v == edge.target)
                    || (!edge.directed && u == edge.target && v == edge.source)
            });
            let on_tda_loop =
                tda_loop_edges.contains(&(edge.source.as_str(), edge.target.as_str()));

            let base_rgb = match self.color_mode {
                ColorSchemeMode::FunctionalNetwork => color_map
                    .get(edge.source.as_str())
                    .copied()
                    .unwrap_or([56, 189, 248]),
                ColorSchemeMode::DtiDirectionalRgb => {
                    let dec = DecColor::from_endpoints(c_a, c_b);
                    [dec.r, dec.g, dec.b]
                }
                ColorSchemeMode::BetweennessHeatmap => {
                    let ba = betw_map.get(edge.source.as_str()).copied().unwrap_or(0.0) as f32;
                    let bb = betw_map.get(edge.target.as_str()).copied().unwrap_or(0.0) as f32;
                    betweenness_rgb(((ba + bb) * 0.5 * 4.5).clamp(0.0, 1.0))
                }
                ColorSchemeMode::HarmonicModePsi => {
                    let ha = harmonic_map.get(edge.source.as_str()).copied().unwrap_or(0.0);
                    let hb = harmonic_map.get(edge.target.as_str()).copied().unwrap_or(0.0);
                    diverging_harmonic_rgb((ha + hb) * 0.5)
                }
                ColorSchemeMode::CriticalitySusceptibility => {
                    let ca = crit_map.get(edge.source.as_str()).copied().unwrap_or(0.0);
                    let cb = crit_map.get(edge.target.as_str()).copied().unwrap_or(0.0);
                    thermal_plasma_rgb((ca + cb) * 0.5)
                }
                ColorSchemeMode::TdaBeta1Loops => {
                    if on_tda_loop {
                        [45, 226, 200]
                    } else {
                        [62, 82, 118]
                    }
                }
            };

            let strands = match self.tract_style {
                TractRenderStyle::BraidedFascicles => {
                    generate_braided_fascicle_bundle(c_a, c_b, edge.modality, edge.weight, 4, 14)
                }
                TractRenderStyle::SingleStreamline => {
                    vec![generate_arched_tract_curve(c_a, c_b, edge.modality, 14)]
                }
            };
            for (strand_i, curve) in strands.iter().enumerate() {
                let n_seg = curve.len().saturating_sub(1);
                for s_idx in 0..n_seg {
                    let (Some((p0, _)), Some((p1, _))) =
                        (to_screen(curve[s_idx]), to_screen(curve[s_idx + 1]))
                    else {
                        continue;
                    };
                    let u_frac = (s_idx as f32) / (n_seg.max(1) as f32);
                    let wave =
                        0.78 + 0.22 * ((u_frac - self.tonic_phase) * std::f32::consts::TAU).sin();
                    let seg_rgb = match self.color_mode {
                        ColorSchemeMode::HarmonicModePsi if !on_shortest_path => {
                            let ha = harmonic_map.get(edge.source.as_str()).copied().unwrap_or(0.0);
                            let hb = harmonic_map.get(edge.target.as_str()).copied().unwrap_or(0.0);
                            diverging_harmonic_rgb(ha + (hb - ha) * u_frac)
                        }
                        ColorSchemeMode::CriticalitySusceptibility if !on_shortest_path => {
                            let ca = crit_map.get(edge.source.as_str()).copied().unwrap_or(0.0);
                            let cb = crit_map.get(edge.target.as_str()).copied().unwrap_or(0.0);
                            thermal_plasma_rgb(ca + (cb - ca) * u_frac)
                        }
                        _ => base_rgb,
                    };
                    let (col, width) = if on_shortest_path {
                        (
                            Color32::from_rgba_unmultiplied(
                                253,
                                (224.0 * wave) as u8,
                                71,
                                242,
                            ),
                            3.6,
                        )
                    } else {
                        let base_alpha = if self.color_mode == ColorSchemeMode::TdaBeta1Loops {
                            if on_tda_loop {
                                (edge.weight as f32 * 135.0 + 95.0) * wave
                            } else {
                                (edge.weight as f32 * 45.0 + 18.0) * wave
                            }
                        } else {
                            (edge.weight as f32 * 145.0 + 42.0) * wave
                        };
                        let alpha = base_alpha.clamp(16.0, 232.0) as u8;
                        let boost = if self.color_mode == ColorSchemeMode::TdaBeta1Loops && on_tda_loop {
                            0.55
                        } else {
                            0.0
                        };
                        (
                            Color32::from_rgba_unmultiplied(
                                seg_rgb[0],
                                seg_rgb[1],
                                seg_rgb[2],
                                alpha,
                            ),
                            if strand_i == 0 {
                                1.25 + (edge.weight as f32) * 1.1 + boost
                            } else {
                                0.70 + (edge.weight as f32) * 0.55 + boost * 0.4
                            },
                        )
                    };
                    painter.line_segment([p0, p1], Stroke::new(width, col));
                }
            }
        }

        // Pass 5: 3D Crossing-Fiber Diffusion Tensor Ellipsoid Glyphs (3-Axis Shaded + Orthogonal Axes)
        if self.show_tensor_glyphs {
            let mut tensor_tris: Vec<ShadedTri> = Vec::new();
            for r in &self.atlas.regions {
                if self.atlas.ablated_regions.contains(&r.id) {
                    continue;
                }
                if let Some(ft) = self.fiber_tensors.iter().find(|t| t.region_id == r.id) {
                    if ft.fractional_anisotropy < 0.14 {
                        continue;
                    }
                    let u = Vec3::from_array(ft.principal_axis).normalize_or_zero();
                    let span = (r.radius_mm * (1.25 + 0.80 * ft.fractional_anisotropy as f32)).min(9.5);
                    let glyph = generate_tensor_ellipsoid_mesh(r.coord, ft, span, 4, 8);
                    let rgb = [
                        (u.x.abs() * 230.0 + 25.0) as u8,
                        (u.y.abs() * 230.0 + 25.0) as u8,
                        (u.z.abs() * 230.0 + 25.0) as u8,
                    ];

                    let g_proj: Vec<Option<(Pos2, f32)>> =
                        glyph.vertices.iter().map(|&v| to_screen(v)).collect();
                    for tri in &glyph.triangles {
                        let (i0, i1, i2) = (tri[0] as usize, tri[1] as usize, tri[2] as usize);
                        let (Some((p0, d0)), Some((p1, d1)), Some((p2, d2))) =
                            (g_proj[i0], g_proj[i1], g_proj[i2])
                        else {
                            continue;
                        };
                        let w0 = Vec3::from_array(glyph.vertices[i0]);
                        let w1 = Vec3::from_array(glyph.vertices[i1]);
                        let w2 = Vec3::from_array(glyph.vertices[i2]);
                        let face_n = (w1 - w0).cross(w2 - w0);
                        let face_c = (w0 + w1 + w2) / 3.0;
                        if face_n.dot(eye - face_c) <= 0.0 {
                            continue;
                        }
                        let shade_g = |wp: Vec3, wn: Vec3| -> Color32 {
                            let n = wn.normalize_or_zero();
                            let v = (eye - wp).normalize_or_zero();
                            let diff = 0.38 + 0.56 * n.dot(light_key).max(0.0) + 0.22 * n.dot(light_rim).max(0.0);
                            let fres = (1.0 - n.dot(v).abs().clamp(0.0, 1.0)).powf(1.8);
                            Color32::from_rgba_unmultiplied(
                                (rgb[0] as f32 * diff + 180.0 * fres * 0.25).clamp(0.0, 255.0) as u8,
                                (rgb[1] as f32 * diff + 210.0 * fres * 0.25).clamp(0.0, 255.0) as u8,
                                (rgb[2] as f32 * diff + 240.0 * fres * 0.25).clamp(0.0, 255.0) as u8,
                                92,
                            )
                        };
                        tensor_tris.push(ShadedTri {
                            depth: (d0 + d1 + d2) / 3.0,
                            p: [p0, p1, p2],
                            c: [
                                shade_g(w0, Vec3::from_array(glyph.normals[i0])),
                                shade_g(w1, Vec3::from_array(glyph.normals[i1])),
                                shade_g(w2, Vec3::from_array(glyph.normals[i2])),
                            ],
                        });
                    }

                    for (ax_idx, &(ea, eb)) in glyph.axis_segments.iter().enumerate() {
                        if let (Some((pa, _)), Some((pb, _))) = (to_screen(ea), to_screen(eb)) {
                            let (ax_col, ax_w) = if ax_idx == 0 {
                                (Color32::from_rgba_unmultiplied(rgb[0], rgb[1], rgb[2], 195), 1.8_f32)
                            } else {
                                (Color32::from_rgba_unmultiplied(186, 230, 253, 120), 1.0_f32)
                            };
                            painter.line_segment([pa, pb], Stroke::new(ax_w, ax_col));
                        }
                    }
                }
            }
            flush_tris(&mut tensor_tris);
        }

        // Pass 6: Depth-Sorted 3D Anatomical Region Parcels / Somas
        let mut projected_nodes = Vec::new();
        for r in &self.atlas.regions {
            if let Some((sp, depth)) = to_screen(r.coord) {
                projected_nodes.push((r, sp, depth));
            }
        }
        projected_nodes.sort_by(|a, b| b.2.total_cmp(&a.2));

        let click_pos = if response.clicked() {
            response.interact_pointer_pos()
        } else {
            None
        };
        let mut clicked_region: Option<String> = None;

        for (r, sp, _depth) in projected_nodes {
            let ablated = self.atlas.ablated_regions.contains(&r.id);
            let is_selected = self.selected_region.as_deref() == Some(r.id.as_str());
            let radius_px = (r.radius_mm * 235.0 / self.camera.distance).clamp(4.2, 13.5);

            if let Some(cp) = click_pos {
                if cp.distance(sp) <= radius_px + 5.0 {
                    clicked_region = Some(r.id.clone());
                }
            }

            let rgb = if ablated {
                [65, 70, 82]
            } else {
                match self.color_mode {
                    ColorSchemeMode::BetweennessHeatmap => {
                        let b = betw_map.get(r.id.as_str()).copied().unwrap_or(0.0) as f32;
                        betweenness_rgb((b * 4.5).clamp(0.0, 1.0))
                    }
                    ColorSchemeMode::HarmonicModePsi => {
                        let h = harmonic_map.get(r.id.as_str()).copied().unwrap_or(0.0);
                        diverging_harmonic_rgb(h)
                    }
                    ColorSchemeMode::CriticalitySusceptibility => {
                        let s = crit_map.get(r.id.as_str()).copied().unwrap_or(0.0);
                        thermal_plasma_rgb(s)
                    }
                    ColorSchemeMode::TdaBeta1Loops => {
                        let (loop_cnt, score) =
                            tda_node_map.get(r.id.as_str()).copied().unwrap_or((0, 0.0));
                        if loop_cnt > 0 {
                            let t = score.clamp(0.25, 1.0);
                            [
                                (38.0 + 205.0 * t) as u8,
                                (212.0 + 28.0 * (1.0 - t)) as u8,
                                (205.0 - 135.0 * t) as u8,
                            ]
                        } else {
                            [68, 92, 130]
                        }
                    }
                    _ => parse_hex_color(&r.color_hex),
                }
            };

            painter.circle_filled(sp, radius_px + 1.2, Color32::from_rgba_unmultiplied(6, 10, 18, 210));
            painter.circle_filled(sp, radius_px, Color32::from_rgb(rgb[0], rgb[1], rgb[2]));
            // Specular sphere highlight
            painter.circle_filled(
                Pos2::new(sp.x - radius_px * 0.28, sp.y - radius_px * 0.28),
                (radius_px * 0.30).max(1.2),
                Color32::from_rgba_unmultiplied(255, 255, 255, 145),
            );
            if ablated {
                painter.circle_stroke(
                    sp,
                    radius_px + 2.5,
                    Stroke::new(1.8_f32, Color32::from_rgb(239, 68, 68)),
                );
            } else if is_selected {
                painter.circle_stroke(sp, radius_px + 3.0, Stroke::new(2.2_f32, Color32::WHITE));
            }

            if is_selected || radius_px >= 6.4 {
                painter.text(
                    Pos2::new(sp.x + 1.0, sp.y - radius_px - 2.0),
                    egui::Align2::CENTER_BOTTOM,
                    &r.id,
                    FontId::proportional(10.5),
                    Color32::from_rgba_unmultiplied(8, 12, 20, 220),
                );
                painter.text(
                    Pos2::new(sp.x, sp.y - radius_px - 3.0),
                    egui::Align2::CENTER_BOTTOM,
                    &r.id,
                    FontId::proportional(10.5),
                    Color32::from_rgba_unmultiplied(241, 245, 249, 235),
                );
            }
        }

        // Pass 7: Front-Facing Fresnel Rim-Lit Cortical / Organism Outer Glass Pass + Cutaway Edge Glow
        flush_tris(&mut shell_front_tris);
        if self.shell_opacity > 0.01 && !cutaway_edges.is_empty() {
            let cut_col = Color32::from_rgba_unmultiplied(103, 232, 249, 175);
            for (pa, pb) in cutaway_edges {
                painter.line_segment([pa, pb], Stroke::new(1.6_f32, cut_col));
            }
        }

        // Pass 8: Species Anatomical Landmark Callouts (inward-aligned so edges never clip)
        if self.show_anatomical_labels {
            let mid_x = rect.center().x;
            for lm in &self.anatomy_scene.landmarks {
                if let Some((sp, _)) = to_screen(lm.coord) {
                    if !rect.contains(sp) {
                        continue;
                    }
                    let on_right = sp.x > mid_x + 35.0;
                    let tick_dx = if on_right { -10.0 } else { 10.0 };
                    let tick_dy = if sp.y < rect.top() + 36.0 { 10.0 } else { -10.0 };
                    let tick_end = Pos2::new(
                        (sp.x + tick_dx).clamp(rect.left() + 6.0, rect.right() - 6.0),
                        (sp.y + tick_dy).clamp(rect.top() + 10.0, rect.bottom() - 10.0),
                    );
                    let align = if on_right {
                        egui::Align2::RIGHT_CENTER
                    } else {
                        egui::Align2::LEFT_CENTER
                    };
                    let text_dx = if on_right { -3.0 } else { 3.0 };
                    let tick_col =
                        Color32::from_rgba_unmultiplied(lm.color[0], lm.color[1], lm.color[2], 160);
                    painter.circle_filled(sp, 2.2, tick_col);
                    painter.line_segment([sp, tick_end], Stroke::new(1.0_f32, tick_col));
                    painter.text(
                        Pos2::new(tick_end.x + text_dx + 1.0, tick_end.y + 1.0),
                        align,
                        lm.label,
                        FontId::proportional(10.5),
                        Color32::from_rgba_unmultiplied(6, 10, 18, 225),
                    );
                    painter.text(
                        Pos2::new(tick_end.x + text_dx, tick_end.y),
                        align,
                        lm.label,
                        FontId::proportional(10.5),
                        Color32::from_rgba_unmultiplied(lm.color[0], lm.color[1], lm.color[2], 225),
                    );
                }
            }
        }

        if let Some(cid) = clicked_region {
            self.select_region(cid.clone(), true);
            match self.probe_mode {
                ProbeMode::Inspect => {
                    self.status_msg = format!("Inspecting region {cid}.");
                }
                ProbeMode::ShortestPath => {
                    if self.path_source.is_none() || self.path_target.is_some() {
                        self.path_source = Some(cid.clone());
                        self.path_target = None;
                        self.active_path = None;
                        self.status_msg = format!("Signaling path source set to {cid} - click target region.");
                    } else {
                        self.path_target = Some(cid.clone());
                        self.recompute();
                        self.status_msg = format!(
                            "Computed shortest signaling path {} -> {}.",
                            self.path_source.as_deref().unwrap_or(""),
                            cid
                        );
                    }
                }
                ProbeMode::LesionAblate => {
                    self.atlas.toggle_region_ablation(&cid);
                    self.recompute();
                    self.status_msg = format!("Toggled structural lesion on {cid}.");
                }
            }
        }
    }
}

impl ConnectomeAtlasApp {
    fn draw_triplanar_2d_strip(&mut self, ui: &mut egui::Ui) {
        let ext = self.anatomy_scene.slice_extent_mm.max(40.0);
        ui.columns(3, |cols| {
            let planes = [
                ("AXIAL (XY)", 0_u8, self.slice_z_mm, ["L", "R", "A", "P"]),
                ("CORONAL (XZ)", 1_u8, self.slice_y_mm, ["L", "R", "S", "I"]),
                ("SAGITTAL (YZ)", 2_u8, self.slice_x_mm, ["P", "A", "S", "I"]),
            ];
            for (idx, (title, mode, slice_v, dirs)) in planes.into_iter().enumerate() {
                let col_ui = &mut cols[idx];
                col_ui.label(egui::RichText::new(format!("{title} @ {slice_v:.0} mm")).small().strong());
                let avail = col_ui.available_size();
                let (resp, painter) = col_ui.allocate_painter(avail, Sense::click_and_drag());
                let r = resp.rect;
                painter.rect_filled(r, 4.0, Color32::from_rgb(6, 9, 17));
                let cx = r.center().x;
                let cy = r.center().y;
                let scale = (r.height() / (ext * 2.15))
                    .min(r.width() / (ext * 2.35))
                    .max(0.45);

                // Species-accurate anatomical cross-section contours
                let contours = generate_2d_slice_contours(&self.atlas.id, mode, slice_v);
                for contour in contours {
                    let screen_pts: Vec<Pos2> = contour
                        .points
                        .iter()
                        .map(|&[u, v]| Pos2::new(cx + u * scale, cy - v * scale))
                        .collect();
                    if contour.closed && screen_pts.len() >= 3 {
                        if contour.fill_rgba[3] > 0 {
                            let fill_col = Color32::from_rgba_unmultiplied(
                                contour.fill_rgba[0],
                                contour.fill_rgba[1],
                                contour.fill_rgba[2],
                                contour.fill_rgba[3],
                            );
                            let mut c_sum = Vec2::ZERO;
                            for &sp in &screen_pts {
                                c_sum += sp.to_vec2();
                            }
                            let centroid = (c_sum / (screen_pts.len() as f32)).to_pos2();
                            let mut em = EguiMesh::default();
                            em.vertices.push(EguiVertex { pos: centroid, uv: WHITE_UV, color: fill_col });
                            for &sp in &screen_pts {
                                em.vertices.push(EguiVertex { pos: sp, uv: WHITE_UV, color: fill_col });
                            }
                            let m = screen_pts.len() as u32;
                            for i in 0..m {
                                em.indices.push(0);
                                em.indices.push(1 + i);
                                em.indices.push(1 + (i + 1) % m);
                            }
                            painter.add(egui::Shape::mesh(em));
                        }
                        if contour.stroke_rgba[3] > 0 {
                            let stroke_col = Color32::from_rgba_unmultiplied(
                                contour.stroke_rgba[0],
                                contour.stroke_rgba[1],
                                contour.stroke_rgba[2],
                                contour.stroke_rgba[3],
                            );
                            painter.add(egui::Shape::closed_line(
                                screen_pts,
                                Stroke::new(contour.stroke_width, stroke_col),
                            ));
                        }
                    } else if screen_pts.len() >= 2 && contour.stroke_rgba[3] > 0 {
                        let stroke_col = Color32::from_rgba_unmultiplied(
                            contour.stroke_rgba[0],
                            contour.stroke_rgba[1],
                            contour.stroke_rgba[2],
                            contour.stroke_rgba[3],
                        );
                        painter.add(egui::Shape::line(
                            screen_pts,
                            Stroke::new(contour.stroke_width, stroke_col),
                        ));
                    }
                }

                let slab_half = if ext > 80.0 { 28.0 } else { 24.0 };
                let hover_pos = resp.hover_pos();
                let click_pos = if resp.clicked() { resp.interact_pointer_pos() } else { None };
                let mut hovered_best: Option<(f32, String)> = None;
                let mut clicked_best: Option<(f32, String)> = None;

                for reg in &self.atlas.regions {
                    let (u, v, ortho) = match mode {
                        0 => (reg.coord[0], reg.coord[1], reg.coord[2]),
                        1 => (reg.coord[0], reg.coord[2], reg.coord[1]),
                        _ => (reg.coord[1], reg.coord[2], reg.coord[0]),
                    };
                    if (ortho - slice_v).abs() <= slab_half {
                        let rgb = parse_hex_color(&reg.color_hex);
                        let pos = Pos2::new(cx + u * scale, cy - v * scale);
                        painter.circle_filled(pos, 3.0, Color32::from_rgb(rgb[0], rgb[1], rgb[2]));
                        if self.atlas.ablated_regions.contains(&reg.id) {
                            painter.circle_stroke(
                                pos,
                                4.8_f32,
                                Stroke::new(1.5_f32, Color32::from_rgb(239, 68, 68)),
                            );
                        }
                        if self.selected_region.as_deref() == Some(reg.id.as_str()) {
                            painter.circle_stroke(
                                pos,
                                5.2_f32,
                                Stroke::new(1.6_f32, Color32::WHITE),
                            );
                        }
                        if let Some(hp) = hover_pos {
                            let d = hp.distance(pos);
                            if d <= 7.5_f32 && hovered_best.as_ref().map_or(true, |(bd, _)| d < *bd) {
                                hovered_best = Some((
                                    d,
                                    format!(
                                        "{} — {}\n{} | ({:.1}, {:.1}, {:.1}) mm",
                                        reg.id, reg.name, reg.network, reg.coord[0], reg.coord[1], reg.coord[2]
                                    ),
                                ));
                            }
                        }
                        if let Some(cp) = click_pos {
                            let d = cp.distance(pos);
                            if d <= 7.5_f32 && clicked_best.as_ref().map_or(true, |(bd, _)| d < *bd) {
                                clicked_best = Some((d, reg.id.clone()));
                            }
                        }
                    }
                }

                if let Some((_, rid)) = clicked_best {
                    self.select_region(rid.clone(), true);
                    self.status_msg = format!("Selected region {rid} from 2D slice.");
                } else if resp.clicked() || resp.dragged() {
                    if let Some(ptr) = resp.interact_pointer_pos() {
                        let u_mm = ((ptr.x - cx) / scale).clamp(-ext, ext);
                        let v_mm = ((cy - ptr.y) / scale).clamp(-ext, ext);
                        match mode {
                            0 => { self.slice_x_mm = u_mm; self.slice_y_mm = v_mm; }
                            1 => { self.slice_x_mm = u_mm; self.slice_z_mm = v_mm; }
                            _ => { self.slice_y_mm = u_mm; self.slice_z_mm = v_mm; }
                        }
                    }
                }

                if let Some((_, tip)) = hovered_best {
                    resp.on_hover_text(tip);
                }

                let (ch_u, ch_v) = match mode {
                    0 => (self.slice_x_mm, self.slice_y_mm),
                    1 => (self.slice_x_mm, self.slice_z_mm),
                    _ => (self.slice_y_mm, self.slice_z_mm),
                };
                let ch_x = (cx + ch_u * scale).clamp(r.left(), r.right());
                let ch_y = (cy - ch_v * scale).clamp(r.top(), r.bottom());
                let line_col = Color32::from_rgba_unmultiplied(56, 189, 248, 120);
                painter.line_segment([Pos2::new(ch_x, r.top()), Pos2::new(ch_x, r.bottom())], Stroke::new(1.0_f32, line_col));
                painter.line_segment([Pos2::new(r.left(), ch_y), Pos2::new(r.right(), ch_y)], Stroke::new(1.0_f32, line_col));

                // Radiological axis orientation letters
                let dir_col = Color32::from_rgba_unmultiplied(148, 163, 184, 175);
                painter.text(Pos2::new(r.left() + 6.0, cy), egui::Align2::LEFT_CENTER, dirs[0], FontId::proportional(9.5), dir_col);
                painter.text(Pos2::new(r.right() - 6.0, cy), egui::Align2::RIGHT_CENTER, dirs[1], FontId::proportional(9.5), dir_col);
                painter.text(Pos2::new(cx, r.top() + 4.0), egui::Align2::CENTER_TOP, dirs[2], FontId::proportional(9.5), dir_col);
                painter.text(Pos2::new(cx, r.bottom() - 4.0), egui::Align2::CENTER_BOTTOM, dirs[3], FontId::proportional(9.5), dir_col);
            }
        });
    }
}

fn parse_hex_color(hex: &str) -> [u8; 3] {
    let s = hex.trim_start_matches('#');
    if s.len() == 6 {
        let r = u8::from_str_radix(&s[0..2], 16).unwrap_or(56);
        let g = u8::from_str_radix(&s[2..4], 16).unwrap_or(189);
        let b = u8::from_str_radix(&s[4..6], 16).unwrap_or(248);
        [r, g, b]
    } else {
        [56, 189, 248]
    }
}

fn lerp_rgb(a: [u8; 3], b: [u8; 3], t: f32) -> [u8; 3] {
    let u = t.clamp(0.0, 1.0);
    [
        (a[0] as f32 + (b[0] as f32 - a[0] as f32) * u) as u8,
        (a[1] as f32 + (b[1] as f32 - a[1] as f32) * u) as u8,
        (a[2] as f32 + (b[2] as f32 - a[2] as f32) * u) as u8,
    ]
}

fn betweenness_rgb(t: f32) -> [u8; 3] {
    lerp_rgb([50, 190, 245], [255, 95, 55], t)
}

fn diverging_harmonic_rgb(val: f32) -> [u8; 3] {
    let v = val.clamp(-1.0, 1.0);
    if v < 0.0 {
        lerp_rgb([100, 116, 139], [34, 211, 238], -v)
    } else {
        lerp_rgb([100, 116, 139], [251, 113, 62], v)
    }
}

fn thermal_plasma_rgb(t: f32) -> [u8; 3] {
    let u = t.clamp(0.0, 1.0);
    if u < 0.5 {
        lerp_rgb([37, 99, 235], [45, 212, 191], u * 2.0)
    } else {
        lerp_rgb([45, 212, 191], [250, 204, 21], (u - 0.5) * 2.0)
    }
}

fn matrix_weight_rgb(w: f32) -> [u8; 3] {
    let u = w.clamp(0.0, 1.0);
    if u < 0.5 {
        lerp_rgb([15, 34, 72], [34, 211, 238], u * 2.0)
    } else {
        lerp_rgb([34, 211, 238], [250, 204, 21], (u - 0.5) * 2.0)
    }
}

struct NoopWaker;

impl std::task::Wake for NoopWaker {
    fn wake(self: std::sync::Arc<Self>) {}
}

fn block_on_future<F: std::future::Future>(fut: F) -> F::Output {
    let mut fut = std::pin::pin!(fut);
    let waker = std::task::Waker::from(std::sync::Arc::new(NoopWaker));
    let mut cx = std::task::Context::from_waker(&waker);
    loop {
        match fut.as_mut().poll(&mut cx) {
            std::task::Poll::Ready(val) => return val,
            std::task::Poll::Pending => std::thread::yield_now(),
        }
    }
}

fn write_rgba_png(
    path: &std::path::Path,
    width: u32,
    height: u32,
    rgba: &[u8],
) -> Result<(), Box<dyn std::error::Error>> {
    let expected = (width as usize) * (height as usize) * 4;
    if rgba.len() != expected {
        return Err(format!(
            "Unexpected RGBA buffer length {} (expected {})",
            rgba.len(),
            expected
        )
        .into());
    }
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)?;
        }
    }
    let file = std::fs::File::create(path)?;
    let writer = std::io::BufWriter::new(file);
    let mut encoder = png::Encoder::new(writer, width, height);
    encoder.set_color(png::ColorType::Rgba);
    encoder.set_depth(png::BitDepth::Eight);
    encoder.set_compression(png::Compression::High);
    encoder.set_filter(png::Filter::Adaptive);
    let mut png_writer = encoder.write_header()?;
    png_writer.write_image_data(rgba)?;
    png_writer.finish()?;
    Ok(())
}

pub fn render_headless_screenshot(
    out_path: &std::path::Path,
    view: u8,
) -> Result<(), Box<dyn std::error::Error>> {
    use eframe::egui_wgpu::{wgpu, Renderer, ScreenDescriptor};

    let width: u32 = 1600;
    let height: u32 = 1000;
    let instance = wgpu::Instance::new(&wgpu::InstanceDescriptor::default());
    let adapter = block_on_future(instance.request_adapter(&wgpu::RequestAdapterOptions {
        power_preference: wgpu::PowerPreference::HighPerformance,
        compatible_surface: None,
        force_fallback_adapter: false,
    }))
    .ok_or("No GPU adapter available for offscreen rendering")?;
    let (device, queue) = block_on_future(adapter.request_device(
        &wgpu::DeviceDescriptor {
            label: Some("connectome-atlas-offscreen-device"),
            ..Default::default()
        },
        None,
    ))?;

    let format = wgpu::TextureFormat::Rgba8Unorm;
    let texture = device.create_texture(&wgpu::TextureDescriptor {
        label: Some("connectome-atlas-offscreen-target"),
        size: wgpu::Extent3d {
            width,
            height,
            depth_or_array_layers: 1,
        },
        mip_level_count: 1,
        sample_count: 1,
        dimension: wgpu::TextureDimension::D2,
        format,
        usage: wgpu::TextureUsages::RENDER_ATTACHMENT | wgpu::TextureUsages::COPY_SRC,
        view_formats: &[],
    });
    let texture_view = texture.create_view(&wgpu::TextureViewDescriptor::default());
    let bytes_per_row = width * 4;
    let staging = device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("connectome-atlas-offscreen-staging"),
        size: (bytes_per_row * height) as u64,
        usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
        mapped_at_creation: false,
    });

    let mut renderer = Renderer::new(&device, format, None, 1, false);
    let egui_ctx = egui::Context::default();
    egui_ctx.set_pixels_per_point(1.0);
    let mut visuals = egui::Visuals::dark();
    visuals.panel_fill = Color32::from_rgb(13, 17, 26);
    visuals.window_fill = Color32::from_rgb(16, 21, 32);
    egui_ctx.set_visuals(visuals);

    let mut app = match view {
        2 => {
            let mut a = ConnectomeAtlasApp::new_headless(2);
            a.right_tab = RightTab::AdjacencyMatrix;
            a.color_mode = ColorSchemeMode::DtiDirectionalRgb;
            a.cutaway_mode = SurgicalCutawayMode::Intact;
            a
        }
        3 => {
            let mut a = ConnectomeAtlasApp::new_headless(3);
            a.right_tab = RightTab::RichClubHubs;
            a.color_mode = ColorSchemeMode::BetweennessHeatmap;
            a.cutaway_mode = SurgicalCutawayMode::Intact;
            a
        }
        4 => {
            let mut a = ConnectomeAtlasApp::new_headless(1);
            a.right_tab = RightTab::Inspector;
            a.color_mode = ColorSchemeMode::FunctionalNetwork;
            a.cutaway_mode = SurgicalCutawayMode::AnteriorSuperiorOctant;
            a
        }
        _ => {
            let mut a = ConnectomeAtlasApp::new_headless(0);
            a.right_tab = RightTab::Inspector;
            a.cutaway_mode = SurgicalCutawayMode::AnteriorSuperiorOctant;
            a
        }
    };

    let warmup_delays_ms = [40u64, 20, 20, 20];
    for (frame_idx, delay_ms) in warmup_delays_ms.into_iter().enumerate() {
        std::thread::sleep(std::time::Duration::from_millis(delay_ms));
        let raw_input = egui::RawInput {
            screen_rect: Some(egui::Rect::from_min_size(
                egui::Pos2::ZERO,
                egui::vec2(width as f32, height as f32),
            )),
            viewports: std::iter::once((
                egui::ViewportId::ROOT,
                egui::ViewportInfo {
                    native_pixels_per_point: Some(1.0),
                    ..Default::default()
                },
            ))
            .collect(),
            time: Some(frame_idx as f64 / 60.0),
            predicted_dt: 1.0 / 60.0,
            ..Default::default()
        };
        let full_output = egui_ctx.run(raw_input, |ctx| {
            app.ui_step(ctx, false);
        });
        for (id, image_delta) in &full_output.textures_delta.set {
            renderer.update_texture(&device, &queue, *id, image_delta);
        }
        let paint_jobs = egui_ctx.tessellate(full_output.shapes, full_output.pixels_per_point);
        let screen_descriptor = ScreenDescriptor {
            size_in_pixels: [width, height],
            pixels_per_point: 1.0,
        };
        let mut encoder = device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("connectome-atlas-offscreen-encoder"),
        });
        let user_cmd_bufs = renderer.update_buffers(
            &device,
            &queue,
            &mut encoder,
            &paint_jobs,
            &screen_descriptor,
        );
        {
            let render_pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("connectome-atlas-offscreen-pass"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: &texture_view,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Clear(wgpu::Color::BLACK),
                        store: wgpu::StoreOp::Store,
                    },
                })],
                depth_stencil_attachment: None,
                occlusion_query_set: None,
                timestamp_writes: None,
            });
            let mut render_pass = render_pass.forget_lifetime();
            renderer.render(&mut render_pass, &paint_jobs, &screen_descriptor);
        }
        if frame_idx == 3 {
            encoder.copy_texture_to_buffer(
                wgpu::TexelCopyTextureInfo {
                    texture: &texture,
                    mip_level: 0,
                    origin: wgpu::Origin3d::ZERO,
                    aspect: wgpu::TextureAspect::All,
                },
                wgpu::TexelCopyBufferInfo {
                    buffer: &staging,
                    layout: wgpu::TexelCopyBufferLayout {
                        offset: 0,
                        bytes_per_row: Some(bytes_per_row),
                        rows_per_image: Some(height),
                    },
                },
                wgpu::Extent3d {
                    width,
                    height,
                    depth_or_array_layers: 1,
                },
            );
        }
        queue.submit(user_cmd_bufs.into_iter().chain(std::iter::once(encoder.finish())));
        device.poll(wgpu::Maintain::Wait);
        for id in &full_output.textures_delta.free {
            renderer.free_texture(id);
        }
    }

    let slice = staging.slice(..);
    let (tx, rx) = std::sync::mpsc::channel();
    slice.map_async(wgpu::MapMode::Read, move |res| {
        let _ = tx.send(res);
    });
    device.poll(wgpu::Maintain::Wait);
    rx.recv()??;
    let mapped = slice.get_mapped_range();
    write_rgba_png(out_path, width, height, &mapped)?;
    drop(mapped);
    staging.unmap();
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn select_region_syncs_slices_and_refreshes_domain_decode() {
        let mut app = ConnectomeAtlasApp::new_headless(0);
        let broca = app
            .atlas
            .regions
            .iter()
            .find(|r| r.id == "LANG.Broca_L")
            .expect("LANG.Broca_L exists in human preset")
            .clone();

        app.select_region("LANG.Broca_L".to_string(), true);
        assert_eq!(app.selected_region.as_deref(), Some("LANG.Broca_L"));
        assert!((app.slice_x_mm - broca.coord[0]).abs() < 1e-4);
        assert!((app.slice_y_mm - broca.coord[1]).abs() < 1e-4);
        assert!((app.slice_z_mm - broca.coord[2]).abs() < 1e-4);
        assert!(!app.domain_decode.is_empty());
        let lang_entry = app
            .domain_decode
            .iter()
            .find(|d| d.domain_label == "Language & Semantic")
            .expect("Language & Semantic domain entry present");
        assert!(
            (lang_entry.score - 0.88).abs() < 1e-6 && !lang_entry.is_lesion_disruption,
            "domain_decode should boost Language & Semantic for LANG.Broca_L, got {:?}",
            lang_entry
        );
    }

    #[test]
    fn switching_receptor_channel_updates_control_harmonics() {
        let mut app = ConnectomeAtlasApp::new_headless(0);
        app.agonist_gain = 1.0;
        for ch in ["ht", "da", "gaba"] {
            app.receptor_channel = ch.to_string();
            app.recompute();
            assert_eq!(app.control_harmonics.receptor_channel, ch);
        }
    }

    #[test]
    fn headless_f1_and_escape_shortcuts_toggle_and_cascade() {
        let ctx = egui::Context::default();
        let mut app = ConnectomeAtlasApp::new_headless(0);
        assert!(!app.show_about_window);

        let run_key = |app: &mut ConnectomeAtlasApp, key: Option<egui::Key>| {
            let mut raw = egui::RawInput {
                screen_rect: Some(egui::Rect::from_min_size(
                    egui::Pos2::ZERO,
                    egui::vec2(1400.0, 900.0),
                )),
                ..Default::default()
            };
            if let Some(k) = key {
                raw.events.push(egui::Event::Key {
                    key: k,
                    physical_key: None,
                    pressed: true,
                    repeat: false,
                    modifiers: egui::Modifiers::NONE,
                });
            }
            let _ = ctx.run(raw, |ctx| {
                app.ui_step(ctx, false);
            });
        };

        run_key(&mut app, None);
        run_key(&mut app, Some(egui::Key::F1));
        assert!(app.show_about_window, "F1 should open About window");
        run_key(&mut app, Some(egui::Key::F1));
        assert!(!app.show_about_window, "F1 should toggle About window closed");

        // Escape cascade: cmd_palette_open -> show_about_window -> fly_ctrl.active -> zen_mode
        app.cmd_palette_open = true;
        app.show_about_window = true;
        app.fly_ctrl.active = true;
        app.zen_mode = true;
        app.show_left_panel = false;
        app.show_right_panel = false;
        app.show_bottom_slices = false;

        run_key(&mut app, Some(egui::Key::Escape));
        assert!(!app.cmd_palette_open);
        assert!(app.show_about_window);

        run_key(&mut app, Some(egui::Key::Escape));
        assert!(!app.show_about_window);
        assert!(app.fly_ctrl.active);

        run_key(&mut app, Some(egui::Key::Escape));
        assert!(!app.fly_ctrl.active);
        assert!(app.zen_mode);

        run_key(&mut app, Some(egui::Key::Escape));
        assert!(!app.zen_mode);
        assert!(app.show_left_panel && app.show_right_panel && app.show_bottom_slices);
    }

    #[test]
    fn hotkeys_palette_enter_and_cached_lesion_delta_behave_accurately() {
        let ctx = egui::Context::default();
        let mut app = ConnectomeAtlasApp::new_headless(0);

        let send_event = |app: &mut ConnectomeAtlasApp, key: egui::Key, mods: egui::Modifiers, text: Option<&str>| {
            let mut raw = egui::RawInput {
                screen_rect: Some(egui::Rect::from_min_size(
                    egui::Pos2::ZERO,
                    egui::vec2(1400.0, 900.0),
                )),
                modifiers: mods,
                ..Default::default()
            };
            if let Some(t) = text {
                raw.events.push(egui::Event::Text(t.to_string()));
            }
            raw.events.push(egui::Event::Key {
                key,
                physical_key: None,
                pressed: true,
                repeat: false,
                modifiers: mods,
            });
            let _ = ctx.run(raw, |ctx| {
                app.ui_step(ctx, false);
            });
        };

        // Warmup frame
        let _ = ctx.run(
            egui::RawInput {
                screen_rect: Some(egui::Rect::from_min_size(
                    egui::Pos2::ZERO,
                    egui::vec2(1400.0, 900.0),
                )),
                ..Default::default()
            },
            |ctx| app.ui_step(ctx, false),
        );

        // Hotkeys 1-4 set presets and cancel fly-through
        app.fly_ctrl.active = true;
        send_event(&mut app, egui::Key::Num3, egui::Modifiers::NONE, None);
        assert!(!app.fly_ctrl.active);
        let (_, ax_pitch) = OrbitCamera::preset_angles(ViewPreset::SuperiorAxial);
        assert!((app.camera.pitch_rad - ax_pitch).abs() < 1e-3);

        // R resets camera to species default framing
        send_event(&mut app, egui::Key::R, egui::Modifiers::NONE, None);
        assert!((app.camera.yaw_rad - app.anatomy_scene.default_yaw).abs() < 1e-3);

        // I / P / L probe modes and C cutaway cycle
        send_event(&mut app, egui::Key::P, egui::Modifiers::NONE, None);
        assert_eq!(app.probe_mode, ProbeMode::ShortestPath);
        send_event(&mut app, egui::Key::L, egui::Modifiers::NONE, None);
        assert_eq!(app.probe_mode, ProbeMode::LesionAblate);
        send_event(&mut app, egui::Key::I, egui::Modifiers::NONE, None);
        assert_eq!(app.probe_mode, ProbeMode::Inspect);
        let before_cut = app.cutaway_mode;
        send_event(&mut app, egui::Key::C, egui::Modifiers::NONE, None);
        assert_ne!(app.cutaway_mode, before_cut);

        // Ctrl+K palette focus blocks single-key viewport hotkeys & Enter executes top match
        send_event(&mut app, egui::Key::K, egui::Modifiers::CTRL, None);
        assert!(app.cmd_palette_open);
        send_event(&mut app, egui::Key::P, egui::Modifiers::NONE, Some("broca"));
        assert_eq!(app.probe_mode, ProbeMode::Inspect, "Typing inside Ctrl+K must not switch probe mode");
        send_event(&mut app, egui::Key::Enter, egui::Modifiers::NONE, None);
        assert!(!app.cmd_palette_open);
        assert_eq!(app.selected_region.as_deref(), Some("LANG.Broca_L"));

        // Cached lesion_delta matches evaluate_lesion_impact after toggling a lesion
        app.atlas.toggle_tract_sever("Corpus Callosum");
        app.recompute();
        let expected = app.atlas.evaluate_lesion_impact(app.min_weight);
        assert!((app.lesion_delta.efficiency_retention_pct - expected.efficiency_retention_pct).abs() < 1e-9);
    }

    #[test]
    fn renders_all_color_modes_matrix_and_criticality_tda_plots_headless() {
        let ctx = egui::Context::default();
        let mut app = ConnectomeAtlasApp::new_headless(0);
        let modes = [
            ColorSchemeMode::FunctionalNetwork,
            ColorSchemeMode::DtiDirectionalRgb,
            ColorSchemeMode::BetweennessHeatmap,
            ColorSchemeMode::HarmonicModePsi,
            ColorSchemeMode::CriticalitySusceptibility,
            ColorSchemeMode::TdaBeta1Loops,
        ];
        let tabs = [
            RightTab::Inspector,
            RightTab::AdjacencyMatrix,
            RightTab::RichClubHubs,
            RightTab::CriticalityTda,
            RightTab::VaultIo,
        ];
        for (i, &mode) in modes.iter().enumerate() {
            app.color_mode = mode;
            app.right_tab = tabs[i % tabs.len()];
            let raw = egui::RawInput {
                screen_rect: Some(egui::Rect::from_min_size(
                    egui::Pos2::ZERO,
                    egui::vec2(1480.0, 920.0),
                )),
                ..Default::default()
            };
            let out = ctx.run(raw, |ctx| app.ui_step(ctx, false));
            assert!(!out.shapes.is_empty());
        }
    }

    #[test]
    fn command_palette_inspect_topology_digest_switches_to_hubs_tab() {
        let mut app = ConnectomeAtlasApp::new_headless(0);
        assert_eq!(app.topology_digest.entries.len(), 4);
        app.right_tab = RightTab::Inspector;
        app.execute_command_action(CommandAction::InspectTopologyDigest, false);
        assert_eq!(app.right_tab, RightTab::RichClubHubs);
        assert!(!app.zen_mode);
        assert!(app.show_right_panel);
    }

    #[test]
    fn recompute_refreshes_active_preset_in_topology_digest_after_lesion() {
        let mut app = ConnectomeAtlasApp::new_headless(0);
        let before_eff = app.topology_digest.entries[0].global_efficiency;
        app.atlas.toggle_tract_sever("Corpus Callosum");
        app.recompute();
        let after_eff = app.topology_digest.entries[0].global_efficiency;
        assert!(after_eff < before_eff);
    }
}
