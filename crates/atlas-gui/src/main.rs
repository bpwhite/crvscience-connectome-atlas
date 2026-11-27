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

#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]
#![deny(unsafe_code)]

mod app;
mod icon;

use app::{render_headless_screenshot, ConnectomeAtlasApp};
use atlas_core::{
    all_presets, compute_all_presets_topology_digest, load_preset_by_id, AtlasEnvelope,
};
use atlas_io::{format_swc_document, generate_region_swc_arbor};
use clap::{Parser, Subcommand};
use eframe::egui;
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(
    name = "connectome-atlas",
    version = "0.8.0",
    about = "Standalone 3D Brain & Multi-Scale Connectome Atlas Workstation"
)]
struct Cli {
    /// Path to local SQLite WAL vault
    #[arg(long, default_value = "connectome_vault.sqlite")]
    db: PathBuf,

    /// Initial connectome preset ID (human_mni152, mouse_ccfv3, drosophila_flywire, celegans_wormatlas)
    #[arg(long)]
    preset: Option<String>,

    /// Render an offscreen 1600x1000 RGBA PNG screenshot without opening a GUI window
    #[arg(long)]
    headless_screenshot: Option<PathBuf>,

    /// Select preset view (1=Human MNI152, 2=Drosophila FlyWire, 3=C. elegans WormAtlas, 4=Mouse CCFv3) for --headless-screenshot
    #[arg(long, default_value_t = 1)]
    screenshot_view: u8,

    #[command(subcommand)]
    command: Option<Commands>,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// List bundled multi-scale connectome presets
    Presets,
    /// Compute structural network summary metrics for a preset
    Metrics {
        #[arg(long, default_value = "human_mni152")]
        preset: String,
    },
    /// Compute brain criticality branching ratio and dynamic range report
    Criticality {
        #[arg(long, default_value = "human_mni152")]
        preset: String,
        #[arg(long, default_value_t = 1.0)]
        kappa: f64,
    },
    /// Compute clique complex and Z2 persistent homology (beta-1) report
    Topology {
        #[arg(long, default_value = "human_mni152")]
        preset: String,
        #[arg(long, default_value_t = 0.0)]
        min_weight: f64,
    },
    /// Output comparative graph-topology digest (sigma, C, L, rich-club Phi, modularity Q)
    TopologyDigest {
        #[arg(long, default_value_t = 0.0)]
        min_weight: f64,
        #[arg(long, default_value_t = 1337)]
        seed: u64,
        /// Output format: json or csv
        #[arg(long, default_value = "json")]
        format: String,
        #[arg(long)]
        out: Option<PathBuf>,
    },
    /// Compute weighted shortest signaling path between two regions
    Path {
        #[arg(long, default_value = "human_mni152")]
        preset: String,
        #[arg(long)]
        from: String,
        #[arg(long)]
        to: String,
        #[arg(long, default_value_t = 9.6)]
        velocity: f64,
    },
    /// Simulate structural lesion impact after severing a tract or ablating a region
    Lesion {
        #[arg(long, default_value = "human_mni152")]
        preset: String,
        #[arg(long)]
        tract: Option<String>,
        #[arg(long)]
        region: Option<String>,
    },
    /// Export a signed SHA-256 .connectome.json envelope
    ExportJson {
        #[arg(long, default_value = "human_mni152")]
        preset: String,
        #[arg(long)]
        out: PathBuf,
    },
    /// Export a sculpted 3D .swc morphology arbor for a region
    ExportSwc {
        #[arg(long, default_value = "human_mni152")]
        preset: String,
        #[arg(long)]
        region: String,
        #[arg(long)]
        out: PathBuf,
    },
}

#[cfg(windows)]
#[allow(unsafe_code)]
fn attach_parent_console_if_cli_args() {
    if std::env::args_os().len() > 1 {
        unsafe {
            let _ = windows_sys::Win32::System::Console::AttachConsole(u32::MAX);
        }
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    #[cfg(windows)]
    attach_parent_console_if_cli_args();

    let cli = Cli::parse();

    if let Some(out_path) = cli.headless_screenshot {
        render_headless_screenshot(&out_path, cli.screenshot_view)?;
        return Ok(());
    }

    if let Some(cmd) = cli.command {
        match cmd {
            Commands::Presets => {
                for p in all_presets() {
                    println!(
                        "{:<22} | {:<18} | {} regions, {} pathways | {}",
                        p.id,
                        p.organism,
                        p.regions.len(),
                        p.edges.len(),
                        p.coordinate_space
                    );
                }
            }
            Commands::Metrics { preset } => {
                let atlas = load_preset_by_id(&preset)
                    .ok_or_else(|| format!("Unknown preset '{preset}'"))?;
                let m = atlas.compute_summary_metrics(0.0);
                println!("{}", serde_json::to_string_pretty(&m)?);
            }
            Commands::Criticality { preset, kappa } => {
                let atlas = load_preset_by_id(&preset)
                    .ok_or_else(|| format!("Unknown preset '{preset}'"))?;
                let rep = atlas.compute_criticality_report(0.0, kappa);
                println!("{}", serde_json::to_string_pretty(&rep)?);
            }
            Commands::Topology { preset, min_weight } => {
                let atlas = load_preset_by_id(&preset)
                    .ok_or_else(|| format!("Unknown preset '{preset}'"))?;
                let tda = atlas.compute_topology_tda_report(min_weight);
                println!("{}", serde_json::to_string_pretty(&tda)?);
            }
            Commands::TopologyDigest {
                min_weight,
                seed,
                format,
                out,
            } => {
                let digest = compute_all_presets_topology_digest(min_weight, seed);
                let text = match format.to_ascii_lowercase().as_str() {
                    "csv" => digest.to_csv(),
                    "json" => digest.to_json_pretty()?,
                    other => {
                        return Err(format!("Unsupported format '{other}'; expected json or csv").into());
                    }
                };
                if let Some(out_path) = out {
                    std::fs::write(&out_path, &text)?;
                    println!("Wrote topology digest ({format}) to {}", out_path.display());
                } else {
                    println!("{text}");
                }
            }
            Commands::Path {
                preset,
                from,
                to,
                velocity,
            } => {
                let mut atlas = load_preset_by_id(&preset)
                    .ok_or_else(|| format!("Unknown preset '{preset}'"))?;
                atlas.conduction_velocity_m_s = velocity;
                match atlas.shortest_signaling_path(&from, &to, 0.0) {
                    Some(res) => println!("{}", serde_json::to_string_pretty(&res)?),
                    None => {
                        return Err(format!(
                            "No intact signaling path between '{from}' and '{to}'."
                        )
                        .into());
                    }
                }
            }
            Commands::Lesion {
                preset,
                tract,
                region,
            } => {
                let mut atlas = load_preset_by_id(&preset)
                    .ok_or_else(|| format!("Unknown preset '{preset}'"))?;
                if let Some(t) = tract {
                    atlas.toggle_tract_sever(&t);
                }
                if let Some(r) = region {
                    atlas.toggle_region_ablation(&r);
                }
                let delta = atlas.evaluate_lesion_impact(0.0);
                println!("{}", serde_json::to_string_pretty(&delta)?);
            }
            Commands::ExportJson { preset, out } => {
                let atlas = load_preset_by_id(&preset)
                    .ok_or_else(|| format!("Unknown preset '{preset}'"))?;
                let env = AtlasEnvelope::wrap(atlas)?;
                std::fs::write(&out, serde_json::to_string_pretty(&env)?)?;
                println!("Wrote {} (SHA-256: {})", out.display(), env.sha256_digest);
            }
            Commands::ExportSwc {
                preset,
                region,
                out,
            } => {
                let atlas = load_preset_by_id(&preset)
                    .ok_or_else(|| format!("Unknown preset '{preset}'"))?;
                let reg = atlas
                    .regions
                    .iter()
                    .find(|r| r.id == region)
                    .ok_or_else(|| format!("Unknown region '{region}' in {preset}"))?;
                let arbor = generate_region_swc_arbor(reg, &atlas);
                let doc = format_swc_document(reg, &atlas.organism, &arbor);
                std::fs::write(&out, doc)?;
                println!(
                    "Exported {} SWC compartments for {} to {}",
                    arbor.len(),
                    reg.id,
                    out.display()
                );
            }
        }
        return Ok(());
    }

    let db_path = cli.db;
    let initial_preset = cli.preset;
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_title("ConnectomeAtlas v0.8.0 - Standalone 3D Brain & Connectome Workstation")
            .with_icon(icon::app_egui_icon())
            .with_inner_size([1480.0, 920.0])
            .with_min_inner_size([1080.0, 700.0]),
        ..Default::default()
    };
    eframe::run_native(
        "ConnectomeAtlas v0.8.0",
        options,
        Box::new(move |cc| {
            Ok(Box::new(ConnectomeAtlasApp::new_with_db(
                cc,
                db_path,
                initial_preset.as_deref(),
            )))
        }),
    )
    .map_err(|e| format!("eframe error: {e}").into())
}
