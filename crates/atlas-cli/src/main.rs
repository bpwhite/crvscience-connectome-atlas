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

use atlas_core::{all_presets, build_command_index, load_preset_by_id, search_command_index, AtlasEnvelope};
use atlas_io::{format_swc_document, generate_region_swc_arbor, parse_nifti1_header, parse_trk_streamlines, AtlasVault};
use atlas_render::compute_atlas_fiber_tensors;
use clap::{Parser, Subcommand};
use std::fs;
use std::path::PathBuf;

#[derive(Parser)]
#[command(
    name = "atlas-cli",
    version,
    about = "Headless CLI for ConnectomeAtlas — Multi-Scale 3D Brain & Connectome Workstation"
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// List all bundled multi-scale reference connectome presets
    Presets,
    /// Quick-jump search across regions, functional roles, and workstation commands
    Search {
        #[arg(long, default_value = "human_mni152")]
        preset: String,
        #[arg(long)]
        query: String,
        #[arg(long, default_value_t = 10)]
        limit: usize,
    },
    /// Compute topological summary & rich-club hubs for a preset
    Inspect {
        #[arg(long, default_value = "human_mni152")]
        preset: String,
        #[arg(long, default_value_t = 0.0)]
        min_weight: f64,
    },
    /// Inspect a single anatomical module's functional role, afferent inputs, and efferent outputs
    InspectRegion {
        #[arg(long, default_value = "human_mni152")]
        preset: String,
        #[arg(long)]
        region: String,
    },
    /// Find the weighted Dijkstra shortest signaling path between two regions
    ShortestPath {
        #[arg(long, default_value = "human_mni152")]
        preset: String,
        #[arg(long)]
        from: String,
        #[arg(long)]
        to: String,
    },
    /// Run receptor-informed network control theory (average/modal controllability & transition energy)
    Control {
        #[arg(long, default_value = "human_mni152")]
        preset: String,
        #[arg(long, default_value = "ht")]
        receptor: String,
        #[arg(long, default_value_t = 1.2)]
        gain: f64,
    },
    /// Compute normalized Graph Laplacian connectome harmonics & Principal Gradient G1/G2
    Harmonics {
        #[arg(long, default_value = "human_mni152")]
        preset: String,
        #[arg(long, default_value_t = 5)]
        top_modes: usize,
    },
    /// Compute Brain Criticality branching ratio sigma, sigmoidal response & dynamic range Delta (dB)
    Criticality {
        #[arg(long, default_value = "human_mni152")]
        preset: String,
        #[arg(long, default_value_t = 1.0)]
        kappa: f64,
    },
    /// Compute Simplicial Clique Complex (3/4-cliques) & Z2 Persistent Homology beta-1 cavities
    Topology {
        #[arg(long, default_value = "human_mni152")]
        preset: String,
        #[arg(long, default_value_t = 0.0)]
        min_weight: f64,
    },
    /// Compute Synergy-Redundancy (PhiID), Metabolic Wiring Economy & Transcriptomic G1 correlation
    Synergy {
        #[arg(long, default_value = "human_mni152")]
        preset: String,
    },
    /// Run Maslov-Sneppen degree-preserving null ensemble & Betzel-Sporns generative wiring benchmark
    NullBenchmark {
        #[arg(long, default_value = "human_mni152")]
        preset: String,
        #[arg(long, default_value_t = 1337)]
        seed: u64,
    },
    /// Decode functional domain recruitment or lesion-network % disruption across cognitive/motor systems
    DecodeLesion {
        #[arg(long, default_value = "human_mni152")]
        preset: String,
        #[arg(long)]
        ablate: Option<String>,
        #[arg(long)]
        sever_tract: Option<String>,
    },
    /// Simulate virtual tractotomy or focal region ablation and report efficiency delta
    Lesion {
        #[arg(long, default_value = "human_mni152")]
        preset: String,
        #[arg(long)]
        tract: Option<String>,
        #[arg(long)]
        region: Option<String>,
    },
    /// Compute 3x3 crossing-fiber orientation tensors (FA, Westin cl/cp/cs & superquadric exponents)
    Tensors {
        #[arg(long, default_value = "human_mni152")]
        preset: String,
    },
    /// Export a standardized 7-column NeuroMorpho/INCF `.swc` 3D neuronal arbor for a module
    ExportSwc {
        #[arg(long, default_value = "drosophila_flywire")]
        preset: String,
        #[arg(long)]
        region: String,
        #[arg(long)]
        out: PathBuf,
    },
    /// Export a synchronized SHA-256 `.connectome.json` envelope
    ExportJson {
        #[arg(long, default_value = "human_mni152")]
        preset: String,
        #[arg(long)]
        out: PathBuf,
    },
    /// Inspect a NIfTI-1 (.nii or .nii.gz) volumetric parcellation header
    InspectNifti {
        #[arg(long)]
        file: PathBuf,
    },
    /// Inspect a TrackVis (.trk) diffusion tractography streamline file
    InspectTrk {
        #[arg(long)]
        file: PathBuf,
    },
    /// Save or list snapshots in the local SQLite WAL connectome vault
    Vault {
        #[arg(long, default_value = "connectome_vault.sqlite")]
        db: PathBuf,
        #[arg(long)]
        save_preset: Option<String>,
        #[arg(long, default_value = "default_slot")]
        slot: String,
    },
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();
    match cli.command {
        Commands::Search {
            preset,
            query,
            limit,
        } => {
            let atlas = load_preset_by_id(&preset)
                .ok_or_else(|| format!("Unknown preset '{preset}'"))?;
            let index = build_command_index(&atlas);
            let hits = search_command_index(&index, &query, limit);
            println!("Quick-Jump Command Search ('{query}' in {preset} — {} hits):", hits.len());
            for h in hits {
                println!("  [{:<10}] {:<38} | {}", h.category, h.title, h.subtitle);
            }
        }
        Commands::Presets => {
            for p in all_presets() {
                println!(
                    "{:<20} | {:<45} | {} ROIs | {} pathways",
                    p.id,
                    p.title,
                    p.regions.len(),
                    p.edges.len()
                );
            }
        }
        Commands::Inspect { preset, min_weight } => {
            let atlas = load_preset_by_id(&preset)
                .ok_or_else(|| format!("Unknown preset '{preset}'"))?;
            let m = atlas.compute_summary_metrics(min_weight);
            println!("Atlas: {} ({})", atlas.title, atlas.organism);
            println!("Space: {}", atlas.coordinate_space);
            println!(
                "Nodes: {}/{} | Pathways: {}/{} | Components: {}",
                m.active_nodes, m.total_nodes, m.active_edges, m.total_edges, m.connected_components
            );
            println!(
                "Global Efficiency: {:.4} | Characteristic Path Cost: {:.4} | Mean Clustering: {:.4} | Rich-Club Phi: {:.4}",
                m.global_efficiency,
                m.characteristic_path_length,
                m.mean_clustering,
                m.rich_club_coeff_top_quartile
            );
            println!("
Top Rich-Club Hubs:");
            let mut hubs: Vec<_> = m.node_metrics.iter().filter(|n| n.is_rich_club_hub).collect();
            hubs.sort_by(|a, b| b.strength_total.total_cmp(&a.strength_total));
            for h in hubs.iter().take(8) {
                println!(
                    "  - {:<18} deg={:<2} strength={:.2} part={:.3} betw={:.3}",
                    h.region_id,
                    h.degree_total,
                    h.strength_total,
                    h.participation_coeff,
                    h.betweenness_centrality
                );
            }
        }
        Commands::InspectRegion { preset, region } => {
            let atlas = load_preset_by_id(&preset)
                .ok_or_else(|| format!("Unknown preset '{preset}'"))?;
            let reg = atlas
                .regions
                .iter()
                .find(|r| r.id == region)
                .ok_or_else(|| format!("Unknown region '{region}' in {preset}"))?;
            let (inc, out, coup) = atlas.incident_io_breakdown(&region, 0.0);
            println!("{} — {} [{}]", reg.id, reg.name, reg.network);
            println!("Coord: ({:.1}, {:.1}, {:.1})", reg.coord[0], reg.coord[1], reg.coord[2]);
            println!("What It Does: {}", reg.what_it_does);
            println!(
                "Receptors: DA={:.2} 5-HT={:.2} ACh={:.2} Glu={:.2} GABA={:.2}",
                reg.receptors.da, reg.receptors.serotonin, reg.receptors.ach, reg.receptors.glutamate, reg.receptors.gaba
            );
            println!("Afferent Inputs ({}):", inc.len());
            for e in inc {
                println!("  <- {:<18} w={:.2} ({})", e.source, e.weight, e.tract_name);
            }
            println!("Efferent Outputs ({}):", out.len());
            for e in out {
                println!("  -> {:<18} w={:.2} ({})", e.target, e.weight, e.tract_name);
            }
            println!("Reciprocal / Coupled Pathways ({}):", coup.len());
            for e in coup {
                let peer = if e.source == region { &e.target } else { &e.source };
                println!("  <-> {:<17} w={:.2} ({})", peer, e.weight, e.tract_name);
            }
        }
        Commands::ShortestPath { preset, from, to } => {
            let atlas = load_preset_by_id(&preset)
                .ok_or_else(|| format!("Unknown preset '{preset}'"))?;
            match atlas.shortest_signaling_path(&from, &to, 0.0) {
                Some(res) => {
                    println!("Signaling Path: {}", res.hops.join(" -> "));
                    println!(
                        "Hops: {} | Latency: {:.2} ms ({:.1} mm @ 9.6 m/s) | Weighted Cost: {:.4} | Bottleneck Weight: {:.3}",
                        res.hops.len().saturating_sub(1),
                        res.conduction_latency_ms,
                        res.arc_length_mm,
                        res.total_distance_cost,
                        res.bottleneck_weight
                    );
                }
                None => println!("No intact signaling path between '{from}' and '{to}'."),
            }
        }
        Commands::Control {
            preset,
            receptor,
            gain,
        } => {
            let atlas = load_preset_by_id(&preset)
                .ok_or_else(|| format!("Unknown preset '{preset}'"))?;
            let rep = atlas.compute_control_and_harmonics(0.0, &receptor, gain);
            println!(
                "Receptor-Informed Control ({preset}, channel={receptor}, alpha={gain:+.2}):"
            );
            println!(
                "Mean Transition Energy: {:.4} (Placebo) -> {:.4} (Modulated) [{:+.2}% Landscape Shift]",
                rep.placebo_mean_transition_energy,
                rep.modulated_mean_transition_energy,
                rep.energy_landscape_delta_pct
            );
            let mut by_avg = rep.node_profiles.clone();
            by_avg.sort_by(|a, b| b.average_controllability.total_cmp(&a.average_controllability));
            println!("Top Average Controllability Drivers:");
            for p in by_avg.iter().take(5) {
                println!(
                    "  - {:<18} avg_ctrl={:.3} modal_ctrl={:.3} G1={:+.3}",
                    p.region_id, p.average_controllability, p.modal_controllability, p.principal_gradient_g1
                );
            }
        }
        Commands::Harmonics { preset, top_modes } => {
            let atlas = load_preset_by_id(&preset)
                .ok_or_else(|| format!("Unknown preset '{preset}'"))?;
            let rep = atlas.compute_control_and_harmonics(0.0, "ht", 0.0);
            println!("Connectome Laplacian Harmonics ({preset}):");
            for (k, &omega) in rep.laplacian_wavenumbers.iter().take(top_modes).enumerate() {
                println!("  Mode psi_{} | wavenumber omega = {:.4}", k + 1, omega);
            }
        }
        Commands::Criticality { preset, kappa } => {
            let atlas = load_preset_by_id(&preset)
                .ok_or_else(|| format!("Unknown preset '{preset}'"))?;
            let rep = atlas.compute_criticality_report(0.0, kappa);
            println!(
                "Brain Criticality Observatory ({preset}, kappa={:.2}):",
                rep.coupling_kappa
            );
            println!(
                "  Branching Ratio sigma = {:.4} [{}] | Dynamic Range Delta = {:.2} dB | Spont. Floor F0 = {:.4}",
                rep.branching_ratio_sigma,
                rep.regime.label(),
                rep.dynamic_range_db,
                rep.spontaneous_floor_f0
            );
            let mut by_sus = rep.node_profiles.clone();
            by_sus.sort_by(|a, b| b.avalanche_susceptibility.total_cmp(&a.avalanche_susceptibility));
            println!("  Top Critical Avalanche Driver Modules:");
            for p in by_sus.iter().take(6) {
                println!(
                    "    - {:<18} susceptibility={:.3} E/I_gain={:.2}",
                    p.region_id, p.avalanche_susceptibility, p.local_ei_gain
                );
            }
        }
        Commands::Topology { preset, min_weight } => {
            let atlas = load_preset_by_id(&preset)
                .ok_or_else(|| format!("Unknown preset '{preset}'"))?;
            let tda = atlas.compute_topology_tda_report(min_weight);
            println!(
                "Algebraic Topology Report ({preset}, min_weight={min_weight:.2}):"
            );
            println!(
                "  3-Cliques (triangles) = {} | 4-Cliques (tetrahedra) = {} | Max Simplex Dim = {}",
                tda.clique_3_total, tda.clique_4_total, tda.max_clique_dimension
            );
            println!(
                "  Persistent beta-1 Cavities = {} ({} essential) | Max Persistence pi = {:.3}",
                tda.persistent_beta1_count, tda.essential_beta1_count, tda.max_persistence
            );
            for cav in tda.cavities.iter().take(6) {
                println!(
                    "    - gamma_{:<2} [{:.2} -> {:.2}{}] pi={:.3} | {}",
                    cav.cavity_index,
                    cav.birth_weight,
                    cav.death_weight,
                    if cav.is_essential { " inf" } else { "" },
                    cav.persistence,
                    cav.cycle_nodes.join(" -> ")
                );
            }
        }
        Commands::Synergy { preset } => {
            let atlas = load_preset_by_id(&preset)
                .ok_or_else(|| format!("Unknown preset '{preset}'"))?;
            let rep = atlas.compute_systems_neuro_report(0.0, 1337);
            println!(
                "Systems Neuroscience Report ({preset}) | Wiring Cost: {:.1} mm*w | Economy E/C: {:.3} | |r(SST-PVALB, G1)| = {:.3}",
                rep.total_wiring_cost_mm,
                rep.economy_efficiency_ratio,
                rep.transcriptomic_g1_pearson_r
            );
            let mut by_syn = rep.node_profiles.clone();
            by_syn.sort_by(|a, b| b.synergy_rank_gradient.total_cmp(&a.synergy_rank_gradient));
            println!("Top Synergistic Workspace Hubs:");
            for p in by_syn.iter().take(6) {
                println!(
                    "  - {:<18} syn_rank={:+.3} AG={:.2} SST={:.2} PVALB={:.2}",
                    p.region_id,
                    p.synergy_rank_gradient,
                    p.aerobic_glycolysis_index,
                    p.transcriptomics.sst,
                    p.transcriptomics.pvalb
                );
            }
        }
        Commands::NullBenchmark { preset, seed } => {
            let atlas = load_preset_by_id(&preset)
                .ok_or_else(|| format!("Unknown preset '{preset}'"))?;
            let nb = atlas.compute_systems_neuro_report(0.0, seed).null_benchmark;
            println!("Maslov-Sneppen & Betzel-Sporns Generative Null Benchmark ({preset}):");
            println!(
                "  Small-Worldness sigma = {:.3}x | Rich-Club Phi_norm = {:.3}x | Wiring Cost vs Null = {:+.2}%",
                nb.small_world_sigma, nb.rich_club_normalized, nb.wiring_cost_saving_pct
            );
            println!(
                "  Generative Fit (E_ij ~ D^-eta * K^gamma): eta* = {:.2}, gamma* = {:.2} (Discrepancy = {:.3})",
                nb.generative_eta_distance, nb.generative_gamma_homophily, nb.generative_discrepancy
            );
        }
        Commands::DecodeLesion {
            preset,
            ablate,
            sever_tract,
        } => {
            let mut atlas = load_preset_by_id(&preset)
                .ok_or_else(|| format!("Unknown preset '{preset}'"))?;
            if let Some(r) = &ablate {
                atlas.toggle_region_ablation(r);
            }
            if let Some(t) = &sever_tract {
                atlas.toggle_tract_sever(t);
            }
            let doms = atlas.decode_functional_domains(ablate.as_deref(), 0.0);
            println!("Functional Domain / Lesion-Network Decoding ({preset}):");
            for d in doms {
                if d.is_lesion_disruption {
                    println!("  - {:<28} -{:.1}% network disruption", d.domain_label, d.score * 100.0);
                } else {
                    println!("  - {:<28} recruitment = {:.2}", d.domain_label, d.score);
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
                println!("Severed Tract: {t}");
            }
            if let Some(r) = region {
                atlas.toggle_region_ablation(&r);
                println!("Ablated Region: {r}");
            }
            let d = atlas.evaluate_lesion_impact(0.0);
            println!(
                "Global Efficiency: {:.4} (Intact) -> {:.4} (Lesioned) [{:.2}% Retained]",
                d.intact_efficiency, d.lesioned_efficiency, d.efficiency_retention_pct
            );
            println!(
                "Connected Components: {} -> {}",
                d.intact_components, d.lesioned_components
            );
        }
        Commands::Tensors { preset } => {
            let atlas = load_preset_by_id(&preset)
                .ok_or_else(|| format!("Unknown preset '{preset}'"))?;
            let mut tensors = compute_atlas_fiber_tensors(&atlas, 0.0);
            tensors.sort_by(|a, b| b.westin_planar_cp.total_cmp(&a.westin_planar_cp));
            println!("Local 3D Crossing-Fiber Orientation Tensors ({preset}):");
            for t in tensors.iter().take(10) {
                println!(
                    "  - {:<18} FA={:.3} cl={:.3} cp={:.3} cs={:.3} sq(a={:.2}, b={:.2})",
                    t.region_id,
                    t.fractional_anisotropy,
                    t.westin_linear_cl,
                    t.westin_planar_cp,
                    t.westin_spherical_cs,
                    t.superquadric_alpha_beta[0],
                    t.superquadric_alpha_beta[1]
                );
            }
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
            fs::write(&out, doc)?;
            println!(
                "Exported {} SWC compartments for {} to {}",
                arbor.len(),
                reg.id,
                out.display()
            );
        }
        Commands::ExportJson { preset, out } => {
            let atlas = load_preset_by_id(&preset)
                .ok_or_else(|| format!("Unknown preset '{preset}'"))?;
            let env = AtlasEnvelope::wrap(atlas)?;
            fs::write(&out, serde_json::to_string_pretty(&env)?)?;
            println!("Wrote {} (SHA-256: {})", out.display(), env.sha256_digest);
        }
        Commands::InspectNifti { file } => {
            let bytes = fs::read(&file)?;
            let hdr = parse_nifti1_header(&bytes)?;
            println!("{}", serde_json::to_string_pretty(&hdr)?);
        }
        Commands::InspectTrk { file } => {
            let bytes = fs::read(&file)?;
            let (hdr, tracks) = parse_trk_streamlines(&bytes, 500)?;
            println!("TrackVis Header: {}", serde_json::to_string_pretty(&hdr)?);
            println!("Sampled Streamlines: {}", tracks.len());
        }
        Commands::Vault {
            db,
            save_preset,
            slot,
        } => {
            let vault = AtlasVault::open(&db)?;
            if let Some(pid) = save_preset {
                let atlas = load_preset_by_id(&pid)
                    .ok_or_else(|| format!("Unknown preset '{pid}'"))?;
                let sha = vault.save_snapshot(&slot, &atlas)?;
                println!("Saved snapshot '{slot}' (SHA-256: {sha}) to {}", db.display());
            }
            for item in vault.list_snapshots()? {
                println!(
                    "#{} [{}] {} ({}) sha256={}.. saved={}",
                    item.id,
                    item.slot_name,
                    item.atlas_id,
                    item.organism,
                    &item.sha256_digest[..12],
                    item.saved_at_utc
                );
            }
        }
    }
    Ok(())
}
