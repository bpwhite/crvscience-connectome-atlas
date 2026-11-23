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

use crate::model::{BrainRegion, ConnectomeAtlas};
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct TranscriptomicProfile {
    pub pvalb: f64,
    pub sst: f64,
    pub mbp: f64,
    pub grin2b: f64,
    pub sst_minus_pvalb_axis: f64,
}

pub fn region_transcriptomic_profile(preset_id: &str, region: &BrainRegion) -> TranscriptomicProfile {
    let id = region.id.as_str();
    let net = region.network.as_str();
    let (pvalb, sst, mbp, grin2b) = if preset_id == "human_mni152" || preset_id == "mouse_ccfv3" {
        if net.contains("Visual") || id.contains("VIS") {
            (0.93, 0.28, 0.89, 0.31)
        } else if net.contains("Sensorimotor") || id.contains("SMN") || id.contains("MOp") || id.contains("SSp") {
            (0.88, 0.34, 0.91, 0.36)
        } else if net.contains("Default") || id.contains("DMN") || id.contains("ACA") || id.contains("RSP") {
            (0.31, 0.91, 0.38, 0.89)
        } else if net.contains("Limbic") || id.contains("Hipp") || id.contains("Amyg") || id.contains("CA1") {
            (0.26, 0.94, 0.32, 0.92)
        } else if net.contains("Frontoparietal") || net.contains("Language") || net.contains("Salience") {
            (0.41, 0.84, 0.47, 0.85)
        } else if net.contains("Cerebell") || id.contains("CB.") {
            (0.95, 0.24, 0.84, 0.29)
        } else {
            (0.64, 0.52, 0.68, 0.55)
        }
    } else if preset_id == "drosophila_flywire" {
        if net.contains("Optic") || id.contains("ME_") || id.contains("LO_") {
            (0.91, 0.29, 0.84, 0.35)
        } else if net.contains("Mushroom") || id.contains("MB") {
            (0.28, 0.93, 0.42, 0.88)
        } else if net.contains("Central") || id == "EB" || id == "FB" || id == "PB" {
            (0.46, 0.84, 0.58, 0.82)
        } else {
            (0.62, 0.60, 0.62, 0.60)
        }
    } else if net.contains("Motor") || id.starts_with("DA0") || id.starts_with("DB0") || id.starts_with("VB0") {
        (0.89, 0.30, 0.82, 0.32)
    } else if id.starts_with("AVA") || id.starts_with("AVB") || id.starts_with("PVC") || id == "RIBL" {
        (0.34, 0.90, 0.52, 0.94)
    } else {
        (0.55, 0.68, 0.54, 0.66)
    };
    TranscriptomicProfile {
        pvalb,
        sst,
        mbp,
        grin2b,
        sst_minus_pvalb_axis: (0.5 + (sst - pvalb) * 0.62).clamp(0.0, 1.0),
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NodeSystemsProfile {
    pub region_id: String,
    pub synergy_rank_gradient: f64,
    pub mean_synergy: f64,
    pub mean_redundancy: f64,
    pub wiring_cost_mm: f64,
    pub aerobic_glycolysis_index: f64,
    pub transcriptomics: TranscriptomicProfile,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NullGenerativeBenchmark {
    pub small_world_sigma: f64,
    pub rich_club_normalized: f64,
    pub wiring_cost_saving_pct: f64,
    pub generative_eta_distance: f64,
    pub generative_gamma_homophily: f64,
    pub generative_discrepancy: f64,
    #[serde(default)]
    pub obs_characteristic_path_length: f64,
    #[serde(default)]
    pub null_characteristic_path_length: f64,
    #[serde(default)]
    pub null_clustering_coeff: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DomainDecodeEntry {
    pub domain_label: String,
    pub score: f64,
    pub is_lesion_disruption: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SystemsNeuroReport {
    pub total_wiring_cost_mm: f64,
    pub economy_efficiency_ratio: f64,
    pub transcriptomic_g1_pearson_r: f64,
    pub null_benchmark: NullGenerativeBenchmark,
    pub node_profiles: Vec<NodeSystemsProfile>,
}

fn dist_mm(a: [f32; 3], b: [f32; 3]) -> f64 {
    let dx = (b[0] - a[0]) as f64;
    let dy = (b[1] - a[1]) as f64;
    let dz = (b[2] - a[2]) as f64;
    ((dx * dx + dy * dy + dz * dz).sqrt() * 1.15).max(2.0)
}

fn mean_weighted_path_length(
    n: usize,
    active_idxs: &[usize],
    edges: &[(usize, usize, f64)],
) -> (f64, bool) {
    let mut adj: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
    for &(u, v, w) in edges {
        if u != v && w > 0.0 {
            adj[u].push((v, w));
            adj[v].push((u, w));
        }
    }
    let mut finite_dist_sum = 0.0_f64;
    let mut finite_pairs = 0_usize;
    for &s in active_idxs {
        let (dist, _, _, _) = crate::metrics::dijkstra_brandes_from(s, n, &adj);
        for &t in active_idxs {
            if s == t {
                continue;
            }
            let d = dist[t];
            if d.is_finite() && d > 0.0 {
                finite_dist_sum += d;
                finite_pairs += 1;
            }
        }
    }
    let expected_pairs = active_idxs.len() * active_idxs.len().saturating_sub(1);
    let l = if finite_pairs > 0 {
        finite_dist_sum / (finite_pairs as f64)
    } else {
        0.0
    };
    (l, expected_pairs > 0 && finite_pairs == expected_pairs)
}

impl ConnectomeAtlas {
    /// Computes the v0.7.0 Systems Neuroscience suite:
    /// 1) Synergy–Redundancy (`ΦID` proxy) rank gradient
    /// 2) Metabolic Wiring Economy & Aerobic Glycolysis index
    /// 3) Normative `PVALB`/`SST` Transcriptomic correlation against Principal Gradient `G1`
    /// 4) Maslov–Sneppen degree-preserving null benchmark & Betzel–Sporns generative parameters
    pub fn compute_systems_neuro_report(&self, min_weight: f64, seed: u64) -> SystemsNeuroReport {
        let n = self.regions.len();
        let summary = self.compute_summary_metrics(min_weight);
        let ctrl_harm = self.compute_control_and_harmonics(min_weight, "ht", 0.0);

        let mut idx_of = HashMap::new();
        for (i, r) in self.regions.iter().enumerate() {
            idx_of.insert(r.id.as_str(), i);
        }

        let mut sym_w = vec![vec![0.0_f64; n]; n];
        let mut node_cost = vec![0.0_f64; n];
        let mut total_cost = 0.0_f64;
        let mut undir_edges: Vec<(usize, usize, f64)> = Vec::new();
        let mut adj_sets: Vec<HashSet<usize>> = vec![HashSet::new(); n];

        for e in &self.edges {
            if !self.is_edge_intact(e, min_weight) {
                continue;
            }
            let (Some(&u), Some(&v)) = (idx_of.get(e.source.as_str()), idx_of.get(e.target.as_str())) else {
                continue;
            };
            if u == v {
                continue;
            }
            let length = dist_mm(self.regions[u].coord, self.regions[v].coord);
            let c_mm = e.weight * length;
            total_cost += c_mm;
            node_cost[u] += c_mm;
            node_cost[v] += c_mm;
            if e.weight > sym_w[u][v] {
                sym_w[u][v] = e.weight;
                sym_w[v][u] = e.weight;
            }
            let (a, b) = if u < v { (u, v) } else { (v, u) };
            if adj_sets[a].insert(b) {
                adj_sets[b].insert(a);
                undir_edges.push((a, b, e.weight));
            }
        }

        let row_norm: Vec<f64> = (0..n)
            .map(|i| sym_w[i].iter().map(|w| w * w).sum::<f64>().sqrt())
            .collect();

        // Synergy vs Redundancy decomposition
        let mut mean_syn = vec![0.0_f64; n];
        let mut mean_red = vec![0.0_f64; n];
        for i in 0..n {
            if self.ablated_regions.contains(&self.regions[i].id) {
                continue;
            }
            let mut s_acc = 0.0;
            let mut r_acc = 0.0;
            let mut cnt = 0_usize;
            for j in 0..n {
                if i == j || self.ablated_regions.contains(&self.regions[j].id) {
                    continue;
                }
                let dot: f64 = (0..n).map(|k| sym_w[i][k] * sym_w[j][k]).sum();
                let cos_sc = if row_norm[i] > 1e-6 && row_norm[j] > 1e-6 {
                    dot / (row_norm[i] * row_norm[j])
                } else {
                    0.0
                };
                let same_net = self.regions[i].network == self.regions[j].network;
                let direct = sym_w[i][j];
                let fc_proxy = (0.45 * direct + 0.55 * cos_sc).clamp(0.05, 1.0);
                let rij = fc_proxy.min(direct * 0.7 + cos_sc * 0.6) * if same_net { 1.22 } else { 0.68 };
                let sij = (fc_proxy * (1.0 - 0.58 * cos_sc)).max(0.0) * if same_net { 0.48 } else { 1.28 };
                s_acc += sij;
                r_acc += rij;
                cnt += 1;
            }
            if cnt > 0 {
                mean_syn[i] = s_acc / (cnt as f64);
                mean_red[i] = r_acc / (cnt as f64);
            }
        }

        let active_idxs: Vec<usize> = (0..n)
            .filter(|&i| !self.ablated_regions.contains(&self.regions[i].id))
            .collect();
        let rank_vec = |vals: &[f64]| -> Vec<f64> {
            let mut order = active_idxs.clone();
            order.sort_by(|&a, &b| vals[a].total_cmp(&vals[b]));
            let mut out = vec![0.0_f64; n];
            let denom = (order.len().saturating_sub(1)).max(1) as f64;
            for (pos, &idx) in order.iter().enumerate() {
                out[idx] = (pos as f64) / denom;
            }
            out
        };
        let syn_r = rank_vec(&mean_syn);
        let red_r = rank_vec(&mean_red);

        let max_node_cost = node_cost.iter().copied().fold(1.0_f64, f64::max);
        let mut max_possible_cost = 1.0_f64;
        for i in 0..n {
            for j in (i + 1)..n {
                max_possible_cost += dist_mm(self.regions[i].coord, self.regions[j].coord);
            }
        }
        let cost_frac = total_cost / max_possible_cost;
        let economy_ratio = summary.global_efficiency / (cost_frac * 6.5).max(0.04);

        let mut node_profiles = Vec::with_capacity(n);
        let mut g_vals = Vec::new();
        let mut g1_vals = Vec::new();

        for i in 0..n {
            let reg = &self.regions[i];
            let trans = region_transcriptomic_profile(&self.id, reg);
            let sr = if self.ablated_regions.contains(&reg.id) {
                0.0
            } else {
                syn_r[i] - red_r[i]
            };
            let c_norm = node_cost[i] / max_node_cost;
            let ag = if self.ablated_regions.contains(&reg.id) {
                0.0
            } else {
                (0.34 * c_norm
                    + 0.36 * (reg.receptors.glutamate as f64) * trans.sst
                    + 0.30 * ((sr + 1.0) * 0.5))
                    .clamp(0.0, 1.0)
            };
            if !self.ablated_regions.contains(&reg.id) {
                if let Some(cp) = ctrl_harm.node_profiles.iter().find(|p| p.region_id == reg.id) {
                    g_vals.push(trans.sst_minus_pvalb_axis);
                    g1_vals.push(cp.principal_gradient_g1);
                }
            }
            node_profiles.push(NodeSystemsProfile {
                region_id: reg.id.clone(),
                synergy_rank_gradient: sr,
                mean_synergy: mean_syn[i],
                mean_redundancy: mean_red[i],
                wiring_cost_mm: node_cost[i],
                aerobic_glycolysis_index: ag,
                transcriptomics: trans,
            });
        }

        // Pearson correlation |r(SST-PVALB, G1)| (aligned with transmodal pole)
        let transcriptomic_g1_pearson_r = if g_vals.len() > 2 {
            let mg = g_vals.iter().sum::<f64>() / (g_vals.len() as f64);
            let mh = g1_vals.iter().sum::<f64>() / (g1_vals.len() as f64);
            let mut num = 0.0;
            let mut dg = 0.0;
            let mut dh = 0.0;
            for (&gv, &hv) in g_vals.iter().zip(g1_vals.iter()) {
                let a = gv - mg;
                let b = hv - mh;
                num += a * b;
                dg += a * a;
                dh += b * b;
            }
            if dg > 1e-9 && dh > 1e-9 {
                (num / (dg * dh).sqrt()).abs()
            } else {
                0.0
            }
        } else {
            0.0
        };

        let null_benchmark = self.run_null_and_generative_benchmark(
            &undir_edges,
            &adj_sets,
            &active_idxs,
            summary.mean_clustering,
            summary.rich_club_coeff_top_quartile,
            total_cost,
            seed,
        );

        SystemsNeuroReport {
            total_wiring_cost_mm: total_cost,
            economy_efficiency_ratio: economy_ratio,
            transcriptomic_g1_pearson_r,
            null_benchmark,
            node_profiles,
        }
    }

    fn run_null_and_generative_benchmark(
        &self,
        undir_edges: &[(usize, usize, f64)],
        adj_sets: &[HashSet<usize>],
        active_idxs: &[usize],
        obs_clust: f64,
        obs_rc: f64,
        obs_cost: f64,
        seed: u64,
    ) -> NullGenerativeBenchmark {
        let n = self.regions.len();
        if undir_edges.len() < 4 || active_idxs.len() < 4 {
            return NullGenerativeBenchmark {
                small_world_sigma: 1.0,
                rich_club_normalized: 1.0,
                wiring_cost_saving_pct: 0.0,
                generative_eta_distance: 1.6,
                generative_gamma_homophily: 0.8,
                generative_discrepancy: 0.18,
                obs_characteristic_path_length: 0.0,
                null_characteristic_path_length: 0.0,
                null_clustering_coeff: 0.0,
            };
        }

        let (obs_l, obs_conn) = mean_weighted_path_length(n, active_idxs, undir_edges);

        let mut rng = seed.max(1);
        let mut next_u64 = || -> u64 {
            rng ^= rng << 13;
            rng ^= rng >> 7;
            rng ^= rng << 17;
            rng
        };

        let mut degs: Vec<usize> = active_idxs.iter().map(|&i| adj_sets[i].len()).collect();
        degs.sort_unstable();
        let rc_thresh = degs[(degs.len() * 3) / 4].max(2);

        let n_surr = 10_usize;
        let swaps = undir_edges.len() * 4;
        let mut sum_c0 = 0.0_f64;
        let mut sum_l0_conn = 0.0_f64;
        let mut cnt_l0_conn = 0_usize;
        let mut sum_l0_all = 0.0_f64;
        let mut sum_rc0 = 0.0_f64;
        let mut sum_cost0 = 0.0_f64;

        for _ in 0..n_surr {
            let mut elist = undir_edges.to_vec();
            let mut adj = adj_sets.to_vec();
            for _ in 0..swaps {
                let i1 = (next_u64() as usize) % elist.len();
                let i2 = (next_u64() as usize) % elist.len();
                if i1 == i2 {
                    continue;
                }
                let (u, v, w1) = elist[i1];
                let (x, y, w2) = elist[i2];
                if u == x || u == y || v == x || v == y {
                    continue;
                }
                if adj[u].contains(&y) || adj[x].contains(&v) {
                    continue;
                }
                adj[u].remove(&v);
                adj[v].remove(&u);
                adj[x].remove(&y);
                adj[y].remove(&x);
                adj[u].insert(y);
                adj[y].insert(u);
                adj[x].insert(v);
                adj[v].insert(x);
                elist[i1] = (u, y, w1);
                elist[i2] = (x, v, w2);
            }

            let mut c_acc = 0.0;
            for &u in active_idxs {
                let nbrs: Vec<usize> = adj[u].iter().copied().collect();
                let k = nbrs.len();
                if k >= 2 {
                    let mut tri = 0_usize;
                    for a in 0..k {
                        for b in (a + 1)..k {
                            if adj[nbrs[a]].contains(&nbrs[b]) {
                                tri += 1;
                            }
                        }
                    }
                    c_acc += (2 * tri) as f64 / ((k * (k - 1)) as f64);
                }
            }
            sum_c0 += c_acc / (active_idxs.len() as f64);

            let mut c_mm = 0.0;
            for &(u, v, w) in &elist {
                c_mm += w * dist_mm(self.regions[u].coord, self.regions[v].coord);
            }
            sum_cost0 += c_mm;

            let hubs: Vec<usize> = active_idxs
                .iter()
                .copied()
                .filter(|&i| adj[i].len() >= rc_thresh)
                .collect();
            if hubs.len() > 1 {
                let mut internal = 0_usize;
                for a in 0..hubs.len() {
                    for b in (a + 1)..hubs.len() {
                        if adj[hubs[a]].contains(&hubs[b]) {
                            internal += 1;
                        }
                    }
                }
                sum_rc0 += (2 * internal) as f64 / ((hubs.len() * (hubs.len() - 1)) as f64);
            }

            let (l0, l0_conn) = mean_weighted_path_length(n, active_idxs, &elist);
            sum_l0_all += l0;
            if l0_conn {
                sum_l0_conn += l0;
                cnt_l0_conn += 1;
            }
        }

        let mean_c0 = (sum_c0 / (n_surr as f64)).max(1e-4);
        let mean_l0 = if obs_conn && cnt_l0_conn > 0 {
            sum_l0_conn / (cnt_l0_conn as f64)
        } else {
            sum_l0_all / (n_surr as f64)
        }
        .max(1e-6);
        let mean_rc0 = (sum_rc0 / (n_surr as f64)).max(0.05);
        let mean_cost0 = (sum_cost0 / (n_surr as f64)).max(1.0);

        let l_ratio = (obs_l.max(1e-6) / mean_l0).max(1e-6);
        let small_world_sigma = (obs_clust / mean_c0) / l_ratio;
        let rich_club_normalized = if obs_rc > 0.0 { obs_rc / mean_rc0 } else { 1.0 };
        let wiring_cost_saving_pct = ((obs_cost - mean_cost0) / mean_cost0) * 100.0;

        // Betzel-Sporns (2016) generative wiring fit: E_ij ~ D_ij^{-eta} * K_ij^{gamma}
        let mut hom = vec![vec![0.08_f64; n]; n];
        for i in 0..n {
            for j in (i + 1)..n {
                let shared = adj_sets[i].intersection(&adj_sets[j]).count();
                let denom = (adj_sets[i].len() + adj_sets[j].len()).max(1) as f64;
                let k_ij = 0.08 + (2.0 * (shared as f64)) / denom;
                hom[i][j] = k_ij;
                hom[j][i] = k_ij;
            }
        }
        let m_edges = undir_edges.len();
        let obs_mean_d = undir_edges
            .iter()
            .map(|&(u, v, _)| dist_mm(self.regions[u].coord, self.regions[v].coord))
            .sum::<f64>()
            / (m_edges as f64);
        let obs_mean_k = undir_edges
            .iter()
            .map(|&(u, v, _)| hom[u][v])
            .sum::<f64>()
            / (m_edges as f64);

        let mut best_eta = 1.6_f64;
        let mut best_gamma = 0.9_f64;
        let mut best_err = f64::MAX;

        for &eta in &[0.6_f64, 1.1, 1.6, 2.1, 2.7] {
            for &gam in &[0.2_f64, 0.6, 1.0, 1.5, 2.0] {
                let mut cands: Vec<(f64, f64, f64)> = Vec::new();
                for a in 0..active_idxs.len() {
                    let u = active_idxs[a];
                    for b in (a + 1)..active_idxs.len() {
                        let v = active_idxs[b];
                        let d = dist_mm(self.regions[u].coord, self.regions[v].coord);
                        let k = hom[u][v];
                        let s = d.powf(-eta) * k.powf(gam);
                        cands.push((s, d, k));
                    }
                }
                cands.sort_by(|a, b| b.0.total_cmp(&a.0));
                let take_n = cands.len().min(m_edges).max(1);
                let sim_d = cands[..take_n].iter().map(|c| c.1).sum::<f64>() / (take_n as f64);
                let sim_k = cands[..take_n].iter().map(|c| c.2).sum::<f64>() / (take_n as f64);
                let err = (sim_d - obs_mean_d).abs() / obs_mean_d.max(10.0)
                    + (sim_k - obs_mean_k).abs() / obs_mean_k.max(0.1);
                if err < best_err {
                    best_err = err;
                    best_eta = eta;
                    best_gamma = gam;
                }
            }
        }

        NullGenerativeBenchmark {
            small_world_sigma,
            rich_club_normalized,
            wiring_cost_saving_pct,
            generative_eta_distance: best_eta,
            generative_gamma_homophily: best_gamma,
            generative_discrepancy: (best_err * 0.55).clamp(0.06, 0.42),
            obs_characteristic_path_length: obs_l,
            null_characteristic_path_length: mean_l0,
            null_clustering_coeff: mean_c0,
        }
    }

    /// Decodes either single-region functional domain recruitment (when intact) or
    /// lesion-network `% domain disruption` (when focal ablations or tractotomies are active).
    pub fn decode_functional_domains(
        &self,
        selected_region_id: Option<&str>,
        min_weight: f64,
    ) -> Vec<DomainDecodeEntry> {
        let domains: Vec<(&str, Vec<&str>)> = if self.id == "drosophila_flywire" {
            vec![
                ("Associative Memory (MB)", vec!["Mushroom", "MB"]),
                ("Compass Navigation (CX)", vec!["Central", "EB", "FB", "PB"]),
                ("Motion Vision (OL)", vec!["Optic", "ME_", "LO_"]),
                ("Olfactory Valence (AL/LH)", vec!["Olfactory", "AL_", "LH_"]),
                ("Motor Descent (GNG)", vec!["Gnathal", "GNG", "LAL"]),
            ]
        } else if self.id == "celegans_wormatlas" {
            vec![
                ("Amphid Chemotaxis", vec!["Sensory", "ASH", "ASE", "AWC"]),
                ("Reversal Escape (AVA)", vec!["AVA", "AVD", "AVE"]),
                ("Forward Run (AVB)", vec!["AVB", "PVC", "DB0", "VB0"]),
                ("Body-Wall Motor", vec!["Motor", "DA0", "VA0", "DD0"]),
                ("Pharyngeal Feeding", vec!["Pharyn", "RIP", "M4_"]),
            ]
        } else {
            vec![
                ("Internal Narrative (DMN)", vec!["Default", "DMN", "ACA", "RSP"]),
                ("Executive Working Memory", vec!["Frontoparietal", "Dorsal", "FPN", "DAN", "Prefrontal"]),
                ("Language & Semantic", vec!["Language", "LANG", "Temp"]),
                ("Sensorimotor & Action", vec!["Sensorimotor", "Cerebell", "SMN", "CB.", "MOp", "SSp", "Put", "Caud"]),
                ("Visual & Salience", vec!["Visual", "Salience", "VIS", "SAL", "AIns", "Amyg"]),
            ]
        };

        let has_lesion = !self.ablated_regions.is_empty() || !self.severed_tracts.is_empty();
        let les_metrics = self.compute_summary_metrics(min_weight);
        let mut intact_clone = self.clone();
        intact_clone.heal_all();
        let intact_metrics = intact_clone.compute_summary_metrics(min_weight);

        domains
            .into_iter()
            .map(|(label, toks)| {
                let members: Vec<&BrainRegion> = self
                    .regions
                    .iter()
                    .filter(|r| toks.iter().any(|tk| r.network.contains(tk) || r.id.contains(tk)))
                    .collect();
                if members.is_empty() {
                    return DomainDecodeEntry {
                        domain_label: label.to_string(),
                        score: 0.15,
                        is_lesion_disruption: has_lesion,
                    };
                }
                if has_lesion {
                    let mut s0 = 0.0;
                    let mut s1 = 0.0;
                    for mr in &members {
                        if let Some(m0) = intact_metrics.node_metrics.iter().find(|n| n.region_id == mr.id) {
                            s0 += m0.communicability_centrality;
                        }
                        if let Some(m1) = les_metrics.node_metrics.iter().find(|n| n.region_id == mr.id) {
                            s1 += m1.communicability_centrality;
                        }
                    }
                    let drop_frac = if s0 > 1e-6 {
                        (1.0 - s1 / s0).clamp(0.0, 1.0)
                    } else {
                        0.0
                    };
                    DomainDecodeEntry {
                        domain_label: label.to_string(),
                        score: drop_frac,
                        is_lesion_disruption: true,
                    }
                } else {
                    let sel = selected_region_id.unwrap_or_else(|| self.regions[0].id.as_str());
                    let is_member = members.iter().any(|r| r.id == sel);
                    let score = if is_member { 0.88 } else { 0.28 };
                    DomainDecodeEntry {
                        domain_label: label.to_string(),
                        score,
                        is_lesion_disruption: false,
                    }
                }
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use crate::presets::load_preset_by_id;

    #[test]
    fn computes_synergy_economy_transcriptomics_and_null_models() {
        let mut human = load_preset_by_id("human_mni152").expect("human preset");
        let rep = human.compute_systems_neuro_report(0.0, 42);
        let sum = human.compute_summary_metrics(0.0);
        assert_eq!(rep.node_profiles.len(), human.regions.len());
        assert!(rep.total_wiring_cost_mm > 100.0);
        assert!(rep.null_benchmark.small_world_sigma > 0.5);
        assert!(rep.null_benchmark.obs_characteristic_path_length > 0.0);
        assert!(rep.null_benchmark.null_characteristic_path_length > 0.0);
        assert!(rep.null_benchmark.null_clustering_coeff > 0.0);

        let expected_sigma = (sum.mean_clustering / rep.null_benchmark.null_clustering_coeff)
            / (rep.null_benchmark.obs_characteristic_path_length
                / rep.null_benchmark.null_characteristic_path_length);
        assert!(
            (rep.null_benchmark.small_world_sigma - expected_sigma).abs() < 1e-9,
            "small_world_sigma must equal (C/C0)/(L/L0): {} vs {}",
            rep.null_benchmark.small_world_sigma,
            expected_sigma
        );
        assert!(rep.transcriptomic_g1_pearson_r >= 0.0 && rep.transcriptomic_g1_pearson_r <= 1.0);

        human.toggle_region_ablation("LANG.Broca_L");
        let decoded = human.decode_functional_domains(Some("LANG.Broca_L"), 0.0);
        assert_eq!(decoded.len(), 5);
        assert!(decoded.iter().all(|d| d.is_lesion_disruption));
    }
}
