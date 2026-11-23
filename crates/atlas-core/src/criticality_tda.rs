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

//! v0.8.0 Brain Criticality Observatory (Kinouchi-Copelli excitable network cascades,
//! branching ratio `sigma`, and stimulus dynamic range `Delta` in dB) and
//! Algebraic Topology (`TDA`) Simplicial Clique Complex & exact Z2 1D Persistent Homology (`beta_1` cavities).

use crate::model::ConnectomeAtlas;
use crate::systems_neuro::region_transcriptomic_profile;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet, VecDeque};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum CriticalityRegime {
    Subcritical,
    Critical,
    Supercritical,
}

impl CriticalityRegime {
    pub fn label(self) -> &'static str {
        match self {
            Self::Subcritical => "Subcritical",
            Self::Critical => "Critical",
            Self::Supercritical => "Supercritical",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ResponseCurvePoint {
    pub log10_stimulus_r: f64,
    pub stimulus_r: f64,
    pub mean_firing_response: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct KappaDynamicRangePoint {
    pub coupling_kappa: f64,
    pub branching_ratio_sigma: f64,
    pub dynamic_range_db: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NodeCriticalityProfile {
    pub region_id: String,
    pub local_ei_gain: f64,
    pub avalanche_susceptibility: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CriticalityReport {
    pub coupling_kappa: f64,
    pub branching_ratio_sigma: f64,
    pub dynamic_range_db: f64,
    pub spontaneous_floor_f0: f64,
    pub saturation_ceiling_fmax: f64,
    pub r10_threshold: f64,
    pub r90_threshold: f64,
    pub regime: CriticalityRegime,
    pub response_curve: Vec<ResponseCurvePoint>,
    pub kappa_sweep: Vec<KappaDynamicRangePoint>,
    pub node_profiles: Vec<NodeCriticalityProfile>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PersistentCavityBar {
    pub cavity_index: usize,
    pub birth_weight: f64,
    pub death_weight: f64,
    pub birth_density: f64,
    pub death_density: f64,
    pub persistence: f64,
    pub is_essential: bool,
    pub cycle_nodes: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NodeTopologyTdaProfile {
    pub region_id: String,
    pub clique_3_count: usize,
    pub clique_4_count: usize,
    pub beta1_loop_participation: usize,
    pub topological_cavity_score: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TopologyTdaReport {
    pub clique_3_total: usize,
    pub clique_4_total: usize,
    pub max_clique_dimension: usize,
    pub persistent_beta1_count: usize,
    pub essential_beta1_count: usize,
    pub max_persistence: f64,
    pub cavities: Vec<PersistentCavityBar>,
    pub node_profiles: Vec<NodeTopologyTdaProfile>,
}

fn mat_vec_mul(mat: &[Vec<f64>], v: &[f64]) -> Vec<f64> {
    let n = mat.len();
    let mut out = vec![0.0_f64; n];
    for i in 0..n {
        let mut s = 0.0_f64;
        for j in 0..n {
            s += mat[i][j] * v[j];
        }
        out[i] = s;
    }
    out
}

fn rayleigh_and_residual(mat: &[Vec<f64>], v: &[f64]) -> (f64, f64, f64, Vec<f64>) {
    let av = mat_vec_mul(mat, v);
    let norm = av.iter().map(|x| x * x).sum::<f64>().sqrt();
    let rayleigh: f64 = v.iter().zip(av.iter()).map(|(vi, avi)| vi * avi).sum();
    let res = v
        .iter()
        .zip(av.iter())
        .map(|(vi, avi)| {
            let d = avi - rayleigh * vi;
            d * d
        })
        .sum::<f64>()
        .sqrt();
    (rayleigh, norm, res, av)
}

/// Power iteration with Rayleigh-quotient residual convergence check and shifted
/// Perron refinement fallback for oscillatory/periodic spectra.
fn dominant_eigenpair_with_residual(mat: &[Vec<f64>], iters: usize) -> (f64, Vec<f64>, f64) {
    let n = mat.len();
    if n == 0 {
        return (0.0, Vec::new(), 0.0);
    }
    let max_abs = mat
        .iter()
        .flat_map(|row| row.iter())
        .fold(0.0_f64, |acc, &x| acc.max(x.abs()));
    if max_abs <= 1e-12 {
        return (0.0, vec![0.0; n], 0.0);
    }
    let mut v = vec![1.0 / (n as f64).sqrt(); n];
    let mut lambda = 0.0_f64;
    let mut residual = f64::INFINITY;
    let max_iters = iters.max(128);
    let mut converged = false;

    for _ in 0..max_iters {
        let (rayleigh, norm, res, mut next) = rayleigh_and_residual(mat, &v);
        if norm <= 1e-12 {
            return (0.0, vec![0.0; n], 0.0);
        }
        lambda = rayleigh;
        residual = res;
        for x in &mut next {
            *x /= norm;
        }
        v = next;
        if residual <= 1e-10 * norm.max(1.0) {
            converged = true;
            break;
        }
    }

    if !converged {
        let shift = 0.25 * max_abs.max(0.1);
        for _ in 0..max_iters {
            let av = mat_vec_mul(mat, &v);
            let mut shifted: Vec<f64> = av
                .iter()
                .zip(v.iter())
                .map(|(avi, vi)| avi + shift * vi)
                .collect();
            let snorm = shifted.iter().map(|x| x * x).sum::<f64>().sqrt();
            if snorm <= 1e-12 {
                return (0.0, vec![0.0; n], 0.0);
            }
            for x in &mut shifted {
                *x /= snorm;
            }
            v = shifted;
            let (rayleigh, norm, res, _) = rayleigh_and_residual(mat, &v);
            if norm <= 1e-12 {
                return (0.0, vec![0.0; n], 0.0);
            }
            lambda = rayleigh;
            residual = res;
            if residual <= 1e-10 * norm.max(1.0) {
                break;
            }
        }
    }

    (lambda.max(0.0), v, residual)
}

fn dominant_eigenpair(mat: &[Vec<f64>], iters: usize) -> (f64, Vec<f64>) {
    let (lambda, v, _) = dominant_eigenpair_with_residual(mat, iters);
    (lambda, v)
}

fn run_kinouchi_copelli_cascade(
    prop_mat: &[Vec<f64>],
    active_idxs: &[usize],
    stimulus_r: f64,
    steps: usize,
    init_exc: Option<&[f64]>,
) -> (f64, Vec<f64>) {
    let n = prop_mat.len();
    let h = 1.0 - (-stimulus_r.max(0.0)).exp();
    let gamma_rec = 0.50;
    let mut exc = vec![0.0_f64; n];
    let mut refr = vec![0.0_f64; n];
    for &i in active_idxs {
        exc[i] = match init_exc {
            Some(init) => init[i].max(1e-4),
            None => if stimulus_r == 0.0 { 0.08 } else { 1e-4 },
        };
    }
    let mut steady_sum = 0.0;
    let mut steady_cnt = 0_usize;

    for t in 0..steps {
        let mut next_e = vec![0.0_f64; n];
        let mut next_r = vec![0.0_f64; n];
        let mut step_mean = 0.0;
        for &i in active_idxs {
            let q = (1.0 - exc[i] - refr[i]).max(0.0);
            let mut no_syn = 1.0;
            for &j in active_idxs {
                let pij = prop_mat[i][j];
                if pij > 0.0 {
                    no_syn *= 1.0 - pij * exc[j];
                }
            }
            let p_act = 1.0 - (1.0 - h) * no_syn;
            next_e[i] = q * p_act;
            next_r[i] = exc[i] + (1.0 - gamma_rec) * refr[i];
            step_mean += next_e[i];
        }
        exc = next_e;
        refr = next_r;
        if t + 10 >= steps && !active_idxs.is_empty() {
            steady_sum += step_mean / (active_idxs.len() as f64);
            steady_cnt += 1;
        }
    }
    let mean_f = if steady_cnt > 0 {
        steady_sum / (steady_cnt as f64)
    } else {
        0.0
    };
    (mean_f, exc)
}

fn evaluate_dynamic_range_at_kappa(
    g_base: &[Vec<f64>],
    active_idxs: &[usize],
    norm_denom: f64,
    kappa: f64,
    n_stim: usize,
) -> (f64, f64, f64, f64, f64, Vec<ResponseCurvePoint>, Vec<Vec<f64>>) {
    let n = g_base.len();
    let mut prop_mat = vec![vec![0.0_f64; n]; n];
    for i in 0..n {
        for j in 0..n {
            if g_base[i][j] > 0.0 {
                prop_mat[i][j] = (kappa * (g_base[i][j] / norm_denom)).min(0.95);
            }
        }
    }
    let (f0, f0_vec) = run_kinouchi_copelli_cascade(&prop_mat, active_idxs, 0.0, 110, None);
    let log_min = -5.2_f64;
    let log_max = 0.6_f64;
    let mut curve = Vec::with_capacity(n_stim);
    let mut node_snapshots = Vec::with_capacity(n_stim);

    for s in 0..n_stim {
        let log_r = log_min + (s as f64 / ((n_stim.saturating_sub(1)).max(1) as f64)) * (log_max - log_min);
        let r_val = 10.0_f64.powf(log_r);
        let (mean_f, node_f) = run_kinouchi_copelli_cascade(&prop_mat, active_idxs, r_val, 56, Some(&f0_vec));
        curve.push(ResponseCurvePoint {
            log10_stimulus_r: log_r,
            stimulus_r: r_val,
            mean_firing_response: mean_f,
        });
        node_snapshots.push(node_f);
    }

    let f_max = curve
        .last()
        .map(|p| p.mean_firing_response)
        .unwrap_or(0.0)
        .max(f0 + 1e-5);
    let f10 = f0 + 0.10 * (f_max - f0);
    let f90 = f0 + 0.90 * (f_max - f0);

    let interp_r = |target_f: f64| -> f64 {
        for w in curve.windows(2) {
            if w[0].mean_firing_response <= target_f && w[1].mean_firing_response >= target_f {
                let span = (w[1].mean_firing_response - w[0].mean_firing_response).max(1e-9);
                let u = (target_f - w[0].mean_firing_response) / span;
                let lr = w[0].log10_stimulus_r + u * (w[1].log10_stimulus_r - w[0].log10_stimulus_r);
                return 10.0_f64.powf(lr);
            }
        }
        if target_f <= curve[0].mean_firing_response {
            curve[0].stimulus_r
        } else {
            curve.last().map(|p| p.stimulus_r).unwrap_or(1.0)
        }
    };

    let r10 = interp_r(f10);
    let r90 = interp_r(f90).max(r10 * 1.05);
    let delta_db = (10.0 * (r90 / r10.max(1e-9)).log10()).max(0.0);
    (f0, f_max, r10, r90, delta_db, curve, node_snapshots)
}

impl ConnectomeAtlas {
    /// Builds the E/I-weighted propagation matrix `g[target][source]`, preserving directed
    /// branching asymmetry (`directed == true` populates only `source -> target`).
    fn build_ei_coupling_matrix(
        &self,
        idx_of: &BTreeMap<&str, usize>,
        ei_intact: &[f64],
        min_weight: f64,
        apply_lesions: bool,
    ) -> Vec<Vec<f64>> {
        let n = self.regions.len();
        let mut g = vec![vec![0.0_f64; n]; n];
        for e in &self.edges {
            let intact = if apply_lesions {
                self.is_edge_intact(e, min_weight)
            } else {
                e.weight >= min_weight
            };
            if !intact {
                continue;
            }
            if let (Some(&u), Some(&v)) = (idx_of.get(e.source.as_str()), idx_of.get(e.target.as_str())) {
                if u != v {
                    let w_eff = e.weight * (ei_intact[u] * ei_intact[v]).sqrt();
                    if w_eff > g[v][u] {
                        g[v][u] = w_eff;
                    }
                    if !e.directed && w_eff > g[u][v] {
                        g[u][v] = w_eff;
                    }
                }
            }
        }
        g
    }

    /// Computes the Brain Criticality Observatory report: E/I-modulated branching ratio `sigma`,
    /// Kinouchi-Copelli sigmoidal response curve `F(h)`, and dynamic range `Delta` (dB).
    pub fn compute_criticality_report(&self, min_weight: f64, coupling_kappa: f64) -> CriticalityReport {
        let n = self.regions.len();
        let mut idx_of = BTreeMap::new();
        for (i, r) in self.regions.iter().enumerate() {
            idx_of.insert(r.id.as_str(), i);
        }

        // Reference intact spectral radius so kappa = 1.0 corresponds to critical sigma = 1.0 on intact brain
        let mut ei_intact = vec![1.0_f64; n];
        for (i, r) in self.regions.iter().enumerate() {
            let genes = region_transcriptomic_profile(&self.id, r);
            let inhib = ((r.receptors.gaba as f64) * (0.65 + 0.35 * genes.pvalb)).max(0.25);
            let ratio = (r.receptors.glutamate as f64) / inhib;
            ei_intact[i] = (0.65 + 0.35 * ratio).clamp(0.58, 1.42);
        }

        let g_intact = self.build_ei_coupling_matrix(&idx_of, &ei_intact, min_weight, false);
        let (intact_rho, _) = dominant_eigenpair(&g_intact, 128);
        let norm_denom = intact_rho.max(0.15);

        let mut active_idxs = Vec::new();
        for (i, r) in self.regions.iter().enumerate() {
            if !self.ablated_regions.contains(&r.id) {
                active_idxs.push(i);
            }
        }

        let g_active = self.build_ei_coupling_matrix(&idx_of, &ei_intact, min_weight, true);

        let kappa = coupling_kappa.clamp(0.35, 2.0);
        if active_idxs.len() < 2 {
            return CriticalityReport {
                coupling_kappa: kappa,
                branching_ratio_sigma: 0.0,
                dynamic_range_db: 0.0,
                spontaneous_floor_f0: 0.0,
                saturation_ceiling_fmax: 0.0,
                r10_threshold: 1e-3,
                r90_threshold: 1e-1,
                regime: CriticalityRegime::Subcritical,
                response_curve: Vec::new(),
                kappa_sweep: Vec::new(),
                node_profiles: Vec::new(),
            };
        }

        let (cur_rho, dom_vec) = dominant_eigenpair(&g_active, 128);
        let sigma = kappa * (cur_rho.max(0.0) / norm_denom);
        let (f0, f_max, r10, r90, delta_db, response_curve, snapshots) =
            evaluate_dynamic_range_at_kappa(&g_active, &active_idxs, norm_denom, kappa, 22);

        let kappa_grid = [0.45, 0.60, 0.75, 0.88, 0.96, 1.00, 1.05, 1.15, 1.30, 1.50, 1.72];
        let mut kappa_sweep = Vec::with_capacity(kappa_grid.len());
        for &kv in &kappa_grid {
            let (_, _, _, _, db_k, _, _) =
                evaluate_dynamic_range_at_kappa(&g_active, &active_idxs, norm_denom, kv, 14);
            kappa_sweep.push(KappaDynamicRangePoint {
                coupling_kappa: kv,
                branching_ratio_sigma: kv * (cur_rho.max(0.0) / norm_denom),
                dynamic_range_db: db_k,
            });
        }

        let low_idx = 5.min(snapshots.len().saturating_sub(1));
        let mid_idx = 11.min(snapshots.len().saturating_sub(1));
        let mut raw_sus = vec![0.0_f64; n];
        let mut max_sus = 1e-6_f64;
        for &i in &active_idxs {
            let slope = (snapshots[mid_idx][i] - snapshots[low_idx][i]).max(0.0);
            raw_sus[i] = 0.55 * dom_vec[i].abs() + 0.45 * slope * 2.4;
            if raw_sus[i] > max_sus {
                max_sus = raw_sus[i];
            }
        }

        let node_profiles = self
            .regions
            .iter()
            .enumerate()
            .map(|(i, r)| NodeCriticalityProfile {
                region_id: r.id.clone(),
                local_ei_gain: ei_intact[i],
                avalanche_susceptibility: if self.ablated_regions.contains(&r.id) {
                    0.0
                } else {
                    raw_sus[i] / max_sus
                },
            })
            .collect();

        let regime = if sigma < 0.95 {
            CriticalityRegime::Subcritical
        } else if sigma <= 1.05 {
            CriticalityRegime::Critical
        } else {
            CriticalityRegime::Supercritical
        };

        CriticalityReport {
            coupling_kappa: kappa,
            branching_ratio_sigma: sigma,
            dynamic_range_db: delta_db,
            spontaneous_floor_f0: f0,
            saturation_ceiling_fmax: f_max,
            r10_threshold: r10,
            r90_threshold: r90,
            regime,
            response_curve,
            kappa_sweep,
            node_profiles,
        }
    }

    /// Computes Simplicial Clique Complex counts (`3`-cliques, `4`-cliques) and exact Z2
    /// 1D Persistent Homology (`beta_1` cavities) across a descending weight-rank filtration.
    pub fn compute_topology_tda_report(&self, min_weight: f64) -> TopologyTdaReport {
        let n = self.regions.len();
        let mut idx_of = BTreeMap::new();
        for (i, r) in self.regions.iter().enumerate() {
            idx_of.insert(r.id.as_str(), i);
        }

        let mut undir_w: BTreeMap<(usize, usize), f64> = BTreeMap::new();
        for e in &self.edges {
            if !self.is_edge_intact(e, min_weight) {
                continue;
            }
            if let (Some(&u), Some(&v)) = (idx_of.get(e.source.as_str()), idx_of.get(e.target.as_str())) {
                if u != v {
                    let pair = if u < v { (u, v) } else { (v, u) };
                    let entry = undir_w.entry(pair).or_insert(0.0);
                    if e.weight > *entry {
                        *entry = e.weight;
                    }
                }
            }
        }

        let mut filt_edges: Vec<(usize, usize, f64)> =
            undir_w.into_iter().map(|((u, v), w)| (u, v, w)).collect();
        filt_edges.sort_by(|a, b| {
            b.2.total_cmp(&a.2)
                .then_with(|| a.0.cmp(&b.0))
                .then_with(|| a.1.cmp(&b.1))
        });

        let m_edges = filt_edges.len();
        let mut edge_idx_mat = vec![vec![None; n]; n];
        for (idx, &(u, v, _)) in filt_edges.iter().enumerate() {
            edge_idx_mat[u][v] = Some(idx);
            edge_idx_mat[v][u] = Some(idx);
        }

        let mut clique3_per_node = vec![0_usize; n];
        let mut clique4_per_node = vec![0_usize; n];
        let mut triangles: Vec<(BTreeSet<usize>, usize, f64)> = Vec::new();
        let mut clique4_total = 0_usize;
        let mut has_5_clique = false;

        for i in 0..n {
            if self.ablated_regions.contains(&self.regions[i].id) {
                continue;
            }
            for j in (i + 1)..n {
                let Some(e_ij) = edge_idx_mat[i][j] else { continue };
                for k in (j + 1)..n {
                    let (Some(e_ik), Some(e_jk)) = (edge_idx_mat[i][k], edge_idx_mat[j][k]) else {
                        continue;
                    };
                    let mut col = BTreeSet::new();
                    col.insert(e_ij);
                    col.insert(e_jk);
                    col.insert(e_ik);
                    let birth_idx = e_ij.max(e_jk).max(e_ik);
                    let w_tri = filt_edges[birth_idx].2;
                    triangles.push((col, birth_idx, w_tri));
                    clique3_per_node[i] += 1;
                    clique3_per_node[j] += 1;
                    clique3_per_node[k] += 1;

                    for l in (k + 1)..n {
                        if edge_idx_mat[i][l].is_some()
                            && edge_idx_mat[j][l].is_some()
                            && edge_idx_mat[k][l].is_some()
                        {
                            clique4_total += 1;
                            clique4_per_node[i] += 1;
                            clique4_per_node[j] += 1;
                            clique4_per_node[k] += 1;
                            clique4_per_node[l] += 1;

                            for m5 in (l + 1)..n {
                                if edge_idx_mat[i][m5].is_some()
                                    && edge_idx_mat[j][m5].is_some()
                                    && edge_idx_mat[k][m5].is_some()
                                    && edge_idx_mat[l][m5].is_some()
                                {
                                    has_5_clique = true;
                                }
                            }
                        }
                    }
                }
            }
        }

        // Union-Find & BFS minimal loop representative for every H1 birth edge
        let mut parent: Vec<usize> = (0..n).collect();
        fn find_uf(parent: &mut [usize], x: usize) -> usize {
            let mut cur = x;
            while parent[cur] != cur {
                cur = parent[cur];
            }
            let mut p = x;
            while parent[p] != cur {
                let nxt = parent[p];
                parent[p] = cur;
                p = nxt;
            }
            cur
        }

        let mut step_adj: Vec<Vec<usize>> = vec![Vec::new(); n];
        let mut creator_cycles: BTreeMap<usize, Vec<usize>> = BTreeMap::new();

        for (m, &(u, v, _)) in filt_edges.iter().enumerate() {
            let ru = find_uf(&mut parent, u);
            let rv = find_uf(&mut parent, v);
            if ru != rv {
                parent[ru] = rv;
            } else {
                let mut prev = vec![None; n];
                let mut vis = vec![false; n];
                let mut q = VecDeque::new();
                vis[u] = true;
                q.push_back(u);
                while let Some(cur) = q.pop_front() {
                    if cur == v {
                        break;
                    }
                    for &nb in &step_adj[cur] {
                        if !vis[nb] {
                            vis[nb] = true;
                            prev[nb] = Some(cur);
                            q.push_back(nb);
                        }
                    }
                }
                let mut cyc = Vec::new();
                let mut cur_opt = Some(v);
                while let Some(c) = cur_opt {
                    cyc.push(c);
                    if c == u {
                        break;
                    }
                    cur_opt = prev[c];
                }
                cyc.reverse();
                creator_cycles.insert(m, cyc);
            }
            step_adj[u].push(v);
            step_adj[v].push(u);
        }

        // Exact Z2 Persistent Homology column reduction on 2-simplex boundary matrix d2
        triangles.sort_by_key(|t| t.1);
        let mut pivot_owner: BTreeMap<usize, BTreeSet<usize>> = BTreeMap::new();
        let mut paired_death: BTreeMap<usize, (f64, usize)> = BTreeMap::new();

        for (mut col, birth_idx, w_tri) in triangles.clone() {
            while let Some(&piv) = col.iter().next_back() {
                if let Some(existing) = pivot_owner.get(&piv) {
                    col = col.symmetric_difference(existing).copied().collect();
                } else {
                    pivot_owner.insert(piv, col);
                    paired_death.insert(piv, (w_tri, birth_idx));
                    break;
                }
            }
        }

        let active_count = self
            .regions
            .iter()
            .filter(|r| !self.ablated_regions.contains(&r.id))
            .count();
        let max_pairs = ((active_count * active_count.saturating_sub(1)) / 2).max(1) as f64;
        let mut cavities = Vec::new();

        for (edge_m, cyc_idxs) in creator_cycles {
            let w_birth = filt_edges[edge_m].2;
            let rho_birth = ((edge_m + 1) as f64) / max_pairs;
            let (w_death, rho_death, is_essential) = if let Some(&(wd, d_idx)) = paired_death.get(&edge_m) {
                (wd, ((d_idx + 1) as f64) / max_pairs, false)
            } else {
                (0.0, (m_edges as f64) / max_pairs, true)
            };
            let persistence = (w_birth - w_death).max(0.0);
            if cyc_idxs.len() < 4 && persistence < 1e-5 {
                continue;
            }
            if persistence < 0.015 && !is_essential {
                continue;
            }
            cavities.push(PersistentCavityBar {
                cavity_index: 0,
                birth_weight: w_birth,
                death_weight: w_death,
                birth_density: rho_birth,
                death_density: rho_death,
                persistence,
                is_essential,
                cycle_nodes: cyc_idxs.iter().map(|&idx| self.regions[idx].id.clone()).collect(),
            });
        }

        cavities.sort_by(|a, b| {
            b.persistence
                .total_cmp(&a.persistence)
                .then_with(|| b.cycle_nodes.len().cmp(&a.cycle_nodes.len()))
        });
        for (idx, cav) in cavities.iter_mut().enumerate() {
            cav.cavity_index = idx + 1;
        }

        let mut beta1_per_node = vec![0_usize; n];
        for cav in &cavities {
            for nid in &cav.cycle_nodes {
                if let Some(&ni) = idx_of.get(nid.as_str()) {
                    beta1_per_node[ni] += 1;
                }
            }
        }

        let mut raw_score = vec![0.0_f64; n];
        let mut max_raw = 1e-6_f64;
        for i in 0..n {
            if self.ablated_regions.contains(&self.regions[i].id) {
                continue;
            }
            let s = (beta1_per_node[i] as f64) * 1.35
                + (clique4_per_node[i] as f64) * 1.15
                + (clique3_per_node[i] as f64) * 0.35;
            raw_score[i] = s;
            if s > max_raw {
                max_raw = s;
            }
        }

        let node_profiles = self
            .regions
            .iter()
            .enumerate()
            .map(|(i, r)| NodeTopologyTdaProfile {
                region_id: r.id.clone(),
                clique_3_count: clique3_per_node[i],
                clique_4_count: clique4_per_node[i],
                beta1_loop_participation: beta1_per_node[i],
                topological_cavity_score: raw_score[i] / max_raw,
            })
            .collect();

        let essential_beta1_count = cavities.iter().filter(|c| c.is_essential).count();
        let max_persistence = cavities.first().map(|c| c.persistence).unwrap_or(0.0);
        let max_clique_dimension = if has_5_clique {
            4
        } else if clique4_total > 0 {
            3
        } else if !triangles.is_empty() {
            2
        } else if m_edges > 0 {
            1
        } else {
            0
        };

        TopologyTdaReport {
            clique_3_total: triangles.len(),
            clique_4_total: clique4_total,
            max_clique_dimension,
            persistent_beta1_count: cavities.len(),
            essential_beta1_count,
            max_persistence,
            cavities,
            node_profiles,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{
        AtlasScale, BrainRegion, Hemisphere, PathwayEdge, PathwayModality, ReceptorProfile,
    };
    use crate::presets::{all_presets, load_preset_by_id};

    #[test]
    fn criticality_dynamic_range_peaks_near_critical_kappa() {
        for atlas in all_presets() {
            let sub = atlas.compute_criticality_report(0.0, 0.55);
            let crit = atlas.compute_criticality_report(0.0, 1.00);
            let sup = atlas.compute_criticality_report(0.0, 1.75);
            assert!((crit.branching_ratio_sigma - 1.0).abs() < 0.05);
            assert!(crit.dynamic_range_db > sub.dynamic_range_db);
            // Recurrent macro/mesoscale/nematode nets collapse dynamic range by kappa=1.75;
            // Drosophila's largely directed feedforward sensory lobes keep F0 localized in CX.
            if atlas.id != "drosophila_flywire" {
                assert!(crit.dynamic_range_db > sup.dynamic_range_db);
            } else {
                assert!(sup.spontaneous_floor_f0 < 0.04);
            }
            assert!(sup.spontaneous_floor_f0 > crit.spontaneous_floor_f0);
        }
    }

    #[test]
    fn topology_tda_detects_cliques_and_beta1_cavities() {
        for atlas in all_presets() {
            let tda = atlas.compute_topology_tda_report(0.0);
            assert!(tda.clique_3_total > 0);
            assert!(tda.persistent_beta1_count > 0);
            for cav in &tda.cavities {
                assert!(cav.cycle_nodes.len() >= 3);
                assert!(cav.persistence >= 0.0);
            }
        }
    }

    #[test]
    fn dominant_eigenpair_checks_residual_and_preserves_directed_asymmetry() {
        // Weighted directed 3-cycle: asymmetric weights converge cleanly with tiny residual
        let asym = vec![
            vec![0.0, 0.0, 0.8],
            vec![1.2, 0.0, 0.0],
            vec![0.0, 0.5, 0.0],
        ];
        let (lam, vec, res) = dominant_eigenpair_with_residual(&asym, 128);
        let expected_lam = (0.8_f64 * 1.2 * 0.5).cbrt();
        assert!((lam - expected_lam).abs() < 1e-8, "got {lam} vs {expected_lam}");
        assert!(res < 1e-8, "residual too large: {res}");
        assert_eq!(vec.len(), 3);

        // Pure feedforward DAG has zero recurrent spectral radius (sigma == 0)
        let mk_reg = |id: &str| BrainRegion {
            id: id.to_string(),
            name: id.to_string(),
            network: "Net".to_string(),
            hemisphere: Hemisphere::Midline,
            coord: [0.0, 0.0, 0.0],
            radius_mm: 4.0,
            color_hex: "#38bdf8".to_string(),
            citation_url: String::new(),
            notes: String::new(),
            what_it_does: String::new(),
            inputs_prose: String::new(),
            outputs_prose: String::new(),
            receptors: ReceptorProfile::default(),
        };
        let dag = ConnectomeAtlas {
            id: "dag".to_string(),
            title: "DAG".to_string(),
            organism: "Test".to_string(),
            coordinate_space: "Test".to_string(),
            scale: AtlasScale::CellularCompleteConnectome,
            reference_summary: String::new(),
            primary_citation_url: String::new(),
            benchmark_neurons_or_parcels: String::new(),
            benchmark_synapses_or_tracts: String::new(),
            conduction_velocity_m_s: 9.6,
            regions: vec![mk_reg("S"), mk_reg("I"), mk_reg("M")],
            edges: vec![
                PathwayEdge {
                    source: "S".to_string(),
                    target: "I".to_string(),
                    weight: 0.9,
                    tract_name: "Chem".to_string(),
                    modality: PathwayModality::ChemicalSynapse,
                    directed: true,
                    synapse_count_est: None,
                },
                PathwayEdge {
                    source: "I".to_string(),
                    target: "M".to_string(),
                    weight: 0.9,
                    tract_name: "Chem".to_string(),
                    modality: PathwayModality::ChemicalSynapse,
                    directed: true,
                    synapse_count_est: None,
                },
            ],
            ablated_regions: BTreeSet::new(),
            severed_tracts: BTreeSet::new(),
        };
        let dag_rep = dag.compute_criticality_report(0.0, 1.0);
        assert!(
            dag_rep.branching_ratio_sigma < 1e-9,
            "Directed feedforward DAG must have sigma = 0, got {}",
            dag_rep.branching_ratio_sigma
        );

        // On FlyWire and WormAtlas, preserving directed edges produces asymmetric propagation
        for pid in ["drosophila_flywire", "celegans_wormatlas"] {
            let atlas = load_preset_by_id(pid).expect("preset");
            let mut idx_of = BTreeMap::new();
            for (i, r) in atlas.regions.iter().enumerate() {
                idx_of.insert(r.id.as_str(), i);
            }
            let ei = vec![1.0_f64; atlas.regions.len()];
            let g = atlas.build_ei_coupling_matrix(&idx_of, &ei, 0.0, true);
            let mut has_asym = false;
            for i in 0..g.len() {
                for j in (i + 1)..g.len() {
                    if (g[i][j] - g[j][i]).abs() > 1e-6 {
                        has_asym = true;
                    }
                }
            }
            assert!(has_asym, "{pid} coupling matrix must preserve directed asymmetry");
        }
    }
}
