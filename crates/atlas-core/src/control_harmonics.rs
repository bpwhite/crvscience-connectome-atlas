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

use crate::model::ConnectomeAtlas;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NodeControlHarmonicMetrics {
    pub region_id: String,
    pub average_controllability: f64,
    pub modal_controllability: f64,
    pub principal_gradient_g1: f64,
    pub secondary_gradient_g2: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ControlHarmonicsReport {
    pub receptor_channel: String,
    pub agonist_gain: f64,
    pub placebo_mean_transition_energy: f64,
    pub modulated_mean_transition_energy: f64,
    pub energy_landscape_delta_pct: f64,
    pub state_names: Vec<String>,
    pub transition_energy_matrix: Vec<Vec<f64>>,
    pub laplacian_wavenumbers: Vec<f64>,
    pub harmonic_modes: Vec<Vec<f64>>,
    pub node_profiles: Vec<NodeControlHarmonicMetrics>,
}

/// Cyclic Jacobi symmetric eigensolver in 100% safe Rust (`A = V * diag(eigenvalues) * V^T`).
/// Eigenvalues and corresponding eigenvector columns are sorted in ascending order.
pub fn jacobi_eigen_symmetric(matrix: &[Vec<f64>], max_sweeps: usize) -> (Vec<f64>, Vec<Vec<f64>>) {
    let n = matrix.len();
    if n == 0 {
        return (Vec::new(), Vec::new());
    }
    let mut a = matrix.to_vec();
    let mut v = vec![vec![0.0_f64; n]; n];
    for (i, row) in v.iter_mut().enumerate() {
        row[i] = 1.0;
    }

    for _sweep in 0..max_sweeps.max(12) {
        let mut max_off = 0.0_f64;
        for p in 0..n {
            for q in (p + 1)..n {
                let apq = a[p][q].abs();
                if apq > max_off {
                    max_off = apq;
                }
                if apq < 1e-12 {
                    continue;
                }
                let app = a[p][p];
                let aqq = a[q][q];
                let tau = (aqq - app) / (2.0 * a[p][q]);
                let t = if tau >= 0.0 {
                    1.0 / (tau + (1.0 + tau * tau).sqrt())
                } else {
                    -1.0 / (-tau + (1.0 + tau * tau).sqrt())
                };
                let c = 1.0 / (1.0 + t * t).sqrt();
                let s = t * c;

                a[p][p] = app - t * a[p][q];
                a[q][q] = aqq + t * a[p][q];
                a[p][q] = 0.0;
                a[q][p] = 0.0;

                for r in 0..n {
                    if r != p && r != q {
                        let arp = a[r][p];
                        let arq = a[r][q];
                        a[r][p] = c * arp - s * arq;
                        a[p][r] = a[r][p];
                        a[r][q] = s * arp + c * arq;
                        a[q][r] = a[r][q];
                    }
                    let vrp = v[r][p];
                    let vrq = v[r][q];
                    v[r][p] = c * vrp - s * vrq;
                    v[r][q] = s * vrp + c * vrq;
                }
            }
        }
        if max_off < 1e-11 {
            break;
        }
    }

    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&i, &j| a[i][i].total_cmp(&a[j][j]));

    let mut evals = vec![0.0_f64; n];
    let mut evecs = vec![vec![0.0_f64; n]; n];
    for (col, &src_col) in order.iter().enumerate() {
        evals[col] = a[src_col][src_col];
        let mut best_idx = 0;
        let mut best_abs = -1.0_f64;
        for (r, row) in v.iter().enumerate() {
            let av = row[src_col].abs();
            if av > best_abs {
                best_abs = av;
                best_idx = r;
            }
        }
        let sign = if v[best_idx][src_col] < 0.0 { -1.0 } else { 1.0 };
        for r in 0..n {
            evecs[r][col] = v[r][src_col] * sign;
        }
    }
    (evals, evecs)
}

/// Gaussian elimination with partial pivoting for `A x = b`.
pub fn solve_linear_system(mat_a: &[Vec<f64>], vec_b: &[f64]) -> Vec<f64> {
    let n = vec_b.len();
    let mut m = vec![vec![0.0_f64; n + 1]; n];
    for i in 0..n {
        m[i][..n].copy_from_slice(&mat_a[i][..n]);
        m[i][n] = vec_b[i];
    }
    for col in 0..n {
        let mut piv = col;
        let mut max_abs = m[col][col].abs();
        for (r, row) in m.iter().enumerate().skip(col + 1) {
            if row[col].abs() > max_abs {
                max_abs = row[col].abs();
                piv = r;
            }
        }
        if piv != col {
            m.swap(col, piv);
        }
        let diag = if m[col][col].abs() < 1e-9 { 1e-9 } else { m[col][col] };
        for r in (col + 1)..n {
            let factor = m[r][col] / diag;
            if factor == 0.0 {
                continue;
            }
            for c in col..=n {
                let sub = factor * m[col][c];
                m[r][c] -= sub;
            }
        }
    }
    let mut x = vec![0.0_f64; n];
    for i in (0..n).rev() {
        let mut sum = m[i][n];
        for j in (i + 1)..n {
            sum -= m[i][j] * x[j];
        }
        let diag = if m[i][i].abs() < 1e-9 { 1e-9 } else { m[i][i] };
        x[i] = sum / diag;
    }
    x
}

impl ConnectomeAtlas {
    /// Computes Receptor-Informed Network Control Theory (Gu et al. 2015; Singleton et al. 2022)
    /// and Normalized Graph Laplacian Connectome Harmonics (Margulies et al. 2016; Atasoy et al. 2016).
    pub fn compute_control_and_harmonics(
        &self,
        min_weight: f64,
        receptor_channel: &str,
        agonist_gain: f64,
    ) -> ControlHarmonicsReport {
        let n = self.regions.len();
        let mut idx_of = HashMap::with_capacity(n);
        for (i, r) in self.regions.iter().enumerate() {
            idx_of.insert(r.id.as_str(), i);
        }

        let mut sym_w = vec![vec![0.0_f64; n]; n];
        for e in &self.edges {
            if !self.is_edge_intact(e, min_weight) {
                continue;
            }
            if let (Some(&u), Some(&v)) = (idx_of.get(e.source.as_str()), idx_of.get(e.target.as_str())) {
                if u != v {
                    sym_w[u][v] = sym_w[u][v].max(e.weight);
                    if !e.directed {
                        sym_w[v][u] = sym_w[v][u].max(e.weight);
                    } else {
                        sym_w[v][u] = sym_w[v][u].max(e.weight * 0.5);
                    }
                }
            }
        }

        let strengths: Vec<f64> = (0..n).map(|i| sym_w[i].iter().sum()).collect();
        let mut lap_sym = vec![vec![0.0_f64; n]; n];
        for i in 0..n {
            if !self.ablated_regions.contains(&self.regions[i].id) && strengths[i] > 0.0 {
                lap_sym[i][i] = 1.0;
            }
            for j in 0..n {
                if i != j && sym_w[i][j] > 0.0 && strengths[i] > 0.0 && strengths[j] > 0.0 {
                    lap_sym[i][j] = -sym_w[i][j] / (strengths[i] * strengths[j]).sqrt();
                }
            }
        }

        let (lap_evals, lap_evecs) = jacobi_eigen_symmetric(&lap_sym, 22);
        let non_zero_cols: Vec<usize> = (0..n).filter(|&c| lap_evals[c] > 1e-4).collect();
        let g1_col = non_zero_cols.first().copied().unwrap_or(0);
        let g2_col = non_zero_cols.get(1).copied().unwrap_or(g1_col);

        let (w_evals, w_evecs) = jacobi_eigen_symmetric(&sym_w, 20);
        let lambda_max = w_evals.last().copied().unwrap_or(1.0).max(0.1);
        let scale_a = 1.0 / (1.0 + lambda_max);
        let mut a_norm = vec![vec![0.0_f64; n]; n];
        for i in 0..n {
            for j in 0..n {
                a_norm[i][j] = sym_w[i][j] * scale_a;
            }
        }

        let (state_names, placebo_mat, placebo_mean, _) =
            self.eval_transition_energy_grid(&a_norm, receptor_channel, 0.0);
        let (_, mod_mat, mod_mean, b_diag) = if agonist_gain.abs() < 1e-9 {
            let b0: Vec<f64> = (0..n)
                .map(|i| if self.ablated_regions.contains(&self.regions[i].id) { 0.0 } else { 1.0 })
                .collect();
            (state_names.clone(), placebo_mat.clone(), placebo_mean, b0)
        } else {
            self.eval_transition_energy_grid(&a_norm, receptor_channel, agonist_gain)
        };

        let delta_pct = if placebo_mean > 1e-9 {
            ((mod_mean - placebo_mean) / placebo_mean) * 100.0
        } else {
            0.0
        };

        let mut node_profiles = Vec::with_capacity(n);
        for i in 0..n {
            if self.ablated_regions.contains(&self.regions[i].id) {
                node_profiles.push(NodeControlHarmonicMetrics {
                    region_id: self.regions[i].id.clone(),
                    average_controllability: 0.0,
                    modal_controllability: 0.0,
                    principal_gradient_g1: 0.0,
                    secondary_gradient_g2: 0.0,
                });
                continue;
            }
            let bi2 = b_diag[i] * b_diag[i];
            let mut avg_c = 0.0_f64;
            let mut mod_c = 0.0_f64;
            for m in 0..n {
                let lam = w_evals[m] * scale_a;
                let vim2 = w_evecs[i][m] * w_evecs[i][m];
                avg_c += bi2 * (vim2 / (1.0 - lam * lam).max(1e-6));
                mod_c += (1.0 - lam * lam) * vim2;
            }
            node_profiles.push(NodeControlHarmonicMetrics {
                region_id: self.regions[i].id.clone(),
                average_controllability: avg_c,
                modal_controllability: mod_c,
                principal_gradient_g1: if n > 0 { lap_evecs[i][g1_col] } else { 0.0 },
                secondary_gradient_g2: if n > 0 { lap_evecs[i][g2_col] } else { 0.0 },
            });
        }

        // Transpose harmonic columns into mode-major vectors `harmonic_modes[k][node]`
        let mut harmonic_modes = Vec::with_capacity(non_zero_cols.len());
        let mut wavenumbers = Vec::with_capacity(non_zero_cols.len());
        for &c in &non_zero_cols {
            wavenumbers.push(lap_evals[c]);
            harmonic_modes.push((0..n).map(|r| lap_evecs[r][c]).collect());
        }

        ControlHarmonicsReport {
            receptor_channel: receptor_channel.to_string(),
            agonist_gain,
            placebo_mean_transition_energy: placebo_mean,
            modulated_mean_transition_energy: mod_mean,
            energy_landscape_delta_pct: delta_pct,
            state_names,
            transition_energy_matrix: mod_mat,
            laplacian_wavenumbers: wavenumbers,
            harmonic_modes,
            node_profiles,
        }
    }

    fn eval_transition_energy_grid(
        &self,
        a_norm: &[Vec<f64>],
        receptor_channel: &str,
        alpha_gain: f64,
    ) -> (Vec<String>, Vec<Vec<f64>>, f64, Vec<f64>) {
        let n = self.regions.len();
        let b_diag: Vec<f64> = (0..n)
            .map(|i| {
                if self.ablated_regions.contains(&self.regions[i].id) {
                    0.0
                } else {
                    let rho = self.regions[i].receptors.channel_value(receptor_channel);
                    (1.0 + alpha_gain * rho).max(0.15)
                }
            })
            .collect();

        let horizon = 8_usize;
        let mut gram = vec![vec![0.0_f64; n]; n];
        let mut cur_m = vec![vec![0.0_f64; n]; n];
        let mut a_pow = vec![vec![0.0_f64; n]; n];
        for i in 0..n {
            cur_m[i][i] = b_diag[i];
            a_pow[i][i] = 1.0;
        }

        for _step in 0..horizon {
            for i in 0..n {
                for j in i..n {
                    let dot: f64 = (0..n).map(|k| cur_m[i][k] * cur_m[j][k]).sum();
                    gram[i][j] += dot;
                    if i != j {
                        gram[j][i] += dot;
                    }
                }
            }
            let mut next_m = vec![vec![0.0_f64; n]; n];
            let mut next_a = vec![vec![0.0_f64; n]; n];
            for i in 0..n {
                for k in 0..n {
                    let aik = a_norm[i][k];
                    if aik == 0.0 {
                        continue;
                    }
                    for j in 0..n {
                        next_m[i][j] += aik * cur_m[k][j];
                        next_a[i][j] += aik * a_pow[k][j];
                    }
                }
            }
            cur_m = next_m;
            a_pow = next_a;
        }
        for (i, row) in gram.iter_mut().enumerate() {
            row[i] += 1e-2;
        }

        let mut groups: BTreeMap<String, Vec<usize>> = BTreeMap::new();
        for (i, r) in self.regions.iter().enumerate() {
            if self.ablated_regions.contains(&r.id) {
                continue;
            }
            groups.entry(r.network.clone()).or_default().push(i);
        }
        let state_names: Vec<String> = groups.keys().take(6).cloned().collect();
        let state_vecs: Vec<Vec<f64>> = state_names
            .iter()
            .map(|name| {
                let idxs = &groups[name];
                let mut v = vec![0.0_f64; n];
                let sc = 1.0 / (idxs.len().max(1) as f64).sqrt();
                for &idx in idxs {
                    v[idx] = sc;
                }
                v
            })
            .collect();

        let m_sys = state_names.len();
        let mut energy_mat = vec![vec![0.0_f64; m_sys]; m_sys];
        let mut sum_e = 0.0_f64;
        let mut cnt_e = 0_usize;
        for a_idx in 0..m_sys {
            let x0 = &state_vecs[a_idx];
            let ax0: Vec<f64> = (0..n)
                .map(|i| (0..n).map(|j| a_pow[i][j] * x0[j]).sum())
                .collect();
            for b_idx in 0..m_sys {
                if a_idx == b_idx {
                    continue;
                }
                let xt = &state_vecs[b_idx];
                let dx: Vec<f64> = (0..n).map(|i| xt[i] - ax0[i]).collect();
                let sol = solve_linear_system(&gram, &dx);
                let ev: f64 = dx.iter().zip(sol.iter()).map(|(a, b)| a * b).sum::<f64>().max(0.0);
                energy_mat[a_idx][b_idx] = ev;
                sum_e += ev;
                cnt_e += 1;
            }
        }
        let mean_e = if cnt_e > 0 { sum_e / (cnt_e as f64) } else { 0.0 };
        (state_names, energy_mat, url_safe_finite(mean_e), b_diag)
    }
}

fn url_safe_finite(v: f64) -> f64 {
    if v.is_finite() { v } else { 0.0 }
}

#[cfg(test)]
mod tests {
    use crate::presets::human_mni152_atlas;

    #[test]
    fn serotonergic_agonism_flattens_control_energy_and_harmonics_are_ordered() {
        let atlas = human_mni152_atlas();
        let rep = atlas.compute_control_and_harmonics(0.0, "ht", 1.2);
        assert!(
            rep.modulated_mean_transition_energy < rep.placebo_mean_transition_energy,
            "Expected 5-HT agonism (+1.2) to lower transition control energy"
        );
        assert!(rep.energy_landscape_delta_pct < -20.0);
        assert!(rep.laplacian_wavenumbers.len() >= 10);
        for pair in rep.laplacian_wavenumbers.windows(2) {
            assert!(pair[0] <= pair[1] + 1e-9);
        }
        assert!(rep.node_profiles.iter().all(|p| p.average_controllability >= 1.0));
        assert!(rep.node_profiles.iter().all(|p| p.modal_controllability > 0.0));
    }
}
