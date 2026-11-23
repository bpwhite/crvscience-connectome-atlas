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
use std::cmp::Ordering;
use std::collections::{BinaryHeap, BTreeMap, HashMap};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NodeTopologyMetrics {
    pub region_id: String,
    pub degree_total: usize,
    pub strength_in: f64,
    pub strength_out: f64,
    pub strength_total: f64,
    pub participation_coeff: f64,
    pub clustering_coeff: f64,
    pub betweenness_centrality: f64,
    pub communicability_centrality: f64,
    pub mean_conduction_delay_ms: f64,
    pub is_rich_club_hub: bool,
    pub ablated: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AtlasSummaryMetrics {
    pub active_nodes: usize,
    pub total_nodes: usize,
    pub active_edges: usize,
    pub total_edges: usize,
    pub global_efficiency: f64,
    pub characteristic_path_length: f64,
    pub mean_clustering: f64,
    pub rich_club_coeff_top_quartile: f64,
    pub mean_communicability: f64,
    pub mean_conduction_delay_ms: f64,
    pub connected_components: usize,
    pub node_metrics: Vec<NodeTopologyMetrics>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PathResult {
    pub source: String,
    pub target: String,
    pub hops: Vec<String>,
    pub total_distance_cost: f64,
    pub bottleneck_weight: f64,
    pub arc_length_mm: f64,
    pub conduction_latency_ms: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LesionDeltaReport {
    pub intact_efficiency: f64,
    pub lesioned_efficiency: f64,
    pub efficiency_retention_pct: f64,
    pub intact_components: usize,
    pub lesioned_components: usize,
}

#[derive(Clone, Copy, Debug)]
struct State {
    cost: f64,
    node: usize,
}

impl PartialEq for State {
    fn eq(&self, other: &Self) -> bool {
        self.cost.total_cmp(&other.cost) == Ordering::Equal && self.node == other.node
    }
}
impl Eq for State {}

impl Ord for State {
    fn cmp(&self, other: &Self) -> Ordering {
        other
            .cost
            .total_cmp(&self.cost)
            .then_with(|| self.node.cmp(&other.node))
    }
}

impl PartialOrd for State {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl ConnectomeAtlas {
    pub fn compute_summary_metrics(&self, min_weight: f64) -> AtlasSummaryMetrics {
        let n = self.regions.len();
        let mut idx_of = HashMap::with_capacity(n);
        for (i, r) in self.regions.iter().enumerate() {
            idx_of.insert(r.id.as_str(), i);
        }

        let mut adj: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
        let mut undir_neighbors: Vec<BTreeMap<usize, f64>> = vec![BTreeMap::new(); n];
        let mut strength_in = vec![0.0_f64; n];
        let mut strength_out = vec![0.0_f64; n];
        let mut active_edges = 0_usize;

        for edge in &self.edges {
            if !self.is_edge_intact(edge, min_weight) {
                continue;
            }
            let (Some(&u), Some(&v)) = (idx_of.get(edge.source.as_str()), idx_of.get(edge.target.as_str())) else {
                continue;
            };
            if u == v || edge.weight <= 0.0 {
                continue;
            }
            active_edges += 1;
            adj[u].push((v, edge.weight));
            strength_out[u] += edge.weight;
            strength_in[v] += edge.weight;
            let entry_uv = undir_neighbors[u].entry(v).or_insert(0.0);
            *entry_uv = entry_uv.max(edge.weight);
            let entry_vu = undir_neighbors[v].entry(u).or_insert(0.0);
            *entry_vu = entry_vu.max(edge.weight);

            if !edge.directed {
                adj[v].push((u, edge.weight));
                strength_out[v] += edge.weight;
                strength_in[u] += edge.weight;
            }
        }

        // Connected components on intact nodes (undirected projection)
        let mut visited = vec![false; n];
        let mut components = 0_usize;
        let mut active_nodes = 0_usize;
        for i in 0..n {
            if self.ablated_regions.contains(&self.regions[i].id) {
                visited[i] = true;
                continue;
            }
            active_nodes += 1;
        }
        for i in 0..n {
            if visited[i] {
                continue;
            }
            components += 1;
            let mut stack = vec![i];
            visited[i] = true;
            while let Some(curr) = stack.pop() {
                for (&nbr, _) in &undir_neighbors[curr] {
                    if !visited[nbr] {
                        visited[nbr] = true;
                        stack.push(nbr);
                    }
                }
            }
        }

        // All-pairs weighted Dijkstra (cost = 1.0 / weight) + Brandes tied-path shortest-path betweenness
        let mut betweenness = vec![0.0_f64; n];
        let mut inv_dist_sum = 0.0_f64;
        let mut finite_dist_sum = 0.0_f64;
        let mut finite_pairs = 0_usize;

        for s in 0..n {
            if self.ablated_regions.contains(&self.regions[s].id) {
                continue;
            }
            let (dist, sigma, preds, order) = dijkstra_brandes_from(s, n, &adj);
            for t in 0..n {
                if s == t || self.ablated_regions.contains(&self.regions[t].id) {
                    continue;
                }
                let d = dist[t];
                if d.is_finite() && d > 0.0 {
                    inv_dist_sum += 1.0 / d;
                    finite_dist_sum += d;
                    finite_pairs += 1;
                }
            }
            let mut delta = vec![0.0_f64; n];
            for &w in order.iter().rev() {
                let coeff = (1.0 + delta[w]) / sigma[w].max(1.0);
                for &p in &preds[w] {
                    delta[p] += sigma[p] * coeff;
                }
                if w != s && !self.ablated_regions.contains(&self.regions[w].id) {
                    betweenness[w] += delta[w];
                }
            }
        }

        let norm_pairs = if active_nodes > 1 {
            (active_nodes * (active_nodes - 1)) as f64
        } else {
            1.0
        };
        let global_efficiency = inv_dist_sum / norm_pairs;
        let characteristic_path_length = if finite_pairs > 0 {
            finite_dist_sum / (finite_pairs as f64)
        } else {
            0.0
        };
        let betw_norm = if active_nodes > 2 {
            ((active_nodes - 1) * (active_nodes - 2)) as f64
        } else {
            1.0
        };

        // Degree & Rich-club threshold (top quartile degree among active nodes)
        let degrees: Vec<usize> = (0..n).map(|i| undir_neighbors[i].len()).collect();
        let mut sorted_deg: Vec<usize> = (0..n)
            .filter(|&i| !self.ablated_regions.contains(&self.regions[i].id))
            .map(|i| degrees[i])
            .collect();
        sorted_deg.sort_unstable();
        let rc_threshold = if sorted_deg.is_empty() {
            usize::MAX
        } else {
            sorted_deg[(sorted_deg.len() * 3) / 4].max(2)
        };

        let hub_indices: Vec<usize> = (0..n)
            .filter(|&i| !self.ablated_regions.contains(&self.regions[i].id) && degrees[i] >= rc_threshold)
            .collect();
        let rich_club_coeff_top_quartile = if hub_indices.len() > 1 {
            let mut internal_edges = 0_usize;
            for &u in &hub_indices {
                for &v in &hub_indices {
                    if u < v && undir_neighbors[u].contains_key(&v) {
                        internal_edges += 1;
                    }
                }
            }
            (2.0 * internal_edges as f64) / ((hub_indices.len() * (hub_indices.len() - 1)) as f64)
        } else {
            0.0
        };

        // Per-node clustering & participation coefficient
        let mut clustering_sum = 0.0_f64;
        let mut node_metrics = Vec::with_capacity(n);
        for i in 0..n {
            let ablated = self.ablated_regions.contains(&self.regions[i].id);
            let nbrs: Vec<usize> = undir_neighbors[i].keys().copied().collect();
            let k = nbrs.len();
            let clust = if k >= 2 {
                let mut triangles = 0_usize;
                for a_idx in 0..k {
                    for b_idx in (a_idx + 1)..k {
                        if undir_neighbors[nbrs[a_idx]].contains_key(&nbrs[b_idx]) {
                            triangles += 1;
                        }
                    }
                }
                (2.0 * triangles as f64) / ((k * (k - 1)) as f64)
            } else {
                0.0
            };
            if !ablated {
                clustering_sum += clust;
            }

            // Participation coefficient P_i = 1 - sum_m (k_{i,m} / k_i)^2 across networks
            let mut comm_weights: BTreeMap<&str, f64> = BTreeMap::new();
            let mut total_undir_w = 0.0_f64;
            for (&nbr, &w) in &undir_neighbors[i] {
                *comm_weights.entry(self.regions[nbr].network.as_str()).or_insert(0.0) += w;
                total_undir_w += w;
            }
            let participation = if total_undir_w > 0.0 {
                let sq_sum: f64 = comm_weights
                    .values()
                    .map(|&cw| {
                        let frac = cw / total_undir_w;
                        frac * frac
                    })
                    .sum();
                (1.0 - sq_sum).clamp(0.0, 1.0)
            } else {
                0.0
            };

            let s_tot = strength_in[i] + strength_out[i];
            node_metrics.push(NodeTopologyMetrics {
                region_id: self.regions[i].id.clone(),
                degree_total: k,
                strength_in: strength_in[i],
                strength_out: strength_out[i],
                strength_total: s_tot,
                participation_coeff: participation,
                clustering_coeff: clust,
                betweenness_centrality: betweenness[i] / betw_norm,
                communicability_centrality: 0.0,
                mean_conduction_delay_ms: 0.0,
                is_rich_club_hub: !ablated && k >= rc_threshold,
                ablated,
            });
        }

        // Matrix-Exponential Weighted Communicability G = exp(D^{-1/2} W D^{-1/2}) & Conduction Delays
        let comm_mat = self.compute_weighted_communicability_matrix(min_weight);
        let vel_m_s = self.configured_velocity_m_s();
        let mut delay_sum = 0.0_f64;
        let mut delay_cnt = 0_usize;
        let mut node_delay_acc = vec![(0.0_f64, 0_usize); n];
        for edge in &self.edges {
            if !self.is_edge_intact(edge, min_weight) {
                continue;
            }
            if let (Some(&u), Some(&v)) = (idx_of.get(edge.source.as_str()), idx_of.get(edge.target.as_str())) {
                let ca = glam::Vec3::from_array(self.regions[u].coord);
                let cb = glam::Vec3::from_array(self.regions[v].coord);
                let arch = if edge.modality == crate::model::PathwayModality::Commissural { 1.28 } else { 1.14 };
                let len_mm = ((ca.distance(cb) as f64) * arch).max(2.0);
                let tau_ms = len_mm / vel_m_s;
                delay_sum += tau_ms;
                delay_cnt += 1;
                node_delay_acc[u].0 += tau_ms;
                node_delay_acc[u].1 += 1;
                node_delay_acc[v].0 += tau_ms;
                node_delay_acc[v].1 += 1;
            }
        }

        let mut comm_sum = 0.0_f64;
        let mut comm_pairs = 0_usize;
        for i in 0..n {
            if node_metrics[i].ablated {
                continue;
            }
            let mut row_c = 0.0_f64;
            let mut row_n = 0_usize;
            for j in 0..n {
                if i != j && !node_metrics[j].ablated {
                    row_c += comm_mat[i][j];
                    row_n += 1;
                    comm_sum += comm_mat[i][j];
                    comm_pairs += 1;
                }
            }
            node_metrics[i].communicability_centrality = if row_n > 0 { row_c / (row_n as f64) } else { 0.0 };
            node_metrics[i].mean_conduction_delay_ms = if node_delay_acc[i].1 > 0 {
                node_delay_acc[i].0 / (node_delay_acc[i].1 as f64)
            } else {
                0.0
            };
        }

        let mean_clustering = if active_nodes > 0 {
            clustering_sum / (active_nodes as f64)
        } else {
            0.0
        };

        AtlasSummaryMetrics {
            active_nodes,
            total_nodes: n,
            active_edges,
            total_edges: self.edges.len(),
            global_efficiency,
            characteristic_path_length,
            mean_clustering,
            rich_club_coeff_top_quartile,
            mean_communicability: if comm_pairs > 0 { comm_sum / (comm_pairs as f64) } else { 0.0 },
            mean_conduction_delay_ms: if delay_cnt > 0 { delay_sum / (delay_cnt as f64) } else { 0.0 },
            connected_components: components,
            node_metrics,
        }
    }

    fn configured_velocity_m_s(&self) -> f64 {
        if self.conduction_velocity_m_s.is_finite() && self.conduction_velocity_m_s > 0.0 {
            self.conduction_velocity_m_s.max(0.1)
        } else {
            9.6
        }
    }

    pub fn shortest_signaling_path(
        &self,
        source_id: &str,
        target_id: &str,
        min_weight: f64,
    ) -> Option<PathResult> {
        self.shortest_signaling_path_with_velocity(
            source_id,
            target_id,
            min_weight,
            self.configured_velocity_m_s(),
        )
    }

    pub fn shortest_signaling_path_with_velocity(
        &self,
        source_id: &str,
        target_id: &str,
        min_weight: f64,
        conduction_velocity_m_s: f64,
    ) -> Option<PathResult> {
        let vel_m_s = if conduction_velocity_m_s.is_finite() && conduction_velocity_m_s > 0.0 {
            conduction_velocity_m_s.max(0.1)
        } else {
            9.6
        };
        let n = self.regions.len();
        let mut idx_of = HashMap::with_capacity(n);
        for (i, r) in self.regions.iter().enumerate() {
            idx_of.insert(r.id.as_str(), i);
        }
        let (&s, &t) = (idx_of.get(source_id)?, idx_of.get(target_id)?);
        if self.ablated_regions.contains(source_id) || self.ablated_regions.contains(target_id) {
            return None;
        }

        let mut adj: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
        for edge in &self.edges {
            if !self.is_edge_intact(edge, min_weight) || edge.weight <= 0.0 {
                continue;
            }
            if let (Some(&u), Some(&v)) = (
                idx_of.get(edge.source.as_str()),
                idx_of.get(edge.target.as_str()),
            ) {
                adj[u].push((v, edge.weight));
                if !edge.directed {
                    adj[v].push((u, edge.weight));
                }
            }
        }

        let (dist, _, preds, _) = dijkstra_brandes_from(s, n, &adj);
        if !dist[t].is_finite() {
            return None;
        }

        let mut rev_indices = vec![t];
        let mut cur = t;
        let mut bottleneck = f64::INFINITY;
        while cur != s {
            let p = *preds[cur]
                .iter()
                .max_by(|&&a, &&b| {
                    let wa = adj[a]
                        .iter()
                        .filter(|&&(nbr, _)| nbr == cur)
                        .map(|&(_, w)| w)
                        .fold(0.0_f64, f64::max);
                    let wb = adj[b]
                        .iter()
                        .filter(|&&(nbr, _)| nbr == cur)
                        .map(|&(_, w)| w)
                        .fold(0.0_f64, f64::max);
                    wa.total_cmp(&wb).then_with(|| b.cmp(&a))
                })?;
            if let Some(w) = adj[p]
                .iter()
                .filter(|&&(nbr, _)| nbr == cur)
                .map(|&(_, w)| w)
                .reduce(f64::max)
            {
                bottleneck = bottleneck.min(w);
            }
            rev_indices.push(p);
            cur = p;
        }
        rev_indices.reverse();
        let hops: Vec<String> = rev_indices
            .into_iter()
            .map(|i| self.regions[i].id.clone())
            .collect();

        let mut arc_len_mm = 0.0_f64;
        for w_pair in hops.windows(2) {
            if let (Some(&u), Some(&v)) = (idx_of.get(w_pair[0].as_str()), idx_of.get(w_pair[1].as_str())) {
                let ca = glam::Vec3::from_array(self.regions[u].coord);
                let cb = glam::Vec3::from_array(self.regions[v].coord);
                arc_len_mm += ((ca.distance(cb) as f64) * 1.15).max(2.0);
            }
        }
        Some(PathResult {
            source: source_id.to_string(),
            target: target_id.to_string(),
            hops,
            total_distance_cost: dist[t],
            bottleneck_weight: if bottleneck.is_finite() { bottleneck } else { 0.0 },
            arc_length_mm: arc_len_mm,
            conduction_latency_ms: arc_len_mm / vel_m_s,
        })
    }

    /// Computes Crofts & Higham (2009) weighted communicability `G = exp(D^{-1/2} W D^{-1/2})` via 16-term Taylor series.
    pub fn compute_weighted_communicability_matrix(&self, min_weight: f64) -> Vec<Vec<f64>> {
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
        let mut w_tilde = vec![vec![0.0_f64; n]; n];
        for i in 0..n {
            for j in 0..n {
                if sym_w[i][j] > 0.0 && strengths[i] > 0.0 && strengths[j] > 0.0 {
                    w_tilde[i][j] = sym_w[i][j] / (strengths[i] * strengths[j]).sqrt();
                }
            }
        }
        let mut comm = vec![vec![0.0_f64; n]; n];
        let mut term = vec![vec![0.0_f64; n]; n];
        for i in 0..n {
            if !self.ablated_regions.contains(&self.regions[i].id) {
                comm[i][i] = 1.0;
                term[i][i] = 1.0;
            }
        }
        for k in 1..=16 {
            let inv_k = 1.0 / (k as f64);
            let mut next = vec![vec![0.0_f64; n]; n];
            for i in 0..n {
                for m in 0..n {
                    let tim = term[i][m];
                    if tim == 0.0 {
                        continue;
                    }
                    for j in 0..n {
                        next[i][j] += tim * w_tilde[m][j] * inv_k;
                    }
                }
            }
            for i in 0..n {
                for j in 0..n {
                    comm[i][j] += next[i][j];
                }
            }
            term = next;
        }
        comm
    }

    pub fn evaluate_lesion_impact(&self, min_weight: f64) -> LesionDeltaReport {
        let mut intact_clone = self.clone();
        intact_clone.heal_all();
        let intact = intact_clone.compute_summary_metrics(min_weight);
        let current = self.compute_summary_metrics(min_weight);
        let retention = if intact.global_efficiency > 0.0 {
            (current.global_efficiency / intact.global_efficiency) * 100.0
        } else {
            100.0
        };
        LesionDeltaReport {
            intact_efficiency: intact.global_efficiency,
            lesioned_efficiency: current.global_efficiency,
            efficiency_retention_pct: retention,
            intact_components: intact.connected_components,
            lesioned_components: current.connected_components,
        }
    }
}

/// Weighted Dijkstra + Brandes shortest-path DAG (`dist`, path counts `sigma`,
/// tied predecessors `preds`, and non-decreasing distance settle `order`).
pub(crate) fn dijkstra_brandes_from(
    start: usize,
    n: usize,
    adj: &[Vec<(usize, f64)>],
) -> (Vec<f64>, Vec<f64>, Vec<Vec<usize>>, Vec<usize>) {
    let mut dist = vec![f64::INFINITY; n];
    let mut sigma = vec![0.0_f64; n];
    let mut preds = vec![Vec::new(); n];
    let mut order = Vec::with_capacity(n);
    let mut settled = vec![false; n];
    let mut heap = BinaryHeap::new();

    dist[start] = 0.0;
    sigma[start] = 1.0;
    heap.push(State {
        cost: 0.0,
        node: start,
    });

    while let Some(State { cost, node }) = heap.pop() {
        if settled[node] || cost > dist[node] + 1e-12 {
            continue;
        }
        settled[node] = true;
        order.push(node);

        for &(next, weight) in &adj[node] {
            let step_cost = 1.0 / weight.max(1e-6);
            let next_cost = cost + step_cost;
            if next_cost + 1e-12 < dist[next] {
                dist[next] = next_cost;
                sigma[next] = sigma[node];
                preds[next].clear();
                preds[next].push(node);
                heap.push(State {
                    cost: next_cost,
                    node: next,
                });
            } else if (next_cost - dist[next]).abs() <= 1e-12 {
                sigma[next] += sigma[node];
                if !preds[next].contains(&node) {
                    preds[next].push(node);
                }
            }
        }
    }

    (dist, sigma, preds, order)
}

#[cfg(test)]
mod tests {
    use crate::model::{
        AtlasScale, BrainRegion, ConnectomeAtlas, Hemisphere, PathwayEdge, PathwayModality,
        ReceptorProfile,
    };
    use crate::presets::human_mni152_atlas;
    use std::collections::BTreeSet;

    #[test]
    fn brandes_tracks_tied_shortest_paths_symmetrically() {
        let mk_reg = |id: &str, coord: [f32; 3]| BrainRegion {
            id: id.to_string(),
            name: id.to_string(),
            network: "Net".to_string(),
            hemisphere: Hemisphere::Midline,
            coord,
            radius_mm: 4.0,
            color_hex: "#38bdf8".to_string(),
            citation_url: String::new(),
            notes: String::new(),
            what_it_does: String::new(),
            inputs_prose: String::new(),
            outputs_prose: String::new(),
            receptors: ReceptorProfile::default(),
        };
        let mk_edge = |s: &str, t: &str| PathwayEdge {
            source: s.to_string(),
            target: t.to_string(),
            weight: 1.0,
            tract_name: "T".to_string(),
            modality: PathwayModality::Association,
            directed: false,
            synapse_count_est: None,
        };
        let diamond = ConnectomeAtlas {
            id: "diamond".to_string(),
            title: "Diamond".to_string(),
            organism: "Test".to_string(),
            coordinate_space: "Test".to_string(),
            scale: AtlasScale::MacroParcellation,
            reference_summary: String::new(),
            primary_citation_url: String::new(),
            benchmark_neurons_or_parcels: String::new(),
            benchmark_synapses_or_tracts: String::new(),
            conduction_velocity_m_s: 9.6,
            regions: vec![
                mk_reg("N0", [0.0, 0.0, 0.0]),
                mk_reg("N1", [-10.0, 10.0, 0.0]),
                mk_reg("N2", [10.0, 10.0, 0.0]),
                mk_reg("N3", [0.0, 20.0, 0.0]),
            ],
            edges: vec![
                mk_edge("N0", "N1"),
                mk_edge("N0", "N2"),
                mk_edge("N1", "N3"),
                mk_edge("N2", "N3"),
            ],
            ablated_regions: BTreeSet::new(),
            severed_tracts: BTreeSet::new(),
        };
        let summary = diamond.compute_summary_metrics(0.0);
        let b0 = summary.node_metrics[0].betweenness_centrality;
        let b1 = summary.node_metrics[1].betweenness_centrality;
        let b2 = summary.node_metrics[2].betweenness_centrality;
        let b3 = summary.node_metrics[3].betweenness_centrality;
        // On a symmetric 4-cycle N0-N1-N3-N2, every vertex splits 2 tied paths (1.0 / 6.0)
        assert!((b1 - b2).abs() < 1e-9, "Tied shortest paths must yield equal betweenness: {b1} vs {b2}");
        assert!((b1 - (1.0 / 6.0)).abs() < 1e-9, "Expected 1/6, got {b1}");
        assert!((b0 - b3).abs() < 1e-9 && (b0 - b1).abs() < 1e-9);
    }

    #[test]
    fn shortest_signaling_path_scales_latency_by_configured_velocity() {
        let mut atlas = human_mni152_atlas();
        let p_base = atlas
            .shortest_signaling_path("LANG.Wernicke_L", "LANG.Broca_L", 0.0)
            .expect("path exists");
        atlas.conduction_velocity_m_s = 19.2;
        let p_fast = atlas
            .shortest_signaling_path("LANG.Wernicke_L", "LANG.Broca_L", 0.0)
            .expect("path exists");
        assert!(
            (p_fast.conduction_latency_ms * 2.0 - p_base.conduction_latency_ms).abs() < 1e-9,
            "Doubling velocity should halve latency: base={} fast={}",
            p_base.conduction_latency_ms,
            p_fast.conduction_latency_ms
        );
        let m_base = human_mni152_atlas().compute_summary_metrics(0.0);
        let m_fast = atlas.compute_summary_metrics(0.0);
        assert!((m_fast.mean_conduction_delay_ms * 2.0 - m_base.mean_conduction_delay_ms).abs() < 1e-9);
    }
}
