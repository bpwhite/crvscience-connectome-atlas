# ConnectomeAtlas Scientific Methods (`v0.8.0`)

This document specifies the mathematical formulations, units, numerical safeguards, and peer-reviewed foundations implemented in `crates/atlas-core` and `crates/atlas-render`.

- **License**: GNU General Public License v3.0 or later (`GPL-3.0-or-later`, see [../LICENSE](../LICENSE))

---

## 1. Coordinate Frames, Arc Lengths & Conduction Latency

- **Human (`human_mni152`)**: MNI152 stereotaxic millimeters (`+X` Right, `+Y` Anterior, `+Z` Superior).
- **Mouse (`mouse_ccfv3`)**: Allen Mouse Brain CCFv3 centered stereotaxic space (`10x mm` scale, Anterior Commissure origin).
- **Drosophila (`drosophila_flywire`)**: FlyWire FAFB v783 neuropil space (scaled `µm`).
- **C. elegans (`celegans_wormatlas`)**: Body-axis ganglion/ventral-cord space (`+Y` head/pharynx, `-Y` tail).

Given two connected regions `u`, `v` at coordinates `p_u`, `p_v`, curved axonal tract length `L_mm(u, v)` and conduction delay `tau_ms(u, v)` at configured conduction velocity `v_m_s` (default `9.6 m/s`, guarded by `max(0.1)` when finite and positive) are computed as:

```text
d_euclid(u, v) = || p_u - p_v ||_2
arch_factor    = 1.28  (Commissural tracts),  1.14  (nodal mean delay),  1.15  (path arc)
L_mm(u, v)     = max(2.0, d_euclid(u, v) * arch_factor)
tau_ms(u, v)   = L_mm(u, v) / v_m_s          (since 1 mm / (1 m/s) = 1 ms)
```

---

## 2. Weighted Graph Topology & Brandes Tied-Path Betweenness

Edge communication cost is inverse connection weight:

```text
c(u, v) = 1.0 / max(w_uv, 1e-6)
```

For each intact source node `s`, `dijkstra_brandes_from` solves weighted single-source shortest paths while tracking all tied predecessors `preds[w]` within numerical tolerance `1e-12` and shortest-path counts `sigma[w]`. Nodal dependencies `delta[w]` propagate in reverse settle order (`Brandes 2001`):

```text
sigma[v] = sum_{p in preds[v]} sigma[p]
delta[p] += (sigma[p] / max(sigma[w], 1.0)) * (1.0 + delta[w])
BC_norm(w) = betweenness[w] / ((N_active - 1) * (N_active - 2))      (for N_active > 2)
```

Global efficiency `E_glob`, characteristic path length `L`, participation coefficient `P_i`, and Crofts–Higham (`2009`) weighted communicability `G` (16-term Taylor series) are:

```text
E_glob = (1 / (N_active * (N_active - 1))) * sum_{s != t} (1 / dist(s, t))
L      = (1 / |{(s,t): 0 < dist(s,t) < inf}|) * sum dist(s, t)
P_i    = 1 - sum_m (k_{i,m} / s_i)^2
W_norm = D^{-1/2} W D^{-1/2}
G      = sum_{k=0..16} (W_norm^k / k!)
```

---

## 3. Maslov–Sneppen Null Ensemble & Small-Worldness

In `atlas-core::systems_neuro`, degree-preserving Maslov–Sneppen rewiring performs `4 * |E|` edge-pair swaps `(u, v, w1), (x, y, w2) -> (u, y, w1), (x, v, w2)` across `n_surr = 10` deterministic surrogates (preventing self-loops and multi-edges while preserving exact degree sequence and edge-weight multiset). True null characteristic path length `L0` and null clustering `C0` are computed directly on each surrogate:

```text
gamma  = C / max(C0, 1e-4)
lambda = max(L, 1e-6) / max(L0, 1e-6)
sigma_SW = gamma / lambda = (C / C0) / (L / L0)
```

When the observed network is fully connected (`obs_conn == true`), `L0` averages over connected surrogates (`l0_conn == true`, falling back to all surrogates if none remain connected). Wiring cost `C_mm = sum_{(u,v)} w_uv * L_mm(u, v)` is benchmarked against null wiring cost `C0_mm`, and Betzel–Sporns (`2016`) generative parameters `(eta*, gamma*)` minimize combined distance and homophily discrepancy under `E_ij ~ D_ij^{-eta} * K_ij^{gamma}`.

---

## 4. Network Control Theory & Connectome Laplacian Harmonics

### Schur-Stable Spectral Normalization
Given the symmetrized structural adjacency `W`, cyclic Jacobi eigendecomposition (`jacobi_eigen_symmetric`) yields real eigenvalues `{lambda_m}` and orthonormal eigenvectors `V`. To guarantee discrete-time Schur stability `rho(A_norm) < 1` even for signed or anti-ferromagnetic matrices where `|lambda_min| > lambda_max`:

```text
rho_max = max(0.1, max_i |lambda_i|)
scale_a = 1.0 / (1.0 + rho_max)
A_norm  = scale_a * W
```

Thus `|lambda_m * scale_a| < 1` for every eigenmode `m`, keeping `(1 - lambda_m^2 * scale_a^2) > 0`.

### Receptor-Modulated Controllability & Transition Energy
Input gain `B_i = max(0.15, 1.0 + alpha * rho_i)` couples regional neurotransmitter receptor density `rho_i` (`5-HT`, `DA`, `ACh`, `Glu`, `GABA`) with agonist/antagonist gain `alpha`:

```text
avg_ctrl(i)   = B_i^2 * sum_m ( V_{i,m}^2 / max(1e-6, 1 - (lambda_m * scale_a)^2) )
modal_ctrl(i) = sum_m ( (1 - (lambda_m * scale_a)^2) * V_{i,m}^2 )
W_T           = sum_{t=0..7} (A_norm^t B) (A_norm^t B)^T + 1e-2 I
E(x_0 -> x_T) = (x_T - A_norm^8 x_0)^T W_T^{-1} (x_T - A_norm^8 x_0)
```

All active functional systems (`groups.keys()`) in each preset—including all 9 systems in `human_mni152`—are evaluated in `transition_energy_matrix` without truncation.

### Connectome Laplacian Harmonics (`psi_k`)
The symmetric normalized Graph Laplacian `L_sym = I - D^{-1/2} W D^{-1/2}` is diagonalized via `jacobi_eigen_symmetric`. Non-trivial eigenvectors (`omega_k > 1e-4`) define spatial harmonics `psi_k` (`G1`, `G2`, ...).

---

## 5. Directed/Undirected Brain Criticality & Kinouchi–Copelli Cascades

In `atlas-core::criticality_tda`, regional E/I gain `g_EI(i) = clamp(0.65 + 0.35 * Glu_i / max(0.25, GABA_i * (0.65 + 0.35 * PVALB_i)), 0.58, 1.42)` weights structural edges while **preserving directed branching asymmetry**:

```text
G[v][u] = max(G[v][u], w_uv * sqrt(g_EI(u) * g_EI(v)))   for directed edge u -> v
G[u][v] = max(G[u][v], w_uv * sqrt(g_EI(u) * g_EI(v)))   additionally when !directed
```

`dominant_eigenpair_with_residual` computes the spectral radius `rho(G)` and Perron eigenvector using power iteration with Rayleigh-quotient residual verification:

```text
rayleigh = v^T (G v)
residual = || G v - rayleigh * v ||_2 <= 1e-10 * max(1.0, ||G v||_2)
```

If an oscillatory/periodic spectrum prevents residual convergence, a shifted Perron refinement `v_{k+1} = normalize(G v_k + shift * v_k)` restores convergence while re-evaluating the true Rayleigh quotient and residual on `G`. Pure feedforward DAGs terminate cleanly at `rho = 0.0`. Branching ratio `sigma = kappa * rho_active / max(rho_intact, 0.15)` classifies regimes into Subcritical (`sigma < 0.95`), Critical (`0.95 <= sigma <= 1.05`), and Supercritical (`sigma > 1.05`). Dynamic range across 5 decades of Poisson drive `r in [10^-5.2, 10^0.6]` is:

```text
Delta_dB = 10 * log10(r_0.9 / r_0.1)
```

---

## 6. Simplicial Clique Complex & Z2 Persistent Homology (`TDA`)

Edges are ordered in descending weight filtration `w_1 >= w_2 >= ... >= w_M`. Triangles (`2`-simplices) and tetrahedra (`3`-simplices) are enumerated alongside Union–Find cycle-birth detection with BFS minimal loop extraction. Column reduction over `Z_2` on the `2`-simplex boundary matrix `d_2` pairs birth edges `w_birth` with death triangles `w_death` (`persistence = w_birth - w_death`), yielding persistent 1D cavities `beta_1` (`PersistentCavityBar`).

---

## 7. Crossing-Fiber Diffusion Tensors & 3D Ellipsoid Glyphs

In `atlas-render::surfaces_and_glyphs`, each region accumulates a symmetric 3x3 orientation tensor from intact incident fiber unit vectors `u_e`:

```text
T_i = 0.04 * I_3 + sum_{e incident i} w_e * (u_e u_e^T)
```

Jacobi diagonalization yields ordered eigenvalues `l1 >= l2 >= l3` and orthonormal axes `(e1, e2, e3)`, Fractional Anisotropy `FA`, Westin metrics (`c_l + c_p + c_s = 1`), and Kindlmann superquadric exponents `(alpha, beta)`. `generate_tensor_ellipsoid_mesh` constructs a shaded 3D ellipsoid with semi-axes `r_k = r_max * sqrt(l_k / l_1)` aligned to `(e1, e2, e3)` and orthogonal principal-axis segments.

---

## 8. Peer-Reviewed References

1. **Brandes, U. (2001).** A faster algorithm for betweenness centrality. *Journal of Mathematical Sociology*, 25(2), 163–177.
2. **Maslov, S., & Sneppen, K. (2002).** Specificity and stability in topology of protein networks. *Science*, 296(5569), 910–913.
3. **Crofts, J. J., & Higham, D. J. (2009).** A weighted communicability measure applied to complex brain networks. *Journal of the Royal Society Interface*, 6(33), 411–414.
4. **Gu, S., et al. (2015).** Controllability of structural brain networks. *Nature Communications*, 6, 8414.
5. **Singleton, S. P., et al. (2022).** Receptor-informed network control theory links 5-HT2A agonism to flattened control energy landscapes. *Nature Communications*, 13, 5812.
6. **Margulies, D. S., et al. (2016).** Situating the default-mode network along a principal gradient of macroscale cortical organization. *PNAS*, 113(44), 12574–12579.
7. **Atasoy, S., Donnelly, I., & Pearson, J. (2016).** Human brain networks function in connectome-specific harmonic waves. *Nature Communications*, 7, 10340.
8. **Betzel, R. F., et al. (2016).** Generative models of the human connectome. *NeuroImage*, 124, 1054–1064.
9. **Kinouchi, O., & Copelli, M. (2006).** Optimal dynamical range of excitable networks at criticality. *Nature Physics*, 2(5), 348–351.
10. **Gollo, L. L., et al. (2016).** Frustrated critical brain dynamics. *PLOS Computational Biology* / PMC4860327.
11. **Sizemore, A. E., et al. (2018).** Cliques and cavities in the human connectome. *Journal of Computational Neuroscience*, 44(1), 115–145.
12. **Kindlmann, G. (2004).** Superquadric tensor glyphs. *IEEE TVCG / VisSym*.
13. **Glasser, M. F., et al. (2016).** A multi-modal parcellation of human cerebral cortex (HCP-MMP1.0). *Nature*, 536, 171–178.
14. **Dorkenwald, S., et al. (2024).** Neuronal wiring diagram of an adult brain (*Drosophila* FlyWire FAFB v783). *Nature*, 634, 124–138.
15. **Cook, S. J., et al. (2019).** Whole-animal connectomes of both *Caenorhabditis elegans* sexes. *Nature*, 571, 63–71.
