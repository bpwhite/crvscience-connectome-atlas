# ConnectomeAtlas (`connectome-atlas` v0.8.0)

Standalone `#![forbid(unsafe_code)]`-core 3D Brain & Multi-Scale Connectome Atlas desktop workstation, headless CLI analyzer, binary neuroimaging parser (`NIfTI-1`, `TrackVis .trk`, `SWC`), and SQLite WAL vault in pure safe Rust.

- **Website**: [CRV Science Network](https://crvscience.com) · [CRV Science Public Apps](https://crvscience.com/apps)
- **License**: GNU General Public License v3.0 or later (`GPL-3.0-or-later`, see [LICENSE](LICENSE))
- **Guides**: [User Guide](docs/USER_GUIDE.md) · [Scientific Methods](docs/SCIENTIFIC_METHODS.md) · [Preprints](Preprints/README.md)

---

## Build & Run Quickstart

```bash
# Run all workspace unit tests
cargo test --workspace

# Launch the interactive 3D desktop workstation
cargo run --release -p atlas-gui

# Launch directly into a specific organism preset
cargo run --release -p atlas-gui -- --preset drosophila_flywire

# Render an offscreen 1600x1000 RGBA PNG screenshot (view: 1=Human, 2=Drosophila, 3=C. elegans, 4=Mouse)
cargo run --release -p atlas-gui -- --headless-screenshot validation/headless_screenshot.png --screenshot-view 1

# Run headless connectome analyses via atlas-cli
cargo run --release -p atlas-cli -- presets
cargo run --release -p atlas-cli -- shortest-path --preset human_mni152 --from LANG.Wernicke_L --to LANG.Broca_L --velocity 9.6
cargo run --release -p atlas-cli -- control --preset human_mni152 --receptor ht --gain 1.2
cargo run --release -p atlas-cli -- harmonics --preset human_mni152 --top-modes 6
cargo run --release -p atlas-cli -- criticality --preset human_mni152 --kappa 1.0
cargo run --release -p atlas-cli -- topology --preset human_mni152 --min-weight 0.0
cargo run --release -p atlas-cli -- null-benchmark --preset human_mni152 --seed 1337
```

---

## Interactive Controls & Hotkeys

Viewport single-key hotkeys (`1–4`, `R`, `I`, `P`, `L`, `C`, `Z`) are guarded by `!ctx.wants_keyboard_input()` so typing in `Ctrl+K` or vault text boxes never triggers camera or probe changes.

| Input / Hotkey | Action |
| :--- | :--- |
| **Left-Drag** (3D Viewport) | Orbit 3D stereotaxic camera around focus target (cancels pathway flight) |
| **Right-Drag / Middle-Drag** | Pan 3D camera target in the current view plane |
| **Mouse Scroll** | Zoom camera distance (`65–520 mm`) |
| **Click Node** (3D / 2D Slices / Matrix) | Execute active Probe Mode (`Inspect`, `Signaling Path`, or `Lesion / Ablate`) |
| **`1` / `2` / `3` / `4`** | Smoothly interpolate camera to `3D Orbit` / `Left Sagittal` / `Superior Axial` / `Anterior Coronal` |
| **`R`** | Smoothly reset camera to active species default framing |
| **`I` / `P` / `L`** | Switch Probe Mode to `Inspect` / `Shortest Signaling Path` / `Lesion / Ablate` |
| **`C`** | Cycle 3D surgical cutaway (`Intact` → `Octant` → `Left Hemisphere` → `Superior Dome`) |
| **`Z`** | Toggle Zen Full-Canvas 3D mode (collapses sidebars and 2D slice dock) |
| **`Ctrl/Cmd + K`** | Open Quick-Jump Command Palette (`ArrowUp`/`ArrowDown` select, `Enter` executes) |
| **`Ctrl/Cmd + O` / `Ctrl/Cmd + S`** | Open `.connectome.json`/`.swc`/`.trk`/`.nii(.gz)` or export `.connectome.json` |
| **`F1` / `Esc`** | Toggle Shortcuts & About sheet / Cascade close (`Palette` → `About` → `Flight` → `Zen`) |
| **Drag & Drop onto Window** | Ingest `.connectome.json`, `.swc`, `.trk`, `.nii`, or `.nii.gz` files directly |

---

## CLI Flags & Subcommands

### `connectome-atlas` (`atlas-gui`) Binary Flags & Headless Subcommands

| Flag / Subcommand | Arguments & Defaults | Purpose |
| :--- | :--- | :--- |
| `--db <PATH>` | Default `connectome_vault.sqlite` | Path to local SQLite WAL snapshot database |
| `--preset <ID>` | `human_mni152`, `mouse_ccfv3`, `drosophila_flywire`, `celegans_wormatlas` | Initial connectome preset loaded on startup |
| `--headless-screenshot <PATH>` | `--screenshot-view <1\|2\|3\|4>` (default `1`) | Render offscreen `1600x1000` RGBA PNG without a GUI window |
| `presets` | — | List bundled multi-scale connectome presets |
| `metrics` | `--preset <ID>` | Output JSON structural network summary metrics |
| `criticality` | `--preset <ID> --kappa <F64>` (default `1.0`) | Output JSON Kinouchi–Copelli branching ratio `sigma` and dynamic range `Delta` |
| `topology` | `--preset <ID> --min-weight <F64>` (default `0.0`) | Output JSON clique complex & Z2 `beta_1` persistent homology report |
| `path` | `--preset <ID> --from <ID> --to <ID> --velocity <F64>` (default `9.6`) | Output JSON shortest weighted signaling path & velocity-scaled latency |
| `lesion` | `--preset <ID> [--tract <NAME>] [--region <ID>]` | Output JSON `LesionDeltaReport` efficiency retention & component count |
| `export-json` | `--preset <ID> --out <PATH>` | Write SHA-256 verified `.connectome.json` envelope |
| `export-swc` | `--preset <ID> --region <ID> --out <PATH>` | Export standardized 7-column NeuroMorpho `.swc` 3D arbor |

### `atlas-cli` Subcommands

| Subcommand | Key Flags | Purpose |
| :--- | :--- | :--- |
| `presets` | — | List all 4 bundled multi-scale reference presets |
| `search` | `--preset --query --limit` | Search the `Ctrl+K` region & command index |
| `inspect` | `--preset --min-weight` | Compute global efficiency, Brandes betweenness, rich-club `Phi`, and delays |
| `inspect-region` | `--preset --region` | Print functional prose, receptor profile, and afferent/efferent pathways |
| `shortest-path` | `--preset --from --to --velocity` | Compute weighted Dijkstra path, arc length, and latency at configured `v` (m/s) |
| `control` | `--preset --receptor --gain` | Compute Schur-normalized receptor-informed controllability and transition energy |
| `harmonics` | `--preset --top-modes` | Print normalized Graph Laplacian harmonic wavenumbers `omega_k` |
| `criticality` | `--preset --kappa` | Compute directed/undirected E/I branching ratio `sigma` and dynamic range `Delta` (dB) |
| `topology` | `--preset --min-weight` | Enumerate `3`/`4`-cliques and exact Z2 persistent `beta_1` cavities |
| `synergy` | `--preset` | Compute PhiID synergy-redundancy rank, metabolic economy, and `|r(SST-PVALB, G1)|` |
| `null-benchmark` | `--preset --seed` | Run Maslov–Sneppen null ensemble `(C/C0)/(L/L0)` & Betzel–Sporns generative fit |
| `decode-lesion` | `--preset [--ablate] [--sever-tract]` | Decode cognitive/motor domain recruitment or `%` lesion disruption |
| `lesion` | `--preset [--tract] [--region]` | Evaluate virtual tractotomy or focal ablation impact on global efficiency |
| `tensors` | `--preset` | Compute 3x3 crossing-fiber orientation tensors (`FA`, Westin `c_l/c_p/c_s`, superquadrics) |
| `export-swc` / `export-json` | `--preset [--region] --out` | Export `.swc` neuronal morphology or SHA-256 `.connectome.json` |
| `inspect-nifti` / `inspect-trk` | `--file` | Inspect `NIfTI-1` (`.nii`/`.nii.gz`) or TrackVis `.trk` binary headers |
| `vault` | `--db [--save-preset] --slot` | Save/list SHA-256 verified snapshots in SQLite WAL |

---

## Scientific Summary & Architecture

| Crate / Preset | Scope & Implementation Highlights |
| :--- | :--- |
| `crates/atlas-core` | `#![forbid(unsafe_code)]`; stereotaxic data model, 4 reference connectomes, Brandes tied-path weighted Dijkstra betweenness, velocity-configured signaling latency (`conduction_velocity_m_s`), Crofts–Higham 16-term matrix-exponential communicability, Schur-normalized (`1 + max_i \|lambda_i\|`) network control across all functional networks, Jacobi symmetric eigensolver for Graph Laplacian harmonics (`psi_k`), residual-checked directed/undirected Kinouchi–Copelli criticality (`sigma`, `Delta` dB), exact Z2 `beta_1` persistent homology, and Maslov–Sneppen true-`L0` small-world & Betzel–Sporns null benchmarks. |
| `crates/atlas-io` | `#![forbid(unsafe_code)]`; binary `NIfTI-1` (`.nii`/`.nii.gz`), TrackVis `.trk`, `.swc` arbor parser/exporter, SHA-256 `AtlasEnvelope` JSON, and bundled `rusqlite` WAL vault. |
| `crates/atlas-render` | `#![forbid(unsafe_code)]`; `OrbitCamera` with right-drag panning, smooth exponential preset interpolation, per-frame precomputed `view_proj_matrix()`, braided multi-strand fascicles, CSG surgical cutaways, sculpted organoid meshes, 3-axis ellipsoidal diffusion-tensor glyphs (`generate_tensor_ellipsoid_mesh`), 6 `ColorSchemeMode` overlays (`FunctionalNetwork`, `DtiDirectionalRgb`, `BetweennessHeatmap`, `HarmonicModePsi`, `CriticalitySusceptibility`, `TdaBeta1Loops`), and `wgpu` 24.0 WGSL shaders. |
| `crates/atlas-cli` | `#![forbid(unsafe_code)]`; 19 headless subcommands for reproducible batch analysis and format inspection. |
| `crates/atlas-gui` | `#![deny(unsafe_code)]` (localized Windows `AttachConsole` helper); `eframe`/`egui` `0.31` 3D Glass Brain viewport, cached `LesionDeltaReport`, `Ctrl+K` fuzzy command palette with `Enter` execution, community-ordered matrix with colorbar legend, interactive `Crit / TDA` 2D plots (`F(r)`, `kappa` sweep, `beta_1` barcode), and triplanar 2D slice dock. |

### Bundled Multi-Scale Reference Connectomes

1. **`human_mni152`** (*Homo sapiens*): 40 stereotaxic MNI152 ROIs (mm) across 9 functional systems (`Default Mode`, `Sensorimotor`, `Visual`, `Salience`, `Dorsal Attention`, `Frontoparietal`, `Language`, `Subcortical & Limbic`, `Cerebellar & Brainstem`) and 62 pathways across 8 tract systems.
2. **`mouse_ccfv3`** (*Mus musculus* C57BL/6J): 28 Allen CCFv3 structures across 9 subdivisions and 39 sensorimotor, thalamocortical, basal ganglia, hippocampal, and olfactory pathways.
3. **`drosophila_flywire`** (*Drosophila melanogaster*): 22 FlyWire FAFB v783 neuropils across 5 systems and 32 directed/commissural pathways.
4. **`celegans_wormatlas`** (*Caenorhabditis elegans*): 26 WormAtlas hermaphrodite neurons across 4 functional classes and 35 chemical/gap-junction pathways.

---

## Validation Metrics

Observed workspace artifacts and test suite status (from `validation/connectome_atlas_validation_summary.json` and `cargo test --workspace`):

- **Packaged Windows Binary**: `20_Connectome-Atlas.exe` (`v0.8.0`, `4` reference presets).
- **Headless Offscreen Render Captures (`1600x1000` RGBA PNG)**:
  - View 1 (`human_mni152`): `validation/headless_screenshot.png` (`720,314` bytes)
  - View 2 (`drosophila_flywire`): `validation/headless_screenshot_2.png` (`397,657` bytes)
  - View 3 (`celegans_wormatlas`): `validation/headless_screenshot_3.png` (`403,864` bytes)
- **Workspace Test Suite (`cargo test --workspace`)**: **`22 passed; 0 failed`** (`atlas_core`: `11 passed; 0 failed`, `atlas_gui`: `6 passed; 0 failed`, `atlas_render`: `3 passed; 0 failed`, `atlas_io`: `2 passed; 0 failed`, `atlas_cli`: `0 tests`; baseline snapshot in `validation/connectome_atlas_validation_summary.json` recorded `12/12` passing prior to these upgrades).

---

## Citation & License

This software is licensed under the **GNU General Public License v3.0 or later (`GPL-3.0-or-later`)**. See [LICENSE](LICENSE) and [CITATION.cff](CITATION.cff) for full citation metadata (`DOI: 10.5281/zenodo.24359109`).
