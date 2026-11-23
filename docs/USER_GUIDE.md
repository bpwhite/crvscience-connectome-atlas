# ConnectomeAtlas User Guide (`v0.8.0`)

ConnectomeAtlas is a standalone pure-Rust 3D brain & multi-scale connectome workstation supporting human macroscale parcellations, mouse mesoscale viral projections, *Drosophila* whole-brain neuropil circuits, and *C. elegans* cellular wiring diagrams.

- **License**: GNU General Public License v3.0 or later (`GPL-3.0-or-later`, see [../LICENSE](../LICENSE))
- **Methods Reference**: [Scientific Methods](SCIENTIFIC_METHODS.md)

---

## 1. Installation & Building

### Prerequisites
- Rust toolchain (`rustc` / `cargo` 1.78+ recommended, Edition `2021`).
- Standard OS graphics drivers capable of Vulkan, DX12, or Metal (used by `eframe`/`wgpu` `24.0`).

### Build & Test Commands
```bash
# Run all workspace unit tests
cargo test --workspace

# Build optimized release binaries (GUI + headless CLI)
cargo build --release --workspace

# Launch the desktop workstation
cargo run --release -p atlas-gui

# Run the headless CLI
cargo run --release -p atlas-cli -- --help
```

---

## 2. Graphical Workstation Layout

The `connectome-atlas` interface is organized into five coordinated docks:

1. **Top Menu Bar**:
   - **Organism Preset Switcher**: Instantaneously switches between *Homo sapiens* (`human_mni152`), *Mus musculus* (`mouse_ccfv3`), *Drosophila melanogaster* (`drosophila_flywire`), and *Caenorhabditis elegans* (`celegans_wormatlas`).
   - **Probe Mode Selector**: `Inspect` (`I`), `Signaling Path` (`P`), and `Lesion / Ablate` (`L`), plus a contextual `Heal All` button when lesions are active.
   - **Dock Toggles & Palette**: `Jump / Search (Ctrl+K)`, `Controls`, `2D Slices`, `Analytics`, and `Zen 3D` (`Z`).

2. **Left Controls Sidebar** (`Controls`):
   - **1. Camera & Palette**: Camera presets (`3D Orbit`, `Sagittal`, `Axial`, `Coronal`, `Reset`, `Fly Path`), 6 coloring modes (`Functional Network`, `DTI Directional RGB`, `Betweenness Heatmap`, `Harmonic Mode psi_k`, `Criticality sigma_i`, `TDA beta-1 Loops`), and sliders for `Min Weight`, `Shell Alpha`, and `River Flow`.
   - **2. 3D Shaders, Cutaway & Slices**: Toggle `Braided Fascicles` vs `Single Tube`, CSG surgical cutaways (`Intact`, `Octant`, `L-Hemi`, `Dome`), `Cutaway Offset`, internal subcortical anatomy, anatomical landmark callouts, 3-axis diffusion-tensor ellipsoids, and stereotaxic `Sagittal X` / `Coronal Y` / `Axial Z` planes.
   - **3. Receptors, Control & Harmonics**: Receptor channel (`5-HT`, `DA`, `ACh`, `Glu`, `GABA`), `Receptor Gain alpha`, `Criticality kappa`, `Harmonic psi_k`, and `Velocity (m/s)` sliders.
   - **4. Virtual Tractotomy Bundles**: Toggle individual white-matter or synaptic tracts (`Corpus Callosum`, `Arcuate Fasciculus`, etc.) to simulate focal disconnections.

3. **Central 3D Glass Brain Viewport**:
   - Renders depth-sorted Fresnel outer cortical/organism glass shells with CSG cutaway contours, sculpted internal anatomy, braided pathways with tonic laminar wave motion, 3-axis ellipsoidal fiber-tensor glyphs, and depth-sorted ROI spheres/somas.

4. **Bottom Triplanar 2D Slice Dock** (`2D Slices`):
   - Synchronized `AXIAL (XY)`, `CORONAL (XZ)`, and `SAGITTAL (YZ)` stereotaxic cross-sections with radiological orientation markers (`L/R/A/P/S/I`). Clicking or dragging crosshairs updates slice planes; clicking a projected node selects it across all docks.

5. **Right Analytics Panel** (`Analytics`):
   - **Inspector**: Shortest signaling path telemetry (hops, velocity-scaled latency in ms, weighted cost, bottleneck weight), selected ROI prose, afferent/efferent/coupled links, receptor densities, controllability, gradients `G1/G2`, fiber tensors, transcriptomics, criticality `sigma_i`, and functional domain/lesion disruption decoding.
   - **Matrix**: Community-sorted structural/synaptic adjacency matrix with network color strips, intra-community block outlines, calibrated weight colorbar legend (`0.0–1.0` + `min_weight` threshold tick), and hover/click inspection.
   - **Hubs**: Maslov–Sneppen small-worldness `sigma = (C/C0)/(L/L0)`, normalized rich-club `Phi_norm`, wiring cost saving, Betzel–Sporns generative parameters (`eta*`, `gamma*`), transcriptomic correlation `|r(mRNA, G1)|`, and degree/strength/participation rankings.
   - **Crit / TDA**: Interactive 2D plots for Kinouchi–Copelli response curve `F(r)`, coupling sweep `Delta(kappa)` (click anywhere on the curve to tune `kappa`), and descending-weight `H1` (`beta_1`) persistence barcode (hover to inspect cycle nodes; click a bar to highlight that topological cavity in 3D).
   - **Vault / JSON**: Save/load SHA-256 verified snapshots to SQLite WAL (`connectome_vault.sqlite`), open external files, or export `.swc` and `.connectome.json`.

---

## 3. Interactive Controls & Hotkeys

| Control / Key | Scope | Behavior |
| :--- | :--- | :--- |
| **Left-Drag** | 3D Viewport | Orbit camera yaw/pitch |
| **Right-Drag / Middle-Drag** | 3D Viewport | Pan camera focus target in the view plane |
| **Scroll Wheel** | 3D Viewport | Zoom camera distance (`65–520 mm`) |
| **Click** | 3D / 2D / Matrix | Select region or execute `Inspect` / `Path` / `Lesion` |
| **`1` / `2` / `3` / `4`** | Global (when not typing) | Smoothly transition to `3D Orbit` / `Left Sagittal` / `Superior Axial` / `Anterior Coronal` |
| **`R`** | Global (when not typing) | Smoothly reset camera to active species default framing |
| **`I` / `P` / `L`** | Global (when not typing) | Select `Inspect` / `Shortest Signaling Path` / `Lesion / Ablate` probe |
| **`C`** | Global (when not typing) | Cycle 3D surgical cutaway mode |
| **`Z`** | Global (when not typing) | Toggle Zen Full-Canvas 3D mode |
| **`Ctrl/Cmd + K`** | Global | Open Quick-Jump Command Palette (`ArrowUp`/`ArrowDown`, `Enter` executes) |
| **`Ctrl/Cmd + O` / `S`** | Global | Open neuroimaging file / Save `.connectome.json` envelope |
| **`F1` / `Esc`** | Global | Toggle Shortcuts & About sheet / Cascade close (`Palette` → `About` → `Flight` → `Zen`) |

---

## 4. Supported File Formats & CLI Reference

### Supported File Formats
- **`.connectome.json`**: SHA-256 verified `connectome-atlas-v1` envelope (`AtlasEnvelope`).
- **`.nii` / `.nii.gz`**: Uncompressed or gzip-compressed `NIfTI-1` header parser (`dimensions`, `voxel_spacing_mm`, `sform` RAS affine).
- **`.trk`**: TrackVis binary diffusion tractography streamline parser (`1000`-byte header + polyline vertices).
- **`.swc`**: 7-column NeuroMorpho/INCF compartmental morphology parser and generator (`atlas-cli export-swc`).
- **`.sqlite` / `.sqlite3` / `.db`**: Bundled `rusqlite` WAL vault storing JSON envelopes and audit logs with SHA-256 verification on load.

### CLI Examples
```bash
# Evaluate shortest path at 14.0 m/s conduction velocity
cargo run --release -p atlas-cli -- shortest-path --preset human_mni152 --from LANG.Wernicke_L --to LANG.Broca_L --velocity 14.0

# Inspect receptor-modulated transition energy across all human functional networks
cargo run --release -p atlas-cli -- control --preset human_mni152 --receptor ht --gain 1.2

# Run Maslov-Sneppen null & Betzel-Sporns generative benchmark
cargo run --release -p atlas-cli -- null-benchmark --preset mouse_ccfv3 --seed 42

# Save a snapshot to a temporary SQLite WAL vault and inspect it
cargo run --release -p atlas-cli -- vault --db my_study.sqlite --save-preset human_mni152 --slot baseline
```

---

## 5. Troubleshooting

- **GPU Adapter / Headless Screenshot Fails**: Ensure Vulkan/DX12/Metal drivers are installed; `render_headless_screenshot` requests a `wgpu` high-performance adapter for offscreen `Rgba8Unorm` rendering.
- **Hotkeys Do Not Respond**: Check that no text field (such as `Ctrl+K` or `Slot:`) currently has keyboard focus (`!ctx.wants_keyboard_input()` protects typing). Press `Esc` to exit modals or `Zen` mode.
- **No Signaling Path Found**: Lower `Min Weight` or click `Heal All` if lesions/tractotomies disconnected the source and target nodes.
- **SQLite Integrity Error**: Snapshots failing `AtlasEnvelope::verify_digest()` raise `VaultError::IntegrityMismatch` if modified outside the app.
