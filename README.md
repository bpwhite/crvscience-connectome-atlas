# ConnectomeAtlas (`connectome-atlas` v0.8.0)

Standalone `#![forbid(unsafe_code)]` core 3D Brain & Multi-Scale Connectome Atlas desktop workstation, CLI analyzer, binary neuroimaging parser (`NIfTI-1`, `TrackVis .trk`, `SWC`), and SQLite WAL vault in pure safe Rust.

- **Website**: [CRV Science Network](https://crvscience.com) · [CRV Science Public Apps](https://crvscience.com/apps)
- **License**: GNU General Public License v3.0 or later (`GPL-3.0-or-later`, see [LICENSE](LICENSE))

## Workspace Architecture

| Crate | Purpose |
| :--- | :--- |
| `crates/atlas-core` | Multi-scale stereotaxic data structures, 4 grounded reference connectomes, Dijkstra signaling paths, Rich-Club & participation metrics, and virtual tractotomy / focal lesion impact modeling. |
| `crates/atlas-io` | Safe binary readers for volumetric [NIfTI-1](https://docs.rs/nifti) (`.nii` / `.nii.gz`), [TrackVis `.trk`](https://pmc.ncbi.nlm.nih.gov/articles/PMC11631367/) diffusion MRI streamlines, `.swc` morphologies, CSV matrices, SHA-256 `.connectome.json` envelopes (synchronized with the HTML5/Three.js preview), and bundled [`rusqlite`](https://docs.rs/rusqlite/) WAL persistence. |
| `crates/atlas-render` | Orbital & stereotaxic slicing camera (`glam`), DTI directional RGB Catmull-Rom tract curve generator with calm non-strobing tonic laminar flow, depth-sorted 3D projector, and [`wgpu`](https://docs.rs/wgpu/) `24.0` WGSL shader module (`shaders/connectome.wgsl`). |
| `crates/atlas-cli` | Headless `clap` CLI for listing presets, computing weighted signaling paths, running virtual tractotomies, inspecting `.nii.gz` / `.trk` headers, and auditing the SQLite WAL vault. |
| `crates/atlas-gui` | Native [`eframe`](https://docs.rs/eframe/latest/eframe/) / `egui` `0.31` 3D desktop workstation with interactive Glass Brain viewport, Axial/Coronal/Sagittal slicing planes, `Inspect` / `Path` / `Lesion` probes, community-ordered adjacency heatmap, and SQLite snapshot manager. |

## Bundled Multi-Scale Reference Connectomes

1. **Human Macro-Connectome (`human_mni152`)** — 35 stereotaxic MNI152 ROIs (mm) across 8 canonical intrinsic networks (`DMN`, `SMN`, `VIS`, `SAL`, `DAN`, `FPN`, `LANG`, `SUBC`) and 7 major white-matter tract bundles (Corpus Callosum, Arcuate Fasciculus, Cingulum, Superior Longitudinal Fasciculus, Uncinate Fasciculus, Inferior Fronto-Occipital Fasciculus, and Corticospinal/Thalamocortical projections), contextualized against [HCP-MMP1.0](https://pmc.ncbi.nlm.nih.gov/articles/PMC4990127/) (360 cortical + 14 subcortical parcels) and [AAL3v1](https://www.oxcns.org/aal3.html) (166 parcellations).
2. **Mammalian Mesoscale Connectome (`mouse_ccfv3`)** — 24 anatomical structures mapped to the [Allen Mouse Brain Common Coordinate Framework v3 (CCFv3)](https://alleninstitute.github.io/CCF-MAP/descriptions/mouse_ccf.html) (1,675 C57BL/6J mice, 10 µm isotropic voxels, 43 isocortical areas, 314 subcortical gray matter structures, 81 fiber tracts, and 8 ventricular structures) with verified sensorimotor and thalamocortical (`VAL -> MOp/MOs`, `VPM -> SSp`, `PF -> CP`) pathways.
3. **Insect Whole-Brain Synaptic Connectome (`drosophila_flywire`)** — 20 neuropil compartments from [FlyWire FAFB v783](https://codex.flywire.ai/faq) ([Dorkenwald et al., *Nature* 2024](https://www.brainfacts.org/neuroscience-in-society/supporting-research/2024/researchers-create-first-adult-fruit-fly-brain-connectome-110724): 139,255 proofread neurons, >50,000,000 synapses, >3,600 cell types across 78 neuropils), spanning Optic Lobes (`ME`, `LO`, `AOTU`), Olfactory/Associative pathways (`AL -> MB_CA -> MB_Lobes` & `LH`), Central Complex (`PB`, `FB`, `EB`, `NO`), and `GNG`.
4. **Nematode Complete Cellular Connectome (`celegans_wormatlas`)** — 26 identifiable sensory, rich-club command interneuron (`AVA`, `AVB`, `AVD`, `AVE`, `PVC`), nerve-ring/pharyngeal bridge (`RIB`, `RIP`), and motor neurons from [WormAtlas](https://www.wormatlas.org/neuronalwiring.html) ([Varshney et al. 2011](https://arxiv.org/html/0907.2373v4); [Cook et al. 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC10983851/): 302 hermaphrodite neurons — 282 somatic + 20 pharyngeal — connected by 6,393 chemical synapses, 890 gap junctions, and 1,410 neuromuscular junctions).

## Quick Start

```bash
# Run all workspace unit tests
cargo test --workspace

# Launch the interactive 3D desktop workstation
cargo run --release -p atlas-gui

# Inspect presets & compute shortest weighted signaling path from CLI
cargo run --release -p atlas-cli -- presets
cargo run --release -p atlas-cli -- shortest-path --preset human_mni152 --from LANG.Wernicke_L --to LANG.Broca_L
cargo run --release -p atlas-cli -- control --preset human_mni152 --receptor ht --gain 1.2
cargo run --release -p atlas-cli -- harmonics --preset human_mni152 --top-modes 6
cargo run --release -p atlas-cli -- lesion --preset human_mni152 --tract "Corpus Callosum"
cargo run --release -p atlas-cli -- criticality --preset human_mni152 --kappa 1.0
cargo run --release -p atlas-cli -- topology --preset human_mni152 --min-weight 0.0
cargo run --release -p atlas-cli -- export-json --preset celegans_wormatlas --out worm.connectome.json
```

## New in v0.7.0 — 3D GPU Shaders, Braided Fascicles, Surgical Cutaways, Organoids/SWC & Tensor Glyphs

1. **Braided Multi-Strand Illuminated Fascicles (`atlas-render::generate_braided_fascicle_bundle` & `connectome.wgsl`)**: Fans white-matter pathways into braided multi-fiber bundles with callosal waist compaction, cortical insertion fanning, depth-darkening cues, and Kajiya-Kay / Banks tangent-aligned anisotropic diffuse + specular shading.
2. **Continuous Cortical Surface Heatmaps & Fresnel Rim-Lit Glass Shader (`project_parcel_scalars_to_surface` & `fs_fresnel_cortex`)**: Interpolates receptor densities, controllability, and Laplacian harmonics (`ψₖ`) smoothly across 3D cortical vertices with Fresnel rim shading so deep tracts remain crisp underneath.
3. **Interactive 3D Surgical Octant & Hemisphere Cutaways (`SurgicalCutawayMode`)**: Real-time CSG shader clipping (`AnteriorSuperiorOctant`, `LeftHemisphere`, `SuperiorAxialDome`) with pial incision margin highlighting that carves only the outer shell while preserving deep nuclei and tracts in situ.
4. **Sculpted 3D Subcortical Organoids & Branched `.SWC` Neuronal Arbors (`OrganoidClass`, `atlas-io::generate_region_swc_arbor`)**: Replaces plain spheres with C-shaped Caudate, lentiform Putamen, coiled Hippocampus, almond Amygdala, pulvinar Thalamus, foliated Cerebellum, and gyral crowns in mammals, plus true branched soma-dendrite-axon `.swc` frustum trees in *Drosophila* and *C. elegans* (`atlas-cli export-swc`).
5. **3D Crossing-Fiber Superquadric Tensor / ODF Glyphs & Pathway Fly-Through (`compute_atlas_fiber_tensors`, `FlyThroughController`)**: Computes local 3x3 symmetric orientation tensors (`FA`, Westin `c_l`, `c_p`, `c_s`, Kindlmann superquadric glyph exponents) and adds a smooth spline-guided chase camera along traced signaling pathways (`atlas-cli tensors`).

## New in v0.8.0 — Brain Criticality Dynamic-Range Observatory & Algebraic Topology (`TDA`)

1. **Brain Criticality & Kinouchi–Copelli Dynamic-Range Observatory (`atlas-core::criticality_tda::CriticalityReport`, `atlas-cli criticality`)**: Couples local Excitation/Inhibition (`Glu / (GABA * PVALB)`) gain with structural connectome spectral branching ratio `sigma(kappa)` ([Zimmern 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7902700/)) and sweeps 5 decades of Poisson stimulus drive `h` to compute the steady-state response curve `F(h)` and dynamic range `Delta = 10 log10(r_0.9 / r_0.1)` dB ([Gollo et al. 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4860327/)).
2. **Algebraic Topology Clique Complex & Exact Z2 Persistent Homology (`atlas-core::criticality_tda::TopologyTdaReport`, `atlas-cli topology`)**: Enumerates simplicial `3`-cliques (triangles) and `4`-cliques (tetrahedra) and performs exact Z2 boundary-matrix reduction across a descending weight-rank filtration to extract persistent 1D `beta_1` topological loop cavities `[w_birth, w_death)` and minimal stereotaxic cycle representatives ([Sizemore et al. 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC5769855/)).
