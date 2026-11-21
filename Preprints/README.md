# Preprints for ConnectomeAtlas (`crvscience-connectome-atlas`)

**Author:** Bryan P. White (ORCID: [0000-0003-4191-3511](https://orcid.org/0000-0003-4191-3511))  
**Hub:** https://crvscience.com/apps

This directory contains the self-contained LaTeX source (`.tex`) and single-column PDF (`.pdf`) editions for both:

1. **Companion Application Preprint (`APP-PREPRINT-21`)**
   - **Title:** ConnectomeAtlas: Standalone 3D Multi-Scale Brain Connectome, Persistent Homology, and Lesion Workstation in Rust
   - **Files:** `connectome-atlas-3d-brain-connectome.tex`, `connectome-atlas-3d-brain-connectome.pdf`
   - **Canonical Zenodo DOI:** https://doi.org/10.5281/zenodo.24359109
   - **Figshare Mirror:** https://figshare.com/articles/preprint/ConnectomeAtlas_Standalone_3D_Multi-Scale_Brain_Connectome_Persistent_Homology_and_Lesion_Workstation_in_Rust/35056468
   - **Web Preprint:** https://crvscience.com/apps/preprints/connectome-atlas-3d-brain-connectome

2. **CRV Science Public Apps Network Overview (`APP-PREPRINT-01`)**
   - **Title:** The CRV Science App Network: Twenty Open-Source Rust Programs for Science, Data, and Simulation
   - **Files:** `crv-science-public-apps-network.tex`, `crv-science-public-apps-network.pdf`
   - **Canonical Zenodo DOI:** https://doi.org/10.5281/zenodo.24222398
   - **Figshare Mirror:** https://figshare.com/articles/preprint/The_CRV_Science_App_Network_Seventeen_Open-Source_Rust_Programs_for_Science_Data_and_Simulation/35013736
   - **Web Preprint:** https://crvscience.com/apps/preprints/crv-science-public-apps-network

## Compiling Standalone LaTeX

Both `.tex` files inline their bibliography so each compiles standalone with standard `pdflatex`:

```bash
pdflatex connectome-atlas-3d-brain-connectome.tex
pdflatex crv-science-public-apps-network.tex
```
