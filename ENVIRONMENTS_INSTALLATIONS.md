# Environments & Installations

Environment documentation for the RV_Atlas figure-generation pipeline (`Figure_*.R`,
`Supplementary_Figure_*.R`, `Generate_Supplementary_*.R`).

**Source of truth:** [`environment.yml`](environment.yml) at the repo root defines the
canonical, reproducible conda environment. This file summarizes it and records what is
actually installed on the development machine (2026-07, Apple Silicon / osx-arm64).

---

## Primary Environment (canonical / reproducible)

- **Manager**: conda (miniforge3) — *not* anaconda3
- **Env name**: `rv_atlas`
- **R version**: r-base 4.4.* (currently **R 4.4.3** in `~/miniforge3/envs/rv_atlas`)
- **Channels**: conda-forge, bioconda, defaults (**flexible** channel priority required)
- **Platform**: osx-arm64
- **Scope**: only the packages needed to render the figure panels. scVelo / cellrank /
  python-docx / Pillow are **excluded by design** — `reticulate` is loaded by some scripts
  but never actually invoked in the current figure pipeline.

### Setup from scratch

```bash
# 1. Create the env (strict channel priority CANNOT resolve libgfortran/r-base/bioconductor
#    cross-deps — flexible is mandatory)
CONDA_CHANNEL_PRIORITY=flexible conda env create -f environment.yml
conda activate rv_atlas

# 2. Post-conda: packages not available / unsolvable on conda channels for r-base 4.4
Rscript -e 'BiocManager::install(c("decoupleR","OmnipathR","WGCNA","rhdf5"), update = FALSE, ask = FALSE)'
Rscript -e 'install.packages(c("colormap","harmony","ashr","ClusterR"), repos = "https://cloud.r-project.org")'
Rscript -e 'remotes::install_github("smorabit/hdWGCNA",       ref = "dev")'
Rscript -e 'remotes::install_github("satijalab/seurat-wrappers")'
Rscript -e 'remotes::install_github("bnprks/BPCells",         subdir = "r")'
Rscript -e 'remotes::install_github("samuel-marsh/scCustomize")'   # also on CRAN
```

Why several packages are installed post-conda rather than in `environment.yml`:
- **harmony** — recent conda builds need libgfortran5 ≥15, which conflicts with the
  bioconductor stack on r-base 4.4; install from CRAN instead.
- **decoupleR + OmnipathR** (Fig 5F TF activity) — no bioconda build aligns with r-base 4.4.
- **WGCNA** — no current bioconda build for r-base 4.4.

See `environment.yml` for the full dependency list (single-cell, bulk DE, enrichment,
plotting, tidyverse, spatial/NN, and build deps). It is the authoritative manifest; this
file does not duplicate it.

### Installed package versions (`rv_atlas`, captured 2026-07)

Run figures with `~/miniforge3/envs/rv_atlas/bin/Rscript` (or `conda activate rv_atlas`):

| Package | Version | | Package | Version |
|---|---|---|---|---|
| Seurat | 5.5.0 | | DESeq2 | 1.46.0 |
| SeuratWrappers | 0.4.0 | | edgeR | 4.4.0 |
| harmony | 2.0.3 | | tximport | 1.34.0 |
| scCustomize | 3.3.0 | | sva | 3.54.0 |
| hdWGCNA | 0.4.11 | | biomaRt | 2.62.0 |
| igraph | 2.0.3 | | EnhancedVolcano | 1.24.0 |
| ggplot2 | 4.0.3 | | enrichR | 3.4 |
| ggpubr | 0.6.3 | | decoupleR | 2.12.0 |
| ggrastr | 1.0.2 | | WGCNA | 1.74 |
| openxlsx | 4.2.8.1 | | reticulate | 1.46.0 |

**Not currently installed in `rv_atlas`** (install on demand if a script needs them):
`monocle3` (trajectory), `JASPAR2020` / `TFBSTools` / `motifmatchr` (TF-motif panels),
`Nebulosa` (density FeaturePlots). These are outside the yml's figure-panel scope; add via
`BiocManager::install(...)` / `remotes::install_github("cole-trapnell-lab/monocle3")` if a
panel that uses them is regenerated.

> A second R install (framework **R 4.3.3** arm64 at `/usr/local/bin/Rscript`,
> lib `/Library/Frameworks/R.framework/Versions/4.3-arm64`) also exists on `PATH` with an
> older package set. **`rv_atlas` is canonical — prefer it.** Don't run figures against the
> bare `Rscript` on PATH expecting reproducibility.

---

## Python / reticulate

- Current **figure scripts do not invoke Python** — `reticulate` is loaded in a few scripts
  but no live `use_python`/scVelo call runs during panel generation (the `use_python(...)`
  line in `Figure_1.R:94` is commented out).
- If RNA-velocity / scVelo is ever needed (legacy `Figure_8.R` workflow), the working
  interpreter on this machine is the miniforge3 **`cmacs_traj`** env:
  `~/miniforge3/envs/cmacs_traj/bin/python` — Python 3.10.20, scvelo 0.3.4, anndata 0.10.8.
- `zellkonverter` manages its own basilisk conda env automatically
  (`~/Library/Caches/org.R-project.R/R/basilisk/.../zellkonverterAnnDataEnv-0.10.2`); no
  manual setup required.

---

## System Dependencies

Installed into the conda env via `environment.yml` (needed by several Bioconductor packages):
`pkg-config`, `hdf5`, `libxml2`, `libpng`, `libtiff`, `libjpeg-turbo`, `gsl`, `cmake`.

---

## Known discrepancies / gotchas

Stale references to fix if they cause confusion — the actual working setup is documented above:

- **`CLAUDE.md`** cites `/Users/ikuz/anaconda3/envs/velocity/bin/python` for scVelo. That env
  **does not exist**, and the manager is **miniforge3**, not anaconda3. scVelo actually lives
  in `cmacs_traj` (see Python section), and current figure scripts don't call it anyway.
- **`environment.yml`** header references a `requirements.txt.osx_arm64_lock` explicit-spec
  snapshot — that file is **not present** in the repo.
- Reticulate configs across scripts are inconsistent and mostly inert:
  `preprocess_pipeline/preprocess_sn.R:91` uses `use_condaenv('r-reticulate', required=FALSE)`
  (no such env → silent fallback); `helper_scripts/annotate_subnames.r:7` points at a pyenv
  `3.11.6` python.
- **R version**: canonical env is r-base 4.4.3; the default `Rscript` on PATH is framework
  R 4.3.3. Pin to one when reproducibility matters.
