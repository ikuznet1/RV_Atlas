# ── RV Atlas — full bioinformatics environment ────────────────────────────────
# R 4.4.3 + RStudio Server + all figure/Shiny dependencies + Python (scVelo)
#
# Build:  docker compose build
# Run:    docker compose up
#         then open http://localhost:8787  (user: rstudio  password: see compose)
#
# Data volumes (set in docker-compose.yml):
#   ./dependencies  → /home/rstudio/RV_Atlas/dependencies   (RDS/CSV files)
#   ./output                    → /home/rstudio/output

FROM rocker/rstudio:4.4.3

LABEL maintainer="RV Atlas"
LABEL description="RV Atlas — R 4.4.3, all figure packages, Python scVelo, RStudio Server"

# ── 1. System libraries ───────────────────────────────────────────────────────
RUN apt-get update && apt-get install -y --no-install-recommends \
        # HDF5 (hdf5r, Seurat)
        libhdf5-dev \
        # Linear programming (igraph, leidenbase)
        libglpk-dev \
        # GNU Scientific Library (Seurat internals)
        libgsl-dev \
        # Network / cryptography
        libcurl4-openssl-dev \
        libssl-dev \
        libssh2-1-dev \
        libxml2-dev \
        # Graphics (ragg, textshaping, systemfonts, Cairo)
        libpng-dev \
        libfreetype6-dev \
        libharfbuzz-dev \
        libfribidi-dev \
        libfontconfig1-dev \
        libtiff-dev \
        libcairo2-dev \
        # Geospatial (sf, sp)
        libgdal-dev \
        libproj-dev \
        libgeos-dev \
        libudunits2-dev \
        # Parallel / Fortran
        libgomp1 \
        gfortran \
        # Utilities
        make \
        git \
        pandoc \
        wget \
        curl \
    && rm -rf /var/lib/apt/lists/*

# ── 2. Python (Miniconda) — scVelo / RNA velocity for Figure 8 ───────────────
# Creates conda env named 'velocity' matching the project's expected env name.
# RETICULATE_PYTHON overrides the hardcoded macOS path in the R scripts.
ENV CONDA_DIR=/opt/conda
ENV PATH="${CONDA_DIR}/bin:${PATH}"
ENV RETICULATE_PYTHON="${CONDA_DIR}/envs/velocity/bin/python"

RUN wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
        -O /tmp/miniconda.sh \
    && bash /tmp/miniconda.sh -b -p "${CONDA_DIR}" \
    && rm /tmp/miniconda.sh \
    && conda clean -afy

RUN conda create -n velocity python=3.8 -y \
    && conda run -n velocity pip install --no-cache-dir \
        scvelo \
        anndata \
        leidenalg \
        python-igraph \
        matplotlib \
        numpy \
        pandas \
        cellrank \
    && conda clean -afy

# ── 3. R package installation ─────────────────────────────────────────────────
# Packages are split into four layers ordered by install time so that
# editing app.R or adding a small package doesn't invalidate the heavy layers.

RUN Rscript -e "install.packages('pak', repos = 'https://r-lib.github.io/p/pak/stable/')"

# Layer A — Seurat ecosystem (heaviest; ~20–30 min)
RUN Rscript -e "pak::pak(c( \
        'Seurat', 'SeuratObject', 'sctransform', 'hdf5r' \
    ), ask = FALSE)"

# Layer B — Bioconductor core packages
RUN Rscript -e "pak::pak(c( \
        'bioc::DESeq2', \
        'bioc::edgeR', \
        'bioc::sva', \
        'bioc::apeglm', \
        'bioc::vsn', \
        'bioc::tximport', \
        'bioc::tximportData', \
        'bioc::biomaRt', \
        'bioc::DOSE', \
        'bioc::EnhancedVolcano', \
        'bioc::GeneOverlap', \
        'bioc::JASPAR2020', \
        'bioc::TFBSTools', \
        'bioc::motifmatchr', \
        'bioc::EnsDb.Hsapiens.v86', \
        'bioc::Nebulosa' \
    ), ask = FALSE)"

# Layer C — CRAN packages
RUN Rscript -e "pak::pak(c( \
        'shiny', 'plotly', \
        'ggplot2', 'ggpubr', 'ggfortify', 'ggrepel', 'ggeasy', \
        'patchwork', 'cowplot', \
        'dplyr', 'tidyverse', 'reshape2', 'stringr', 'forcats', \
        'RColorBrewer', 'viridis', 'colormap', \
        'igraph', 'enrichR', 'scCustomize', 'harmony', \
        'arrow', 'readxl', 'matrixStats', 'pracma', 'sf', \
        'reticulate', 'gplots', 'ashr' \
    ), ask = FALSE)"

# Layer D — GitHub packages (seurat-wrappers, hdWGCNA, monocle3)
RUN Rscript -e "pak::pak(c( \
        'satijalab/seurat-wrappers', \
        'smorabit/hdWGCNA', \
        'cole-trapnell-lab/monocle3' \
    ), ask = FALSE)"

# ── 4. Directory structure ────────────────────────────────────────────────────
# Output directory (scripts write PDFs here via ./output)
RUN mkdir -p /home/rstudio/output

# Mount points for external data drives (populated via docker-compose volumes)
#RUN mkdir -p 

# ── 5. Project source files ───────────────────────────────────────────────────
# Heavy data (dependencies/, output/, TOM/) are excluded by .dockerignore
# and mounted as volumes at runtime instead.
COPY . /home/rstudio/RV_Atlas/
RUN chown -R rstudio:rstudio /home/rstudio/RV_Atlas

# ── 6. Runtime ────────────────────────────────────────────────────────────────
# RStudio Server: http://localhost:8787  (user: rstudio)
# Set PASSWORD in docker-compose.yml or via -e PASSWORD=...
# ROOT=TRUE allows sudo inside RStudio (needed for some package installs)
ENV ROOT=TRUE

EXPOSE 8787
# rocker/rstudio's default CMD starts RStudio Server; no override needed.
