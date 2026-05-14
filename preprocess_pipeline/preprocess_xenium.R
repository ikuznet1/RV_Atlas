###############################################################################
## preprocess_xenium.R — Xenium spatial transcriptomics preprocessing (v52)
##
## Reproduces the full Xenium preprocessing pipeline: Proseg resegmentation →
## RCTD cell type annotation → SPLIT ambient RNA correction → SpaGE gene
## imputation → BPCells merge → UMAP projection → CellNEST export.
##
## Steps:
##   STEP 1 — Convert Proseg SpatialData zarr stores to Seurat RDS
##             (spatialdata2seurat.py; run as shell commands)
##   STEP 2 — Merge region sub-sections per patient
##             (transfer_seg_idents.py; run as shell commands)
##   STEP 3 — R: Load per-patient corrected RDS → Harmony integration
##             Output: Xenium_resegmented_corrected.rds
##   STEP 4 — RCTD cell type annotation (doublet mode)
##             Uses snRNA-seq Post_R3_FINAL.rds as reference
##   STEP 5 — BPCells on-disk merge of imputed RDS → Harmony integration
##             Output: Xenium_resegmented_imputed.rds
##                     Xenium_resegmented_imputed_final.rds
##   STEP 6 — MapQuery: project imputed cells onto corrected UMAP
##   STEP 7 — CellNEST export (seurat2cellnest.py; run as shell command)
##
## Working directory: /Users/ikuz/Documents/RV_Atlas/additional_scripts/
## (Proseg zarr stores assumed to be under ~/.xenium_explorer_cache/)
## Estimated RAM: 256 GB recommended (STEP 5 BPCells merge)
##
## Patient manifest:
##   1467 → XETG00217__0038291__Region_1__20241212__142808 (RVF)
##   1632 → XETG00217__0038291__Region_2__20241212__142808 (RVF)
##   1561 → XETG00217__0038290__Region_1__20241212__142808 (NF)
##   1343 → XETG00217__0038213__Region_1__20241206__182124 (RVF) [3 sections]
##   1697 → XETG00217__0038213__Region_1__20241206__182124 (NF)  [top section]
##   1691 → XETG00217__0038213__Region_1__20241206__182124 (NF)  [middle section]
##   1618 → XETG00217__0038216__Region_1__20241206__182124 (pRV) [top section]
##   1567 → XETG00217__0038216__Region_1__20241206__182124 (pRV) [bottom section]
##   1692 → XETG00217__0038290__Region_2__20241212__142808 (pRV)
###############################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(harmony)
  library(spacexr)   # RCTD — devtools::install_github('dmcable/spacexr')
  library(BPCells)
  library(dplyr)
  library(Matrix)
})

CACHE_DIR <- path.expand('~/.xenium_explorer_cache')
FUNC_DIR  <- '/Users/ikuz/Documents/XeniumWorkflow/functions'
WD        <- '/Users/ikuz/Documents/RV_Atlas/additional_scripts'

## ─────────────────────────────────────────────────────────────────────────────
## Sample manifest
## ─────────────────────────────────────────────────────────────────────────────
## Each entry: patient ID, Proseg param tag, output prefix, group
SAMPLES_CORRECTED <- list(
  list(patient = '1467', group = 'RVF',
       zarr_dir = 'proseg_proseg-output-XETG00217__0038291__Region_1__20241212__142808_2dec91fb',
       out_prefix = 'proseg_proseg-output-XETG00217__0038291__Region_1__20241212__142808_2dec91fb'),
  list(patient = '1632', group = 'RVF',
       zarr_dir = 'proseg_proseg-output-XETG00217__0038291__Region_2__20241212__142808_2dec91fb',
       out_prefix = 'proseg_proseg-output-XETG00217__0038291__Region_2__20241212__142808_2dec91fb'),
  list(patient = '1561', group = 'NF',
       zarr_dir = 'proseg_proseg-output-XETG00217__0038290__Region_1__20241212__142808_2dec91fb',
       out_prefix = 'proseg_proseg-output-XETG00217__0038290__Region_1__20241212__142808_2dec91fb'),
  list(patient = '1692', group = 'pRV',
       zarr_dir = 'proseg_proseg-output-XETG00217__0038290__Region_2__20241212__142808_2dec91fb',
       out_prefix = 'proseg_proseg-output-XETG00217__0038290__Region_2__20241212__142808_2dec91fb'),
  ## 1343 slide: three sections merged via transfer_seg_idents.py (see STEP 2)
  list(patient = '1343', group = 'RVF',
       zarr_dir = 'proseg_proseg-output-XETG00217__0038213__Region_1__20241206__182124_2dec91fb',
       out_prefix = 'proseg_proseg-output-XETG00217__0038213__Region_1__20241206__182124_2dec91fb_corrected'),
  ## 1567/1618 slide: two sections (top = 1618 pRV, bottom = 1567 pRV)
  list(patient = '1567', group = 'pRV',
       zarr_dir = 'proseg_proseg-output-XETG00217__0038216__Region_1__20241206__182124_2dec91fb',
       out_prefix = 'proseg_proseg-output-XETG00217__0038216__Region_1__20241206__182124_2dec91fb_corrected')
)

###############################################################################
## STEP 1 — Convert Proseg SpatialData zarr stores → Seurat RDS
##
## The spatialdata2seurat.py script reads a SpatialData zarr store produced
## by the Proseg resegmentation (via the Xenium Explorer app) and exports a
## Seurat-compatible RDS with the corrected cell×gene count matrix (RCTD
## and SPLIT have already been applied inside the zarr if run via the app).
##
## Two modes per patient:
##   --corrected : exports the SPLIT-corrected count layer  (X_corrected)
##   --imputed   : exports the SpaGE-imputed gene matrix    (all genes)
##
## Run these commands in a terminal from FUNC_DIR before running the R steps.
###############################################################################

cat('
# ─────────────────────────────────────────────────────────────────────────────
# STEP 1: Run spatialdata2seurat.py for each patient (shell commands)
# ─────────────────────────────────────────────────────────────────────────────
# Execute from:', FUNC_DIR, '
# Requires: Python environment with spatialdata, anndata, rpy2, scipy

cd', FUNC_DIR, '

# ── Corrected (panel-gene only, SPLIT-corrected counts) ──────────────────────
python spatialdata2seurat.py \\
  ~/.xenium_explorer_cache/proseg_proseg-output-XETG00217__0038291__Region_1__20241212__142808_2dec91fb/spatialdata_proseg.zarr \\
  --corrected --output proseg_proseg-output-XETG00217__0038291__Region_1__20241212__142808_2dec91fb_corrected.rds

python spatialdata2seurat.py \\
  ~/.xenium_explorer_cache/proseg_proseg-output-XETG00217__0038291__Region_2__20241212__142808_2dec91fb/spatialdata_proseg.zarr \\
  --corrected --output proseg_proseg-output-XETG00217__0038291__Region_2__20241212__142808_2dec91fb_corrected.rds

python spatialdata2seurat.py \\
  ~/.xenium_explorer_cache/proseg_proseg-output-XETG00217__0038290__Region_1__20241212__142808_2dec91fb/spatialdata_proseg.zarr \\
  --corrected --output proseg_proseg-output-XETG00217__0038290__Region_1__20241212__142808_2dec91fb_corrected.rds

python spatialdata2seurat.py \\
  ~/.xenium_explorer_cache/proseg_proseg-output-XETG00217__0038290__Region_2__20241212__142808_2dec91fb/spatialdata_proseg.zarr \\
  --corrected --output proseg_proseg-output-XETG00217__0038290__Region_2__20241212__142808_2dec91fb_corrected.rds

# Slide 0038213: three sections (top=1697 NF, middle=1691 NF, bottom=1343 RVF)
python spatialdata2seurat.py \\
  ~/.xenium_explorer_cache/proseg_proseg-output-XETG00217__0038213__Region_1__20241206__182124_2dec91fb/spatialdata_proseg.zarr \\
  --corrected --section top --output proseg_proseg-output-XETG00217__0038213__Region_1__20241206__182124_2dec91fb_top_corrected.rds

python spatialdata2seurat.py \\
  ~/.xenium_explorer_cache/proseg_proseg-output-XETG00217__0038213__Region_1__20241206__182124_2dec91fb/spatialdata_proseg.zarr \\
  --corrected --section middle --output proseg_proseg-output-XETG00217__0038213__Region_1__20241206__182124_2dec91fb_middle_corrected.rds

python spatialdata2seurat.py \\
  ~/.xenium_explorer_cache/proseg_proseg-output-XETG00217__0038213__Region_1__20241206__182124_2dec91fb/spatialdata_proseg.zarr \\
  --corrected --section bottom --output proseg_proseg-output-XETG00217__0038213__Region_1__20241206__182124_2dec91fb_bottom_corrected.rds

# Slide 0038216: two sections (top=1618 pRV, bottom=1567 pRV)
python spatialdata2seurat.py \\
  ~/.xenium_explorer_cache/proseg_proseg-output-XETG00217__0038216__Region_1__20241206__182124_2dec91fb/spatialdata_proseg.zarr \\
  --corrected --section top --output proseg_proseg-output-XETG00217__0038216__Region_1__20241206__182124_2dec91fb_top_corrected.rds

python spatialdata2seurat.py \\
  ~/.xenium_explorer_cache/proseg_proseg-output-XETG00217__0038216__Region_1__20241206__182124_2dec91fb/spatialdata_proseg.zarr \\
  --corrected --section bottom --output proseg_proseg-output-XETG00217__0038216__Region_1__20241206__182124_2dec91fb_bottom_corrected.rds

# ── Imputed (all 477 panel + SpaGE imputed genes) ─────────────────────────────
python spatialdata2seurat.py \\
  ~/.xenium_explorer_cache/proseg_proseg-output-XETG00217__0038291__Region_1__20241212__142808_2dec91fb/spatialdata_proseg.zarr \\
  --imputed --output proseg_proseg-output-XETG00217__0038291__Region_1__20241212__142808_2dec91fb_imputed.rds

python spatialdata2seurat.py \\
  ~/.xenium_explorer_cache/proseg_proseg-output-XETG00217__0038291__Region_2__20241212__142808_2dec91fb/spatialdata_proseg.zarr \\
  --imputed --output proseg_proseg-output-XETG00217__0038291__Region_2__20241212__142808_2dec91fb_imputed.rds

python spatialdata2seurat.py \\
  ~/.xenium_explorer_cache/proseg_proseg-output-XETG00217__0038290__Region_1__20241212__142808_2dec91fb/spatialdata_proseg.zarr \\
  --imputed --output proseg_proseg-output-XETG00217__0038290__Region_1__20241212__142808_2dec91fb_imputed.rds

python spatialdata2seurat.py \\
  ~/.xenium_explorer_cache/proseg_proseg-output-XETG00217__0038290__Region_2__20241212__142808_2dec91fb/spatialdata_proseg.zarr \\
  --imputed --output proseg_proseg-output-XETG00217__0038290__Region_2__20241212__142808_2dec91fb_imputed.rds

python spatialdata2seurat.py \\
  ~/.xenium_explorer_cache/proseg_proseg-output-XETG00217__0038213__Region_1__20241206__182124_2dec91fb/spatialdata_proseg.zarr \\
  --imputed --section top --output proseg_proseg-output-XETG00217__0038213__Region_1__20241206__182124_2dec91fb_top_imputed.rds

python spatialdata2seurat.py \\
  ~/.xenium_explorer_cache/proseg_proseg-output-XETG00217__0038213__Region_1__20241206__182124_2dec91fb/spatialdata_proseg.zarr \\
  --imputed --section middle --output proseg_proseg-output-XETG00217__0038213__Region_1__20241206__182124_2dec91fb_middle_imputed.rds

python spatialdata2seurat.py \\
  ~/.xenium_explorer_cache/proseg_proseg-output-XETG00217__0038213__Region_1__20241206__182124_2dec91fb/spatialdata_proseg.zarr \\
  --imputed --section bottom --output proseg_proseg-output-XETG00217__0038213__Region_1__20241206__182124_2dec91fb_bottom_imputed.rds

python spatialdata2seurat.py \\
  ~/.xenium_explorer_cache/proseg_proseg-output-XETG00217__0038216__Region_1__20241206__182124_2dec91fb/spatialdata_proseg.zarr \\
  --imputed --section top --output proseg_proseg-output-XETG00217__0038216__Region_1__20241206__182124_2dec91fb_top_imputed.rds

python spatialdata2seurat.py \\
  ~/.xenium_explorer_cache/proseg_proseg-output-XETG00217__0038216__Region_1__20241206__182124_2dec91fb/spatialdata_proseg.zarr \\
  --imputed --section bottom --output proseg_proseg-output-XETG00217__0038216__Region_1__20241206__182124_2dec91fb_bottom_imputed.rds
', sep = '\n')

###############################################################################
## STEP 2 — Merge multi-section slides via transfer_seg_idents.py
##
## Slides 0038213 (patients 1343/1697/1691) and 0038216 (patients 1618/1567)
## were captured as multiple tissue sections on a single slide. Each section
## was segmented independently; transfer_seg_idents.py assigns the per-section
## patient/group labels to cells by y-coordinate cutoffs derived from the
## manual nuclei assignment CSVs.
##
## Output: per-patient assignment CSV written to additional_scripts/output/
###############################################################################

cat('
# ─────────────────────────────────────────────────────────────────────────────
# STEP 2: Merge multi-section slides (shell commands; run from FUNC_DIR)
# ─────────────────────────────────────────────────────────────────────────────
cd', FUNC_DIR, '

# Slide 0038213 — three sections: top (1697 NF), middle (1691 NF), bottom (1343 RVF)
python transfer_seg_idents.py \\
  ~/.xenium_explorer_cache/proseg_proseg-output-XETG00217__0038213__Region_1__20241206__182124_2dec91fb/spatialdata_proseg.zarr \\
  output-XETG00217__0038213__Region_1__20241206__182124_top_nuclei_assignments.csv    --section top
python transfer_seg_idents.py \\
  ~/.xenium_explorer_cache/proseg_proseg-output-XETG00217__0038213__Region_1__20241206__182124_2dec91fb/spatialdata_proseg.zarr \\
  output-XETG00217__0038213__Region_1__20241206__182124_middle_nuclei_assignments.csv --section middle
python transfer_seg_idents.py \\
  ~/.xenium_explorer_cache/proseg_proseg-output-XETG00217__0038213__Region_1__20241206__182124_2dec91fb/spatialdata_proseg.zarr \\
  output-XETG00217__0038213__Region_1__20241206__182124_bottom_nuclei_assignments.csv --section bottom

# Slide 0038216 — two sections: top (1618 pRV), bottom (1567 pRV)
python transfer_seg_idents.py \\
  ~/.xenium_explorer_cache/proseg_proseg-output-XETG00217__0038216__Region_1__20241206__182124_2dec91fb/spatialdata_proseg.zarr \\
  output-XETG00217__0038216__Region_1__20241206__182124_top_nuclei_assignments.csv    --section top
python transfer_seg_idents.py \\
  ~/.xenium_explorer_cache/proseg_proseg-output-XETG00217__0038216__Region_1__20241206__182124_2dec91fb/spatialdata_proseg.zarr \\
  output-XETG00217__0038216__Region_1__20241206__182124_bottom_nuclei_assignments.csv --section bottom
', sep = '\n')

###############################################################################
## STEP 3 — Load per-patient corrected RDS → Harmony integration
##
## All 9 corrected per-patient Seurat objects are loaded from FUNC_DIR,
## labelled with patient/group metadata, subset to the 477 panel genes,
## and merged. A Harmony-corrected UMAP is computed for the corrected data.
##
## Input:  {FUNC_DIR}/{patient}_corrected.rds  (9 files from STEP 1)
## Output: Xenium_resegmented_corrected.rds
###############################################################################
message('\n\n====== STEP 3: Corrected data — Harmony integration ======')

setwd(FUNC_DIR)

CORRECTED_FILES <- c(
  'proseg_proseg-output-XETG00217__0038291__Region_1__20241212__142808_2dec91fb_corrected.rds',
  'proseg_proseg-output-XETG00217__0038291__Region_2__20241212__142808_2dec91fb_corrected.rds',
  'proseg_proseg-output-XETG00217__0038290__Region_1__20241212__142808_2dec91fb_corrected.rds',
  'proseg_proseg-output-XETG00217__0038290__Region_2__20241212__142808_2dec91fb_corrected.rds',
  'proseg_proseg-output-XETG00217__0038213__Region_1__20241206__182124_2dec91fb_top_corrected.rds',
  'proseg_proseg-output-XETG00217__0038213__Region_1__20241206__182124_2dec91fb_middle_corrected.rds',
  'proseg_proseg-output-XETG00217__0038213__Region_1__20241206__182124_2dec91fb_bottom_corrected.rds',
  'proseg_proseg-output-XETG00217__0038216__Region_1__20241206__182124_2dec91fb_top_corrected.rds',
  'proseg_proseg-output-XETG00217__0038216__Region_1__20241206__182124_2dec91fb_bottom_corrected.rds'
)
PATIENTS_CORR  <- c('1467','1632','1561','1692','1697','1691','1343','1618','1567')
GROUPS_CORR    <- c('RVF','RVF','NF','pRV','NF','NF','RVF','pRV','pRV')

## Load and label each per-patient object
obj_list <- vector('list', length(CORRECTED_FILES))
for (i in seq_along(CORRECTED_FILES)) {
  f <- CORRECTED_FILES[i]
  if (!file.exists(f)) { warning('Missing: ', f); next }
  message('Loading patient ', PATIENTS_CORR[i], ' (', GROUPS_CORR[i], ')')
  obj <- readRDS(f)
  obj$patient <- PATIENTS_CORR[i]
  obj$group   <- GROUPS_CORR[i]

  ## Retain only panel genes (first 477 features)
  panel_features <- rownames(obj)[seq_len(min(477L, nrow(obj)))]
  obj <- subset(obj, features = panel_features)

  obj_list[[i]] <- obj
  rm(obj); gc(verbose = FALSE)
}
obj_list <- Filter(Negate(is.null), obj_list)

## Merge sequentially (memory-efficient)
message('Merging ', length(obj_list), ' corrected objects ...')
xenium.obj <- Reduce(function(a, b) merge(a, b), obj_list)
rm(obj_list); gc(verbose = FALSE)

## Harmony integration
message('NormalizeData → PCA → Harmony ...')
xenium.obj <- NormalizeData(xenium.obj, verbose = FALSE)
xenium.obj <- FindVariableFeatures(xenium.obj, verbose = FALSE)
xenium.obj <- ScaleData(xenium.obj, split.by = 'patient', verbose = FALSE)
xenium.obj <- RunPCA(xenium.obj, npcs = 30, verbose = FALSE)
xenium.obj <- JoinLayers(xenium.obj)
xenium.obj <- RunHarmony(xenium.obj, group.by.vars = 'patient',
                         kmeans_init_nstart = 20, kmeans_init_iter_max = 100,
                         verbose = FALSE)
xenium.obj <- FindNeighbors(xenium.obj, dims = 1:30, reduction = 'harmony',
                             verbose = FALSE)
xenium.obj <- FindClusters(xenium.obj, resolution = 1, verbose = FALSE)
xenium.obj <- RunUMAP(xenium.obj, reduction = 'pca', dims = 1:30,
                      reduction.name = 'umap_orig', verbose = FALSE)
xenium.obj <- RunUMAP(xenium.obj, reduction = 'harmony', dims = 1:30,
                      reduction.name = 'umap', verbose = FALSE)

saveRDS(xenium.obj, './Xenium_resegmented_corrected.rds')
message('Saved Xenium_resegmented_corrected.rds  (', ncol(xenium.obj), ' cells)')


###############################################################################
## STEP 4 — RCTD cell type annotation (doublet mode)
##
## spacexr RCTD is run in doublet mode, allowing each Xenium cell to be
## assigned to up to two cell types (for heterotypic doublets / edge cells).
## The snRNA-seq Post_R3_FINAL.rds provides the reference transcriptome.
## Cell types with <25 reference cells are excluded from the reference.
##
## Annotation stored in: xenium.obj$cell_type_rctd_doublet
##                        xenium.obj$cell_type_seurat  (from snRNA label transfer)
###############################################################################
message('\n\n====== STEP 4: RCTD annotation ======')

## snRNA reference
sn_ref_path <- '/Users/ikuz/Documents/RV_Atlas/dependencies/shared/Post_R3_FINAL_with_counts.rds'
if (!file.exists(sn_ref_path)) {
  stop('snRNA reference not found: ', sn_ref_path)
}
message('Loading snRNA reference ...')
sn_ref <- readRDS(sn_ref_path)

## Build RCTD reference from snRNA (lineage-level annotation)
Idents(sn_ref) <- 'Names'
ref_counts <- GetAssayData(sn_ref, assay = 'RNA', layer = 'counts')
ref_meta   <- sn_ref@meta.data

## Remove cell types with < 25 cells
type_counts <- table(ref_meta$Names)
valid_types <- names(type_counts[type_counts >= 25])
keep_cells  <- rownames(ref_meta[ref_meta$Names %in% valid_types, ])
ref_counts  <- ref_counts[, keep_cells]
ref_labels  <- ref_meta[keep_cells, 'Names']

rctd_reference <- Reference(ref_counts,
                             cell_types = setNames(ref_labels, keep_cells),
                             nUMI       = colSums(ref_counts))
rm(sn_ref, ref_counts); gc(verbose = FALSE)
message('RCTD reference: ', length(valid_types), ' cell types, ',
        length(keep_cells), ' reference cells')

## Identify unassigned cells (no label from panel-only subclustering)
xenium.obj <- readRDS('./Xenium_resegmented_corrected.rds')

## Build RCTD query spatial object from Xenium corrected data
query_counts <- GetAssayData(xenium.obj, layer = 'counts')
coords <- data.frame(
  x = xenium.obj$x_centroid,
  y = xenium.obj$y_centroid,
  row.names = colnames(xenium.obj)
)
query_numi  <- colSums(query_counts)
query_spatial <- SpatialRNA(coords, query_counts, query_numi)

## Run RCTD (doublet mode; max_cores = number of available CPUs)
message('Running RCTD doublet mode (may take 30-60 min) ...')
rctd_query <- create.RCTD(query_spatial, rctd_reference,
                           max_cores = 16, test_mode = FALSE, UMI_min = 10)
rctd_query <- run.RCTD(rctd_query, doublet_mode = 'doublet')

## Extract results
rctd_results <- rctd_query@results
rctd_cells   <- names(rctd_results)

spot_class   <- sapply(rctd_cells, function(b)
  as.character(rctd_results[[b]]$spot_class))
first_type   <- sapply(rctd_cells, function(b)
  as.character(rctd_results[[b]]$first_type))

## Assign to Seurat metadata; reject 'reject' spot class
xenium.obj$cell_type_rctd_doublet <- 'Unassigned'
assigned_mask <- spot_class != 'reject'
xenium.obj$cell_type_rctd_doublet[rctd_cells[assigned_mask]] <-
  first_type[assigned_mask]

message('RCTD spot class distribution:')
print(table(spot_class))
message('Assigned: ', sum(assigned_mask), ' / ', length(rctd_cells))

## Store RCTD also as Seurat ident for downstream use
Idents(xenium.obj) <- 'cell_type_rctd_doublet'

saveRDS(xenium.obj, './Xenium_resegmented_corrected.rds')
message('Updated Xenium_resegmented_corrected.rds with RCTD annotations')
rm(rctd_query, query_spatial); gc(verbose = FALSE)


###############################################################################
## STEP 5 — BPCells on-disk merge of imputed RDS → Harmony integration
##
## Each imputed Seurat object contains panel genes + SpaGE-imputed genes
## (~20 GB each; 9 × 20 GB = 180 GB if all loaded simultaneously).
## BPCells converts each count matrix to a memory-mapped on-disk format
## so the merge operates on lightweight pointers rather than dense matrices.
##
## Input:  9 × {patient}_imputed.rds  (from STEP 1, in FUNC_DIR)
## Output: Xenium_resegmented_imputed.rds
##         Xenium_resegmented_imputed_final.rds  (>10 tx, non-Unassigned cells)
###############################################################################
message('\n\n====== STEP 5: BPCells imputed merge ======')

IMPUTED_FILES <- c(
  'proseg_proseg-output-XETG00217__0038291__Region_1__20241212__142808_2dec91fb_imputed.rds',
  'proseg_proseg-output-XETG00217__0038291__Region_2__20241212__142808_2dec91fb_imputed.rds',
  'proseg_proseg-output-XETG00217__0038290__Region_1__20241212__142808_2dec91fb_imputed.rds',
  'proseg_proseg-output-XETG00217__0038290__Region_2__20241212__142808_2dec91fb_imputed.rds',
  'proseg_proseg-output-XETG00217__0038213__Region_1__20241206__182124_2dec91fb_top_imputed.rds',
  'proseg_proseg-output-XETG00217__0038213__Region_1__20241206__182124_2dec91fb_middle_imputed.rds',
  'proseg_proseg-output-XETG00217__0038213__Region_1__20241206__182124_2dec91fb_bottom_imputed.rds',
  'proseg_proseg-output-XETG00217__0038216__Region_1__20241206__182124_2dec91fb_top_imputed.rds',
  'proseg_proseg-output-XETG00217__0038216__Region_1__20241206__182124_2dec91fb_bottom_imputed.rds'
)
PATIENTS_IMP <- c('1467','1632','1561','1692','1697','1691','1343','1618','1567')
GROUPS_IMP   <- c('RVF','RVF','NF','pRV','NF','NF','RVF','pRV','pRV')

bp_dir <- './bpcells_imputed'
dir.create(bp_dir, showWarnings = FALSE, recursive = TRUE)

## Pass 1: discover shared gene set across all 9 imputed files
message('=== Pass 1: discovering shared gene set ===')
gene_sets <- list()
for (i in seq_along(IMPUTED_FILES)) {
  f <- IMPUTED_FILES[i]
  if (!file.exists(f)) { warning('Missing: ', f); gene_sets[[i]] <- character(0); next }
  message('  Genes from sample ', i, ': ', PATIENTS_IMP[i])
  obj <- readRDS(f)
  gene_sets[[i]] <- rownames(obj)
  rm(obj); gc(verbose = FALSE)
}
gene_sets  <- Filter(function(x) length(x) > 0, gene_sets)
shared_genes <- Reduce(intersect, gene_sets)
message('Shared genes across all imputed samples: ', length(shared_genes))
rm(gene_sets); gc(verbose = FALSE)

## Pass 2: convert each to on-disk BPCells and retain lightweight Seurat shell
message('=== Pass 2: converting to BPCells on-disk format ===')
imp_list <- vector('list', length(IMPUTED_FILES))
for (i in seq_along(IMPUTED_FILES)) {
  f <- IMPUTED_FILES[i]
  if (!file.exists(f)) next
  message('  Processing sample ', i, ' / ', length(IMPUTED_FILES),
          ': patient ', PATIENTS_IMP[i])
  obj <- readRDS(f)
  obj$patient <- PATIENTS_IMP[i]
  obj$group   <- GROUPS_IMP[i]

  ## Transfer RCTD annotations from corrected object
  corrected <- readRDS('./Xenium_resegmented_corrected.rds')
  shared_bc  <- intersect(colnames(obj), colnames(corrected))
  obj$cell_type_rctd_doublet <- 'Unassigned'
  obj$cell_type_rctd_doublet[shared_bc] <-
    corrected$cell_type_rctd_doublet[shared_bc]
  rm(corrected); gc(verbose = FALSE)

  ## Subset to shared gene set
  obj <- subset(obj, features = shared_genes)

  ## Write counts to disk; replace in-memory matrix with on-disk pointer
  assay_name  <- DefaultAssay(obj)
  counts_path <- file.path(bp_dir, paste0('sample_', i, '_counts'))
  mat <- GetAssayData(obj, assay = assay_name, layer = 'counts')
  BPCells::write_matrix_dir(mat, dir = counts_path, overwrite = TRUE)
  obj[[assay_name]]$counts <- BPCells::open_matrix_dir(counts_path)

  imp_list[[i]] <- obj
  rm(obj, mat); gc(verbose = FALSE)
  message('  On-disk at ', counts_path)
}
imp_list <- Filter(Negate(is.null), imp_list)

## Merge (matrices are on-disk pointers — only metadata held in RAM)
message('=== Merging imputed objects ===')
xenium.imp <- merge(imp_list[[1]], imp_list[seq(2, length(imp_list))])
rm(imp_list); gc(verbose = FALSE)

## Standard Seurat v5 workflow (BPCells-backed)
message('NormalizeData → PCA → Harmony ...')
xenium.imp <- NormalizeData(xenium.imp, verbose = FALSE)
xenium.imp <- FindVariableFeatures(xenium.imp, verbose = FALSE)
xenium.imp <- ScaleData(xenium.imp, split.by = 'patient', verbose = FALSE)
xenium.imp <- RunPCA(xenium.imp, npcs = 30, verbose = FALSE)
xenium.imp <- JoinLayers(xenium.imp)
xenium.imp <- RunHarmony(xenium.imp, group.by.vars = 'patient',
                         reduction.save = 'harmony', verbose = FALSE)
xenium.imp <- FindNeighbors(xenium.imp, dims = 1:30, reduction = 'harmony',
                             verbose = FALSE)
xenium.imp <- FindClusters(xenium.imp, resolution = 1, verbose = FALSE)
xenium.imp <- RunUMAP(xenium.imp, reduction = 'pca', dims = 1:30,
                      reduction.name = 'umap_orig', verbose = FALSE)
xenium.imp <- RunUMAP(xenium.imp, reduction = 'harmony', dims = 1:30,
                      reduction.name = 'umap', verbose = FALSE)

saveRDS(xenium.imp, './Xenium_resegmented_imputed.rds')
message('Saved Xenium_resegmented_imputed.rds  (', ncol(xenium.imp), ' cells)')

## Filter: >10 panel-gene transcripts, exclude Unassigned / None / Unknown
xenium.obj <- readRDS('./Xenium_resegmented_corrected.rds')
panel_genes <- rownames(xenium.obj)[seq_len(min(477L, nrow(xenium.obj)))]
rm(xenium.obj)
high_tx  <- colSums(GetAssayData(xenium.imp, layer = 'counts')[panel_genes, ]) > 10
good_ann <- !xenium.imp$cell_type_rctd_doublet %in% c('None','Unknown','Unassigned')
xenium.imp <- xenium.imp[, high_tx & good_ann]
saveRDS(xenium.imp, './Xenium_resegmented_imputed_final.rds')
message('Saved Xenium_resegmented_imputed_final.rds  (', ncol(xenium.imp), ' cells)')


###############################################################################
## STEP 6 — MapQuery: project imputed cells onto corrected UMAP
##
## FindTransferAnchors uses the Harmony embedding of the corrected object as
## the reference; MapQuery projects imputed cells onto the reference UMAP and
## transfers cell_type_rctd_doublet labels as a validation check.
###############################################################################
message('\n\n====== STEP 6: MapQuery projection ======')

options(future.globals.maxSize = Inf)
plan('sequential')

xenium.obj <- readRDS('./Xenium_resegmented_corrected.rds')
xenium.imp <- readRDS('./Xenium_resegmented_imputed_final.rds')

## Build reference UMAP model on corrected object
xenium.obj <- RunUMAP(xenium.obj, reduction = 'pca', dims = 1:30,
                      reduction.name = 'umap', return.model = TRUE,
                      verbose = FALSE)

message('Finding transfer anchors (imputed → corrected) ...')
anchors <- FindTransferAnchors(reference = xenium.obj,
                               query     = xenium.imp,
                               dims      = 1:30,
                               reference.reduction = 'harmony',
                               verbose   = FALSE)

xenium.imp <- MapQuery(
  anchorset         = anchors,
  reference         = xenium.obj,
  query             = xenium.imp,
  refdata           = list(cell_type_rctd_doublet_pred = 'cell_type_rctd_doublet'),
  reference.reduction = 'harmony',
  reduction.model   = 'umap',
  verbose           = FALSE
)

saveRDS(xenium.imp, './Xenium_resegmented_imputed_final.rds')
message('Updated Xenium_resegmented_imputed_final.rds with MapQuery UMAP')
rm(anchors, xenium.obj, xenium.imp); gc(verbose = FALSE)


###############################################################################
## STEP 7 — CellNEST export
##
## Converts Xenium_resegmented_imputed_final.rds to an h5ad file suitable
## for CellNEST LR-pair inference. FOVs are split by patient; cell type
## annotation is sourced from cell_type_rctd_doublet.
##
## CellNEST is run externally (PyTorch model); outputs per-patient
## *_histogram_byFrequency_table_top1500.csv  and  *_allCCC.csv  files.
###############################################################################

cat('
# ─────────────────────────────────────────────────────────────────────────────
# STEP 7: CellNEST export (shell command; run from FUNC_DIR)
# ─────────────────────────────────────────────────────────────────────────────
cd', FUNC_DIR, '

# List available metadata columns (run first to confirm column names)
python seurat2cellnest.py Xenium_resegmented_imputed_final.rds --list_columns

# Export to CellNEST h5ad format (one h5ad per patient FOV)
python seurat2cellnest.py Xenium_resegmented_imputed_final.rds \\
  --output_dir cellnest_input/ \\
  --fov_col patient \\
  --cell_type_col cell_type_rctd_doublet

# Then run CellNEST per patient (requires GPU):
# for PATIENT in 1343 1392 1467 1561 1567 1618 1632 1681 1691 1692 1697; do
#   python cellnest_train.py \\
#     --input cellnest_input/Xenium_resegmented_imputed_final_${PATIENT}.h5ad \\
#     --output_dir cellnest_output/Xenium_resegmented_imputed_final_${PATIENT}/
# done
', sep = '\n')

message('\nAll Xenium preprocessing steps complete.')
message('Output files:')
message('  Xenium_resegmented_corrected.rds    — panel-gene SPLIT-corrected, RCTD annotated')
message('  Xenium_resegmented_imputed.rds      — all genes, pre-filter')
message('  Xenium_resegmented_imputed_final.rds — all genes, QC filtered, MapQuery UMAP')
message('  cellnest_input/                     — h5ad files for CellNEST')
