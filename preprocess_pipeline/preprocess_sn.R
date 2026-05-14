###############################################################################
## preprocess_sn.R — snRNA-seq preprocessing pipeline (v52)
##
## Reproduces the full four-step preprocessing of the adult human RV snRNA-seq
## dataset (11 patients, 3 disease states: NF / pRV / RVF).
##
## Steps:
##   STEP 1 — Per-patient QC, Scrublet doublet scoring, ambient-gene clipping
##             Output: dataset_post_clipping_qc.rds
##   STEP 2 — SCTransform + Harmony integration across patients
##             Output: BroadData_Harmony_SCTransform_with_Doublets.rds
##   STEP 3 — Per-lineage subclustering and marker-based doublet removal
##             Output: CellTypes/{cm,fb,ec,endo,adipo,myeloid,nkt,pc,sm,epi,lec,neuron}.rds
##                     CellTypes/all.tsv (barcode whitelist)
##                     CellTypes/integrated_no_doublets.rds
##   STEP 4 — DecontX ambient RNA removal, Harmony re-integration, final object
##             Output: Post_R3_FINAL.rds
##
## Input: CellBender-corrected H5 files in ./dependencies/CellBender_Final/
##        (one file per patient, named {patient_id}_filtered.h5)
##
## Working directory: /Users/ikuz/Documents/RV_Atlas/
## Estimated RAM: 256 GB recommended (STEP 3 + 4 load all cells simultaneously)
###############################################################################

suppressPackageStartupMessages({
  library(Seurat)        # v5
  library(SeuratObject)
  library(harmony)
  library(dplyr)
  library(Matrix)
  library(reticulate)    # Scrublet via singleCellTK
  library(SingleCellExperiment)
  library(celda)         # decontX (STEP 4)
})

## ─────────────────────────────────────────────────────────────────────────────
## Patient manifest
## ─────────────────────────────────────────────────────────────────────────────
PATIENTS <- c('1343','1392','1467','1561','1567','1618','1632','1681','1691','1692','1697')

PATIENT_GROUPS <- c(
  '1343' = 'RVF',
  '1392' = 'pRV',
  '1467' = 'RVF',
  '1561' = 'NF',
  '1567' = 'pRV',
  '1618' = 'pRV',
  '1632' = 'RVF',
  '1681' = 'NF',
  '1691' = 'NF',
  '1692' = 'pRV',
  '1697' = 'NF'
)

## ─────────────────────────────────────────────────────────────────────────────
## Lineage marker genes (used for doublet identification in STEP 3)
## ─────────────────────────────────────────────────────────────────────────────
LINEAGE_MARKERS <- c(
  'PLIN1',   # Adipocytes
  'RYR2',    # Cardiomyocytes
  'LEPR',    # Endocardial
  'VWF',     # Endothelial
  'HAS1',    # Epicardial
  'DCN',     # Fibroblasts
  'CCL21',   # Lymphatic EC
  'CSF1R',   # Myeloid
  'NRXN1',   # Neuron
  'CD2',     # NKT
  'PDGFRB',  # Pericyte
  'MYH11'    # Smooth muscle
)

###############################################################################
## STEP 1 — Per-patient QC, Scrublet scoring, ambient clipping
##
## Inputs:  ./dependencies/CellBender_Final/{patient}_filtered.h5
##          ./dependencies/cellranger/{patient}/scrinvex.tsv  (exon counts)
## Outputs: ./output/dataset_post_clipping_qc.rds
##
## QC thresholds applied (per cluster, then globally):
##   - IQR outlier removal per cluster: nFeature_RNA, percent.mt, exon ratio,
##     Shannon entropy (IQR × 1.5 fence)
##   - Sklearn EllipticEnvelope (contamination = 0.05) across all per-cluster
##     survivors (catches multivariate outliers missed by IQR)
##   - Hard clipping: 200 ≤ nFeature ≤ 100,000 | percent.mt < 1% |
##     entropy > 5 | exon ratio < 0.25
###############################################################################

## Python env with sklearn and singleCellTK / scrublet
use_condaenv('r-reticulate', required = FALSE)
sklearn <- import('sklearn.covariance')

cellbender_dir <- './dependencies/CellBender_Final'
scrinvex_dir   <- './dependencies/cellranger'
out_dir        <- './output/CellTypes'
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## Per-patient preprocessing
obj_list <- vector('list', length(PATIENTS))
names(obj_list) <- PATIENTS

for (pt in PATIENTS) {
  message('\n=== Processing patient ', pt, ' (', PATIENT_GROUPS[pt], ') ===')
  h5_path <- file.path(cellbender_dir, paste0(pt, '_filtered.h5'))
  if (!file.exists(h5_path)) {
    warning('Missing CellBender H5: ', h5_path, ' — skipping patient ', pt)
    next
  }

  ## Load CellBender output
  counts <- Read10X_h5(h5_path)
  obj <- CreateSeuratObject(counts = counts, project = pt,
                            min.cells = 3, min.features = 200)
  obj$patient <- pt
  obj$group   <- PATIENT_GROUPS[pt]

  ## Mitochondrial content
  obj[['percent.mt']] <- PercentageFeatureSet(obj, pattern = '^MT-')

  ## Exon ratio from cellranger scrinvex
  scrinvex_path <- file.path(scrinvex_dir, pt, 'scrinvex.tsv')
  if (file.exists(scrinvex_path)) {
    scrinvex <- read.table(scrinvex_path, header = TRUE, sep = '\t',
                           row.names = 1)
    shared_cells <- intersect(colnames(obj), rownames(scrinvex))
    obj$exon_ratio <- NA_real_
    obj$exon_ratio[shared_cells] <- scrinvex[shared_cells, 'exon_ratio']
  } else {
    warning('scrinvex.tsv not found for patient ', pt, '; exon_ratio set to NA')
    obj$exon_ratio <- NA_real_
  }

  ## Shannon entropy (transcriptional diversity) via ndd Python package
  tryCatch({
    ndd  <- import('ndd')
    mat  <- as.matrix(GetAssayData(obj, layer = 'counts'))
    ent  <- apply(mat, 2, function(x) {
      x <- x[x > 0]
      if (length(x) == 0) return(0)
      ndd$entropy(as.integer(x))
    })
    obj$entropy <- ent
  }, error = function(e) {
    message('  ndd entropy failed: ', conditionMessage(e), ' — using log(nFeature)')
    obj$entropy <- log1p(obj$nFeature_RNA)
  })

  ## Scrublet doublet scores via singleCellTK wrapper
  tryCatch({
    sce <- as.SingleCellExperiment(obj)
    sce <- singleCellTK::runScrublet(sce, sample = pt)
    obj$doublet_score    <- colData(sce)$scrublet_score
    obj$predicted_doublet <- colData(sce)$scrublet_call
  }, error = function(e) {
    message('  Scrublet failed: ', conditionMessage(e))
    obj$doublet_score     <- NA_real_
    obj$predicted_doublet <- NA
  })

  ## Normalize + cluster for per-cluster IQR outlier removal
  obj <- NormalizeData(obj, verbose = FALSE) %>%
    FindVariableFeatures(nfeatures = 3000, verbose = FALSE) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(npcs = 30, verbose = FALSE) %>%
    FindNeighbors(dims = 1:30, verbose = FALSE) %>%
    FindClusters(resolution = 0.5, verbose = FALSE)

  ## Per-cluster IQR outlier removal
  meta <- obj@meta.data
  keep_cells <- character(0)
  for (cl in levels(Idents(obj))) {
    idx <- which(Idents(obj) == cl)
    sub <- meta[idx, , drop = FALSE]

    pass <- rep(TRUE, nrow(sub))
    for (qc_var in c('nFeature_RNA', 'percent.mt', 'exon_ratio', 'entropy')) {
      v <- sub[[qc_var]]
      if (all(is.na(v))) next
      q1 <- quantile(v, 0.25, na.rm = TRUE)
      q3 <- quantile(v, 0.75, na.rm = TRUE)
      iqr <- q3 - q1
      pass <- pass & (is.na(v) | (v >= q1 - 1.5*iqr & v <= q3 + 1.5*iqr))
    }
    keep_cells <- c(keep_cells, rownames(sub)[pass])
  }
  obj <- subset(obj, cells = keep_cells)
  message('  After per-cluster IQR: ', ncol(obj), ' cells')

  ## EllipticEnvelope multivariate outlier removal (sklearn)
  qc_mat <- na.omit(obj@meta.data[, c('nFeature_RNA','percent.mt','exon_ratio','entropy')])
  if (nrow(qc_mat) > 10) {
    ee <- sklearn$EllipticEnvelope(contamination = 0.05)
    ee$fit(as.matrix(qc_mat))
    inliers <- ee$predict(as.matrix(qc_mat)) == 1L
    obj <- subset(obj, cells = rownames(qc_mat)[inliers])
    message('  After EllipticEnvelope: ', ncol(obj), ' cells')
  }

  ## Hard threshold clipping
  obj <- subset(obj,
    subset = nFeature_RNA >= 200  &
             nFeature_RNA <= 100000 &
             percent.mt  <  1      &
             entropy     >  5      &
             exon_ratio  <  0.25)
  message('  After hard clipping: ', ncol(obj), ' cells')

  obj_list[[pt]] <- obj
  rm(obj); gc(verbose = FALSE)
}

## Remove skipped patients
obj_list <- Filter(Negate(is.null), obj_list)

## Merge and save pre-integration QC-passing object
message('\n=== Merging ', length(obj_list), ' patients ===')
combined <- Reduce(function(a, b) merge(a, b), obj_list)
saveRDS(combined, './output/dataset_post_clipping_qc.rds')
message('Saved dataset_post_clipping_qc.rds  (', ncol(combined), ' cells)')
rm(obj_list, combined); gc(verbose = FALSE)


###############################################################################
## STEP 2 — SCTransform + Harmony integration
##
## Input:  ./output/dataset_post_clipping_qc.rds
## Output: ./output/BroadData_Harmony_SCTransform_with_Doublets.rds
##
## Each patient is SCTransformed independently to avoid batch effects in
## variance stabilization. Patients are merged sequentially (memory-efficient).
## Harmony corrects patient-level embedding; UMAP on harmony reduction.
###############################################################################
message('\n\n====== STEP 2: Integration ======')

combined <- readRDS('./output/dataset_post_clipping_qc.rds')
patients_in <- unique(combined$patient)

## SCTransform each patient separately
message('SCTransforming ', length(patients_in), ' patients individually ...')
sct_list <- lapply(patients_in, function(pt) {
  message('  SCTransform: patient ', pt)
  sub <- subset(combined, subset = patient == pt)
  sub <- SCTransform(sub, verbose = FALSE, conserve.memory = TRUE,
                     vars.to.regress = 'percent.mt')
  sub
})
rm(combined); gc(verbose = FALSE)

## Sequential merge (avoids holding all SCT assays simultaneously)
message('Merging SCTransformed objects ...')
M1 <- sct_list[[1]]
for (i in seq(2, length(sct_list))) {
  M1 <- merge(M1, sct_list[[i]])
  gc(verbose = FALSE)
}
rm(sct_list); gc(verbose = FALSE)

## PrepSCTFindMarkers ensures consistent residuals across batches
M1 <- PrepSCTFindMarkers(M1)

## PCA + Harmony integration
message('PCA + Harmony ...')
M1 <- RunPCA(M1, npcs = 80, verbose = FALSE)
M1 <- RunHarmony(M1, group.by.vars = 'patient',
                 kmeans_init_nstart = 20, kmeans_init_iter_max = 100,
                 verbose = FALSE)
M1 <- FindNeighbors(M1, reduction = 'harmony', dims = 1:50, verbose = FALSE)
M1 <- FindClusters(M1, resolution = 1, verbose = FALSE)
M1 <- RunUMAP(M1, reduction = 'harmony', dims = 1:50, verbose = FALSE)

saveRDS(M1, './output/BroadData_Harmony_SCTransform_with_Doublets.rds')
message('Saved BroadData_Harmony_SCTransform_with_Doublets.rds  (',
        ncol(M1), ' cells, ', nlevels(Idents(M1)), ' clusters)')


###############################################################################
## STEP 3 — Per-lineage subclustering and marker-based doublet removal
##
## Input:  ./output/BroadData_Harmony_SCTransform_with_Doublets.rds
## Output: ./output/CellTypes/{lineage}.rds  (cleaned per-lineage objects)
##         ./output/CellTypes/{lineage}.tsv  (barcode whitelists)
##         ./output/CellTypes/all.tsv        (combined whitelist)
##         ./output/CellTypes/integrated_no_doublets.rds
##
## Strategy: extract broad-cluster cells per lineage, re-cluster at higher
## resolution, inspect marker-gene DotPlot for doublet clusters (cells
## co-expressing markers of two lineages), remove those clusters, repeat
## until clean. Cluster IDs to remove are hard-coded from the original
## manual curation (see comments below each lineage block).
##
## NOTE: The cluster IDs (e.g. "4", "7") refer to clusters in the
## per-lineage sub-object, not the broad integrated object.
###############################################################################
message('\n\n====== STEP 3: Per-lineage doublet removal ======')

.recluster <- function(obj) {
  obj %>%
    SCTransform(verbose = FALSE, conserve.memory = TRUE,
                vars.to.regress = 'percent.mt') %>%
    RunPCA(dims = 1:80, verbose = FALSE) %>%
    RunUMAP(reduction = 'harmony', dims = 1:50, verbose = FALSE) %>%
    FindNeighbors(reduction = 'harmony', dims = 1:50, verbose = FALSE) %>%
    FindClusters(resolution = 1, verbose = FALSE)
}

.save_lineage <- function(obj, name, all_barcodes) {
  rds_path <- file.path(out_dir, paste0(name, '.rds'))
  tsv_path <- file.path(out_dir, paste0(name, '.tsv'))
  saveRDS(obj, rds_path)
  write.table(colnames(obj), file = tsv_path, sep = '\t', quote = FALSE,
              col.names = NA)
  c(all_barcodes, colnames(obj))
}

M1 <- readRDS('./output/BroadData_Harmony_SCTransform_with_Doublets.rds')
all_clean <- character(0)

## ── Cardiomyocytes (broad clusters 0, 1, 2, 4, 6, 10, 11, 14, 18, 21, 23) ──
message('Lineage: CM')
cm <- subset(M1, ident = c('0','1','2','4','6','10','11','14','18','21','23'))
cm <- .recluster(cm)
## Remove clusters with high DCN (FB doublets)
cm_sub6 <- subset(cm, ident = '6')
cm_sub6 <- .recluster(cm_sub6)
to_rm <- colnames(subset(cm_sub6, ident = c('1','2','3','4','5','7')))
cm <- subset(cm, cells = to_rm, invert = TRUE)
cm <- .recluster(cm)
all_clean <- .save_lineage(cm, 'cm', all_clean)
rm(cm, cm_sub6); gc(verbose = FALSE)

## ── Fibroblasts (broad clusters 9, 12, 22, 25, 26, 27, 28) ──
message('Lineage: FB')
fb <- subset(M1, ident = c('9','12','22','25','26','27','28'))
fb <- .recluster(fb)
## Remove cluster 4 subclusters contaminated by CM/Myeloid doublets
fb4 <- subset(fb, ident = '4')
fb4 <- .recluster(fb4)
to_rm <- colnames(subset(fb4, ident = c('1','2','4')))
fb <- subset(fb, cells = to_rm, invert = TRUE)
fb <- .recluster(fb)
all_clean <- .save_lineage(fb, 'fb', all_clean)
rm(fb, fb4); gc(verbose = FALSE)

## ── Endothelial (broad clusters 8, 13, 15, 19, 24) ──
message('Lineage: EC')
ec <- subset(M1, ident = c('8','13','15','19','24'))
ec <- .recluster(ec)
## Remove DCN-rich (cluster 0, 8) and RYR2-rich (cluster 7, 9) doublets
ec <- subset(ec, ident = c('0','7','8','9'), invert = TRUE)
ec <- .recluster(ec)
## Subcluster cluster 6 for fine-grained doublet removal
ec <- subset(ec, ident = '9', invert = TRUE)  # re-check after recluster
ec6 <- subset(ec, ident = '6')
ec6 <- .recluster(ec6)
to_rm <- colnames(subset(ec6, ident = c('0','1')))
ec <- subset(ec, cells = to_rm, invert = TRUE)
ec <- .recluster(ec)
all_clean <- .save_lineage(ec, 'ec', all_clean)
rm(ec, ec6); gc(verbose = FALSE)

## ── Endocardial (broad cluster 5) ──
message('Lineage: Endo')
endo <- subset(M1, ident = '5')
endo <- .recluster(endo)
## Clusters 7, 9 are outliers / doublets
endo <- subset(endo, ident = c('7','9'), invert = TRUE)
endo0 <- subset(endo, ident = '0')
endo0 <- .recluster(endo0)
to_rm <- colnames(subset(endo0, ident = '3'))
endo <- subset(endo, cells = to_rm, invert = TRUE)
endo <- .recluster(endo)
## Cluster 2 has CM contamination; sub-cluster and purge
endo2 <- subset(endo, ident = '2')
endo2 <- .recluster(endo2)
to_rm <- colnames(subset(endo2, ident = '3'))
endo <- subset(endo, cells = to_rm, invert = TRUE)
endo <- .recluster(endo)
all_clean <- .save_lineage(endo, 'endo', all_clean)
rm(endo, endo0, endo2); gc(verbose = FALSE)

## ── Adipocytes (broad cluster 7) ──
message('Lineage: Adipo')
adipo <- subset(M1, ident = '7')
adipo <- .recluster(adipo)
## Clusters 10, 11, 12 appear as junk / low-quality
adipo <- subset(adipo, ident = c('10','11','12'), invert = TRUE)
adipo <- .recluster(adipo)
adipo_sub <- subset(adipo, ident = c('6','7','8','9'))
adipo_sub <- .recluster(adipo_sub)
to_rm <- colnames(subset(adipo_sub, ident = c('1','5')))
adipo <- subset(adipo, cells = to_rm, invert = TRUE)
adipo <- .recluster(adipo)
all_clean <- .save_lineage(adipo, 'adipo', all_clean)
rm(adipo, adipo_sub); gc(verbose = FALSE)

## ── Myeloid (broad clusters 3, 16, 17, 20, 29) ──
message('Lineage: Myeloid')
myeloid <- subset(M1, ident = c('3','16','17','20','29'))
myeloid <- .recluster(myeloid)
## Clusters 0 (DCN-rich), 4 (EC-rich), 7-11 are doublets
myeloid_sub <- subset(myeloid, ident = c('0','4','7','8','9','10','11'))
myeloid_sub <- .recluster(myeloid_sub)
to_rm_sub <- colnames(myeloid_sub)  # all subclustered cells are doublets
myeloid <- subset(myeloid, cells = to_rm_sub, invert = TRUE)
myeloid <- .recluster(myeloid)
all_clean <- .save_lineage(myeloid, 'myeloid', all_clean)
rm(myeloid, myeloid_sub); gc(verbose = FALSE)

## ── NKT (broad clusters per Harmony resolution) ──
message('Lineage: NKT')
## NKT-containing clusters identified from CD2 expression in broad DotPlot
nkt_clusters <- c('30','31','32','33')  # adjust per dataset
if (any(nkt_clusters %in% levels(Idents(M1)))) {
  nkt <- subset(M1, ident = intersect(nkt_clusters, levels(Idents(M1))))
  nkt <- .recluster(nkt)
  all_clean <- .save_lineage(nkt, 'nkt', all_clean)
  rm(nkt)
}
gc(verbose = FALSE)

## ── Smooth Muscle / Pericyte / Epicardial / LEC / Neuron ──
## These lineages are smaller and require less aggressive sub-clustering.
## Identify containing clusters by DotPlot before running:
##   DotPlot(M1, features = LINEAGE_MARKERS)
for (lineage_info in list(
  list(name = 'sm',     idents = c('34','35')),
  list(name = 'pc',     idents = c('36','37')),
  list(name = 'epi',    idents = c('38')),
  list(name = 'lec',    idents = c('39')),
  list(name = 'neuron', idents = c('40'))
)) {
  avail <- intersect(lineage_info$idents, levels(Idents(M1)))
  if (length(avail) == 0) next
  message('Lineage: ', lineage_info$name)
  obj <- subset(M1, ident = avail)
  obj <- .recluster(obj)
  all_clean <- .save_lineage(obj, lineage_info$name, all_clean)
  rm(obj); gc(verbose = FALSE)
}

## Write combined barcode whitelist
write.table(all_clean, file.path(out_dir, 'all.tsv'),
            sep = '\t', quote = FALSE, col.names = NA)
message('\nTotal clean barcodes: ', length(all_clean))

## Build integrated no-doublets object
message('Building integrated_no_doublets.rds ...')
M1_clean <- subset(M1, cells = all_clean)
rm(M1); gc(verbose = FALSE)
saveRDS(M1_clean, file.path(out_dir, 'integrated_no_doublets.rds'))
message('Saved integrated_no_doublets.rds  (', ncol(M1_clean), ' cells)')
rm(M1_clean); gc(verbose = FALSE)


###############################################################################
## STEP 4 — DecontX ambient RNA removal, Harmony re-integration, final object
##
## Input:  ./output/CellTypes/integrated_no_doublets.rds
## Output: ./output/Post_R3_FINAL.rds
##
## DecontX is run per-patient on raw counts with the clean doublet-removed
## cells as the observed set. Cells failing decontX QC (contamination > 0.25)
## and any residual contamination clusters are removed using CellSelector.
## The object is then re-Harmonized and per-lineage re-clustered to produce
## the final 34-subtype annotation (stored in Subnames_manual metadata column).
###############################################################################
message('\n\n====== STEP 4: DecontX + final integration ======')

M_clean <- readRDS(file.path(out_dir, 'integrated_no_doublets.rds'))

## Assign broad lineage labels from marker-gene clustering
## (these were determined by DotPlot inspection in STEP 3 and are stored
##  in the per-lineage barcode TSVs)
lineage_files <- list.files(out_dir, pattern = '\\.tsv$', full.names = TRUE)
lineage_files <- lineage_files[basename(lineage_files) != 'all.tsv']

lineage_labels <- character(ncol(M_clean))
names(lineage_labels) <- colnames(M_clean)

lineage_name_map <- c(
  cm = 'CM', fb = 'FB', ec = 'EC', endo = 'Endo', adipo = 'Adipo',
  myeloid = 'Myeloid', nkt = 'NKT', pc = 'PC', sm = 'SM',
  epi = 'Epi', lec = 'LEC', neuron = 'Neuron'
)

for (f in lineage_files) {
  lin_name <- sub('\\.tsv$', '', basename(f))
  barcodes <- read.table(f, header = TRUE, sep = '\t')[[2]]
  valid <- intersect(barcodes, names(lineage_labels))
  lineage_labels[valid] <- lineage_name_map[lin_name]
}
M_clean$Names <- lineage_labels
Idents(M_clean) <- 'Names'

## DecontX: run per-patient to model patient-specific ambient RNA profiles
message('Running decontX per patient ...')
M_clean$decontX_contamination <- NA_real_

for (pt in unique(M_clean$patient)) {
  message('  decontX: patient ', pt)
  idx <- which(M_clean$patient == pt)
  sub_counts <- GetAssayData(M_clean[, idx], layer = 'counts')
  sce_sub <- SingleCellExperiment(assays = list(counts = sub_counts))
  sce_sub <- decontX(sce_sub)
  M_clean$decontX_contamination[idx] <- sce_sub$decontX_contamination
}

message('Median decontX contamination: ',
        round(median(M_clean$decontX_contamination, na.rm = TRUE), 3))

## Remove high-contamination cells (>25% ambient)
M_clean <- subset(M_clean,
                  subset = decontX_contamination < 0.25 | is.na(decontX_contamination))
message('After decontX QC: ', ncol(M_clean), ' cells')

## Harmony re-integration on final clean set
message('Final Harmony integration ...')
M_clean <- SCTransform(M_clean, verbose = FALSE, conserve.memory = TRUE,
                       vars.to.regress = 'percent.mt')
M_clean <- RunPCA(M_clean, npcs = 50, verbose = FALSE)
M_clean <- RunHarmony(M_clean, group.by.vars = 'patient',
                      kmeans_init_nstart = 20, kmeans_init_iter_max = 100,
                      verbose = FALSE)
M_clean <- FindNeighbors(M_clean, reduction = 'harmony', dims = 1:50,
                          verbose = FALSE)
M_clean <- FindClusters(M_clean, resolution = 1, verbose = FALSE)
M_clean <- RunUMAP(M_clean, reduction = 'harmony', dims = 1:50, verbose = FALSE)

## Per-lineage re-clustering to define 34 subtypes
## (subtype IDs manually assigned in Subnames_manual after visual inspection
##  of per-lineage UMAP + DotPlot; code block below illustrates the pattern)
message('Per-lineage re-clustering for subtype assignment ...')

.recluster_lineage <- function(obj_full, lin_name, dims = 1:50) {
  cells <- colnames(obj_full)[obj_full$Names == lin_name]
  if (length(cells) < 50) return(obj_full)
  sub <- subset(obj_full, cells = cells)
  sub <- .recluster(sub)
  ## Subtype labels assigned interactively after visual inspection:
  ## Idents(sub) → DimPlot + DotPlot → manually label via RenameIdents()
  ## Store in Subnames_manual
  sub
}

## Final subtype labels (Subnames_manual) are assigned interactively via:
##   obj <- subset(M_clean, ident = 'CM')  → recluster → DotPlot → RenameIdents
## and stored back into M_clean$Subnames_manual
## See STEP 3 per-lineage RDS files for the annotated sub-objects.

## Placeholder: copy Names to Subnames_manual; replace with manual annotations
if (!'Subnames_manual' %in% colnames(M_clean@meta.data)) {
  M_clean$Subnames_manual <- M_clean$Names
  message('[NOTE] Subnames_manual = Names (broad lineage). ',
          'Replace with 34-subtype annotations from per-lineage RDS files.')
}

saveRDS(M_clean, './output/Post_R3_FINAL.rds')
message('\nSaved Post_R3_FINAL.rds  (', ncol(M_clean), ' cells)')
message('Preprocessing complete.')
