###############################################################################
## Supplementary Figure 2 (v53 layout) — Cardiomyocyte subclustering and
## spatial validation. Supports the 3-axis stress-reactivation narrative.
##
## Panels:
##   (A) UMAP of CM subclustering (10 subtypes: CM_Baseline, CM_BetaMHC,
##       CM_XIRP, CM_RORA, CM_HAND2, CM_HTR4, CM_HMGCS2, CM_KCNJ3, CM_MYH6,
##       CM_NPP) — colored by Subnames, patient, and disease group
##   (B) Marker-gene dot plot per subtype (12 genes; default Seurat blue)
##   (C) CM subtype proportions across disease states (NF / pRV / RVF) —
##       three subtypes go up (CM_BetaMHC, CM_HAND2, CM_NPP) and several
##       specialized subtypes go down
##   (D) GO BP enrichment per subtype (per-cluster + summary dot plot)
##   (E) Bulk RV-WGCNA module eigengenes projected onto the 10 CM subtypes
##   (F) MitoCarta3.0 module score by CM subcluster × disease group
##   (G) CM_HMGCS2: HMGCS2 violin + metabolic gene-set scores
##       (ketogenesis / FAO / glycolysis / OXPHOS) — anchors the metabolic
##       loss axis
##   (H) CM_HAND2 + EDNRA dual dot plot — RV-specifying TF + endothelin-1
##       receptor (clinically targeted by ERAs) across CM subtype × group
##   (I) SpaGE-imputed sn_subtype labels projected onto Xenium spatial
##       (loads cached predictions from _xenium_sn_subtype_transfer.R)
##   (J) NPPA / NPPB spatial visualization on the clean Xenium CM object —
##       anchors CM_NPP cluster spatially
##   (supp) k=20 kNN neighborhood enrichment heatmap of SpaGE-imputed CM
##          subtypes — quantifies same-subtype clustering and cross-subtype
##          co-localization, stratified by NF / pRV / RVF
##   (supp) per-subtype spatial overlays — for each of the 10 CM_* subtypes,
##          a 9-FOV composite with the focal subtype in red and all others
##          in gray (one PDF per subtype)
##   (supp) Xenium-projected CM subtype proportions across NF/pRV/RVF —
##          analog of S2C using override-based assignments on Xenium cells
##          (patient-level proportions; Kruskal–Wallis per subtype)
##
## Source:
##   - cm_subclust_new_new.rds (snRNA CM subset, post-merge labeled Subnames)
##   - bulk_heart_modules.csv (bulk WGCNA module assignments — for E)
##   - Human.MitoCarta3.0.xls (gene set for F)
##   - Xenium clean object: ./dependencies/shared/Xenium/cm_clean_clean.rds
##     (12 cm_subtype labels; from additional_scripts/XeniumReanalysis.r)
##   - Xenium SpaGE-imputed: ./dependencies/shared/Xenium/xenium_obj_subclustered.rds
##   - SpaGE label cache: ./output/Supplementary_Figure_2/Xenium/Xenium_cm_sn_subtype_labels.rds
##     (produced by _xenium_sn_subtype_transfer.R)
##
## Output:
##   - per-panel PDFs in ./output/Supplementary_Figure_2/ and ./output/Supplementary_Figure_2/Xenium/
##   - composite ./output/Supplementary_Figure_2/v52_figures/SupplementaryFigure_2.pdf via save_figure()
###############################################################################

source('./helper_scripts/_shared_helpers.R')

## Per-figure output directory (introduced for consistent output paths)
V52_FIG_DIR <- './output/Supplementary_Figure_2'
dir.create(V52_FIG_DIR, showWarnings = FALSE, recursive = TRUE)


## Suppress R's default Rplots.pdf in cwd when Rscript hits a plot call
## that's outside an explicit pdf() ... dev.off() envelope.
pdf(NULL)
dir.create(file.path(V52_FIG_DIR, 'Xenium'), showWarnings = FALSE, recursive = TRUE)
COMP_W <- 5
COMP_H <- 5

## Publication geom scales (linewidths, point sizes, dot ranges, bracket widths)
PS <- pub_scales(COMP_W)

library(Seurat)
library(harmony)
library(ggplot2)
library(dplyr)
library(hdWGCNA)
library(ggpubr)

source('./helper_scripts/spatial_functions.R')

# SCTransform / FindTransferAnchors ship large objects to future workers via
# globals; default 500 MiB ceiling is too small for the full atlas. Bump to 16 GiB.
options(future.globals.maxSize = 16 * 1024^3)

## ── CM subtype rename helper ─────────────────────────────────────────────────
## Applied after every load of cm_new_subclust.rds; does not modify the saved
## RDS. Underscored names are preserved by writing directly to @active.ident
## (Seurat's Idents<- silently rewrites underscores to dashes).
.cm_subname_map <- c(
  Cm1  = 'CM_Stressed',
  Cm2  = 'CM_MYH6',
  Cm3  = 'CM_Lipogenic',
  Cm4  = 'CM_Ketogenic',
  Cm5  = 'CM_GJA5',
  Cm6  = 'CM_CAMK2D',
  Cm7  = 'CM_LMCD1',
  Cm8  = 'CM_Endocardial',
  Cm9  = 'CM_KCNJ3',
  Cm10 = 'CM_OXPHOS'
)

.rename_cm_subnames <- function(obj) {
  new_levels <- unname(.cm_subname_map)
  new_subs   <- factor(unname(.cm_subname_map[as.character(obj$Subnames)]),
                       levels = new_levels)
  obj$Subnames    <- new_subs
  names(new_subs) <- colnames(obj)
  obj@active.ident <- new_subs   # direct slot assignment
  obj
}

## ── Materialize BPCells assay layers to in-memory dgCMatrix ──────────────────
## The clean Xenium object uses BPCells-backed layers; Seurat v5
## AverageExpression / GetAssayData['gene', ] error with
## "data.use.i == Inf : comparison (==) is possible only for atomic and list types"
## on those layers. Replacing the assay with a dgCMatrix-backed one fixes it.
.materialize_xenium_assay <- function(obj) {
  a  <- DefaultAssay(obj)
  cm <- as(GetAssayData(obj, assay = a, layer = 'counts'), 'dgCMatrix')
  obj[[a]] <- CreateAssayObject(counts = cm)
  obj <- NormalizeData(obj, assay = a, verbose = FALSE)
  obj
}

## ── Figure 3D-style multi-FOV canvas ────────────────────────────────────────
## Tiles per-patient ggplots into NF / pRV / RVF column groups at a single
## physical scale (default 3000 µm/inch, matching Figure 3D and Figure 5F).
## Each tissue is an independent rasterised grob placed via cowplot::draw_plot
## so it stays movable in Illustrator. Mirrors new_scripts/Figure_3.R lines
## ~2099-2218 but is generic over the colour aesthetic.
##
## Args
##   .df           data.frame with columns: patient, group, x, y (and whatever
##                 column the plot_fn uses for colour). Coords in raw µm
##                 (function handles Y-flip internally if flip_y = TRUE).
##   plot_fn       function(d) -> ggplot for one patient (geom_point + scale +
##                 theme_void + NoLegend; coord_fixed will be appended)
##   legend_grob   pre-built grob from cowplot::get_legend(...) (or NULL)
##   out_path      output PDF
##   ums_per_inch  physical scale, default 3000
##   tissue_bboxes optional named list (by patient) of full-tissue bboxes:
##                 list("1561" = list(x_min,x_max,y_min,y_max), ...). When
##                 provided, each tile is sized to the full tissue extent
##                 (not just the data subset) so cross-tissue scale is
##                 visually comparable in Figure 3D style.
##   flip_y        flip Y so tissues are right-side up (default TRUE)
.tile_3D_canvas <- function(.df, plot_fn, legend_grob, out_path,
                            ums_per_inch = 3000,
                            tissue_bboxes = NULL,
                            flip_y = TRUE) {
  suppressPackageStartupMessages({ library(cowplot); library(ggrastr) })

  if (flip_y) {
    .df$y <- -.df$y
    if (!is.null(tissue_bboxes)) {
      tissue_bboxes <- lapply(tissue_bboxes, function(bb) {
        list(x_min = bb$x_min, x_max = bb$x_max,
             y_min = -bb$y_max, y_max = -bb$y_min)
      })
    }
  }

  ## Origin-normalize each tissue. If full-tissue bboxes provided, shift by
  ## bbox origin and use bbox dims (so each tile occupies the full tissue
  ## extent regardless of where the plotted points fall).
  .pt_names <- unique(.df$patient)
  .pt_info  <- list(); .pt_plots <- list()
  for (.pt in .pt_names) {
    d <- .df[.df$patient == .pt, ]
    if (!is.null(tissue_bboxes) && .pt %in% names(tissue_bboxes)) {
      bb   <- tissue_bboxes[[.pt]]
      d$x  <- d$x - bb$x_min
      d$y  <- d$y - bb$y_min
      w_um <- bb$x_max - bb$x_min
      h_um <- bb$y_max - bb$y_min
    } else {
      d$x  <- d$x - min(d$x)
      d$y  <- d$y - min(d$y)
      w_um <- diff(range(d$x))
      h_um <- diff(range(d$y))
    }
    .pt_info[[.pt]] <- list(group = as.character(d$group[1]),
                            w_um = w_um, h_um = h_um)
    ## Append coord_fixed with explicit limits — overrides any coord_fixed
    ## already in plot_fn so the device matches the bbox exactly.
    .pt_plots[[.pt]] <- plot_fn(d) +
      coord_fixed(xlim = c(0, w_um), ylim = c(0, h_um), expand = FALSE)
  }

  ## Per-group layout: handles 1, 2, or 3 patients.
  ## 3 patients: 2 largest side-by-side bottom, smallest centred above.
  .gap_um     <- 500
  .row_gap_um <- 1500
  .grp_gap_um <- 5000
  .grp_levels <- intersect(c('NF','pRV','RVF'),
                           unique(sapply(.pt_info, function(x) x$group)))

  .layout_group <- function(.g) {
    pts <- names(.pt_info)[sapply(.pt_info, function(x) x$group == .g)]
    bbs <- .pt_info[pts]
    pts <- pts[order(-sapply(bbs, function(b) b$w_um * b$h_um))]
    bbs <- bbs[pts]
    if (length(pts) == 1) {
      tiles <- list(list(x = 0, y = 0, w = bbs[[1]]$w_um, h = bbs[[1]]$h_um))
      names(tiles) <- pts
      return(list(tiles = tiles,
                  canvas = list(w = bbs[[1]]$w_um, h = bbs[[1]]$h_um),
                  top_pt = pts[1], bot_h = bbs[[1]]$h_um))
    }
    if (length(pts) == 2) {
      tiles <- list(
        list(x = 0,                       y = 0, w = bbs[[1]]$w_um, h = bbs[[1]]$h_um),
        list(x = bbs[[1]]$w_um + .gap_um, y = 0, w = bbs[[2]]$w_um, h = bbs[[2]]$h_um)
      )
      names(tiles) <- pts
      bot_w <- bbs[[1]]$w_um + .gap_um + bbs[[2]]$w_um
      bot_h <- max(bbs[[1]]$h_um, bbs[[2]]$h_um)
      return(list(tiles = tiles,
                  canvas = list(w = bot_w, h = bot_h),
                  top_pt = pts[1], bot_h = bot_h))
    }
    b1 <- bbs[[1]]; b2 <- bbs[[2]]; b3 <- bbs[[3]]
    bot_w <- b1$w_um + .gap_um + b2$w_um
    bot_h <- max(b1$h_um, b2$h_um)
    tiles <- list(
      list(x = 0,                       y = 0,                   w = b1$w_um, h = b1$h_um),
      list(x = b1$w_um + .gap_um,       y = 0,                   w = b2$w_um, h = b2$h_um),
      list(x = (bot_w - b3$w_um) / 2,   y = bot_h + .row_gap_um, w = b3$w_um, h = b3$h_um)
    )
    names(tiles) <- pts
    list(tiles  = tiles,
         canvas = list(w = bot_w, h = bot_h + .row_gap_um + b3$h_um),
         top_pt = pts[3], bot_h = bot_h)
  }
  .layouts <- setNames(lapply(.grp_levels, .layout_group), .grp_levels)

  .canvas_w_um <- max(sapply(.layouts, function(l) l$canvas$w))
  .canvas_h_um <- sum(sapply(.layouts, function(l) l$canvas$h)) +
                  (length(.grp_levels) - 1) * .grp_gap_um

  ## Group Y origins (RVF at y=0, NF on top).
  .grp_y <- numeric(length(.grp_levels)); names(.grp_y) <- .grp_levels
  .cur_y <- 0
  for (.g in rev(.grp_levels)) {
    .grp_y[.g] <- .cur_y
    .cur_y     <- .cur_y + .layouts[[.g]]$canvas$h + .grp_gap_um
  }

  ## Build composite.
  .canvas <- ggdraw()
  for (.g in .grp_levels) {
    lay    <- .layouts[[.g]]
    grp_yo <- .grp_y[.g]
    grp_xo <- (.canvas_w_um - lay$canvas$w) / 2
    for (.pt in names(lay$tiles)) {
      t <- lay$tiles[[.pt]]
      .canvas <- .canvas + draw_plot(
        .pt_plots[[.pt]],
        x      = (grp_xo + t$x)        / .canvas_w_um,
        y      = (grp_yo + t$y)        / .canvas_h_um,
        width  = t$w                   / .canvas_w_um,
        height = t$h                   / .canvas_h_um
      )
    }
    .canvas <- .canvas + draw_label(
      .g, x = 0.01,
      y = (grp_yo + lay$canvas$h) / .canvas_h_um,
      hjust = 0, vjust = 1, fontface = 'bold', size = 22
    )
  }

  ## Legend in RVF top-right whitespace (only available with 3-patient layout).
  if (!is.null(legend_grob) && 'RVF' %in% .grp_levels) {
    .rvf_lay <- .layouts[['RVF']]
    if (length(.rvf_lay$tiles) >= 3) {
      .rvf_yo  <- .grp_y[['RVF']]
      .rvf_xo  <- (.canvas_w_um - .rvf_lay$canvas$w) / 2
      .top_t   <- .rvf_lay$tiles[[.rvf_lay$top_pt]]
      .leg_x_um <- .rvf_xo + .top_t$x + .top_t$w + 1000
      .leg_y_um <- .rvf_yo + .rvf_lay$bot_h + .row_gap_um
      .leg_w_um <- (.rvf_xo + .rvf_lay$canvas$w) - .leg_x_um
      .leg_h_um <- .top_t$h * 0.70
      if (.leg_w_um > 1000 && .leg_h_um > 1000) {
        .canvas <- .canvas + draw_plot(
          legend_grob,
          x      = .leg_x_um / .canvas_w_um,
          y      = (.leg_y_um + .top_t$h - .leg_h_um) / .canvas_h_um,
          width  = .leg_w_um / .canvas_w_um,
          height = .leg_h_um / .canvas_h_um
        )
      }
    }
  }

  ## 2-mm scale bar in bottom-left.
  .sb_x0 <- 0.02
  .sb_x1 <- .sb_x0 + (2000 / .canvas_w_um)
  .sb_y  <- 0.02
  .canvas <- .canvas +
    draw_line(x = c(.sb_x0, .sb_x1), y = c(.sb_y, .sb_y),
              size = 0.8, color = 'black') +
    draw_label('2 mm',
               x = (.sb_x0 + .sb_x1) / 2,
               y = .sb_y - 0.012,
               hjust = 0.5, vjust = 1, size = 12)

  .fig_w_in <- .canvas_w_um / ums_per_inch
  .fig_h_in <- .canvas_h_um / ums_per_inch
  ggsave(out_path, plot = .canvas,
         width = .fig_w_in, height = .fig_h_in,
         device = cairo_pdf, limitsize = FALSE)
  cat('  wrote ', out_path, ' (', round(.fig_w_in, 2),
      ' x ', round(.fig_h_in, 2), ' in)\n', sep = '')
  invisible(NULL)
}

## ── Full-tissue bbox cache (for Figure 3D-style consistent scale) ──────────
## Loads the SpaGE-imputed full Xenium ONCE to compute per-patient tissue
## bboxes (xmin/xmax/ymin/ymax in µm). After that the small RDS is reused.
## Without this, Section I/J tiles use bbox-of-CMs which under-counts tissue
## extent for tissues where CMs cluster in one region.
.xenium_bbox_cache <- './dependencies/shared/Xenium/Xenium_tissue_bboxes.rds'

.get_xenium_tissue_bboxes <- function() {
  if (file.exists(.xenium_bbox_cache)) {
    return(readRDS(.xenium_bbox_cache))
  }
  cat('Building Xenium tissue-bbox cache from full SpaGE-imputed object...\n')
  .xen_full_path <- path.expand(
    './dependencies/shared/Xenium/xenium_obj_subclustered.rds')
  if (!file.exists(.xen_full_path)) {
    cat('  full Xenium not found at', .xen_full_path,
        '— bbox cache unavailable\n')
    return(NULL)
  }
  .xen_full <- readRDS(.xen_full_path)
  .xen_full <- UpdateSeuratObject(.xen_full)
  .md <- .xen_full@meta.data
  .pt_col <- intersect(c('patient','sample','orig.ident'), colnames(.md))[1]
  .x_col  <- intersect(c('x_centroid','x'), colnames(.md))[1]
  .y_col  <- intersect(c('y_centroid','y'), colnames(.md))[1]
  if (is.na(.x_col) || is.na(.y_col)) {
    ## Centroids not in metadata — pull from FOV objects.
    .imgs <- Images(.xen_full)
    .all_coords <- do.call(rbind, lapply(.imgs, function(im) {
      cc <- tryCatch(GetTissueCoordinates(.xen_full, which = 'centroids',
                                          image = im),
                     error = function(e) GetTissueCoordinates(.xen_full[[im]]))
      if (!('cell' %in% colnames(cc))) cc$cell <- rownames(cc)
      cc
    }))
    .xc <- intersect(c('x','x_centroid'), colnames(.all_coords))[1]
    .yc <- intersect(c('y','y_centroid'), colnames(.all_coords))[1]
    .md_coords <- data.frame(cell = .all_coords[['cell']],
                             x = .all_coords[[.xc]],
                             y = .all_coords[[.yc]])
    .md_coords$patient <- .md[.md_coords$cell, .pt_col]
  } else {
    .md_coords <- data.frame(cell = rownames(.md),
                             x = .md[[.x_col]],
                             y = .md[[.y_col]],
                             patient = as.character(.md[[.pt_col]]))
  }
  rm(.xen_full); gc(verbose = FALSE)

  .bb_df <- aggregate(cbind(x, y) ~ patient, data = .md_coords,
                      FUN = function(v) c(min(v), max(v)))
  .bboxes <- list()
  for (i in seq_len(nrow(.bb_df))) {
    .bboxes[[as.character(.bb_df$patient[i])]] <- list(
      x_min = .bb_df$x[i, 1], x_max = .bb_df$x[i, 2],
      y_min = .bb_df$y[i, 1], y_max = .bb_df$y[i, 2]
    )
  }
  saveRDS(.bboxes, .xenium_bbox_cache)
  cat('  wrote', .xenium_bbox_cache, '(', length(.bboxes), 'patients)\n')
  .bboxes
}

##############################################
##############################################
####### Pre-A: Re-cluster CMs (patient-balanced, ~10 final clusters)
##############################################
##############################################
# Goal: build a CM clustering where no cluster is dominated by a single
# patient. Heavy compute (~30-60 min, 5 GB RDS load + SCTransform). Cached
# via file.exists() so it only runs once. Output: cm_subclust_new_new.rds.
#
# Approach:
#   1. Load full snRNA atlas, subset to CMs.
#   2. Aggressive Harmony (theta=4) to push cross-patient integration.
#   3. Pass 1: high-res cluster (res=1.0), drop any cluster >50% from a
#      single donor (donor-driven artifacts).
#   4. Pass 2: re-cluster the cleaned cells at high res.
#   5. Hierarchical merge of cluster centroids → exactly 10 final Cm groups.

.cm_recluster_path  <- './output/Supplementary_Figure_2/cm_subclust_new_new.rds'
.cm_stage2_path     <- './output/Supplementary_Figure_2/cm_subclust_stage2.rds'
.cm_high_res_path   <- './output/Supplementary_Figure_2/cm_subclust_high_res.rds'

## Skip ALL 3 regeneration stages if the canonical cm_subclust_new_new.rds is
## already provided in dependencies/shared/ (the common case). Stages 1+2 are
## intermediates; we only need them to produce stage 3.
.cm_canonical_dep <- './dependencies/shared/cm_subclust_new_new.rds'
.have_canonical_cm <- file.exists(.cm_canonical_dep)
if (.have_canonical_cm) {
  message('Canonical cm_subclust_new_new.rds found in dependencies — skipping CM rebuild stages 1-3.')
}
.cm_high_res_marks  <- './output/Supplementary_Figure_2/cm_recluster_high_res_markers.csv'

# Manual cluster identity assignment (from marker review of cm_recluster_cleaned_markers.csv).
# After merging cluster 1 → 0 (weak markers, transcriptionally closest to 0):
#   0 (was 0+1) → CM_CORIN          atrial-secretory baseline (CORIN/PCSK6/NCALD)
#   1 (was 2)   → CM_BetaMHC        fetal/failure program (MYH7/TNNT1/ACTA1)
#   2 (was 3)   → CM_XIRP           Z-disc remodeling (XIRP1/XIRP2/ANKRD1/CNN1)
#   3 (was 4)   → CM_RORA           circadian/synaptic-like (RORA/SPOCK1)
#   4 (was 5)   → CM_HAND2          RV-canonical (HAND2/EDNRA/TENM3)
#   5 (was 6)   → CM_HTR4           atrial-specialized (HTR4/CACNA1E)
#   6 (was 7)   → CM_HMGCS2         circadian-metabolic (HMGCS2/ARNTL/SOX5)
#   7 (was 8)   → CM_KCNJ3          atrial conduction (KCNJ3/HS6ST3/NEB)
#   8 (was 9)   → CM_MYH6           active baseline (MYH6/STAT3/CAMK1D)
#   9 (was 10)  → CM_NPP            stress/endocardial (NPPA/NPPB/NRXN3)
.cm_label_map <- c(
  '0' = 'CM_Baseline',
  '1' = 'CM_BetaMHC',
  '2' = 'CM_XIRP',
  '3' = 'CM_RORA',
  '4' = 'CM_HAND2',
  '5' = 'CM_HTR4',
  '6' = 'CM_HMGCS2',
  '7' = 'CM_KCNJ3',
  '8' = 'CM_MYH6',
  '9' = 'CM_NPP'
)

# Doublet cluster IDs (high-res cluster numbers from FindAllMarkers output).
# Populate this AFTER inspecting cm_recluster_high_res_markers.csv: clusters
# with mixed-lineage markers (e.g., CM + EC, CM + FB, CM + immune) are likely
# doublets and should be dropped before the hierarchical merge.
.doublet_cluster_ids <- c('10','11','18','19','20')   # 10=fibroblast/Schwann, 11=endothelial (VWF/PECAM1), 18=macrophage (CD163/MRC1/F13A1), 19=pericyte/SMC (RGS5/PDGFRB/CARMN), 20=pericyte/SMC (CARMN/LDB2/PLA2G5/EBF1)

# ---------- Stage 1: heavy SCT / PCA / Harmony / high-res clustering -------
if (!.have_canonical_cm && !file.exists(.cm_high_res_path)) {
  cat('=== Stage 1: building', .cm_high_res_path, '(heavy: 30-60 min) ===\n')

  RV_data <- readRDS('./dependencies/shared/RV_data.rds')
  RV_data <- subset(RV_data, Names == 'CM')
  cat('  CM cells:', ncol(RV_data), '\n')

  # Drop sex-chromosome / lncRNA / mito / metallothionein genes that bias
  # integration or drive clustering on metabolic / stress state rather than
  # CM identity.
  feats <- setdiff(rownames(RV_data), c('XIST','TTTY14','TTTY10','UTY'))
  feats <- setdiff(feats, grep('^LINC', feats, value = TRUE))
  # Mitochondrial-encoded (MT-*), mito-rRNA pseudogenes (MTRNR*),
  # and metallothionein family (MT1*/MT2*/MT3/MT4):
  feats <- setdiff(feats, grep('^MT[-0-9]|^MTRNR', feats, value = TRUE))
  feats <- grep('-',  feats, invert = TRUE, value = TRUE)
  feats <- grep('\\.', feats, invert = TRUE, value = TRUE)
  RV_data <- RV_data[feats, ]

  RV_data <- SCTransform(RV_data, vst.flavor = 'v2',
                         vars.to.regress = c('nFeature_RNA','percent.mt'),
                         verbose = FALSE)
  RV_data <- RunPCA(RV_data, npcs = 30, verbose = FALSE)
  RV_data <- RunHarmony(RV_data, 'patient', theta = 4, max_iter = 30,
                        verbose = FALSE)
  RV_data <- RunUMAP(RV_data, reduction = 'harmony', dims = 1:30, verbose = FALSE)
  RV_data <- FindNeighbors(RV_data, reduction = 'harmony', dims = 1:30,
                           verbose = FALSE)

  ## High-resolution clustering. Donor-dominated micro-clusters are NOT
  ## dropped; the hierarchical merge below will fold them into transcriptionally
  ## similar parents.
  RV_data <- FindClusters(RV_data, resolution = 1.0, verbose = FALSE)
  pass1_share <- sapply(levels(Idents(RV_data)), function(cl) {
    cells <- WhichCells(RV_data, idents = cl)
    pt_counts <- table(RV_data$patient[cells])
    max(pt_counts) / sum(pt_counts)
  })
  cat('High-res clustering — top-patient share per cluster:\n')
  print(round(sort(pass1_share, decreasing = TRUE), 2))

  ## Per-cluster markers — review for doublet identification.
  cat('Running FindAllMarkers on high-res clusters...\n')
  high_res_markers <- FindAllMarkers(RV_data, only.pos = TRUE,
                                     min.pct = 0.25, logfc.threshold = 0.25,
                                     verbose = FALSE)
  dir.create(dirname(.cm_high_res_marks), showWarnings = FALSE, recursive = TRUE)
  write.csv(high_res_markers, .cm_high_res_marks, row.names = FALSE)
  cat('  wrote', .cm_high_res_marks, '\n')

  saveRDS(RV_data, .cm_high_res_path)
  cat('Saved', .cm_high_res_path, '\n')
  cat('  → Inspect', .cm_high_res_marks, 'for doublet clusters\n')
  cat('  → Populate .doublet_cluster_ids near the top of S2 and re-run\n\n')
}

# ---------- Stage 2: drop doublets, re-integrate, FindAllMarkers ----------
# No hclust merging — final cluster unification is done manually in Stage 3.
if (!.have_canonical_cm && !file.exists(.cm_stage2_path)) {
  cat('=== Stage 2: re-integrate cleaned cells ===\n')
  RV_data <- readRDS(.cm_high_res_path)

  if (length(.doublet_cluster_ids) > 0) {
    cat('Dropping doublet clusters (per manual review):',
        paste(.doublet_cluster_ids, collapse = ','), '\n')
    RV_data <- subset(RV_data,
                      idents = setdiff(levels(Idents(RV_data)),
                                       .doublet_cluster_ids))
    cat('  retained cells:', ncol(RV_data), '\n')
  }

  # Re-run the full integration pipeline on the cleaned cells. SCT residuals
  # were computed on the original (incl. doublet) cell set, so recompute.
  DefaultAssay(RV_data) <- 'decontXcounts'
  RV_data <- SCTransform(RV_data, vst.flavor = 'v2',
                         vars.to.regress = c('nFeature_RNA','percent.mt'),
                         verbose = FALSE, assay = 'decontXcounts')
  RV_data <- RunPCA(RV_data, npcs = 30, verbose = FALSE)
  RV_data <- RunHarmony(RV_data, 'patient', theta = 4, max_iter = 30,
                        verbose = FALSE)
  RV_data <- RunUMAP(RV_data, reduction = 'harmony', dims = 1:30, verbose = FALSE)
  RV_data <- FindNeighbors(RV_data, reduction = 'harmony', dims = 1:30,
                           verbose = FALSE)
  RV_data <- FindClusters(RV_data, resolution = 0.5, verbose = FALSE)

  pt_share <- sapply(levels(Idents(RV_data)), function(cl) {
    cells <- WhichCells(RV_data, idents = cl)
    pt_counts <- table(RV_data$patient[cells])
    max(pt_counts) / sum(pt_counts)
  })
  cat('Cleaned re-clustering (res=0.5) — top-patient share per cluster:\n')
  print(round(sort(pt_share, decreasing = TRUE), 2))
  cat('Cleaned cluster count:', length(unique(Idents(RV_data))), '\n')

  # FindAllMarkers on the cleaned clusters — for manual merge review.
  cat('Running FindAllMarkers on cleaned clusters...\n')
  cleaned_marks <- FindAllMarkers(RV_data, only.pos = TRUE,
                                  min.pct = 0.25, logfc.threshold = 0.25,
                                  verbose = FALSE)
  cleaned_marks_csv <- './output/Supplementary_Figure_2/cm_recluster_cleaned_markers.csv'
  write.csv(cleaned_marks, cleaned_marks_csv, row.names = FALSE)
  cat('  wrote', cleaned_marks_csv, '\n')

  # Persist the FindClusters labels in $Subnames so Stage 3 has a stable
  # reference column to merge from.
  RV_data$Subnames <- Idents(RV_data)
  saveRDS(RV_data, .cm_stage2_path)
  cat('Saved', .cm_stage2_path, '\n\n')
} else {
  cat('Stage 2 cache exists at', .cm_stage2_path, '— skipping.\n')
}

# ---------- Stage 3: merge cluster 1 → 0, apply biological labels ---------
# Cluster 1 from Stage 2 had only 23 weak markers (max log2FC=0.44) and is
# transcriptionally closest to cluster 0 (Pearson 0.996, lowest Euclidean).
# Merge it in, renumber 0..9, run FindAllMarkers on the final 10 clusters,
# and apply CM_* biological labels to $Subnames.
if (!.have_canonical_cm && !file.exists(.cm_recluster_path)) {
  cat('=== Stage 3: merge 1 → 0, label, FindAllMarkers ===\n')
  RV_data <- readRDS(.cm_stage2_path)

  # Merge cluster 1 into 0
  ids <- as.character(Idents(RV_data))
  ids[ids == '1'] <- '0'
  # Remap remaining IDs (0,2,3,4,5,6,7,8,9,10) → (0,1,2,3,4,5,6,7,8,9)
  remap <- setNames(as.character(0:9),
                    c('0','2','3','4','5','6','7','8','9','10'))
  new_num_ids <- remap[ids]

  # Apply biological labels (numeric IDs → CM_* names) using direct slot
  # assignment to preserve underscores (Seurat's Idents<- would rewrite them).
  bio_subs <- factor(unname(.cm_label_map[new_num_ids]),
                     levels = unname(.cm_label_map))
  names(bio_subs) <- colnames(RV_data)
  RV_data@active.ident <- bio_subs
  RV_data$Subnames    <- bio_subs

  cat('Running FindAllMarkers on labeled 10 clusters...\n')
  merged_marks <- FindAllMarkers(RV_data, only.pos = TRUE,
                                 min.pct = 0.25, logfc.threshold = 0.25,
                                 verbose = FALSE)
  merged_marks_csv <- './output/Supplementary_Figure_2/cm_recluster_merged_markers.csv'
  write.csv(merged_marks, merged_marks_csv, row.names = FALSE)
  cat('  wrote', merged_marks_csv, '\n')

  pt_share <- sapply(levels(Idents(RV_data)), function(cl) {
    cells <- WhichCells(RV_data, idents = cl)
    pt_counts <- table(RV_data$patient[cells])
    max(pt_counts) / sum(pt_counts)
  })
  cat('Final 10 labeled clusters — top-patient share:\n')
  print(round(sort(pt_share, decreasing = TRUE), 2))

  saveRDS(RV_data, .cm_recluster_path)
  cat('Saved', .cm_recluster_path, '\n\n')
} else {
  cat('Stage 3 cache exists at', .cm_recluster_path, '— skipping.\n')
}

##############################################
##############################################
####### Figure S2A inset — CM subclusters coloured by patient (light-weight,
####### reuses cached cm_subclust_new_new.rds; runs in seconds vs the full
####### S2A rerun below).
##############################################
##############################################

## Prefer the canonical dep; fall back to S2's regenerated output if present.
.cm_cache_for_inset <- if (file.exists('./dependencies/shared/cm_subclust_new_new.rds')) './dependencies/shared/cm_subclust_new_new.rds' else './output/Supplementary_Figure_2/cm_subclust_new_new.rds'
if (file.exists(.cm_cache_for_inset)) {
  message('Figure S2A inset: loading cached CM subcluster object')
  .cm_for_inset <- readRDS(.cm_cache_for_inset)
  if ('patient' %in% colnames(.cm_for_inset@meta.data)) {
    p_S2A_patient <- DimPlot(.cm_for_inset, reduction = 'umap',
                             group.by = 'patient',
                             raster = TRUE, raster.dpi = c(400, 400),
                             pt.size = 1) +
      ggtitle(NULL) + NoLegend() +
      theme_void() +
      theme(plot.margin = margin(0, 0, 0, 0))
    save_figure(p_S2A_patient,
                'Figure_S2_panel_A_inset_patient.pdf',
                width = 3.5, height = 3.5)
    message('  wrote Figure_S2_panel_A_inset_patient.pdf')

    p_S2A_patient_lg <- DimPlot(.cm_for_inset, reduction = 'umap',
                                group.by = 'patient',
                                raster = TRUE, raster.dpi = c(400, 400),
                                pt.size = 1) +
      ggtitle('CM subclustering coloured by patient')
    save_figure(p_S2A_patient_lg,
                'Figure_S2_panel_A_inset_patient_with_legend.pdf',
                width = 5.5, height = 4.0)
    message('  wrote Figure_S2_panel_A_inset_patient_with_legend.pdf')
  }
  rm(.cm_for_inset); gc(verbose = FALSE)
} else {
  message('Figure S2A inset: skipping (cm_subclust_new_new.rds not on disk)')
}


##############################################
##############################################
####### Figure S2A
##############################################
##############################################

RV_data <- readRDS('./dependencies/shared/RV_data.rds')
RV_data <- subset(RV_data, Names == 'CM')

RV_data <- SCTransform(RV_data, vst.flavor = 'v2', assay = 'RNA',
                       vars.to.regress = c('nFeature_RNA', 'percent.mt'))
RV_data <- RunPCA(RV_data, npcs = 30)
RV_data <- RunHarmony(RV_data, 'patient')
RV_data <- RunUMAP(RV_data, reduction = 'harmony', dims = 1:30)
RV_data <- FindNeighbors(RV_data, reduction = 'harmony', dims = 1:30)
RV_data <- FindClusters(RV_data, resolution = 0.3)

## Only remove clusters that actually exist (Louvain cluster count varies by run)
.bad_clusters <- intersect(c('4', '8'), levels(Idents(RV_data)))
if (length(.bad_clusters) > 0) {
  RV_data <- subset(RV_data, idents = .bad_clusters, invert = TRUE)
} else {
  message('S2 cleanup pass 1: clusters 4/8 not present in current Louvain output — skipping.')
}

RV_data <- SCTransform(RV_data, vst.flavor = 'v2', assay = 'decontXcounts',
                       vars.to.regress = c('nFeature_RNA', 'percent.mt'))
RV_data <- RunPCA(RV_data, npcs = 30, assay = 'SCT')
RV_data <- RunHarmony(RV_data, 'patient')
RV_data <- FindNeighbors(RV_data, reduction = 'harmony', dims = 1:30)
RV_data <- RunUMAP(RV_data, reduction = 'harmony', dims = 1:30)
RV_data <- FindClusters(RV_data, resolution = 0.3)

.bad_clusters <- intersect(c('8'), levels(Idents(RV_data)))
if (length(.bad_clusters) > 0) {
  RV_data <- subset(RV_data, idents = .bad_clusters, invert = TRUE)
} else {
  message('S2 cleanup pass 2: cluster 8 not present in current Louvain output — skipping.')
}

feats <- setdiff(rownames(RV_data), c('XIST','TTTY14','TTTY10','UTY'))
feats <- setdiff(feats, grep('^LINC', feats, value = TRUE))
feats <- grep('-',  feats, invert = TRUE, value = TRUE)
feats <- grep('\\.', feats, invert = TRUE, value = TRUE)

RV_data <- SCTransform(RV_data[feats,], vst.flavor = 'v2', assay = 'decontXcounts',
                       vars.to.regress = c('nFeature_RNA','percent.mt'))
RV_data <- RunPCA(RV_data, npcs = 30, assay = 'SCT')
RV_data <- RunHarmony(RV_data, 'patient')
RV_data <- FindNeighbors(RV_data, reduction = 'harmony', dims = 1:30)
RV_data <- RunUMAP(RV_data, reduction = 'harmony', dims = 1:30)
RV_data <- FindClusters(RV_data, resolution = 0.3)

RV.cm.marks <- FindAllMarkers(RV_data)
marker.genes <- unique(subset(RV.cm.marks, p_val_adj < 0.05)$gene)

new.cluster.ids <- c('Cm1','Cm2','Cm3','Cm4','Cm5','Cm6','Cm7','Cm8','Cm9','Cm10')
names(new.cluster.ids) <- levels(RV_data)
RV_data <- RenameIdents(RV_data, new.cluster.ids)

RV_data$Subnames <- RV_data@active.ident
RV_data$SubNames_Groups <- paste(RV_data$Subnames, RV_data$group, sep = '_')

# Apply the canonical CM_* rename so the UMAP labels match all downstream panels.
RV_data <- .rename_cm_subnames(RV_data)

pdf(paste0('./output/Supplementary_Figure_2/', 'CM_New_snUMAP.pdf'), width = 5, height = 5)
print(PlotEmbedding(RV_data, group.by = 'Subnames', point_size = 1, plot_under = TRUE,
                    plot_theme = umap_theme() + NoLegend(),
                    raster_dpi = 400, raster_scale = 0.5))
dev.off()

pdf(paste0('./output/Supplementary_Figure_2/', 'CM_New_snUMAP_patient.pdf'), width = 5, height = 5)
print(PlotEmbedding(RV_data, group.by = 'patient', point_size = 1, plot_under = TRUE,
                    plot_theme = umap_theme() + NoLegend(),
                    raster_dpi = 400, raster_scale = 0.5))
dev.off()

pdf(paste0('./output/Supplementary_Figure_2/', 'CM_New_snUMAP_group.pdf'), width = 5, height = 5)
print(PlotEmbedding(RV_data, group.by = 'group', point_size = 1, plot_under = TRUE,
                    plot_theme = umap_theme() + NoLegend(),
                    raster_dpi = 400, raster_scale = 0.5))
dev.off()

# saveRDS(RV_data, './dependencies/shared/cm_new_subclust.rds')

##############################################
##############################################
####### Figure S2B
##############################################
##############################################

RV_data <- readRDS('./dependencies/shared/cm_subclust_new_new.rds')

# Markers in cm_label_map order (CM_Baseline..CM_NPP). Reduced to one per
# subtype where the second marker was non-specific (per review).
cm_markers <- c(
  'CORIN',  'PCSK6',     # CM_Baseline  — atrial NP convertase + proprotein convertase
  'MYH7',                # CM_BetaMHC   — β-MHC fetal/failure program
  'XIRP1',               # CM_XIRP      — Xin actin-binding repeat (Z-disc/IC-disc)
  'RORA',                # CM_RORA      — circadian TF
  'HAND2',  'EDNRA',     # CM_HAND2     — RV-specifying cardiac TF + endothelin-1 receptor
  'HTR4',                # CM_HTR4      — 5-HT4 receptor
  'HMGCS2',              # CM_HMGCS2    — ketogenesis rate-limiting enzyme
  'KCNJ3',               # CM_KCNJ3     — GIRK1 atrial K+ channel
  'MYH6',                # CM_MYH6      — α-MHC mature CM
  'NPPA',   'NPPB'       # CM_NPP       — paired natriuretic peptides
)

pdf(paste0('./output/Supplementary_Figure_2/', 'CM_new_sn_Dot.pdf'), width = 11, height = 5)
print(DotPlot(RV_data, features = cm_markers,
              group.by = 'Subnames', col.min = 0, col.max = 2) +
      RotatedAxis() +
      ylab('CM subtype') + xlab('Marker gene') +
      theme(panel.border = element_rect(size = 1, fill = NA, color = 'black'),
            axis.line.x = element_blank(),
            axis.line.y = element_blank()))
dev.off()


##############################################
##############################################
####### Figure S2C
##############################################
##############################################

cm.patient <- table(RV_data$Subnames, RV_data$patient)
cm.patient <- t(t(cm.patient) / colSums(cm.patient))

disease  <- c('RVF','pRV','RVF','NF','pRV','pRV','RVF','NF','NF','pRV','NF')
disease  <- c(t(replicate(10, disease)))

cm.patient <- data.frame(disease = disease, cm.patient)
cm.patient$disease <- factor(cm.patient$disease, levels = c('NF','pRV','RVF'))

pdf('./output/Supplementary_Figure_2/CM_clust_counts.pdf', width = 10, height = 3)
print(ggplot(cm.patient, aes(Var1, Freq, color = disease)) +
        geom_boxplot() + theme_classic())
dev.off()

pdf('./output/Supplementary_Figure_2/CM_clust_freq_stats.pdf', width = 12.5, height = 15)
p <- ggboxplot(cm.patient, x = 'disease', y = 'Freq',
               fill = 'disease', group = 'disease') +
  theme_classic() +
  theme(axis.text.x  = element_text(size = 16),
        axis.text.y  = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text  = element_text(size = 16),
        text         = element_text(color = 'black'),
        axis.text    = element_text(color = 'black')) +
  facet_wrap(~Var1, ncol = 5) +
  stat_compare_means(aes(group = disease), method = 'kruskal.test')
print(p)
dev.off()


##############################################
##############################################
####### Figure S2D
##############################################
##############################################

library(enrichR)

M1 <- readRDS(file = './dependencies/shared/cm_subclust_new_new.rds')
a  <- FindAllMarkers(M1)

dbs <- c('GO_Biological_Process_2023','GO_Cellular_Component_2021',
         'GO_Molecular_Function_2021','LINCS_L1000_Chem_Pert_up',
         'LINCS_L1000_Chem_Pert_down','WikiPathway_2021_Human','KEGG_2021_Human')

combined_output <- data.frame()
for (i in 1:length(unique(M1$Subnames))) {
  cur_mod  <- unique(M1$Subnames)[i]
  cur_info <- subset(a, cluster == cur_mod & avg_log2FC > 0 & p_val_adj < 0.05)
  cur_info <- cur_info[, c('gene')]

  enriched <- enrichR::enrichr(cur_info, dbs)
  Sys.sleep(5)
  for (db in names(enriched)) {
    cur_df <- enriched[[db]]
    if (nrow(cur_df) > 1) {
      cur_df$db     <- db
      cur_df$module <- cur_mod
      combined_output <- rbind(combined_output, cur_df)
    }
  }
}

M1 <- SetupForWGCNA(M1, wgcna_name = 'temp')
M1 <- SetEnrichrTable(M1, combined_output)

outdir <- './output/Supplementary_Figure_2/scCM_subclust_enrichr_plot'
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

wrapText <- function(x, len) {
  sapply(x, function(y) paste(strwrap(y, len), collapse = '\n'),
         USE.NAMES = FALSE)
}

enrichr_df <- GetEnrichrTable(M1)

for (i in 1:length(unique(M1$Subnames))) {
  cur_mod   <- unique(M1$Subnames)[i]
  cur_terms <- subset(enrichr_df, module == cur_mod)
  if (nrow(cur_terms) == 0) next
  cur_terms$wrap <- wrapText(cur_terms$Term, 45)
  plot_list <- list()
  for (cur_db in dbs) {
    plot_df <- subset(cur_terms, db == cur_db) %>% top_n(5, wt = Combined.Score)
    plot_df$Combined.Score <- log(plot_df$Combined.Score)
    plot_list[[cur_db]] <- ggplot(plot_df, aes(x = Combined.Score,
                                               y = reorder(wrap, Combined.Score))) +
      geom_bar(stat = 'identity', position = 'identity', color = 'white', fill = 'green') +
      geom_text(aes(label = wrap), x = 0.2, color = 'black', size = 3.5, hjust = 'left') +
      ylab('Term') + xlab('Enrichment log(combined score)') + ggtitle(cur_db) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.title     = element_blank(),
            axis.ticks.y     = element_blank(),
            axis.text.y      = element_blank(),
            plot.title       = element_text(hjust = 0.5))
  }
  pdf(paste0(outdir, '/', cur_mod, '.pdf'), width = 5, height = 4)
  for (plot in plot_list) print(plot)
  dev.off()
}

# Summary GO BP dot plot
selected_terms <- subset(combined_output, db == 'GO_Biological_Process_2023')
selected_terms <- subset(selected_terms, Adjusted.P.value < 0.05)

idx_top_1 <- match(unique(selected_terms$module), selected_terms$module)
idx_top_3 <- sort(c(idx_top_1, idx_top_1 + 1, idx_top_1 + 2))
selected_terms <- selected_terms[idx_top_3, ]

selected_terms$group <- factor(as.character(selected_terms$module),
                               levels = levels(selected_terms$module))
selected_terms$logp  <- pmin(-log(selected_terms$Adjusted.P.value), 10)

library(stringr)
selected_terms$Term <- str_replace(selected_terms$Term, '\\(GO.*', '')
selected_terms <- selected_terms %>% arrange(group)
selected_terms$wrap <- wrapText(selected_terms$Term, 30)
selected_terms$Term <- factor(as.character(selected_terms$Term),
                              levels = rev(unique(as.character(selected_terms$Term))))
selected_terms$wrap <- factor(as.character(selected_terms$wrap),
                              levels = rev(unique(as.character(selected_terms$wrap))))

library(viridis)
p <- selected_terms %>%
  ggplot(aes(x = group, y = wrap, color = logp, size = log(Combined.Score))) +
  geom_point() +
  scale_color_stepsn(colors = rev(magma(256))) +
  RotatedAxis() + xlab('') + ylab('') +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_rect(size = 1, color = 'black', fill = NA),
        axis.line.x  = element_blank(),
        axis.line.y  = element_blank(),
        plot.margin  = margin(0,0,0,0),
        panel.grid   = element_line(size = 0.25, color = 'lightgrey'))

pdf(paste0('./output/Supplementary_Figure_2/', 'CM_by_cluster_New_GO.pdf'), width = 5.25, height = 8)
print(p)
dev.off()


##############################################
##############################################
######## Figure S2E
##############################################
##############################################

seurat_ref <- readRDS(file = './dependencies/shared/cm_subclust_new_new.rds')
seurat_ref <- SetIdent(seurat_ref, value = 'Subnames')

modules_up   <- c('M2','M12','M28')
modules_down <- c('M10','M25','M26')

mapping <- labels2colors(1:100)

consensus_modules <- read.csv('./dependencies/shared/bulk_heart_modules.csv')
consensus_modules <- consensus_modules[, 1:3]
consensus_modules <- subset(consensus_modules, gene_name %in% rownames(seurat_ref))
consensus_modules <- consensus_modules[
  match(unique(consensus_modules$gene_name), consensus_modules$gene_name), ]

score_calc    <- consensus_modules %>% group_by(module) %>% group_split()
module_colors <- unique(unlist(lapply(score_calc, '[[', 'module')))
module_colors <- paste0('M', match(module_colors, mapping))

seurat_ref <- AddModuleScore(seurat_ref, lapply(score_calc, '[[', 'gene_name'),
                             name = 'module_score')

cols_current <- colnames(seurat_ref@meta.data)
cols_current[startsWith(colnames(seurat_ref@meta.data), 'module_score')] <-
  paste0('module_', module_colors)
colnames(seurat_ref@meta.data) <- cols_current

modules_all <- c('M2','M12','M10','M25','M26','M28')

pdf(paste0('./output/Supplementary_Figure_2/', 'CM_cluster_dot_new_subclust_mods.pdf'),
    width = 6, height = 3.125)
p <- DotPlot(seurat_ref, paste0('module_', modules_all),
             group.by = 'Subnames', dot.min = 0, col.min = 0) +
  RotatedAxis() + ylab('') + xlab('') +
  scale_color_gradient2(high = 'red', mid = 'grey95', low = 'blue') +
  coord_flip() +
  theme(panel.border = element_rect(size = 1, fill = NA, color = 'black'),
        axis.line.x  = element_blank(),
        axis.line.y  = element_blank())
print(p)
dev.off()


##############################################
##############################################
####### Figure S2G — MitoCarta3.0 module score by CM subcluster × disease group
####### (was labeled S2F in v52/v54; v57 caption maps this to panel G)
##############################################
##############################################
# X = CM subtype ($Subnames), Y = group (NF/pRV/RVF). Color = z(mean MitoCarta
# module score), size = % cells with score > 0. Mirrors Figure_5.R Panel 5F.

suppressPackageStartupMessages(library(readxl))

RV_data <- readRDS('./dependencies/shared/cm_subclust_new_new.rds')

.mc_path <- './dependencies/shared/Human.MitoCarta3.0.xls'
mito_genes_all <- tryCatch({
  sheet <- readxl::read_excel(.mc_path, sheet = 'A Human MitoCarta3.0')
  sym_col <- intersect(c('Symbol','Gene Symbol','gene_symbol','Gene'),
                       colnames(sheet))[1]
  unique(toupper(as.character(sheet[[sym_col]])))
}, error = function(e) {
  message('  MitoCarta load failed (', e$message,
          '); falling back to ETC/mitoribosome regex.')
  grep('^(NDUF|COX[0-9]|SDH|UQCR|ATP5|MRPS|MRPL|COA|TIMM|TOMM)',
       rownames(RV_data), value = TRUE)
})

mito_inscope <- intersect(mito_genes_all, toupper(rownames(RV_data)))
cat('Section F (MitoCarta):', length(mito_inscope), '/', length(mito_genes_all),
    'MitoCarta3.0 genes present in CM object.\n')

RV_data <- AddModuleScore(RV_data, features = list(mito_inscope),
                          name = 'mitoAll_', assay = 'SCT')

md <- RV_data@meta.data[, c('Subnames', 'group', 'mitoAll_1')]
md$group    <- factor(md$group,
                      levels = intersect(c('NF','pRV','RVF'),
                                         unique(as.character(md$group))))
md$Subnames <- factor(md$Subnames, levels = unname(.cm_label_map))

agg_mito <- md %>%
  dplyr::group_by(group, Subnames) %>%
  dplyr::summarise(mean_score = mean(mitoAll_1),
                   pct_pos    = 100 * mean(mitoAll_1 > 0),
                   n_cells    = dplyr::n(),
                   .groups    = 'drop') %>%
  dplyr::mutate(z = as.numeric(scale(mean_score)))

p_F <- ggplot(agg_mito,
              aes(x = Subnames, y = group, colour = z, size = pct_pos)) +
  geom_point() +
  scale_colour_gradient2(low = '#2166AC', mid = 'grey92', high = '#B2182B',
                         midpoint = 0, name = 'z(mean score)') +
  scale_size_area(max_size = 6, name = '% cells > 0') +
  labs(x = 'CM subcluster', y = NULL,
       title = 'MitoCarta3.0 module score (all genes) by CM subcluster × group') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

pdf('./output/Supplementary_Figure_2/CM_MitoCarta_dotplot.pdf', width = 8, height = 3)
print(p_F)
dev.off()


##############################################
##############################################
####### Figure S2G — CM_HMGCS2 / metabolic gene-set scores
##############################################
##############################################
# HMGCS2 violin + ketogenesis / FAO / glycolysis / OXPHOS gene-set scores
# across the 10 snRNA CM subtypes. Anchors the metabolic-loss axis.

RV_data <- readRDS('./dependencies/shared/cm_subclust_new_new.rds')

ketogenesis_genes <- c('HMGCS2','BDH1','BDH2','OXCT1','OXCT2','ACAT1','HMGCL','ACAT2','ACSS2')
fao_genes         <- c('CPT1A','CPT1B','CPT2','ACADM','ACADL','ACADVL','ACADS','HADHA','HADHB','ECH1','ECHS1','ACAA2','ACSL1','ACSL3','SLC25A20','PPARA')
glycolysis_genes  <- c('HK1','HK2','GPI','PFKM','PFKL','ALDOA','TPI1','GAPDH','PGK1','PGM1','ENO1','ENO3','PKM','LDHA','LDHB','SLC2A1','SLC2A4')
oxphos_genes      <- c('NDUFA1','NDUFA2','NDUFB1','NDUFB7','NDUFS1','NDUFS3','SDHA','SDHB','UQCRC1','UQCRC2','UQCRH','COX4I1','COX5A','COX6A2','COX7A1','ATP5F1A','ATP5F1B','ATP5F1C','ATP5MC1')

ketogenesis_genes <- ketogenesis_genes[ketogenesis_genes %in% rownames(RV_data)]
fao_genes         <- fao_genes[fao_genes                 %in% rownames(RV_data)]
glycolysis_genes  <- glycolysis_genes[glycolysis_genes   %in% rownames(RV_data)]
oxphos_genes      <- oxphos_genes[oxphos_genes           %in% rownames(RV_data)]

RV_data <- AddModuleScore(RV_data, features = list(ketogenesis_genes), name = 'Ketogenesis', assay = 'SCT')
RV_data <- AddModuleScore(RV_data, features = list(fao_genes),         name = 'FAO',          assay = 'SCT')
RV_data <- AddModuleScore(RV_data, features = list(glycolysis_genes),  name = 'Glycolysis',   assay = 'SCT')
RV_data <- AddModuleScore(RV_data, features = list(oxphos_genes),      name = 'OXPHOS',       assay = 'SCT')

pdf('./output/Supplementary_Figure_2/CM_HMGCS2_violin.pdf', width = 6, height = 3.5)
print(VlnPlot(RV_data, features = 'HMGCS2', group.by = 'Subnames', pt.size = 0) + NoLegend())
dev.off()

pdf('./output/Supplementary_Figure_2/CM_HMGCS2_featureplot.pdf', width = 5, height = 5)
print(FeaturePlot(RV_data, features = 'HMGCS2', order = TRUE) + NoAxes())
dev.off()

pdf('./output/Supplementary_Figure_2/CM_metabolic_dot.pdf', width = 6.5, height = 3)
print(DotPlot(RV_data, features = c('Ketogenesis1','FAO1','Glycolysis1','OXPHOS1'),
              group.by = 'Subnames', col.min = 0, dot.min = 0) +
      RotatedAxis() + ylab('') + xlab('') +
      scale_color_gradient2(high = 'red', mid = 'grey95', low = 'blue') + coord_flip() +
      theme(panel.border = element_rect(size = 1, fill = NA, color = 'black'),
            axis.line.x  = element_blank(),
            axis.line.y  = element_blank()))
dev.off()


##############################################
##############################################
####### Figure S2H — CM_HAND2 + EDNRA dual dot plot
##############################################
##############################################
# HAND2 = RV-specifying cardiac transcription factor (developmental).
# EDNRA = endothelin-1 receptor (clinically targeted by ERAs in PAH/RVF).
# Both surface in CM_HAND2 and expand with disease — one of the three
# disease-induced stress-reactivation axes.

RV_data <- readRDS('./dependencies/shared/cm_subclust_new_new.rds')

# Composite Subnames × group identity for a 30-row × 2-column dot plot.
RV_data$Subnames_group <- factor(
  paste(as.character(RV_data$Subnames), as.character(RV_data$group), sep = '_'),
  levels = paste(rep(unname(.cm_label_map), each = 3),
                 c('NF','pRV','RVF'), sep = '_')
)

pdf('./output/Supplementary_Figure_2/CM_HAND2_EDNRA_dot.pdf', width = 5, height = 8)
print(DotPlot(RV_data, features = c('HAND2','EDNRA'),
              group.by = 'Subnames_group',
              col.min = 0, col.max = 2, dot.min = 0) +
      RotatedAxis() +
      ylab('CM subtype × disease group') + xlab('Gene') +
      ggtitle('HAND2 / EDNRA expression by CM subcluster × group'))
dev.off()


##############################################
##############################################
## Supplementary: CM sn_subtype spatial tiles (Figure 3D-style composite).
## This block depends on `.df` / `.df_nb` / `.legend_grob` / `sn_cols` which
## are constructed in the NPPA/NPPB section further below — earlier cleanup
## passes left this render block ahead of its inputs. Wrapped in tryCatch so
## the ordering issue can't halt S2's lettered panels (A-J already written).
tryCatch({
  ## Figure 3D-style composite (NF/pRV/RVF columns, single physical scale,
  ## tiles sized to FULL tissue bboxes so cross-tissue scale is comparable).
  .tissue_bboxes <- .get_xenium_tissue_bboxes()
  .plot_fn <- function(d) {
    d <- d[order(d$.zorder), ]
    ggplot(d, aes(x = x, y = y, color = sn_subtype)) +
      ggrastr::rasterise(geom_point(size = 0.7, stroke = 0), dpi = 400) +
      scale_color_manual(values = sn_cols, name = NULL, drop = FALSE) +
      theme_void() +
      theme(legend.position = 'none', plot.margin = margin(0, 0, 0, 0))
  }
  .tile_3D_canvas(.df, .plot_fn, .legend_grob,
                  './output/Supplementary_Figure_2/Xenium/CM_sn_subtype_spatial_tiled.pdf',
                  ums_per_inch = 3000,
                  tissue_bboxes = .tissue_bboxes,
                  flip_y = TRUE)

  ## Non-baseline composite (drops Baseline pool so disease-axis pops).
  .df_nb <- .df[.df$sn_subtype != 'CM_Baseline', ]
  if (nrow(.df_nb) > 0) {
    .tile_3D_canvas(.df_nb, .plot_fn, .legend_grob,
                    './output/Supplementary_Figure_2/Xenium/CM_sn_subtype_spatial_tiled_nonbase.pdf',
                    ums_per_inch = 3000,
                    tissue_bboxes = .tissue_bboxes,
                    flip_y = TRUE)
  }

  cat('  per-patient sn_subtype tiles in', .tile_dir, '\n')
}, error = function(e) {
  message('Supplementary CM sn_subtype spatial tiles skipped: ', conditionMessage(e))
})


##############################################
##############################################
####### Figure S2I — WGA CM cross-sectional area (CM hypertrophy)
##############################################
##############################################
# WGA-stained FFPE human RV sections; CellProfiler pipeline measures per-CM
# AreaShape_MinFeretDiameter (min Feret diameter, proxy for CM short-axis
# size). Linear mixed model (cell-level data, random intercept by HHTB_ID)
# tests Group effect; BH-adjusted pairwise contrasts via emmeans.
# Per-patient median min Feret used for visualisation (n=4 NF / 4 pRV /
# 4 RVF, ~18K CMs total).

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2)
  library(lme4); library(lmerTest); library(emmeans)
})

.wga_path <- './dependencies/shared/wga_minferet_human_rv.csv'
if (file.exists(.wga_path)) {
  message('Figure S2 panel I: building WGA CM hypertrophy plot')
  df_wga <- read.csv(.wga_path, check.names = FALSE)
  ## CellProfiler exports a leading unnamed row-index column; with
  ## check.names=FALSE its name stays "" and dplyr cannot build a data
  ## mask over it ("Can't transform a data frame with NA or '' names").
  .blank <- is.na(names(df_wga)) | names(df_wga) == ''
  if (any(.blank)) df_wga <- df_wga[, !.blank, drop = FALSE]
  df_wga <- df_wga %>%
    filter(!is.na(AreaShape_MinFeretDiameter),
           !is.na(Group), !is.na(HHTB_ID),
           Group %in% c('NF', 'pRV', 'RVF')) %>%
    mutate(Group = factor(Group, levels = c('NF', 'pRV', 'RVF')),
           HHTB_ID = factor(HHTB_ID))

  # Cell-level LMM with random intercept by patient
  .lmm <- lmer(AreaShape_MinFeretDiameter ~ Group + (1 | HHTB_ID),
               data = df_wga, REML = FALSE)
  .emm <- emmeans(.lmm, ~ Group)
  .pw  <- as.data.frame(summary(pairs(.emm, adjust = 'BH')))
  write.csv(.pw,
            './output/Supplementary_Figure_2/Figure_S2_panel_I_WGA_LMM_pairwise_BH.csv',
            row.names = FALSE)

  # Per-patient summary for plot
  df_patient <- df_wga %>%
    group_by(HHTB_ID, Group) %>%
    summarize(median_MinFeret = median(AreaShape_MinFeretDiameter,
                                       na.rm = TRUE),
              n_cells = n(), .groups = 'drop')
  write.csv(df_patient,
            './output/Supplementary_Figure_2/Figure_S2_panel_I_WGA_per_patient_summary.csv',
            row.names = FALSE)

  # Significance brackets from BH-adjusted pairwise LMM contrasts
  .y_max   <- max(df_patient$median_MinFeret, na.rm = TRUE)
  .y_range <- max(diff(range(df_patient$median_MinFeret, na.rm = TRUE)), 1)

  .brk <- data.frame(
    x1 = c(1, 1, 2),
    x2 = c(2, 3, 3),
    contrast = c('NF - pRV', 'NF - RVF', 'pRV - RVF'))
  .brk$p <- .pw$p.value[match(.brk$contrast, .pw$contrast)]
  .brk$label <- ifelse(.brk$p < 0.001, 'p<0.001',
                       paste0('p=', signif(.brk$p, 2)))
  .brk$y <- .y_max + .y_range * c(0.12, 0.26, 0.40)

  p_S2I <- ggplot(df_patient,
                  aes(x = Group, y = median_MinFeret, fill = Group)) +
    geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 0.4, linewidth = 0.7) +
    geom_jitter(width = 0.12, size = 4.8, alpha = 0.9, shape = 21,
                colour = 'black') +
    scale_fill_manual(values = c(NF = '#4DAF4A', pRV = '#FF7F00',
                                  RVF = '#E41A1C')) +
    geom_segment(data = .brk,
                 aes(x = x1, xend = x2, y = y, yend = y),
                 inherit.aes = FALSE) +
    geom_segment(data = .brk,
                 aes(x = x1, xend = x1, y = y, yend = y - 0.02 * .y_range),
                 inherit.aes = FALSE) +
    geom_segment(data = .brk,
                 aes(x = x2, xend = x2, y = y, yend = y - 0.02 * .y_range),
                 inherit.aes = FALSE) +
    geom_text(data = .brk,
              aes(x = (x1 + x2) / 2, y = y + 0.04 * .y_range,
                  label = label),
              inherit.aes = FALSE, size = 6.4) +
    labs(x = NULL,
         y = expression('CM min Feret diameter ('*mu*'m)'),
         title = 'WGA CM cross-sectional area') +
    theme_classic(base_size = 24) +
    theme(legend.position = 'none')

  save_figure(p_S2I, 'Figure_S2_panel_I_WGA_CSA.pdf',
              width = 4.5, height = 4.5)
  message('  wrote Figure_S2_panel_I_WGA_CSA.pdf')
  message('  BH-adjusted pairwise LMM contrasts:')
  print(.pw)
} else {
  warning(paste('WGA min Feret CSV not found at', .wga_path,
                '\u2014 skipping Figure S2 panel I.'))
}


##############################################
##############################################
####### Figure S2J — NPPA / NPPB spatial tiles (per Figure 5 Panel F approach)
##############################################
##############################################
# Per-patient rasterised tiles + composite faceted overview, continuous
# gray→red color scale on the combined NPP score (mean of NPPA + NPPB log-
# norm expression). Mirrors new_scripts/Figure_5.R lines ~1241-1361 but
# emits to S2 paths and reuses the cached CM-only Xenium object.

suppressPackageStartupMessages({
  library(ggrastr); library(cowplot)
})

.xen_cm_cache <- './dependencies/shared/Xenium/Xenium_cm_with_sn_subtype.rds'
.xen_slim_cache <- './dependencies/shared/Xenium/Xenium_cm_plot_data.rds'
.xen_full_path <- path.expand(
  './dependencies/shared/Xenium/xenium_obj_subclustered.rds')

## Prefer slim plot-data cache (~5 MB) when available — avoids loading the
## multi-GB Seurat object when only NPPA/NPPB + coords are needed.
.df <- NULL
if (file.exists(.xen_slim_cache)) {
  cat('Section J: loading slim plot-data cache\n')
  .df <- readRDS(.xen_slim_cache)
  .df$NPP_score <- rowMeans(cbind(.df$NPPA, .df$NPPB), na.rm = TRUE)
} else if (file.exists(.xen_cm_cache)) {
  cat('Section J: loading heavy Xenium CM cache (slim cache missing)\n')
  .xen_cm <- readRDS(.xen_cm_cache)
} else if (file.exists(.xen_full_path)) {
  cat('Section J: building CM-only Xenium from full subclustered RDS\n')
  .xen_full <- readRDS(.xen_full_path)
  .xen_full <- UpdateSeuratObject(.xen_full)
  .lin_col <- intersect(c('cell_type_rctd_doublet','cell_type_rctd_full',
                          'cell_type','Names','Lineage','lineage'),
                        colnames(.xen_full@meta.data))[1]
  .cm_cells <- rownames(.xen_full@meta.data)[
    grepl('^CM$|Cardiom', .xen_full@meta.data[[.lin_col]],
          ignore.case = TRUE)]
  .xen_cm <- subset(.xen_full, cells = .cm_cells)
  rm(.xen_full); gc(verbose = FALSE)
} else {
  cat('Section J: no Xenium object available — skipping\n')
}

if (is.null(.df) && exists('.xen_cm') && !is.null(.xen_cm)) {
  ## Build .df from heavy cache (only path when slim cache absent).
  .imgs <- Images(.xen_cm)
  cat('  FOVs:', paste(.imgs, collapse = ', '), '\n')
  .coords <- do.call(rbind, lapply(.imgs, function(im) {
    cc <- tryCatch(GetTissueCoordinates(.xen_cm, which = 'centroids', image = im),
                   error = function(e) GetTissueCoordinates(.xen_cm[[im]]))
    if (!('cell' %in% colnames(cc))) cc$cell <- rownames(cc)
    cc
  }))
  .x_col <- intersect(c('x','x_centroid'), colnames(.coords))[1]
  .y_col <- intersect(c('y','y_centroid'), colnames(.coords))[1]

  .md     <- .xen_cm@meta.data
  .pt_col <- intersect(c('patient','sample','orig.ident'), colnames(.md))[1]
  .grp_col <- intersect(c('group','disease','disease_state'), colnames(.md))[1]

  .feat_df <- tryCatch(
    Seurat::FetchData(.xen_cm, vars = c('NPPA','NPPB'), layer = 'data'),
    error = function(e) Seurat::FetchData(.xen_cm, vars = c('NPPA','NPPB')))
  .nppa <- if ('NPPA' %in% colnames(.feat_df)) .feat_df[, 'NPPA'] else NA_real_
  .nppb <- if ('NPPB' %in% colnames(.feat_df)) .feat_df[, 'NPPB'] else NA_real_

  .df <- data.frame(
    cell    = colnames(.xen_cm),
    patient = as.character(.md[[.pt_col]]),
    group   = if (!is.na(.grp_col)) as.character(.md[[.grp_col]]) else NA_character_,
    NPPA    = .nppa,
    NPPB    = .nppb,
    stringsAsFactors = FALSE
  )
  .df$NPP_score <- rowMeans(cbind(.df$NPPA, .df$NPPB), na.rm = TRUE)

  .coords_keep <- data.frame(cell = .coords[['cell']],
                              x = .coords[[.x_col]],
                              y = .coords[[.y_col]])
  .df <- merge(.df, .coords_keep, by = 'cell', all.x = FALSE)
  rm(.xen_cm); gc(verbose = FALSE)
}

if (!is.null(.df)) {

  .ums_per_inch <- 3000  # match Fig 3D / Fig 5F physical scale
  .tile_dir <- './output/Supplementary_Figure_2/Xenium/CM_NPP_tiles'
  dir.create(.tile_dir, showWarnings = FALSE, recursive = TRUE)

  ## Per-patient rasterised tiles (Illustrator-friendly, no legend).
  ## Tile bbox = full tissue bbox when available so all per-patient PDFs
  ## share a true cross-tissue physical scale.
  .tile_bboxes <- .get_xenium_tissue_bboxes()
  for (.pt in unique(.df$patient)) {
    .sub <- .df[.df$patient == .pt, ]
    if (all(is.na(.sub$NPP_score))) next
    .grp <- .sub$group[1]
    if (!is.null(.tile_bboxes) && .pt %in% names(.tile_bboxes)) {
      bb  <- .tile_bboxes[[.pt]]
      .xr <- c(bb$x_min, bb$x_max)
      .yr <- c(bb$y_min, bb$y_max)
    } else {
      .xr <- range(.sub$x); .yr <- range(.sub$y)
    }
    .w_in <- diff(.xr) / .ums_per_inch
    .h_in <- diff(.yr) / .ums_per_inch
    .sub2 <- .sub[order(.sub$NPP_score), ]
    .p <- ggplot(.sub2, aes(x = x, y = y, color = NPP_score)) +
      ggrastr::rasterise(geom_point(size = 0.45, stroke = 0), dpi = 400) +
      scale_color_gradient(low = 'grey90', high = 'red', name = 'NPPA/NPPB') +
      coord_fixed(xlim = .xr, ylim = .yr, expand = FALSE) +
      theme_void() +
      theme(legend.position = 'none', plot.margin = margin(0, 0, 0, 0))
    ggsave(file.path(.tile_dir,
                     sprintf('CM_NPP_%s_%s.pdf', .grp, .pt)),
           plot = .p, width = .w_in, height = .h_in, device = cairo_pdf)
  }

  ## Figure 3D-style composite + standalone legend.
  if (!all(is.na(.df$NPP_score))) {
    .npp_max <- max(.df$NPP_score, na.rm = TRUE)
    .scale_lo <- 'grey90'; .scale_hi <- 'red'

    ## Standalone legend (continuous gradient).
    .p_legend_src <- ggplot(.df[1:min(2000, nrow(.df)), ],
                             aes(x = x, y = y, color = NPP_score)) +
      geom_point(size = 0.5) +
      scale_color_gradient(low = .scale_lo, high = .scale_hi,
                           limits = c(0, .npp_max), name = 'NPPA/NPPB') +
      theme(legend.position = 'right')
    .legend_grob <- cowplot::get_legend(.p_legend_src)
    pdf('./output/Supplementary_Figure_2/Xenium/CM_NPP_spatial_legend.pdf', width = 2.5, height = 4)
    grid::grid.newpage()
    grid::grid.draw(.legend_grob)
    dev.off()
    cat('  wrote ./output/Supplementary_Figure_2/Xenium/CM_NPP_spatial_legend.pdf\n')

    ## Per-patient ggplot factory: low-NPP first so high-expressing cells
    ## are drawn on top of the gray pool.
    .plot_fn <- function(d) {
      d <- d[order(d$NPP_score), ]
      ggplot(d, aes(x = x, y = y, color = NPP_score)) +
        ggrastr::rasterise(geom_point(size = 0.7, stroke = 0), dpi = 400) +
        scale_color_gradient(low = .scale_lo, high = .scale_hi,
                             limits = c(0, .npp_max), name = NULL) +
        theme_void() +
        theme(legend.position = 'none', plot.margin = margin(0, 0, 0, 0))
    }
    .tissue_bboxes <- .get_xenium_tissue_bboxes()
    .tile_3D_canvas(.df, .plot_fn, .legend_grob,
                    './output/Supplementary_Figure_2/Xenium/CM_NPP_spatial_tiled.pdf',
                    ums_per_inch = 3000,
                    tissue_bboxes = .tissue_bboxes,
                    flip_y = TRUE)
  }
  cat('  per-patient NPP tiles in', .tile_dir, '\n')
}


##############################################
##############################################
