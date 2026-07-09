#!/usr/bin/env Rscript

# =============================================================================
# RV_Atlas integration note (added on archiving into the repo):
#   Replicability source for Supplementary Data S14 (CellChat spatial niche
#   communication). Author: Yonathan (spatial_cellchat_FINAL pipeline).
#   Paths below are repointed to ./dependencies/shared/Xenium; inputs + BPCells
#   backing are present there. NOT run at supp-data build time — Generate_
#   Supplementary_Data.R ships the pre-computed table
#   (dependencies/shared/Xenium/xenium_cellchat_niche_communication.csv). Rerun
#   only to regenerate: heavy (per-patient x per-niche CellChat, ~hours) and
#   needs CellChat (GitHub) + BPCells (C compiler). Writes to output/spatial_cellchat.
# =============================================================================
# COMBINED CELLCHAT PIPELINE v3
#
# Edwards Lab, CHOP — RV Failure spatial transcriptomics CellChat analysis
# CellChat v3 | R | spatial_cellchat_FINAL
#
# Four sequential stages run in a single script:
#
# STAGE 1 — METADATA ASSEMBLY
#   Loads the resegmented Xenium Seurat object, adds niche assignments from
#   the kmeans CSV, and transfers subclustering metadata from the subclustered
#   object. Enriched object held in memory (not saved to disk).
#
# STAGE 2 — CELLCHAT PER-PATIENT PER-NICHE ANALYSIS
#   Runs CellChat (CellChatDB.human) on each valid patient x niche combination
#   using truncatedMean (trim=0.25) to handle sparse Xenium data. Checkpoint
#   logic skips already-completed combinations. Saves per-combo .rds objects
#   and interaction CSVs to cellchat_per_patient_per_niche/.
#
# STAGE 3 — DIFFERENTIAL ANALYSIS + SUPPLEMENTARY TABLE
#   Wilcoxon tests at cell-type-pair and LR x niche aggregated levels across
#   NF vs pRV and pRV vs RVF comparisons. Outputs: per-comparison Wilcoxon
#   tables, tier tables, master comparison table, and a formatted supplementary
#   table of niche-specific communication patterns for publication.
#   Results saved to cellchat_lr_differential/.
#
# STAGE 4 — FIGURE GENERATION
#   fig1: Top 15 pathways, Myocardium niche, normalized to sum=1 per group
#   fig2: Chord diagram, Myocardium, normalized by cell type proportion
#   fig3:  SEMA3C circle plots, NF | pRV | RVF
#   fig4:  Custom hexbin spatial network plots, SEMA3 pathway, per patient
#          (full + zoomed versions; FOV artifact mitigated by 100um binning)
#   fig5:  SEMA3C quantification: SEMA3C-PLXND1 LR prob, Neural->EC, CM->EC
#   Figures saved to cellchat_lr_differential/figures/.
# =============================================================================

# =============================================================================
# LIBRARIES
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(BPCells)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(pheatmap)
  library(circlize)
  library(CellChat)
  library(future)
})

plan("sequential")
options(future.globals.maxSize = 8 * 1024^3)

# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================

# Repointed for the RV_Atlas repo (originally /mnt/isilon/edwardslab/Yonathan/...).
# All three raw inputs (Xenium_resegmented_imputed_final.rds, xenium_obj_subclustered.rds,
# Xenium_niche_labels.csv) plus the BPCells `bpcells_imputed` backing dir live under
# dependencies/shared/Xenium, so the file.path(RAW_DIR, ...) inputs below resolve as-is.
RAW_DIR     <- '/Users/ikuz/Documents/RV_Atlas/dependencies/shared/Xenium'
# DATA_DIR removed — standard CellChatDB.human used, no custom database files needed
RESULTS_DIR <- '/Users/ikuz/Documents/RV_Atlas/output/spatial_cellchat'
set.seed(1)  # CellChat computeCommunProb permutation p-values + scPalette are otherwise nondeterministic

SEURAT_RAW_PATH        <- file.path(RAW_DIR, 'Xenium_resegmented_imputed_final.rds')
SEURAT_SUBCLUST_PATH   <- file.path(RAW_DIR, 'xenium_obj_subclustered.rds')
NICHE_LABELS_PATH      <- file.path(RAW_DIR, 'Xenium_niche_labels.csv')
# SEURAT_OUT_PATH removed — enriched object used in-memory only, not saved to disk
# CellChatDB: use standard human database (CellChatDB.human)

PP_OUT_DIR    <- file.path(RESULTS_DIR, 'cellchat_per_patient_per_niche')
DIFF_OUT_DIR  <- file.path(RESULTS_DIR, 'cellchat_lr_differential')
FIG_DIR       <- file.path(DIFF_OUT_DIR, 'figures')

patient_groups <- list(
  NF  = c("1561", "1697", "1691"),
  pRV = c("1618", "1567", "1692"),
  RVF = c("1467", "1632", "1343")
)
group_order  <- c("NF", "pRV", "RVF")
group_colors <- c(NF = "#4393C3", pRV = "#F4A582", RVF = "#D6604D")
TARGET_NICHE <- "Myocardium"

# Differential analysis thresholds
HI_THRESH <- 5e-3
LO_THRESH <- 5e-4

# Chord diagram visualization thresholds
VIZ_NF_HI  <- 1e-5
VIZ_RVF_HI <- 1e-3
VIZ_LO     <- 1e-6

# CellChat filtering thresholds
CELLCHAT_TRIM     <- 0.25
CELLCHAT_MIN_CELL <- 5
COMBO_MIN_CELLS   <- 100
COMBO_MIN_CT      <- 10
COMBO_MIN_CT_N    <- 2

# Figure LR pairs and cell subsets
lr_pairs_to_plot <- list(
  list(name = "SEMA3C_PLXND1", label = "SEMA3C - PLXND1"),
  list(name = "PTPRM_PTPRM",   label = "PTPRM - PTPRM")
)
lr_cell_subsets <- list(
  PTPRM_PTPRM   = c("CM", "EC", "FB", "Pericyte", "VSMC"),
  SEMA3C_PLXND1 = c("Neural", "CM", "FB", "EC", "Macrophage", "Pericyte")
)
lr_primary_niche <- list(
  SEMA3C_PLXND1 = "Myocardium",
  PTPRM_PTPRM   = "Fibrotic remodeling stroma"
)
pathway_primary_niche <- list(
  SEMA3 = "Myocardium",
  PTPRM = "Fibrotic remodeling stroma"
)
pathway_receivers <- list(
  SEMA3 = c("EC", "Macrophage", "Pericyte", "LEC", "Myeloid"),
  PTPRM = c("EC", "CM", "Pericyte", "VSMC", "FB")
)

# =============================================================================
# SHARED HELPER FUNCTIONS
# =============================================================================

# ── Cell type grouping ────────────────────────────────────────────────────────

create_cell_type_groups <- function(x) {
  if (grepl("Arterial_EC|Venous_EC|Capillary_EC|Activated_EC|Endocardial", x)) return("EC")
  if (grepl("^CM_", x))                                    return("CM")
  if (grepl("^FB_", x))                                    return("FB")
  if (grepl("^Macrophage_", x))                            return("Macrophage")
  if (grepl("^Pericyte", x))                               return("Pericyte")
  if (grepl("^VSMC", x))                                   return("VSMC")
  if (grepl("^CD4_T|^CD8_T|^NK$", x))                     return("Lymphocyte")
  if (grepl("^Monocyte|^Dendritic_Cell|^Myeloid_Proliferating", x)) return("Myeloid")
  if (grepl("^Adipo_", x))                                 return("Adipo")
  if (grepl("^Epi_", x))                                   return("Epi")
  if (grepl("^Neuron_|^Schwann_", x))                      return("Neural")
  if (x == "LEC")                                          return("LEC")
  if (x == "Mast_Cell")                                    return("Mast")
  return(x)
}

# ── Matrix utilities ──────────────────────────────────────────────────────────

pad_matrix <- function(mat, all_types) {
  missing_r <- setdiff(all_types, rownames(mat))
  if (length(missing_r) > 0)
    mat <- rbind(mat, matrix(0, length(missing_r), ncol(mat),
                             dimnames = list(missing_r, colnames(mat))))
  missing_c <- setdiff(all_types, colnames(mat))
  if (length(missing_c) > 0)
    mat <- cbind(mat, matrix(0, nrow(mat), length(missing_c),
                             dimnames = list(rownames(mat), missing_c)))
  mat[all_types, all_types]
}

make_source_color_matrix <- function(mat, ct_colors) {
  col_mat <- matrix(NA_character_, nrow(mat), ncol(mat), dimnames = dimnames(mat))
  for (src in rownames(mat))
    col_mat[src, ] <- if (src %in% names(ct_colors))
      adjustcolor(ct_colors[src], alpha.f = 0.7) else "#AAAAAA"
  col_mat
}

# ── Safe CellChat array accessors ─────────────────────────────────────────────
# Guard against NULL, wrong dimensions, or missing pathway/LR name.
# These can occur when computeCommunProbPathway() found no signal or when
# only 1-2 cell types were active, collapsing the 3D probability array.

safe_get_netP <- function(obj, pathway) {
  if (is.null(obj@netP$prob))                              return(NULL)
  if (length(dim(obj@netP$prob)) < 3)                      return(NULL)
  if (!pathway %in% dimnames(obj@netP$prob)[[3]])          return(NULL)
  obj@netP$prob[, , pathway]
}

safe_get_net <- function(obj, lr_name) {
  if (is.null(obj@net$prob))                               return(NULL)
  if (length(dim(obj@net$prob)) < 3)                       return(NULL)
  if (!lr_name %in% dimnames(obj@net$prob)[[3]])           return(NULL)
  obj@net$prob[, , lr_name]
}

# ── Object loaders ────────────────────────────────────────────────────────────

load_obj <- function(fname, base_dir = PP_OUT_DIR) {
  tryCatch(readRDS(file.path(base_dir, 'cellchat_objects', fname)),
           error = function(e) NULL)
}

extract_lr_df <- function(obj, patient_id, niche_label) {
  tryCatch({
    df <- subsetCommunication(obj)
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df$patient <- patient_id
    df$niche   <- niche_label
    df
  }, error = function(e) NULL)
}

# ── Matrix averaging across patients/objects ──────────────────────────────────
# Average netP$prob[,,pathway] matrices across all patients in a group
# for a specific niche. Each patient may have different cell types so matrices
# are padded to their union before averaging.

load_niche_pathway_mats <- function(niche_label, pathway, fdf = file_df) {
  nf_df <- fdf[fdf$niche == niche_label, ]
  lapply(setNames(group_order, group_order), function(grp) {
    mats <- list()
    for (fn in nf_df$fname[nf_df$group == grp]) {
      obj <- load_obj(fn); if (is.null(obj)) next
      mat <- safe_get_netP(obj, pathway)
      if (!is.null(mat)) mats[[length(mats) + 1]] <- mat
      rm(obj); gc()
    }
    if (length(mats) == 0) return(NULL)
    all_ct <- sort(unique(unlist(lapply(mats, rownames))))
    Reduce("+", lapply(mats, pad_matrix, all_types = all_ct)) / length(mats)
  })
}

load_niche_lr_mats <- function(niche_label, lr_name, fdf = file_df) {
  nf_df <- fdf[fdf$niche == niche_label, ]
  lapply(setNames(group_order, group_order), function(grp) {
    mats <- list()
    for (fn in nf_df$fname[nf_df$group == grp]) {
      obj <- load_obj(fn); if (is.null(obj)) next
      mat <- safe_get_net(obj, lr_name)
      if (!is.null(mat)) mats[[length(mats) + 1]] <- mat
      rm(obj); gc()
    }
    if (length(mats) == 0) return(NULL)
    all_ct <- sort(unique(unlist(lapply(mats, rownames))))
    Reduce("+", lapply(mats, pad_matrix, all_types = all_ct)) / length(mats)
  })
}

get_group_avg_net <- function(group_label, fdf = file_df) {
  mats <- list()
  for (fn in fdf$fname[fdf$group == group_label]) {
    obj <- load_obj(fn); if (is.null(obj)) next
    if (!is.null(obj@net$weight)) mats[[length(mats) + 1]] <- obj@net$weight
    rm(obj); gc()
  }
  if (length(mats) == 0) return(NULL)
  all_ct <- sort(unique(unlist(lapply(mats, rownames))))
  Reduce("+", lapply(mats, pad_matrix, all_types = all_ct)) / length(mats)
}

# ── Averaged CellChat object builder ─────────────────────────────────────────
# Constructs a representative CellChat object for a group x niche by averaging
# net$count, net$weight, and netP$prob across all patients. The patient with
# the most interactions is used as the structural template. Used by
# netVisual_aggregate for circle/hierarchy plots.

get_avg_cellchat_obj <- function(group_label, niche_label, fdf = file_df) {
  niche_fdf <- fdf[fdf$niche == niche_label & fdf$group == group_label, ]
  if (nrow(niche_fdf) == 0) return(NULL)
  objs <- Filter(Negate(is.null), lapply(niche_fdf$fname, load_obj))
  if (length(objs) == 0) return(NULL)
  n_int    <- sapply(objs, function(o) tryCatch(sum(o@net$count), error = function(e) 0))
  template <- objs[[which.max(n_int)]]
  if (length(objs) == 1) return(template)
  all_ct     <- sort(unique(unlist(lapply(objs, function(o) rownames(o@net$count)))))
  avg_count  <- Reduce("+", lapply(objs, function(o) pad_matrix(o@net$count,  all_ct))) / length(objs)
  avg_weight <- Reduce("+", lapply(objs, function(o) pad_matrix(o@net$weight, all_ct))) / length(objs)
  if (!is.null(template@netP$prob) && length(dim(template@netP$prob)) == 3) {
    netP_arrays <- Filter(Negate(is.null),
                          lapply(objs, function(o) {
                            if (!is.null(o@netP$prob) &&
                                length(dim(o@netP$prob)) == 3) o@netP$prob else NULL
                          }))
    if (length(netP_arrays) > 0) {
      all_pw <- sort(unique(unlist(lapply(netP_arrays, function(a) dimnames(a)[[3]]))))
      pad_array <- function(arr) {
        ct <- dimnames(arr)[[1]]; pw <- dimnames(arr)[[3]]
        out <- array(0, dim = c(length(all_ct), length(all_ct), length(all_pw)),
                     dimnames = list(all_ct, all_ct, all_pw))
        out[ct, ct, pw] <- arr; out
      }
      template@netP$prob     <- Reduce("+", lapply(netP_arrays, pad_array)) / length(netP_arrays)
      template@netP$pathways <- all_pw
    }
  }
  shared_ct           <- intersect(rownames(template@net$count), all_ct)
  template@net$count  <- avg_count[shared_ct,  shared_ct]
  template@net$weight <- avg_weight[shared_ct, shared_ct]
  template
}

# ── Statistical helpers ───────────────────────────────────────────────────────

assign_tier <- function(r, p) {
  tier <- rep("3_exploratory", length(r))
  tier[abs(r) >= 0.67 & p < 0.2] <- "2_large_effect"
  tier[abs(r) == 1    & p < 0.1] <- "1_top"
  tier
}

assign_pattern <- function(prob_g1, prob_g2, group2_label) {
  dplyr::case_when(
    prob_g2 >= HI_THRESH & prob_g1 < LO_THRESH ~ paste0("gained_in_", group2_label),
    prob_g1 >= HI_THRESH & prob_g2 < LO_THRESH ~ paste0("lost_in_",   group2_label),
    TRUE ~ "quantitative_change"
  )
}

safe_wilcox <- function(x, y) {
  if (length(x) == 0 || length(y) == 0) return(list(p = NA_real_, r = NA_real_))
  tryCatch({
    wt <- suppressWarnings(wilcox.test(x, y, exact = FALSE))
    W  <- as.numeric(wt$statistic)
    list(p = wt$p.value, r = 1 - (2 * W) / (length(x) * length(y)))
  }, error = function(e) list(p = NA_real_, r = NA_real_))
}

# ── Chord diagram helpers ─────────────────────────────────────────────────────

# Prepare NF and RVF display matrices with asymmetric thresholds.
# Returns a list with pre-filtered matrices, union cell types, and global max.
prep_chord_mats <- function(avg_mats, nf_hi, rvf_hi, lo) {
  mat_nf  <- avg_mats[["NF"]]
  mat_rvf <- avg_mats[["RVF"]]
  all_ct  <- sort(unique(c(
    if (!is.null(mat_nf))  rownames(mat_nf)  else character(0),
    if (!is.null(mat_rvf)) rownames(mat_rvf) else character(0)
  )))
  if (length(all_ct) == 0) return(NULL)
  m_nf  <- if (!is.null(mat_nf))  pad_matrix(mat_nf,  all_ct) else
    matrix(0, length(all_ct), length(all_ct), dimnames = list(all_ct, all_ct))
  m_rvf <- if (!is.null(mat_rvf)) pad_matrix(mat_rvf, all_ct) else
    matrix(0, length(all_ct), length(all_ct), dimnames = list(all_ct, all_ct))
  nf_active  <- m_nf  >= nf_hi  & !(m_nf < lo & m_rvf < lo)
  rvf_active <- m_rvf >= rvf_hi & !(m_nf < lo & m_rvf < lo)
  m_nf_d  <- m_nf;  m_nf_d[!nf_active]  <- 0
  m_rvf_d <- m_rvf; m_rvf_d[!rvf_active] <- 0
  keep_ct <- all_ct[sapply(all_ct, function(ct)
    any(m_nf_d[ct,] > 0 | m_nf_d[,ct] > 0 | m_rvf_d[ct,] > 0 | m_rvf_d[,ct] > 0))]
  if (length(keep_ct) < 2) return(NULL)
  list(m_nf_d     = m_nf_d[keep_ct, keep_ct],
       m_rvf_d    = m_rvf_d[keep_ct, keep_ct],
       keep_ct    = keep_ct,
       global_max = max(c(m_nf[keep_ct, keep_ct], m_rvf[keep_ct, keep_ct]), na.rm = TRUE))
}

# Draw one chord panel into the active graphics device
draw_chord <- function(mat_display, keep_ct, global_max, grp,
                       label, niche_label, panel_title, min_frac = 0.02) {
  ct_colors <- scPalette(length(keep_ct)) |> setNames(keep_ct)
  mat_plot  <- mat_display
  mat_plot[mat_plot < global_max * min_frac] <- 0
  if (sum(mat_plot) == 0) {
    plot.new()
    title(main = sprintf("%s\n%s | %s\n(no active edges)", panel_title, label, grp),
          col.main = group_colors[grp])
    return(invisible(NULL))
  }
  present <- keep_ct[keep_ct %in% rownames(mat_plot)]
  use_col <- ct_colors[intersect(present, names(ct_colors))]
  col_mat <- make_source_color_matrix(mat_plot[present, present, drop = FALSE], use_col)
  tryCatch({
    circos.clear()
    chordDiagram(mat_plot[present, present, drop = FALSE],
                 grid.col = use_col, col = col_mat, transparency = 0.3,
                 annotationTrack = c("grid", "name"),
                 preAllocateTracks = list(track.height = 0.08),
                 directional = 1, direction.type = "arrows",
                 link.arr.type = "big.arrow", self.link = 1, scale = FALSE)
    circos.clear()
    title(main = sprintf("%s\n%s | %s\nmax prob: %.2e",
                         panel_title, label, grp, max(mat_plot)),
          col.main = group_colors[grp], cex.main = 1.1, line = -1)
  }, error = function(e) {
    circos.clear(); plot.new()
    title(main = sprintf("%s | %s\nError: %s", label, grp, e$message))
    cat(sprintf("      Chord error %s %s: %s\n", label, grp, e$message))
  })
}

# Two-panel NF vs RVF chord saved to a single PDF
plot_nf_rvf_chord <- function(avg_mats, label, niche_label, out_path,
                               ct_colors = NULL, cell_subset = NULL,
                               binary = FALSE, min_frac = 0.02,
                               nf_hi = VIZ_NF_HI, rvf_hi = VIZ_RVF_HI,
                               lo = VIZ_LO) {
  mat_nf  <- avg_mats[["NF"]]
  mat_rvf <- avg_mats[["RVF"]]
  if (!is.null(cell_subset)) {
    subset_mat <- function(m) {
      if (is.null(m)) return(NULL)
      m_pad <- pad_matrix(m, union(rownames(m), cell_subset))
      keep  <- intersect(cell_subset, rownames(m_pad))
      m_pad[keep, keep]
    }
    mat_nf <- subset_mat(mat_nf); mat_rvf <- subset_mat(mat_rvf)
  }
  all_ct <- sort(unique(c(
    if (!is.null(mat_nf))  rownames(mat_nf)  else character(0),
    if (!is.null(mat_rvf)) rownames(mat_rvf) else character(0)
  )))
  if (length(all_ct) == 0) { cat(sprintf("      No cell types for %s\n", label)); return(invisible(NULL)) }
  m_nf  <- if (!is.null(mat_nf))  pad_matrix(mat_nf,  all_ct) else
    matrix(0, length(all_ct), length(all_ct), dimnames = list(all_ct, all_ct))
  m_rvf <- if (!is.null(mat_rvf)) pad_matrix(mat_rvf, all_ct) else
    matrix(0, length(all_ct), length(all_ct), dimnames = list(all_ct, all_ct))
  nf_active  <- m_nf  >= nf_hi  & !(m_nf < lo & m_rvf < lo)
  rvf_active <- m_rvf >= rvf_hi & !(m_nf < lo & m_rvf < lo)
  m_nf_d  <- m_nf;  m_nf_d[!nf_active]  <- 0
  m_rvf_d <- m_rvf; m_rvf_d[!rvf_active] <- 0
  keep_ct <- all_ct[sapply(all_ct, function(ct)
    any(m_nf_d[ct,] > 0 | m_nf_d[,ct] > 0 | m_rvf_d[ct,] > 0 | m_rvf_d[,ct] > 0))]
  if (length(keep_ct) < 2) { cat(sprintf("      < 2 active cell types for %s\n", label)); return(invisible(NULL)) }
  cat(sprintf("      Active: %s\n", paste(keep_ct, collapse = ", ")))
  m_nf_d    <- m_nf_d[keep_ct, keep_ct]
  m_rvf_d   <- m_rvf_d[keep_ct, keep_ct]
  global_max <- max(c(m_nf[keep_ct, keep_ct], m_rvf[keep_ct, keep_ct]), na.rm = TRUE)
  if (global_max == 0) return(invisible(NULL))
  if (is.null(ct_colors)) ct_colors <- scPalette(length(keep_ct)) |> setNames(keep_ct)
  panels <- list(
    list(mat = m_nf_d,  grp = "NF",  hi = nf_hi,
         title = "Protective signaling in NF\n(lost in disease progression)"),
    list(mat = m_rvf_d, grp = "RVF", hi = rvf_hi,
         title = "Maladaptive signaling in RVF\n(gained in disease progression)")
  )
  pdf(out_path, width = 12, height = 6)
  par(mfrow = c(1, 2), mar = c(1, 1, 5, 1))
  for (panel in panels) {
    mat_plot <- if (binary) (panel$mat >= panel$hi) * 1.0 else {
      tmp <- panel$mat; tmp[tmp < global_max * min_frac] <- 0; tmp
    }
    if (sum(mat_plot) == 0) {
      plot.new()
      title(main = sprintf("%s\n%s | %s\n(no active edges)", panel$title, label, panel$grp),
            col.main = group_colors[panel$grp])
      next
    }
    present <- keep_ct[keep_ct %in% rownames(mat_plot)]
    use_col <- ct_colors[intersect(present, names(ct_colors))]
    col_mat <- make_source_color_matrix(mat_plot[present, present, drop = FALSE], use_col)
    n_label <- if (binary) sprintf("n edges: %d", sum(mat_plot > 0)) else
      sprintf("max prob: %.4f", max(mat_plot))
    tryCatch({
      circos.clear()
      chordDiagram(mat_plot[present, present, drop = FALSE],
                   grid.col = use_col, col = col_mat, transparency = 0.3,
                   annotationTrack = c("grid", "name"),
                   preAllocateTracks = list(track.height = 0.08),
                   directional = 1, direction.type = "arrows",
                   link.arr.type = "big.arrow", self.link = 1, scale = FALSE)
      circos.clear()
      title(main = sprintf("%s\n%s | %s\n%s\n%s",
                           panel$title, label, panel$grp, niche_label, n_label),
            col.main = group_colors[panel$grp], cex.main = 1.1)
    }, error = function(e) {
      circos.clear(); plot.new()
      title(main = sprintf("%s | %s\nError: %s", label, panel$grp, e$message))
    })
  }
  dev.off()
  cat(sprintf("      Saved: %s\n", basename(out_path)))
}

# Draw one netVisual_aggregate circle panel
draw_circle <- function(group_label, niche_label, pathway) {
  obj <- get_avg_cellchat_obj(group_label, niche_label)
  if (is.null(obj)) {
    plot.new()
    title(main = sprintf("%s | %s\n(no data)", pathway, group_label),
          col.main = group_colors[group_label])
    return(invisible(NULL))
  }
  if (!pathway %in% obj@netP$pathways) {
    plot.new()
    title(main = sprintf("%s | %s\n(not detected)", pathway, group_label),
          col.main = group_colors[group_label])
    return(invisible(NULL))
  }
  tryCatch({
    netVisual_aggregate(obj, signaling = pathway, layout = "circle",
                        color.use = scPalette(length(unique(obj@idents))) |>
                          setNames(levels(obj@idents)),
                        vertex.label.cex = 0.9, arrow.size = 0.3)
    title(main = sprintf("%s | %s", pathway, group_label),
          col.main = group_colors[group_label], cex.main = 1.2, line = -1)
  }, error = function(e) {
    plot.new()
    title(main = sprintf("%s | %s\nError: %s", pathway, group_label, e$message))
    cat(sprintf("      Circle error %s %s: %s\n", pathway, group_label, e$message))
  })
}

# ── Wilcoxon test functions ───────────────────────────────────────────────────

run_wilcox_lr <- function(lr_data, group1, group2, patients_g1, patients_g2,
                          presence_map) {
  cat(sprintf("   Testing %s vs %s...\n", group1, group2))
  key_cols <- intersect(c("interaction_name", "pathway_name", "ligand", "receptor",
                           "niche", "source", "target"), colnames(lr_data))
  sub    <- lr_data[lr_data$group %in% c(group1, group2), ]
  combos <- unique(sub[, key_cols, drop = FALSE])
  pt_df  <- data.frame(patient = c(patients_g1, patients_g2),
                        group   = c(rep(group1, length(patients_g1)),
                                    rep(group2, length(patients_g2))),
                        stringsAsFactors = FALSE)
  grid <- merge(combos, pt_df, by = character())
  grid <- merge(grid, sub[, c(key_cols, "patient", "prob")],
                by = c(key_cols, "patient"), all.x = TRUE)
  grid <- merge(grid, presence_map[, c("patient", "niche", "niche_present")],
                by = c("patient", "niche"), all.x = TRUE)
  grid$niche_present[is.na(grid$niche_present)] <- FALSE
  grid$zero_type <- dplyr::case_when(
    !is.na(grid$prob)   ~ "observed",
    grid$niche_present  ~ "tested_zero",
    !grid$niche_present ~ "absent_niche"
  )
  grid$prob[is.na(grid$prob)] <- 0
  cat(sprintf("      Grid: %s rows | observed=%d tested_zero=%d absent_niche=%d\n",
              format(nrow(grid), big.mark=","),
              sum(grid$zero_type=="observed"),
              sum(grid$zero_type=="tested_zero"),
              sum(grid$zero_type=="absent_niche")))
  combo_key  <- do.call(paste, c(grid[, key_cols, drop = FALSE], sep = "|||"))
  results_list <- lapply(split(grid, combo_key), function(d) {
    meta  <- d[1, key_cols, drop = FALSE]
    x     <- d$prob[d$group == group1]; y <- d$prob[d$group == group2]
    wt    <- safe_wilcox(x, y)
    zt_g1 <- d$zero_type[d$group == group1]; zt_g2 <- d$zero_type[d$group == group2]
    data.frame(meta,
               mean_prob_g1 = mean(x, na.rm=TRUE), mean_prob_g2 = mean(y, na.rm=TRUE),
               mean_diff    = mean(y, na.rm=TRUE) - mean(x, na.rm=TRUE),
               log2fc       = log2((mean(y,na.rm=TRUE)+1e-6)/(mean(x,na.rm=TRUE)+1e-6)),
               wilcox_p = wt$p, rank_biserial_r = wt$r,
               n_patients_g1=length(x), n_patients_g2=length(y),
               n_observed_g1=sum(zt_g1=="observed"), n_tested_zero_g1=sum(zt_g1=="tested_zero"),
               n_absent_niche_g1=sum(zt_g1=="absent_niche"),
               n_observed_g2=sum(zt_g2=="observed"), n_tested_zero_g2=sum(zt_g2=="tested_zero"),
               n_absent_niche_g2=sum(zt_g2=="absent_niche"), stringsAsFactors=FALSE)
  })
  result_df <- do.call(rbind, results_list)
  rownames(result_df) <- NULL
  names(result_df)[names(result_df)=="mean_prob_g1"] <- paste0("mean_prob_", group1)
  names(result_df)[names(result_df)=="mean_prob_g2"] <- paste0("mean_prob_", group2)
  prob_g1 <- result_df[[paste0("mean_prob_", group1)]]
  prob_g2 <- result_df[[paste0("mean_prob_", group2)]]
  result_df$tier            <- assign_tier(result_df$rank_biserial_r, result_df$wilcox_p)
  result_df$pattern         <- assign_pattern(prob_g1, prob_g2, group2)
  result_df$max_prob        <- pmax(prob_g1, prob_g2)
  result_df$direction       <- ifelse(result_df$rank_biserial_r > 0,
                                      paste0("higher_in_", group2), paste0("lower_in_", group2))
  result_df$has_absent_niche <- result_df$n_absent_niche_g1 > 0 | result_df$n_absent_niche_g2 > 0
  result_df[order(result_df$tier, -abs(result_df$rank_biserial_r), result_df$wilcox_p), ]
}

aggregate_lr_niche <- function(lr_data, group1, group2, patients_g1, patients_g2,
                               presence_map) {
  key_cols <- intersect(c("interaction_name", "pathway_name", "ligand", "receptor", "niche"),
                        colnames(lr_data))
  pt_df <- data.frame(patient = c(patients_g1, patients_g2),
                       group   = c(rep(group1, length(patients_g1)),
                                   rep(group2, length(patients_g2))),
                       stringsAsFactors = FALSE)
  sub <- lr_data[lr_data$group %in% c(group1, group2), ]
  agg <- sub |>
    group_by(across(all_of(c(key_cols, "patient", "group")))) |>
    summarise(prob_sum = sum(prob, na.rm = TRUE), .groups = "drop") |>
    as.data.frame()
  grid <- merge(unique(agg[, key_cols, drop = FALSE]), pt_df, by = character())
  grid <- merge(grid, agg[, c(key_cols, "patient", "prob_sum")],
                by = c(key_cols, "patient"), all.x = TRUE)
  grid <- merge(grid, presence_map[, c("patient", "niche", "niche_present")],
                by = c("patient", "niche"), all.x = TRUE)
  grid$niche_present[is.na(grid$niche_present)] <- FALSE
  grid$zero_type <- dplyr::case_when(
    !is.na(grid$prob_sum) ~ "observed",
    grid$niche_present    ~ "tested_zero",
    !grid$niche_present   ~ "absent_niche"
  )
  grid$prob_sum[is.na(grid$prob_sum)] <- 0
  cat(sprintf("      Aggregated: %s rows | %d combos\n",
              format(nrow(grid), big.mark=","),
              length(unique(do.call(paste, c(grid[,key_cols,drop=FALSE], sep="|||"))))))
  combo_key  <- do.call(paste, c(grid[, key_cols, drop = FALSE], sep = "|||"))
  results_list <- lapply(split(grid, combo_key), function(d) {
    meta  <- d[1, key_cols, drop = FALSE]
    x     <- d$prob_sum[d$group == group1]; y <- d$prob_sum[d$group == group2]
    wt    <- safe_wilcox(x, y)
    zt_g1 <- d$zero_type[d$group == group1]; zt_g2 <- d$zero_type[d$group == group2]
    data.frame(meta,
               mean_prob_sum_g1=mean(x,na.rm=TRUE), mean_prob_sum_g2=mean(y,na.rm=TRUE),
               mean_diff=mean(y,na.rm=TRUE)-mean(x,na.rm=TRUE),
               log2fc=log2((mean(y,na.rm=TRUE)+1e-6)/(mean(x,na.rm=TRUE)+1e-6)),
               wilcox_p=wt$p, rank_biserial_r=wt$r,
               n_patients_g1=length(x), n_patients_g2=length(y),
               n_observed_g1=sum(zt_g1=="observed"), n_tested_zero_g1=sum(zt_g1=="tested_zero"),
               n_absent_niche_g1=sum(zt_g1=="absent_niche"),
               n_observed_g2=sum(zt_g2=="observed"), n_tested_zero_g2=sum(zt_g2=="tested_zero"),
               n_absent_niche_g2=sum(zt_g2=="absent_niche"), stringsAsFactors=FALSE)
  })
  result_df <- do.call(rbind, results_list)
  rownames(result_df) <- NULL
  names(result_df)[names(result_df)=="mean_prob_sum_g1"] <- paste0("mean_prob_sum_", group1)
  names(result_df)[names(result_df)=="mean_prob_sum_g2"] <- paste0("mean_prob_sum_", group2)
  result_df$tier            <- assign_tier(result_df$rank_biserial_r, result_df$wilcox_p)
  result_df$direction       <- ifelse(result_df$rank_biserial_r > 0,
                                      paste0("higher_in_", group2), paste0("lower_in_", group2))
  result_df$has_absent_niche <- result_df$n_absent_niche_g1 > 0 | result_df$n_absent_niche_g2 > 0
  result_df[order(result_df$tier, -abs(result_df$rank_biserial_r), result_df$wilcox_p), ]
}

write_tier_tables <- function(ct_df, agg_df, group1, group2, comp_label, out_dir) {
  for (level in c("celltype_pair", "lr_niche_aggregated")) {
    df <- if (level == "celltype_pair") ct_df else agg_df
    t1 <- df[df$tier == "1_top", ]
    t2 <- df[df$tier %in% c("1_top", "2_large_effect"), ]
    if ("pattern" %in% names(t1)) {
      t1_gl <- t1[t1$pattern != "quantitative_change", ]
      write.csv(t1_gl, file.path(out_dir, sprintf("%s_%s_tier1_gained_lost.csv",
                                                   comp_label, level)), row.names=FALSE)
    }
    write.csv(t1, file.path(out_dir, sprintf("%s_%s_tier1.csv",   comp_label, level)), row.names=FALSE)
    write.csv(t2, file.path(out_dir, sprintf("%s_%s_tier1_2.csv", comp_label, level)), row.names=FALSE)
    cat(sprintf("      [%s] T1: %d | T1+2: %d\n", level, nrow(t1), nrow(t2)))
  }
}

# ── SEMA pathway recoding (F10 only) ─────────────────────────────────────────
# Collapse SEMA3 and SEMA6 into a single "SEMA" label for pathway-level
# figures. Applied only to the F10 summary — raw data is unchanged.
recode_sema <- function(df) {
  df$pathway_name <- ifelse(grepl("^SEMA3|^SEMA6", df$pathway_name),
                            "SEMA", df$pathway_name)
  df
}

# ── Stage 4 helpers ───────────────────────────────────────────────────────────

# Average net$weight across patients, optionally filtered to a single niche.
# When niche_filter=NULL, averages across all niches for the group.
get_avg_weight <- function(group_label, niche_filter = NULL, fdf = file_df) {
  sub <- if (!is.null(niche_filter))
    fdf[fdf$group == group_label & fdf$niche == niche_filter, ]
  else fdf[fdf$group == group_label, ]
  mats <- list()
  for (fn in sub$fname) {
    obj <- load_obj(fn); if (is.null(obj)) next
    if (!is.null(obj@net$weight)) mats[[length(mats) + 1]] <- obj@net$weight
    rm(obj); gc()
  }
  if (length(mats) == 0) return(NULL)
  all_ct <- sort(unique(unlist(lapply(mats, rownames))))
  Reduce("+", lapply(mats, pad_matrix, all_types = all_ct)) / length(mats)
}

# Normalize a weight matrix by sender x receiver cell type proportions so
# rare cell types (Neural, LEC) are comparable to abundant ones (CM, EC, FB).
# Divides each edge weight by (prop_sender * prop_receiver).
normalize_by_proportion <- function(weight_mat, group_label,
                                     niche_filter = NULL, fdf = file_df) {
  if (is.null(weight_mat)) return(NULL)
  sub <- if (!is.null(niche_filter))
    fdf[fdf$group == group_label & fdf$niche == niche_filter, ]
  else fdf[fdf$group == group_label, ]
  ct_counts <- list()
  for (fn in sub$fname) {
    obj <- load_obj(fn); if (is.null(obj)) next
    counts <- table(obj@idents)
    ct_counts[[length(ct_counts) + 1]] <- setNames(as.numeric(counts), names(counts))
    rm(obj); gc()
  }
  if (length(ct_counts) == 0) return(weight_mat)
  all_ct <- sort(unique(unlist(lapply(ct_counts, names))))
  pooled <- setNames(numeric(length(all_ct)), all_ct)
  for (v in ct_counts) for (ct in names(v)) pooled[ct] <- pooled[ct] + v[ct]
  total <- sum(pooled); if (total == 0) return(weight_mat)
  props     <- pooled / total
  ct_in     <- rownames(weight_mat)
  props_sub <- props[ct_in]; props_sub[is.na(props_sub)] <- 1e-6
  norm_mat  <- weight_mat
  for (i in seq_along(ct_in))
    for (j in seq_along(ct_in))
      norm_mat[i, j] <- weight_mat[i, j] / max(props_sub[i] * props_sub[j], 1e-10)
  norm_mat
}

# Draw a single chord into an open device for 3-panel chord PDFs (Stage 4 Fig 2)
draw_one_chord <- function(mat, group_label, title_suffix = "",
                            min_frac = 0.02, ct_colors = NULL) {
  if (is.null(mat) || sum(mat, na.rm = TRUE) == 0) {
    plot.new()
    title(main = sprintf("%s\n(no signal)", group_label),
          col.main = group_colors[group_label])
    return(invisible(NULL))
  }
  keep_ct <- rownames(mat)[rowSums(mat) > 0 | colSums(mat) > 0]
  if (length(keep_ct) < 2) {
    plot.new()
    title(main = sprintf("%s\n(< 2 active cell types)", group_label),
          col.main = group_colors[group_label])
    return(invisible(NULL))
  }
  mat_sub <- mat[keep_ct, keep_ct]
  mat_sub[mat_sub < max(mat_sub, na.rm = TRUE) * min_frac] <- 0
  if (sum(mat_sub) == 0) {
    plot.new()
    title(main = sprintf("%s\n(all below threshold)", group_label),
          col.main = group_colors[group_label])
    return(invisible(NULL))
  }
  if (is.null(ct_colors))
    ct_colors <- scPalette(length(keep_ct)) |> setNames(keep_ct)
  use_col <- ct_colors[intersect(keep_ct, names(ct_colors))]
  col_mat <- make_source_color_matrix(mat_sub, use_col)
  tryCatch({
    circos.clear()
    chordDiagram(mat_sub, grid.col = use_col, col = col_mat, transparency = 0.3,
                 annotationTrack = c("grid", "name"),
                 preAllocateTracks = list(track.height = 0.08),
                 directional = 1, direction.type = "arrows",
                 link.arr.type = "big.arrow", self.link = 1, scale = FALSE)
    circos.clear()
    title(main = sprintf("%s%s", group_label, title_suffix),
          col.main = group_colors[group_label], cex.main = 1.2, line = -1)
  }, error = function(e) {
    circos.clear(); plot.new()
    title(main = sprintf("%s\nError: %s", group_label, e$message))
    cat(sprintf("   Chord error %s: %s\n", group_label, e$message))
  })
}

# Also define get_avg_cellchat_obj variant that accepts niche_filter=NULL
# for all-niche averaging (Stage 4 Fig 3 uses niche_filter=TARGET_NICHE)
# The main get_avg_cellchat_obj uses niche_label (required) — this wrapper
# makes niche_filter optional for Stage 4 compatibility.
get_avg_cellchat_obj_filtered <- function(group_label, niche_filter = NULL,
                                           fdf = file_df) {
  sub <- if (!is.null(niche_filter))
    fdf[fdf$niche == niche_filter & fdf$group == group_label, ]
  else fdf[fdf$group == group_label, ]
  if (nrow(sub) == 0) return(NULL)
  objs <- Filter(Negate(is.null), lapply(sub$fname, load_obj))
  if (length(objs) == 0) return(NULL)
  n_int    <- sapply(objs, function(o) tryCatch(sum(o@net$count), error = function(e) 0))
  template <- objs[[which.max(n_int)]]
  if (length(objs) == 1) return(template)
  all_ct     <- sort(unique(unlist(lapply(objs, function(o) rownames(o@net$count)))))
  avg_count  <- Reduce("+", lapply(objs, function(o) pad_matrix(o@net$count,  all_ct))) / length(objs)
  avg_weight <- Reduce("+", lapply(objs, function(o) pad_matrix(o@net$weight, all_ct))) / length(objs)
  if (!is.null(template@netP$prob) && length(dim(template@netP$prob)) == 3) {
    netP_arrays <- Filter(Negate(is.null), lapply(objs, function(o)
      if (!is.null(o@netP$prob) && length(dim(o@netP$prob)) == 3) o@netP$prob else NULL))
    if (length(netP_arrays) > 0) {
      all_pw <- sort(unique(unlist(lapply(netP_arrays, function(a) dimnames(a)[[3]]))))
      pad_array <- function(arr) {
        ct <- dimnames(arr)[[1]]; pw <- dimnames(arr)[[3]]
        out <- array(0, dim = c(length(all_ct), length(all_ct), length(all_pw)),
                     dimnames = list(all_ct, all_ct, all_pw))
        out[ct, ct, pw] <- arr; out
      }
      template@netP$prob     <- Reduce("+", lapply(netP_arrays, pad_array)) / length(netP_arrays)
      template@netP$pathways <- all_pw
    }
  }
  shared_ct           <- intersect(rownames(template@net$count), all_ct)
  template@net$count  <- avg_count[shared_ct, shared_ct]
  template@net$weight <- avg_weight[shared_ct, shared_ct]
  template
}

# =============================================================================
# CREATE OUTPUT DIRECTORIES
# =============================================================================

for (d in c(PP_OUT_DIR,
            file.path(PP_OUT_DIR, c('cellchat_objects','summary_tables','plots')),
            DIFF_OUT_DIR,
            file.path(DIFF_OUT_DIR, c('tables','plots', file.path('plots','bubbles'))),
            FIG_DIR)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# =============================================================================
# STAGE 1 — METADATA ASSEMBLY
# =============================================================================

cat("\n=== STAGE 1: METADATA ASSEMBLY ===\n\n")

cat("Loading objects...\n")
setwd(RAW_DIR)  # BPCells stores matrix paths relative to this directory
seurat_obj  <- readRDS(SEURAT_RAW_PATH)
subclust    <- readRDS(SEURAT_SUBCLUST_PATH)
niche_labels <- read.csv(NICHE_LABELS_PATH)

cat(sprintf("   Resegmented: %d cells\n", ncol(seurat_obj)))
cat(sprintf("   Subclustered: %d cells\n", ncol(subclust)))
cat(sprintf("   Niche CSV: %d rows\n", nrow(niche_labels)))

# Check barcode overlap
overlap_niche     <- intersect(colnames(seurat_obj), niche_labels$cell_id)
overlap_subcluster <- intersect(colnames(seurat_obj), colnames(subclust))
cat(sprintf("   Niche overlap: %d cells (%.1f%%)\n",
            length(overlap_niche), 100*length(overlap_niche)/ncol(seurat_obj)))
cat(sprintf("   Subcluster overlap: %d cells (%.1f%%)\n",
            length(overlap_subcluster), 100*length(overlap_subcluster)/ncol(seurat_obj)))

if (length(overlap_niche) == 0)
  stop("No barcode overlap between Seurat object and niche CSV — check data files.")

# Add niche assignments
seurat_obj$niche_kmeans_15 <- setNames(niche_labels$kmeans_15,
                                        niche_labels$cell_id)[colnames(seurat_obj)]
seurat_obj$niche_label     <- setNames(niche_labels$niche_label,
                                        niche_labels$cell_id)[colnames(seurat_obj)]

cat(sprintf("   niche_kmeans_15: %d non-NA\n", sum(!is.na(seurat_obj$niche_kmeans_15))))
cat(sprintf("   niche_label:     %d non-NA\n", sum(!is.na(seurat_obj$niche_label))))

# Transfer new metadata columns from subclustered object
new_cols <- setdiff(colnames(subclust@meta.data), colnames(seurat_obj@meta.data))
cat(sprintf("   Transferring %d new metadata columns from subclustered object\n",
            length(new_cols)))
for (col in new_cols) {
  vals <- setNames(subclust@meta.data[[col]], colnames(subclust))
  seurat_obj@meta.data[[col]] <- vals[colnames(seurat_obj)]
  cat(sprintf("     %-30s %d non-NA\n", col,
              sum(!is.na(seurat_obj@meta.data[[col]]))))
}

# Niche distribution summary
cat("\nNiche distribution:\n")
niche_counts <- sort(table(seurat_obj$niche_label[!is.na(seurat_obj$niche_label)]),
                     decreasing = TRUE)
for (i in seq_along(niche_counts))
  cat(sprintf("  %2d. %-35s: %6d cells (%.1f%%)\n", i, names(niche_counts)[i],
              niche_counts[i], 100*niche_counts[i]/sum(niche_counts)))

seurat_obj <- UpdateSeuratObject(seurat_obj)
cat("\nEnriched Seurat object ready (not saved to disk — used in-memory only).\n")

rm(subclust, niche_labels); gc()

# =============================================================================
# STAGE 2 — CELLCHAT PER-PATIENT PER-NICHE
# =============================================================================

cat("\n=== STAGE 2: CELLCHAT PER-PATIENT PER-NICHE ===\n\n")
setwd(RAW_DIR)  # BPCells needs this working directory for relative matrix paths

cat("Loading CellChat database...\n")
# Use standard CellChatDB.human — no custom additions
cellchat_db <- CellChatDB.human
cat(sprintf("   %d interactions, %d pathways\n",
            nrow(cellchat_db$interaction),
            length(unique(cellchat_db$interaction$pathway_name))))

# UpdateSeuratObject before subset — subset() validates the FOV slot
# during the call itself, so update must happen first
seurat_obj <- UpdateSeuratObject(seurat_obj)
# Filter to high-confidence cells
high_conf  <- !is.na(seurat_obj$niche_kmeans_15) &
               seurat_obj$cell_types_manual != "Unassigned"
seurat_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[high_conf])
cat(sprintf("   After filtering: %s cells\n", format(ncol(seurat_obj), big.mark=",")))

# Apply cell type grouping (vectorised with sapply)
seurat_obj$cell_type_grouped <- sapply(seurat_obj$cell_types_manual,
                                        create_cell_type_groups)

cat("   Cell type groups:\n")
for (ct in sort(unique(seurat_obj$cell_type_grouped)))
  cat(sprintf("     %-20s %s\n", ct,
              format(sum(seurat_obj$cell_type_grouped == ct), big.mark=",")))

# Materialise expression matrix once — reused for all patient x niche combos
# to avoid repeated GetAssayData calls which re-read from BPCells on disk
cat("\nMaterialising expression matrix...\n")
full_mat <- as(GetAssayData(seurat_obj, assay = "Xenium", layer = "data"), "dgCMatrix")
cat(sprintf("   %d genes x %d cells | sparsity %.1f%%\n",
            nrow(full_mat), ncol(full_mat),
            100*(1 - Matrix::nnzero(full_mat)/prod(dim(full_mat)))))

# Build valid combos
patients <- sort(unique(seurat_obj$patient))
niches   <- sort(unique(seurat_obj$niche_label[!is.na(seurat_obj$niche_label)]))
cat(sprintf("\nPatients: %d | Niches: %d\n", length(patients), length(niches)))

combos <- expand.grid(patient = patients, niche = niches, stringsAsFactors = FALSE)
combos$keep <- FALSE; combos$n_cells <- 0
combos$n_cell_types <- 0; combos$filter_reason <- ""

for (i in seq_len(nrow(combos))) {
  mask    <- seurat_obj$patient == combos$patient[i] &
             seurat_obj$niche_label == combos$niche[i]
  n_cells <- sum(mask)
  if (n_cells < COMBO_MIN_CELLS) {
    combos$filter_reason[i] <- "too_few_cells"; combos$n_cells[i] <- n_cells; next
  }
  ct_counts <- table(seurat_obj$cell_type_grouped[mask])
  n_valid   <- sum(ct_counts >= COMBO_MIN_CT)
  if (n_valid < COMBO_MIN_CT_N) {
    combos$filter_reason[i] <- "insufficient_cell_types"; combos$n_cells[i] <- n_cells; next
  }
  combos$keep[i] <- TRUE; combos$n_cells[i] <- n_cells
  combos$n_cell_types[i] <- n_valid; combos$filter_reason[i] <- "pass"
}

valid_combos <- combos[combos$keep, ]
cat(sprintf("Valid combos: %d / %d\n", nrow(valid_combos), nrow(combos)))
write.csv(combos, file.path(PP_OUT_DIR,'summary_tables','filtering_summary.csv'),
          row.names = FALSE)

# Checkpoint: skip already-completed combos with interactions
existing_files <- list.files(file.path(PP_OUT_DIR,'cellchat_objects'),
                              pattern="\\.rds$", full.names=TRUE)
completed_names <- character(0)
if (length(existing_files) > 0) {
  has_interactions <- vapply(existing_files, function(path) {
    cc <- readRDS(path); !is.null(cc@net$count) && sum(cc@net$count) > 0
  }, logical(1))
  completed_names <- basename(existing_files)[has_interactions]
  cat(sprintf("Existing with interactions: %d\n", sum(has_interactions)))
}

valid_combos$completed <- vapply(seq_len(nrow(valid_combos)), function(i) {
  sprintf("patient_%s_niche_%s.rds",
          valid_combos$patient[i], gsub(" ","_",valid_combos$niche[i])) %in% completed_names
}, logical(1))

remaining <- valid_combos[!valid_combos$completed, ]
cat(sprintf("Remaining: %d / %d\n", nrow(remaining), nrow(valid_combos)))

if (nrow(remaining) > 0) {
  cat(sprintf("\nRunning CellChat on %d combos...\n", nrow(remaining)))
  start_time <- Sys.time()

  for (i in seq_len(nrow(remaining))) {
    patient_id  <- remaining$patient[i]
    niche_label <- remaining$niche[i]
    out_name    <- sprintf("patient_%s_niche_%s", patient_id, gsub(" ","_",niche_label))
    cat(sprintf("   [%d/%d] Patient %s | Niche: %s\n",
                i, nrow(remaining), patient_id, niche_label))
    tryCatch({
      mask        <- seurat_obj$patient == patient_id & seurat_obj$niche_label == niche_label
      cells_sub   <- colnames(seurat_obj)[mask]
      meta_sub    <- seurat_obj@meta.data[cells_sub, ]
      valid_types <- names(table(meta_sub$cell_type_grouped)[table(meta_sub$cell_type_grouped) >= COMBO_MIN_CT])
      keep_cells  <- cells_sub[meta_sub$cell_type_grouped %in% valid_types]
      meta_keep   <- seurat_obj@meta.data[keep_cells, ]
      meta_keep$cell_type_grouped <- droplevels(as.factor(meta_keep$cell_type_grouped))
      cat(sprintf("      %d cells, %d types: %s\n", length(keep_cells), length(valid_types),
                  paste(sort(valid_types), collapse=", ")))
      cc <- createCellChat(object=full_mat[, keep_cells], meta=meta_keep,
                           group.by="cell_type_grouped")
      cc@DB <- cellchat_db
      cc     <- subsetData(cc)
      cc     <- identifyOverExpressedGenes(cc)
      cc     <- identifyOverExpressedInteractions(cc)
      cc     <- computeCommunProb(cc, type="truncatedMean", trim=CELLCHAT_TRIM,
                                  population.size=TRUE)
      cc     <- filterCommunication(cc, min.cells=CELLCHAT_MIN_CELL)
      cc     <- computeCommunProbPathway(cc)
      cc     <- aggregateNet(cc)
      elapsed <- as.numeric(difftime(Sys.time(), start_time, units="mins"))
      if (!is.null(cc@net$count) && sum(cc@net$count) > 0) {
        saveRDS(cc, file.path(PP_OUT_DIR,'cellchat_objects', paste0(out_name,'.rds')))
        df_net <- subsetCommunication(cc)
        df_net$patient <- patient_id; df_net$niche <- niche_label
        write.csv(df_net, file.path(PP_OUT_DIR,'summary_tables',
                                    paste0(out_name,'_interactions.csv')), row.names=FALSE)
        tryCatch({
          pdf(file.path(PP_OUT_DIR,'plots', paste0(out_name,'_summary.pdf')),
              width=12, height=10)
          par(mfrow=c(1,2))
          netVisual_circle(cc@net$count, weight.scale=TRUE, label.edge=FALSE,
                           title.name="Number of interactions")
          netVisual_circle(cc@net$weight, weight.scale=TRUE, label.edge=FALSE,
                           title.name="Interaction weights")
          par(mfrow=c(1,1))
          netVisual_heatmap(cc, color.heatmap="Reds", title.name="Interaction strength")
          if (length(cc@netP$pathways) > 0)
            netAnalysis_signalingRole_scatter(cc, title="Signaling roles")
          dev.off()
        }, error=function(e) { try(dev.off(),silent=TRUE) })
        cat(sprintf("      ✓ %d interactions | %.1f min\n", sum(cc@net$count), elapsed))
      } else {
        cat(sprintf("      no interactions | %.1f min\n", elapsed))
      }
      rm(cc); gc()
    }, error=function(e) cat(sprintf("      FAILED: %s\n", e$message)))
  }
  cat(sprintf("Stage 2 total: %.1f min\n",
              as.numeric(difftime(Sys.time(), start_time, units="mins"))))
}

# Stage 2 summary
all_rds <- list.files(file.path(PP_OUT_DIR,'cellchat_objects'), pattern="\\.rds$", full.names=TRUE)
results_df <- do.call(rbind, Filter(Negate(is.null), lapply(all_rds, function(f) {
  tryCatch({
    cc    <- readRDS(f); fname <- basename(f)
    parts <- strsplit(gsub("\\.rds$","",fname),"_niche_")[[1]]
    pid   <- gsub("patient_","",parts[1]); niche <- gsub("_"," ",parts[2])
    n_int <- if (!is.null(cc@net$count) && sum(cc@net$count)>0)
      nrow(subsetCommunication(cc)) else 0
    data.frame(patient=pid, niche=niche, n_interactions=n_int,
               n_cells=length(cc@idents), n_cell_types=length(unique(cc@idents)),
               stringsAsFactors=FALSE)
  }, error=function(e) NULL)
})))

if (!is.null(results_df) && nrow(results_df) > 0) {
  pt_grp <- unique(data.frame(patient=as.character(seurat_obj$patient),
                               group=as.character(seurat_obj$group)))
  results_df <- merge(results_df, pt_grp, by="patient", all.x=TRUE)
  write.csv(results_df, file.path(PP_OUT_DIR,'summary_tables','analysis_results_summary.csv'),
            row.names=FALSE)
  cat(sprintf("Analyses: %d | With interactions: %d (%.1f%%)\n",
              nrow(results_df), sum(results_df$n_interactions>0),
              100*mean(results_df$n_interactions>0)))
}

rm(full_mat); gc()

# =============================================================================
# STAGE 3 — DIFFERENTIAL ANALYSIS + FIGURES
# =============================================================================

cat("\n=== STAGE 3: DIFFERENTIAL ANALYSIS + FIGURES ===\n\n")

# Build file catalogues
all_files <- list.files(file.path(PP_OUT_DIR,'cellchat_objects'),
                         pattern="\\.rds$", full.names=FALSE)
file_df <- data.frame(
  fname   = all_files,
  patient = gsub("patient_([0-9]+)_niche_.*","\\1", all_files),
  niche   = gsub("_"," ", gsub(".*_niche_(.*)\\.rds","\\1", all_files)),
  stringsAsFactors = FALSE)
file_df$group <- NA
for (grp in names(patient_groups))
  file_df$group[file_df$patient %in% patient_groups[[grp]]] <- grp

niches <- sort(unique(file_df$niche))
cat(sprintf("Per-patient: %d objects\n", nrow(file_df)))

# Patient-niche presence map
presence_map <- file_df[, c("patient","niche","group")]
presence_map$niche_present <- TRUE
all_patients <- unlist(patient_groups)
full_grid    <- expand.grid(patient=all_patients, niche=niches, stringsAsFactors=FALSE)
full_grid$group <- NA
for (grp in names(patient_groups))
  full_grid$group[full_grid$patient %in% patient_groups[[grp]]] <- grp
presence_map_full <- merge(full_grid, presence_map[,c("patient","niche","niche_present")],
                            by=c("patient","niche"), all.x=TRUE)
presence_map_full$niche_present[is.na(presence_map_full$niche_present)] <- FALSE

# Load all interaction data
cat("\nLoading interaction CSVs...\n")
all_lr <- do.call(rbind, Filter(Negate(is.null), lapply(seq_len(nrow(file_df)), function(i) {
  obj <- load_obj(file_df$fname[i]); if (is.null(obj)) return(NULL)
  df  <- extract_lr_df(obj, file_df$patient[i], file_df$niche[i])
  if (!is.null(df)) df$group <- file_df$group[i]
  rm(obj); gc(); df
})))
cat(sprintf("   Per-patient rows: %s\n", format(nrow(all_lr), big.mark=",")))

# Ensure required columns exist
for (col in c("interaction_name","pathway_name","ligand","receptor")) {
  if (!col %in% colnames(all_lr))
    all_lr[[col]] <- if (col=="interaction_name")
      paste(all_lr$ligand, all_lr$receptor, sep="_") else NA_character_
}

# Differential comparisons
comparisons <- list(
  list(group1="NF",  group2="pRV", label="NF_vs_pRV"),
  list(group1="pRV", group2="RVF", label="pRV_vs_RVF")
)

cat("\nRunning differential comparisons...\n")
for (comp in comparisons) {
  g1 <- comp$group1; g2 <- comp$group2; cl <- comp$label
  cat(sprintf("\n--- %s ---\n", cl))
  out_bubbles <- file.path(DIFF_OUT_DIR,'plots','bubbles',cl)
  dir.create(out_bubbles, showWarnings=FALSE, recursive=TRUE)

  ct_results  <- run_wilcox_lr(all_lr, g1, g2, patient_groups[[g1]], patient_groups[[g2]],
                                presence_map_full)
  agg_results <- aggregate_lr_niche(all_lr, g1, g2, patient_groups[[g1]], patient_groups[[g2]],
                                    presence_map_full)
  write.csv(ct_results,
            file.path(DIFF_OUT_DIR,'tables',sprintf('%s_wilcox_celltype_pair_level.csv',cl)),
            row.names=FALSE)
  write.csv(agg_results,
            file.path(DIFF_OUT_DIR,'tables',sprintf('%s_wilcox_lr_niche_aggregated.csv',cl)),
            row.names=FALSE)
  write_tier_tables(ct_results, agg_results, g1, g2, cl,
                    file.path(DIFF_OUT_DIR,'tables'))
  cat(sprintf("   CT rows: %d | T1: %d | T1+2: %d\n", nrow(ct_results),
              sum(ct_results$tier=="1_top"),
              sum(ct_results$tier %in% c("1_top","2_large_effect"))))

  # Bubble plots per niche
  for (niche in niches) {
    plot_df <- ct_results[ct_results$niche == niche, ]
    if (nrow(plot_df) == 0 || !"source" %in% colnames(plot_df)) next
    plot_df <- head(plot_df[order(plot_df$tier,-abs(plot_df$rank_biserial_r),
                                   plot_df$wilcox_p),], 30)
    plot_df$lr_pair   <- paste0(plot_df$interaction_name,"\n",plot_df$source,"->",plot_df$target)
    plot_df$direction <- ifelse(plot_df$rank_biserial_r > 0,
                                paste0("^ in ",g2), paste0("^ in ",g1))
    plot_df$tier_label <- ifelse(plot_df$tier=="1_top","T1",
                                  ifelse(plot_df$tier=="2_large_effect","T2",""))
    p_b <- ggplot(plot_df, aes(x=reorder(lr_pair,rank_biserial_r),
                                y=rank_biserial_r, fill=direction)) +
      geom_col(width=0.65,alpha=0.85) +
      geom_hline(yintercept=0,linetype="dashed",color="grey50") +
      coord_flip() +
      scale_fill_manual(values=setNames(c(group_colors[g2],group_colors[g1]),
                                         c(paste0("^ in ",g2),paste0("^ in ",g1)))) +
      labs(title=sprintf("%s vs %s - %s",g1,g2,niche), x=NULL,
           y=sprintf("Rank-biserial r (positive = higher in %s)",g2), fill="Direction") +
      theme_bw(base_size=11) +
      theme(plot.title=element_text(face="bold",size=12), axis.text.y=element_text(size=7))
    pdf(file.path(out_bubbles,sprintf("bubble_%s_%s.pdf",cl,gsub(" ","_",niche))),
        width=11, height=max(6,nrow(plot_df)*0.3+3))
    print(p_b); dev.off()
  }

  # Summary heatmaps
  tryCatch({
    t12 <- agg_results[agg_results$tier %in% c("1_top","2_large_effect"),]
    use_df <- if (nrow(t12)>=5) t12 else agg_results
    agg_wide <- use_df |> select(interaction_name,niche,mean_diff) |>
      pivot_wider(names_from=niche, values_from=mean_diff, values_fill=0)
    mat_A <- as.matrix(agg_wide[,-1]); rownames(mat_A) <- agg_wide$interaction_name
    if (nrow(mat_A) > 50) mat_A <- mat_A[order(apply(mat_A,1,var),decreasing=TRUE)[1:50],,drop=FALSE]
    max_abs <- max(abs(mat_A),na.rm=TRUE); if(max_abs==0) max_abs <- 1
    pdf(file.path(DIFF_OUT_DIR,'plots',sprintf("heatmap_LR_x_niche_%s.pdf",cl)),
        width=max(8,ncol(mat_A)*0.6+3), height=max(8,nrow(mat_A)*0.22+3))
    pheatmap::pheatmap(mat_A,
      color=colorRampPalette(c(group_colors[g1],"white",group_colors[g2]))(100),
      breaks=seq(-max_abs,max_abs,length.out=101),
      main=sprintf("LR x Niches | mean prob diff\n%s",cl),
      fontsize_row=7,fontsize_col=9,border_color=NA,
      cluster_rows=TRUE,cluster_cols=TRUE,angle_col=45)
    dev.off()
  }, error=function(e) { try(dev.off(),silent=TRUE) })

  tryCatch({
    pw_niche <- ct_results |>
      filter(tier %in% c("1_top","2_large_effect"), !is.na(pathway_name)) |>
      group_by(pathway_name,niche) |> summarise(n=n(),.groups="drop") |>
      pivot_wider(names_from=niche,values_from=n,values_fill=0) |> as.data.frame()
    mat_B <- as.matrix(pw_niche[,-1]); rownames(mat_B) <- pw_niche$pathway_name
    if (nrow(mat_B) > 0) {
      pdf(file.path(DIFF_OUT_DIR,'plots',sprintf("heatmap_pathway_x_niche_%s.pdf",cl)),
          width=max(8,ncol(mat_B)*0.6+3), height=max(6,nrow(mat_B)*0.28+3))
      pheatmap::pheatmap(mat_B,
        color=colorRampPalette(c("white","#FDBB84","#D7301F"))(50),
        main=sprintf("Pathways x Niches | N T1+2 LR pairs\n%s",cl),
        fontsize_row=8,fontsize_col=9,border_color=NA,
        cluster_rows=TRUE,cluster_cols=TRUE,angle_col=45,
        display_numbers=TRUE,number_color="grey20")
      dev.off()
    }
  }, error=function(e) { try(dev.off(),silent=TRUE) })
}

# Master table
master_list <- lapply(comparisons, function(comp) {
  f <- file.path(DIFF_OUT_DIR,'tables',
                 sprintf('%s_wilcox_celltype_pair_level.csv',comp$label))
  if (!file.exists(f)) return(NULL)
  df <- read.csv(f, stringsAsFactors=FALSE)
  df$comparison <- comp$label; df$group1_label <- comp$group1; df$group2_label <- comp$group2
  prob_cols <- grep("^mean_prob_",names(df),value=TRUE)
  if (length(prob_cols)==2) {
    names(df)[names(df)==prob_cols[1]] <- "mean_prob_g1"
    names(df)[names(df)==prob_cols[2]] <- "mean_prob_g2"
  }
  df
})
master_df <- dplyr::bind_rows(Filter(Negate(is.null), master_list))
master_df <- master_df[order(master_df$tier,-abs(master_df$rank_biserial_r),master_df$wilcox_p),]
write.csv(master_df, file.path(DIFF_OUT_DIR,'tables','master_all_comparisons_celltype_level.csv'),
          row.names=FALSE)
cat(sprintf("Master table: %s rows\n", format(nrow(master_df),big.mark=",")))

# =============================================================================
# SUPPLEMENTARY TABLE — Niche-specific communication patterns
# =============================================================================
# One row per ligand-receptor pair x niche x condition, reporting:
#   - Core identifiers (niche, condition, ligand, receptor, interaction,
#     pathway, source cell type, target cell type)
#   - Mean communication probability per condition
#   - Number of patients with detected signal
#   - Wilcoxon effect size and p-value for NF vs pRV and pRV vs RVF
#   - Tier and directional pattern
# This is formatted from master_df joined back to per-condition mean probs
# and the per-comparison Wilcoxon outputs.

cat("\nBuilding supplementary table...\n")

# Compute per-condition mean probability from all_lr (per-patient LR data)
per_cond_mean <- all_lr |>
  group_by(interaction_name, pathway_name, ligand, receptor,
           source, target, niche, group) |>
  summarise(
    mean_prob       = mean(prob, na.rm = TRUE),
    n_patients_det  = sum(prob > 0, na.rm = TRUE),
    .groups = "drop"
  )

# Pivot conditions wide so each row has NF/pRV/RVF side by side
per_cond_wide <- per_cond_mean |>
  pivot_wider(
    names_from  = group,
    values_from = c(mean_prob, n_patients_det),
    values_fill = 0
  )

# Load both comparison Wilcoxon tables and join
comp_tables <- lapply(list(
  list(label = "NF_vs_pRV",  g1 = "NF",  g2 = "pRV"),
  list(label = "pRV_vs_RVF", g1 = "pRV", g2 = "RVF")
), function(comp) {
  f <- file.path(DIFF_OUT_DIR, "tables",
                 sprintf("%s_wilcox_celltype_pair_level.csv", comp$label))
  if (!file.exists(f)) return(NULL)
  df <- read.csv(f, stringsAsFactors = FALSE)
  # Rename group-specific prob cols to generic names
  prob_cols <- grep("^mean_prob_", names(df), value = TRUE)
  if (length(prob_cols) == 2) {
    names(df)[names(df) == prob_cols[1]] <- "mean_prob_g1"
    names(df)[names(df) == prob_cols[2]] <- "mean_prob_g2"
  }
  df$comparison <- comp$label
  # Keep only the stats columns we want
  keep <- c("interaction_name", "pathway_name", "ligand", "receptor",
            "source", "target", "niche",
            "wilcox_p", "rank_biserial_r", "tier", "pattern",
            "direction", "has_absent_niche", "comparison")
  df[, intersect(keep, names(df))]
})

stats_df <- bind_rows(Filter(Negate(is.null), comp_tables))

# Pivot stats wide by comparison
stats_wide <- stats_df |>
  pivot_wider(
    names_from  = comparison,
    values_from = c(wilcox_p, rank_biserial_r, tier, pattern, direction),
    names_sep   = "__"
  )

# Join everything together
supp_table <- per_cond_wide |>
  left_join(stats_wide,
            by = c("interaction_name", "pathway_name", "ligand", "receptor",
                   "source", "target", "niche")) |>
  # Clean up column order
  select(
    niche, interaction_name, pathway_name, ligand, receptor,
    source, target,
    any_of(paste0("mean_prob_", group_order)),
    any_of(paste0("n_patients_det_", group_order)),
    any_of(c("wilcox_p__NF_vs_pRV", "rank_biserial_r__NF_vs_pRV",
             "tier__NF_vs_pRV", "pattern__NF_vs_pRV",
             "wilcox_p__pRV_vs_RVF", "rank_biserial_r__pRV_vs_RVF",
             "tier__pRV_vs_RVF", "pattern__pRV_vs_RVF")),
    any_of("has_absent_niche")
  ) |>
  # Rename for readability
  rename_with(~ gsub("mean_prob_",     "mean_prob_",       .x)) |>
  rename_with(~ gsub("n_patients_det_", "n_patients_",     .x)) |>
  rename_with(~ gsub("wilcox_p__",      "wilcox_p_",       .x)) |>
  rename_with(~ gsub("rank_biserial_r__","effect_size_r_", .x)) |>
  rename_with(~ gsub("tier__",          "tier_",           .x)) |>
  rename_with(~ gsub("pattern__",       "pattern_",        .x)) |>
  arrange(niche, pathway_name, interaction_name, source, target)

# Round numeric columns for readability
numeric_cols <- sapply(supp_table, is.numeric)
supp_table[numeric_cols] <- lapply(supp_table[numeric_cols],
                                    function(x) round(x, 8))

supp_table_path <- file.path(DIFF_OUT_DIR, "tables",
                              "supplementary_table_niche_communication.csv")
write.csv(supp_table, supp_table_path, row.names = FALSE)
cat(sprintf("Supplementary table: %s rows x %s columns\n",
            format(nrow(supp_table), big.mark = ","), ncol(supp_table)))
cat(sprintf("Saved: %s\n", supp_table_path))

# =============================================================================
# FIGURES
# =============================================================================


# =============================================================================
# STAGE 4 — FIGURE GENERATION
# =============================================================================

cat("\n=== STAGE 4: FIGURE GENERATION ===\n\n")

cat("\nLoading per-patient interaction CSVs...\n")
pp_csvs <- list.files(file.path(PP_OUT_DIR, 'summary_tables'),
                       pattern = '_interactions\\.csv$', full.names = TRUE)
all_lr <- bind_rows(lapply(pp_csvs, read.csv))
all_lr$group <- NA
for (grp in names(patient_groups))
  all_lr$group[all_lr$patient %in% patient_groups[[grp]]] <- grp
cat(sprintf("   %s rows loaded\n", format(nrow(all_lr), big.mark = ",")))

# Collapse SEMA3/SEMA6 into SEMA for LR-pair display
recode_sema_lr <- function(df) {
  df$pathway_name <- ifelse(grepl("^SEMA3|^SEMA6", df$pathway_name),
                            "SEMA", df$pathway_name)
  df
}

# Checkpoint: skip Figures 1-3 if all outputs already exist (cheaper figures,
# rebuilding wouldn't take long, but skipping avoids redundant CellChat
# object averaging across all patients for fig2/fig3)
fig_outputs <- file.path(FIG_DIR, c(
  "fig1_top15_pathways_myocardium.pdf",
  "fig2_chord_myocardium_normalized.pdf",
  "fig3_sema3c_circles.pdf",
  "fig5_sema3c_quantification.pdf"
))
if (all(file.exists(fig_outputs))) {
  cat("\nFigures 1, 2, 3, 5: all outputs already exist — skipping.\n")
} else {

# =============================================================================
# FIGURE 1 — TOP 15 LR PAIRS (per-patient, SEMA collapsed)
# =============================================================================

cat("\nFigure 1: Top 15 LR pairs...\n")

# Summarise at the pathway level
# Ranked by summed communication probability across all LR pairs within the
# pathway. SEMA highlighted in red.
make_top_pathway_bar <- function(lr_data, niche_filter = NULL, top_n = 15,
                                  fig_title = "Top 15 pathways") {
  df <- if (!is.null(niche_filter))
    lr_data[lr_data$niche == niche_filter, ] else lr_data
  df <- recode_sema_lr(df)

  summ <- df |>
    group_by(group, pathway_name) |>
    summarise(
      pooled_prob = sum(prob, na.rm = TRUE),
      n_lr        = n_distinct(interaction_name),
      .groups = "drop"
    ) |>
    filter(!is.na(pathway_name), pathway_name != "") |>
    # Normalize so all pathways within each group sum to 1 —
    # each bar = fraction of total communication probability for that pathway
    group_by(group) |>
    mutate(pooled_prob = pooled_prob / sum(pooled_prob)) |>
    ungroup() |>
    mutate(group = factor(group, levels = group_order))

  # Shared x-axis limit: slightly above the max pooled_prob across all groups
  x_max <- max(summ$pooled_prob, na.rm = TRUE) * 1.08

  plots <- lapply(group_order, function(grp) {
    sub <- summ[summ$group == grp, ]
    if (nrow(sub) == 0) return(ggplot() + theme_void())
    sub <- head(sub[order(-sub$pooled_prob), ], top_n)
    sub$pathway_name <- factor(sub$pathway_name,
                                levels = rev(unique(sub$pathway_name)))
    y_colors <- ifelse(levels(sub$pathway_name) == "SEMA", "#D6604D", "black")
    y_faces  <- ifelse(levels(sub$pathway_name) == "SEMA", "bold", "plain")
    ggplot(sub, aes(x = pathway_name, y = pooled_prob, fill = pooled_prob)) +
      geom_col(width = 0.75, alpha = 0.9) +
      coord_flip() +
      scale_fill_gradient(low  = adjustcolor(group_colors[grp], 0.3),
                          high = group_colors[grp], guide = "none") +
      scale_y_continuous(limits = c(0, x_max), expand = expansion(mult = c(0, 0))) +
      labs(title = grp, x = NULL, y = "Proportion of total communication probability") +
      theme_bw(base_size = 11) +
      theme(
        plot.title  = element_text(face = "bold", size = 13,
                                    color = group_colors[grp]),
        axis.text.y = element_text(size = 8, color = y_colors, face = y_faces)
      )
  })

  wrap_plots(plots, nrow = 1)
}

p1a <- make_top_pathway_bar(all_lr, niche_filter = TARGET_NICHE,
                             fig_title = "Top 15 pathways — Myocardium")
pdf(file.path(FIG_DIR, "fig1_top15_pathways_myocardium.pdf"), width = 18, height = 6)
print(p1a); dev.off()
cat("   Saved fig1\n")

# =============================================================================
# FIGURE 2 — MYOCARDIUM CHORD DIAGRAM (NORMALIZED)
# =============================================================================


cat("\nFigure 2: Myocardium chord diagram (normalized)...\n")

cat("   Building averaged weight matrices (Myocardium)...\n")
mats_myocardium <- setNames(lapply(group_order, get_avg_weight,
                                    niche_filter = TARGET_NICHE), group_order)

cat("   Normalizing by cell type proportion...\n")
mats_myocardium_norm <- setNames(
  lapply(group_order, function(grp)
    normalize_by_proportion(mats_myocardium[[grp]], grp,
                             niche_filter = TARGET_NICHE)),
  group_order)

all_ct_global <- sort(unique(unlist(lapply(
  mats_myocardium,
  function(m) if (!is.null(m)) rownames(m) else character(0)
))))
ct_colors_global <- scPalette(length(all_ct_global)) |> setNames(all_ct_global)

write_3panel_chord <- function(mat_list, out_path, main_title,
                                subtitle = "", min_frac = 0.02) {
  has_title <- nchar(trimws(paste(main_title, subtitle))) > 0
  top_mar   <- if (has_title) 5 else 1
  pdf(out_path, width = 18, height = 6)
  par(mfrow = c(1, 3), mar = c(2, 2, top_mar, 2),
      oma = c(0, 0, if (has_title) 2 else 0, 0))
  for (grp in group_order)
    draw_one_chord(mat_list[[grp]], grp, ct_colors = ct_colors_global,
                   min_frac = min_frac)
  if (has_title)
    mtext(paste(main_title, subtitle, sep = "\n"), outer = TRUE, side = 3,
          line = 0, cex = 1.1, font = 2)
  dev.off()
  cat(sprintf("   Saved: %s\n", basename(out_path)))
}

write_3panel_chord(
  mats_myocardium_norm,
  file.path(FIG_DIR, "fig2_chord_myocardium_normalized.pdf"),
  main_title = "",
  subtitle   = ""
)

# =============================================================================
# FIGURE 3 — SEMA3C CIRCLE PLOTS (panel3_sema3c_circles equivalent)
# =============================================================================

cat("\nFigure 3: SEMA3C circle plots...\n")

pdf(file.path(FIG_DIR, "fig3_sema3c_circles.pdf"), width = 18, height = 6)
par(mfrow = c(1, 3), mar = c(2, 2, 2, 2), oma = c(0, 0, 0, 0))

for (grp in group_order) {
  cat(sprintf("   Building averaged object for %s...\n", grp))
  obj <- get_avg_cellchat_obj_filtered(grp, niche_filter = TARGET_NICHE)
  if (is.null(obj)) {
    plot.new()
    title(main = sprintf("SEMA3 | %s\n(no data)", grp), col.main = group_colors[grp])
    next
  }
  if (!"SEMA3" %in% obj@netP$pathways) {
    plot.new()
    title(main = sprintf("SEMA3 | %s\n(not detected)", grp), col.main = group_colors[grp])
    next
  }
  tryCatch({
    netVisual_aggregate(obj, signaling = "SEMA3", layout = "circle",
                        color.use = scPalette(length(unique(obj@idents))) |>
                          setNames(levels(obj@idents)),
                        vertex.label.cex = 0.9, arrow.size = 0.3)
    title(main = sprintf("SEMA3 | %s\n%s", grp, TARGET_NICHE),
          col.main = group_colors[grp], cex.main = 1.2, line = -1)
  }, error = function(e) {
    plot.new()
    title(main = sprintf("SEMA3 | %s\nError: %s", grp, e$message))
    cat(sprintf("   Circle error %s: %s\n", grp, e$message))
  })
}

dev.off()
cat("   Saved fig3_sema3c_circles.pdf\n")

# =============================================================================
# FIGURE 5 — SEMA3C SIGNALING QUANTIFICATION (group mean bars + patient dots)
# =============================================================================
# Quantify amount of SEMA3C signalling with n=3 per group, individual points on
# group mean bars are more informative than box-and-whisker plots.
#
# Three panels:
#   A: Total SEMA3 pathway probability (netP sum) — overall signal strength
#   B: Number of active sender->receiver cell type pairs — "connections"
#   C: SEMA3C-PLXND1 LR pair probability — most specific metric

cat("\nFigure 5: SEMA3C quantification...\n")

metrics_list_fig5 <- list()
for (grp_fig5 in group_order) {
  for (pid_fig5 in patient_groups[[grp_fig5]]) {
    fname_fig5 <- sprintf("patient_%s_niche_%s.rds", pid_fig5,
                           gsub(" ", "_", TARGET_NICHE))
    fpath_fig5 <- file.path(PP_OUT_DIR, "cellchat_objects", fname_fig5)
    if (!file.exists(fpath_fig5)) next
    cc_fig5 <- readRDS(fpath_fig5)

    sema3c_lr <- tryCatch({
      if (!is.null(cc_fig5@net$prob) && length(dim(cc_fig5@net$prob)) == 3 &&
          "SEMA3C_PLXND1" %in% dimnames(cc_fig5@net$prob)[[3]])
        sum(cc_fig5@net$prob[, , "SEMA3C_PLXND1"]) else 0
    }, error = function(e) 0)

    # Neural->EC: dominant sender pair in NF
    neural_ec <- tryCatch({
      if ("SEMA3" %in% cc_fig5@netP$pathways &&
          !is.null(cc_fig5@netP$prob) && length(dim(cc_fig5@netP$prob)) == 3 &&
          "Neural" %in% dimnames(cc_fig5@netP$prob)[[1]] &&
          "EC"     %in% dimnames(cc_fig5@netP$prob)[[2]])
        cc_fig5@netP$prob["Neural", "EC", "SEMA3"] else 0
    }, error = function(e) 0)

    # CM->EC: dominant sender pair in RVF
    cm_ec <- tryCatch({
      if ("SEMA3" %in% cc_fig5@netP$pathways &&
          !is.null(cc_fig5@netP$prob) && length(dim(cc_fig5@netP$prob)) == 3 &&
          "CM" %in% dimnames(cc_fig5@netP$prob)[[1]] &&
          "EC" %in% dimnames(cc_fig5@netP$prob)[[2]])
        cc_fig5@netP$prob["CM", "EC", "SEMA3"] else 0
    }, error = function(e) 0)

    metrics_list_fig5[[pid_fig5]] <- data.frame(
      patient = pid_fig5, group = grp_fig5,
      sema3c_plxnd1 = sema3c_lr,
      neural_ec     = neural_ec,
      cm_ec         = cm_ec,
      stringsAsFactors = FALSE)
    rm(cc_fig5); gc()
  }
}

metrics_fig5 <- bind_rows(metrics_list_fig5)
metrics_fig5$group <- factor(metrics_fig5$group, levels = group_order)
# Lock x-axis order to NF -> pRV -> RVF in all panels
metrics_fig5 <- metrics_fig5[order(metrics_fig5$group), ]

make_quant_panel <- function(df, y_col, y_label, y_log = TRUE,
                              floor_val = 1e-8) {
  # Split into detected (> floor) and undetected (== floor or 0) patients
  df$detected <- df[[y_col]] > floor_val

  # Group mean uses only detected patients, so pRV bar reflects true signal
  # for 1692 only (not dragged up by floor-substituted zeros)
  grp_means <- df[df$detected, ] |>
    group_by(group) |>
    summarise(mean_val = mean(.data[[y_col]]), .groups = "drop") |>
    # pRV may have 0 detected patients — ensure all groups present with NA
    right_join(data.frame(group = factor(group_order, levels = group_order)),
               by = "group")

  # Axis floor for "not detected" symbols: one log unit below lowest detected
  detected_min <- min(df[[y_col]][df$detected], na.rm = TRUE)
  nd_y <- detected_min / 20  # position for "not detected" triangles

  # Detected patient points
  df_det    <- df[df$detected,  ]
  df_nd     <- df[!df$detected, ]
  df_nd[[y_col]] <- nd_y  # place undetected at floor position

  # Enforce group order on all data slices going into the plot
  df$group           <- factor(df$group,           levels = group_order)
  grp_means$group    <- factor(grp_means$group,    levels = group_order)
  if (nrow(df_det) > 0) df_det$group <- factor(df_det$group, levels = group_order)
  if (nrow(df_nd)  > 0) df_nd$group  <- factor(df_nd$group,  levels = group_order)

  p <- ggplot(mapping = aes(x = group)) +
    # Group mean bars (only for groups with detected signal)
    geom_col(data = grp_means[!is.na(grp_means$mean_val), ],
             aes(y = mean_val, fill = group), width = 0.6, alpha = 0.7) +
    # Detected patients: jittered filled circles
    # Use position_jitter with fixed seed so points and labels stay aligned
    geom_point(data = df_det,
               aes(y = .data[[y_col]], fill = group),
               shape = 21, size = 4, color = "black", stroke = 0.8,
               alpha = 0.95,
               position = position_jitter(width = 0.12, height = 0, seed = 42)) +
    geom_text(data = df_det,
              aes(y = .data[[y_col]], label = patient),
              vjust = -0.8, size = 2.8, color = "grey30",
              position = position_jitter(width = 0.12, height = 0, seed = 42)) +
    # Undetected patients: jittered open downward triangles at floor
    geom_point(data = df_nd,
               aes(y = .data[[y_col]]),
               shape = 25, size = 3, color = "grey50", fill = "white",
               stroke = 0.8, alpha = 0.9,
               position = position_jitter(width = 0.12, height = 0, seed = 42)) +
    geom_text(data = df_nd,
              aes(y = .data[[y_col]], label = patient),
              vjust = 1.8, size = 2.5, color = "grey50",
              position = position_jitter(width = 0.12, height = 0, seed = 42)) +
    # Dashed line marking the "not detected" threshold
    geom_hline(yintercept = nd_y * 3, linetype = "dashed",
               color = "grey70", linewidth = 0.4) +
    annotate("text", x = 0.4, y = nd_y * 3,
             label = "not detected", hjust = 1, vjust = 0.5,
             size = 2.5, color = "grey60") +
    scale_fill_manual(values = group_colors, guide = "none") +
    scale_x_discrete(limits = group_order) +
    scale_y_log10(labels = scales::label_log()) +
    labs(x = NULL, y = y_label) +
    theme_bw(base_size = 12) +
    theme(panel.grid.major.x = element_blank(),
          axis.text.x = element_text(face = "bold", size = 12))
  p
}

# make_quant_panel handles undetected patients internally:
# detected patients shown as filled circles with group mean bars,
# undetected shown as open downward triangles below a "not detected" line.
# Flag zeros as floor so the function can distinguish detected vs not.
FLOOR_VAL <- 1e-8
metrics_fig5$sema3c_plxnd1[metrics_fig5$sema3c_plxnd1 == 0] <- FLOOR_VAL
metrics_fig5$neural_ec[metrics_fig5$neural_ec           == 0] <- FLOOR_VAL
metrics_fig5$cm_ec[metrics_fig5$cm_ec                   == 0] <- FLOOR_VAL

p_q_lr <- make_quant_panel(
  metrics_fig5,
  "sema3c_plxnd1", "SEMA3C-PLXND1 LR\nprobability (log10)") +
  ggtitle("SEMA3C-PLXND1 signal")

p_q_neural <- make_quant_panel(
  metrics_fig5,
  "neural_ec", "Neural -> EC\nSEMA3 prob (log10)") +
  ggtitle("Neural sender (NF-dominant)")

p_q_cm <- make_quant_panel(
  metrics_fig5,
  "cm_ec", "CM -> EC\nSEMA3 prob (log10)") +
  ggtitle("CM sender (RVF-dominant)")

p_quant <- p_q_lr + p_q_neural + p_q_cm +
  plot_annotation(
    title    = "SEMA3C signaling: overall activity and sender-identity shift",
    subtitle = sprintf("%s niche | n=3 per group | Bar = group mean (detected only) | Triangles = not detected",
                       TARGET_NICHE),
    theme = theme(plot.title    = element_text(face = "bold", size = 14),
                  plot.subtitle = element_text(size = 9, color = "grey40"))
  )

pdf(file.path(FIG_DIR, "fig5_sema3c_quantification.pdf"), width = 12, height = 5)
print(p_quant)
dev.off()
cat("   Saved fig5_sema3c_quantification.pdf\n")

}  # end Figures 1-3+5 checkpoint block

# =============================================================================
# FIGURE 4 — SEMA3C SPATIAL HEXBIN NETWORK PLOTS
# =============================================================================
# A custom hexbin-smoothed spatial network plot that:
#   1. Bins cells into a square grid (HEX_BINWIDTH microns)
#   2. Computes dominant cell type per bin
#   3. For each nonzero SEMA3 (source_type, target_type) pair from
#      cc@netP$prob[,,"SEMA3"], draws edges between nearby bins whose
#      dominant types match the source/target pair
#   4. Edge linewidth scaled by sqrt(prob) to compress the 130-fold dynamic
#      range across groups while preserving the correct rank ordering
#
# sqrt scaling rationale: SEMA3 signal spans ~5 orders of magnitude
# (1e-6 to 1e-1 across patients). Linear scaling makes NF edges invisible
# next to RVF; log scaling loses edge width meaning entirely. Square-root
# is standard practice for compressing dynamic range while maintaining
# proportionality — analogous to sqrt-transforming count data.
#
# pRV patients with no detected SEMA3 (1618, 1567) are shown with nodes
# only (no edges) to preserve spatial geometry context and make the
# signal absence visible rather than leaving them out entirely.

cat("\nFigure 4: SEMA3C spatial hexbin network plots...\n")

HEX_BINWIDTH  <- 100              # microns — smooths FOV grid artifact
EDGE_MAX_DIST <- HEX_BINWIDTH * 3 # only connect bins within this distance

# Reuse seurat_obj from Stage 2 (already filtered + cell_type_grouped)
setwd(RAW_DIR)
full_mat_spatial <- as(GetAssayData(seurat_obj, assay = "Xenium", layer = "data"),
                       "dgCMatrix")

get_spatial_coords_fig4 <- function(seurat_subset) {
  tryCatch({
    coords <- GetTissueCoordinates(seurat_subset, which = "centroids")
    m <- as.matrix(coords[, c("x", "y")]); rownames(m) <- coords$cell; m
  }, error = function(e) { cat(sprintf("   Coord error: %s\n", e$message)); NULL })
}

make_hexbin_spatial_plot <- function(patient_id, grp, cc_obj, coords_mat,
                                      meta_keep, has_sema3) {
  cell_df <- data.frame(
    x = coords_mat[, "x"], y = coords_mat[, "y"],
    cell_type = meta_keep$cell_type_grouped
  ) |> mutate(
    bin_x = round(x / HEX_BINWIDTH) * HEX_BINWIDTH,
    bin_y = round(y / HEX_BINWIDTH) * HEX_BINWIDTH
  )

  bin_summary <- cell_df |>
    group_by(bin_x, bin_y) |>
    summarise(n_cells = n(),
              dominant_type = names(sort(table(cell_type), decreasing=TRUE))[1],
              .groups = "drop")

  edge_df <- NULL
  if (has_sema3 && !is.null(cc_obj)) {
    prob_mat <- tryCatch(cc_obj@netP$prob[,,"SEMA3"], error=function(e) NULL)
    if (!is.null(prob_mat)) {
      prob_df <- as.data.frame(as.table(prob_mat))
      colnames(prob_df) <- c("source_type","target_type","prob")
      prob_df <- prob_df[prob_df$prob > 0, ]

      edge_list <- list()
      for (i in seq_len(nrow(prob_df))) {
        src <- as.character(prob_df$source_type[i])
        tgt <- as.character(prob_df$target_type[i])
        ep  <- prob_df$prob[i]
        sb  <- bin_summary[bin_summary$dominant_type == src, ]
        tb  <- bin_summary[bin_summary$dominant_type == tgt, ]
        if (nrow(sb)==0 || nrow(tb)==0) next
        for (si in seq_len(nrow(sb))) {
          dx <- tb$bin_x - sb$bin_x[si]; dy <- tb$bin_y - sb$bin_y[si]
          dist <- sqrt(dx^2 + dy^2)
          close <- if (src==tgt) which(dist<=EDGE_MAX_DIST)
                   else which(dist<=EDGE_MAX_DIST & dist>0)
          for (ti in close)
            edge_list[[length(edge_list)+1]] <- data.frame(
              x1=sb$bin_x[si], y1=sb$bin_y[si],
              x2=tb$bin_x[ti], y2=tb$bin_y[ti],
              prob=ep, sqrt_prob=sqrt(ep))
        }
      }
      if (length(edge_list) > 0) edge_df <- bind_rows(edge_list)
    }
  }

  # Bins that participate in edges (for colored highlighting)
  if (!is.null(edge_df) && nrow(edge_df) > 0) {
    edge_coords <- unique(rbind(
      setNames(edge_df[,c("x1","y1")], c("bin_x","bin_y")),
      setNames(edge_df[,c("x2","y2")], c("bin_x","bin_y"))
    ))
    bin_summary$in_edge <- paste(bin_summary$bin_x, bin_summary$bin_y) %in%
                           paste(edge_coords$bin_x, edge_coords$bin_y)
  } else {
    bin_summary$in_edge <- FALSE
  }

  subtitle_txt <- if (!has_sema3) "No SEMA3 signal detected" else
    sprintf("Edges within %dum | sqrt(prob) linewidth", EDGE_MAX_DIST)

  p <- ggplot() +
    geom_point(data = bin_summary,
               aes(x=bin_x, y=bin_y, size=n_cells),
               color="grey85", alpha=0.6)

  if (!is.null(edge_df) && nrow(edge_df) > 0) {
    p <- p +
      geom_segment(data=edge_df,
                   aes(x=x1, y=y1, xend=x2, yend=y2, linewidth=prob),
                   color="red", alpha=0.5, lineend="round") +
      scale_linewidth_continuous(range=c(0.1, 1.5), name="SEMA3 prob",
                                 guide="none")
  }

  # Only color bins that are edge endpoints — for no-signal patients,
  # show nothing colored so the plot reads as an empty tissue map
  if (!is.null(edge_df) && nrow(edge_df) > 0) {
    p <- p +
      geom_point(data=bin_summary[bin_summary$in_edge, ],
                 aes(x=bin_x, y=bin_y, fill=dominant_type, size=n_cells),
                 shape=21, color="black", alpha=0.9, stroke=0.4)
  }

  p <- p +
    scale_size_continuous(range=c(0.5, 4), guide="none") +
    coord_fixed() +
    theme_minimal(base_size=10) +
    theme(legend.position="right",
          plot.title=element_text(face="bold", size=10)) +
    labs(fill="Cell type", x=NULL, y=NULL)

  p
}

for (patient_id in unlist(patient_groups)) {
  grp <- names(which(sapply(patient_groups, function(pts) patient_id %in% pts)))

  out_full   <- file.path(FIG_DIR,
    sprintf("fig4_sema3c_hexbin_%s_%s_full.pdf",   patient_id, grp))
  out_zoomed <- file.path(FIG_DIR,
    sprintf("fig4_sema3c_hexbin_%s_%s_zoomed.pdf", patient_id, grp))

  if (file.exists(out_full) && file.exists(out_zoomed)) {
    cat(sprintf("   Patient %s (%s): already done — skipping\n", patient_id, grp))
    next
  }
  cat(sprintf("   Patient %s (%s)...\n", patient_id, grp))

  tryCatch({
    mask <- seurat_obj$patient == patient_id &
            seurat_obj$niche_label == TARGET_NICHE &
            !is.na(seurat_obj$niche_label)
    if (sum(mask) < 50) {
      cat(sprintf("      Skipping: only %d cells\n", sum(mask))); next
    }

    cells_sub   <- colnames(seurat_obj)[mask]
    meta_sub    <- seurat_obj@meta.data[cells_sub, ]
    ct_counts   <- table(meta_sub$cell_type_grouped)
    valid_types <- names(ct_counts[ct_counts >= 5])
    if (length(valid_types) < 2) {
      cat("      Skipping: insufficient cell types\n"); next
    }

    keep_cells <- cells_sub[meta_sub$cell_type_grouped %in% valid_types]
    meta_keep  <- seurat_obj@meta.data[keep_cells, ]
    meta_keep$cell_type_grouped <- droplevels(as.factor(meta_keep$cell_type_grouped))
    meta_keep$samples <- patient_id

    seurat_sub <- subset(seurat_obj, cells = keep_cells)
    seurat_sub <- UpdateSeuratObject(seurat_sub)
    coords     <- get_spatial_coords_fig4(seurat_sub); rm(seurat_sub)
    if (is.null(coords)) next

    keep_cells <- keep_cells[keep_cells %in% rownames(coords)]
    meta_keep  <- meta_keep[keep_cells, ]
    meta_keep$cell_type_grouped <- droplevels(as.factor(meta_keep$cell_type_grouped))

    # Build spatial CellChat to get cell-type x cell-type prob matrix
    cc_spatial <- tryCatch({
      cc <- createCellChat(
        object=full_mat_spatial[, keep_cells], meta=meta_keep,
        group.by="cell_type_grouped", datatype="spatial",
        coordinates=coords[keep_cells, ],
        spatial.factors=data.frame(ratio=0.5, tol=0.08)
      )
      cc@DB <- cellchat_db
      cc <- subsetData(cc)
      cc <- identifyOverExpressedGenes(cc)
      cc <- identifyOverExpressedInteractions(cc)
      cc <- computeCommunProb(cc, type="truncatedMean", trim=0.25,
                               population.size=TRUE, distance.use=FALSE,
                               contact.knn.k=10)
      cc <- filterCommunication(cc, min.cells=5)
      cc <- computeCommunProbPathway(cc)
      cc <- aggregateNet(cc)
      cc
    }, error=function(e) {
      cat(sprintf("      CellChat error: %s\n", e$message)); NULL
    })

    has_sema3 <- !is.null(cc_spatial) && "SEMA3" %in% cc_spatial@netP$pathways
    sema3_sum <- if (has_sema3)
      tryCatch(sum(cc_spatial@netP$prob[,,"SEMA3"]), error=function(e) 0) else 0
    cat(sprintf("      SEMA3 detected: %s | netP sum: %.2e\n", has_sema3, sema3_sum))

    p_full <- make_hexbin_spatial_plot(patient_id, grp, cc_spatial,
                                        coords[keep_cells,], meta_keep, has_sema3)

    pdf(out_full, width=10, height=8)
    print(p_full); dev.off()
    cat(sprintf("      Saved: %s\n", basename(out_full)))

    # Zoomed version: crop to the central 60% of the tissue bounding box,
    # centered on the tissue centroid — removes sparse edge cells while
    # keeping the main myocardial section where signaling occurs.
    xvals <- coords[keep_cells, "x"]; yvals <- coords[keep_cells, "y"]
    xrange <- range(xvals); yrange <- range(yvals)
    xpad <- diff(xrange) * 0.2;  ypad <- diff(yrange) * 0.2
    xlim_z <- c(xrange[1] + xpad, xrange[2] - xpad)
    ylim_z <- c(yrange[1] + ypad, yrange[2] - ypad)
    p_zoom <- p_full + coord_fixed(xlim=xlim_z, ylim=ylim_z)
    pdf(out_zoomed, width=10, height=8)
    print(p_zoom); dev.off()
    cat(sprintf("      Saved: %s\n", basename(out_zoomed)))

    if (!is.null(cc_spatial)) { rm(cc_spatial); gc() }

  }, error=function(e) cat(sprintf("      ERROR: %s\n", e$message)))
}

rm(full_mat_spatial); gc()

cat("\n=== PIPELINE v3 COMPLETE ===\n")
cat(sprintf("Results: %s\n\n", RESULTS_DIR))
cat(sprintf("  Per-patient: %s\n", PP_OUT_DIR))
cat(sprintf("  Differential: %s\n", DIFF_OUT_DIR))
cat(sprintf("  Final figures: %s\n", FIG_DIR))
cat("Files generated:\n")
cat("  Stage 3 tables/: master_all_comparisons_celltype_level.csv\n")
cat("                   supplementary_table_niche_communication.csv\n")
cat("                   NF_vs_pRV and pRV_vs_RVF Wilcoxon tables\n")
cat("  Stage 4 figures/: fig1_top15_pathways_myocardium.pdf\n")
cat("                            fig2_chord_myocardium_normalized.pdf\n")
cat("                            fig3_sema3c_circles.pdf\n")
cat("                            fig4_sema3c_hexbin_<patient>_<group>_full.pdf\n")
cat("                            fig4_sema3c_hexbin_<patient>_<group>_zoomed.pdf\n")
cat("                            fig5_sema3c_quantification.pdf\n")
