###############################################################################
## Supplementary Figure 7 (v54) — Pediatric myeloid and endothelial comparison
##
## Panels (per v53 manuscript caption):
##   (A) Myeloid cluster annotation concordance pediatric vs adult
##       (DotPlot of new 6-cluster adult-RV module scores in pediatric myeloid)
##   (B) Myeloid-specific bulk WGCNA modules (M1, M3, M4, M8) in pediatric
##   (C-D) M8 / MHCII component induction in ped-RVF (DotPlots, RVF vs NF & RVF vs pRV)
##   (E) ChEA: NR3C1 target induction and IRF8 target downregulation in pediatric
##   (F) Adult EC hdWGCNA modules (ecM1–ecM7) projected onto pediatric EC
##   (G) Pediatric EC subtype prevalence — venous-driven expansion in ped-RVF
##   (H) Angiogenesis-related GO scores in pediatric ECs
##   (I) MECOM expression: reduced in ped-HLHS/RVF arterial ECs
##   (J) Arterialization TF score across pediatric ECs
##   (K) Notch target score across pediatric ECs
##   (L) SMAD1 expression in pediatric capillary EC
##   (M) NR2F2 expression in pediatric venous EC
##
## Source code base: ~/Downloads/hdWGCNA_TOM/Manuscripts/Supplementary_Figure_7.R
## Key change vs source: ModuleScore gene lists rewritten to reflect the new
## 6-cluster RV myeloid taxonomy (per new_scripts/Figure_6.R) instead of the
## old 7-cluster taxonomy (CCR2+/CCR2- rMac1/rMac2 split).
##
## ── Statistical methodology disclosure ─────────────────────────────────────
##   Panels D and E (per-celltype × per-module ChEA enrichment dotplots) use
##   single-cell `FindMarkers` (Wilcoxon, p_val_adj < 0.05) on individual
##   cells, exactly matching v51 source layout. enrichr ChEA_2022 is run on
##   the up/down DEG lists with cutoff P.value < 0.05 (raw, not FDR), top-1
##   per (module × celltype). This treats each cell as an independent
##   observation, which inflates significance vs pseudobulk; some specific
##   IRF8/CIITA driver-gene enrichments do not retain FDR significance
##   under per-patient pseudobulk DESeq2 (n=4 ped-NF / n=5 ped-HLHS / n=5
##   ped-RVF). Panels are kept single-cell for visual consistency with the
##   originally published v51 figure; manuscript text reports the more
##   conservative pseudobulk interpretation for the IRF8 specific claim.
##
##   Panel "E2" (S7E_program_scores_box.pdf) — GR_homeostatic /
##   HIF_vascular / NFkB_MHCII / IFNg_AP module scores — IS plotted at
##   per-patient resolution (cell-level Seurat AddModuleScore aggregated
##   to per-patient means; ANOVA across 14 patients). Each dot = one
##   patient, so this panel is pseudobulk-style aggregation.
##
## Pediatric input objects (large; expensive to load):
##   ./dependencies/Figure_8/myeloid_annotated.rds       ~3.5 GB
##   ./dependencies/Figure_8/endothelium_annotated.rds   ~9.3 GB
## Output PDFs land in ./output/v52_figures/ with prefix S7*_
###############################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratWrappers)
  library(harmony)
  library(ggplot2)
  library(ggpubr)
  library(ggrepel)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(WGCNA)
  library(hdWGCNA)
  library(EnhancedVolcano)
  library(enrichR)
  library(forcats)
  library(viridis)
})

source('./helper_scripts/_shared_helpers.R')
## Working-repo adaptation: per-figure output dir.
V52_FIG_DIR <- './output/Supplementary_Figure_7'

COMP_W <- FINAL_WIDTH_IN
PS     <- pub_scales(COMP_W)

OUT <- V52_FIG_DIR
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

CACHE <- './output/S7_cache'
dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)

## peds object label convention: condition = {Donor, NF, SystolicHF}
## remap to v53 disease-state labels:
##   Donor      → ped-NF    (non-failing donor)
##   NF         → ped-HLHS  (compensated single-ventricle / HLHS)
##   SystolicHF → ped-RVF   (failing HLHS)
relabel_condition <- function(x) {
  x <- as.character(x)
  out <- rep(NA_character_, length(x))
  out[x == 'Donor']      <- 'ped-NF'
  out[x == 'NF']         <- 'ped-HLHS'
  out[x == 'SystolicHF'] <- 'ped-RVF'
  factor(out, levels = c('ped-NF','ped-HLHS','ped-RVF'))
}

###############################################################################
## ── Pediatric MYELOID — panels A-E ───────────────────────────────────────────
###############################################################################

cat('\n=== Loading pediatric myeloid object ===\n')
M1 <- readRDS('./dependencies/Figure_8/myeloid_annotated.rds')
DefaultAssay(M1) <- 'SCT'
M1$group <- relabel_condition(M1$condition)
Idents(M1) <- M1$sub.type

## ─── New adult-RV 6-cluster myeloid module scores ──────────────────────────
## (replaces v51 7-cluster CCR2+/CCR2-rMac1/rMac2 split)
## Marker genes harvested from new_scripts/Figure_6.R cluster annotations.
new_myeloid_sets <- list(
  ResMac          = c('LYVE1','FOLR2','VSIG4','MARCO','F13A1','MRC1',
                      'TIMD4','HPGDS','RBPJ','CD163L1'),
  iMac            = c('NR4A3','NR4A2','CD83','NR4A1','IL1B','CCL3',
                      'CCL4','CXCL3','CXCL8'),
  Mono            = c('FCN1','VCAN','CD300E','S100A8','S100A9','CCR2','SELL'),
  TREM2_Mac       = c('TREM2','GPNMB','SPP1','MITF','PLA2G7','LPL'),
  DC              = c('FLT3','CIITA','HLA-DRA','HLA-DQA1','CLEC10A','CLEC9A','CD1C'),
  Proliferating   = c('MKI67','TOP2A','CDK1')
)
new_myeloid_sets <- lapply(new_myeloid_sets, function(g) intersect(g, rownames(M1)))

for (nm in names(new_myeloid_sets)) {
  M1 <- AddModuleScore(M1, features = list(new_myeloid_sets[[nm]]),
                       assay = 'SCT', name = paste0(nm, '_Score'))
}

###############################################################################
## (A) Myeloid cluster annotation concordance pediatric vs adult
###############################################################################
.score_cols  <- paste0(names(new_myeloid_sets), '_Score1')
.score_label <- names(new_myeloid_sets)

Idents(M1) <- 'sub.type'
pdf(file.path(OUT, 'S7A_peds_myeloid_concordance_dot.pdf'),
    width = 7, height = 5)
print(
  DotPlot(M1, features = .score_cols, dot.min = 0,
          col.min = 0, col.max = 2) +
    RotatedAxis() +
    scale_x_discrete(labels = .score_label) +
    scale_color_gradient2(high = 'red', mid = 'grey95', low = 'blue') +
    xlab('New 6-cluster RV myeloid module score') + ylab('Pediatric subtype') +
    theme(panel.border = element_rect(linewidth = 1, fill = NA, color = 'black'))
)
dev.off()

## same view collapsed by disease group
pdf(file.path(OUT, 'S7A_peds_myeloid_concordance_dot_bygroup.pdf'),
    width = 6, height = 3)
print(
  DotPlot(M1, features = .score_cols, group.by = 'group',
          dot.min = 0, col.min = 0, col.max = 2) +
    RotatedAxis() +
    scale_x_discrete(labels = .score_label) +
    scale_color_gradient2(high = 'red', mid = 'grey95', low = 'blue') +
    xlab('New 6-cluster RV myeloid module score') + ylab('') +
    theme(panel.border = element_rect(linewidth = 1, fill = NA, color = 'black'))
)
dev.off()

## supporting volcanoes (RVF vs pRV, pRV vs NF) — kept compact
Idents(M1) <- 'group'
.cache_vol <- file.path(CACHE, 'S7_peds_myeloid_volcano_DEGs.rds')
if (file.exists(.cache_vol)) {
  .vol <- readRDS(.cache_vol)
} else {
  .vol <- list(
    RVF_pRV = FindMarkers(M1, ident.1 = 'ped-RVF', ident.2 = 'ped-HLHS'),
    pRV_NF  = FindMarkers(M1, ident.1 = 'ped-HLHS', ident.2 = 'ped-NF')
  )
  saveRDS(.vol, .cache_vol)
}
pdf(file.path(OUT, 'S7A_peds_myeloid_volcano.pdf'), width = 12, height = 5)
print(
  EnhancedVolcano(.vol$pRV_NF, lab = rownames(.vol$pRV_NF),
                  x = 'avg_log2FC', y = 'p_val_adj', FCcutoff = 0.1,
                  title = 'pRV vs NF (peds)') + coord_flip()
  +
  EnhancedVolcano(.vol$RVF_pRV, lab = rownames(.vol$RVF_pRV),
                  x = 'avg_log2FC', y = 'p_val_adj', FCcutoff = 0.1,
                  title = 'RVF vs pRV (peds)') + coord_flip()
)
dev.off()

###############################################################################
## (B) Pediatric expression of bulk WGCNA myeloid-specific modules M1/M3/M4/M8
###############################################################################
.cache_consensus <- file.path(CACHE, 'S7_consensus_modules.rds')
if (file.exists(.cache_consensus)) {
  consensus_modules <- readRDS(.cache_consensus)
} else {
  mapping <- WGCNA::labels2colors(1:100)
  cm <- read.csv('./dependencies/shared/bulk_heart_modules.csv',
                 stringsAsFactors = FALSE)[, 1:3]
  cm <- subset(cm, gene_name %in% rownames(M1))
  cm <- cm[match(unique(cm$gene_name), cm$gene_name), ]
  consensus_modules <- list(cm = cm, mapping = mapping)
  saveRDS(consensus_modules, .cache_consensus)
}

cm      <- consensus_modules$cm
mapping <- consensus_modules$mapping

score_calc      <- cm %>% dplyr::group_by(module) %>% dplyr::group_split()
module_colors   <- unique(unlist(lapply(score_calc, '[[', 'module')))
module_M_labels <- paste0('M', match(module_colors, mapping))

M1 <- AddModuleScore(M1, lapply(score_calc, '[[', 'gene_name'),
                     name = 'module_score')
.cols <- colnames(M1@meta.data)
.cols[startsWith(.cols, 'module_score')] <- paste0('module_', module_M_labels)
colnames(M1@meta.data) <- .cols

myeloid_M_panels <- c('M1','M3','M4','M8')
myeloid_M_avail  <- intersect(paste0('module_', myeloid_M_panels), colnames(M1@meta.data))

Idents(M1) <- 'sub.type'
pdf(file.path(OUT, 'S7B_peds_myeloid_modules_subtype.pdf'),
    width = 5, height = 3)
print(
  DotPlot(M1, features = myeloid_M_avail,
          dot.min = 0, col.min = 0, col.max = 2) +
    RotatedAxis() + xlab('') + ylab('') +
    scale_x_discrete(labels = sub('module_', '', myeloid_M_avail)) +
    scale_color_gradient2(high = 'red', mid = 'grey95', low = 'blue') +
    theme(panel.border = element_rect(linewidth = 1, fill = NA, color = 'black'))
)
dev.off()

Idents(M1) <- 'group'
Idents(M1) <- factor(Idents(M1), levels = c('ped-NF','ped-HLHS','ped-RVF'))
pdf(file.path(OUT, 'S7B_peds_myeloid_modules_group.pdf'),
    width = 5, height = 2.5)
print(
  DotPlot(M1, features = myeloid_M_avail,
          dot.min = 0, col.min = 0, col.max = 2) +
    RotatedAxis() + xlab('') + ylab('') +
    scale_x_discrete(labels = sub('module_', '', myeloid_M_avail)) +
    scale_color_gradient2(high = 'red', mid = 'grey95', low = 'blue') +
    theme(panel.border = element_rect(linewidth = 1, fill = NA, color = 'black'))
)
dev.off()

## S7B composite (v57: "by cell type AND disease state") — subtype dotplot
## (top) stacked over disease-group dotplot (bottom).  Rebuilt as objects.
.s7b_dot <- function(idv, lv) {
  Idents(M1) <- idv
  if (!is.null(lv)) Idents(M1) <- factor(Idents(M1), levels = lv)
  DotPlot(M1, features = myeloid_M_avail, dot.min = 0,
          col.min = 0, col.max = 2) +
    RotatedAxis() + xlab('') + ylab('') +
    scale_x_discrete(labels = sub('module_', '', myeloid_M_avail)) +
    scale_color_gradient2(high = 'red', mid = 'grey95', low = 'blue') +
    theme(panel.border = element_rect(linewidth = 1, fill = NA, color = 'black'))
}
pdf(file.path(OUT, 'S7B_peds_myeloid_modules_combined.pdf'),
    width = 5, height = 5.6)
print(.s7b_dot('sub.type', NULL) /
      .s7b_dot('group', c('ped-NF','ped-HLHS','ped-RVF')) +
      plot_layout(heights = c(3, 2.5)))
dev.off()
cat('Wrote: S7B_peds_myeloid_modules_combined.pdf\n')

###############################################################################
## (C-D) Per-celltype × per-module ChEA enrichment dotplots — pediatric myeloid
##
## Verbatim port of v51 Supplementary_Figure_7.R lines 540-720 with two fixes:
##   - Reactome_2016 → Reactome_2022 (older DB returns 0 terms here)
##   - Per-celltype × per-module FindMarkers loop, NOT a single global rbind
##     (avoids a row-name prefix bug that would break enrichr lookups)
##
## For each module i ∈ {M1, M3, M4, M8}, peds sub.type j, and comparison
## k ∈ {(RVF, NF), (RVF, pRV)}:
##   FindMarkers within the bulk M_i gene set →
##   enrichr(up genes) and enrichr(down genes) →
##   accumulate into combined_output, then plot ChEA UP terms as a dotplot
##   coloured by -log10(P-adj) and sized by log(Combined.Score), with a
##   bottom colour bar showing module identity per column.
##
## Result is cached under output/S7_cache/S7_celltype_module_enrichr.rds —
## subsequent runs skip the (slow) FindMarkers + enrichr loop.
###############################################################################
mapping_wgcna <- WGCNA::labels2colors(1:100)
mod_color_lookup <- c(M1 = mapping_wgcna[1], M3 = mapping_wgcna[3],
                      M4 = mapping_wgcna[4], M8 = mapping_wgcna[8])
dbs <- c('ChEA_2022','Reactome_2022','GO_Biological_Process_2023')

## _rawp suffix: invalidates the stale FDR-filtered cache so combined_output
## is rebuilt with the corrected raw-P.value DEG selection (fixes blank S7E).
.cache_celltype_enr <- file.path(CACHE, 'S7_celltype_module_enrichr_adjp.rds')
if (file.exists(.cache_celltype_enr)) {
  combined_output <- readRDS(.cache_celltype_enr)
} else {
  M1$Names_group <- paste0(as.character(M1$sub.type), '_', as.character(M1$group))
  Idents(M1)     <- 'Names_group'

  bulk_modules <- read.csv('./dependencies/shared/bulk_heart_modules.csv',
                           stringsAsFactors = FALSE)
  bulk_modules$module_idx <- match(bulk_modules$module, mapping_wgcna)

  mods_idx   <- c(1, 3, 4, 8)
  cell_types <- unique(as.character(M1$sub.type))
  ## M1$group is relabel_condition() output (ped-NF/ped-HLHS/ped-RVF), so the
  ## FindMarkers idents (Names_group) carry ped-* values.  Match on those but
  ## keep the short stored label (RVF_NF / pRV_NF / RVF_pRV) that the panels
  ## filter on.  Previously this used raw RVF/NF/pRV -> id1 never matched ->
  ## combined_output stayed EMPTY -> panels C/D/E were all blank.  Also add
  ## the pRV_NF (ped-HLHS vs ped-NF) pair that Panel D requires.
  comparison <- list(
    list(g = c('ped-RVF','ped-NF'),   lab = 'RVF_NF'),
    list(g = c('ped-HLHS','ped-NF'),  lab = 'pRV_NF'),
    list(g = c('ped-RVF','ped-HLHS'), lab = 'RVF_pRV'))

  combined_output <- data.frame(
    Term = character(), Overlap = character(), P.value = numeric(),
    Adjusted.P.value = numeric(), Combined.Score = numeric(),
    db = character(), module = character(), celltype = character(),
    comparison = character(), color = character(), direction = character(),
    stringsAsFactors = FALSE)

  for (i in mods_idx) {
    key_genes <- intersect(subset(bulk_modules, module_idx == i)$gene_name,
                           rownames(M1))
    for (j in cell_types) {
      for (k in comparison) {
        id1 <- paste0(j, '_', k$g[1]); id2 <- paste0(j, '_', k$g[2])
        if (!(id1 %in% Idents(M1)) || !(id2 %in% Idents(M1))) next
        gs <- tryCatch(FindMarkers(M1, ident.1 = id1, ident.2 = id2,
                                   features = key_genes, verbose = FALSE),
                       error = function(e) NULL)
        if (is.null(gs) || nrow(gs) == 0) next
        ## Legacy DEG filter (Manuscripts/Supplementary_Figure_7.R:141):
        ## FDR p_val_adj < 0.05.  (The earlier raw-p workaround was added
        ## to un-blank Panel E, but the true blank cause was the ped-*
        ## ident mismatch — now fixed — so restore the legacy threshold
        ## for an exact match to the published TF set.)
        gs <- subset(gs, p_val_adj < 0.05)
        if (nrow(gs) == 0) next
        for (dir_label in c('up','down')) {
          syms <- rownames(if (dir_label == 'up') subset(gs, avg_log2FC > 0)
                                                  else subset(gs, avg_log2FC < 0))
          if (length(syms) < 3) next
          enr <- tryCatch(enrichR::enrichr(syms, dbs),
                          error = function(e) NULL)
          Sys.sleep(2)
          if (is.null(enr)) next
          for (db in names(enr)) {
            cur <- enr[[db]]
            if (is.null(cur) || nrow(cur) < 1) next
            cur$db         <- db
            cur$module     <- paste0('M', i)
            cur$celltype   <- j
            cur$comparison <- k$lab
            cur$color      <- mapping_wgcna[i]
            cur$direction  <- dir_label
            combined_output <- rbind(combined_output, cur)
          }
        }
      }
    }
  }
  saveRDS(combined_output, .cache_celltype_enr)
}

## Selection logic ports v51 Supplementary_Figure_7.R lines 1300-1395:
##   filter by raw P.value < 0.05, take TOP-1 per (module × celltype),
##   colour by -log(P.value) clipped at max_p = 10. (Adjusted.P.value
##   is intentionally NOT used here — it over-filters and drops most
##   modules/celltypes, contradicting the published figure.)
build_chea_dotplot <- function(co, comparison_label, direction_label, title_lab,
                               modules_keep = c('M1','M3','M4','M8'),
                               pval_cut = 0.05, max_p = 10,
                               tf_keep = NULL, mc_order = NULL) {
  sel <- subset(co, db == 'ChEA_2022' & direction == direction_label &
                comparison == comparison_label & P.value < pval_cut &
                module %in% modules_keep)
  if (nrow(sel) == 0)
    return(list(plot = ggplot() + theme_void() +
                  ggtitle(paste(title_lab, '- no terms')),
                colorbar = ggplot() + theme_void(), n_terms = 0))
  sel$module_celltype <- paste0(sel$module, '_', sel$celltype)
  sel$wrap <- sub(' .*$', '', sel$Term)
  sel$logp <- -log(sel$P.value)
  sel$logp <- ifelse(sel$logp > max_p, max_p, sel$logp)

  if (!is.null(tf_keep) || !is.null(mc_order)) {
    ## Curated pixel-match mode: keep exactly the requested TFs and
    ## module_celltype columns, in the requested order; one dot per
    ## (TF, module_celltype) = best (lowest) P.value.
    if (!is.null(tf_keep))  sel <- sel[sel$wrap %in% tf_keep, , drop = FALSE]
    if (!is.null(mc_order)) sel <- sel[sel$module_celltype %in% mc_order, ,
                                       drop = FALSE]
    if (nrow(sel) == 0)
      return(list(plot = ggplot() + theme_void() +
                    ggtitle(paste(title_lab, '- no terms')),
                  colorbar = ggplot() + theme_void(), n_terms = 0))
    sel <- sel[order(sel$P.value), ]
    sel <- sel[!duplicated(paste(sel$wrap, sel$module_celltype)), ]
    .mcl <- if (!is.null(mc_order)) mc_order else unique(sel$module_celltype)
    .tfl <- if (!is.null(tf_keep))  tf_keep  else unique(sel$wrap)
    sel$module          <- factor(sel$module, levels = modules_keep)
    sel$module_celltype <- factor(sel$module_celltype, levels = .mcl)
    sel$wrap            <- factor(sel$wrap, levels = rev(.tfl))
    .drop <- FALSE
  } else {
    ## v51 default: TOP-1 per module_celltype.  v51 RELIED on enrichr
    ## returning rows pre-sorted by Combined.Score desc, but the rebuilt
    ## combined_output stores enrichR output in P.value order — so "first
    ## per module_celltype" was picking the min-P TF (a mix of immune +
    ## generic TFs) instead of the max-Combined.Score TF.  Explicitly sort
    ## by Combined.Score desc first so the biologically meaningful
    ## regulators (CIITA, IRF8, NR3C1, RELB, MYB, ...) surface naturally.
    sel <- sel[order(-sel$Combined.Score), ]
    idx_top_1 <- match(unique(sel$module_celltype), sel$module_celltype)
    sel <- sel[idx_top_1, ]
    sel$module <- factor(sel$module, levels = modules_keep)
    sel <- sel %>% arrange(module, celltype)
    sel$module_celltype <- factor(sel$module_celltype,
                                  levels = unique(sel$module_celltype))
    sel$wrap <- factor(sel$wrap, levels = rev(unique(sel$wrap)))
    .mcl  <- levels(sel$module_celltype)
    .drop <- TRUE
  }

  p <- ggplot(sel, aes(x = module_celltype, y = wrap,
                       color = logp, size = log(Combined.Score))) +
    geom_point() +
    scale_x_discrete(drop = .drop) +
    scale_color_stepsn(colors = rev(viridis::magma(256)),
                       name = 'logP') +
    scale_size_continuous(name = 'log(Score)') +
    RotatedAxis() + xlab('') + ylab('') + ggtitle(title_lab) +
    theme(panel.border = element_rect(linewidth = 1, color = 'black', fill = NA),
          axis.text.x  = element_blank(), axis.ticks.x = element_blank(),
          panel.grid   = element_line(linewidth = 0.25, color = 'lightgrey'),
          plot.margin  = margin(0, 0, 0, 0))

  .cbdf <- data.frame(module_celltype = factor(.mcl, levels = .mcl),
                       stringsAsFactors = FALSE)
  .cbdf$module <- sub('_.*$', '', as.character(.cbdf$module_celltype))
  .cbdf$colour <- mod_color_lookup[.cbdf$module]
  .cbdf$y      <- 1
  colorbar <- ggplot(.cbdf, aes(x = module_celltype, y = y,
                                fill = module_celltype)) +
    geom_tile() +
    scale_x_discrete(drop = .drop) +
    scale_fill_manual(values = setNames(.cbdf$colour,
                                        as.character(.cbdf$module_celltype))) +
    NoLegend() + RotatedAxis() +
    theme(plot.title = element_blank(), axis.line = element_blank(),
          axis.ticks.y = element_blank(), axis.text.y = element_blank(),
          axis.title = element_blank(), plot.margin = margin(0, 0, 0, 0))
  list(plot = p, colorbar = colorbar, n_terms = nrow(sel))
}

## Panel D — ped-HLHS vs ped-NF Up — all 4 modules.
## Comparison = pRV_NF (NF condition vs Donor) under Option A label mapping:
##   ped-NF   = Donor       (true non-failing biventricular)
##   ped-HLHS = NF condition (compensated single-ventricle baseline)
##   ped-RVF  = SystolicHF  (failed single-ventricle)
## combined_output is now correctly populated (ped-* label fix), so the
## DEFAULT top-1-TF-per-module_celltype logic reproduces the published
## diagonal staircase.  (Hardcoded tf_keep/mc_order produced nonsense.)
panel_D <- build_chea_dotplot(combined_output, 'pRV_NF', 'up',
                              'ped-HLHS vs ped-NF Up',
                              modules_keep = c('M1','M3','M4','M8'))
pdf(file.path(OUT, 'S7D_chea_celltype_HLHS_vs_NF_up.pdf'),
    width = 10, height = 5)
print(panel_D$plot / panel_D$colorbar + plot_layout(heights = c(8, 1)))
dev.off()

## MHCII component DotPlot — anchors the C/D narrative as a complementary view.
.mhc_genes <- intersect(c('CIITA','HLA-DRA','HLA-DRB1','HLA-DPA1','HLA-DPB1',
                          'HLA-DQA1','HLA-DQB1','CD74'),
                        rownames(M1))
Idents(M1) <- 'group'
Idents(M1) <- factor(Idents(M1), levels = c('ped-NF','ped-HLHS','ped-RVF'))
pdf(file.path(OUT, 'S7CD_peds_myeloid_MHCII_dot.pdf'),
    width = 6, height = 2.5)
print(
  DotPlot(M1, features = .mhc_genes, dot.min = 0, col.min = 0, col.max = 2) +
    RotatedAxis() +
    scale_color_gradient2(high = 'red', mid = 'grey95', low = 'blue') +
    xlab('MHC class II components') + ylab('') +
    theme(panel.border = element_rect(linewidth = 1, fill = NA, color = 'black'))
)
dev.off()

###############################################################################
## Panel S7C (v57) — GO-BP enrichment of UPregulated M8 DEGs for failing
## single-ventricle (ped-RVF) vs NF bi-ventricle (ped-NF). Built from the
## same combined_output (db = GO_Biological_Process_2023); raw P.value
## selection, ranked by Combined.Score, GO-id suffix stripped. This panel
## had no code previously (emission map flagged it MISSING).
###############################################################################
.s7c <- subset(combined_output,
               db == 'GO_Biological_Process_2023' &
               comparison == 'RVF_NF' & direction == 'up' &
               module == 'M8' & P.value < 0.05)
if (nrow(.s7c) > 0) {
  .s7c <- .s7c[order(.s7c$Combined.Score, decreasing = TRUE), ]
  .s7c <- .s7c[!duplicated(.s7c$Term), ]
  .s7c <- head(.s7c, 10)
  .s7c$Term  <- sub(' \\(GO:[0-9]+\\)', '', .s7c$Term)
  .s7c$Term  <- factor(.s7c$Term, levels = rev(.s7c$Term))
  .s7c$logp  <- pmin(-log10(.s7c$P.value), 10)
  p_S7C <- ggplot(.s7c, aes(x = Combined.Score, y = Term,
                            color = logp, size = log(Combined.Score))) +
    geom_point() +
    scale_color_gradient(low = 'grey80', high = '#B2182B',
                         name = expression(-log[10] ~ P)) +
    scale_size(range = c(2, 7), name = 'log(Comb. Score)') +
    xlab('Combined Score') + ylab('') +
    ggtitle('M8 up: failing-SV vs NF-biventricle (GO-BP)') +
    theme_classic() +
    theme(axis.text = element_text(colour = 'black'),
          panel.border = element_rect(linewidth = 0.6, fill = NA,
                                      color = 'black'),
          plot.title = element_text(size = 9, face = 'bold'))
} else {
  message('Panel S7C: no GO-BP rows for M8 RVF_NF up — placeholder emitted')
  p_S7C <- ggplot() +
    annotate('text', x = 0.5, y = 0.5, size = 4,
             label = 'S7C: no enriched GO-BP terms\n(M8 up, failing-SV vs NF-biventricle)') +
    theme_void()
}
pdf(file.path(OUT, 'S7C_peds_myeloid_M8_GO_up.pdf'), width = 10, height = 4)
print(p_S7C)
dev.off()
cat('Wrote: S7C_peds_myeloid_M8_GO_up.pdf\n')

###############################################################################
## Panel E — ped-RVF vs ped-NF Up + Down — restricted to M1 and M8 only.
## Comparison = RVF_NF (SystolicHF vs Donor):
##   E up : note CIITA enrichment in M8, NR3C1/RELB/MYB/FOXM1 in M1.
##   E dn : note IRF8 (and PPARG, CIITA, NR1H3, ESR1) enrichment.
###############################################################################
panel_E_up <- build_chea_dotplot(combined_output, 'RVF_NF', 'up',
                                 'ped-RVF vs ped-NF Up',
                                 modules_keep = c('M1','M8'))
pdf(file.path(OUT, 'S7E_chea_celltype_RVF_vs_NF_up.pdf'),
    width = 7, height = 4)
print(panel_E_up$plot / panel_E_up$colorbar + plot_layout(heights = c(8, 1)))
dev.off()

panel_E_dn <- build_chea_dotplot(combined_output, 'RVF_NF', 'down',
                                 'ped-RVF vs ped-NF Down',
                                 modules_keep = c('M1','M8'))
pdf(file.path(OUT, 'S7E_chea_celltype_RVF_vs_NF_down.pdf'),
    width = 7, height = 4)
print(panel_E_dn$plot / panel_E_dn$colorbar + plot_layout(heights = c(8, 1)))
dev.off()

## S7E (v57): ChEA enriched TFs for ped-RVF (failing single-ventricle) vs
## ped-NF (NF bi-ventricle), M1+M8 genes only — UPregulated DEGs (left,
## expect NR3C1/CIITA) beside DOWNregulated DEGs (right, expect IRF8).
.s7e_L <- panel_E_up$plot / panel_E_up$colorbar + plot_layout(heights = c(8, 1))
.s7e_R <- panel_E_dn$plot / panel_E_dn$colorbar + plot_layout(heights = c(8, 1))
pdf(file.path(OUT, 'S7E_chea_RVF_vs_NF_up_down.pdf'), width = 13, height = 4.5)
print(.s7e_L | .s7e_R)
dev.off()
cat('Wrote: S7E_chea_RVF_vs_NF_up_down.pdf\n')

###############################################################################
## (E2) Curated myeloid program module scores — boxplots per disease state
##      Gene lists pulled from new_scripts/Figure_6.R lines 915-925:
##        GR_homeostatic, HIF_vascular, NFkB_MHCII, IFNg_AP
##      Rendered as a horizontally-concatenated panel of box+whisker plots
##      with NF/pRV/RVF on x-axis, jittered points overlaid, ggpubr ANOVA.
###############################################################################
GR_homeostatic <- c('TSC22D3','FKBP5','KLF9','KLF13','MERTK','CD163','MARCO',
                    'MT2A','GLUL','ZBTB16','BCL6','CEBPB','TFCP2L1','MAFB',
                    'GADD45B','ANGPTL4')
HIF_vascular   <- c('EPAS1','RAMP3','NAMPT','TIMP3','PECAM1','MYH9','LDLR',
                    'PODXL','ID1','MAF','CHD7','ADIPOR2')
NFkB_MHCII     <- c('CD74','HLA-DRB1','HLA-DRA','HLA-DPA1','HLA-DPB1',
                    'NLRP1','NLRP3','AIM2','CIITA','RUNX1','FOSL2','SRGN',
                    'FOXO3','IRF1','CD83')
IFNg_AP        <- c('STAT1','IRF1','IRF8','JAK2',
                    'B2M','HLA-A','HLA-B','HLA-C',
                    'PSMB8','PSMB9','TAP1','IFI30',
                    'GBP1','GBP2','GBP4','GBP5','ICAM1')

prog_list <- list(
  GR_homeostatic = intersect(GR_homeostatic, rownames(M1)),
  HIF_vascular   = intersect(HIF_vascular,   rownames(M1)),
  NFkB_MHCII     = intersect(NFkB_MHCII,     rownames(M1)),
  IFNg_AP        = intersect(IFNg_AP,        rownames(M1))
)
prog_label <- c(GR_homeostatic = 'GR homeostatic',
                HIF_vascular   = 'HIF / vascular drift',
                NFkB_MHCII     = 'NF-kB MHCII / inflammasome',
                IFNg_AP        = 'IFNg antigen presentation')

for (nm in names(prog_list)) {
  M1 <- AddModuleScore(M1, list(prog_list[[nm]]),
                       assay = 'SCT', name = paste0(nm, '_'))
}

prog_score_cols <- paste0(names(prog_list), '_1')

## Aggregate to per-patient (sample) means so each patient is one dot.
## Disease labels: ped-NF (Donor), ped-HLHS (NF), ped-RVF (SystolicHF).
## M1$group is already relabel_condition() output (ped-NF/ped-HLHS/ped-RVF).
## The old NF/pRV/RVF -> ped-* remap produced all-NA (1-group / blank S7F).
prog_cell_df <- data.frame(
  patient = M1$sample,
  group   = as.character(M1$group),
  M1@meta.data[, prog_score_cols, drop = FALSE],
  check.names = FALSE
)
prog_pat_df <- prog_cell_df %>%
  dplyr::group_by(patient, group) %>%
  dplyr::summarise(across(all_of(prog_score_cols), mean), .groups = 'drop')

prog_long_pat <- prog_pat_df %>%
  tidyr::pivot_longer(all_of(prog_score_cols),
                      names_to = 'program', values_to = 'score') %>%
  dplyr::mutate(program = sub('_1$', '', program),
                program = factor(program,
                                 levels = names(prog_label),
                                 labels = prog_label),
                group   = factor(group,
                                 levels = c('ped-NF','ped-HLHS','ped-RVF')))

## ggpubr ANOVA expects matching colour scale; recolour with disease_pal keys
ped_pal <- setNames(unname(disease_pal[c('NF','pRV','RVF')]),
                    c('ped-NF','ped-HLHS','ped-RVF'))

## Reverse factor levels so the topmost row after coord_flip() is ped-NF
prog_long_pat$group <- factor(prog_long_pat$group,
                              levels = rev(c('ped-NF','ped-HLHS','ped-RVF')))

p_box <- ggplot(prog_long_pat, aes(x = group, y = score, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.85, linewidth = 0.4) +
  geom_jitter(width = 0.18, size = 2.4, shape = 21,
              colour = 'black', stroke = 0.3, alpha = 0.9) +
  facet_wrap(~ program, nrow = 1, scales = 'free_x') +
  coord_flip() +
  scale_fill_manual(values = ped_pal) +
  ggpubr::stat_compare_means(method = 'anova',
                             label = 'p.format', size = 2.6) +
  theme_classic() +
  xlab('') + ylab('Module score (per-patient mean)') +
  theme(strip.text   = element_text(face = 'bold'),
        legend.position = 'none',
        panel.border = element_rect(linewidth = 0.6, fill = NA, color = 'black'))

pdf(file.path(OUT, 'S7E_program_scores_box.pdf'),
    width = 11, height = 2.38)
print(p_box)
dev.off()

rm(M1); invisible(gc())

###############################################################################
## ── Pediatric ENDOTHELIUM — panels F-M ──────────────────────────────────────
###############################################################################

cat('\n=== Loading pediatric endothelial object ===\n')
M1 <- readRDS('./dependencies/Figure_8/endothelium_annotated.rds')
DefaultAssay(M1) <- 'SCT'
M1$group <- relabel_condition(M1$condition)

## Keep peds-native EC subtype labels (Artery1, Artery2, Capillary, Kit+,
## Vein, Endocardium, Lymphatic, ...) — do NOT canonicalise to adult names.
M1$ec_subtype <- as.character(M1$sub.type)
Idents(M1) <- 'ec_subtype'

## Local disease palette keyed on ped-NF / ped-HLHS / ped-RVF
## (shadows the global disease_pal which uses NF/pRV/RVF keys)
disease_pal <- c('ped-NF'   = unname(disease_pal[['NF']]),
                 'ped-HLHS' = unname(disease_pal[['pRV']]),
                 'ped-RVF'  = unname(disease_pal[['RVF']]))

## ── Per-patient pseudobulk helpers ──────────────────────────────────────────
## Aggregate cell-level scores or expression to per-patient means and plot a
## ggboxplot + ANOVA across n≈14 patients (mirrors S7E_program_scores_box.pdf).
.peds_pb_box <- function(obj, score_cols, label_map = NULL,
                         cell_filter = NULL, out_path,
                         width = 10, height = 2.8,
                         ylab = 'Module score (per-patient mean)',
                         patient_col = 'sample',
                         group_levels = c('ped-NF','ped-HLHS','ped-RVF'),
                         palette = NULL,
                         data_cache_path = NULL) {
  if (is.null(palette)) palette <- disease_pal
  meta <- obj@meta.data
  if (!is.null(cell_filter)) meta <- meta[cell_filter, , drop = FALSE]
  cell_df <- data.frame(
    patient = meta[[patient_col]],
    group   = as.character(meta$group),
    meta[, score_cols, drop = FALSE],
    check.names = FALSE
  )
  pat_df <- cell_df %>%
    dplyr::group_by(patient, group) %>%
    dplyr::summarise(dplyr::across(dplyr::all_of(score_cols), mean),
                     .groups = 'drop')
  long_pat <- pat_df %>%
    tidyr::pivot_longer(dplyr::all_of(score_cols),
                        names_to = 'program', values_to = 'score') %>%
    dplyr::mutate(program = sub('_1$', '', program))
  if (!is.null(label_map)) {
    long_pat$program <- factor(long_pat$program,
                               levels = names(label_map),
                               labels = label_map)
  }
  long_pat$group <- factor(long_pat$group, levels = group_levels)
  if (!is.null(data_cache_path)) saveRDS(long_pat, data_cache_path)
  p <- ggplot(long_pat, aes(x = group, y = score, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.85, linewidth = 0.4) +
    geom_jitter(width = 0.18, size = 2.4, shape = 21,
                colour = 'black', stroke = 0.3, alpha = 0.9) +
    facet_wrap(~ program, nrow = 1, scales = 'free_y') +
    scale_fill_manual(values = palette) +
    ggpubr::stat_compare_means(method = 'anova',
                               label = 'p.format', size = 2.6) +
    theme_classic() +
    xlab('') + ylab(ylab) +
    theme(strip.text      = element_text(face = 'bold'),
          legend.position = 'none',
          panel.border    = element_rect(linewidth = 0.6, fill = NA, color = 'black'))
  pdf(out_path, width = width, height = height); print(p); dev.off()
}

.peds_pb_gene_box <- function(obj, gene, cell_filter = NULL, out_path,
                              width = 3.4, height = 3, title = NULL,
                              slot = 'data',
                              patient_col = 'sample',
                              group_levels = c('ped-NF','ped-HLHS','ped-RVF'),
                              palette = NULL,
                              data_cache_path = NULL) {
  if (is.null(palette)) palette <- disease_pal
  expr <- Seurat::FetchData(obj, vars = gene, layer = slot)
  meta <- obj@meta.data
  if (!is.null(cell_filter)) {
    expr <- expr[cell_filter, , drop = FALSE]
    meta <- meta[cell_filter, , drop = FALSE]
  }
  df <- data.frame(patient = meta[[patient_col]],
                   group   = as.character(meta$group),
                   expr    = expr[, gene])
  pat_df <- df %>%
    dplyr::group_by(patient, group) %>%
    dplyr::summarise(expr = mean(expr), .groups = 'drop') %>%
    dplyr::mutate(group = factor(group, levels = group_levels))
  if (!is.null(data_cache_path)) saveRDS(pat_df, data_cache_path)
  p <- ggplot(pat_df, aes(x = group, y = expr, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.85, linewidth = 0.4) +
    geom_jitter(width = 0.18, size = 2.4, shape = 21,
                colour = 'black', stroke = 0.3, alpha = 0.9) +
    scale_fill_manual(values = palette) +
    ggpubr::stat_compare_means(method = 'anova',
                               label = 'p.format', size = 2.6) +
    theme_classic() +
    xlab('') + ylab(paste0(gene, ' expression (per-patient mean)')) +
    theme(legend.position = 'none',
          panel.border = element_rect(linewidth = 0.6, fill = NA, color = 'black'))
  if (!is.null(title)) p <- p + ggtitle(title)
  pdf(out_path, width = width, height = height); print(p); dev.off()
}

## ── Fig 7 Phase-2 framework: gene categories + DESeq2 helpers ────────────────
## Mirrors Figure_7.R lines 1019-1037. Used by panels H (volcano), I
## (peds-vs-adult sn concordance), J (NRG1 endo), K (MECOM arterial).
phase2_cats <- list(
  Mesenchymal   = c('SNAI1'),
  Angiogenic    = c('ANGPT2','MYC','ADAMTS1'),
  EC_activation = c('PLA2G2A','ACKR1','PECAM1','VWF'),
  EC_stress     = c('STC1','NR4A1','NR4A3'),
  Identity      = c('MECOM','CA4')
)
phase2_pal <- c(Mesenchymal   = '#B2182B', Angiogenic    = '#D6604D',
                EC_activation = '#2166AC', EC_stress     = '#1B9E77',
                Identity      = '#762A83', Other         = 'grey80')

## ── Phase-1 EC framework: vasoprotective program erosion ────────────────────
## Empirically anchored to sn + Xenium pRV_vs_NF EC pseudobulk (Phase 1).
## ACKR1 included for biphasic display; excluded from vasoprot composite score.
phase1_cats <- list(
  Antioxidant       = c('GPX3','TXNRD1','HMOX1','NQO1'),
  RAAS_vasoreactive = c('ACE2','EDN1'),
  Anti_inflammatory = c('IL1RL1','CD163','VSIG4','DNASE1L3','F13A1','AQP3','ACKR1'),
  ECM_stabilization = c('TIMP3','TIMP4'),
  EC_quiescence     = c('GATA2','BMP6','TMEM100')
)
phase1_pal <- c(Antioxidant       = '#FF7F00',
                RAAS_vasoreactive = '#1F78B4',
                Anti_inflammatory = '#33A02C',
                ECM_stabilization = '#E31A1C',
                EC_quiescence     = '#6A3D9A',
                Other             = 'grey80')
phase1_lookup <- unlist(lapply(names(phase1_cats), function(k)
  setNames(rep(k, length(phase1_cats[[k]])), phase1_cats[[k]])))
phase1_genes  <- names(phase1_lookup)

phase2_lookup <- unlist(lapply(names(phase2_cats), function(k)
  setNames(rep(k, length(phase2_cats[[k]])), phase2_cats[[k]])))
phase2_genes  <- names(phase2_lookup)

## Pan-EC pseudobulk DESeq2 (lfcShrink ashr) — returns df with gene, log2FoldChange, padj
.ec_pan_pseudobulk_de <- function(obj, contrast, patient_col, cache_path) {
  if (!is.null(cache_path) && file.exists(cache_path)) {
    return(readRDS(cache_path))
  }
  suppressPackageStartupMessages({ library(DESeq2) })
  DefaultAssay(obj) <- 'RNA'
  agg <- Seurat::AggregateExpression(obj, assays = 'RNA', slot = 'counts',
                                     group.by = c(patient_col, 'group'),
                                     return.seurat = FALSE)$RNA
  cn <- colnames(agg)
  cd <- data.frame(sample  = cn,
                   patient = sub('_[^_]+$', '', cn),
                   group   = sub('.*_', '', cn),
                   stringsAsFactors = FALSE)
  cd$group <- factor(cd$group, levels = unique(c(contrast[3], contrast[2],
                                                  setdiff(unique(cd$group),
                                                          c(contrast[2], contrast[3])))))
  dds <- DESeq2::DESeqDataSetFromMatrix(round(as.matrix(agg)),
                                        colData = cd, design = ~ group)
  dds <- DESeq2::DESeq(dds)
  res <- DESeq2::lfcShrink(dds, contrast = c('group', contrast[2], contrast[3]),
                            type = 'ashr')
  out <- as.data.frame(res); out$gene <- rownames(out)
  if (!is.null(cache_path)) saveRDS(out, cache_path)
  rm(agg, dds, res); gc(verbose = FALSE)
  out
}

## Subtype pseudobulk DESeq2-VST → per-donor expression of one gene
## (mirrors Figure_7.R lines 691-709 but parameterized for subtype-column +
##  patient-column + contrast levels).
.ec_subtype_pseudobulk_vst <- function(obj, subtype_col, subtype_values, gene,
                                        patient_col, group_levels, cache_path) {
  if (!is.null(cache_path) && file.exists(cache_path)) {
    return(readRDS(cache_path))
  }
  suppressPackageStartupMessages({ library(DESeq2); library(SummarizedExperiment) })
  keep <- as.character(obj@meta.data[[subtype_col]]) %in% subtype_values
  if (sum(keep) == 0) {
    cat('No cells matched subtype_values=', paste(subtype_values, collapse=','),
        ' in column ', subtype_col,
        '. Available values: ',
        paste(unique(as.character(obj@meta.data[[subtype_col]])), collapse=', '),
        '\n', sep = '')
    return(NULL)
  }
  sub  <- subset(obj, cells = colnames(obj)[keep])
  DefaultAssay(sub) <- 'RNA'
  agg <- Seurat::AggregateExpression(sub, assays = 'RNA', slot = 'counts',
                                     group.by = c(patient_col, 'group'),
                                     return.seurat = FALSE)$RNA
  cn <- colnames(agg)
  cd <- data.frame(sample  = cn,
                   patient = sub('_[^_]+$', '', cn),
                   group   = sub('.*_', '', cn),
                   stringsAsFactors = FALSE)
  cd$group <- factor(cd$group, levels = group_levels)
  dds <- DESeq2::DESeqDataSetFromMatrix(round(as.matrix(agg)),
                                        colData = cd, design = ~ group)
  vsd <- DESeq2::vst(dds, blind = FALSE)
  vmat <- as.matrix(SummarizedExperiment::assay(vsd))
  if (!gene %in% rownames(vmat)) return(NULL)
  out <- data.frame(sample = cd$sample, patient = cd$patient,
                    group  = cd$group,  expr = vmat[gene, ])
  if (!is.null(cache_path)) saveRDS(out, cache_path)
  out
}

## Volcano with Phase-2 genes coloured by program (Figure_7.R-style, coord-flipped)
.ec_volcano <- function(df, label_text, fc_xlim, sig_ylim, show_legend = TRUE,
                         subtitle_text = 'EC pseudobulk (Phase 2)') {
  suppressPackageStartupMessages(library(EnhancedVolcano))
  df <- df[!is.na(df$padj) & !is.na(df$log2FoldChange), ]
  rownames(df) <- make.unique(as.character(df$gene))
  df$category  <- ifelse(df$gene %in% phase2_genes,
                         phase2_lookup[df$gene], 'Other')
  keyvals <- setNames(phase2_pal[df$category], df$category)
  sel     <- intersect(phase2_genes, df$gene)
  EnhancedVolcano(df, lab = rownames(df),
    x = 'log2FoldChange', y = 'padj',
    FCcutoff = 0.5, pCutoff = 0.05,
    xlim = fc_xlim, ylim = sig_ylim,
    selectLab = sel, max.overlaps = Inf, drawConnectors = TRUE,
    arrowheads = FALSE, widthConnectors = 0.3, caption = '',
    title    = label_text, subtitle = subtitle_text,
    legendPosition = if (show_legend) 'bottom' else 'none',
    colCustom = keyvals, labFace = 'italic') + coord_flip()
}

## Subtype-pseudobulk-VST violin + pairwise Wilcoxon (Fig 7 J/K style)
.ec_vst_violin <- function(df, gene, subtype_lbl, group_levels,
                            comparisons, palette) {
  if (is.null(df)) return(NULL)
  df$group <- factor(df$group, levels = group_levels)
  ggplot(df, aes(x = group, y = expr, fill = group)) +
    geom_violin(scale = 'width', alpha = 0.6, linewidth = 0.4) +
    geom_jitter(width = 0.08, size = 1.4, alpha = 0.9) +
    scale_fill_manual(values = palette) +
    ggpubr::stat_compare_means(comparisons = comparisons,
                               method = 'wilcox.test', size = 2.6) +
    labs(x = NULL,
         y = bquote(italic(.(gene)) ~ 'expression (VST)'),
         title = paste0(subtype_lbl, ' (pseudobulk)')) +
    theme_classic() +
    theme(legend.position = 'none',
          panel.border = element_rect(linewidth = 0.6, fill = NA, color = 'black'))
}

###############################################################################
## (F) Adult EC hdWGCNA modules (ecM1–ecM7) projected onto pediatric EC
##     Use module gene lists from EC_hdWGCNA_by_celltype.rds, score in pediatric.
###############################################################################
.cache_ec_mods <- file.path(CACHE, 'S7_adult_ec_module_genelists.rds')
if (file.exists(.cache_ec_mods)) {
  ec_module_genes <- readRDS(.cache_ec_mods)
} else {
  ec_ref      <- readRDS('./dependencies/shared/EC_hdWGCNA_by_celltype.rds')
  ec_modules  <- hdWGCNA::GetModules(ec_ref)
  ec_modules  <- subset(ec_modules, module != 'grey')
  ## Match Figure_7.R numbering: M_n where n = match(color, labels2colors(1:N)).
  ## (M1 = turquoise / capillary / Phase 1, M4 = yellow / arterial-angiogenic,
  ##  M7 = black / activated / Phase 2.)
  mapping_ec      <- WGCNA::labels2colors(1:100)
  ec_mod_colors   <- as.character(unique(ec_modules$module))
  ec_color2M      <- setNames(paste0('ecM', match(ec_mod_colors, mapping_ec)),
                              ec_mod_colors)
  ec_modules$ec_M <- ec_color2M[as.character(ec_modules$module)]
  ec_module_genes <- split(ec_modules$gene_name, ec_modules$ec_M)
  saveRDS(ec_module_genes, .cache_ec_mods)
  rm(ec_ref); invisible(gc())
}

ec_panel_modules <- intersect(paste0('ecM', 1:7), names(ec_module_genes))
for (m in ec_panel_modules) {
  g <- intersect(ec_module_genes[[m]], rownames(M1))
  if (length(g) >= 5) {
    M1 <- AddModuleScore(M1, list(g), assay = 'SCT', name = paste0(m, '_'))
  }
}

ec_score_cols <- paste0(ec_panel_modules, '_1')
ec_score_cols <- intersect(ec_score_cols, colnames(M1@meta.data))

Idents(M1) <- 'ec_subtype'
pdf(file.path(OUT, 'S7F_peds_EC_adult_modules_subtype.pdf'),
    width = 7, height = 3.5)
print(
  DotPlot(M1, features = ec_score_cols, dot.min = 0,
          col.min = 0, col.max = 2) +
    RotatedAxis() +
    scale_x_discrete(labels = sub('_1$', '', ec_score_cols)) +
    scale_color_gradient2(high = 'red', mid = 'grey95', low = 'blue') +
    xlab('Adult EC hdWGCNA modules') + ylab('Pediatric EC subtype') +
    theme(panel.border = element_rect(linewidth = 1, fill = NA, color = 'black'))
)
dev.off()

Idents(M1) <- 'group'
Idents(M1) <- factor(Idents(M1), levels = c('ped-NF','ped-HLHS','ped-RVF'))
pdf(file.path(OUT, 'S7F_peds_EC_adult_modules_group.pdf'),
    width = 7, height = 2.5)
print(
  DotPlot(M1, features = ec_score_cols, dot.min = 0,
          col.min = 0, col.max = 2) +
    RotatedAxis() +
    scale_x_discrete(labels = sub('_1$', '', ec_score_cols)) +
    scale_color_gradient2(high = 'red', mid = 'grey95', low = 'blue') +
    xlab('Adult EC hdWGCNA modules') + ylab('') +
    theme(panel.border = element_rect(linewidth = 1, fill = NA, color = 'black'))
)
dev.off()

## S7G composite (v57 S7G = 2-part EC module dotplot): pediatric-EC-subtype
## (top) stacked over disease-group (bottom).  Rebuilt as objects.
.s7g_dot <- function(idv, lv, ylab_) {
  Idents(M1) <- idv
  if (!is.null(lv)) Idents(M1) <- factor(Idents(M1), levels = lv)
  DotPlot(M1, features = ec_score_cols, dot.min = 0,
          col.min = 0, col.max = 2) +
    RotatedAxis() +
    scale_x_discrete(labels = sub('_1$', '', ec_score_cols)) +
    scale_color_gradient2(high = 'red', mid = 'grey95', low = 'blue') +
    xlab('Adult EC hdWGCNA modules') + ylab(ylab_) +
    theme(panel.border = element_rect(linewidth = 1, fill = NA, color = 'black'))
}
pdf(file.path(OUT, 'S7G_peds_EC_adult_modules_combined.pdf'),
    width = 7, height = 6)
print(.s7g_dot('ec_subtype', NULL, 'Pediatric EC subtype') /
      .s7g_dot('group', c('ped-NF','ped-HLHS','ped-RVF'), '') +
      plot_layout(heights = c(3.5, 2.5)))
dev.off()
cat('Wrote: S7G_peds_EC_adult_modules_combined.pdf\n')

###############################################################################
## (G) Pediatric EC subtype prevalence by disease group
###############################################################################
.prop <- as.data.frame(prop.table(
  table(M1$ec_subtype, M1$group), margin = 2))
colnames(.prop) <- c('subtype','group','freq')
.prop$group     <- factor(.prop$group, levels = c('ped-NF','ped-HLHS','ped-RVF'))

pdf(file.path(OUT, 'S7G_peds_EC_subtype_prev_stacked.pdf'),
    width = 4, height = 4)
print(
  ggplot(.prop, aes(x = group, y = freq, fill = subtype)) +
    geom_bar(stat = 'identity', width = 0.7) +
    theme_classic() +
    scale_y_continuous(expand = c(0,0)) +
    xlab('Disease state') + ylab('Fraction of EC') +
    theme(panel.border = element_rect(linewidth = 1, fill = NA, color = 'black'))
)
dev.off()

## per-patient boxplot view (KW / ANOVA)
.pt <- as.data.frame(prop.table(
  table(M1$ec_subtype, M1$sample), margin = 2))
colnames(.pt) <- c('subtype','sample','freq')
.pt$group <- M1$group[match(.pt$sample, M1$sample)]
.pt$group <- factor(.pt$group, levels = c('ped-NF','ped-HLHS','ped-RVF'))

pdf(file.path(OUT, 'S7G_peds_EC_subtype_freq_box.pdf'),
    width = 8, height = 4)
print(
  ggboxplot(.pt, x = 'group', y = 'freq', fill = 'group') +
    facet_wrap(~ subtype, ncol = 5, scales = 'free_y') +
    stat_compare_means(method = 'anova', label = 'p.format') +
    scale_fill_manual(values = disease_pal) +
    theme(strip.text = element_text(face = 'bold'))
)
dev.off()

###############################################################################
## (H) Pediatric EC angiogenesis-related GO scores
###############################################################################
.cache_angio <- file.path(CACHE, 'S7_angio_genelists.rds')
if (file.exists(.cache_angio)) {
  angio_sets <- readRDS(.cache_angio)
} else {
  ## biomaRt is slow / can fail offline — guard with a fallback gene list.
  angio_sets <- tryCatch({
    suppressPackageStartupMessages(library(biomaRt))
    em <- biomaRt::useMart('ensembl', dataset = 'hsapiens_gene_ensembl')
    gset <- function(go) {
      d <- biomaRt::getBM(
        attributes = c('hgnc_symbol','ensembl_transcript_id','go_id'),
        filters    = 'go', values = go, mart = em)
      unique(d$hgnc_symbol)
    }
    list(
      AllAngio       = gset('GO:0001525'),
      RegAngio       = gset('GO:0045765'),
      PositiveAngio  = gset('GO:0045766'),
      NegativeAngio  = gset('GO:0016525'),
      CoronaryAngio  = gset('GO:0060978')
    )
  }, error = function(e) {
    cat('biomaRt unavailable — using minimal fallback gene set.\n')
    list(
      AllAngio       = c('VEGFA','VEGFB','VEGFC','KDR','FLT1','FLT4',
                          'ANGPT1','ANGPT2','TEK','TIE1','PECAM1','CDH5'),
      RegAngio       = c('THBS1','THBS2','TIMP3','SERPINF1','SPARC',
                          'NRP1','NRP2','EFNB2','EPHB4'),
      PositiveAngio  = c('VEGFA','VEGFB','VEGFC','FGF2','PDGFB',
                          'ANGPT1','HGF','MMP2'),
      NegativeAngio  = c('THBS1','THBS2','TIMP3','SERPINF1','ANGPTL4'),
      CoronaryAngio  = c('VEGFA','APLN','APLNR','CXCL12','CXCR4')
    )
  })
  saveRDS(angio_sets, .cache_angio)
}
for (nm in names(angio_sets)) {
  g <- intersect(angio_sets[[nm]], rownames(M1))
  if (length(g) >= 3) {
    M1 <- AddModuleScore(M1, list(g), assay = 'SCT', name = paste0(nm, 'Score_'))
  }
}
.angio_cols <- paste0(names(angio_sets), 'Score_1')
.angio_cols <- intersect(.angio_cols, colnames(M1@meta.data))

pdf(file.path(OUT, 'S7H_peds_EC_angio_combined.pdf'),
    width = 10, height = 4)
print(
  VlnPlot(M1, .angio_cols, group.by = 'group', pt.size = 0,
          ncol = length(.angio_cols)) &
    scale_fill_manual(values = disease_pal)
)
dev.off()

## per-patient pseudobulk companion (n≈14, ANOVA)
.peds_pb_box(M1, .angio_cols,
             label_map = setNames(sub('Score$', '', sub('_1$','',.angio_cols)),
                                  sub('_1$','',.angio_cols)),
             out_path = file.path(OUT, 'S7H_peds_EC_angio_combined_box.pdf'),
             width = 11, height = 2.6,
             data_cache_path = file.path(CACHE, 'S7H_peds_pat.rds'))

###############################################################################
## (I) MECOM expression — reduced in ped-HLHS/RVF (cf. adult induction)
###############################################################################
if ('MECOM' %in% rownames(M1)) {
  pdf(file.path(OUT, 'S7J_subpanel_MECOM_peds_artery_violin.pdf'),
      width = 3, height = 3)
  print(
    VlnPlot(subset(M1, ec_subtype %in% c('Artery1','Artery2')), 'MECOM',
            group.by = 'group', pt.size = 0) +
      scale_fill_manual(values = disease_pal) +
      ggtitle('MECOM - pediatric arterial EC')
  )
  dev.off()

  pdf(file.path(OUT, 'S7J_subpanel_MECOM_peds_dotplot.pdf'),
      width = 4, height = 3)
  print(
    DotPlot(M1, features = 'MECOM', group.by = 'ec_subtype', col.min = 0) +
      RotatedAxis() +
      scale_color_gradient2(high = 'red', mid = 'grey95', low = 'blue') +
      theme(panel.border = element_rect(linewidth = 1, fill = NA, color = 'black'))
  )
  dev.off()

  .peds_pb_gene_box(M1, 'MECOM',
                    cell_filter = M1$ec_subtype %in% c('Artery1','Artery2'),
                    out_path = file.path(OUT, 'S7J_subpanel_MECOM_peds_artery_box.pdf'),
                    title = 'MECOM - peds arterial',
                    data_cache_path = file.path(CACHE, 'S7J_subpanel_MECOM_peds_pat.rds'))
}

###############################################################################
## (J) Arterialization TF score (Wythe/De Val 2013 arterial-EC TF code)
###############################################################################
art_tf <- intersect(c('AFF3','HEY2','SOX17','MSX1','EMX2','TOX2','PRDM16'),
                    rownames(M1))
M1 <- AddModuleScore(M1, list(art_tf), name = 'ArtTF_')
pdf(file.path(OUT, 'S7J_subpanel_ArtTF_peds_violin.pdf'), width = 3.2, height = 3)
print(
  VlnPlot(M1, 'ArtTF_1', group.by = 'group', pt.size = 0) +
    scale_fill_manual(values = disease_pal) +
    ggtitle('Arterialization TF score (peds EC)')
)
dev.off()

.peds_pb_box(M1, 'ArtTF_1',
             label_map = c(ArtTF = 'Arterialization TF'),
             out_path = file.path(OUT, 'S7J_subpanel_ArtTF_peds_box.pdf'),
             width = 3.4, height = 2.8,
             ylab = 'Score (per-patient mean)',
             data_cache_path = file.path(CACHE, 'S7J_subpanel_ArtTF_peds_pat.rds'))

###############################################################################
## (K) Notch target score (canonical arterial Notch program)
###############################################################################
notch_set <- intersect(
  c('RND1','RASSF10','HEY2','EFNB2','PRICKLE2','SLC46A3','DLL4',
    'RAPGEF5','PLPP3','HEY1','LRRC32','SAT1','HES4','ANKRD33B','FLT1'),
  rownames(M1))
M1 <- AddModuleScore(M1, list(notch_set), name = 'Notch_')
pdf(file.path(OUT, 'S7J_subpanel_Notch_peds_violin.pdf'), width = 3.2, height = 3)
print(
  VlnPlot(M1, 'Notch_1', group.by = 'group', pt.size = 0) +
    scale_fill_manual(values = disease_pal) +
    ggtitle('Notch target score (peds EC)')
)
dev.off()

.peds_pb_box(M1, 'Notch_1',
             label_map = c(Notch = 'Notch target'),
             out_path = file.path(OUT, 'S7J_subpanel_Notch_peds_box.pdf'),
             width = 3.4, height = 2.8,
             ylab = 'Score (per-patient mean)',
             data_cache_path = file.path(CACHE, 'S7J_subpanel_Notch_peds_pat.rds'))

###############################################################################
## (L) SMAD1 — capillary-restricted view
###############################################################################
if ('SMAD1' %in% rownames(M1)) {
  pdf(file.path(OUT, 'S7J_subpanel_SMAD1_peds_capillary_violin.pdf'),
      width = 3, height = 3)
  print(
    VlnPlot(subset(M1, ec_subtype == 'Capillary'), 'SMAD1',
            group.by = 'group', pt.size = 0) +
      scale_fill_manual(values = disease_pal) +
      ggtitle('SMAD1 - pediatric capillary EC')
  )
  dev.off()

  .peds_pb_gene_box(M1, 'SMAD1',
                    cell_filter = M1$ec_subtype == 'Capillary',
                    out_path = file.path(OUT, 'S7J_subpanel_SMAD1_peds_capillary_box.pdf'),
                    title = 'SMAD1 - capillary',
                    data_cache_path = file.path(CACHE, 'S7J_subpanel_SMAD1_peds_pat.rds'))
}

###############################################################################
## (M) NR2F2 — venous-restricted view
###############################################################################
if ('NR2F2' %in% rownames(M1)) {
  pdf(file.path(OUT, 'S7J_subpanel_NR2F2_peds_vein_violin.pdf'),
      width = 3, height = 3)
  print(
    VlnPlot(subset(M1, ec_subtype == 'Vein'), 'NR2F2',
            group.by = 'group', pt.size = 0) +
      scale_fill_manual(values = disease_pal) +
      ggtitle('NR2F2 - pediatric venous EC')
  )
  dev.off()

  .peds_pb_gene_box(M1, 'NR2F2',
                    cell_filter = M1$ec_subtype == 'Vein',
                    out_path = file.path(OUT, 'S7J_subpanel_NR2F2_peds_vein_box.pdf'),
                    title = 'NR2F2 - vein',
                    data_cache_path = file.path(CACHE, 'S7J_subpanel_NR2F2_peds_pat.rds'))
}

###############################################################################
## (N) Phase 1 EC: vasoprotective erosion + IFN engagement (data-anchored)
##     Vasoprotective composite: pooled antioxidant, RAAS/vasoreactive,
##     anti-inflammatory, ECM stabilization, and EC quiescence categories.
##     ACKR1 is biphasic and gets its own panel.
###############################################################################
phase1_vasoprot <- intersect(
  c(phase1_cats$Antioxidant,
    phase1_cats$RAAS_vasoreactive,
    setdiff(phase1_cats$Anti_inflammatory, 'ACKR1'),
    phase1_cats$ECM_stabilization,
    phase1_cats$EC_quiescence),
  rownames(M1))
phase1_ifn <- intersect(
  c('IFI44L','IFI44','RSAD2','MX1','MX2','OAS1','OAS2','OAS3',
    'IRF7','STAT1','HERC6','EPSTI1','XAF1','TRIM22','ISG15',
    'IFIT1','IFIT3','IFITM1','IFITM3'),
  rownames(M1))
M1 <- AddModuleScore(M1, list(phase1_vasoprot), name = 'Vasoprot_')
M1 <- AddModuleScore(M1, list(phase1_ifn),       name = 'IFN_')

## (N1) Vasoprotective + IFN composite scores (per-patient pseudobulk)
.peds_pb_box(M1, c('Vasoprot_1','IFN_1'),
             label_map = c(Vasoprot = 'Vasoprotective program',
                           IFN       = 'IFN engagement'),
             out_path = file.path(OUT, 'S7I_peds_EC_phase1_box.pdf'),
             width = 6, height = 2.8,
             ylab = 'Score (per-patient mean)',
             data_cache_path = file.path(CACHE, 'S7I_peds_EC_phase1_pat.rds'))

## (N2) GPX3 single-gene per-patient pseudobulk box (paragraph headline)
if ('GPX3' %in% rownames(M1)) {
  .peds_pb_gene_box(M1, 'GPX3',
                    out_path = file.path(OUT, 'S7I_supp_peds_EC_GPX3_box.pdf'),
                    title = 'GPX3 - peds EC',
                    data_cache_path = file.path(CACHE, 'S7I_supp_GPX3_peds_pat.rds'))
}

## (N3) ACKR1 biphasic single-gene per-patient pseudobulk box
if ('ACKR1' %in% rownames(M1)) {
  .peds_pb_gene_box(M1, 'ACKR1',
                    out_path = file.path(OUT, 'S7I_supp_peds_EC_ACKR1_box.pdf'),
                    title = 'ACKR1 - peds EC (biphasic)',
                    data_cache_path = file.path(CACHE, 'S7I_supp_ACKR1_peds_pat.rds'))
}

## (N4) Phase 1 pan-EC pseudobulk DESeq2 (peds: ped-HLHS vs ped-NF)
.cache_peds_phase1 <- file.path(CACHE, 'S7_phase1_peds_pan_DE.rds')
peds_phase1_de <- .ec_pan_pseudobulk_de(
  M1,
  contrast    = c('group', 'ped-HLHS', 'ped-NF'),
  patient_col = 'sample',
  cache_path  = .cache_peds_phase1
)

###############################################################################
## (O) NRG1 in Endocardial EC — pseudobulk DESeq2-VST per donor (Fig 7J rep)
###############################################################################
.cache_O_peds_NRG1 <- file.path(CACHE, 'S7_supp_peds_endocardial_NRG1_vst.rds')
.peds_NRG1_endo_df <- .ec_subtype_pseudobulk_vst(
  M1,
  subtype_col    = 'ec_subtype',
  subtype_values = c('Endocardium','Endocardial'),
  gene           = 'NRG1',
  patient_col    = 'sample',
  group_levels   = c('ped-NF','ped-HLHS','ped-RVF'),
  cache_path     = .cache_O_peds_NRG1
)
if (!is.null(.peds_NRG1_endo_df)) {
  pdf(file.path(OUT, 'S7_supp_peds_EC_NRG1_endocardial.pdf'),
      width = 3.4, height = 3.4)
  print(.ec_vst_violin(.peds_NRG1_endo_df, 'NRG1', 'Pediatric Endocardial',
                       group_levels = c('ped-NF','ped-HLHS','ped-RVF'),
                       comparisons  = list(c('ped-NF','ped-HLHS'),
                                           c('ped-HLHS','ped-RVF'),
                                           c('ped-NF','ped-RVF')),
                       palette      = disease_pal))
  dev.off()
}

###############################################################################
## ── Phase-2 framework: peds pan-EC pseudobulk DESeq2 + volcano (Fig 7H rep) ─
###############################################################################
.cache_peds_phase2 <- file.path(CACHE, 'S7_phase2_peds_pan_DE.rds')
peds_pan_de <- .ec_pan_pseudobulk_de(
  M1,
  contrast    = c('group', 'ped-RVF', 'ped-HLHS'),
  patient_col = 'sample',
  cache_path  = .cache_peds_phase2
)

pdf(file.path(OUT, 'S7H_supp_peds_EC_phase2_volcano.pdf'),
    width = 5, height = 5)
print(.ec_volcano(peds_pan_de,
                  label_text   = 'Pediatric pan-EC (snRNA-seq)',
                  fc_xlim      = c(-3, 3), sig_ylim = c(0, 5),
                  show_legend  = TRUE,
                  subtitle_text = 'pan-EC pseudobulk - ped-RVF vs ped-HLHS (Phase 2)'))
dev.off()

###############################################################################
## ── Adult RV ENDOTHELIUM — comparison panels (H–M equivalents) ─────────────
##   Generates _adult versions of S7{H,I,J,K,L,M} cell-level violin + per-
##   patient pseudobulk box plots so the peds vs adult readouts sit side-by-
##   side. Adult metadata: $patient (per-sample), $Names (EC subtype:
##   Arterial/Capillary/Venous/Lymph/Endocardial), $group (NF/pRV/RVF).
###############################################################################
rm(M1); invisible(gc())
cat('\n=== Loading adult EC object for comparison ===\n')
A1 <- readRDS('./dependencies/shared/EC_hdWGCNA_by_celltype.rds')
DefaultAssay(A1) <- 'SCT'

## adult palette uses the global NF/pRV/RVF keys
adult_pal       <- c(NF = '#4DAF4A', pRV = '#377EB8', RVF = '#E41A1C')
adult_levels    <- c('NF','pRV','RVF')
A1$ec_subtype   <- as.character(A1$Names)

## angio gene sets — reuse cached lists
angio_sets_a <- readRDS(file.path(CACHE, 'S7_angio_genelists.rds'))
for (nm in names(angio_sets_a)) {
  g <- intersect(angio_sets_a[[nm]], rownames(A1))
  if (length(g) >= 3) {
    A1 <- AddModuleScore(A1, list(g), assay = 'SCT', name = paste0(nm, 'Score_'))
  }
}
.angio_cols_a <- intersect(paste0(names(angio_sets_a), 'Score_1'),
                           colnames(A1@meta.data))

## ArtTF + Notch
art_tf_a <- intersect(c('AFF3','HEY2','SOX17','MSX1','EMX2','TOX2','PRDM16'),
                      rownames(A1))
A1 <- AddModuleScore(A1, list(art_tf_a), name = 'ArtTF_')
notch_set_a <- intersect(
  c('RND1','RASSF10','HEY2','EFNB2','PRICKLE2','SLC46A3','DLL4',
    'RAPGEF5','PLPP3','HEY1','LRRC32','SAT1','HES4','ANKRD33B','FLT1'),
  rownames(A1))
A1 <- AddModuleScore(A1, list(notch_set_a), name = 'Notch_')

## (H_adult) angiogenesis violin + per-patient box
pdf(file.path(OUT, 'S7H_adult_EC_angio_combined.pdf'), width = 10, height = 4)
print(
  VlnPlot(A1, .angio_cols_a, group.by = 'group', pt.size = 0,
          ncol = length(.angio_cols_a)) &
    scale_fill_manual(values = adult_pal)
)
dev.off()
.peds_pb_box(A1, .angio_cols_a,
             label_map = setNames(sub('Score$', '', sub('_1$','',.angio_cols_a)),
                                  sub('_1$','',.angio_cols_a)),
             out_path = file.path(OUT, 'S7H_adult_EC_angio_combined_box.pdf'),
             width = 11, height = 2.6,
             patient_col = 'patient',
             group_levels = adult_levels, palette = adult_pal,
             data_cache_path = file.path(CACHE, 'S7H_adult_pat.rds'))

## (I_adult) MECOM
if ('MECOM' %in% rownames(A1)) {
  pdf(file.path(OUT, 'S7J_subpanel_MECOM_adult_artery_violin.pdf'),
      width = 3, height = 3)
  print(
    VlnPlot(subset(A1, ec_subtype == 'Arterial'), 'MECOM',
            group.by = 'group', pt.size = 0) +
      scale_fill_manual(values = adult_pal) +
      ggtitle('MECOM - adult arterial EC')
  )
  dev.off()
  .peds_pb_gene_box(A1, 'MECOM',
                    cell_filter = A1$ec_subtype == 'Arterial',
                    out_path = file.path(OUT, 'S7J_subpanel_MECOM_adult_artery_box.pdf'),
                    title = 'MECOM - adult arterial',
                    patient_col = 'patient',
                    group_levels = adult_levels, palette = adult_pal,
                    data_cache_path = file.path(CACHE, 'S7J_subpanel_MECOM_adult_pat.rds'))
}

## (J_adult) ArtTF score
pdf(file.path(OUT, 'S7J_subpanel_ArtTF_adult_violin.pdf'), width = 3.2, height = 3)
print(
  VlnPlot(A1, 'ArtTF_1', group.by = 'group', pt.size = 0) +
    scale_fill_manual(values = adult_pal) +
    ggtitle('Arterialization TF score (adult EC)')
)
dev.off()
.peds_pb_box(A1, 'ArtTF_1',
             label_map = c(ArtTF = 'Arterialization TF'),
             out_path = file.path(OUT, 'S7J_subpanel_ArtTF_adult_box.pdf'),
             width = 3.4, height = 2.8,
             ylab = 'Score (per-patient mean)',
             patient_col = 'patient',
             group_levels = adult_levels, palette = adult_pal,
             data_cache_path = file.path(CACHE, 'S7J_subpanel_ArtTF_adult_pat.rds'))

## (K_adult) Notch target score
pdf(file.path(OUT, 'S7J_subpanel_Notch_adult_violin.pdf'), width = 3.2, height = 3)
print(
  VlnPlot(A1, 'Notch_1', group.by = 'group', pt.size = 0) +
    scale_fill_manual(values = adult_pal) +
    ggtitle('Notch target score (adult EC)')
)
dev.off()
.peds_pb_box(A1, 'Notch_1',
             label_map = c(Notch = 'Notch target'),
             out_path = file.path(OUT, 'S7J_subpanel_Notch_adult_box.pdf'),
             width = 3.4, height = 2.8,
             ylab = 'Score (per-patient mean)',
             patient_col = 'patient',
             group_levels = adult_levels, palette = adult_pal,
             data_cache_path = file.path(CACHE, 'S7J_subpanel_Notch_adult_pat.rds'))

## (L_adult) SMAD1 capillary
if ('SMAD1' %in% rownames(A1)) {
  pdf(file.path(OUT, 'S7J_subpanel_SMAD1_adult_capillary_violin.pdf'),
      width = 3, height = 3)
  print(
    VlnPlot(subset(A1, ec_subtype == 'Capillary'), 'SMAD1',
            group.by = 'group', pt.size = 0) +
      scale_fill_manual(values = adult_pal) +
      ggtitle('SMAD1 - adult capillary EC')
  )
  dev.off()
  .peds_pb_gene_box(A1, 'SMAD1',
                    cell_filter = A1$ec_subtype == 'Capillary',
                    out_path = file.path(OUT, 'S7J_subpanel_SMAD1_adult_capillary_box.pdf'),
                    title = 'SMAD1 - adult capillary',
                    patient_col = 'patient',
                    group_levels = adult_levels, palette = adult_pal,
                    data_cache_path = file.path(CACHE, 'S7J_subpanel_SMAD1_adult_pat.rds'))
}

## (M_adult) NR2F2 vein
if ('NR2F2' %in% rownames(A1)) {
  pdf(file.path(OUT, 'S7J_subpanel_NR2F2_adult_vein_violin.pdf'),
      width = 3, height = 3)
  print(
    VlnPlot(subset(A1, ec_subtype == 'Venous'), 'NR2F2',
            group.by = 'group', pt.size = 0) +
      scale_fill_manual(values = adult_pal) +
      ggtitle('NR2F2 - adult venous EC')
  )
  dev.off()
  .peds_pb_gene_box(A1, 'NR2F2',
                    cell_filter = A1$ec_subtype == 'Venous',
                    out_path = file.path(OUT, 'S7J_subpanel_NR2F2_adult_vein_box.pdf'),
                    title = 'NR2F2 - adult venous',
                    patient_col = 'patient',
                    group_levels = adult_levels, palette = adult_pal,
                    data_cache_path = file.path(CACHE, 'S7J_subpanel_NR2F2_adult_pat.rds'))
}

###############################################################################
## (N_adult) Phase 1 EC: vasoprotective erosion + IFN engagement (data-anchored)
###############################################################################
phase1_vasoprot_a <- intersect(
  c(phase1_cats$Antioxidant,
    phase1_cats$RAAS_vasoreactive,
    setdiff(phase1_cats$Anti_inflammatory, 'ACKR1'),
    phase1_cats$ECM_stabilization,
    phase1_cats$EC_quiescence),
  rownames(A1))
phase1_ifn_a <- intersect(
  c('IFI44L','IFI44','RSAD2','MX1','MX2','OAS1','OAS2','OAS3',
    'IRF7','STAT1','HERC6','EPSTI1','XAF1','TRIM22','ISG15',
    'IFIT1','IFIT3','IFITM1','IFITM3'),
  rownames(A1))
A1 <- AddModuleScore(A1, list(phase1_vasoprot_a), name = 'Vasoprot_')
A1 <- AddModuleScore(A1, list(phase1_ifn_a),       name = 'IFN_')

## (N1) Vasoprotective + IFN composite scores
.peds_pb_box(A1, c('Vasoprot_1','IFN_1'),
             label_map = c(Vasoprot = 'Vasoprotective program',
                           IFN       = 'IFN engagement'),
             out_path = file.path(OUT, 'S7N_adult_EC_phase1_box.pdf'),
             width = 6, height = 2.8,
             ylab = 'Score (per-patient mean)',
             patient_col = 'patient',
             group_levels = adult_levels, palette = adult_pal,
             data_cache_path = file.path(CACHE, 'S7N_adult_pat.rds'))

## (N2) GPX3 single-gene per-patient pseudobulk box
if ('GPX3' %in% rownames(A1)) {
  .peds_pb_gene_box(A1, 'GPX3',
                    out_path = file.path(OUT, 'S7N_adult_EC_GPX3_box.pdf'),
                    title = 'GPX3 - adult EC',
                    patient_col = 'patient',
                    group_levels = adult_levels, palette = adult_pal,
                    data_cache_path = file.path(CACHE, 'S7N_GPX3_adult_pat.rds'))
}

## (N3) ACKR1 biphasic single-gene per-patient pseudobulk box
if ('ACKR1' %in% rownames(A1)) {
  .peds_pb_gene_box(A1, 'ACKR1',
                    out_path = file.path(OUT, 'S7N_adult_EC_ACKR1_box.pdf'),
                    title = 'ACKR1 - adult EC (biphasic)',
                    patient_col = 'patient',
                    group_levels = adult_levels, palette = adult_pal,
                    data_cache_path = file.path(CACHE, 'S7N_ACKR1_adult_pat.rds'))
}

## (N4) Phase 1 pan-EC pseudobulk DESeq2 (adult: pRV vs NF)
.cache_adult_phase1 <- file.path(CACHE, 'S7_phase1_adult_pan_DE.rds')
adult_phase1_de <- .ec_pan_pseudobulk_de(
  A1,
  contrast    = c('group', 'pRV', 'NF'),
  patient_col = 'patient',
  cache_path  = .cache_adult_phase1
)

###############################################################################
## (O_adult) NRG1 in Endocardial EC — pseudobulk DESeq2-VST (Fig 7J)
###############################################################################
.cache_O_adult_NRG1 <- file.path(CACHE, 'S7_supp_adult_endocardial_NRG1_vst.rds')
.adult_NRG1_endo_df <- .ec_subtype_pseudobulk_vst(
  A1,
  subtype_col    = 'ec_subtype',
  subtype_values = c('Endocardial'),
  gene           = 'NRG1',
  patient_col    = 'patient',
  group_levels   = adult_levels,
  cache_path     = .cache_O_adult_NRG1
)
if (!is.null(.adult_NRG1_endo_df)) {
  pdf(file.path(OUT, 'S7_supp_adult_EC_NRG1_endocardial.pdf'),
      width = 3.4, height = 3.4)
  print(.ec_vst_violin(.adult_NRG1_endo_df, 'NRG1', 'Adult Endocardial',
                       group_levels = adult_levels,
                       comparisons  = list(c('NF','pRV'),
                                           c('pRV','RVF'),
                                           c('NF','RVF')),
                       palette      = adult_pal))
  dev.off()
}

###############################################################################
## ── Phase-2 framework: adult pan-EC pseudobulk DESeq2 + volcano (Fig 7H rep) ─
###############################################################################
.cache_adult_phase2 <- file.path(CACHE, 'S7_phase2_adult_pan_DE.rds')
adult_pan_de <- .ec_pan_pseudobulk_de(
  A1,
  contrast    = c('group', 'RVF', 'pRV'),
  patient_col = 'patient',
  cache_path  = .cache_adult_phase2
)

pdf(file.path(OUT, 'S7H_supp_adult_EC_phase2_volcano.pdf'),
    width = 5, height = 5)
print(.ec_volcano(adult_pan_de,
                  label_text   = 'Adult pan-EC (snRNA-seq)',
                  fc_xlim      = c(-3, 3), sig_ylim = c(0, 5),
                  show_legend  = TRUE,
                  subtitle_text = 'pan-EC pseudobulk - RVF vs pRV (Phase 2)'))
dev.off()

###############################################################################
## ── Combined comparison panels (peds top / adult bottom) ───────────────────
##   Reads cached per-patient long-format dataframes from each cohort and
##   stacks them into a 2-row pseudobulk panel per analysis. Separate
##   peds_*/adult_* PDFs above are kept as fallback views.
###############################################################################
peds_pal_ec <- c('ped-NF'   = '#4DAF4A',
                 'ped-HLHS' = '#377EB8',
                 'ped-RVF'  = '#E41A1C')
adult_pal_ec <- c(NF = '#4DAF4A', pRV = '#377EB8', RVF = '#E41A1C')

.combined_panel_box <- function(peds_path, adult_path, out_path,
                                  width, height, ylab,
                                  is_gene = FALSE, gene_label = NULL) {
  if (!file.exists(peds_path) || !file.exists(adult_path)) {
    cat('Skipping ', basename(out_path), ' — missing cache\n', sep = '')
    return(invisible(NULL))
  }
  peds_df  <- readRDS(peds_path)
  adult_df <- readRDS(adult_path)
  if (is_gene) {
    peds_df$program  <- gene_label
    adult_df$program <- gene_label
    peds_df$score    <- peds_df$expr
    adult_df$score   <- adult_df$expr
    peds_df$expr <- NULL; adult_df$expr <- NULL
  }
  peds_df$cohort  <- 'Pediatric'
  adult_df$cohort <- 'Adult'
  combined <- dplyr::bind_rows(peds_df, adult_df)
  combined$cohort <- factor(combined$cohort, levels = c('Pediatric','Adult'))
  full_pal <- c(peds_pal_ec, adult_pal_ec)
  p <- ggplot(combined, aes(x = group, y = score, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.85, linewidth = 0.4) +
    geom_jitter(width = 0.18, size = 2.4, shape = 21,
                colour = 'black', stroke = 0.3, alpha = 0.9) +
    ## facet_wrap (not facet_grid): independent x AND y per panel so the
    ## peds row shows ped-NF/ped-HLHS/ped-RVF, the adult row shows
    ## NF/pRV/RVF, and the y-axis is NOT shared across panels. nrow=2 →
    ## row1 = Pediatric, row2 = Adult (cohort-major fill).
    facet_wrap(vars(cohort, program), scales = 'free', nrow = 2) +
    scale_fill_manual(values = full_pal) +
    ggpubr::stat_compare_means(method = 'anova',
                               label = 'p.format', size = 2.6) +
    theme_classic() +
    xlab('') + ylab(ylab) +
    theme(strip.text      = element_text(face = 'bold'),
          legend.position = 'none',
          panel.border    = element_rect(linewidth = 0.6, fill = NA, color = 'black'),
          axis.text.x     = element_text(angle = 35, hjust = 1))
  pdf(out_path, width = width, height = height); print(p); dev.off()
}

## H – angiogenesis scores (5 programs)
.combined_panel_box(file.path(CACHE, 'S7H_peds_pat.rds'),
                    file.path(CACHE, 'S7H_adult_pat.rds'),
                    file.path(OUT,   'S7H_combined_EC_angio_box.pdf'),
                    width = 11, height = 4.6,
                    ylab = 'Module score (per-patient mean)')

## I – MECOM in arterial EC
.combined_panel_box(file.path(CACHE, 'S7J_subpanel_MECOM_peds_pat.rds'),
                    file.path(CACHE, 'S7J_subpanel_MECOM_adult_pat.rds'),
                    file.path(OUT,   'S7I_combined_EC_MECOM_artery_box.pdf'),
                    width = 3.6, height = 5,
                    ylab = 'MECOM expression (per-patient mean)',
                    is_gene = TRUE, gene_label = 'MECOM')

## J – Arterialization TF score
.combined_panel_box(file.path(CACHE, 'S7J_subpanel_ArtTF_peds_pat.rds'),
                    file.path(CACHE, 'S7J_subpanel_ArtTF_adult_pat.rds'),
                    file.path(OUT,   'S7J_combined_EC_ArtTF_box.pdf'),
                    width = 3.6, height = 5,
                    ylab = 'Score (per-patient mean)')

## K – Notch target score
.combined_panel_box(file.path(CACHE, 'S7J_subpanel_Notch_peds_pat.rds'),
                    file.path(CACHE, 'S7J_subpanel_Notch_adult_pat.rds'),
                    file.path(OUT,   'S7K_combined_EC_Notch_box.pdf'),
                    width = 3.6, height = 5,
                    ylab = 'Score (per-patient mean)')

## L – SMAD1 in capillary EC
.combined_panel_box(file.path(CACHE, 'S7J_subpanel_SMAD1_peds_pat.rds'),
                    file.path(CACHE, 'S7J_subpanel_SMAD1_adult_pat.rds'),
                    file.path(OUT,   'S7L_combined_EC_SMAD1_capillary_box.pdf'),
                    width = 3.6, height = 5,
                    ylab = 'SMAD1 expression (per-patient mean)',
                    is_gene = TRUE, gene_label = 'SMAD1')

## M – NR2F2 in vein EC
.combined_panel_box(file.path(CACHE, 'S7J_subpanel_NR2F2_peds_pat.rds'),
                    file.path(CACHE, 'S7J_subpanel_NR2F2_adult_pat.rds'),
                    file.path(OUT,   'S7M_combined_EC_NR2F2_vein_box.pdf'),
                    width = 3.6, height = 5,
                    ylab = 'NR2F2 expression (per-patient mean)',
                    is_gene = TRUE, gene_label = 'NR2F2')

## S7I (v57): vasoprotective + IFN composite, 2x2 (cohort x program).
## FIX: peds_path was wrongly the MECOM cache; use the peds vasoprot+IFN
## cache (S7I_peds_EC_phase1_pat.rds) paired with the adult one.
.combined_panel_box(file.path(CACHE, 'S7I_peds_EC_phase1_pat.rds'),
                    file.path(CACHE, 'S7N_adult_pat.rds'),
                    file.path(OUT,   'S7N_combined_EC_phase1_box.pdf'),
                    width = 6, height = 5,
                    ylab = 'Score (per-patient mean)')

## S7J (v57): one 5x2 grid — MECOM | ArtTF | Notch | SMAD1 | NR2F2 columns,
## Pediatric (row 1) over Adult (row 2). Assembled from the cached
## per-patient _pat.rds for each readout x cohort.
.s7j_specs <- list(
  list(p = 'S7J_subpanel_MECOM_peds_pat.rds',  a = 'S7J_subpanel_MECOM_adult_pat.rds',  lab = 'MECOM',  gene = TRUE),
  list(p = 'S7J_subpanel_ArtTF_peds_pat.rds',  a = 'S7J_subpanel_ArtTF_adult_pat.rds',  lab = 'Arterialization TF', gene = FALSE),
  list(p = 'S7J_subpanel_Notch_peds_pat.rds',  a = 'S7J_subpanel_Notch_adult_pat.rds',  lab = 'Notch target',       gene = FALSE),
  list(p = 'S7J_subpanel_SMAD1_peds_pat.rds',  a = 'S7J_subpanel_SMAD1_adult_pat.rds',  lab = 'SMAD1',  gene = TRUE),
  list(p = 'S7J_subpanel_NR2F2_peds_pat.rds',  a = 'S7J_subpanel_NR2F2_adult_pat.rds',  lab = 'NR2F2',  gene = TRUE))
.s7j_rows <- list()
for (.sp in .s7j_specs) {
  .pf <- file.path(CACHE, .sp$p); .af <- file.path(CACHE, .sp$a)
  if (!file.exists(.pf) || !file.exists(.af)) next
  .pd <- readRDS(.pf); .ad <- readRDS(.af)
  .norm <- function(d, ch) {
    if (.sp$gene) d$score <- d$expr else
      d$score <- d[[setdiff(intersect(c('score','expr'), names(d)), 'expr')[1]]]
    if (is.null(d$score) && 'expr' %in% names(d)) d$score <- d$expr
    data.frame(patient = d$patient, group = as.character(d$group),
               score = d$score, program = .sp$lab, cohort = ch,
               stringsAsFactors = FALSE)
  }
  .s7j_rows[[length(.s7j_rows) + 1]] <-
    rbind(.norm(.pd, 'Pediatric'), .norm(.ad, 'Adult'))
}
if (length(.s7j_rows) > 0) {
  .s7j_df <- do.call(rbind, .s7j_rows)
  .s7j_df$program <- factor(.s7j_df$program,
    levels = c('MECOM','Arterialization TF','Notch target','SMAD1','NR2F2'))
  .s7j_df$cohort  <- factor(.s7j_df$cohort, levels = c('Pediatric','Adult'))
  .s7j_df$group   <- factor(.s7j_df$group,
    levels = c('ped-NF','ped-HLHS','ped-RVF','NF','pRV','RVF'))
  p_S7J <- ggplot(.s7j_df, aes(x = group, y = score, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.85, linewidth = 0.4) +
    geom_jitter(width = 0.18, size = 2.0, shape = 21,
                colour = 'black', stroke = 0.3, alpha = 0.9) +
    ## facet_wrap (not facet_grid): independent x AND y per panel so the
    ## peds row shows ped-NF/ped-HLHS/ped-RVF, the adult row shows
    ## NF/pRV/RVF, and the y-axis is NOT shared across panels. nrow=2 →
    ## row1 = Pediatric, row2 = Adult (cohort-major fill).
    facet_wrap(vars(cohort, program), scales = 'free', nrow = 2) +
    scale_fill_manual(values = c(peds_pal_ec, adult_pal_ec)) +
    ggpubr::stat_compare_means(method = 'anova', label = 'p.format',
                               size = 2.4) +
    theme_classic() + xlab('') + ylab('Per-patient mean (score / expr)') +
    theme(strip.text = element_text(face = 'bold'),
          legend.position = 'none',
          panel.border = element_rect(linewidth = 0.6, fill = NA,
                                      color = 'black'),
          axis.text.x = element_text(angle = 35, hjust = 1))
  pdf(file.path(OUT, 'S7J_combined_5x2.pdf'), width = 13, height = 5.4)
  print(p_S7J)
  dev.off()
  cat('Wrote: S7J_combined_5x2.pdf (5 readouts x peds/adult)\n')
}

## N2 – GPX3 single-gene combined
.combined_panel_box(file.path(CACHE, 'S7I_supp_GPX3_peds_pat.rds'),
                    file.path(CACHE, 'S7N_GPX3_adult_pat.rds'),
                    file.path(OUT,   'S7N_combined_EC_GPX3_box.pdf'),
                    width = 3.6, height = 5,
                    ylab = 'GPX3 expression (per-patient mean)',
                    is_gene = TRUE, gene_label = 'GPX3')

## N3 – ACKR1 biphasic single-gene combined
.combined_panel_box(file.path(CACHE, 'S7I_supp_ACKR1_peds_pat.rds'),
                    file.path(CACHE, 'S7N_ACKR1_adult_pat.rds'),
                    file.path(OUT,   'S7N_combined_EC_ACKR1_box.pdf'),
                    width = 3.6, height = 5,
                    ylab = 'ACKR1 expression (per-patient mean)',
                    is_gene = TRUE, gene_label = 'ACKR1 (biphasic)')

###############################################################################
## ── Combined Phase-1 volcano (peds + adult) + concordance scatter ──────────
###############################################################################
if (file.exists(.cache_peds_phase1) && file.exists(.cache_adult_phase1)) {
  peds_p1  <- readRDS(.cache_peds_phase1)
  adult_p1 <- readRDS(.cache_adult_phase1)

  ## Volcano helper using phase1 categories/palette
  .ec_volcano_phase1 <- function(df, label_text, subtitle_text,
                                  fc_xlim = c(-6, 6), sig_ylim = c(0, 8),
                                  show_legend = TRUE) {
    suppressPackageStartupMessages(library(EnhancedVolcano))
    df <- df[!is.na(df$padj) & !is.na(df$log2FoldChange), ]
    rownames(df) <- make.unique(as.character(df$gene))
    df$category  <- ifelse(df$gene %in% phase1_genes,
                           phase1_lookup[df$gene], 'Other')
    keyvals <- setNames(phase1_pal[df$category], df$category)
    sel     <- intersect(phase1_genes, df$gene)
    EnhancedVolcano(df, lab = rownames(df),
      x = 'log2FoldChange', y = 'padj',
      FCcutoff = 0.5, pCutoff = 0.05,
      xlim = fc_xlim, ylim = sig_ylim,
      selectLab = sel, max.overlaps = Inf,
      drawConnectors = TRUE, arrowheads = FALSE,
      widthConnectors = 0.3, caption = '',
      title    = label_text, subtitle = subtitle_text,
      legendPosition = if (show_legend) 'right' else 'none',
      colCustom = keyvals, labFace = 'italic')
  }

  v1_peds  <- .ec_volcano_phase1(peds_p1,
                                  'Pediatric pan-EC',
                                  'ped-HLHS vs ped-NF (Phase 1)',
                                  show_legend = FALSE)
  v1_adult <- .ec_volcano_phase1(adult_p1,
                                  'Adult pan-EC',
                                  'pRV vs NF (Phase 1)',
                                  show_legend = TRUE)

  pdf(file.path(OUT, 'S7N_combined_EC_phase1_volcano.pdf'),
      width = 14, height = 6)
  tryCatch(print(cowplot::plot_grid(v1_peds, v1_adult, ncol = 2,
                                     align = 'h', axis = 'tb')),
           error = function(e) cat('phase1 volcano render error:',
                                    conditionMessage(e), '\n'))
  dev.off()

  ## Concordance scatter — peds Phase 1 log2FC vs adult Phase 1 log2FC
  peds_sub  <- peds_p1[peds_p1$gene  %in% phase1_genes,
                       c('gene','log2FoldChange','padj')]
  adult_sub <- adult_p1[adult_p1$gene %in% phase1_genes,
                         c('gene','log2FoldChange','padj')]
  conc1 <- merge(peds_sub, adult_sub, by = 'gene',
                 suffixes = c('_peds','_adult'))
  conc1$category <- factor(phase1_lookup[conc1$gene],
                           levels = setdiff(names(phase1_pal), 'Other'))
  ct1 <- suppressWarnings(cor.test(conc1$log2FoldChange_peds,
                                    conc1$log2FoldChange_adult,
                                    method = 'spearman', exact = FALSE))
  p_conc1 <- ggplot(conc1, aes(x = log2FoldChange_peds,
                                y = log2FoldChange_adult,
                                fill = category, label = gene)) +
    geom_hline(yintercept = 0, linetype = 'dashed', colour = 'grey60') +
    geom_vline(xintercept = 0, linetype = 'dashed', colour = 'grey60') +
    geom_abline(slope = 1, intercept = 0, linetype = 'dotted',
                colour = 'grey60') +
    geom_point(size = 3.4, shape = 21, colour = 'black', stroke = 0.4) +
    ggrepel::geom_text_repel(size = 2.8, fontface = 'italic',
                              max.overlaps = Inf, box.padding = 0.4,
                              segment.size = 0.3) +
    scale_fill_manual(values = phase1_pal) +
    annotate('text', x = -Inf, y = Inf,
             label = sprintf('Spearman rho = %.2f\np = %.2g',
                             ct1$estimate, ct1$p.value),
             hjust = -0.1, vjust = 1.4, size = 3) +
    xlab('Pediatric log2FC (ped-HLHS vs ped-NF)') +
    ylab('Adult log2FC (pRV vs NF)') +
    ggtitle('Phase-1 vasoprotective genes - peds vs adult sn concordance') +
    theme_classic() +
    theme(panel.border = element_rect(linewidth = 0.6, fill = NA, color = 'black'),
          legend.title = element_blank())
  pdf(file.path(OUT, 'S7N_combined_EC_phase1_concordance.pdf'),
      width = 6, height = 5)
  print(p_conc1); dev.off()
  cat('wrote Phase 1 volcano + concordance\n')
}

## O – NRG1 in Endocardial EC (pseudobulk DESeq2-VST per donor)
##   Built from .ec_subtype_pseudobulk_vst caches (different schema from
##   .peds_pb_* boxes — has expr column, no program). Use a custom builder.
{
  pp <- file.path(CACHE, 'S7_supp_peds_endocardial_NRG1_vst.rds')
  ap <- file.path(CACHE, 'S7_supp_adult_endocardial_NRG1_vst.rds')
  if (file.exists(pp) && file.exists(ap)) {
    p_df <- readRDS(pp); a_df <- readRDS(ap)
    p_df$cohort  <- 'Pediatric'; a_df$cohort <- 'Adult'
    library(patchwork)
    plot_one <- function(df, gl, pal) {
      df$group <- factor(as.character(df$group), levels = gl)
      ggplot(df, aes(x = group, y = expr, fill = group)) +
        geom_violin(scale = 'width', alpha = 0.6, linewidth = 0.4) +
        geom_jitter(width = 0.08, size = 1.6, alpha = 0.9) +
        scale_fill_manual(values = pal) +
        ggpubr::stat_compare_means(
          comparisons = list(c(gl[1],gl[2]), c(gl[2],gl[3]), c(gl[1],gl[3])),
          method = 'wilcox.test', size = 2.6) +
        labs(x = NULL, y = expression(italic('NRG1') ~ '(VST)')) +
        theme_classic() +
        theme(legend.position = 'none',
              panel.border = element_rect(linewidth = 0.6, fill = NA, color = 'black'))
    }
    p_pe <- plot_one(p_df, c('ped-NF','ped-HLHS','ped-RVF'), peds_pal_ec) +
      ggtitle('Pediatric Endocardial')
    p_ad <- plot_one(a_df, c('NF','pRV','RVF'), adult_pal_ec) +
      ggtitle('Adult Endocardial')
    pdf(file.path(OUT, 'S7O_combined_EC_NRG1_endocardial.pdf'),
        width = 6.5, height = 4)
    print(p_pe + p_ad + plot_layout(ncol = 2))
    dev.off()
  }
}

###############################################################################
## ── Combined Phase-2 volcano (peds top / adult bottom) + concordance scatter ─
###############################################################################
suppressPackageStartupMessages(library(patchwork))

if (file.exists(.cache_peds_phase2) && file.exists(.cache_adult_phase2)) {
  peds_de  <- readRDS(.cache_peds_phase2)
  adult_de <- readRDS(.cache_adult_phase2)

  v_peds  <- .ec_volcano(peds_de,
                         label_text    = 'Pediatric pan-EC',
                         fc_xlim = c(-3, 3), sig_ylim = c(0, 5),
                         show_legend   = FALSE,
                         subtitle_text = 'ped-RVF vs ped-HLHS')
  v_adult <- .ec_volcano(adult_de,
                         label_text    = 'Adult pan-EC',
                         fc_xlim = c(-3, 3), sig_ylim = c(0, 5),
                         show_legend   = TRUE,
                         subtitle_text = 'RVF vs pRV')
  pdf(file.path(OUT, 'S7H_combined_EC_phase2_volcano.pdf'),
      width = 14, height = 6)
  tryCatch(print(cowplot::plot_grid(v_peds, v_adult, ncol = 2,
                                     align = 'h', axis = 'tb')),
           error = function(e) cat('combined volcano render error:',
                                    conditionMessage(e), '\n'))
  dev.off()
  cat('wrote S7H_combined_EC_phase2_volcano.pdf\n')

  ## ── Peds-sn vs adult-sn concordance scatter for the 13 Phase-2 genes ──────
  peds_sub  <- peds_de[peds_de$gene  %in% phase2_genes,
                       c('gene','log2FoldChange','padj')]
  adult_sub <- adult_de[adult_de$gene %in% phase2_genes,
                         c('gene','log2FoldChange','padj')]
  conc <- merge(peds_sub, adult_sub, by = 'gene',
                suffixes = c('_peds','_adult'))
  conc$category <- phase2_lookup[conc$gene]
  rho <- suppressWarnings(cor(conc$log2FoldChange_peds,
                              conc$log2FoldChange_adult,
                              method = 'spearman'))
  pval <- suppressWarnings(cor.test(conc$log2FoldChange_peds,
                                     conc$log2FoldChange_adult,
                                     method = 'spearman',
                                     exact = FALSE)$p.value)

  p_conc <- ggplot(conc, aes(x = log2FoldChange_peds,
                             y = log2FoldChange_adult,
                             fill = category, label = gene)) +
    geom_hline(yintercept = 0, linetype = 'dashed', colour = 'grey60') +
    geom_vline(xintercept = 0, linetype = 'dashed', colour = 'grey60') +
    geom_abline(slope = 1, intercept = 0, linetype = 'dotted',
                colour = 'grey60') +
    geom_point(size = 3.4, shape = 21, colour = 'black', stroke = 0.4) +
    ggrepel::geom_text_repel(size = 2.6, fontface = 'italic',
                              max.overlaps = Inf,
                              box.padding = 0.4, segment.size = 0.3) +
    scale_fill_manual(values = phase2_pal) +
    annotate('text', x = -Inf, y = Inf,
             label = sprintf('Spearman rho = %.2f\np = %.2g',
                             rho, pval),
             hjust = -0.1, vjust = 1.4, size = 3) +
    xlab('Pediatric log2FC (ped-RVF vs ped-HLHS)') +
    ylab('Adult log2FC (RVF vs pRV)') +
    ggtitle('Phase-2 candidate genes - peds vs adult sn concordance') +
    theme_classic() +
    theme(panel.border = element_rect(linewidth = 0.6, fill = NA, color = 'black'),
          legend.title = element_blank())

  pdf(file.path(OUT, 'S7I_combined_EC_phase2_concordance.pdf'),
      width = 6, height = 5)
  print(p_conc); dev.off()
  cat('wrote S7I_combined_EC_phase2_concordance.pdf\n')
} else {
  cat('Skipping Phase-2 volcano/concordance - missing DE caches\n')
}

cat('\n=== Supplementary Figure 7 generation complete ===\n')
cat('Outputs in: ', OUT, '\n')

###############################################################################
## v57 standardized per-panel emission — keyed to the v57 S7 legend
## (.figure_run_logs/v57_figure_legends.md), NOT the reference's internal
## letters. v57 S7 is A–J. Reference writes S7<X>_*.pdf to OUT; copy each
## to Figure_S7_panel_<L>.pdf. v57 S7:
##   A myeloid concordance dot | B bulk-WGCNA module dot (2 stacked) |
##   C GO-BP up M8 DEGs | D ChEA up ped-HLHS-vs-ped-NF (CIITA) |
##   E ChEA up+down failing-SV (NR3C1/CIITA/IRF8) |
##   F per-patient myeloid program scores |
##   G adult-EC-module dot in HLHS EC | H EC-proportion stacked bar |
##   I Phase-1 EC vasoprotective + IFN score | J subtype-restricted
##   Phase-2 readouts (MECOM/ArtTF/Notch/SMAD1/NR2F2).
###############################################################################
.s7_v57 <- list(
  A = 'S7A_peds_myeloid_concordance_dot.pdf',
  B = 'S7B_peds_myeloid_modules_combined.pdf',
  C = 'S7C_peds_myeloid_M8_GO_up.pdf',          # may be absent → flagged
  D = 'S7D_chea_celltype_HLHS_vs_NF_up.pdf',
  E = 'S7E_chea_RVF_vs_NF_up_down.pdf',
  F = 'S7E_program_scores_box.pdf',
  G = 'S7G_peds_EC_adult_modules_combined.pdf',
  H = 'S7G_peds_EC_subtype_prev_stacked.pdf',
  I = 'S7N_combined_EC_phase1_box.pdf',
  J = 'S7J_combined_5x2.pdf')
message('Supplementary Figure 7 v57-keyed standardized panels:')
for (.L in names(.s7_v57)) {
  .src <- file.path(OUT, .s7_v57[[.L]])
  if (file.exists(.src)) {
    file.copy(.src, file.path(OUT, sprintf('Figure_S7_panel_%s.pdf', .L)),
              overwrite = TRUE)
    message('  ', .L, ': OK')
  } else message('  ', .L, ': MISSING (', .s7_v57[[.L]], ')')
}
