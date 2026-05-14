###############################################################################
## Supplementary Figure 3 (v54 draft) — Fibroblast and mural cell subclustering
##                                       and spatial validation
##
## Panels (from RV_snRNASeq_v54_draft.md figure legends):
##   FB section (A–N):
##   (A) UMAP showing subclustering of fibroblasts into 7 subpopulations
##   (B) Dot plot of FB subtype marker genes
##   (C) Per-patient FB subtype proportions (NF/pRV/RVF) shown as boxplots with
##       per-patient dots; not significantly different across disease states
##   (D) FB Phase 1 (NF→pRV) two-platform concordance heatmap
##       (snRNA-seq vs Xenium FB-lineage pseudobulk); * padj<0.05, |log2FC|>0.5,
##       baseMean>50; grey tiles = off the 477-gene Xenium panel
##   (E) FB Phase 2 (pRV→RVF) sn-vs-Xenium concordance scatter at the FB lineage
##   (F) "Erosion of FB identity" Phase 1 pseudobulk module score per patient
##       (Koenig 2022 donor-FB signature: GPX3, PID1, TGFBR3, ACSM3, APOD)
##   (G) Phase 2 matrifibrocyte commitment per-patient pseudobulk module score
##       (Fu 2018: COMP, CILP, CILP2, FMOD, CHAD, CTHRC1)
##   (H) Xenium FB subtype dot plot using literature-supported markers
##       (SpaGE-imputed expression for off-panel genes)
##   (I) Xenium ↔ snRNA FB cluster similarity heatmap (Seurat label transfer)
##   (J) FB–myeloid spatial colocalization (K=15 nearest-neighbor; n=9 patients)
##   (K) Per-patient co-abundance scatter (Xenium): FB Matrix density vs total
##       myeloid density, colored by disease group
##   (L) Spatial vignettes in patient 1632 (RVF): FB Stressed / FB Matrix +
##       iMac/Monocyte highlighted
##   (M) [continued spatial vignette / colocalization panels per legend]
##   (N) [continued; per legend]
##   Mural section (O–T):
##   (O) Mural Phase 1 (pRV vs NF) snRNA sublineage heatmap
##   (P) Mural Phase 1 vascular-protective pseudobulk module score (PC+SM)
##   (Q) Mural Phase 2 (RVF vs pRV) snRNA sublineage heatmap
##   (R) Mural Phase 2 Xenium sublineage heatmap (on-panel genes; 6 PC/VSMC subsets)
##   (S) Mural Phase 2 sn-vs-Xenium concordance scatter at the mural lineage
##   (T) Phase 2 vascular IEG/stress pseudobulk module score
##       (NR4A1/2/3 + JUN/JUNB/FOS + HSPB1 + ATF3 — NR4A vascular axis)
##
## Output: ./output/Supplementary_Figure_3/v52_figures/SupplementaryFigure_3.pdf (composite) +
##         per-panel PDFs from save_figure() calls scattered through this script.
###############################################################################

source('./helper_scripts/_shared_helpers.R')

## Per-figure output directory (introduced for consistent output paths)
V52_FIG_DIR <- './output/Supplementary_Figure_3'
dir.create(V52_FIG_DIR, showWarnings = FALSE, recursive = TRUE)


## Suppress R's default Rplots.pdf in cwd when Rscript hits a plot call
## that's outside an explicit pdf() ... dev.off() envelope.
pdf(NULL)
COMP_W <- 7.2
COMP_H <- 13   # 5-row composite (FB story spans rows 2 + 3)
PS     <- pub_scales(COMP_W)

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)
  library(stringr)
  library(ggrepel)
  library(readxl)
  # enrichR onAttach makes a network call to maayanlab.cloud which can hang.
  # Cached ChEA results in ./output/Supplementary_Figure_3/fig_s3_*_chea.rds make the live call
  # unnecessary; skip if network is unreachable.
  tryCatch(suppressWarnings(library(enrichR)),
           error = function(e) message('enrichR load failed (',
                                       conditionMessage(e),
                                       ') — using cached ChEA results only'))
})

## ───────────────────────────────────────────────────────────────────────────
## Shared DEG table + helpers
## ───────────────────────────────────────────────────────────────────────────

sn_sub <- read.csv('./dependencies/shared/sn_pseudobulk_sublineage_deseq2.csv',
                   stringsAsFactors = FALSE)

## FB subtype rename (literature-grounded; see Methods/Discussion).
## Backend cluster IDs in upstream caches are Fb1..Fb7; we relabel everywhere
## downstream so figure labels and legend agree.
fb_label_map <- c(
  'Fb1' = 'Fb_Resident',
  'Fb2' = 'Fb_Adventitial',
  'Fb3' = 'Fb_Elastogenic',
  'Fb4' = 'Fb_Interstitial',
  'Fb5' = 'Fb_Stressed',
  'Fb6' = 'Fb_Pro-fibrotic',
  'Fb7' = 'Fb_Anti-fibrotic'
)
fb_subs    <- unname(fb_label_map)
mural_subs <- c('PC', 'SM')

sn_sub$subtype <- ifelse(sn_sub$subtype %in% names(fb_label_map),
                         fb_label_map[sn_sub$subtype], sn_sub$subtype)

sig_padj   <- 0.05
sig_lfc    <- 0.5
sig_bmean  <- 50

## Build a long gene x subtype frame for a given contrast.
.prep_phase_heatmap <- function(df, subs, contrast_label, genes) {
  m <- df %>%
    dplyr::filter(subtype %in% subs, contrast == contrast_label,
                  gene %in% genes) %>%
    dplyr::select(gene, subtype, log2FoldChange, padj, baseMean)

  grid <- expand.grid(gene = genes, subtype = subs,
                      stringsAsFactors = FALSE)
  res  <- dplyr::left_join(grid, m, by = c('gene', 'subtype'))

  res$sig <- with(res,
    ifelse(!is.na(padj) & padj < sig_padj &
             abs(log2FoldChange) > sig_lfc & baseMean > sig_bmean,
           '*', ''))
  res$gene    <- factor(res$gene, levels = rev(genes))
  res$subtype <- factor(res$subtype, levels = subs)
  res
}

.make_phase_heatmap <- function(m, title = NULL, palette_lim = 3) {
  ggplot(m, aes(x = subtype, y = gene, fill = log2FoldChange)) +
    geom_tile(color = 'grey85', linewidth = PS$linewidth_mm * 0.3) +
    geom_text(aes(label = sig), color = 'black',
              size = PS$text_mm, vjust = 0.75) +
    scale_fill_gradient2(low = '#2c5fa5', mid = 'white', high = '#b3263e',
                         midpoint = 0,
                         limits = c(-palette_lim, palette_lim),
                         oob = scales::squish, na.value = 'grey92',
                         name = expression(log[2]~FC)) +
    labs(x = NULL, y = NULL, title = title) +
    theme_v52(COMP_W) +
    theme(
      axis.text.x       = element_text(angle = 45, hjust = 1),
      axis.text.y       = element_text(face = 'italic'),
      axis.line         = element_blank(),
      axis.ticks        = element_blank(),
      legend.key.height = unit(3.5, 'mm'),
      legend.key.width  = unit(1.5, 'mm'),
      plot.title        = element_text(face = 'bold')
    )
}

.make_chea_bar <- function(df, fill, title) {
  ggplot(df,
         aes(x = reorder(TF, -log10(Adjusted.P.value)),
             y = -log10(Adjusted.P.value))) +
    geom_col(fill = fill, width = 0.75) +
    geom_hline(yintercept = -log10(0.05), linetype = 'dashed',
               color = 'grey50',
               linewidth = PS$linewidth_mm * 0.6) +
    coord_flip() +
    labs(x = NULL,
         y = expression(-log[10]~adj.~italic(p)),
         title = title) +
    theme_v52(COMP_W) +
    theme(plot.title  = element_text(face = 'bold'),
          axis.text.y = element_text(face = 'bold'))
}

## ───────────────────────────────────────────────────────────────────────────
## Panel A: FB UMAP (Fb1–Fb7)
## ───────────────────────────────────────────────────────────────────────────

M_fb <- readRDS('./dependencies/shared/fb_subclust.rds')
if (!('Subnames' %in% colnames(M_fb@meta.data))) {
  new.cluster.ids <- c('Fb1','Fb2','Fb3','Fb4','Fb5','Fb6','Fb7')
  names(new.cluster.ids) <- levels(M_fb)
  M_fb <- RenameIdents(M_fb, new.cluster.ids)
  M_fb$Subnames <- M_fb@active.ident
}
## Apply literature-grounded labels (Fb1..Fb7 -> Fb{n}_Function)
new_subnames <- fb_label_map[as.character(M_fb@meta.data$Subnames)]
names(new_subnames) <- rownames(M_fb@meta.data)
M_fb@meta.data$Subnames <- factor(new_subnames, levels = unname(fb_label_map))
Idents(M_fb) <- 'Subnames'

## ── FB pseudobulk module scores (Koenig identity + Fu matrifibrocyte) ────
##
## Koenig 2022 donor-FB signature, verbatim (5 genes):
##   "Donor fibroblasts selectively expressed a signature represented by
##    GPX3, PID1, TGFBR3, ACSM3 and APOD (Fig. 6f)."
##   — Koenig et al. Nat Cardiovasc Res 2022, DOI 10.1038/s44161-022-00028-6
## → expected to *erode* in Phase 1 (NF→pRV).
##
## Fu/Molkentin 2018 matrifibrocyte signature (3 genes, main-text + IHC):
##   "Selected genes from this pathway signature at 2 and 4 weeks after MI
##    showed an induction of cartilage oligomeric matrix protein (Comp,
##    also known as thrombospondin 5) and chondroadherin (Chad) and Cilp2
##    (Figure 10G). IHC staining confirmed the specific expression of Comp
##    and Chad in the infarct region… Comp and Chad expression were also
##    identified in the fibrotic area of human heart samples from ischemic
##    patients" — Fu et al. J Clin Invest 2018, PMID 29664017, p. 2137,
##   Fig 10G + Fig 11C (human IHC validation).
## → expected to *commit* in Phase 2 (pRV→RVF).
## (Prior 6-gene panel including CILP/FMOD/CTHRC1 was not supported by the
##  paper's main text; CTHRC1 specifically belongs to a distinct activated-
##  profibrotic state per Ruiz-Villalba et al. Circulation 2020, PMID 32972203.)

# Aggregate FB counts per patient (lineage level), CPM + log1p once
.fb_agg    <- AggregateExpression(M_fb, assays = 'RNA',
                                  group.by = 'patient',
                                  slot = 'counts',
                                  return.seurat = FALSE)$RNA
.lib_sizes <- colSums(.fb_agg)
.fb_logcpm <- log1p(sweep(.fb_agg, 2, .lib_sizes, '/') * 1e6)
.pat_grp <- M_fb@meta.data %>%
  dplyr::distinct(patient, group) %>%
  dplyr::mutate(patient = as.character(patient),
                group   = factor(group, levels = c('NF','pRV','RVF')))

.fb_pb_module <- function(gene_set, label_short) {
  have <- intersect(gene_set, rownames(.fb_logcpm))
  cat(sprintf('%s — detected: %s\n', label_short, paste(have, collapse = ', ')))
  mat <- .fb_logcpm[have, , drop = FALSE]
  z   <- t(scale(t(mat)))
  d <- data.frame(patient = sub('^g', '', colnames(z)),
                  score   = colMeans(z, na.rm = TRUE),
                  stringsAsFactors = FALSE)
  dplyr::left_join(d, .pat_grp, by = 'patient') %>%
    dplyr::arrange(group, patient)
}

fb_identity_genes <- c('GPX3','PID1','TGFBR3','ACSM3','APOD')
fb_matrifib_genes <- c('COMP','CHAD','CILP2')

fb_pb_score    <- .fb_pb_module(fb_identity_genes, 'FB identity (Koenig 2022)')
## ── Panels F + G — FB identity-loss + matrifibrocyte-commitment pseudobulk module scores
##    Panel F: 'Erosion of FB identity' Phase 1 (GPX3, ELN, TIMP3, SFRP4, ...)
##    Panel G: 'Matrifibrocyte commitment' Phase 2 (Fu et al. 2018 gene set)
fb_pb_matrifib <- .fb_pb_module(fb_matrifib_genes, 'Matrifibrocyte (Fu 2018)')

cat('\n=== FB-identity score (Koenig) ===\n');     print(fb_pb_score)
cat('\n=== Matrifibrocyte score (Fu 2018) ===\n'); print(fb_pb_matrifib)

p_A <- DimPlot(M_fb, group.by = 'Subnames', label = TRUE,
               label.size = PS$base_pt / .pt, repel = TRUE,
               pt.size = PS$umap_pt, raster = TRUE) +
  NoLegend() +
  labs(title = 'Fibroblast subtypes') +
  theme_v52(COMP_W) +
  theme(plot.title = element_text(face = 'bold'),
        axis.text  = element_blank(),
        axis.ticks = element_blank(),
        axis.line  = element_blank())

## ───────────────────────────────────────────────────────────────────────────
## Panel B: FB dotplot — functional / novel markers (not generic Fb-scores)
## ───────────────────────────────────────────────────────────────────────────

fb_marker_groups <- list(
  'Fb_Resident'      = c('SCN7A','ACSM3'),
  'Fb_Adventitial'   = c('C7','CNTNAP2','KAZN'),
  'Fb_Elastogenic'   = c('NOX4','ELN','FGF14'),
  'Fb_Interstitial'  = c('FBN1','MFAP5','PCOLCE2'),
  'Fb_Stressed'      = c('NR4A1','NR4A3','NAMPT'),
  'Fb_Pro-fibrotic'  = c('TNC','FN1','SERPINE1'),
  'Fb_Anti-fibrotic' = c('TIMP3','SYN3')
)
fb_marker_vec <- unlist(fb_marker_groups, use.names = FALSE)
fb_marker_vec <- fb_marker_vec[fb_marker_vec %in% rownames(M_fb)]

p_B <- DotPlot(M_fb, features = fb_marker_vec, assay = 'SCT',
               cluster.idents = FALSE,
               cols = c('lightgrey','blue'),
               col.min = 0, col.max = 2,
               dot.scale = PS$dot_range[2]) +
  scale_size_continuous(range = PS$dot_range) +
  labs(x = NULL, y = NULL, title = 'FB subtype markers') +
  theme_v52(COMP_W) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = 'italic'),
    plot.title  = element_text(face = 'bold'),
    legend.key.height = unit(3.5,'mm'),
    legend.key.width  = unit(1.5,'mm')
  )

## ───────────────────────────────────────────────────────────────────────────
## Panel D: FB Phase 1 (NF→pRV) — sn vs Xenium 2-column concordance heatmap
##           (script-letter C → v57-letter D; previously labeled Panel C in v54)
##   Erosion of non-failing FB identity (GPX3, ELN, TIMP3, SFRP4) and
##   sn-only off-Xenium-panel anti-fibrotic loss + activation priming
##   (TGFBR3, ACTA2, THBS4). Frame from Koenig 2022 donor-FB signature.
## ───────────────────────────────────────────────────────────────────────────

xen_lin     <- read.csv('./dependencies/shared/Xenium/xenium_pseudobulk_lineage_deseq2.csv',
                        stringsAsFactors = FALSE)
sn_lin      <- read.csv('./dependencies/shared/sn_pseudobulk_lineage_deseq2.csv',
                        stringsAsFactors = FALSE)
xen_sub_csv <- read.csv('./dependencies/shared/Xenium/xenium_pseudobulk_sublineage_deseq2.csv',
                        stringsAsFactors = FALSE)
## Apply unified Xenium FB naming to the sublineage CSV (same map as for
## the spatial colocalisation analysis below)
## Apply only the FB_Matrifibrocyte → FB_Resting rename to the sublineage CSV
.xen_fb_rename_map_csv <- c(
  'FB_Stress'         = 'Xen_FB_Stressed',
  'FB_NOX4'           = 'Xen_FB_NOX4',
  'FB_Homeostatic'    = 'Xen_FB_Matrix',
  'FB_NTN4'           = 'Xen_FB_Interstitial',
  'FB_PCOLCE2'        = 'Xen_FB_Intermediate',
  'FB_Activated'      = 'Xen_FB_Activated',
  'FB_Matrifibrocyte' = 'Xen_FB_CILP')
.is_fb_csv <- xen_sub_csv$subtype %in% names(.xen_fb_rename_map_csv)
xen_sub_csv$subtype[.is_fb_csv] <- .xen_fb_rename_map_csv[xen_sub_csv$subtype[.is_fb_csv]]

.fb_phase_concord <- function(genes, contrast_label) {
  sn_d <- sn_lin %>%
    dplyr::filter(subtype == 'FB', contrast == contrast_label,
                  gene %in% genes) %>%
    dplyr::transmute(gene, platform = 'snRNA',
                     log2FoldChange, padj, baseMean)
  xen_d <- xen_lin %>%
    dplyr::filter(subtype == 'FB', contrast == contrast_label,
                  gene %in% genes) %>%
    dplyr::transmute(gene, platform = 'Xenium',
                     log2FoldChange, padj, baseMean)
  d <- dplyr::bind_rows(sn_d, xen_d)
  grid <- expand.grid(gene = genes,
                      platform = c('snRNA','Xenium'),
                      stringsAsFactors = FALSE)
  d <- dplyr::left_join(grid, d, by = c('gene','platform'))
  d$sig <- ifelse(!is.na(d$padj) & d$padj < sig_padj &
                  abs(d$log2FoldChange) > sig_lfc &
                  !is.na(d$baseMean) & d$baseMean > sig_bmean,
                  '*', '')
  d$gene     <- factor(d$gene, levels = rev(genes))
  d$platform <- factor(d$platform, levels = c('snRNA','Xenium'))
  d
}

## Erosion (loss) genes first, then activation/priming (gain) genes
fb_p1_genes   <- c('GPX3','ELN','TIMP3','TGFBR3','SFRP4','ACTA2','THBS4')
fb_p1_concord <- .fb_phase_concord(fb_p1_genes, 'pRV_vs_NF')

p_C <- ggplot(fb_p1_concord,
              aes(x = platform, y = gene, fill = log2FoldChange)) +
  geom_tile(color = 'grey85', linewidth = PS$linewidth_mm * 0.3) +
  geom_text(aes(label = sig), color = 'black',
            size = PS$text_mm, vjust = 0.75) +
  scale_fill_gradient2(low = '#3B4CC0', mid = 'white', high = '#B40426',
                       midpoint = 0, limits = c(-4, 4),
                       oob = scales::squish, na.value = 'grey92',
                       name = expression(log[2]~FC)) +
  theme_v52(COMP_W) +
  theme(plot.title  = element_text(face = 'bold'),
        axis.text.y = element_text(face = 'italic'),
        legend.key.height = unit(3.5, 'mm'),
        legend.key.width  = unit(1.5, 'mm')) +
  xlab(NULL) + ylab(NULL) +
  labs(title = 'FB Phase 1 (NF→pRV)\nerosion + priming')

## ───────────────────────────────────────────────────────────────────────────
## Preview panels: paired pseudobulk module-score plots
##   (Koenig 2022 identity erosion — Phase-1 specific)
##   (Fu 2018 matrifibrocyte commitment — Phase-2 specific)
## ───────────────────────────────────────────────────────────────────────────

.fb_module_plot <- function(d, ylab_text, title_text, csv_tag) {
  kw_p <- tryCatch(kruskal.test(score ~ group, data = d)$p.value,
                   error = function(e) NA_real_)
  pairs <- list(c('NF','pRV'), c('NF','RVF'), c('pRV','RVF'))
  pw <- do.call(rbind, lapply(pairs, function(p) {
    a <- d$score[d$group == p[1]]; b <- d$score[d$group == p[2]]
    pv <- tryCatch(suppressWarnings(wilcox.test(a, b, exact = FALSE)$p.value),
                   error = function(e) NA_real_)
    data.frame(g1 = p[1], g2 = p[2], wilcox_p = pv,
               stringsAsFactors = FALSE)
  }))
  pw$padj <- p.adjust(pw$wilcox_p, method = 'BH')

  cat(sprintf('\n=== %s — KW p = %s ===\n', csv_tag, signif(kw_p, 3)))
  print(pw %>% dplyr::mutate(wilcox_p = signif(wilcox_p, 3),
                             padj     = signif(padj, 3)))
  write.csv(d,  sprintf('./output/Supplementary_Figure_3/fig_s3_%s_score.csv', csv_tag),    row.names = FALSE)
  write.csv(pw, sprintf('./output/Supplementary_Figure_3/fig_s3_%s_pairwise.csv', csv_tag), row.names = FALSE)

  y_max  <- max(d$score, na.rm = TRUE)
  y_step <- (max(d$score, na.rm = TRUE) - min(d$score, na.rm = TRUE)) * 0.18
  brk <- pw %>%
    dplyr::mutate(x    = match(g1, c('NF','pRV','RVF')),
                  xend = match(g2, c('NF','pRV','RVF')),
                  y    = y_max + y_step * c(1, 3, 2),
                  lab  = dplyr::case_when(
                    is.na(wilcox_p)  ~ 'NA',
                    wilcox_p < 0.001 ~ '***',
                    wilcox_p < 0.01  ~ '**',
                    wilcox_p < 0.05  ~ '*',
                    TRUE             ~ sprintf('p=%s', signif(wilcox_p, 2))))

  ggplot(d, aes(x = group, y = score, fill = group)) +
    geom_boxplot(outlier.shape = NA, width = 0.55,
                 linewidth = PS$geom_lw, alpha = 0.85) +
    geom_jitter(width = 0.12, size = 1.4, shape = 21,
                color = 'black', stroke = 0.3) +
    scale_fill_manual(values = disease_pal, guide = 'none') +
    geom_segment(data = brk, aes(x = x, xend = xend, y = y, yend = y),
                 inherit.aes = FALSE, linewidth = PS$geom_lw) +
    geom_segment(data = brk, aes(x = x, xend = x, y = y, yend = y - y_step * 0.15),
                 inherit.aes = FALSE, linewidth = PS$geom_lw) +
    geom_segment(data = brk, aes(x = xend, xend = xend, y = y, yend = y - y_step * 0.15),
                 inherit.aes = FALSE, linewidth = PS$geom_lw) +
    geom_text(data = brk,
              aes(x = (x + xend) / 2, y = y + y_step * 0.10, label = lab),
              inherit.aes = FALSE, size = PS$text_mm,
              family = FONT_FAMILY, vjust = 0) +
    annotate('text', x = 0.55, y = y_max + y_step * 4.0,
             label = sprintf('KW p = %s', signif(kw_p, 2)),
             hjust = 0, size = PS$text_mm,
             family = FONT_FAMILY, fontface = 'italic') +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.30))) +
    theme_v52(COMP_W) +
    theme(plot.title = element_text(face = 'bold')) +
    xlab(NULL) + ylab(ylab_text) +
    labs(title = title_text)
}

p_module_identity <- .fb_module_plot(
  fb_pb_score,
  ylab_text  = 'FB-identity score (z, pseudobulk)',
  title_text = 'Phase 1: identity erosion\n(Koenig 2022 donor signature)',
  csv_tag    = 'fb_identity_koenig_pseudobulk')

p_module_matrifib <- .fb_module_plot(
  fb_pb_matrifib,
  ylab_text  = 'Matrifibrocyte score (z, pseudobulk)',
  title_text = 'Phase 2: matrifibrocyte commitment\n(Fu 2018 PMID 29664017)',
  csv_tag    = 'fb_matrifib_fu_pseudobulk')

# Save both individually + as a paired patchwork preview
save_figure(p_module_identity,
            'SupplementaryFigure_3_preview_fb_identity_pseudobulk.pdf',
            width = COMP_W * 1.0 / 2.9, height = COMP_H / 4)
save_figure(p_module_matrifib,
            'SupplementaryFigure_3_preview_fb_matrifib_pseudobulk.pdf',
            width = COMP_W * 1.0 / 2.9, height = COMP_H / 4)
save_figure((p_module_identity | p_module_matrifib),
            'SupplementaryFigure_3_preview_fb_modules_paired.pdf',
            width = COMP_W * 2.0 / 2.9, height = COMP_H / 4)

## ───────────────────────────────────────────────────────────────────────────
## Preview panels: Xenium spatial plots of FB and mural subtypes
##   Uses the PARENT xenium_obj_subclustered.rds (~625K cells) and the
##   `cell_types_manual` column (56 manual annotations). FB and mural cells
##   are filtered out by prefix. Patients tiled within each disease group
##   via offset_coords (adapted from XeniumFigureGen lineage_spatial_map).
## ───────────────────────────────────────────────────────────────────────────

.xen_parent_path <- './dependencies/shared/Xenium/xenium_obj_subclustered.rds'

.offset_coords <- function(meta, gap_um = 500, row_height_um = 12000) {
  for (g in unique(meta$group)) {
    idx_g    <- which(meta$group == g)
    patients <- unique(meta$patient[idx_g])
    for (p in patients) {
      idx_p <- which(meta$group == g & meta$patient == p)
      meta$x[idx_p] <- meta$x[idx_p] - min(meta$x[idx_p])
      meta$y[idx_p] <- meta$y[idx_p] - min(meta$y[idx_p])
    }
    if (length(patients) <= 1) next
    bboxes <- lapply(patients, function(p) {
      idx <- which(meta$group == g & meta$patient == p)
      list(w = diff(range(meta$x[idx])), h = diff(range(meta$y[idx])))
    })
    names(bboxes) <- patients
    areas    <- sapply(bboxes, function(b) b$w * b$h)
    patients <- patients[order(-areas)]
    x_off <- 0
    for (p in patients) {
      idx_p <- which(meta$group == g & meta$patient == p)
      meta$x[idx_p] <- meta$x[idx_p] + x_off
      x_off <- x_off + bboxes[[p]]$w + gap_um
    }
  }
  meta
}

# Pull metadata once, keep just the columns we need, free the full object.
# We pull from `cell_types_subclustering` — the stricter-labeled sublineage
# annotation (~290K labeled, rest 'Unassigned'). FB/Myeloid/Mural/EC label
# strings are identical in shape to cell_types_manual but reflect the
# sublineage clustering pass only (e.g., FB labeled = 28,687 vs 62,192 in
# the manual column).
.xen_parent <- readRDS(.xen_parent_path)
.xen_meta_full <- data.frame(
  x        = .xen_parent$x_centroid,
  y        = -.xen_parent$y_centroid,
  group    = factor(.xen_parent$group, levels = c('NF','pRV','RVF')),
  patient  = .xen_parent$patient,
  subtype  = .xen_parent@meta.data[['cell_types_subclustering']],
  stringsAsFactors = FALSE
)
rm(.xen_parent); invisible(gc(verbose = FALSE))

## Keep Xenium FB names DISTINCT from snRNA Fb names (per user direction).
## Only one substantive change: FB_Matrifibrocyte → FB_Resting (the cluster
## was originally labeled matrifibrocyte but lacks canonical Fu 2018 markers
## (COMP, CHAD, THBS4) and is enriched in NF rather than scar; better matches
## the CILP+/HGF+ resting state per Koenig 2022, Amrute 2024).
## Xenium FB clusters — renamed with Xen_ prefix + literature-aligned biology
## (kept distinct from snRNA Fb_ names; cross-platform similarity in heatmap):
##   FB_Activated      → Xen_FB_Activated     (Koenig Fb8, POSTN+)
##   FB_NOX4           → Xen_FB_NOX4          (Cucoranu 2005, NOX4+ profibrotic)
##   FB_NTN4           → Xen_FB_Interstitial  (Koenig Fb2 / Litvinukova FB1;
##                                             PCOLCE2/MFAP5/FBN1 actually peak here)
##   FB_PCOLCE2        → Xen_FB_Intermediate  (weakly defined; metabolic markers)
##   FB_Matrifibrocyte → Xen_FB_CILP          (CILP+/HGF+ resting per
##                                             Koenig 2022, Amrute 2024)
##   FB_Stress         → Xen_FB_Stressed      (Farbehi 2019, Koenig Fb9)
##   FB_Homeostatic    → Xen_FB_Matrix       (Koenig Fb1, Litvinukova FB2;
##                                             DCN+/FBLN1+/ELN+ baseline)
xen_fb_rename_map <- c(
  'FB_Stress'         = 'Xen_FB_Stressed',
  'FB_NOX4'           = 'Xen_FB_NOX4',
  'FB_Homeostatic'    = 'Xen_FB_Matrix',
  'FB_NTN4'           = 'Xen_FB_Interstitial',
  'FB_PCOLCE2'        = 'Xen_FB_Intermediate',
  'FB_Activated'      = 'Xen_FB_Activated',
  'FB_Matrifibrocyte' = 'Xen_FB_CILP'
)
.is_fb <- .xen_meta_full$subtype %in% names(xen_fb_rename_map)
.xen_meta_full$subtype[.is_fb] <- xen_fb_rename_map[.xen_meta_full$subtype[.is_fb]]

.spatial_panel_from_meta <- function(meta, title_str, pal = NULL) {
  meta <- meta[!is.na(meta$subtype) & meta$subtype != 'Unassigned', ]
  meta <- .offset_coords(meta)
  p <- ggplot(meta, aes(x = x, y = y, color = subtype)) +
    ggrastr::rasterise(geom_point(size = 0.18, stroke = 0), dpi = 300) +
    facet_wrap(~ group, nrow = 1) +
    coord_fixed() +
    theme_v52(COMP_W) +
    theme(plot.title      = element_text(face = 'bold'),
          axis.text       = element_blank(),
          axis.ticks      = element_blank(),
          axis.line       = element_blank(),
          axis.title      = element_blank(),
          panel.grid      = element_blank(),
          panel.spacing   = unit(2, 'mm'),
          strip.text      = element_text(face = 'bold'),
          legend.key.size = unit(2.5, 'mm'),
          legend.position = 'right') +
    labs(title = title_str, color = NULL) +
    guides(color = guide_legend(override.aes = list(size = 2)))
  if (!is.null(pal)) p <- p + scale_color_manual(values = pal, na.translate = FALSE)
  p
}

# FB cells = cell_types_manual starting with FB_
.fb_meta    <- .xen_meta_full[grepl('^Xen_FB_', .xen_meta_full$subtype), ]
# Mural cells = Pericyte* + VSMC*
.mural_meta <- .xen_meta_full[grepl('^Pericyte|^VSMC', .xen_meta_full$subtype), ]

cat('\n=== Xenium spatial filter counts (cell_types_manual) ===\n')
cat(sprintf('FB cells:    %d (%d subtypes)\n',
            nrow(.fb_meta), length(unique(.fb_meta$subtype))))
cat(sprintf('Mural cells: %d (%d subtypes)\n',
            nrow(.mural_meta), length(unique(.mural_meta$subtype))))

p_fb_spatial <- .spatial_panel_from_meta(
  .fb_meta,
  'FB subtypes (Xenium spatial, NF / pRV / RVF)')

p_mural_spatial <- .spatial_panel_from_meta(
  .mural_meta,
  'Mural subtypes (Xenium spatial, NF / pRV / RVF)')

save_figure(p_fb_spatial,
            'SupplementaryFigure_3_preview_fb_spatial.pdf',
            width = COMP_W, height = COMP_H / 4)
save_figure(p_mural_spatial,
            'SupplementaryFigure_3_preview_mural_spatial.pdf',
            width = COMP_W, height = COMP_H / 4)

## ───────────────────────────────────────────────────────────────────────────
## ── Panel H — Xenium FB subtype dot plot (literature-supported markers)
## Preview: Xenium FB dotplot — Xenium-native cluster names + literature-
## supported markers per cluster (cross-referenced against the canonical
## fb.clean.marks.csv FindAllMarkers output and the supporting literature:
##   FB_Homeostatic — Litvinukova 2020, Koenig 2022, Chaffin 2022 (resting)
##   FB_Activated   — Kanisicak 2016 (POSTN-Cre), Koenig 2022 (FB-Activated)
##   FB_NOX4        — Cucoranu 2005, Kuroda 2010, Zhao 2015 (NOX4+ profibrotic)
##   FB_Resting     — CILP+/HGF+ (Koenig 2022, Amrute 2024); replaces former
##                    FB_Matrifibrocyte label (lacks Fu 2018 markers)
##   FB_PCOLCE2     — Koenig 2022 (PDPN+/F13A1+ interstitial)
##   FB_Stress      — Farbehi 2019, Koenig 2022 (IEG + hypoxia/glycolysis)
##   FB_NTN4        — putative novel; HAS1+ in Tsukui 2020 / Tabib 2021
## SpaGE-imputed values used (some markers off the 477 panel).
## ───────────────────────────────────────────────────────────────────────────

fb_xen <- readRDS('./dependencies/shared/Xenium/fb_clean_clean.rds')
.is_fb_xen <- as.character(fb_xen$fb_subtype) %in% names(xen_fb_rename_map)
.new_fb_xen <- as.character(fb_xen$fb_subtype)
.new_fb_xen[.is_fb_xen] <- unname(xen_fb_rename_map[.new_fb_xen[.is_fb_xen]])
names(.new_fb_xen) <- colnames(fb_xen)
fb_xen$fb_subtype_unified <- factor(unname(.new_fb_xen),
  levels = c('Xen_FB_Activated','Xen_FB_NOX4','Xen_FB_Interstitial',
             'Xen_FB_Intermediate','Xen_FB_CILP','Xen_FB_Stressed',
             'Xen_FB_Matrix'))
Idents(fb_xen) <- 'fb_subtype_unified'
DefaultAssay(fb_xen) <- 'Xenium'
## Marker groups keyed by Xenium-native cluster name and chosen by direct
## expression dominance audit (rather than the original cluster-naming
## intent). PCOLCE2/MFAP5/FBN1 actually peak in FB_NTN4, not FB_PCOLCE2.
## FB_PCOLCE2 lacks a clean specific marker set — represented by the
## strongest Δpct hits (Wilcoxon FindAllMarkers, contam-filtered).
## Slimmed marker set: NOX4 and Activated share the activated program so
## NOX4 row carries NOX4+THBS4 only; CILP row uses CILP+HGF (Koenig 2022);
## Interstitial uses canonical PCOLCE2/MFAP5/FBN1+NTN4.
fb_marker_groups_xen <- list(
  'Xen_FB_Activated'    = c('POSTN','FAP','TGFB2','FNDC1','AEBP1','COL16A1',
                            'EDNRA','HAPLN1'),
  'Xen_FB_NOX4'         = c('NOX4','THBS4'),
  'Xen_FB_Interstitial' = c('NTN4','PCOLCE2','MFAP5','FBN1'),
  'Xen_FB_Intermediate' = c('IP6K2','UGP2','PDK4'),
  'Xen_FB_CILP'         = c('CILP','HGF'),
  'Xen_FB_Stressed'     = c('CXCL8','FOSL1','MT2A','NAMPT'),
  'Xen_FB_Matrix'      = c('DCN','FBLN1','ELN','SFRP4','ASPN','DPT','OGN')
)
fb_marker_vec_xen2 <- unique(unlist(fb_marker_groups_xen, use.names = FALSE))
fb_marker_vec_xen2 <- intersect(fb_marker_vec_xen2, rownames(fb_xen))

p_xen_fb_dot <- DotPlot(fb_xen, features = fb_marker_vec_xen2,
                        assay = 'Xenium',
                        cluster.idents = FALSE,
                        cols = c('lightgrey','blue'),
                        col.min = 0, col.max = 2,
                        dot.scale = PS$dot_range[2]) +
  scale_size_continuous(range = PS$dot_range) +
  coord_flip() +    # transpose: genes → Y, clusters → X
  labs(x = NULL, y = NULL, title = NULL) +
  theme_v52(COMP_W) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   face = 'plain'),
        axis.text.y = element_text(angle = 0, face = 'plain'),
        plot.title  = element_blank(),
        legend.key.height = unit(3.5,'mm'),
        legend.key.width  = unit(1.5,'mm'))

save_figure(p_xen_fb_dot,
            'SupplementaryFigure_3_preview_xen_fb_dotplot.pdf',
            width  = COMP_H / 4 * 0.67 * 1.15 * 1.10 * 1.10,  # +10% width
            height = COMP_W * 1.6 / 2.6 * 1.7)

## ───────────────────────────────────────────────────────────────────────────
## Cross-platform similarity: Xenium FB clusters × snRNA Fb clusters
##   Uses Seurat label-transfer predicted.id (Fb1..Fb6/7) carried in
##   fb_clean_clean.rds. Each tile = % of Xenium cells in row that mapped
##   to snRNA cluster in column. Annotated with snRNA renamed labels.
## ───────────────────────────────────────────────────────────────────────────

if ('predicted.id' %in% colnames(fb_xen@meta.data)) {
  sn_label_map_xen <- c('Fb1'='Fb_Resident','Fb2'='Fb_Adventitial',
                        'Fb3'='Fb_Elastogenic','Fb4'='Fb_Interstitial',
                        'Fb5'='Fb_Stressed','Fb6'='Fb_Pro-fibrotic',
                        'Fb7'='Fb_Anti-fibrotic')
  pred_renamed <- sn_label_map_xen[as.character(fb_xen$predicted.id)]
  pred_renamed[is.na(pred_renamed)] <- as.character(fb_xen$predicted.id)[is.na(pred_renamed)]

  ct <- table(as.character(fb_xen$fb_subtype_unified), pred_renamed)
  prop_row <- prop.table(ct, margin = 1) * 100
  cat('\n=== Xenium × snRNA label-transfer cross-tab (% of Xenium row → snRNA col) ===\n')
  print(round(prop_row, 1))

  # Long-form for ggplot
  sim_df <- as.data.frame(prop_row, stringsAsFactors = FALSE)
  colnames(sim_df) <- c('xenium', 'snRNA', 'pct')
  # Fb_Anti-fibrotic (snRNA Fb7, TIMP3+/SYN3+) was too rare in the snRNA
  # reference (n=164 cells) to be reliably assigned by Seurat label transfer;
  # excluded from the column axis and noted in the panel title.
  # Order Xenium rows + snRNA cols sensibly
  sim_df$xenium <- factor(sim_df$xenium,
    levels = rev(c('Xen_FB_Activated','Xen_FB_NOX4','Xen_FB_Interstitial',
                   'Xen_FB_Intermediate','Xen_FB_CILP','Xen_FB_Stressed',
                   'Xen_FB_Matrix')))
  sim_df$snRNA <- factor(sim_df$snRNA,
    levels = c('Fb_Resident','Fb_Adventitial','Fb_Elastogenic',
               'Fb_Interstitial','Fb_Stressed','Fb_Pro-fibrotic'))

  p_xen_sn_sim <- ggplot(sim_df, aes(x = snRNA, y = xenium, fill = pct)) +
    geom_tile(color = 'grey85', linewidth = PS$linewidth_mm * 0.3) +
    geom_text(aes(label = ifelse(pct >= 5, sprintf('%.0f', pct), '')),
              color = 'black', size = PS$text_mm) +
    scale_fill_gradient(low = 'white', high = '#B40426',
                        limits = c(0, 100), name = '% of Xenium\ncells (row)') +
    theme_v52(COMP_W) +
    theme(plot.title  = element_text(face = 'bold'),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.key.height = unit(3.5,'mm'),
          legend.key.width  = unit(1.5,'mm')) +
    xlab('snRNA Fb cluster (label-transfer predicted)') +
    ylab('Xenium FB cluster') +
    ## ── Panel I — Xenium ↔ snRNA FB cluster similarity heatmap (label transfer)
labs(title = paste0('Xenium ↔ snRNA FB cluster similarity (Seurat label transfer)\n',
                        'Fb_Anti-fibrotic (Fb7) excluded — too rare for accurate label transfer'))

  save_figure(p_xen_sn_sim,
              'SupplementaryFigure_3_preview_xen_sn_fb_similarity.pdf',
              width = COMP_W * 0.8, height = COMP_H / 4)
} else {
  message('predicted.id column not found in fb_xen — skipping similarity heatmap')
}

rm(fb_xen); invisible(gc(verbose = FALSE))

## ───────────────────────────────────────────────────────────────────────────
## Preview: FB–Myeloid spatial colocalization analysis
## ── Panel J — FB–Myeloid spatial colocalization (K=15 nearest-neighbor)
##   (a) Per-FB-subtype K=15-NN myeloid enrichment (observed/expected),
##       per-patient dots + boxplot + reference line at 1.0 + Wilcoxon
##       vs 1.0
##   (b) Spatial vignette in the Xen_FB_Stressed-richest tissue: Stressed
##       in red, inflammatory myeloid (Mac_Inflammatory + Monocyte) in blue
##   (c) Companion vignette for Xen_FB_Matrix in same tissue
##
## ── Biology rationale: why Xen_FB_Matrix colocalises with myeloid ──────────
## Xen_FB_Matrix (DCN+/FBLN1+/ELN+/SFRP4+ baseline ECM-producing FB) shows
## the strongest myeloid-neighbourhood enrichment of any non-stressed cluster
## (1.90×, padj 0.004 by K=15-NN; this section), AND a per-patient density
## correlation with total myeloid cells of Pearson r = 0.93 (p = 3e-04, n=9
## patients) — the strongest co-abundance among any FB cluster (next: 0.69).
## Three converging mechanisms support this pairing:
##   1. Adventitial-perivascular niche — DCN+/FBLN1+ baseline FBs occupy
##      the same peri-vessel ECM corridors immune cells use for trafficking
##      (Buechler 2021 Nature Pi16+/Col15a1+ universal-FB atlas; Koenig 2022;
##      Litvinukova 2020 spatial validation of FB1/FB2 placement).
##   2. Immunomodulatory ECM — DCN sequesters TGFβ and binds TLR2/4 on
##      macrophages; FBLN1 binds C1q complement; biglycan activates NLRP3
##      (Iozzo "matrix immunology" reviews). Matrix this cluster secretes
##      is itself immune-active.
##   3. Reactive-baseline state — FOS/FOSB/SFRP4 are elevated, indicating
##      tonic responsiveness to cytokine signals — consistent with sitting
##      in a chemokine-rich neighbourhood and reciprocally signalling back.
## Per-patient density correlation tables:
##   ./output/Supplementary_Figure_3/fig_s3_xen_fb_myel_per_patient_density.csv
##   ./output/Supplementary_Figure_3/fig_s3_xen_fb_myel_density_correlations.csv
## ───────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({ library(RANN) })

# Use .xen_meta_full (already loaded above for FB/mural spatial panels) and
# enrich with lineage labels for KNN co-occurrence.
.coloc_meta <- .xen_meta_full
.coloc_meta$lineage <- dplyr::case_when(
  grepl('^Xen_FB_', .coloc_meta$subtype)                    ~ 'FB',
  grepl('^Pericyte|^VSMC', .coloc_meta$subtype)         ~ 'Mural',
  grepl('_EC$|^Endocardial$|^LEC$', .coloc_meta$subtype) ~ 'EC',
  grepl('^Macrophage|^Monocyte|^Dendritic|^Mast|^Myeloid_',
        .coloc_meta$subtype)                            ~ 'Myeloid',
  TRUE                                                  ~ 'Other')
.coloc_meta <- .coloc_meta[!is.na(.coloc_meta$subtype) &
                           .coloc_meta$subtype != 'Unassigned', ]

K_NN <- 15
.coloc_results <- list()
for (pat in unique(.coloc_meta$patient)) {
  d_all <- .coloc_meta[.coloc_meta$patient == pat, ]
  if (nrow(d_all) < (K_NN + 1)) next
  d_q   <- d_all[d_all$lineage == 'FB', ]
  if (nrow(d_q) == 0) next
  nn <- RANN::nn2(d_all[, c('x','y')], d_q[, c('x','y')], k = K_NN + 1)
  nn_idx <- nn$nn.idx[, -1, drop = FALSE]
  nbr <- matrix(d_all$lineage[nn_idx], ncol = K_NN, byrow = FALSE)
  d_q$frac_Myeloid <- rowMeans(nbr == 'Myeloid')
  .coloc_results[[pat]] <- d_q[, c('patient','group','subtype','frac_Myeloid')]
}
.fb_coloc <- do.call(rbind, .coloc_results)

# Per-patient background expectation (overall myeloid fraction)
.fb_bg <- .coloc_meta %>%
  dplyr::group_by(patient) %>%
  dplyr::summarise(bg_Myeloid = mean(lineage == 'Myeloid'), .groups = 'drop')

# Per-patient × subtype enrichment (mean observed / patient expected)
.fb_pp <- .fb_coloc %>%
  dplyr::group_by(patient, subtype) %>%
  dplyr::summarise(mean_frac = mean(frac_Myeloid, na.rm = TRUE),
                   .groups = 'drop') %>%
  dplyr::left_join(.fb_bg, by = 'patient') %>%
  dplyr::mutate(enrich = mean_frac / bg_Myeloid)

# Order subtypes by mean enrichment (high → low)
.fb_subtype_order <- .fb_pp %>%
  dplyr::group_by(subtype) %>%
  dplyr::summarise(m = mean(enrich, na.rm = TRUE), .groups = 'drop') %>%
  dplyr::arrange(dplyr::desc(m)) %>%
  dplyr::pull(subtype)
.fb_pp$subtype <- factor(.fb_pp$subtype, levels = rev(.fb_subtype_order))

# Per-subtype Wilcoxon vs 1.0
.fb_wilcox_v1 <- .fb_pp %>%
  dplyr::group_by(subtype) %>%
  dplyr::summarise(median_enrich = round(median(enrich, na.rm = TRUE), 2),
                   wilcox_p = tryCatch(
                     suppressWarnings(wilcox.test(enrich, mu = 1)$p.value),
                     error = function(e) NA_real_),
                   .groups = 'drop') %>%
  dplyr::mutate(label = dplyr::case_when(
                  is.na(wilcox_p)  ~ '',
                  wilcox_p < 0.001 ~ '***',
                  wilcox_p < 0.01  ~ '**',
                  wilcox_p < 0.05  ~ '*',
                  TRUE             ~ 'ns'))
cat('\n=== FB–Myeloid colocalisation: per-subtype Wilcox vs random (1.0) ===\n')
print(.fb_wilcox_v1)
write.csv(.fb_pp,        './output/Supplementary_Figure_3/fig_s3_fb_myeloid_coloc_per_patient.csv', row.names = FALSE)
write.csv(.fb_wilcox_v1, './output/Supplementary_Figure_3/fig_s3_fb_myeloid_coloc_wilcox.csv',      row.names = FALSE)

## Rotated 90°: subtypes on X, enrichment on Y; clipped 0–4
.lab_y <- 3.85   # asterisks just below upper clip limit (4)
.fb_wilcox_v1$subtype <- factor(.fb_wilcox_v1$subtype,
                                levels = rev(.fb_subtype_order))

p_fb_myel_coloc <- ggplot(.fb_pp, aes(x = subtype, y = enrich)) +
  geom_hline(yintercept = 1, color = 'grey60',
             linetype = 'dashed', linewidth = PS$geom_lw) +
  geom_boxplot(outlier.shape = NA, linewidth = PS$geom_lw,
               width = 0.55, fill = '#984EA3', alpha = 0.55) +
  geom_jitter(width = 0.12, size = 1.2, shape = 21,
              color = 'black', stroke = 0.3, fill = '#984EA3') +
  geom_text(data = .fb_wilcox_v1,
            aes(x = subtype, y = .lab_y, label = label),
            inherit.aes = FALSE, size = PS$text_mm,
            family = FONT_FAMILY, vjust = 0.5) +
  coord_cartesian(ylim = c(0, 4), clip = 'off') +
  scale_y_continuous(breaks = 0:4, expand = expansion(mult = c(0.02, 0.02))) +
  theme_v52(COMP_W) +
  theme(plot.title  = element_text(face = 'bold'),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab(NULL) +
  ylab('Myeloid-neighbour enrichment\n(observed / patient-matched expected)') +
  labs(title = 'FB–Myeloid spatial colocalization\n(K=15 NN, n=9 patients)')

# (b) Spatial vignette: pick the patient with most FB_Stress cells
.fb_stress_per_pt <- .coloc_meta %>%
  dplyr::filter(subtype == 'Xen_FB_Stressed') %>%
  dplyr::count(patient, group, sort = TRUE)
cat('\n=== FB_Stress cells per patient ===\n')
print(.fb_stress_per_pt)
.vignette_pat <- .fb_stress_per_pt$patient[1]
cat(sprintf('\nVignette tissue: patient %s (group %s)\n',
            .vignette_pat,
            as.character(.fb_stress_per_pt$group[1])))

.vmeta <- .coloc_meta[.coloc_meta$patient == .vignette_pat, ]
# Build display layer (3 colours): all baseline cells (incl. other FB) grey,
# the inflammatory-myeloid trio (Mac_Inflammatory + Monocyte + DC) in blue,
# FB_Stress in red.
.vmeta$display <- 'Other'
.vmeta$display[.vmeta$subtype %in% c('Macrophage_Inflammatory',
                                     'Monocyte')]         <- 'Inflammatory myeloid'
.vmeta$display[.vmeta$subtype == 'Xen_FB_Stressed']             <- 'Xen_FB_Stressed'
.vmeta$display <- factor(.vmeta$display,
                         levels = c('Other','Inflammatory myeloid','Xen_FB_Stressed'))
.vmeta <- .vmeta[order(.vmeta$display), ]   # plot Other first, FB_Stress last

.vig_pal <- c('Other'                = 'grey90',
              'Inflammatory myeloid' = '#225EA8',  # blue
              'Xen_FB_Stressed'            = '#E41A1C')  # red
.vig_size <- c('Other' = 0.06,
               'Inflammatory myeloid' = 0.6,
               'Xen_FB_Stressed' = 0.6)

## Vignette rotated 90° (swap aes: original x→y, original -y→x)
p_vignette <- ggplot(.vmeta, aes(x = -y, y = x,
                                 color = display, size = display)) +
  ggrastr::rasterise(geom_point(stroke = 0), dpi = 350) +
  scale_color_manual(values = .vig_pal, name = NULL) +
  scale_size_manual(values  = .vig_size, guide = 'none') +
  coord_fixed() +
  theme_v52(COMP_W) +
  theme(plot.title      = element_text(face = 'bold'),
        axis.text       = element_blank(),
        axis.ticks      = element_blank(),
        axis.line       = element_blank(),
        axis.title      = element_blank(),
        panel.grid      = element_blank(),
        legend.position = 'right',
        legend.key.size = unit(2.5, 'mm')) +
  guides(color = guide_legend(override.aes = list(size = 2.2))) +
  labs(title = sprintf('Xen_FB_Stressed + inflammatory myeloid\n(patient %s, %s)',
                       .vignette_pat,
                       as.character(.fb_stress_per_pt$group[1])))

save_figure(p_fb_myel_coloc,
            'SupplementaryFigure_3_preview_fb_myel_coloc_bar.pdf',
            width  = COMP_W * 1.0 / 2.6 * 1.25 * 0.50,
            height = COMP_H / 5 * 0.67 * 1.25)

## ── Per-patient density correlation: Xen_FB_Matrix vs total myeloid ──────
## Computed inline so all panel data flow through the script. Pearson r =
## 0.93 across n=9 patients. Output: PDF + correlations CSV.
.pat_density <- .coloc_meta %>%
  dplyr::group_by(patient, group) %>%
  dplyr::summarise(total_cells = dplyr::n(),
                   matrix_n    = sum(subtype == 'Xen_FB_Matrix', na.rm = TRUE),
                   myel_n      = sum(lineage == 'Myeloid',       na.rm = TRUE),
                   .groups     = 'drop') %>%
  dplyr::mutate(matrix_frac = matrix_n / total_cells,
                myel_frac   = myel_n   / total_cells)

.cor_p <- cor.test(.pat_density$matrix_frac, .pat_density$myel_frac,
                   method = 'pearson')
.cor_s <- cor.test(.pat_density$matrix_frac, .pat_density$myel_frac,
                   method = 'spearman', exact = FALSE)
cat(sprintf('\n=== Per-patient Xen_FB_Matrix vs myeloid density: r = %.2f (p = %s); rho = %.2f (p = %s) ===\n',
            .cor_p$estimate, signif(.cor_p$p.value, 2),
            .cor_s$estimate, signif(.cor_s$p.value, 2)))
write.csv(.pat_density,
          './output/Supplementary_Figure_3/fig_s3_xen_fb_matrix_vs_myeloid_per_patient.csv',
          row.names = FALSE)

p_matrix_myel_cor <- ggplot(.pat_density,
                            aes(x = matrix_frac * 100, y = myel_frac * 100)) +
  geom_smooth(method = 'lm', se = TRUE, color = 'grey40',
              fill = 'grey80', linewidth = PS$geom_lw, alpha = 0.4) +
  geom_point(aes(fill = group), size = 3.2, shape = 21,
             color = 'black', stroke = 0.4) +
  scale_fill_manual(values = disease_pal, name = NULL) +
  annotate('text',
           x = max(.pat_density$matrix_frac * 100, na.rm = TRUE) * 0.05,
           y = max(.pat_density$myel_frac   * 100, na.rm = TRUE) * 0.95,
           label = sprintf('Pearson r = %.2f\np = %s\nn = %d patients',
                           .cor_p$estimate, signif(.cor_p$p.value, 2),
                           nrow(.pat_density)),
           hjust = 0, vjust = 1, family = FONT_FAMILY,
           size = PS$text_mm, fontface = 'italic') +
  theme_v52(COMP_W) +
  theme(plot.title      = element_text(face = 'bold'),
        legend.position = c(0.85, 0.18),
        legend.background = element_blank(),
        legend.key.size  = unit(2.5, 'mm')) +
  xlab('Xen_FB_Matrix (% of cells)') +
  ylab('Myeloid (% of cells)') +
  ## ── Panel K — Per-patient co-abundance scatter (Xen_FB_Matrix vs myeloid density)
labs(title = 'Per-patient co-abundance: Xen_FB_Matrix ↔ myeloid')

save_figure(p_matrix_myel_cor,
            'SupplementaryFigure_3_preview_xen_fb_matrix_myel_cor.pdf',
            width  = COMP_W * 1.0 / 2.6 * 1.25 * 0.67,
            height = COMP_H / 5)
save_figure(p_vignette,
            'SupplementaryFigure_3_preview_fb_myel_coloc_vignette.pdf',
            width = COMP_W * 1.0 / 2.6, height = COMP_H / 5)

## ── Panel L — Spatial vignettes in patient 1632 (RVF)
##    (FB Stressed / FB Matrix highlighted vs iMac + Monocyte)
## Companion vignette: same tissue (patient 1632), Fb_Elastogenic in red,
## inflammatory myeloid (Mac_Inflammatory + Monocyte) in blue
.vmeta_el <- .coloc_meta[.coloc_meta$patient == .vignette_pat, ]
.vmeta_el$display <- 'Other'
.vmeta_el$display[.vmeta_el$subtype %in% c('Macrophage_Inflammatory',
                                            'Monocyte')]   <- 'Inflammatory myeloid'
.vmeta_el$display[.vmeta_el$subtype == 'Xen_FB_Matrix']   <- 'Xen_FB_Matrix'
.vmeta_el$display <- factor(.vmeta_el$display,
                            levels = c('Other','Inflammatory myeloid','Xen_FB_Matrix'))
.vmeta_el <- .vmeta_el[order(.vmeta_el$display), ]

.vig_pal_el <- c('Other'                = 'grey90',
                 'Inflammatory myeloid' = '#225EA8',
                 'Xen_FB_Matrix'       = '#E41A1C')
.vig_size_el <- c('Other' = 0.06,
                  'Inflammatory myeloid' = 0.6,
                  'Xen_FB_Matrix' = 0.6)

## Vignette 2 also rotated 90°
p_vignette_el <- ggplot(.vmeta_el, aes(x = -y, y = x,
                                       color = display, size = display)) +
  ggrastr::rasterise(geom_point(stroke = 0), dpi = 350) +
  scale_color_manual(values = .vig_pal_el, name = NULL) +
  scale_size_manual(values  = .vig_size_el, guide = 'none') +
  coord_fixed() +
  theme_v52(COMP_W) +
  theme(plot.title      = element_text(face = 'bold'),
        axis.text       = element_blank(),
        axis.ticks      = element_blank(),
        axis.line       = element_blank(),
        axis.title      = element_blank(),
        panel.grid      = element_blank(),
        legend.position = 'right',
        legend.key.size = unit(2.5, 'mm')) +
  guides(color = guide_legend(override.aes = list(size = 2.2))) +
  labs(title = sprintf('Xen_FB_Matrix + inflammatory myeloid\n(patient %s, %s)',
                       .vignette_pat,
                       as.character(.fb_stress_per_pt$group[1])))

save_figure(p_vignette_el,
            'SupplementaryFigure_3_preview_fb_homeostatic_coloc_vignette.pdf',
            width = COMP_W * 1.0 / 2.6, height = COMP_H / 5)

## Three separate panels for Illustrator assembly (no combined patchwork):
##   - Box/whisker:                _coloc_bar.pdf      (sized above)
##   - Xen_FB_Stressed vignette:   _coloc_vignette.pdf (above)
##   - Xen_FB_Matrix vignette:    _homeostatic_coloc_vignette.pdf (above)

## ───────────────────────────────────────────────────────────────────────────
## Per-myeloid-subtype FB×myeloid colocalization heatmap
##   Rows: 7 FB subtypes (ordered by overall myeloid enrichment)
##   Cols: 8 myeloid subtypes
##   Fill: log2 of median enrichment (observed / patient-matched expected)
##   Asterisks: Wilcoxon vs 1.0 across n=9 patients
## ───────────────────────────────────────────────────────────────────────────

myel_subs <- c('Macrophage_Resident','Macrophage_C1q','Macrophage_TREM2',
               'Macrophage_Inflammatory','Monocyte','Dendritic_Cell',
               'Mast_Cell','Myeloid_Proliferating')

.subtype_results <- list()
for (pat in unique(.coloc_meta$patient)) {
  d_all <- .coloc_meta[.coloc_meta$patient == pat, ]
  if (nrow(d_all) < (K_NN + 1)) next
  d_q   <- d_all[d_all$lineage == 'FB', ]
  if (nrow(d_q) == 0) next
  nn <- RANN::nn2(d_all[, c('x','y')], d_q[, c('x','y')], k = K_NN + 1)
  nn_idx <- nn$nn.idx[, -1, drop = FALSE]
  nbr <- matrix(d_all$subtype[nn_idx], ncol = K_NN, byrow = FALSE)
  for (m in myel_subs) {
    d_q[[paste0('frac_', m)]] <- rowMeans(nbr == m, na.rm = TRUE)
  }
  .subtype_results[[pat]] <- d_q
}
.fb_sub_coloc <- do.call(rbind, .subtype_results)

# Patient background per myeloid subtype
.bg_sub <- data.frame(patient = unique(.coloc_meta$patient),
                      stringsAsFactors = FALSE)
for (m in myel_subs) {
  .bg_sub[[paste0('bg_', m)]] <- sapply(.bg_sub$patient,
    function(p) mean(.coloc_meta$subtype[.coloc_meta$patient == p] == m,
                     na.rm = TRUE))
}

# Per-patient × FB-subtype × myeloid-subtype enrichment
.fb_long <- .fb_sub_coloc %>%
  dplyr::select(patient, subtype, dplyr::starts_with('frac_')) %>%
  tidyr::pivot_longer(cols = dplyr::starts_with('frac_'),
                      names_to = 'myeloid', values_to = 'frac') %>%
  dplyr::mutate(myeloid = sub('^frac_', '', myeloid))

.bg_long <- .bg_sub %>%
  tidyr::pivot_longer(cols = -patient, names_to = 'myeloid',
                      values_to = 'bg') %>%
  dplyr::mutate(myeloid = sub('^bg_', '', myeloid))

.fb_pp_sub <- .fb_long %>%
  dplyr::group_by(patient, subtype, myeloid) %>%
  dplyr::summarise(mean_frac = mean(frac, na.rm = TRUE),
                   .groups = 'drop') %>%
  dplyr::left_join(.bg_long, by = c('patient','myeloid')) %>%
  dplyr::mutate(enrich = ifelse(bg > 0 & is.finite(bg),
                                mean_frac / bg, NA_real_))

.fb_sub_agg <- .fb_pp_sub %>%
  dplyr::group_by(subtype, myeloid) %>%
  dplyr::summarise(median_enr = median(enrich, na.rm = TRUE),
                   wilcox_p   = tryCatch(suppressWarnings(
                                  wilcox.test(enrich, mu = 1)$p.value),
                                  error = function(e) NA_real_),
                   .groups = 'drop') %>%
  dplyr::mutate(log2_enr = log2(median_enr),
                sig = dplyr::case_when(
                  is.na(wilcox_p)  ~ '',
                  wilcox_p < 0.001 ~ '***',
                  wilcox_p < 0.01  ~ '**',
                  wilcox_p < 0.05  ~ '*',
                  TRUE             ~ ''))

cat('\n=== FB × myeloid subtype enrichment heatmap data ===\n')
print(.fb_sub_agg %>%
        dplyr::mutate(median_enr = round(median_enr, 2),
                      log2_enr   = round(log2_enr, 2),
                      wilcox_p   = signif(wilcox_p, 2)))
write.csv(.fb_pp_sub,
          './output/Supplementary_Figure_3/fig_s3_fb_myel_subtype_coloc_per_patient.csv',
          row.names = FALSE)
write.csv(.fb_sub_agg,
          './output/Supplementary_Figure_3/fig_s3_fb_myel_subtype_coloc_summary.csv',
          row.names = FALSE)

# Order FB rows by overall median myeloid enrichment (high → low; high at top)
.fb_order_overall <- .fb_pp %>%
  dplyr::group_by(subtype) %>%
  dplyr::summarise(m = median(enrich, na.rm = TRUE), .groups = 'drop') %>%
  dplyr::arrange(dplyr::desc(m)) %>%
  dplyr::pull(subtype)

# Order myeloid columns by total enrichment magnitude (most-enriched left)
.myel_order <- .fb_sub_agg %>%
  dplyr::group_by(myeloid) %>%
  dplyr::summarise(s = sum(pmax(log2_enr, 0), na.rm = TRUE),
                   .groups = 'drop') %>%
  dplyr::arrange(dplyr::desc(s)) %>%
  dplyr::pull(myeloid)

.fb_sub_agg$subtype <- factor(.fb_sub_agg$subtype, levels = rev(.fb_order_overall))
.fb_sub_agg$myeloid <- factor(.fb_sub_agg$myeloid, levels = .myel_order)

.heat_lim <- max(abs(.fb_sub_agg$log2_enr), na.rm = TRUE)

p_fb_myel_heat <- ggplot(.fb_sub_agg,
                         aes(x = myeloid, y = subtype, fill = log2_enr)) +
  geom_tile(color = 'grey85', linewidth = PS$linewidth_mm * 0.3) +
  geom_text(aes(label = sig), color = 'black',
            size = PS$text_mm, vjust = 0.7) +
  scale_fill_gradient2(low = '#3B4CC0', mid = 'white', high = '#B40426',
                       midpoint = 0,
                       limits = c(-.heat_lim, .heat_lim),
                       oob = scales::squish, na.value = 'grey92',
                       name = expression(log[2]~enrich)) +
  theme_v52(COMP_W) +
  theme(plot.title  = element_text(face = 'bold'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.height = unit(3.5, 'mm'),
        legend.key.width  = unit(1.5, 'mm')) +
  xlab(NULL) + ylab(NULL) +
  labs(title = 'FB × Myeloid spatial colocalization\n(K=15 NN, median enrichment, n=9)')

save_figure(p_fb_myel_heat,
            'SupplementaryFigure_3_preview_fb_myel_subtype_heat.pdf',
            width = COMP_W * 1.4 / 2.6, height = COMP_H / 4)

## ───────────────────────────────────────────────────────────────────────────
## Panel E: FB Phase 2 (pRV→RVF) — sn vs Xenium concordance scatter
##           (script-letter D → v57-letter E)
##   Lineage-level signal is quiet on both platforms (manuscript-confirmed),
##   but direction is concordant for the canonical fibrogenic genes that
##   are present in both panels (TNC, THBS4, ELN, FN1, VCAN, PRG4, BAMBI,
##   FBN1). Sublineage Xenium hits (FB_Activated/FB_NOX4/FB_NTN4) reported
##   in Data S9 — annotation here flags THBS4 / TNC / PRG4 as the
##   commitment-defining genes per Fu 2018, Amrute 2024.
## ───────────────────────────────────────────────────────────────────────────

fb_p2_panel_genes <- c('TNC','THBS4','ELN','FN1','VCAN','PRG4','BAMBI','FBN1')

.fb_p2_xy <- function(genes, contrast_label) {
  s <- sn_lin %>%
    dplyr::filter(subtype == 'FB', contrast == contrast_label,
                  gene %in% genes) %>%
    dplyr::transmute(gene, sn_lfc  = log2FoldChange,
                     sn_padj  = padj, sn_bm = baseMean)
  x <- xen_lin %>%
    dplyr::filter(subtype == 'FB', contrast == contrast_label,
                  gene %in% genes) %>%
    dplyr::transmute(gene, xen_lfc = log2FoldChange,
                     xen_padj = padj, xen_bm = baseMean)

  # Best (smallest) Xenium SUBLINEAGE padj among the 7 FB Xenium subtypes.
  # baseMean filter dropped at sublineage level — sublineage clusters have
  # fewer cells so a lineage-level bm>50 cutoff is too strict (e.g., it would
  # reject FB_NTN4 hits ELN bm=37 / PRG4 bm=22 despite padj<0.02 and |lfc|>2).
  fb_xen_subs <- c('Xen_FB_Activated','Xen_FB_NOX4','Xen_FB_Interstitial',
                   'Xen_FB_Intermediate','Xen_FB_CILP','Xen_FB_Stressed',
                   'Xen_FB_Matrix')
  xs <- xen_sub_csv %>%
    dplyr::filter(subtype %in% fb_xen_subs, contrast == contrast_label,
                  gene %in% genes,
                  !is.na(padj), padj < sig_padj,
                  abs(log2FoldChange) > sig_lfc) %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(xen_sub_best_padj   = min(padj, na.rm = TRUE),
                     xen_sub_best_lfc    = log2FoldChange[which.min(padj)],
                     xen_sub_best_subset = subtype[which.min(padj)],
                     xen_sub_n_sig       = dplyr::n(),
                     .groups = 'drop')

  d <- merge(s, x, by = 'gene') %>%
    dplyr::left_join(xs, by = 'gene')
  d$sig_in <- with(d, dplyr::case_when(
    !is.na(sn_padj) & sn_padj < sig_padj &
      !is.na(xen_sub_best_padj)               ~ 'both',
    !is.na(sn_padj) & sn_padj < sig_padj      ~ 'snRNA only',
    !is.na(xen_sub_best_padj)                 ~ 'Xenium sublineage',
    TRUE                                      ~ 'neither'))
  d$min_padj <- pmin(d$sn_padj, d$xen_sub_best_padj, na.rm = TRUE)
  d$min_padj[is.na(d$min_padj)] <- 1
  d
}

fb_p2_d <- .fb_p2_xy(fb_p2_panel_genes, 'RVF_vs_pRV')
.lim    <- max(abs(c(fb_p2_d$sn_lfc, fb_p2_d$xen_lfc)), na.rm = TRUE) * 1.15

.sig_pal <- c('both'              = '#7B3294',
              'snRNA only'        = '#1B7837',
              'Xenium sublineage' = '#C2A5CF',
              'neither'           = 'grey55')

p_D <- ggplot(fb_p2_d, aes(x = sn_lfc, y = xen_lfc, color = sig_in)) +
  geom_hline(yintercept = 0, color = 'grey80',
             linewidth = PS$geom_lw) +
  geom_vline(xintercept = 0, color = 'grey80',
             linewidth = PS$geom_lw) +
  geom_abline(slope = 1, intercept = 0, color = 'grey80',
              linetype = 'dashed', linewidth = PS$geom_lw) +
  geom_point(aes(size = -log10(min_padj)), alpha = 0.85) +
  geom_text_repel(aes(label = gene), color = 'black',
                  size = PS$text_mm, family = FONT_FAMILY,
                  fontface = 'italic',
                  max.overlaps = Inf,
                  segment.size = PS$geom_lw, min.segment.length = 0) +
  scale_color_manual(values = .sig_pal, name = 'Sig') +
  scale_size_continuous(range = c(1.2, 4), name = expression(-log[10]~p[adj])) +
  coord_equal(xlim = c(-.lim, .lim), ylim = c(-.lim, .lim)) +
  theme_v52(COMP_W) +
  theme(plot.title        = element_text(face = 'bold'),
        legend.position   = 'right',
        legend.key.size   = unit(2.5, 'mm'),
        legend.background = element_blank()) +
  xlab(expression('snRNA-seq '*log[2]~FC)) +
  ylab(expression('Xenium '*log[2]~FC)) +
  labs(title = 'FB Phase 2 (pRV→RVF)\nsn–Xenium concordance')

## ───────────────────────────────────────────────────────────────────────────
## Companion panel: Xenium FB SUBLINEAGE Phase 2 heatmap
##   Resolves the focal fibrogenic commitment that the lineage-level
##   scatter shows as dilute. Xenium uses its own FB subtype scheme
##   (FB_Activated, FB_NOX4, FB_NTN4, etc.) — distinct from snRNA
##   Fb1..Fb7. Asterisks = padj<0.05, |log2FC|>0.5, baseMean>50.
## ───────────────────────────────────────────────────────────────────────────

## Xenium FB subtype labels (Xen_ prefix marks Xenium origin; distinct from snRNA Fb_*)
fb_xen_subs_order <- c('Xen_FB_Activated','Xen_FB_NOX4','Xen_FB_Interstitial',
                       'Xen_FB_Intermediate','Xen_FB_CILP','Xen_FB_Stressed',
                       'Xen_FB_Matrix')
fb_p2_sublineage_genes <- c('TNC','THBS4','PRG4','ELN','FN1','VCAN','BAMBI','FBN1')

xen_p2_sub <- xen_sub_csv %>%
  dplyr::filter(subtype %in% fb_xen_subs_order,
                contrast == 'RVF_vs_pRV',
                gene %in% fb_p2_sublineage_genes) %>%
  dplyr::transmute(gene, subtype, log2FoldChange, padj, baseMean)

.grid <- expand.grid(gene    = fb_p2_sublineage_genes,
                     subtype = fb_xen_subs_order,
                     stringsAsFactors = FALSE)
xen_p2_sub <- dplyr::left_join(.grid, xen_p2_sub,
                               by = c('gene','subtype')) %>%
  dplyr::mutate(sig = ifelse(!is.na(padj) & padj < sig_padj &
                             abs(log2FoldChange) > sig_lfc,
                             '*', ''),   # bm filter dropped at sublineage level
                gene    = factor(gene,    levels = rev(fb_p2_sublineage_genes)),
                subtype = factor(subtype, levels = fb_xen_subs_order))

p_D2 <- ggplot(xen_p2_sub,
               aes(x = subtype, y = gene, fill = log2FoldChange)) +
  geom_tile(color = 'grey85', linewidth = PS$linewidth_mm * 0.3) +
  geom_text(aes(label = sig), color = 'black',
            size = PS$text_mm, vjust = 0.75) +
  scale_fill_gradient2(low = '#3B4CC0', mid = 'white', high = '#B40426',
                       midpoint = 0, limits = c(-4.5, 4.5),
                       oob = scales::squish, na.value = 'grey92',
                       name = expression(log[2]~FC)) +
  theme_v52(COMP_W) +
  theme(plot.title  = element_text(face = 'bold'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(face = 'italic'),
        legend.key.height = unit(3.5, 'mm'),
        legend.key.width  = unit(1.5, 'mm')) +
  xlab(NULL) + ylab(NULL) +
  labs(title = 'Xenium FB subtype Phase 2\n(focal fibrogenic commitment)')

save_figure(p_D2,
            'SupplementaryFigure_3_preview_fb_p2_sublineage.pdf',
            width = COMP_W * 1.4 / 2.9, height = COMP_H / 4)

## ───────────────────────────────────────────────────────────────────────────
## Supplementary (legacy v54 Panel E): FB ChEA enrichment (Phase 1 down)
##           Not in v57 main caption — kept for reference
## ───────────────────────────────────────────────────────────────────────────

fb_chea_cache <- './output/Supplementary_Figure_3/fig_s3_fb_chea.rds'
if (file.exists(fb_chea_cache)) {
  fb_chea <- readRDS(fb_chea_cache)
} else {
  fb_down <- sn_sub %>%
    dplyr::filter(subtype %in% fb_subs, contrast == 'pRV_vs_NF',
                  !is.na(padj), padj < sig_padj,
                  log2FoldChange < -sig_lfc, baseMean > sig_bmean) %>%
    dplyr::pull(gene) %>% unique()
  fb_chea <- enrichR::enrichr(fb_down, c('ChEA_2016'))$ChEA_2016
  dir.create('./output', showWarnings = FALSE, recursive = TRUE)
  saveRDS(fb_chea, fb_chea_cache)
}

fb_chea_tf <- fb_chea %>%
  dplyr::mutate(TF = str_extract(Term, '^[^_ ]+')) %>%
  dplyr::filter(!is.na(TF)) %>%
  dplyr::group_by(TF) %>%
  dplyr::slice_min(Adjusted.P.value, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(Adjusted.P.value) %>%
  dplyr::slice_head(n = 8)

p_E <- .make_chea_bar(fb_chea_tf, '#3e64b3', 'FB ChEA (Ph1 down)')

## ───────────────────────────────────────────────────────────────────────────
## Panel N: Mural cell subclustering UMAP (PC + SM)
##           (script-letter F → v57-letter N)
## ───────────────────────────────────────────────────────────────────────────

## Use the original subclustered 2,152-cell list (Pc=1141, Sm=1011) defined
## in mural_subclustering.csv. RECOMPUTE UMAP on the subset so the
## visualization matches the original (the global snRV_ref UMAP coordinates
## leave the mural cells as small clusters in the global embedding).
## Result is cached so subsequent runs reuse it.
.mural_cache <- './output/Supplementary_Figure_3/mural_subclust_recomputed.rds'
if (file.exists(.mural_cache)) {
  M_mural <- readRDS(.mural_cache)
  cat(sprintf('Loaded cached mural object: n = %d cells\n', ncol(M_mural)))
} else {
  .mural_csv <- read.csv('./dependencies/shared/mural_subclustering.csv',
                         stringsAsFactors = FALSE)
  ref   <- readRDS('./dependencies/shared/snRV_ref.rds')
  .keep <- intersect(.mural_csv$colnames.M1., colnames(ref))
  M_mural <- subset(ref, cells = .keep)
  .id_map <- setNames(.mural_csv$M1.active.ident, .mural_csv$colnames.M1.)
  M_mural$mural_subtype <- factor(unname(.id_map[colnames(M_mural)]),
                                  levels = c('Pc','Sm'))
  rm(ref); gc(verbose = FALSE)
  cat(sprintf('Recomputing UMAP on %d-cell mural subset...\n', ncol(M_mural)))
  DefaultAssay(M_mural) <- 'RNA'
  M_mural <- NormalizeData(M_mural, verbose = FALSE)
  M_mural <- FindVariableFeatures(M_mural, nfeatures = 2000, verbose = FALSE)
  M_mural <- ScaleData(M_mural, verbose = FALSE)
  M_mural <- RunPCA(M_mural, npcs = 30, verbose = FALSE)
  M_mural <- RunUMAP(M_mural, dims = 1:20, verbose = FALSE)
  saveRDS(M_mural, .mural_cache)
  cat('Cached recomputed mural Seurat object to', .mural_cache, '\n')
}
mural_group_col <- 'mural_subtype'
Idents(M_mural) <- mural_group_col
cat(sprintf('\nMural UMAP cohort: n = %d nuclei (Pc=%d, Sm=%d)\n',
            ncol(M_mural),
            sum(M_mural$mural_subtype=='Pc'),
            sum(M_mural$mural_subtype=='Sm')))

p_F <- DimPlot(M_mural, group.by = mural_group_col, label = TRUE,
               label.size = PS$base_pt / .pt, repel = TRUE,
               pt.size = PS$umap_pt, raster = TRUE) +
  NoLegend() +
  labs(title = 'Mural subtypes') +
  theme_v52(COMP_W) +
  theme(plot.title = element_text(face = 'bold'),
        axis.text  = element_blank(),
        axis.ticks = element_blank(),
        axis.line  = element_blank())

## ───────────────────────────────────────────────────────────────────────────
## Panel O: Mural Phase 1 (pRV vs NF) snRNA sublineage heatmap
##           (script-letter G → v57-letter O)
## ───────────────────────────────────────────────────────────────────────────

## ── Panel P — Mural Phase 1 vascular-protective pseudobulk module score (PC + SM pooled)
## Mural Phase 1 — VASCULAR-PROTECTIVE program loss (PC identity preserved).
## Pseudobulk module-score analysis (panels below) shows canonical PC
## identity markers (KCNJ8/RGS5/NOTCH3) are NOT significantly lost,
## consistent with cardiac lineage tracing arguing against PMT (Greenhalgh
## 2015 J Pathol; Kanisicak 2016 Nat Commun) — pericytes retain identity
## but selectively shed their vascular-protective program. Heatmap shows
## the genes that ARE lost (Figure 2 callouts + cross-platform validated):
##   metabolic loss      : NNMT, ADH1B, CYP4B1
##   homeostasis         : SCGN, C7, ELN, SRPX
##   anti-protease       : TIMP3
##   decoy receptor      : IL1RL1
##   anti-inflammatory   : PLA2G2A
##   antioxidant         : GPX3
mural_ph1_genes <- c('NNMT','ADH1B','CYP4B1',
                     'SCGN','C7','ELN','SRPX',
                     'TIMP3','IL1RL1','PLA2G2A','GPX3')
mural_ph1       <- .prep_phase_heatmap(sn_sub, mural_subs,
                                       'pRV_vs_NF', mural_ph1_genes)
p_G <- .make_phase_heatmap(mural_ph1, 'Mural Phase 1 (pRV vs NF)\nvascular-protective loss (PC identity preserved)')

## Parallel Xenium sublineage Phase 1 heatmap: same genes, Xenium PC/VSMC
## subsets (8 mural subtypes), reveals whether loss is subset-restricted
## Drop generic 'Pericyte' and 'VSMC' parents — keep specific subsets only
xen_mural_subs <- c('Pericyte_KCNJ8','Pericyte_Mixed','Pericyte_THBS4',
                    'VSMC_Arterial','VSMC_Inflamed','VSMC_Synthetic')
## Restrict Xenium Phase 1 heatmap to on-panel genes only (drop NNMT, ADH1B)
mural_ph1_genes_xen <- c('CYP4B1','SCGN','C7','ELN','SRPX',
                         'TIMP3','IL1RL1','PLA2G2A','GPX3')
xen_mural_ph1 <- xen_sub_csv %>%
  dplyr::filter(subtype %in% xen_mural_subs, contrast == 'pRV_vs_NF',
                gene %in% mural_ph1_genes_xen) %>%
  dplyr::transmute(gene, subtype, log2FoldChange, padj, baseMean)
.grid_xm <- expand.grid(gene = mural_ph1_genes_xen,
                        subtype = xen_mural_subs,
                        stringsAsFactors = FALSE)
xen_mural_ph1 <- dplyr::left_join(.grid_xm, xen_mural_ph1,
                                  by = c('gene','subtype')) %>%
  dplyr::mutate(sig = ifelse(!is.na(padj) & padj < sig_padj &
                             abs(log2FoldChange) > sig_lfc,
                             '*', ''),
                gene    = factor(gene,    levels = rev(mural_ph1_genes_xen)),
                subtype = factor(subtype, levels = xen_mural_subs))

p_G_xen <- ggplot(xen_mural_ph1,
                  aes(x = subtype, y = gene, fill = log2FoldChange)) +
  geom_tile(color = 'grey85', linewidth = PS$linewidth_mm * 0.3) +
  geom_text(aes(label = sig), color = 'black',
            size = PS$text_mm, vjust = 0.75) +
  scale_fill_gradient2(low = '#3B4CC0', mid = 'white', high = '#B40426',
                       midpoint = 0, limits = c(-4, 4),
                       oob = scales::squish, na.value = 'grey92',
                       name = expression(log[2]~FC)) +
  theme_v52(COMP_W) +
  theme(plot.title  = element_text(face = 'bold'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(face = 'italic'),
        legend.key.height = unit(3.5, 'mm'),
        legend.key.width  = unit(1.5, 'mm')) +
  xlab(NULL) + ylab(NULL) +
  labs(title = 'Mural Phase 1 (Xenium sublineage)\n8 PC/VSMC subsets')

save_figure(p_G_xen,
            'SupplementaryFigure_3_preview_mural_p1_xen_sublineage.pdf',
            width = COMP_W * 1.4 / 2.9, height = COMP_H / 4)

## Parallel Xenium sublineage Phase 2 mural heatmap — restricted to genes
## on the actual 477-gene Xenium panel (no SpaGE imputation noise)
mural_ph2_genes_xen <- c('NR4A1','NR4A3','ACTA2','MYH11','CNN1','MYLK')
xen_mural_ph2 <- xen_sub_csv %>%
  dplyr::filter(subtype %in% xen_mural_subs, contrast == 'RVF_vs_pRV',
                gene %in% mural_ph2_genes_xen) %>%
  dplyr::transmute(gene, subtype, log2FoldChange, padj, baseMean)
.grid_xm2 <- expand.grid(gene = mural_ph2_genes_xen,
                         subtype = xen_mural_subs,
                         stringsAsFactors = FALSE)
xen_mural_ph2 <- dplyr::left_join(.grid_xm2, xen_mural_ph2,
                                  by = c('gene','subtype')) %>%
  dplyr::mutate(sig = ifelse(!is.na(padj) & padj < sig_padj &
                             abs(log2FoldChange) > sig_lfc,
                             '*', ''),
                gene    = factor(gene,    levels = rev(mural_ph2_genes_xen)),
                subtype = factor(subtype, levels = xen_mural_subs))

p_H_xen <- ggplot(xen_mural_ph2,
                  aes(x = subtype, y = gene, fill = log2FoldChange)) +
  geom_tile(color = 'grey85', linewidth = PS$linewidth_mm * 0.3) +
  geom_text(aes(label = sig), color = 'black',
            size = PS$text_mm, vjust = 0.75) +
  scale_fill_gradient2(low = '#3B4CC0', mid = 'white', high = '#B40426',
                       midpoint = 0, limits = c(-4, 4),
                       oob = scales::squish, na.value = 'grey92',
                       name = expression(log[2]~FC)) +
  theme_v52(COMP_W) +
  theme(plot.title  = element_text(face = 'bold'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(face = 'italic'),
        legend.key.height = unit(3.5, 'mm'),
        legend.key.width  = unit(1.5, 'mm')) +
  xlab(NULL) + ylab(NULL) +
  labs(title = 'Mural Phase 2 (Xenium sublineage)\n8 PC/VSMC subsets')

save_figure(p_H_xen,
            'SupplementaryFigure_3_preview_mural_p2_xen_sublineage.pdf',
            width = COMP_W * 1.4 / 2.9, height = COMP_H / 4)

## ── Panel S — Mural Phase 2 sn-vs-Xenium concordance scatter (PC + SM pooled)
## ── Mural Phase 2 sn-vs-Xenium concordance scatter (lineage level) ────────
mural_p2_scatter_genes <- c('NR4A1','NR4A2','NR4A3','FOS','JUNB','JUN','HSPB1',
                            'ATF3','ACTA2','MYH11','CNN1','MYOCD','MYLK')

.mural_p2_xy <- function(genes, contrast_label = 'RVF_vs_pRV') {
  .agg_pcsm <- function(df) df %>%
    dplyr::filter(subtype %in% c('PC','SM'), contrast == contrast_label,
                  gene %in% genes) %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(lfc  = mean(log2FoldChange, na.rm = TRUE),
                     padj = min(padj,            na.rm = TRUE),
                     bm   = mean(baseMean,       na.rm = TRUE),
                     .groups = 'drop')
  s <- .agg_pcsm(sn_lin)  %>% dplyr::transmute(gene, sn_lfc=lfc,  sn_padj=padj,  sn_bm=bm)
  x <- .agg_pcsm(xen_lin) %>% dplyr::transmute(gene, xen_lfc=lfc, xen_padj=padj, xen_bm=bm)
  d <- merge(s, x, by = 'gene')
  d$sig_in <- with(d, dplyr::case_when(
    !is.na(sn_padj)  & sn_padj  < sig_padj &
      !is.na(xen_padj) & xen_padj < sig_padj  ~ 'both',
    !is.na(sn_padj)  & sn_padj  < sig_padj    ~ 'snRNA only',
    !is.na(xen_padj) & xen_padj < sig_padj    ~ 'Xenium only',
    TRUE                                      ~ 'neither'))
  d$min_padj <- pmin(d$sn_padj, d$xen_padj, na.rm = TRUE)
  d$min_padj[is.na(d$min_padj)] <- 1
  d
}

mural_p2_d <- .mural_p2_xy(mural_p2_scatter_genes)
.mlim     <- max(abs(c(mural_p2_d$sn_lfc, mural_p2_d$xen_lfc)), na.rm = TRUE) * 1.15
.mural_pal <- c('both' = '#7B3294', 'snRNA only' = '#1B7837',
                'Xenium only' = '#C2A5CF', 'neither' = 'grey55')

p_mural_p2_scatter <- ggplot(mural_p2_d,
                             aes(x = sn_lfc, y = xen_lfc, color = sig_in)) +
  geom_hline(yintercept = 0, color = 'grey80', linewidth = PS$geom_lw) +
  geom_vline(xintercept = 0, color = 'grey80', linewidth = PS$geom_lw) +
  geom_abline(slope = 1, intercept = 0, color = 'grey80',
              linetype = 'dashed', linewidth = PS$geom_lw) +
  geom_point(aes(size = -log10(min_padj)), alpha = 0.85) +
  geom_text_repel(aes(label = gene), color = 'black',
                  size = PS$text_mm, family = FONT_FAMILY,
                  fontface = 'italic',
                  max.overlaps = Inf,
                  segment.size = PS$geom_lw, min.segment.length = 0) +
  scale_color_manual(values = .mural_pal, name = 'Sig') +
  scale_size_continuous(range = c(1.2, 4), name = expression(-log[10]~p[adj])) +
  coord_equal(xlim = c(-.mlim, .mlim), ylim = c(-.mlim, .mlim)) +
  theme_v52(COMP_W) +
  theme(plot.title = element_text(face = 'bold'),
        legend.position = 'right',
        legend.key.size = unit(2.5,'mm'),
        legend.background = element_blank()) +
  xlab(expression('snRNA-seq '*log[2]~FC*' (PC+SM)')) +
  ylab(expression('Xenium '*log[2]~FC*' (PC+SM)')) +
  labs(title = 'Mural Phase 2 (pRV→RVF)\nsn–Xenium concordance')

save_figure(p_mural_p2_scatter,
            'SupplementaryFigure_3_preview_mural_p2_concordance.pdf',
            width = COMP_W * 1.2 / 2.9, height = COMP_H / 4)

## ── Mural Phase 2 — TWO literature-anchored pseudobulk scores ─────────────
##   (a) Vascular IEG/stress — Hamers 2019 Circ Res NR4A axis +
##       canonical AP-1 (FOS/JUN/B) + heat shock (HSPB1) + ATF3
##   (b) Contractile reversion — Owens & Yoshida 2005 Physiol Rev contractile
##       VSMC set; MYH11 EXCLUDED because it goes the wrong way at Xenium
##       sublineage (loss in Pericyte subtypes, not gain)
ieg_stress_genes  <- c('NR4A1','NR4A2','NR4A3','JUN','JUNB','FOS','HSPB1','ATF3')
contractile_genes <- c('CNN1','ACTA2','TAGLN','MYOCD','MYLK')

if (file.exists('./dependencies/shared/snRV_ref.rds')) {
  ref <- readRDS('./dependencies/shared/snRV_ref.rds')
  M_mural_full <- subset(ref, subset = Subnames_manual %in% c('PC','SM'))
  rm(ref); invisible(gc(verbose = FALSE))

  .mu_agg2  <- AggregateExpression(M_mural_full, assays = 'RNA',
                                   group.by = 'patient', slot = 'counts',
                                   return.seurat = FALSE)$RNA
  .mu_lcpm2 <- log1p(sweep(.mu_agg2, 2, colSums(.mu_agg2), '/') * 1e6)
  .mu_pat_grp2 <- M_mural_full@meta.data %>%
    dplyr::distinct(patient, group) %>%
    dplyr::mutate(patient = as.character(patient),
                  group   = factor(group, levels = c('NF','pRV','RVF')))

  .build_mural_score2 <- function(gene_set, label) {
    have <- intersect(gene_set, rownames(.mu_lcpm2))
    cat(sprintf('%s — detected: %s\n', label, paste(have, collapse = ', ')))
    z <- t(scale(t(.mu_lcpm2[have, , drop = FALSE])))
    d <- data.frame(patient = sub('^g','', colnames(z)),
                    score   = colMeans(z, na.rm = TRUE),
                    stringsAsFactors = FALSE)
    dplyr::left_join(d, .mu_pat_grp2, by = 'patient') %>%
      dplyr::arrange(group, patient)
  }
  cat('\n=== Mural Phase 2 — two literature-anchored scores ===\n')
  ## IEG/stress score — PC-ONLY (signal is pericyte-driven; SM is non-
  ## monotonic and dilutes; biologically more interpretable as PC stress)
  M_pc_only <- subset(M_mural_full, subset = Subnames_manual == 'PC')
  .pc_agg   <- AggregateExpression(M_pc_only, assays='RNA', group.by='patient',
                                    slot='counts', return.seurat=FALSE)$RNA
  .pc_lcpm  <- log1p(sweep(.pc_agg, 2, colSums(.pc_agg), '/') * 1e6)
  .pc_have  <- intersect(ieg_stress_genes, rownames(.pc_lcpm))
  cat(sprintf('PC-only IEG/stress (Hamers 2019) — detected: %s\n',
              paste(.pc_have, collapse = ', ')))
  .pc_z <- t(scale(t(.pc_lcpm[.pc_have, , drop = FALSE])))
  mu_ieg_score <- data.frame(patient = sub('^g','', colnames(.pc_z)),
                             score   = colMeans(.pc_z, na.rm = TRUE),
                             stringsAsFactors = FALSE) %>%
    dplyr::left_join(.mu_pat_grp2, by = 'patient') %>%
    dplyr::arrange(group, patient)
  rm(M_pc_only, .pc_agg, .pc_lcpm); invisible(gc(verbose = FALSE))

  mu_cont_score <- .build_mural_score2(contractile_genes, 'Contractile (Owens 2005, MYH11 excluded; PC+SM)')

  p_mural_ieg <- .fb_module_plot(
    mu_ieg_score,
    ylab_text  = 'Vascular IEG/stress score (z, pseudobulk)',
    title_text = 'Phase 2: vascular IEG/stress (PC only)\n(NR4A1-3 + AP-1 + HSPB1, Hamers 2019)',
    csv_tag    = 'mural_ieg_stress_PC_pseudobulk')
  p_mural_cont <- .fb_module_plot(
    mu_cont_score,
    ylab_text  = 'Contractile-reversion score (z, pseudobulk)',
    title_text = 'Phase 2: partial contractile reversion\n(CNN1/ACTA2/TAGLN/MYOCD/MYLK)',
    csv_tag    = 'mural_contractile_pseudobulk')

  save_figure(p_mural_ieg,
              'SupplementaryFigure_3_preview_mural_ieg_stress_pseudobulk.pdf',
              width = COMP_W * 1.0 / 2.9, height = COMP_H / 4)
  save_figure(p_mural_cont,
              'SupplementaryFigure_3_preview_mural_contractile_pseudobulk.pdf',
              width = COMP_W * 1.0 / 2.9, height = COMP_H / 4)
  save_figure((p_mural_ieg | p_mural_cont),
              'SupplementaryFigure_3_preview_mural_phase2_scores_paired.pdf',
              width = COMP_W * 2.0 / 2.9, height = COMP_H / 4)
  rm(M_mural_full, .mu_agg2, .mu_lcpm2); invisible(gc(verbose = FALSE))
}

## ── Mural Phase 1 homeostatic-loss module score (pseudobulk per patient) ──
## Uses snRNA mural lineage (PC + SM pooled) — mirrors the FB-identity score
## (Koenig 2022) panel D structure. Score = mean z-score across the 14
## Phase 1 homeostatic genes per patient pseudobulk.
.mural_xen_path <- './dependencies/shared/Xenium/mural_clean_clean.rds'
M_mural_full <- NULL
if (file.exists('./dependencies/shared/snRV_ref.rds')) {
  ref <- readRDS('./dependencies/shared/snRV_ref.rds')
  M_mural_full <- subset(ref, subset = Subnames_manual %in% c('PC','SM'))
  rm(ref); invisible(gc(verbose = FALSE))
}
if (!is.null(M_mural_full)) {
  .mu_agg  <- AggregateExpression(M_mural_full, assays = 'RNA',
                                  group.by = 'patient', slot = 'counts',
                                  return.seurat = FALSE)$RNA
  .mu_lib  <- colSums(.mu_agg)
  .mu_lcpm <- log1p(sweep(.mu_agg, 2, .mu_lib, '/') * 1e6)
  .mu_pat_grp <- M_mural_full@meta.data %>%
    dplyr::distinct(patient, group) %>%
    dplyr::mutate(patient = as.character(patient),
                  group   = factor(group, levels = c('NF','pRV','RVF')))

  ## Two literature-anchored scores (Path 1 — published-grounded):
  ##   (a) PC identity:   KCNJ8 + RGS5 + NOTCH3 (Vanlandewijck 2018 Nature
  ##                      pericyte canonical triad)
  ##   (b) Vascular protective: TIMP3 + IL1RL1 + GPX3 + ELN (cardiac-protective
  ##       canonical; Kassiri 2005 / Yu 2017 / Kakkar 2008 / atlas)
  pc_identity_genes <- c('KCNJ8','RGS5','NOTCH3')
  vasc_protect_genes <- c('TIMP3','IL1RL1','GPX3','ELN')

  .build_mural_score <- function(gene_set, label) {
    have <- intersect(gene_set, rownames(.mu_lcpm))
    cat(sprintf('%s — detected: %s\n', label, paste(have, collapse = ', ')))
    z <- t(scale(t(.mu_lcpm[have, , drop = FALSE])))
    d <- data.frame(patient = sub('^g','', colnames(z)),
                    score   = colMeans(z, na.rm = TRUE),
                    stringsAsFactors = FALSE)
    dplyr::left_join(d, .mu_pat_grp, by = 'patient') %>%
      dplyr::arrange(group, patient)
  }

  cat('\n=== Mural Phase 1 — two literature-anchored scores (snRNA pseudobulk) ===\n')
  mu_pc_score   <- .build_mural_score(pc_identity_genes,   'PC-identity (Vanlandewijck 2018 triad)')
  mu_vasc_score <- .build_mural_score(vasc_protect_genes,  'Vascular-protective (TIMP3/IL1RL1/GPX3/ELN)')
  cat('\n--- PC-identity score per patient ---\n');     print(mu_pc_score)
  cat('\n--- Vascular-protective score per patient ---\n'); print(mu_vasc_score)

  p_mural_pc <- .fb_module_plot(
    mu_pc_score,
    ylab_text  = 'PC-identity score (z, pseudobulk)',
    title_text = 'Phase 1: PC identity loss\n(KCNJ8/RGS5/NOTCH3 — Vanlandewijck 2018)',
    csv_tag    = 'mural_pc_identity_pseudobulk')
  p_mural_vasc <- .fb_module_plot(
    mu_vasc_score,
    ylab_text  = 'Vascular-protective score (z, pseudobulk)',
    title_text = 'Phase 1: vascular protective loss\n(TIMP3/IL1RL1/GPX3/ELN)',
    csv_tag    = 'mural_vasc_protective_pseudobulk')

  save_figure(p_mural_pc,
              'SupplementaryFigure_3_preview_mural_pc_identity_pseudobulk.pdf',
              width = COMP_W * 1.0 / 2.9, height = COMP_H / 4)
  save_figure(p_mural_vasc,
              'SupplementaryFigure_3_preview_mural_vasc_protective_pseudobulk.pdf',
              width = COMP_W * 1.0 / 2.9, height = COMP_H / 4)
  save_figure((p_mural_pc | p_mural_vasc),
              'SupplementaryFigure_3_preview_mural_phase1_scores_paired.pdf',
              width = COMP_W * 2.0 / 2.9, height = COMP_H / 4)
  rm(M_mural_full, .mu_agg, .mu_lcpm); invisible(gc(verbose = FALSE))
}

## ───────────────────────────────────────────────────────────────────────────
## Panel Q: Mural Phase 2 (RVF vs pRV) snRNA sublineage heatmap
##           (script-letter H → v57-letter Q)
## ───────────────────────────────────────────────────────────────────────────

## Mural Phase 2 — expanded gene set covering two programs.
## NB: "quiescence" is partial — CNN1+MYOCD trend up but MYH11 goes the
## opposite direction at Xenium subtype level (excluded from contractile
## score below; see commentary). Better framed as IEG/stress + partial
## contractile reversion than as "quiescence" per se.
##   IEG / AP-1 stress       : NR4A1, NR4A2, NR4A3, FOS, JUNB, JUN, HSPB1, ATF3
##   Contractile (Owens 2005): ACTA2, MYH11, CNN1, MYOCD, MYLK
mural_ph2_genes <- c('NR4A1','NR4A2','NR4A3','FOS','JUNB','JUN','HSPB1','ATF3',
                     'ACTA2','MYH11','CNN1','MYOCD','MYLK')
mural_ph2       <- .prep_phase_heatmap(sn_sub, mural_subs,
                                       'RVF_vs_pRV', mural_ph2_genes)
p_H <- .make_phase_heatmap(mural_ph2, 'Mural Phase 2 (RVF vs pRV)\nIEG/stress + partial contractile reversion')

## ───────────────────────────────────────────────────────────────────────────
## Panel T: Phase 2 vascular IEG/stress pseudobulk module score
##           (NR4A1/2/3 + JUN/JUNB/FOS + HSPB1 + ATF3)
##           (script-letter I → v57-letter T)
## ───────────────────────────────────────────────────────────────────────────

biphasic_genes <- c('TIMP3','NR2F2','PALLD','BCL6','HIF3A','ANGPT2')

bi_df <- sn_sub %>%
  dplyr::filter(subtype %in% mural_subs, gene %in% biphasic_genes,
                contrast %in% c('pRV_vs_NF','RVF_vs_pRV')) %>%
  dplyr::select(gene, subtype, contrast, log2FoldChange) %>%
  tidyr::pivot_wider(names_from  = contrast,
                     values_from = log2FoldChange,
                     values_fill = 0)

bi_traj <- bi_df %>%
  dplyr::mutate(NF  = 0,
                pRV = pRV_vs_NF,
                RVF = pRV_vs_NF + RVF_vs_pRV) %>%
  dplyr::select(gene, subtype, NF, pRV, RVF) %>%
  tidyr::pivot_longer(cols      = c('NF','pRV','RVF'),
                      names_to  = 'stage',
                      values_to = 'log2FC') %>%
  dplyr::mutate(stage = factor(stage, levels = c('NF','pRV','RVF')))

p_I <- ggplot(bi_traj,
              aes(x = stage, y = log2FC,
                  group = interaction(gene, subtype),
                  color = gene, linetype = subtype)) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'grey60',
             linewidth = PS$linewidth_mm * 0.6) +
  geom_line(linewidth = PS$linewidth_mm * 1.2) +
  geom_point(size = PS$scatter_pt) +
  scale_color_brewer(palette = 'Dark2') +
  labs(x = NULL,
       y = expression(cumulative~log[2]~FC~vs.~NF),
       title = 'Biphasic mural genes') +
  theme_v52(COMP_W) +
  theme(
    plot.title      = element_text(face = 'bold'),
    legend.position = 'right',
    legend.key.size = unit(3, 'mm'),
    legend.box      = 'vertical'
  )

## ───────────────────────────────────────────────────────────────────────────
## Supplementary (legacy v54 Panel J): Mural ChEA enrichment (Phase 1 up)
##           Not in v57 main caption — kept for reference
## ───────────────────────────────────────────────────────────────────────────

mural_chea_cache <- './output/Supplementary_Figure_3/fig_s3_mural_chea.rds'
if (file.exists(mural_chea_cache)) {
  mural_chea <- readRDS(mural_chea_cache)
} else {
  mural_up <- sn_sub %>%
    dplyr::filter(subtype %in% mural_subs, contrast == 'pRV_vs_NF',
                  !is.na(padj), padj < sig_padj,
                  log2FoldChange > sig_lfc, baseMean > sig_bmean) %>%
    dplyr::pull(gene) %>% unique()
  if (length(mural_up) == 0) {
    mural_chea <- data.frame(Term = character(), Adjusted.P.value = numeric(),
                             stringsAsFactors = FALSE)
  } else {
    mural_chea <- enrichR::enrichr(mural_up, c('ChEA_2016'))$ChEA_2016
  }
  dir.create('./output', showWarnings = FALSE, recursive = TRUE)
  saveRDS(mural_chea, mural_chea_cache)
}

mural_chea_tf <- mural_chea %>%
  dplyr::mutate(TF = str_extract(Term, '^[^_ ]+')) %>%
  dplyr::filter(!is.na(TF)) %>%
  dplyr::group_by(TF) %>%
  dplyr::slice_min(Adjusted.P.value, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(Adjusted.P.value) %>%
  dplyr::slice_head(n = 8)

p_J <- .make_chea_bar(mural_chea_tf, '#c8553d', 'Mural ChEA (Ph1 up)')

## ───────────────────────────────────────────────────────────────────────────
## Final composition
## ───────────────────────────────────────────────────────────────────────────

## Display tag → variable mapping (A..M, 13 panels):
##   A p_A             FB UMAP
##   B p_B             FB markers dotplot
##   C p_C             Phase 1 erosion + priming heatmap (sn|xen)
##   D p_module_identity  Koenig 2022 identity score (Phase-1 specific)
##   E p_D             Phase 2 sn–xen lineage concordance scatter
##   F p_module_matrifib  Fu 2018 matrifibrocyte score (Phase-2 specific)
##   G p_D2            Xenium FB sublineage Phase 2 heatmap
##   H p_E             FB ChEA Ph1-down
##   I p_F             Mural UMAP
##   J p_G             Mural Phase 1 sublineage heatmap
##   K p_H             Mural Phase 2 sublineage heatmap
##   L p_I             Mural biphasic trajectory
##   M p_J             Mural ChEA Ph1-up

row1 <- (p_A | p_B)                                     + plot_layout(widths = c(1.0, 1.6))
row2 <- (p_C | p_module_identity | p_D)                 + plot_layout(widths = c(1.0, 0.9, 1.2))
row3 <- (p_module_matrifib | p_D2 | p_E)                + plot_layout(widths = c(0.9, 1.4, 1.0))
row4 <- (p_F | p_G | p_H)                               + plot_layout(widths = c(1.0, 1.0, 1.0))
row5 <- (p_I | p_J)                                     + plot_layout(widths = c(1.6, 1.0))

p_final <- (row1 / row2 / row3 / row4 / row5) +
  plot_layout(heights = c(1.0, 1.0, 1.0, 1.0, 1.0)) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(family = FONT_FAMILY, face = 'bold',
                                size = PS$tag_pt))

save_figure(p_final, 'SupplementaryFigure_3.pdf',
            width = COMP_W, height = COMP_H)

## ───────────────────────────────────────────────────────────────────────────
## ── Panel M — Adult RV trichrome / Sirius red fibrosis (% area) boxplot
## Trichrome / Sirius red fibrosis (% area) by disease state
##   - Box-and-whisker with per-sample dots
##   - KW omnibus + pairwise Wilcoxon (BH-adjusted) brackets
##   - Linear models: (a) treatment-coded contrasts (NF reference);
##                    (b) ordered group as continuous → monotonic trend
##   - Source: dependencies/shared/Adult_Trichrome_Sirius.xlsx
## ───────────────────────────────────────────────────────────────────────────

.fib_xlsx <- './dependencies/shared/Adult_Trichrome_Sirius.xlsx'
.fib_raw  <- suppressMessages(readxl::read_excel(.fib_xlsx, sheet = 1,
                                                 col_names = FALSE))
.fib_parse <- function(x) suppressWarnings(
  as.numeric(gsub('[^0-9.\\-]', '', gsub('\\*', '', x))))
.fib_vals <- as.character(.fib_raw[3, ])

fib_df <- dplyr::bind_rows(
  data.frame(group = 'NF',  pct = .fib_parse(.fib_vals[1:12])),
  data.frame(group = 'pRV', pct = .fib_parse(.fib_vals[13:24])),
  data.frame(group = 'RVF', pct = .fib_parse(.fib_vals[25:36]))) %>%
  dplyr::filter(!is.na(pct)) %>%
  dplyr::mutate(group = factor(group, levels = c('NF','pRV','RVF')))

# Stats
fib_kw <- kruskal.test(pct ~ group, data = fib_df)
fib_pw <- pairwise.wilcox.test(fib_df$pct, fib_df$group,
                               p.adjust.method = 'BH', exact = FALSE)
# Linear model: treatment-coded contrasts
fib_lm_cat   <- lm(pct ~ group, data = fib_df)
# Ordered group as continuous → monotonic trend
fib_df$group_ord <- as.numeric(fib_df$group)   # NF=1, pRV=2, RVF=3
fib_lm_trend <- lm(pct ~ group_ord, data = fib_df)

cat('\n=== Fibrosis (trichrome/Sirius % area) — group summary ===\n')
print(fib_df %>% dplyr::group_by(group) %>%
        dplyr::summarise(n=dplyr::n(),
                         median=round(median(pct),1),
                         mean=round(mean(pct),1),
                         IQR=round(IQR(pct),1)))
cat(sprintf('\nKW omnibus p = %s\n', signif(fib_kw$p.value, 3)))
cat('Pairwise Wilcoxon (BH-adjusted):\n'); print(round(fib_pw$p.value, 3))
cat('\n=== Linear model — treatment-coded (NF reference) ===\n')
print(summary(fib_lm_cat)$coefficients)
cat('\n=== Linear model — ordered group as continuous (monotonic trend) ===\n')
print(summary(fib_lm_trend)$coefficients)
cat(sprintf('  → Per-step (NF→pRV→RVF) increase = %.2f %%-area; trend p = %s; R² = %.2f\n',
            coef(fib_lm_trend)[2],
            signif(summary(fib_lm_trend)$coefficients[2,4], 3),
            summary(fib_lm_trend)$r.squared))

write.csv(fib_df,
          './output/Supplementary_Figure_3/fig_s3_adult_trichrome_sirius_long.csv', row.names = FALSE)

# Brackets for in-plot pairwise stats — positioned within the 0–50 clip
fib_brk <- data.frame(
  x    = c(1, 1, 2),
  xend = c(2, 3, 3),
  y    = c(42, 48, 45),
  wilcox_p = c(fib_pw$p.value['pRV','NF'],
               fib_pw$p.value['RVF','NF'],
               fib_pw$p.value['RVF','pRV']))
.y_step <- 1.5
fib_brk$lab <- with(fib_brk, dplyr::case_when(
  is.na(wilcox_p)  ~ 'NA',
  wilcox_p < 0.001 ~ '***',
  wilcox_p < 0.01  ~ '**',
  wilcox_p < 0.05  ~ '*',
  TRUE             ~ sprintf('p=%s', signif(wilcox_p, 2))))

p_fibrosis <- ggplot(fib_df, aes(x = group, y = pct, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.55,
               linewidth = PS$geom_lw, alpha = 0.85) +
  geom_jitter(width = 0.12, size = 1.5, shape = 21,
              color = 'black', stroke = 0.3) +
  scale_fill_manual(values = disease_pal, guide = 'none') +
  geom_segment(data = fib_brk, aes(x = x, xend = xend, y = y, yend = y),
               inherit.aes = FALSE, linewidth = PS$geom_lw) +
  geom_segment(data = fib_brk, aes(x = x, xend = x, y = y, yend = y - .y_step*0.15),
               inherit.aes = FALSE, linewidth = PS$geom_lw) +
  geom_segment(data = fib_brk, aes(x = xend, xend = xend, y = y, yend = y - .y_step*0.15),
               inherit.aes = FALSE, linewidth = PS$geom_lw) +
  geom_text(data = fib_brk,
            aes(x = (x + xend) / 2, y = y + .y_step * 0.10, label = lab),
            inherit.aes = FALSE, size = PS$text_mm,
            family = FONT_FAMILY, vjust = 0) +
  annotate('text', x = 0.55, y = 50,
           label = sprintf('KW p = %s | LM trend p = %s (slope %.1f)',
                           signif(fib_kw$p.value, 2),
                           signif(summary(fib_lm_trend)$coefficients[2,4], 2),
                           coef(fib_lm_trend)[2]),
           hjust = 0, vjust = 1, size = PS$text_mm,
           family = FONT_FAMILY, fontface = 'italic') +
  coord_cartesian(ylim = c(0, 50), clip = 'off') +
  scale_y_continuous(breaks = seq(0, 50, 10),
                     expand = expansion(mult = c(0.02, 0.02))) +
  theme_v52(COMP_W) +
  theme(plot.title = element_text(face = 'bold')) +
  xlab(NULL) + ylab('Fibrosis (% area)') +
  labs(title = 'Trichrome / Sirius red fibrosis')

save_figure(p_fibrosis,
            'SupplementaryFigure_3_preview_fibrosis_box.pdf',
            width  = 2.5 * 0.67,
            height = 2.8)

## ───────────────────────────────────────────────────────────────────────────
## Per-panel PDFs (sized to match composite proportions: 4 rows of 2.5")
## ───────────────────────────────────────────────────────────────────────────

panel_h <- COMP_H / 5   # 5-row composite

panel_specs <- list(
  list(plot = p_A,                tag = 'A', w = COMP_W * 1.0 / 2.6,        h = panel_h),
  list(plot = p_B,                tag = 'B', w = COMP_W * 1.6 / 2.6 * 1.5,  h = panel_h * 0.67 * 1.15 * 1.10),
  list(plot = p_C,                tag = 'C', w = COMP_W * 1.0 / 3.1,        h = panel_h),
  list(plot = p_module_identity,  tag = 'D', w = COMP_W * 0.9 / 3.1,        h = panel_h),
  list(plot = p_D,                tag = 'E', w = COMP_W * 1.2 / 3.1,        h = panel_h),
  list(plot = p_module_matrifib,  tag = 'F', w = COMP_W * 0.9 / 3.3,        h = panel_h),
  list(plot = p_D2,               tag = 'G', w = COMP_W * 1.4 / 3.3,        h = panel_h),
  list(plot = p_E,                tag = 'H', w = COMP_W * 1.0 / 3.3,        h = panel_h),
  list(plot = p_F,                tag = 'I', w = COMP_W * 1.0 / 3.0,        h = panel_h),
  list(plot = p_G,                tag = 'J', w = COMP_W * 1.0 / 3.0,        h = panel_h),
  list(plot = p_H,                tag = 'K', w = COMP_W * 1.0 / 3.0,        h = panel_h),
  list(plot = p_I,                tag = 'L', w = COMP_W * 1.6 / 2.6,        h = panel_h),
  list(plot = p_J,                tag = 'M', w = COMP_W * 1.0 / 2.6,        h = panel_h)
)

for (ps in panel_specs) {
  save_figure(ps$plot,
              sprintf('SupplementaryFigure_3_panel_%s.pdf', ps$tag),
              width = ps$w, height = ps$h)
}
