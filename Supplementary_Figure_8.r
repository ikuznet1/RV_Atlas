###############################################################################
## Supplementary Figure 8 (v54 draft) -- Adult vs pediatric NF RV comparison
##
## Panels (from RV_snRNASeq_v54_draft.{md,docx} S8 legend):
##   (A) UMAP co-embedding of healthy pediatric + adult RV snRNA-seq datasets
##   (B) Cell-type abundance by dataset origin (CM/FB/EC/Myeloid box plots,
##       one-way ANOVA per cell type)
##   (C) Bulk WGCNA module expression by cell type and dataset origin
##       (DotPlot of module scores split on CombinedNames x origin)
##   (D) Volcano of CM pseudobulk DEGs (adult vs ped NF, DESeq2)
##   (E) PCA embedding of CM pseudobulk counts (PC1/PC2 by origin)
##   (F) GO Biological Process enrichment for top-ranked
##       CM-pediatric-up and CM-adult-up gene sets
##
## Source: ported from legacy ./Supplementary_Figure_8.R (repo root) into the
##         new_scripts v54 house style (pub_scales / theme_v52 / save_figure)
##         on 2026-05-10. The legacy v52 Xenium QC content that previously
##         occupied this slot has been moved to new_scripts/_legacy_xenium_qc.R.
##
## Output: ./output/Supplementary_Figure_8/v52_figures/SupplementaryFigure_8.pdf
###############################################################################

source('./helper_scripts/_shared_helpers.R')

## Per-figure output directory (introduced for consistent output paths)
V52_FIG_DIR <- './output/Supplementary_Figure_8'
dir.create(V52_FIG_DIR, showWarnings = FALSE, recursive = TRUE)


## Suppress R's default Rplots.pdf in cwd when Rscript hits a plot call
## that's outside an explicit pdf() ... dev.off() envelope.
pdf(NULL)
suppressPackageStartupMessages({
  library(Seurat)
  library(hdWGCNA)
  library(ggeasy)
  library(dplyr)
  library(EnhancedVolcano)
  library(ggpubr)
  library(DESeq2)
  library(stringr)
  library(enrichR)
  library(forcats)
  library(viridis)
  library(ggplot2)
  library(patchwork)
})

source('./helper_scripts/spatial_functions.R')

COMP_W <- 14
COMP_H <- 28

## Publication scaling constants (geom linewidths, point sizes, text mm, etc.)
PS <- pub_scales(COMP_W)

## ---- Local helpers reused by panels D/F -----------------------------------
parse_ratio <- function(ratio) {
  ratio <- sub("^\\s*", "", as.character(ratio))
  ratio <- sub("\\s*$", "", ratio)
  numerator   <- as.numeric(sub("/\\d+$", "", ratio))
  denominator <- as.numeric(sub("^\\d+/", "", ratio))
  numerator / denominator
}

wrapText <- function(x, len) {
  sapply(x, function(y) paste(strwrap(y, len), collapse = "\n"), USE.NAMES = FALSE)
}

###############################################################################
## Load merged adult+pediatric RV snRNA-seq object (NF subset)
###############################################################################
M2 <- readRDS('./dependencies/Figure_8/RV_Peds_merge_shrunk.rds')
M2$condition[M2$condition == "NF"]         <- "pRV"
M2$condition[M2$condition == "Donor"]      <- "NF"
M2$condition[M2$condition == "SystolicHF"] <- "RVF"

M2$group[is.na(M2$group)] <- M2$condition[is.na(M2$group)]

## NF-only subset used for panels A-C
M1 <- subset(M2, group == 'NF')
rm(M2); invisible(gc())

## origin in this object is logical (TRUE = Peds, FALSE = Adult).
## Coerce to character labels up-front so all panels (A/B/C/D/E/F) see the
## same Adult/Peds factor (otherwise PCA/scale_colour_manual mismatch).
M1$origin <- ifelse(as.logical(M1$origin), 'Peds', 'Adult')

###############################################################################
## Panel A -- UMAP co-embedding (adult vs pediatric NF)
###############################################################################
p_S8A <- PlotEmbedding(
  M1, group.by = 'origin',
  point_size = PS$umap_pt,
  plot_under = TRUE,
  plot_theme = umap_theme() + NoLegend(),
  raster_dpi = 400, raster_scale = 0.5
) +
  ggtitle('Adult vs pediatric NF (co-embedded)') +
  theme_v52(COMP_W) +
  theme(plot.title = element_text(family = FONT_FAMILY, size = PS$base_pt))

###############################################################################
## Panel B -- Cell-type abundance by dataset origin (per-patient frequencies,
##            one-way ANOVA per cell type)
###############################################################################
cells <- table(M1$CombinedNames, M1$patient)
cells <- cells[c('CM', 'FB', 'EC', 'Myeloid'), ]
cells <- sweep(cells, 2, colSums(cells), '/')
cells <- data.frame(cells)
cells$origin <- M1$origin[match(cells$Var2, M1$patient)]
cells$origin <- as.factor(cells$origin)

p_S8B <- ggboxplot(
  cells[length(cells$origin):1, ],
  x = "origin", y = "Freq",
  fill = "origin", group = "origin"
) +
  facet_wrap(~ Var1, ncol = 7) +
  stat_compare_means(aes(group = origin), method = "anova") +
  labs(color = 'Group', x = "Dataset origin", y = 'Frequency') +
  theme_v52(COMP_W) +
  theme(legend.key.size = PS$legend_key)

###############################################################################
## Panel C -- Bulk WGCNA module expression by cell type x origin
###############################################################################
M2 <- M1  ## reuse alias for ProjectModules call below

DefaultAssay(M2) <- 'RNA'
M2 <- FindVariableFeatures(M2)
M2 <- ScaleData(M2, block.size = 1000)

consensus_modules <- read.csv("./dependencies/shared/bulk_heart_modules.csv")
consensus_modules <- consensus_modules[, 1:3]
consensus_modules <- subset(consensus_modules, gene_name %in% rownames(M2))
consensus_modules <- consensus_modules[
  match(unique(consensus_modules$gene_name), consensus_modules$gene_name), ]

M2 <- SetupForWGCNA(
  M2,
  gene_select = "fraction",
  fraction = 0.05,
  wgcna_name = "Cardiomyocyte"
)

## NOTE: hdWGCNA + harmony 2.x are incompatible — ProjectModules passes
## `assay.use=...` to RunHarmony which harmony 2.x rejects (same as F7).
## Drop group.by.vars to skip the internal harmony re-integration; M2 is
## already harmonised when read in.
M2 <- ProjectModules(
  M2,
  modules = consensus_modules,
  seurat_ref = M2,
  wgcna_name = "Cardiomyocyte",
  wgcna_name_proj = 'bulk2sn'
)

M2 <- SetActiveWGCNA(M2, 'bulk2sn')
mapping <- labels2colors(1:100)
MEs <- GetMEs(M2, harmonized = TRUE)
mods <- colnames(MEs); mods <- mods[mods != 'grey']
all_signif <- c('M1', 'M2', 'M3', 'M4', 'M5', 'M8', 'M10',
                'M11', 'M12', 'M14', 'M20', 'M25', 'M26', 'M28')

colnames(MEs) <- paste0('M', match(colnames(MEs), mapping))
M2@meta.data  <- cbind(M2@meta.data, MEs)
M2 <- SetIdent(M2, value = "CombinedNames")

score_calc    <- consensus_modules %>% group_by(module) %>% group_split()
module_colors <- unique(unlist(lapply(score_calc, '[[', 'module')))
module_colors <- paste0('M', match(module_colors, mapping))

DefaultAssay(M2) <- 'SCT'
M2 <- AddModuleScore(M2, lapply(score_calc, '[[', 'gene_name'),
                     name = "module_score")

cols_current <- colnames(M2@meta.data)
cols_current[startsWith(colnames(M2@meta.data), 'module_score')] <-
  paste0('module_', module_colors)
colnames(M2@meta.data) <- cols_current

## origin was ALREADY coerced to 'Peds'/'Adult' on M1 above, and M2 is an alias
## of M1 — do NOT re-run as.logical() here: as.logical('Peds')/as.logical('Adult')
## are both NA, which previously made every label '<lineage>_NA'.
.origin_short <- ifelse(M2$origin == 'Peds', 'ped', 'adult')
M2$CombinedNamesSplit <- paste0(.origin_short, '-', as.character(M2$CombinedNames))
## Order rows to match the published Panel C (lineage-grouped, ped above adult).
## Seurat DotPlot puts the first factor level at the bottom, so list bottom->top.
.cns_lv <- c('adult-CM','ped-CM','adult-EC','ped-EC','adult-FB','ped-FB',
             'adult-Myeloid','ped-Myeloid','adult-PC','ped-PC','adult-SM','ped-SM')
M2$CombinedNamesSplit <- factor(
  M2$CombinedNamesSplit,
  levels = intersect(.cns_lv, unique(M2$CombinedNamesSplit)))

p_S8C <- DotPlot(
  M2,
  paste0('module_',
         c('M20','M5','M1','M3','M4','M8','M2','M12','M25','M26','M10','M28','M14','M11')),
  group.by = 'CombinedNamesSplit',
  dot.min = 0, col.min = 0, col.max = 2,
  idents = c("CM", "EC", "FB", "Myeloid", "PC", "SM")
) +
  RotatedAxis() + ylab('') + xlab('') +
  scale_color_gradient2(high = 'red', mid = 'grey95', low = 'blue') +
  theme_v52(COMP_W) +
  theme(
    panel.border = element_rect(linewidth = PS$geom_lw, fill = NA, color = 'black'),
    axis.line.x  = element_blank(),
    axis.line.y  = element_blank(),
    legend.key.size = PS$legend_key
  )

###############################################################################
## Panel D -- CM pseudobulk DESeq2 (adult vs ped NF) volcano
###############################################################################
M1 <- SetIdent(M1, value = "origin")
slot(M1$SCT@SCTModel.list[[1]], 'median_umi') <-
  median(M1$SCT@SCTModel.list[[1]]@cell.attributes$umi)

M_cm <- subset(M1, CombinedNames == 'CM')

M_pb <- AggregateExpression(M_cm, assays = "RNA",
                            return.seurat = FALSE,
                            group.by = c("origin", "patient"))
pseudo_counts <- M_pb$RNA

col_data <- data.frame(colnames(pseudo_counts))
col_data$origin <- str_split_i(col_data[, 1], '_', 1)

dds <- DESeqDataSetFromMatrix(
  countData = pseudo_counts,
  colData   = col_data,
  design    = ~ origin
)
dds <- DESeq(dds)
res <- results(dds)
res$padj[res$padj < 1e-50] <- 1e-50

p_S8D <- EnhancedVolcano(
  res, lab = rownames(res),
  x = 'log2FoldChange', y = 'padj',
  FCcutoff = 1, xlim = c(-10, 10),
  title = 'CM pseudobulk: adult vs pediatric NF',
  subtitle = NULL
) +
  theme_v52(COMP_W) +
  theme(legend.position = 'top', legend.key.size = PS$legend_key)

###############################################################################
## Panel E -- PCA of CM pseudobulk (vst transform)
###############################################################################
vsdata <- vst(dds, blind = FALSE)

pca_obj <- plotPCA(vsdata, intgroup = "origin", returnData = TRUE)
pct_var <- round(100 * attr(pca_obj, "percentVar"))

p_S8E <- ggplot(pca_obj, aes(x = PC1, y = PC2, color = origin)) +
  geom_point(size = PS$scatter_pt * 3) +   # PS$pt_size doesn't exist; use scatter_pt
  scale_colour_manual(values = c(Adult = '#E41A1C', Peds = '#377EB8'),
                      name = NULL) +
  labs(x = paste0('PC1 (', pct_var[1], '%)'),
       y = paste0('PC2 (', pct_var[2], '%)'),
       title = 'PCA: CM pseudobulk (adult vs ped NF)') +
  theme_v52(COMP_W) +
  theme(legend.key.size = PS$legend_key)

###############################################################################
## Panel F -- GO BP enrichment for CM-pediatric-up and CM-adult-up gene sets
###############################################################################
g <- read.csv('./dependencies/shared/human2mouse.csv', header = FALSE)$V1

peds_cm_up  <- rownames(subset(res, log2FoldChange >  1 & padj < 0.05))
adult_cm_up <- rownames(subset(res, log2FoldChange < -1 & padj < 0.05))

peds_cm_up  <- peds_cm_up[!startsWith(peds_cm_up,  'MT-')]
peds_cm_up  <- peds_cm_up[peds_cm_up  %in% g]
adult_cm_up <- adult_cm_up[!startsWith(adult_cm_up, 'MT-')]
adult_cm_up <- adult_cm_up[adult_cm_up %in% g]

dbs <- c("ChEA_2022", "WikiPathway_2023_Human",
         "Reactome_2016", "GO_Biological_Process_2025")

.go_dot <- function(gene_set, title) {
  enriched <- enrichr(gene_set, dbs)
  go_bp <- subset(enriched[[4]], Adjusted.P.value < 0.05)
  if (nrow(go_bp) == 0) {
    return(ggplot() + ggtitle(paste0(title, ' (no GO BP hits)')) +
             theme_v52(COMP_W))
  }
  go_bp <- go_bp[order(go_bp$Combined.Score, decreasing = TRUE), ]
  top3  <- go_bp[rev(seq_len(min(3, nrow(go_bp)))), ]
  ggplot(top3,
         aes(x = Combined.Score,
             y = fct_inorder(Term),
             color = as.numeric(Adjusted.P.value),
             size  = parse_ratio(Overlap))) +
    geom_point() +
    xlab('Combined Score') + ylab('Term') +
    labs(color = "P value", size = "Overlap", title = title) +
    scale_y_discrete(labels = fct_inorder(wrapText(
      sapply(strsplit(top3$Term, " \\(GO"), `[`, 1), 35))) +
    scale_color_stepsn(colors = rev(magma(256))) +
    theme_v52(COMP_W) +
    theme(axis.text = element_text(colour = "black"),
          legend.key.size = PS$legend_key)
}

p_S8F_peds  <- .go_dot(peds_cm_up,  'CM pediatric up')
p_S8F_adult <- .go_dot(adult_cm_up, 'CM adult up')
p_S8F <- p_S8F_peds | p_S8F_adult

###############################################################################
## Assemble Supplementary Figure 8
###############################################################################
fig_S8 <- (p_S8A / p_S8B / p_S8C / p_S8D / p_S8E / p_S8F) +
  plot_annotation(
    tag_levels = 'A',
    theme = theme(plot.tag = element_text(
      family = FONT_FAMILY, size = PS$tag_pt, face = "bold"))
  )

save_figure(fig_S8, 'SupplementaryFigure_8.pdf', width = 14, height = 28)

## Individual subpanels (same plot objects, saved separately for assembly)
save_figure(p_S8A,       'Supplementary_Figure_8_panel_A_umap.pdf',           width =  7, height = 6)
save_figure(p_S8B,       'Supplementary_Figure_8_panel_B_abundance.pdf',      width = 14, height = 4)
save_figure(p_S8C,       'Supplementary_Figure_8_panel_C_module_dotplot.pdf', width = 10, height = 6)
save_figure(p_S8D,       'Supplementary_Figure_8_panel_D_volcano.pdf',        width =  8, height = 8)
save_figure(p_S8E,       'Supplementary_Figure_8_panel_E_pca.pdf',            width =  6, height = 6)
save_figure(p_S8F,       'Supplementary_Figure_8_panel_F_GO.pdf',             width = 14, height = 5)
save_figure(p_S8F_peds,  'Supplementary_Figure_8_panel_F_GO_peds.pdf',        width =  7, height = 5)
save_figure(p_S8F_adult, 'Supplementary_Figure_8_panel_F_GO_adult.pdf',       width =  7, height = 5)
message('Supplementary Figure 8 (v54: adult vs pediatric NF) complete — combined + 8 subpanels.')
