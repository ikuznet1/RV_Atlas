###############################################################################
## _shared_helpers.R — v52 figure generation helpers
##
## Source once at the top of each new_scripts/Figure_*.R script via:
##   source('./new_scripts/_shared_helpers.R')
##
## Provides:
##   - Consistent disease-state palette (NF / pRV / RVF)
##   - Lineage palette (12 major snRNA-seq lineages)
##   - Xenium spatial layout helper (facet_wrap by patient, rasterised scatter)
##   - Panel labelling theme (for patchwork tag_levels = 'A')
##   - save_figure() helper that writes to ./output/v52_figures/
###############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(ggrastr)
  library(RColorBrewer)
  library(dplyr)
})

## ── Publication typography constants ────────────────────────────────────────
## Final figure width in Illustrator is 17 cm.
FINAL_WIDTH_IN <- 17 / 2.54   # = 6.6929 inches

## ── Arial font registration (guarded, cross-platform) ──────────────────────
.have_showtext <- requireNamespace("showtext", quietly = TRUE) &&
                  requireNamespace("sysfonts",  quietly = TRUE)
if (.have_showtext) {
  if (!("Arial" %in% sysfonts::font_families())) {
    arial_candidates <- c(
      "/Library/Fonts/Arial.ttf",
      "/System/Library/Fonts/Supplemental/Arial.ttf",
      "/Users/ikuz/Library/Fonts/Arial.ttf",
      Sys.glob("/usr/share/fonts/*/*[Aa]rial*.ttf")
    )
    arial_path <- arial_candidates[file.exists(arial_candidates)][1]
    if (!is.na(arial_path)) sysfonts::font_add("Arial", regular = arial_path)
  }
  showtext::showtext_auto()
  showtext::showtext_opts(dpi = 300)
}
FONT_FAMILY <- if (.have_showtext && "Arial" %in% sysfonts::font_families()) "Arial" else "sans"

## ── Output directory ────────────────────────────────────────────────────────
## Default fallback only; each Figure_*.R / Supplementary_Figure_*.R script
## should override V52_FIG_DIR right after sourcing this file (e.g.
## V52_FIG_DIR <- './output/Figure_1') and call dir.create() itself.
## The directory is NOT auto-created here so that a stray './output/v52_figures/'
## isn't left behind by every script that sources the helper.
V52_FIG_DIR <- './output/v52_figures'

## ── Palettes ────────────────────────────────────────────────────────────────
## Disease state colours (NF / pRV / RVF). Keep consistent across every figure.
disease_pal <- c(
  NF  = '#4DAF4A',
  pRV = '#377EB8',
  RVF = '#E41A1C'
)

## 12-lineage palette (matches XeniumFigureGen.r lineage ordering)
lineage_pal <- c(
  CM         = '#E41A1C',
  FB         = '#377EB8',
  EC         = '#4DAF4A',
  Myeloid    = '#984EA3',
  NKT        = '#FF7F00',
  Mural      = '#A65628',
  B          = '#F781BF',
  Mast       = '#999999',
  Neuronal   = '#66C2A5',
  Adipocyte  = '#FC8D62',
  Epicardial = '#8DA0CB',
  Lymphatic  = '#E78AC3'
)

## Phase 1 / Phase 2 gene palettes for Figure 3 + Figure 6
phase1_genes <- c('FKBP5', 'TGFBR3', 'GPX3', 'VSIG4', 'CD163', 'LYVE1')
phase2_genes <- c('ACTA2', 'ANGPT2', 'SNAI1', 'GZMB')  # + BCL6 targets, added dynamically

## ── Shared ggplot theme (publication-grade, scaled for 17 cm display) ───────
## comp_width: the composite PDF width in inches (each script sets this).
## At 17 cm display, all body text renders at 7 pt Arial,
## panel tags at 16 pt Arial bold, axes/ticks at 1 pt black lines.
theme_v52 <- function(comp_width) {
  fs         <- comp_width / FINAL_WIDTH_IN
  base_pt    <- 7  * fs
  tag_pt     <- 16 * fs
  lw_mm      <- 1  * fs * 0.352778   # 1 pt line at display → mm for ggplot

  theme_bw(base_size = base_pt, base_family = FONT_FAMILY) +
  theme(
    text              = element_text(family = FONT_FAMILY, size = base_pt, color = "black"),
    plot.title        = element_text(family = FONT_FAMILY, size = base_pt, face = "plain",
                                     hjust = 0.5, color = "black"),
    plot.subtitle     = element_text(family = FONT_FAMILY, size = base_pt, color = "black"),
    axis.title        = element_text(family = FONT_FAMILY, size = base_pt, color = "black"),
    axis.text         = element_text(family = FONT_FAMILY, size = base_pt, color = "black"),
    axis.line         = element_line(linewidth = lw_mm, color = "black"),
    axis.ticks        = element_line(linewidth = lw_mm, color = "black"),
    axis.ticks.length = unit(2 * fs, "pt"),
    legend.title      = element_text(family = FONT_FAMILY, size = base_pt, color = "black"),
    legend.text       = element_text(family = FONT_FAMILY, size = base_pt, color = "black"),
    legend.key        = element_rect(fill = "white", color = NA),
    legend.key.size   = unit(0.4 * fs, 'cm'),
    strip.text        = element_text(family = FONT_FAMILY, face = 'bold',
                                     size = base_pt, color = "black"),
    strip.background  = element_rect(fill = 'white', color = "black", linewidth = lw_mm),
    plot.tag          = element_text(family = FONT_FAMILY, face = 'bold',
                                     size = tag_pt, color = "black"),
    panel.border      = element_blank(),
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_blank(),
    panel.background  = element_rect(fill = "white", color = NA),
    plot.background   = element_rect(fill = "white", color = NA)
  )
}

## ── Standardised geom scaling ──────────────────────────────────────────────
## Call pub_scales(COMP_W) once per script to get a list of scaled constants.
pub_scales <- function(comp_width) {
  fs    <- comp_width / FINAL_WIDTH_IN
  lw_mm <- fs * 0.352778
  list(
    font_scale   = fs,
    base_pt      = 7  * fs,
    tag_pt       = 16 * fs,
    linewidth_mm = lw_mm,
    geom_lw      = lw_mm,                     # geom_boxplot/violin/bar outline
    umap_pt      = 0.3 * fs,                  # UMAP / dense scatter point size
    scatter_pt   = 1.0 * fs,                  # sparse scatter point size
    xenium_pt    = 0.05 * fs,                 # Xenium spatial point size
    text_mm      = (7 * fs) / ggplot2::.pt,   # geom_text/annotate size (mm)
    bracket_lw   = lw_mm,                     # ggsignif / stat_compare_means
    dot_range    = c(1, 5) * fs,              # DotPlot dot size range
    legend_key   = unit(0.4 * fs, 'cm')
  )
}

## ── Disease palette enforcement ───────────────────────────────────────────
scale_fill_disease   <- function(...) scale_fill_manual(values = disease_pal, ...)
scale_colour_disease <- function(...) scale_colour_manual(values = disease_pal, ...)
scale_color_disease  <- scale_colour_disease

## ── Lineage palette enforcement ───────────────────────────────────────────
scale_fill_lineage   <- function(...) scale_fill_manual(values = lineage_pal, ...)
scale_colour_lineage <- function(...) scale_colour_manual(values = lineage_pal, ...)
scale_color_lineage  <- scale_colour_lineage

## ── Significance symbol convention ────────────────────────────────────────
## Threshold map used everywhere ggsignif / stat_compare_means / stat_pvalue_manual
## or manual annotations place stars on plots.
signif_cutpoints <- c(0, 1e-4, 1e-3, 1e-2, 5e-2, 1)
signif_symbols   <- c('****', '***', '**', '*', 'ns')
pub_signif_args  <- list(cutpoints = signif_cutpoints, symbols = signif_symbols)

## Convert a numeric vector of p-values to the standardised stars.
sig_stars <- function(p) {
  as.character(cut(p,
                   breaks = c(-Inf, 1e-4, 1e-3, 1e-2, 5e-2, Inf),
                   labels = signif_symbols,
                   right  = TRUE))
}

## ── Zero-touch bar/boxplot y-axis ─────────────────────────────────────────
## Replace scale_y_continuous() on bar/boxplot panels so bars sit on the axis.
## Keeps ~5% top headroom by default; pass `expand_top` to override.
scale_y_touch <- function(..., expand_top = 0.05) {
  scale_y_continuous(expand = expansion(mult = c(0, expand_top)), ...)
}

## ── Italic gene-axis helper ───────────────────────────────────────────────
## Drop-in partial theme for axes whose tick labels are gene symbols.
## Usage:   p + italic_gene_axis('y')   # y-axis labels are genes
##          p + italic_gene_axis('x')   # x-axis labels are genes (e.g. DotPlot)
italic_gene_axis <- function(axis = c('y', 'x')) {
  axis <- match.arg(axis)
  if (axis == 'y')
    theme(axis.text.y = element_text(face = 'italic', family = FONT_FAMILY))
  else
    theme(axis.text.x = element_text(face = 'italic', family = FONT_FAMILY))
}

## ── Asset directory for cropped cartoons / external panels ────────────────
ASSET_DIR <- './new_scripts/assets'
dir.create(ASSET_DIR, showWarnings = FALSE, recursive = TRUE)

## Insert a cropped asset image as a ggplot grob for patchwork composition.
## asset_file: filename inside ASSET_DIR (e.g., "Figure_1_heart_cartoon.png")
## Returns a ggdraw() object suitable for patchwork `+` or `|` operators.
insert_asset <- function(asset_file, label = NULL) {
  suppressPackageStartupMessages(library(cowplot))
  path <- file.path(ASSET_DIR, asset_file)
  if (file.exists(path)) {
    p <- ggdraw() + draw_image(path)
    if (!is.null(label)) {
      sc <- if (exists("COMP_W")) COMP_W / FINAL_WIDTH_IN else 1
      p <- p + draw_label(label, x = 0.5, y = 0.02, vjust = 0, hjust = 0.5,
                           fontfamily = FONT_FAMILY, size = 7 * sc, color = "grey40")
    }
  } else {
    sc <- if (exists("COMP_W")) COMP_W / FINAL_WIDTH_IN else 1
    p <- ggdraw() +
      draw_label(paste0("[Missing: ", asset_file, "]"),
                 fontfamily = FONT_FAMILY, size = 7 * sc, color = "grey50") +
      theme(plot.background = element_rect(fill = "grey95", color = "grey70",
                                            linewidth = 0.5, linetype = "dashed"))
  }
  p
}

## ── Spatial layout helper (Xenium tissue scatter, faceted by patient) ──────
## df must contain: patient, x_centroid, y_centroid, and a `value` column
## mapped to colour. Rasterised via ggrastr for small PDF size.
xenium_spatial_panel <- function(df, value_col,
                                 palette = NULL,
                                 title = NULL,
                                 ncol = 3,
                                 point_size = 0.05,
                                 dpi = 300) {
  stopifnot(all(c('patient', 'x_centroid', 'y_centroid', value_col) %in% colnames(df)))
  p <- ggplot(df, aes(x = x_centroid, y = -y_centroid,
                      colour = .data[[value_col]])) +
    ggrastr::rasterise(geom_point(size = point_size), dpi = dpi) +
    facet_wrap(~ patient, scales = 'free', ncol = ncol) +
    labs(title = title, x = NULL, y = NULL, colour = value_col) +
    theme_v52(comp_width = 18) +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank())
  if (!is.null(palette)) {
    p <- p + scale_colour_manual(values = palette)
  }
  p
}

## ── Figure save helper ─────────────────────────────────────────────────────
save_figure <- function(plot, name, width, height, device = cairo_pdf) {
  path <- file.path(V52_FIG_DIR, name)
  ggsave(path, plot = plot, width = width, height = height, device = device)
  message('Saved: ', path)
  invisible(path)
}

## ── Pseudobulk DESeq2 helper (in-script DEG computation) ───────────────────
## Use on any Seurat object to get fresh DESeq2 results per subtype per contrast
## without relying on pre-saved CSVs.
##   seurat_obj: Seurat with meta.data$patient and meta.data$group
##   ident_col:  metadata column defining cell identity (e.g. 'subtype')
##   contrasts:  list of pairs, e.g. list(c('RVF','NF'), c('pRV','NF'))
##   min_cells:  skip subtypes with fewer cells in any group
run_pseudobulk_deseq <- function(seurat_obj, ident_col, contrasts,
                                 min_cells = 10, assay = 'RNA') {
  suppressPackageStartupMessages({
    library(Seurat)
    library(DESeq2)
  })
  ids <- as.character(sort(unique(seurat_obj@meta.data[[ident_col]])))
  out <- list()
  for (id in ids) {
    sub <- subset(seurat_obj, cells = rownames(seurat_obj@meta.data)[
      as.character(seurat_obj@meta.data[[ident_col]]) == id
    ])
    if (ncol(sub) < min_cells) next
    agg <- AggregateExpression(sub, assays = assay,
                               group.by = c('patient', 'group'),
                               return.seurat = FALSE)[[assay]]
    meta <- data.frame(sample = colnames(agg))
    meta$patient <- sub('_.*$', '', meta$sample)
    meta$group   <- sub('^.*_', '',  meta$sample)
    meta$group   <- factor(meta$group, levels = c('NF', 'pRV', 'RVF'))
    rownames(meta) <- meta$sample
    for (cc in contrasts) {
      keep <- meta$group %in% cc
      if (sum(keep) < 4) next
      grp_counts <- table(droplevels(meta$group[keep]))
      if (length(grp_counts) < 2 || any(grp_counts < 2)) next
      md <- meta[keep, , drop = FALSE]
      md$group <- droplevels(md$group)
      ct_counts <- round(agg[, keep])
      keep_genes <- rowSums(ct_counts >= 5) >= 2
      ct_counts <- ct_counts[keep_genes, , drop = FALSE]
      if (nrow(ct_counts) < 10) next
      dds <- DESeqDataSetFromMatrix(countData = ct_counts,
                                    colData   = md,
                                    design    = ~ group)
      dds <- tryCatch(DESeq(dds, quiet = TRUE), error = function(e) NULL)
      if (is.null(dds)) next
      r <- results(dds, contrast = c('group', cc[1], cc[2]))
      res <- data.frame(
        gene           = rownames(r),
        baseMean       = as.numeric(r$baseMean),
        log2FoldChange = as.numeric(r$log2FoldChange),
        lfcSE          = as.numeric(r$lfcSE),
        stat           = as.numeric(r$stat),
        pvalue         = as.numeric(r$pvalue),
        padj           = as.numeric(r$padj),
        subtype        = as.character(id),
        contrast       = paste(cc, collapse = '_vs_'),
        stringsAsFactors = FALSE
      )
      out[[paste(id, paste(cc, collapse = '_'), sep = '__')]] <- res
    }
  }
  dplyr::bind_rows(out)
}

## ── Pan-lineage pseudobulk DESeq2 helper ───────────────────────────────────
## Like run_pseudobulk_deseq but aggregates across ALL cells (no ident split) —
## i.e. donor-level pseudobulk. Use for Pan-lineage genes (e.g. FKBP5, GPX3).
run_pseudobulk_deseq_pan <- function(seurat_obj, contrasts, assay = NULL) {
  suppressPackageStartupMessages({
    library(Seurat)
    library(DESeq2)
  })
  if (is.null(assay)) {
    av <- names(seurat_obj@assays)
    assay <- if ('RNA' %in% av) 'RNA' else av[1]
  }
  agg <- AggregateExpression(seurat_obj, assays = assay,
                             group.by = c('patient', 'group'),
                             return.seurat = FALSE)[[assay]]
  meta <- data.frame(sample = colnames(agg))
  meta$patient <- sub('_.*$', '', meta$sample)
  meta$group   <- sub('^.*_', '',  meta$sample)
  meta$group   <- factor(meta$group, levels = c('NF', 'pRV', 'RVF'))
  rownames(meta) <- meta$sample
  out <- list()
  for (cc in contrasts) {
    keep <- meta$group %in% cc
    if (sum(keep) < 4) next
    grp_counts <- table(droplevels(meta$group[keep]))
    if (length(grp_counts) < 2 || any(grp_counts < 2)) next
    md <- meta[keep, , drop = FALSE]
    md$group <- droplevels(md$group)
    ct_counts <- round(agg[, keep])
    keep_genes <- rowSums(ct_counts >= 5) >= 2
    ct_counts <- ct_counts[keep_genes, , drop = FALSE]
    if (nrow(ct_counts) < 10) next
    dds <- DESeqDataSetFromMatrix(countData = ct_counts,
                                  colData   = md,
                                  design    = ~ group)
    dds <- tryCatch(DESeq(dds, quiet = TRUE), error = function(e) NULL)
    if (is.null(dds)) next
    r <- results(dds, contrast = c('group', cc[1], cc[2]))
    res <- data.frame(
      gene           = rownames(r),
      baseMean       = as.numeric(r$baseMean),
      log2FoldChange = as.numeric(r$log2FoldChange),
      lfcSE          = as.numeric(r$lfcSE),
      stat           = as.numeric(r$stat),
      pvalue         = as.numeric(r$pvalue),
      padj           = as.numeric(r$padj),
      contrast       = paste(cc, collapse = '_vs_'),
      stringsAsFactors = FALSE
    )
    out[[paste(cc, collapse = '_')]] <- res
  }
  dplyr::bind_rows(out)
}
