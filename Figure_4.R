###############################################################################
## Figure 4 (v53 final) — Two-phase model of RV disease progression
##
## Panels (final, derived from new_scripts/Figure_4.png):
##   (A) Two-phase schematic (Illustrator asset): Phase 1 (NF→pRV) erosion of
##       protective programs vs Phase 2 (pRV→RVF) engagement of pathological
##       programs, with per-lineage gene call-outs
##       (Pan, CM, Myeloid, EC, FB, Mural)
##   (B) Cross-platform DEG heatmap — key Phase 1 + Phase 2 genes × lineages ×
##       (Bulk, sn, Xenium) × (pRV vs NF, RVF vs NF); log₂FC colour with
##       significance stars
##   (C) Bulk WGCNA module × phase × platform × lineage Z-heatmap
##       (AddModuleScore), facetted by Pan / CM / Myeloid / EC / FB / Mural
##   (D) Phase 1 vs Phase 2 DEG-count bar chart by platform (Bulk, sn, Xenium)
##   (E) CellChat top-15 signaling pathways, Myocardium niche (generated below
##       from the Supplementary Data S14 table)
##   (F) CellChat normalized cell-type communication chord, Myocardium (generated
##       below from the S14 table + cached Myocardium cell-type proportions)
##
## Outputs (./output/Figure_4/v52_figures/):
##   Figure_4_panel_B.pdf, Figure_4_panel_C.pdf, Figure_4_panel_D.pdf,
##   Figure_4_panel_E.pdf, Figure_4_panel_F.pdf, Figure_4.pdf (composite)
##   (Panel A is an Illustrator asset, no PDF written.)
###############################################################################

source('./helper_scripts/_shared_helpers.R')

## Per-figure output directory (introduced for consistent output paths)
V52_FIG_DIR <- './output/Figure_4'
dir.create(V52_FIG_DIR, showWarnings = FALSE, recursive = TRUE)


## Suppress R's default Rplots.pdf in cwd when Rscript hits a plot call
## that's outside an explicit pdf() ... dev.off() envelope.
pdf(NULL)
COMP_W <- 14
COMP_H <- 18
PS <- pub_scales(COMP_W)

## Panel A is an Illustrator-authored two-phase schematic. The PNG asset is
## loaded via insert_asset('Figure_4_panel_A_twophase.png') below.

suppressPackageStartupMessages({
  library(Seurat)
  library(DESeq2)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(scales)
  library(ggrastr)
  library(tidyr)
  library(tibble)
  library(ggpattern)
})

###############################################################################
## Load snRNA-seq Seurat object (post annotate_subnames.r, has Subnames_manual)
###############################################################################
sn <- readRDS('./dependencies/shared/snRV_ref.rds')

## Build lineage → subtype lookup (Level1 = Names, Level2 = Subnames_manual)
meta <- sn@meta.data
lineage_map <- setNames(meta$Names, meta$Subnames_manual)

## Phase gene sets (legacy, kept for reference)
phase1_genes_fig4 <- c('FKBP5', 'GPX3', 'MYH6', 'TGFBR3', 'VSIG4', 'CD163', 'SFRP4', 'THBS4')
phase2_genes_fig4 <- c('ANGPT2', 'SNAI1', 'POSTN', 'COL1A1', 'THBS4', 'GNLY', 'S100A1', 'MYBPC3')

## Unified heatmap gene groups (Panel B/C combined) — rows grouped by category
## Pan genes query donor-level pseudobulk; lineage genes query only their lineage.
fig4_gene_groups <- list(
  Pan     = c('FKBP5', 'GPX3'),
  CM      = c('MYH6', 'HMGCS2', 'CREB5', 'HIF1A', 'NAMPT', 'BMPR2'),
  Myeloid = c('VSIG4', 'CD163', 'LYVE1', 'MARCO', 'CD83', 'CX3CR1', 'FCER1A', 'CD1C'),
  EC      = c('ACE2', 'IL1RL1', 'ANGPT2', 'SNAI1', 'MYC', 'CXCL9', 'IFI44L', 'RSAD2'),
  FB      = c('TGFBR3', 'ACTA2', 'SFRP4', 'TNC', 'ELN', 'THBS4', 'KIT'),
  PC      = c('NNMT', 'ADH1B', 'SCGN', 'CNN1', 'NR4A3', 'TIMP3', 'IL1RL1', 'PLA2G2A')
)

## Flat gene → category map (preserves group order)
gene_primary_lineage <- unlist(lapply(names(fig4_gene_groups), function(k)
  setNames(rep(k, length(fig4_gene_groups[[k]])), fig4_gene_groups[[k]])))
gene_primary_lineage <- as.list(gene_primary_lineage)

## Legacy filter (retained for subtype-level heatmaps if needed elsewhere)
gene_lineage_filter <- list(
  FKBP5  = NULL,
  TGFBR3 = c('FB'),
  GPX3   = NULL,
  VSIG4  = c('Myeloid'),
  CD163  = c('Myeloid'),
  MYH6   = c('CM'),
  SFRP4  = c('FB'),
  THBS4  = c('FB'),
  ANGPT2 = c('EC'),
  SNAI1  = c('EC'),
  POSTN  = c('FB'),
  COL1A1 = c('FB'),
  GNLY   = c('NKT'),
  S100A1 = c('CM'),
  MYBPC3 = c('CM')
)

## Xenium panel membership 
xenium_panel_genes <- read.csv('./dependencies/shared/panel_genes.csv')[,2]

###############################################################################
## Load consolidated lineage- and sublineage-level DESeq2 outputs (v53 unified)
###############################################################################
contrasts_fig4 <- list(c('pRV', 'NF'), c('RVF', 'NF'), c('RVF', 'pRV'))

## Lineage-level (12 cell types) — used by Panels B, C, D, G
## Xenium DE was regenerated from panel-only object (477 genes, not imputed 7,504).
snRNA_L1   <- read.csv('./dependencies/shared/sn_pseudobulk_lineage_deseq2.csv')
xenium_L1  <- read.csv('./dependencies/shared/Xenium/xenium_pseudobulk_lineage_deseq2.csv')

## Pan-lineage (donor-level, no lineage split) — used for Pan genes in unified B/C
.cache_sn_pan  <- './output/Figure_4/fig4_sn_pan_deseq.rds'
.cache_xen_pan <- './output/Figure_4/fig4_xen_pan_deseq.rds'

if (file.exists(.cache_sn_pan)) {
  message('Loading cached snRNA pan-lineage pseudobulk...')
  snRNA_pan <- readRDS(.cache_sn_pan)
} else {
  message('Running snRNA pan-lineage pseudobulk DESeq2...')
  snRNA_pan <- run_pseudobulk_deseq_pan(sn, contrasts = contrasts_fig4)
  dir.create(dirname(.cache_sn_pan), showWarnings = FALSE, recursive = TRUE)
  saveRDS(snRNA_pan, .cache_sn_pan)
}

xen_obj_file <- './dependencies/shared/Xenium/xenium_obj_subclustered.rds'
if (file.exists(.cache_xen_pan)) {
  message('Loading cached Xenium pan-lineage pseudobulk...')
  xenium_pan <- readRDS(.cache_xen_pan)
} else if (file.exists(xen_obj_file)) {
  message('Loading Xenium object for pan-lineage pseudobulk...')
  xen <- readRDS(xen_obj_file)
  xenium_pan <- run_pseudobulk_deseq_pan(xen, contrasts = contrasts_fig4)
  saveRDS(xenium_pan, .cache_xen_pan)
  rm(xen); gc()
} else {
  message('Xenium obj not found; Xenium pan lookups will be NA.')
  xenium_pan <- data.frame(gene = character(), log2FoldChange = numeric(),
                           padj = numeric(), contrast = character())
}

## Pan lookups also restricted to xenium panel
xenium_pan <- xenium_pan[xenium_pan$gene %in% xenium_panel_genes, ]

## Sublineage (subtypes) — used by Panel E forest plot (TGFBR3 across FB subtypes)
.cache_deseq_L2 <- './output/Figure_4/fig4_deseq_L2_cache.rds'
if (file.exists(.cache_deseq_L2)) {
  message('Loading cached deseq_L2...')
  deseq_L2 <- readRDS(.cache_deseq_L2)
} else {
  message('Running Level-2 pseudobulk DESeq2 (subtypes)...')
  deseq_L2 <- run_pseudobulk_deseq(sn, ident_col = 'Subnames_manual',
                                    contrasts = contrasts_fig4, min_cells = 10)
  saveRDS(deseq_L2, .cache_deseq_L2)
}
deseq_L2$lineage <- lineage_map[deseq_L2$subtype]

## Bulk RNA-seq DEGs (Phase 1 = NF_vs_pRV pos lfc = down in pRV; Phase 2 = pRV_vs_RVF)
bulk_ph1 <- read.csv('./dependencies/shared/NF_vs_pRV_deseq.csv')
bulk_ph2 <- read.csv('./dependencies/shared/pRV_vs_RVF_deseq.csv')

## Normalize bulk column names so gene is in 'gene' and lfc/padj are consistent
.norm_bulk <- function(d) {
  if (!'gene' %in% colnames(d) && 'X' %in% colnames(d)) d$gene <- d$X
  if (!'gene' %in% colnames(d) && 'symbol' %in% colnames(d)) d$gene <- d$symbol
  if (!'log2FoldChange' %in% colnames(d) && 'logFC' %in% colnames(d))
    d$log2FoldChange <- d$logFC
  d
}
bulk_ph1 <- .norm_bulk(bulk_ph1)
bulk_ph2 <- .norm_bulk(bulk_ph2)

###############################################################################
## Panel A — Two-phase model schematic (Illustrator-authored)
###############################################################################
p_4A <- if (file.exists(file.path(ASSET_DIR, 'Figure_4_panel_A_twophase.png'))) {
  insert_asset('Figure_4_panel_A_twophase.png')
} else {
  ggplot() +
    annotate('text', x = 0.5, y = 0.5,
             label = '[Figure 4A: Two-phase model schematic\n(designed in Illustrator)]',
             size = PS$text_mm, family = FONT_FAMILY,
             colour = 'grey50', hjust = 0.5, vjust = 0.5) +
    theme_void() +
    theme(panel.border = element_rect(colour = 'grey80', fill = NA,
                                      linewidth = PS$geom_lw))
}

###############################################################################
## Panel B — Lineage-level key-gene heatmap (Phase 1 + Phase 2)
## Rows = genes grouped by category (Pan / CM / Myeloid / EC / FB / Mural / NKT)
## Cols = Bulk / snRNA / Xenium × (pRV vs NF, RVF vs pRV)
## Pan genes use donor-level pseudobulk (snRNA_pan, xenium_pan).
## Lineage-specific genes use only their lineage row in snRNA_L1 / xenium_L1.
## Xenium cells grey if gene off-panel.
###############################################################################
.lookup_sn <- function(gene, category, contrast_label) {
  if (category == 'Pan') {
    r <- snRNA_pan %>% dplyr::filter(gene == !!gene, contrast == !!contrast_label)
  } else {
    r <- snRNA_L1 %>% dplyr::filter(gene == !!gene, subtype == !!category,
                                     contrast == !!contrast_label)
  }
  if (nrow(r) == 0) return(c(lfc = NA_real_, padj = NA_real_))
  c(lfc = r$log2FoldChange[1], padj = r$padj[1])
}
.lookup_xen <- function(gene, category, contrast_label) {
  if (length(xenium_panel_genes) > 0 && !(gene %in% xenium_panel_genes))
    return(c(lfc = NA_real_, padj = NA_real_, off = 1))
  if (category == 'Pan') {
    r <- xenium_pan %>% dplyr::filter(gene == !!gene, contrast == !!contrast_label)
  } else {
    r <- xenium_L1 %>% dplyr::filter(gene == !!gene, subtype == !!category,
                                      contrast == !!contrast_label)
  }
  if (nrow(r) == 0) return(c(lfc = NA_real_, padj = NA_real_, off = 0))
  c(lfc = r$log2FoldChange[1], padj = r$padj[1], off = 0)
}
.lookup_bulk <- function(gene, bulk_df, invert = FALSE) {
  r <- bulk_df[bulk_df$gene == gene, ]
  if (nrow(r) == 0) return(c(lfc = NA_real_, padj = NA_real_))
  lfc <- r$log2FoldChange[1]
  if (invert) lfc <- -lfc  # NF_vs_pRV pos = down in pRV → invert for pRV_vs_NF
  c(lfc = lfc, padj = r$padj[1])
}

## Build long-format unified table: gene × category × phase × platform → lfc/padj
.build_unified <- function(gene_groups) {
  rows <- list()
  for (cat in names(gene_groups)) {
    for (g in gene_groups[[cat]]) {
      for (ph in c(1, 2)) {
        sn_contrast  <- if (ph == 1) 'pRV_vs_NF' else 'RVF_vs_pRV'
        xen_contrast <- sn_contrast
        phase_label  <- if (ph == 1) 'Phase 1 (pRV vs NF)' else 'Phase 2 (RVF vs pRV)'

        s <- .lookup_sn(g, cat, sn_contrast)
        x <- .lookup_xen(g, cat, xen_contrast)
        b <- if (ph == 1) .lookup_bulk(g, bulk_ph1, invert = TRUE)
             else          .lookup_bulk(g, bulk_ph2, invert = FALSE)

        rows[[length(rows)+1]] <- data.frame(
          gene = g, category = cat, phase = phase_label, platform = 'Bulk',
          lfc = b['lfc'], padj = b['padj'], off_panel = FALSE,
          stringsAsFactors = FALSE)
        rows[[length(rows)+1]] <- data.frame(
          gene = g, category = cat, phase = phase_label, platform = 'snRNA',
          lfc = s['lfc'], padj = s['padj'], off_panel = FALSE,
          stringsAsFactors = FALSE)
        rows[[length(rows)+1]] <- data.frame(
          gene = g, category = cat, phase = phase_label, platform = 'Xenium',
          lfc = x['lfc'], padj = x['padj'],
          off_panel = isTRUE(x['off'] == 1),
          stringsAsFactors = FALSE)
      }
    }
  }
  do.call(rbind, rows)
}

.make_unified_heatmap <- function(df, gene_groups, sig_threshold = 0.05) {
  df <- df %>% dplyr::mutate(
    l2fc_capped = pmin(pmax(lfc, -4), 4),
    sig         = !is.na(padj) & padj < sig_threshold,
    note        = ifelse(sig, '\u2217', ''),
    label       = ifelse(off_panel, 'n/a', note)
  )
  gene_order <- unique(unlist(gene_groups, use.names = FALSE))
  df$gene     <- factor(df$gene, levels = rev(gene_order))
  df$category <- factor(df$category, levels = names(gene_groups))
  df$platform <- factor(df$platform, levels = c('Bulk', 'snRNA', 'Xenium'))
  df$phase    <- factor(df$phase,
                        levels = c('Phase 1 (pRV vs NF)', 'Phase 2 (RVF vs pRV)'))

  ggplot(df, aes(x = platform, y = gene)) +
    geom_tile(aes(fill = l2fc_capped), colour = 'grey70',
              linewidth = PS$geom_lw) +
    geom_tile(data = df[df$off_panel, ],
              fill = 'grey80', colour = 'grey70',
              linewidth = PS$geom_lw) +
    geom_text(aes(label = label), size = PS$text_mm,
              hjust = 0.5, vjust = 0.5,
              colour = 'black', family = FONT_FAMILY) +
    scale_fill_gradientn(
      colours  = c('#2166AC', 'white', '#B2182B'),
      limits   = c(-4, 4),
      na.value = 'grey90',
      name     = 'log\u2082FC',
      guide    = guide_colorbar(
        direction      = 'vertical',
        barwidth       = unit(0.3 * PS$font_scale, 'cm'),
        barheight      = unit(2.0 * PS$font_scale, 'cm'),
        title.position = 'top',
        ticks.colour   = 'black'
      )
    ) +
    scale_x_discrete(expand = c(0, 0), position = 'top') +
    scale_y_discrete(expand = c(0, 0)) +
    facet_grid(category ~ phase, scales = 'free_y', space = 'free_y',
               switch = 'y',
               labeller = labeller(category = c(Pan = 'Pan', CM = 'CM',
                                                Myeloid = 'Myeloid', EC = 'EC',
                                                FB = 'FB', PC = 'Mural'))) +
    labs(x = NULL, y = NULL) +
    theme_v52(COMP_W) +
    theme(
      axis.text.x      = element_text(angle = 0, hjust = 0.5, size = PS$base_pt),
      axis.text.y      = element_text(face = 'italic', size = PS$base_pt),
      strip.text.y.left = element_text(angle = 0, size = PS$base_pt, face = 'bold'),
      strip.text.x     = element_text(size = PS$base_pt, face = 'bold'),
      strip.placement  = 'outside',
      strip.background = element_rect(fill = 'grey95', colour = NA),
      legend.position  = 'right',
      panel.border     = element_rect(colour = 'black', fill = NA,
                                      linewidth = PS$geom_lw),
      panel.spacing.x  = unit(0.4, 'lines'),
      panel.spacing.y  = unit(0.15, 'lines'),
      axis.ticks       = element_blank()
    )
}

###############################################################################
## Panels B+C — Unified cross-platform × two-phase heatmap
## Rows grouped by category (Pan / CM / Myeloid / EC / FB / Mural / NKT)
## Cols: Bulk / snRNA / Xenium × Phase 1 (pRV vs NF) / Phase 2 (RVF vs pRV)
###############################################################################
## Panel B — validated combined two-phase cross-platform heatmap.
## Rows = genes by lineage (Pan/CM/Myeloid/EC/FB/Mural); cols = Bulk/snRNA/Xenium
## × Phase 1 (NF→pRV) | Phase 2 (pRV→RVF). log2FC fill, FDR star (∗), off-panel
## Xenium = grey, low-count (baseMean<20) = blank grey. Bulk files NF_vs_pRV &
## pRV_vs_RVF are inverted (×−1) so + = up in the later stage; sn/Xenium use
## contrasts pRV_vs_NF / RVF_vs_pRV. Pan row sn/Xen = mean across 5 non-CM
## lineages. Mirrors Figure_4_panel_B_crossplatform.R (edit there + here in sync).
## Wrapped in local() with its own +20% font scale (PS_B) so it does not touch
## the global COMP_W/PS used by panels C–F. Reads ONLY ./dependencies/shared.
p_4BC <- local({
  PS_B  <- pub_scales(16.8)                          # 14 × 1.2 → fonts +20%
  DEP   <- './dependencies/shared'; SIG <- 0.05; MINBM <- 20
  NONCM <- c('Myeloid','EC','FB','PC','SM')

  gene_groups <- list(
    Pan     = c('GPX3','TXNRD1','GCLM','GSR','FKBP5','KLF9'),
    CM      = c('MYH6','MLIP','NPPA','NPPB','ACTA1',  'BMPR2','BMPR1A','BAMBI','FSTL3'),
    Myeloid = c('CD163','VSIG4','MERTK','CD1C','FCER1A','CD83'),
    EC      = c('IL1RL1','ACE2','TIMP3','TIMP4','CXCL9','IFI44L','RSAD2','MX1',  'ANGPT2','SNAI1','ADAMTS1'),
    FB      = c('SFRP4','TGFBR3','THBS4','ACTA2',  'TNC','ELN'),
    PC      = c('IL1RL1','ELN','DES','TIMP3','TIMP4','ACE2','SMAD3','ADH1B')   # facet 'Mural'
  )
  cat_lineage  <- c(CM='CM', Myeloid='Myeloid', EC='EC', FB='FB', PC='PC')
  union_genes  <- unique(unlist(gene_groups, use.names = FALSE))

  bd <- function(file) {                              # bulk: invert -> later-stage orientation
    d <- read.csv(file.path(DEP, file)); names(d)[1] <- 'gene'
    d <- d[d$gene %in% union_genes, ]
    setNames(Map(function(l,p,m) c(lfc=-l, padj=p, bm=m), d$log2FoldChange, d$padj, d$baseMean), d$gene)
  }
  bd1 <- bd('NF_vs_pRV_deseq.csv'); bd2 <- bd('pRV_vs_RVF_deseq.csv')
  loadc <- function(file, con) {
    d <- read.csv(file.path(DEP, file))
    d[d$contrast == con & d$gene %in% union_genes, c('gene','subtype','log2FoldChange','padj','baseMean')]
  }
  sn1 <- loadc('sn_pseudobulk_lineage_deseq2.csv','pRV_vs_NF')
  sn2 <- loadc('sn_pseudobulk_lineage_deseq2.csv','RVF_vs_pRV')
  xe1 <- loadc('xenium_pseudobulk_lineage_deseq2.csv','pRV_vs_NF')
  xe2 <- loadc('xenium_pseudobulk_lineage_deseq2.csv','RVF_vs_pRV')
  panel_genes <- read.csv(file.path(DEP,'panel_genes.csv'))[, 2]

  lk_lin <- function(d, g, lin) { r <- d[d$gene==g & d$subtype==lin, ]
    if (!nrow(r)) c(lfc=NA, padj=NA, bm=NA) else c(lfc=r$log2FoldChange[1], padj=r$padj[1], bm=r$baseMean[1]) }
  lk_pan <- function(d, g) { r <- d[d$gene==g & d$subtype %in% NONCM, ]
    if (!nrow(r)) c(lfc=NA, sig=0, bm=NA) else
      c(lfc=mean(r$log2FoldChange, na.rm=TRUE), sig=as.numeric(sum(!is.na(r$padj) & r$padj<SIG)>=2),
        bm=mean(r$baseMean, na.rm=TRUE)) }

  phases <- list(
    list(lab='Phase 1\n(NF→pRV)',  bd=bd1, sn=sn1, xe=xe1),
    list(lab='Phase 2\n(pRV→RVF)', bd=bd2, sn=sn2, xe=xe2)
  )

  rows <- list()
  add  <- function(g,cat,ph,plat,lfc,sig,bm,off=FALSE)
    rows[[length(rows)+1]] <<- data.frame(gene=g, category=cat, phase=ph,
      platform=plat, lfc=unname(lfc), sig=sig, bm=unname(bm), off_panel=off,
      stringsAsFactors=FALSE)

  for (cat in names(gene_groups)) for (g in gene_groups[[cat]]) for (P in phases) {
    b <- P$bd[[g]]; if (is.null(b)) b <- c(lfc=NA, padj=NA, bm=NA)
    add(g, cat, P$lab, 'Bulk', b['lfc'], !is.na(b['padj']) && b['padj']<SIG, b['bm'])
    off <- !(g %in% panel_genes)
    if (cat == 'Pan') {
      s <- lk_pan(P$sn, g); add(g, cat, P$lab, 'snRNA', s['lfc'], s['sig']==1, s['bm'])
      x <- lk_pan(P$xe, g); add(g, cat, P$lab, 'Xenium', if (off) NA else x['lfc'], (!off) && x['sig']==1,
                                if (off) NA else x['bm'], off)
    } else {
      lin <- cat_lineage[[cat]]
      s <- lk_lin(P$sn, g, lin); add(g, cat, P$lab, 'snRNA', s['lfc'], !is.na(s['padj']) && s['padj']<SIG, s['bm'])
      x <- lk_lin(P$xe, g, lin); add(g, cat, P$lab, 'Xenium', if (off) NA else x['lfc'],
          (!off) && !is.na(x['padj']) && x['padj']<SIG, if (off) NA else x['bm'], off)
    }
  }
  df <- do.call(rbind, rows)

  df <- df %>% dplyr::mutate(                          # baseMean<20 -> blank grey
    low         = !off_panel & (is.na(bm) | bm < MINBM),
    sig         = sig & !low,
    l2fc_capped = ifelse(low, NA_real_, pmin(pmax(lfc, -4), 4)),
    label       = ifelse(off_panel, '', ifelse(sig, '∗', ''))
  )

  key_levels  <- unlist(lapply(names(gene_groups), function(c) paste(c, gene_groups[[c]], sep='|')), use.names=FALSE)
  df$rowkey   <- factor(paste(df$category, df$gene, sep='|'), levels = rev(key_levels))  # per-facet order
  df$category <- factor(df$category, levels = names(gene_groups))
  df$platform <- factor(df$platform, levels = c('Bulk','snRNA','Xenium'))
  df$phase    <- factor(df$phase, levels = c('Phase 1\n(NF→pRV)','Phase 2\n(pRV→RVF)'))

  ggplot(df, aes(platform, rowkey)) +
    geom_tile(aes(fill = l2fc_capped), colour='grey70', linewidth=PS_B$geom_lw) +
    geom_tile(data = df[df$off_panel, ], fill='grey80', colour='grey70', linewidth=PS_B$geom_lw) +
    geom_text(aes(label = label), size=PS_B$text_mm, colour='black', family=FONT_FAMILY,
              hjust=0.5, vjust=0.5) +
    scale_fill_gradientn(colours=c('#2166AC','white','#B2182B'), limits=c(-4,4),
      na.value='grey90', name='log₂FC',
      guide=guide_colorbar(direction='horizontal', barwidth=unit(2.4*PS_B$font_scale,'cm'),
        barheight=unit(0.3*PS_B$font_scale,'cm'), title.position='top', ticks.colour='black')) +
    scale_x_discrete(expand=c(0,0), position='top') +
    scale_y_discrete(expand=c(0,0), labels=function(x) sub('.*\\|','',x)) +
    facet_grid(category ~ phase, scales='free_y', space='free_y', switch='y',
      labeller=labeller(category=c(Pan='Pan',CM='CM',Myeloid='Myeloid',EC='EC',FB='FB',PC='Mural'))) +
    labs(x=NULL, y=NULL) +
    theme_v52(16.8) +
    theme(
      axis.text.x        = element_text(angle=0, hjust=0.5, size=PS_B$base_pt),
      axis.text.y        = element_text(size=PS_B$base_pt),
      strip.text.y.left  = element_text(angle=0, size=PS_B$base_pt, face='bold'),
      strip.text.x       = element_text(size=PS_B$base_pt, face='bold'),
      strip.placement    = 'outside',
      strip.background    = element_rect(fill='grey95', colour=NA),
      legend.position     = 'bottom',
      panel.border        = element_rect(colour='black', fill=NA, linewidth=PS_B$geom_lw),
      panel.spacing.x     = unit(0.4,'lines'),
      panel.spacing.y     = unit(0.15,'lines'),
      axis.ticks          = element_blank()
    )
})
save_figure(p_4BC, 'Figure_4_panel_B.pdf', width = 7.2, height = 11.8)
p_4B <- p_4BC  # Panel B = validated combined Phase 1 | Phase 2 cross-platform heatmap
## p_4C is assigned later from the WGCNA panel (was Panel G → now Panel C).

###############################################################################
## Panel D — DEG burden by phase × platform (|lfc|>0.5, padj<0.05, baseMean>50)
###############################################################################
.count_degs <- function(d, phase_label, platform_label, invert = FALSE) {
  base_col <- intersect(c('baseMean'), colnames(d))
  if (length(base_col) == 0) d$baseMean <- Inf
  tested <- d %>% dplyr::filter(!is.na(padj), baseMean > 50)
  sig <- tested %>% dplyr::filter(padj < 0.05, abs(log2FoldChange) > 0.5)
  ## Count UNIQUE gene symbols (deduplicated across cell types) so the figure
  ## matches its caption and the body text. Bulk has one row per gene; sn/xen
  ## are lineage-stratified, so we collapse on the gene column.
  gene_col <- intersect(c('gene', 'symbol', 'X'), colnames(sig))[1]
  n_unique <- if (!is.na(gene_col)) length(unique(sig[[gene_col]])) else nrow(sig)
  data.frame(platform = platform_label, phase = phase_label,
             n = n_unique, universe = nrow(tested),
             stringsAsFactors = FALSE)
}

burden <- dplyr::bind_rows(
  .count_degs(snRNA_L1  %>% dplyr::filter(contrast == 'pRV_vs_NF'),  'Phase 1', 'snRNA'),
  .count_degs(snRNA_L1  %>% dplyr::filter(contrast == 'RVF_vs_pRV'), 'Phase 2', 'snRNA'),
  .count_degs(xenium_L1 %>% dplyr::filter(contrast == 'pRV_vs_NF'),  'Phase 1', 'Xenium'),
  .count_degs(xenium_L1 %>% dplyr::filter(contrast == 'RVF_vs_pRV'), 'Phase 2', 'Xenium'),
  .count_degs(bulk_ph1, 'Phase 1', 'Bulk'),
  .count_degs(bulk_ph2, 'Phase 2', 'Bulk')
)
burden$platform <- factor(burden$platform, levels = c('Bulk', 'snRNA', 'Xenium'))
burden$pct <- 100 * burden$n / burden$universe

# Long form: show raw counts and % of tested gene-lineage universe side-by-side
burden_long <- dplyr::bind_rows(
  transform(burden, metric = 'Count',       value = n,
            label = formatC(n, format = 'd', big.mark = ',')),
  transform(burden, metric = '% of tested', value = pct,
            label = sprintf('%.1f%%', pct))
)
burden_long$metric <- factor(burden_long$metric,
                             levels = c('Count', '% of tested'))

p_4D <- ggplot(burden_long, aes(x = platform, y = value, fill = phase)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7,
           linewidth = PS$geom_lw, colour = 'black') +
  geom_text(aes(label = label), position = position_dodge(width = 0.75),
            vjust = -0.3, size = PS$text_mm, family = FONT_FAMILY) +
  scale_fill_manual(values = c('Phase 1' = '#1B9E77', 'Phase 2' = '#D95F02'),
                    name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  facet_wrap(~ metric, scales = 'free_y') +
  labs(x = NULL, y = NULL,
       caption = paste('% = (unique significant gene symbols, deduplicated across lineages) /',
                       '(sum of gene-lineage tests with padj != NA and baseMean > 50 across all 12 lineages).',
                       'Denominator grows with #lineages tested, so it is NOT the 477-gene Xenium panel size.',
                       sep = '\n')) +
  theme_v52(COMP_W) +
  theme(legend.position = 'top',
        legend.key.size = PS$legend_key,
        plot.caption = element_text(size = 6, hjust = 0))

save_figure(p_4D, 'Figure_4_panel_D.pdf', width = 8, height = 3.5)

###############################################################################
## (Former Panels E + F — TGFBR3 forest plot and FKBP5 violin) moved to the
## bottom of this file and commented out. Current Figure 4 renumbering:
##   4A = schematic          (was 4A)
##   4B = Phase1/Phase2 gene heatmaps (was 4B + 4C, combined label "B")
##   4C = WGCNA modules      (was 4G)
##   4D = DEG burden         (was 4D)
##   4E = CellNEST LR bars   (was 4H-J)
##   4F = Chord diagrams     (was 4K)
###############################################################################

###############################################################################
## Panel C — Bulk WGCNA module × phase × platform × LINEAGE (was Panel G)
## X axis  = WGCNA module numbers (1..N).
## Y axis  = platform × phase (Bulk_Ph1, Bulk_Ph2, snRNA_Ph1, ..., Xenium_Ph2).
## Facet columns = lineage (Pan / CM / Myeloid / EC / FB / Mural). Pan = all cells.
## snRNA & Xenium scores computed per-lineage subset (AddModuleScore).
## Bulk scores are lineage-agnostic and replicated across facets.
###############################################################################
consensus_modules <- read.csv('./dependencies/shared/bulk_heart_modules.csv')

.gene_col <- intersect(c('gene_name', 'gene', 'symbol'), colnames(consensus_modules))[1]
.mod_col  <- intersect(c('module', 'Module'), colnames(consensus_modules))[1]
mod_gene_lists <- split(consensus_modules[[.gene_col]],
                        consensus_modules[[.mod_col]])
mod_names   <- names(mod_gene_lists)
## M-number follows the paper's labels2colors(1:N) order (Figure_1.R L481, L624),
## so M-numbers are consistent between Figure 1 Panel E and Figure 4 Panel G.
mod_numbers <- setNames(
  as.character(match(mod_names, WGCNA::labels2colors(1:100))),
  mod_names
)

panel_g_lineages <- c('Pan', 'CM', 'Myeloid', 'EC', 'FB', 'PC')
panel_g_labels   <- c('Pan', 'CM', 'Myeloid', 'EC', 'FB', 'Mural')
.sn_lineage_col  <- intersect(c('Names','lineage','Lineage','cell_type'),
                               colnames(sn@meta.data))[1]

## --- Module activities + per-SAMPLE pseudobulk Welch's t-test per lineage ---
## Cells → per-sample means → Welch's t-test between disease groups.
## This respects sample as the biological replicate (avoids per-cell
## pseudoreplication inflation). Nominal p used; no BH correction because
## n=3–4 samples per group makes FDR too punitive.
## Returns list(means = df[module×group×lineage], pvalues = df[module×lineage×(p_ph1, p_ph2)])
.score_and_test_by_lineage <- function(obj, mod_lists, lineages, lineage_col,
                                         sample_col = NULL) {
  if (is.null(sample_col)) {
    sample_col <- intersect(c('patient','orig.ident','Sample','sample','donor'),
                             colnames(obj@meta.data))[1]
  }
  if (is.null(sample_col) || is.na(sample_col))
    stop('could not identify per-sample column for pseudobulk testing')
  if (length(mod_lists) == 0)
    return(list(means = data.frame(), pvalues = data.frame()))

  ## Score ONCE on the full object (avoids repeated subset + AddModuleScore,
  ## which OOMs on large Xenium objects).
  message('  scoring all modules on full object (', ncol(obj), ' cells)...')
  obj <- AddModuleScore(obj, features = mod_lists, name = 'WGCNAmod_')
  score_cols <- paste0('WGCNAmod_', seq_along(mod_lists))
  md_full <- obj@meta.data[, c('group', sample_col, score_cols), drop = FALSE]
  colnames(md_full) <- c('group', 'sample', names(mod_lists))
  md_full$.lineage <- if (!is.null(lineage_col) && !is.na(lineage_col))
    as.character(obj@meta.data[[lineage_col]]) else NA_character_
  ## Free the Seurat object as early as possible
  rm(obj); gc(verbose = FALSE)

  means_out <- list()
  pval_out  <- list()
  for (lin in lineages) {
    message('  aggregating lineage: ', lin)
    if (lin == 'Pan') {
      mdf <- md_full
    } else {
      if (is.null(lineage_col) || is.na(lineage_col)) {
        warning('lineage column not found; only Pan will be scored'); break
      }
      mdf <- md_full[md_full$.lineage == lin & !is.na(md_full$.lineage), , drop = FALSE]
      if (nrow(mdf) < 50) {
        message('    too few cells (', nrow(mdf), '); skipping'); next
      }
    }
    mdf$.lineage <- NULL
    long <- mdf %>% tidyr::pivot_longer(cols = -c(group, sample),
                                         names_to  = 'module',
                                         values_to = 'score')
    mns <- long %>% dplyr::group_by(module, group) %>%
      dplyr::summarise(mean_score = mean(score, na.rm = TRUE), .groups = 'drop')
    mns$lineage <- lin
    means_out[[lin]] <- mns
    ## pseudobulk: per-sample module-score means → Welch's t-test across groups
    samp_means <- long %>% dplyr::group_by(module, group, sample) %>%
      dplyr::summarise(sample_mean = mean(score, na.rm = TRUE), .groups = 'drop')
    pv <- samp_means %>% dplyr::group_by(module) %>%
      dplyr::summarise(
        p_ph1 = tryCatch(
          t.test(sample_mean[group == 'pRV'],
                 sample_mean[group == 'NF'])$p.value,
          error = function(e) NA_real_,
          warning = function(w) NA_real_),
        p_ph2 = tryCatch(
          t.test(sample_mean[group == 'RVF'],
                 sample_mean[group == 'pRV'])$p.value,
          error = function(e) NA_real_,
          warning = function(w) NA_real_),
        .groups = 'drop')
    pv$lineage <- lin
    pval_out[[lin]] <- pv
  }
  list(means   = dplyr::bind_rows(means_out),
       pvalues = dplyr::bind_rows(pval_out))
}

.cache_sn_mod_v2  <- './output/Figure_4/fig4_sn_mod_lineage_v2.rds'
.cache_xen_mod_v2 <- './output/Figure_4/fig4_xen_mod_lineage_v2.rds'
xen_mod_file <- xen_obj_file  # already defined above for pan pseudobulk

if (file.exists(.cache_sn_mod_v2)) {
  message('Loading cached snRNA per-lineage module scores + p-values...')
  sn_v2 <- readRDS(.cache_sn_mod_v2)
} else {
  message('Scoring snRNA modules per lineage (means + Wilcoxon p-values)...')
  sn_v2 <- .score_and_test_by_lineage(sn, mod_gene_lists, panel_g_lineages,
                                       .sn_lineage_col)
  saveRDS(sn_v2, .cache_sn_mod_v2)
}
sn_mod_lin <- sn_v2$means
sn_pval    <- sn_v2$pvalues

if (file.exists(.cache_xen_mod_v2)) {
  message('Loading cached Xenium per-lineage module scores + p-values...')
  xen_v2 <- readRDS(.cache_xen_mod_v2)
} else if (file.exists(xen_mod_file)) {
  message('Loading Xenium object for per-lineage module scoring + p-values...')
  xen <- readRDS(xen_mod_file)
  .xen_lineage_col <- intersect(c('cell_type_rctd_doublet','Names','names','lineage',
                                   'Lineage','cell_type'),
                                 colnames(xen@meta.data))[1]
  panel_mod_lists <- lapply(mod_gene_lists, function(g) intersect(g, rownames(xen)))
  panel_mod_lists <- panel_mod_lists[sapply(panel_mod_lists, length) >= 3]
  xen_v2 <- .score_and_test_by_lineage(xen, panel_mod_lists, panel_g_lineages,
                                        .xen_lineage_col)
  saveRDS(xen_v2, .cache_xen_mod_v2)
  rm(xen); gc()
} else {
  message('Xenium obj not found; Xenium slots will be empty.')
  xen_v2 <- list(
    means   = data.frame(module = character(), group = character(),
                         mean_score = numeric(), lineage = character()),
    pvalues = data.frame(module = character(), lineage = character(),
                         p_ph1 = numeric(), p_ph2 = numeric()))
}
xen_mod_lin <- xen_v2$means
xen_pval    <- xen_v2$pvalues

## --- Within-platform z-score per module across {NF, pRV, RVF} ---------------
## Rationale: bulk emits log2FC (units = log2 ratio); sn/Xenium emit AddModuleScore
## means (arbitrary units). Z-scoring per module across the 3 disease groups puts
## all platforms on the same unitless scale. Ph1 = z(pRV)-z(NF); Ph2 = z(RVF)-z(pRV).
## Direction is cross-platform comparable; magnitudes are only within-platform.
.z3 <- function(v) {
  s <- sd(v, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(NA_real_, length(v)))
  (v - mean(v, na.rm = TRUE)) / s
}

## --- Bulk: derive implied per-group mean from phase log2FC vectors ----------
## Bulk files are NF_vs_pRV and pRV_vs_RVF: DESeq2 numerator is the EARLIER stage,
## so positive log2FC = up in the earlier stage. To get the "up in later stage"
## Phase convention (pRV vs NF, RVF vs pRV), negate both vectors.
## NF is reference (0); pRV = -mean(log2FC_NF_vs_pRV); RVF = pRV - mean(log2FC_pRV_vs_RVF).
.bulk_module_delta <- function(mod_lists, bulk_ph1_df, bulk_ph2_df) {
  do.call(rbind, lapply(names(mod_lists), function(m) {
    g   <- mod_lists[[m]]
    v1  <- bulk_ph1_df$log2FoldChange[bulk_ph1_df$gene %in% g]
    v2  <- bulk_ph2_df$log2FoldChange[bulk_ph2_df$gene %in% g]
    v1  <- v1[!is.na(v1)]; v2 <- v2[!is.na(v2)]
    if (length(v1) == 0 || length(v2) == 0)
      return(data.frame(module = m, ph1 = NA_real_, ph2 = NA_real_))
    nf_m  <- 0
    prv_m <- -mean(v1)           # flip NF_vs_pRV → pRV delta from NF
    rvf_m <- prv_m - mean(v2)    # flip pRV_vs_RVF → RVF delta from pRV
    z     <- .z3(c(nf_m, prv_m, rvf_m))
    data.frame(module = m, ph1 = z[2] - z[1], ph2 = z[3] - z[2])
  }))
}
bulk_mod <- .bulk_module_delta(mod_gene_lists, bulk_ph1, bulk_ph2)

## --- Bulk module significance: Wilcoxon of in-module vs out-of-module log2FC --
.bulk_module_pvalue <- function(mod_lists, bulk_ph1_df, bulk_ph2_df) {
  do.call(rbind, lapply(names(mod_lists), function(m) {
    g <- mod_lists[[m]]
    in1  <- bulk_ph1_df$log2FoldChange[bulk_ph1_df$gene %in% g]
    out1 <- bulk_ph1_df$log2FoldChange[!(bulk_ph1_df$gene %in% g)]
    in2  <- bulk_ph2_df$log2FoldChange[bulk_ph2_df$gene %in% g]
    out2 <- bulk_ph2_df$log2FoldChange[!(bulk_ph2_df$gene %in% g)]
    in1 <- in1[!is.na(in1)]; out1 <- out1[!is.na(out1)]
    in2 <- in2[!is.na(in2)]; out2 <- out2[!is.na(out2)]
    p1 <- if (length(in1) >= 3 && length(out1) >= 3)
      suppressWarnings(wilcox.test(in1, out1)$p.value) else NA_real_
    p2 <- if (length(in2) >= 3 && length(out2) >= 3)
      suppressWarnings(wilcox.test(in2, out2)$p.value) else NA_real_
    data.frame(module = m, p_ph1 = p1, p_ph2 = p2)
  }))
}
bulk_pval <- .bulk_module_pvalue(mod_gene_lists, bulk_ph1, bulk_ph2)

## --- sn/Xenium: z-score mean_score across {NF, pRV, RVF} per module×lineage --
.phase_delta_lin <- function(df) {
  if (nrow(df) == 0) return(df)
  w <- tidyr::pivot_wider(df, id_cols = c(module, lineage),
                          names_from = group, values_from = mean_score)
  for (col in c('NF','pRV','RVF'))
    if (!col %in% colnames(w)) w[[col]] <- NA_real_
  z <- t(apply(as.matrix(w[, c('NF','pRV','RVF')]), 1, .z3))
  data.frame(module  = w$module,
             lineage = w$lineage,
             ph1 = z[, 2] - z[, 1],
             ph2 = z[, 3] - z[, 2])
}
sn_delta_lin  <- .phase_delta_lin(sn_mod_lin)
xen_delta_lin <- .phase_delta_lin(xen_mod_lin)

## --- Assemble long-form: module × lineage × platform × phase -------------
.build_row_g <- function(m, lin, plat, ph1, ph2, p1 = NA_real_, p2 = NA_real_) {
  data.frame(
    module   = rep(m, 2),
    lineage  = rep(lin, 2),
    platform_phase = c(paste0(plat, '_Ph1'), paste0(plat, '_Ph2')),
    delta  = c(ph1, ph2),
    pvalue = c(p1,  p2))
}
mod_long <- do.call(rbind, lapply(panel_g_lineages, function(lin) {
  do.call(rbind, lapply(mod_names, function(m) {
    bi <- match(m, bulk_mod$module)
    b1 <- bulk_mod$ph1[bi]; b2 <- bulk_mod$ph2[bi]
    bp <- bulk_pval[bulk_pval$module == m, ]
    bp1 <- if (nrow(bp)) bp$p_ph1[1] else NA_real_
    bp2 <- if (nrow(bp)) bp$p_ph2[1] else NA_real_
    s  <- sn_delta_lin[sn_delta_lin$module == m & sn_delta_lin$lineage == lin, ]
    s1 <- if (nrow(s)) s$ph1[1] else NA_real_
    s2 <- if (nrow(s)) s$ph2[1] else NA_real_
    sp <- sn_pval[sn_pval$module == m & sn_pval$lineage == lin, ]
    sp1 <- if (nrow(sp)) sp$p_ph1[1] else NA_real_
    sp2 <- if (nrow(sp)) sp$p_ph2[1] else NA_real_
    x  <- xen_delta_lin[xen_delta_lin$module == m & xen_delta_lin$lineage == lin, ]
    x1 <- if (nrow(x)) x$ph1[1] else NA_real_
    x2 <- if (nrow(x)) x$ph2[1] else NA_real_
    xp <- xen_pval[xen_pval$module == m & xen_pval$lineage == lin, ]
    xp1 <- if (nrow(xp)) xp$p_ph1[1] else NA_real_
    xp2 <- if (nrow(xp)) xp$p_ph2[1] else NA_real_
    rbind(
      .build_row_g(m, lin, 'Bulk',   b1, b2, bp1, bp2),
      .build_row_g(m, lin, 'snRNA',  s1, s2, sp1, sp2),
      .build_row_g(m, lin, 'Xenium', x1, x2, xp1, xp2)
    )
  }))
}))

mod_long$delta_capped <- pmin(pmax(mod_long$delta, -1.5), 1.5)
mod_long$module_num   <- factor(mod_numbers[mod_long$module],
                                levels = as.character(seq_along(mod_names)))
mod_long$lineage      <- factor(mod_long$lineage, levels = panel_g_lineages)
mod_long$platform_phase <- factor(mod_long$platform_phase,
  levels = rev(c('Bulk_Ph1','Bulk_Ph2','snRNA_Ph1','snRNA_Ph2','Xenium_Ph1','Xenium_Ph2')))
mod_long$platform <- sub('_Ph[12]$', '', as.character(mod_long$platform_phase))
mod_long$phase    <- sub('^.*_',      '', as.character(mod_long$platform_phase))

## --- Significance at nominal p < 0.05 (no BH correction) ---------------------
## Rationale: only n=3–4 samples per group, so BH on the module × lineage × phase
## grid is too punitive and drops biologically-supported signal. The test itself
## is already conservative (per-sample pseudobulk Welch's t-test for sn/Xenium;
## Wilcoxon in-module vs out-of-module log2FC for bulk).
mod_long$significant <- !is.na(mod_long$pvalue) & mod_long$pvalue < 0.05

## --- Cross-platform concordance per (module, lineage, phase): all 3 platforms
##     non-NA and same sign of delta.
concord_tbl <- mod_long %>%
  dplyr::group_by(module_num, lineage, phase) %>%
  dplyr::summarise(
    n_valid   = sum(!is.na(delta)),
    n_up      = sum(delta > 0, na.rm = TRUE),
    n_dn      = sum(delta < 0, na.rm = TRUE),
    concordant = n_valid == 3 & (n_up == 3 | n_dn == 3),
    .groups   = 'drop') %>%
  dplyr::select(module_num, lineage, phase, concordant)
mod_long <- dplyr::left_join(mod_long, concord_tbl,
                              by = c('module_num','lineage','phase'))
mod_long$concordant[is.na(mod_long$concordant)] <- FALSE

## Export panel G source data: module × lineage × platform × phase z-score deltas
.panel_g_out <- data.frame(
  module         = mod_long$module,
  module_num     = as.character(mod_long$module_num),
  lineage        = as.character(mod_long$lineage),
  platform_phase = as.character(mod_long$platform_phase),
  platform       = mod_long$platform,
  phase          = mod_long$phase,
  delta_zscore   = mod_long$delta,
  delta_capped   = mod_long$delta_capped,
  pvalue         = mod_long$pvalue,
  significant    = mod_long$significant,
  concordant     = mod_long$concordant,
  stringsAsFactors = FALSE
)
.panel_c_csv <- './output/Figure_4/Figure_4_panel_C_data.csv'
dir.create(dirname(.panel_c_csv), showWarnings = FALSE, recursive = TRUE)
write.csv(.panel_g_out, .panel_c_csv, row.names = FALSE)
message('Panel C data written: ', .panel_c_csv)

.tile_border_lw <- 0.1   # thinner inter-cell grid

## Manual 45°-in-display stripe segments for non-significant tiles.
## ggpattern's geom_tile_pattern consistently bleeds past tile edges regardless of
## parameters, so we draw segments geometrically clipped to [-0.5, 0.5]^2.
##
## Tiles render taller than they are wide (tall narrow columns, multiple rows),
## so to appear at 45° in display the segment slope in DATA coords must equal
## the tile width:height ratio in display. .stripe_slope is tunable; 0.30 matches
## the current composite layout. Lower = shallower (closer to horizontal) if tiles
## get taller; higher = steeper if tiles become wider.
.stripe_slope <- 0.30
.b_max        <- abs(.stripe_slope) * 0.5 + 0.5
.stripe_offsets <- seq(-0.85 * .b_max, 0.85 * .b_max, length.out = 9)

.clip_stripe <- function(slope, b) {
  pts <- list()
  y_L <- slope * (-0.5) + b
  if (y_L >= -0.5 && y_L <=  0.5) pts[[length(pts)+1]] <- c(-0.5, y_L)
  y_R <- slope *   0.5  + b
  if (y_R >= -0.5 && y_R <=  0.5) pts[[length(pts)+1]] <- c( 0.5, y_R)
  if (slope != 0) {
    x_B <- (-0.5 - b) / slope
    if (x_B > -0.5 && x_B <  0.5) pts[[length(pts)+1]] <- c(x_B, -0.5)
    x_T <- ( 0.5 - b) / slope
    if (x_T > -0.5 && x_T <  0.5) pts[[length(pts)+1]] <- c(x_T,  0.5)
  }
  if (length(pts) < 2) return(NULL)
  m <- do.call(rbind, pts)
  m <- m[!duplicated(round(m, 8)), , drop = FALSE]
  if (nrow(m) < 2) return(NULL)
  c(x = m[1, 1], y = m[1, 2], xend = m[2, 1], yend = m[2, 2])
}

.non_sig <- mod_long[!mod_long$significant & !is.na(mod_long$delta), ]
.pp_levels <- levels(factor(mod_long$platform_phase))
.non_sig$y_num <- match(as.character(.non_sig$platform_phase), .pp_levels)
.non_sig$x_num <- NA_real_
for (.lin in unique(.non_sig$lineage)) {
  .sub_lvls <- levels(factor(mod_long$module_num[mod_long$lineage == .lin]))
  .idx <- .non_sig$lineage == .lin
  .non_sig$x_num[.idx] <- match(as.character(.non_sig$module_num[.idx]), .sub_lvls)
}

.stripe_segs <- do.call(rbind, lapply(seq_len(nrow(.non_sig)), function(i) {
  xc <- .non_sig$x_num[i]; yc <- .non_sig$y_num[i]
  cell_segs <- do.call(rbind, lapply(.stripe_offsets, function(b) {
    s <- .clip_stripe(.stripe_slope, b)
    if (is.null(s)) return(NULL)
    data.frame(x = xc + s['x'], y = yc + s['y'],
               xend = xc + s['xend'], yend = yc + s['yend'])
  }))
  if (is.null(cell_segs) || nrow(cell_segs) == 0) return(NULL)
  cell_segs$lineage <- .non_sig$lineage[i]
  cell_segs
}))

p_4C <- ggplot(mod_long, aes(x = module_num, y = platform_phase,
                              fill = delta_capped)) +
  ## Base fill for all cells (including NA via na.value = 'grey90')
  geom_tile(linewidth = 0) +
  ## Clipped diagonal stripes drawn BEFORE grid lines (cannot exit tile bounds)
  geom_segment(data = .stripe_segs,
               aes(x = x, xend = xend, y = y, yend = yend),
               colour = '#444444', linewidth = 0.15, inherit.aes = FALSE) +
  ## Inter-cell grid lines drawn ON TOP so striations sit behind them
  geom_tile(fill = NA, colour = 'grey90', linewidth = .tile_border_lw) +
  ## Red border for cross-platform concordant cells (sign agreement across 3 platforms)
  #geom_tile(
  #  data = mod_long[mod_long$concordant & !is.na(mod_long$delta), ],
  #  fill = NA, colour = '#D7191C', linewidth = 0.45) +
  scale_fill_gradientn(
    colours  = c('#2166AC','white','#B2182B'),
    limits   = c(-2, 2),
    na.value = 'grey90',
    name     = '\u0394 score',
    guide    = guide_colorbar(direction = 'vertical',
      barwidth  = unit(0.3 * PS$font_scale, 'cm'),
      barheight = unit(2.0 * PS$font_scale, 'cm'),
      title.position = 'top', ticks.colour = 'black')) +
  facet_grid(. ~ lineage, scales = 'free_x', space = 'free_x',
             labeller = labeller(lineage = setNames(panel_g_labels, panel_g_lineages))) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  labs(x = 'WGCNA module', y = NULL,
       title = 'Bulk WGCNA modules \u00d7 platform \u00d7 phase \u00d7 lineage (hatched = nominal p \u2265 0.05, pseudobulk Welch t-test; red border = 3-platform sign concordance)') +
  theme_v52(COMP_W) +
  theme(axis.text.x      = element_text(size = PS$base_pt),
        axis.text.y      = element_text(size = PS$base_pt),
        strip.text.x     = element_text(size = PS$base_pt, face = 'bold'),
        strip.background = element_rect(fill = 'grey95', colour = NA),
        legend.position  = 'right',
        panel.border     = element_rect(colour='black', fill=NA, linewidth = PS$geom_lw),
        panel.spacing.x  = unit(0.4, 'lines'),
        axis.ticks       = element_blank())

save_figure(p_4C, 'Figure_4_panel_C.pdf', width = 14, height = 3.5)

###############################################################################
## Panels E & F — CellNEST analysis RETIRED (2026-07): CellNEST is no longer
## used. Fig 4E is now the CellChat top-15 pathways (Myocardium niche) and 4F
## the CellChat normalized myocardium chord, produced by additional_scripts/
## spatial_cellchat_niche_communication.R (fig1 / fig2) and assembled in
## Illustrator. The original CellNEST code is disabled with `if (FALSE)`; the
## p_4E / p_4F placeholders below keep the composite buildable.
###############################################################################
## Panel E — CellNEST top LR pairs per disease state (was Panels H-J)
## Panel F — Chord diagrams: cell type communication by disease state (was Panel K)
###############################################################################
if (FALSE) {  ## <<< retired CellNEST E/F — begin (disabled, kept for reference)
if (!exists('xen_meta')) {
  xen_meta      <- read.csv('./dependencies/shared/Xenium/Xenium_metadata.csv', row.names = 1)
  patient_group <- xen_meta %>%
    dplyr::select(patient, group) %>%
    dplyr::distinct()
  patient_group$patient <- as.character(patient_group$patient)
}

cellnest_base  <- './dependencies/shared/Xenium/cellnest_output/output'
cellnest_pats  <- list.dirs(cellnest_base, recursive = FALSE, full.names = FALSE)
cellnest_pats  <- cellnest_pats[grepl('Xenium_resegmented', cellnest_pats)]

.extract_pat <- function(d) sub('.*_(\\d+)$', '\\1', d)

## --- validate CellNEST inputs per patient (the FILES, not just the dir) ---
## A patient feeds 4E via its *_histogram_byFrequency_table_top1500.csv and
## 4F via its *_allCCC.csv. An output dir can exist while these files are
## absent (a failed/partial run leaves an empty dir), so check the files and
## restrict the ensemble + n report to patients with usable inputs.
.dir_for  <- function(pat) paste0('Xenium_resegmented_imputed_final_', pat)
.freq_for <- function(d) file.path(cellnest_base, d,
  paste0('CellNEST_', d, '_histogram_byFrequency_table_top1500.csv'))
.ccc_for  <- function(d) file.path(cellnest_base, d,
  paste0('CellNEST_', d, '_allCCC.csv'))

.cellnest_expected <- patient_group$patient[
  patient_group$group %in% c('NF', 'pRV', 'RVF')]
.usable <- character(0)
for (pat in .cellnest_expected) {
  d  <- .dir_for(pat)
  gp <- patient_group$group[patient_group$patient == pat][1]
  has_dir  <- dir.exists(file.path(cellnest_base, d))
  has_freq <- file.exists(.freq_for(d)) && file.size(.freq_for(d)) > 0
  has_ccc  <- file.exists(.ccc_for(d))  && file.size(.ccc_for(d))  > 0
  if (has_freq && has_ccc) { .usable <- c(.usable, pat); next }
  reason <- if (!has_dir) 'no output dir'
            else if (!has_freq && !has_ccc) 'output dir empty (no allCCC / freq table)'
            else if (!has_ccc) 'missing allCCC.csv (4F input)'
            else 'missing freq table (4E input)'
  warning(sprintf(
    'CellNEST: patient %s (group %s) %s under %s -- DROPPED from Panels 4E/4F ensemble.',
    pat, gp, reason, cellnest_base), call. = FALSE, immediate. = TRUE)
}
## restrict downstream ensemble to patients with usable inputs
cellnest_pats <- cellnest_pats[.extract_pat(cellnest_pats) %in% .usable]
.found_grp <- patient_group$group[match(.usable, patient_group$patient)]
message(sprintf('CellNEST 4E/4F ensemble: %d/%d patients (%s).',
  length(.usable), length(.cellnest_expected),
  paste(names(table(.found_grp)), table(.found_grp),
        sep = '=', collapse = ', ')))

lr_all <- lapply(cellnest_pats, function(d) {
  pat   <- .extract_pat(d)
  grp   <- patient_group$group[patient_group$patient == pat]
  if (length(grp) == 0) return(NULL)
  freq_file <- file.path(cellnest_base, d,
    paste0('CellNEST_', d, '_histogram_byFrequency_table_top1500.csv'))
  if (!file.exists(freq_file)) return(NULL)
  df <- read.csv(freq_file)
  colnames(df) <- c('lr_pair', 'count')
  df$patient <- pat; df$group <- grp
  df
})
lr_all <- dplyr::bind_rows(Filter(Negate(is.null), lr_all))

## --- per-patient normalization (replaces pooled sum) ---
## Within each patient, express every LR pair as a fraction of THAT patient's
## total CellNEST edge count, then average across the group's patients. Patients
## that never list a pair contribute 0 for it (zero-fill within group), so the
## mean is over ALL patients in the group, not just those expressing the pair.
## This removes the section-size / patient-count bias of pooled summation.
lr_pat <- lr_all %>%
  dplyr::group_by(group, patient) %>%
  dplyr::mutate(frac = count / sum(count)) %>%
  dplyr::ungroup() %>%
  dplyr::select(group, patient, lr_pair, frac) %>%
  dplyr::group_by(group) %>%
  tidyr::complete(patient, lr_pair, fill = list(frac = 0)) %>%
  dplyr::ungroup()

lr_agg <- lr_pat %>%
  dplyr::group_by(group, lr_pair) %>%
  dplyr::summarise(
    n_pat     = dplyr::n(),
    mean_frac = mean(frac),
    se_frac   = if (dplyr::n() > 1) stats::sd(frac) / sqrt(dplyr::n())
                else NA_real_,
    .groups = 'drop') %>%
  dplyr::group_by(group) %>%
  dplyr::slice_max(mean_frac, n = 15) %>%
  dplyr::ungroup()

.make_lr_bar <- function(grp_label, fill_col, title) {
  df <- dplyr::filter(lr_agg, group == grp_label) %>%
    dplyr::arrange(mean_frac) %>%
    dplyr::mutate(lr_pair = factor(lr_pair, levels = lr_pair))
  if (nrow(df) == 0)
    return(ggplot() + annotate('text', x=0.5, y=0.5,
      label=paste0('[', grp_label, ' LR — no data]'),
      size = PS$text_mm, family = FONT_FAMILY, colour='grey50', hjust=0.5, vjust=0.5) + theme_void())
  ggplot(df, aes(x = mean_frac * 100, y = lr_pair)) +
    geom_col(fill = fill_col, width = 0.7, linewidth = PS$geom_lw) +
    geom_errorbarh(aes(xmin = pmax((mean_frac - se_frac) * 100, 0),
                       xmax = (mean_frac + se_frac) * 100),
                   height = 0.3, linewidth = PS$geom_lw, na.rm = TRUE) +
    labs(x = 'Mean per-patient % of CellNEST edges', y = NULL, title = title) +
    theme_v52(COMP_W) +
    theme(panel.grid.major.y = element_blank())
}

p_4H <- .make_lr_bar('NF',  disease_pal['NF'],  'NF — Top CellNEST LR pairs')
p_4I <- .make_lr_bar('pRV', disease_pal['pRV'], 'pRV — Top CellNEST LR pairs')
p_4J <- .make_lr_bar('RVF', disease_pal['RVF'], 'RVF — Top CellNEST LR pairs')

p_4E <- p_4H | p_4I | p_4J  # renumbered Panel E = old Panels H-J combined
save_figure(p_4E, 'Figure_4_panel_E.pdf', width = 14, height = 5)

cellnest_meta_base <- './dependencies/shared/Xenium/cellnest_output/metadata'

## Cell-type lookup keyed by Proseg cell barcode "<cellid>_<patient_idx>"
## (see /Users/ikuz/Documents/XeniumWorkflow/xenium_explorer.py:_load_cellnest_data —
##  CCC from_id/to_id are 0-based row indices into cell_barcode_*.csv; resulting
##  "<id>_<patient_idx>" strings match colnames of the resegmented Xenium RDS.)
.cache_xen_ct_lookup <- './output/Figure_4/fig4_xen_ct_lookup.rds'
xen_reseg_rds <- './dependencies/shared/Xenium/Xenium_resegmented_imputed_final.rds'

if (file.exists(.cache_xen_ct_lookup)) {
  message('Loading cached Xenium reseg cell-type lookup...')
  ct_lookup <- readRDS(.cache_xen_ct_lookup)
} else {
  message('Building Xenium reseg cell-type lookup from combined RDS...')
  .xen_reseg <- readRDS(xen_reseg_rds)
  ct_lookup <- as.character(.xen_reseg$cell_type_rctd_doublet)
  names(ct_lookup) <- colnames(.xen_reseg)
  saveRDS(ct_lookup, .cache_xen_ct_lookup)
  rm(.xen_reseg); gc(verbose = FALSE)
}

.ccc_per_pat <- lapply(cellnest_pats, function(d) {
  pat   <- .extract_pat(d)
  grp   <- patient_group$group[patient_group$patient == pat]
  if (length(grp) == 0) return(NULL)
  ccc_file <- file.path(cellnest_base, d,
    paste0('CellNEST_', d, '_allCCC.csv'))
  bar_file <- file.path(cellnest_meta_base, d,
    paste0('cell_barcode_', d, '.csv'))
  if (!file.exists(ccc_file) || !file.exists(bar_file)) return(NULL)
  ccc  <- read.csv(ccc_file)
  bars <- read.csv(bar_file, header = FALSE, col.names = c('barcode'))$barcode
  ## from_id/to_id are 0-based row indices into bars; the resulting barcodes
  ## (e.g. "22750_4") are direct keys into ct_lookup.
  ## NOTE: the regenerated allCCC.csv also carries explicit from_cell/to_cell
  ## barcode columns; verified that bars[from_id + 1L] == from_cell exactly
  ## (R 1-based index of the 0-based id), so this reconstruction is correct.
  ccc$from_barcode <- bars[ccc$from_id + 1L]
  ccc$to_barcode   <- bars[ccc$to_id   + 1L]
  ccc$from_type    <- ct_lookup[ccc$from_barcode]
  ccc$to_type      <- ct_lookup[ccc$to_barcode]
  ccc$group <- grp; ccc$patient <- pat
  ccc <- ccc[!is.na(ccc$from_type) & !is.na(ccc$to_type), ]
  edge_df <- ccc %>%
    dplyr::count(group, patient, from_type, to_type, name = 'n_edges')
  ## Per-patient cell counts per type (for abundance normalization)
  bar_types <- ct_lookup[bars]
  bar_types <- bar_types[!is.na(bar_types)]
  cell_df <- data.frame(group = grp, patient = pat,
                        type = unname(bar_types)) %>%
    dplyr::count(group, patient, type, name = 'n_cells')
  list(edges = edge_df, cells = cell_df)
})
.ccc_per_pat <- Filter(Negate(is.null), .ccc_per_pat)
ccc_all   <- dplyr::bind_rows(lapply(.ccc_per_pat, `[[`, 'edges'))
cells_all <- dplyr::bind_rows(lapply(.ccc_per_pat, `[[`, 'cells'))

## Cell types to display + log2(O/E) cap for cross-panel comparability
chord_keep_types <- c('CM', 'Myeloid', 'EC', 'FB', 'PC', 'SM', 'Endo')
chord_oe_cap     <- Inf  # uncapped — let extreme values render at full magnitude

.make_chord <- function(grp_label, title, min_obs = 10,
                         keep_types = chord_keep_types, oe_cap = chord_oe_cap) {
  df <- dplyr::filter(ccc_all, group == grp_label)
  if (nrow(df) == 0)
    return(ggplot() + annotate('text', x=0.5, y=0.5,
      label=paste0('[', grp_label, ' chord — no data]'),
      size = PS$text_mm, family = FONT_FAMILY, colour='grey50', hjust=0.5, vjust=0.5) + theme_void())

  ## --- per-patient O/E, then averaged across the group's patients ---
  ## Replaces pooled-sum normalization. Each patient contributes equally
  ## regardless of section size: abundance-aware log2(O/E) is computed WITHIN
  ## each patient (expected = that patient's total edges x p_A x p_B from that
  ## patient's own cell-type abundances), then the MEAN is taken across the
  ## group's patients. Pairs absent in a patient count as 0 edges for that
  ## patient (zero-filled over the displayed type x type grid), so a pair seen
  ## in only some patients is correctly down-weighted. Self-loops dropped
  ## (within-type spatial autocorrelation is structural niche confinement,
  ## not signaling). A pair is shown only if its MEAN per-patient edge count
  ## clears the min_obs floor; weight clipped to [0, oe_cap].
  kt <- keep_types[keep_types %in%
                   cells_all$type[cells_all$group == grp_label]]
  combos <- expand.grid(from_type = kt, to_type = kt,
                         stringsAsFactors = FALSE)
  combos <- combos[combos$from_type != combos$to_type, , drop = FALSE]

  per_pat <- lapply(unique(df$patient), function(pt) {
    npt <- setNames(rep(0, length(kt)), kt)
    cg  <- cells_all %>%
      dplyr::filter(group == grp_label, patient == pt, type %in% kt)
    npt[cg$type] <- cg$n_cells
    Ntot <- sum(npt)
    if (Ntot == 0) return(NULL)
    e <- dplyr::filter(df, patient == pt,
                       from_type %in% kt, to_type %in% kt,
                       from_type != to_type)
    m <- merge(combos, e[, c('from_type', 'to_type', 'n_edges')],
               by = c('from_type', 'to_type'), all.x = TRUE)
    m$n_edges[is.na(m$n_edges)] <- 0
    tot_e <- sum(m$n_edges)
    m$expected <- tot_e * (npt[m$from_type] / Ntot) * (npt[m$to_type] / Ntot)
    m$log2OE   <- log2((m$n_edges + 1) / (m$expected + 1))
    m$patient  <- pt
    m
  })
  per_pat <- dplyr::bind_rows(Filter(Negate(is.null), per_pat))
  if (nrow(per_pat) == 0)
    return(ggplot() + annotate('text', x=0.5, y=0.5,
      label=paste0('[', grp_label, ' chord — no data]'),
      size = PS$text_mm, family = FONT_FAMILY, colour='grey50') + theme_void())

  df <- per_pat %>%
    dplyr::group_by(from_type, to_type) %>%
    dplyr::summarise(mean_log2OE = mean(log2OE),
                     mean_edges  = mean(n_edges),
                     .groups = 'drop') %>%
    dplyr::mutate(
      weight = ifelse(mean_edges >= min_obs,
                      pmin(pmax(mean_log2OE, 0), oe_cap), 0)
    ) %>%
    dplyr::select(from_type, to_type, weight)

  mat_df <- df %>%
    tidyr::pivot_wider(names_from = to_type, values_from = weight, values_fill = 0) %>%
    tibble::column_to_rownames('from_type')
  mat <- as.matrix(mat_df)
  ## drop rows/cols that are entirely zero (no enriched edges)
  mat <- mat[rowSums(mat) > 0 | colSums(mat) > 0,
             rowSums(mat) > 0 | colSums(mat) > 0, drop = FALSE]
  if (nrow(mat) < 2)
    return(ggplot() + annotate('text', x=0.5, y=0.5,
      label=paste0('[', grp_label, ' chord — no enriched pairs]'),
      size = PS$text_mm, family = FONT_FAMILY, colour='grey50') + theme_void())

  all_types <- union(rownames(mat), colnames(mat))
  ## Distinct color per displayed Xenium cell type (PC≠SM, Endo≠EC, etc.).
  ## Base hues come from lineage_pal; near-related types get a perceptually
  ## distinct shifted hue.
  chord_pal <- c(
    CM      = unname(lineage_pal['CM']),       # red
    Myeloid = unname(lineage_pal['Myeloid']),  # purple
    EC      = unname(lineage_pal['EC']),       # green
    FB      = unname(lineage_pal['FB']),       # blue
    PC      = unname(lineage_pal['Mural']),    # brown
    SM      = '#D2691E',                       # chocolate (mural-adjacent, distinct from PC)
    Endo    = '#1B9E77',                       # teal (vascular-adjacent, distinct from EC)
    Adipo   = unname(lineage_pal['Adipocyte']),
    Neuron  = unname(lineage_pal['Neuronal']),
    Epi     = unname(lineage_pal['Epicardial']),
    LEC     = unname(lineage_pal['Lymphatic']),
    NKT     = unname(lineage_pal['NKT'])
  )
  grid_cols <- chord_pal[all_types]
  grid_cols[is.na(grid_cols)] <- 'grey70'
  names(grid_cols) <- all_types

  suppressPackageStartupMessages({
    library(circlize)
    library(ggplotify)
  })

  ggplotify::as.ggplot(function() {
    circlize::circos.clear()
    circlize::circos.par(gap.after = 2, start.degree = 90,
                         cell.padding = c(0, 0, 0, 0))
    circlize::chordDiagram(
      mat,
      grid.col = grid_cols,
      transparency = 0.25,
      directional = 1,
      direction.type = c('arrows', 'diffHeight'),
      diffHeight = -0.04,
      annotationTrack = 'grid',
      preAllocateTracks = list(track.height = 0.04),
      link.arr.type = 'big.arrow',
      link.arr.length = 0.15
    )
    circlize::circos.trackPlotRegion(
      track.index = 1,
      panel.fun = function(x, y) {
        xlim <- circlize::CELL_META$xlim
        sector.name <- circlize::CELL_META$sector.index
        circlize::circos.text(mean(xlim), 0.4, sector.name,
                              facing = 'clockwise', niceFacing = TRUE,
                              adj = c(0, 0.5), cex = 0.65)
      },
      bg.border = NA)
    circlize::circos.clear()
    graphics::title(title, cex.main = 0.9, line = -0.5)
  })
}

p_4K_nf  <- .make_chord('NF',  'NF: mean per-patient log2(O/E)')
p_4K_prv <- .make_chord('pRV', 'pRV: mean per-patient log2(O/E)')
p_4K_rvf <- .make_chord('RVF', 'RVF: mean per-patient log2(O/E)')
p_4K     <- p_4K_nf | p_4K_prv | p_4K_rvf

p_4F <- p_4K  # renumbered Panel F = old Panel K (chord diagrams)
save_figure(p_4F, 'Figure_4_panel_F.pdf', width = 14, height = 5)
}  ## <<< retired CellNEST E/F — end

###############################################################################
## Panels E & F — CellChat spatial niche communication (Myocardium niche)
##
## Full generating pipeline (per-patient x per-niche CellChat over the resegmented
## Xenium object): additional_scripts/spatial_cellchat_niche_communication.R. That
## pipeline is HEAVY — per-patient x niche CellChat on the 9GB BPCells object;
## needs CellChat (GitHub jinworks/CellChat) + BPCells; ~hours — and writes the
## aggregated table shipped as Supplementary Data S14:
##   dependencies/shared/Xenium/xenium_cellchat_niche_communication.csv
##
## Panels E and F are rebuilt natively below from that shipped table, so Figure 4
## renders WITHOUT CellChat/BPCells at figure time:
##   E (fig1) = top-15 signaling pathways in Myocardium, normalized to sum=1 per
##       group. Reconstructed EXACTLY: the pipeline's per-patient pooled probability
##       equals mean_prob x n_patients per L-R (mean over present patients x n =
##       the summed per-patient probability), so summing that per pathway reproduces
##       make_top_pathway_bar()'s pooled_prob.
##   F (fig2) = cell-type communication chord in Myocardium, normalized by
##       sender x receiver cell-type proportion. The source->target weight matrix is
##       the per-L-R probability summed per cell-type pair (= CellChat net$weight);
##       proportions come from the cached metadata file
##       dependencies/shared/Xenium/xenium_myocardium_celltype_proportions.csv
##       (cell_type_grouped counts within Myocardium per group — the same idents
##       CellChat used; built by helper_scripts/make_myocardium_proportions.R).
##       Link colours use a fixed palette (the pipeline's CellChat::scPalette is
##       unavailable here) — cosmetic only.
##
## Set RUN_CELLCHAT_FROM_RAW <- TRUE to regenerate the underlying table (and the
## pipeline's own object-exact PDF panels) from raw before building E/F.
###############################################################################
RUN_CELLCHAT_FROM_RAW <- FALSE

.cc_group_order  <- c('NF', 'pRV', 'RVF')
.cc_group_colors <- c(NF = '#4393C3', pRV = '#F4A582', RVF = '#D6604D')
.cc_target_niche <- 'Myocardium'
.cc_s14_csv  <- './dependencies/shared/Xenium/xenium_cellchat_niche_communication.csv'
.cc_prop_csv <- './dependencies/shared/Xenium/xenium_myocardium_celltype_proportions.csv'

if (RUN_CELLCHAT_FROM_RAW) {
  message('RUN_CELLCHAT_FROM_RAW=TRUE — regenerating S14 table from raw ',
          '(heavy: per-patient x niche CellChat; needs CellChat + BPCells)...')
  source('./additional_scripts/spatial_cellchat_niche_communication.R')
  ## The pipeline writes output/spatial_cellchat/.../supplementary_table_niche_communication.csv;
  ## point .cc_s14_csv there to build E/F from the freshly regenerated table.
}

## Fallback so the composite still assembles if inputs/packages are unavailable.
.panel_placeholder <- function(txt)
  ggplot2::ggplot() +
  ggplot2::annotate('text', x = 0.5, y = 0.5, label = txt, size = 4) +
  ggplot2::theme_void()

## ---- Panel E: top-15 pathways, Myocardium (native, exact) ----
p_4E <- tryCatch({
  stopifnot(file.exists(.cc_s14_csv))
  cc  <- utils::read.csv(.cc_s14_csv, stringsAsFactors = FALSE)
  myo <- cc[cc$niche == .cc_target_niche, ]
  ## per-patient pooled probability per L-R = mean_prob x n_patients (exact)
  long <- do.call(rbind, lapply(.cc_group_order, function(g)
    data.frame(group   = g,
               pathway = myo$pathway_name,
               pooled  = myo[[paste0('mean_prob_', g)]] * myo[[paste0('n_patients_', g)]],
               stringsAsFactors = FALSE)))
  ## recode SEMA3*/SEMA6* -> SEMA (matches pipeline recode_sema_lr)
  long$pathway <- ifelse(grepl('^SEMA3|^SEMA6', long$pathway), 'SEMA', long$pathway)
  long <- long[!is.na(long$pathway) & long$pathway != '', ]
  summ <- aggregate(pooled ~ group + pathway, long, sum)
  summ$pooled <- ave(summ$pooled, summ$group, FUN = function(v) v / sum(v))
  x_max <- max(summ$pooled) * 1.08
  plots <- lapply(.cc_group_order, function(g) {
    sg  <- summ[summ$group == g, ]
    sub <- utils::head(sg[order(-sg$pooled), ], 15)
    sub$pathway <- factor(sub$pathway, levels = rev(unique(sub$pathway)))
    ycol <- ifelse(levels(sub$pathway) == 'SEMA', '#D6604D', 'black')
    yfc  <- ifelse(levels(sub$pathway) == 'SEMA', 'bold', 'plain')
    ggplot(sub, aes(pathway, pooled, fill = pooled)) +
      geom_col(width = 0.75, alpha = 0.9) +
      coord_flip() +
      scale_fill_gradient(low  = grDevices::adjustcolor(.cc_group_colors[g], 0.3),
                          high = .cc_group_colors[g], guide = 'none') +
      scale_y_continuous(limits = c(0, x_max), expand = expansion(mult = c(0, 0))) +
      labs(title = g, x = NULL, y = 'Proportion of total communication probability') +
      theme_bw(base_size = 11) +
      theme(plot.title  = element_text(face = 'bold', size = 13, color = .cc_group_colors[g]),
            axis.text.y = element_text(size = 8, color = ycol, face = yfc))
  })
  pE <- patchwork::wrap_plots(plots, nrow = 1)
  save_figure(pE, 'Figure_4_panel_E.pdf', width = 18, height = 6)
  pE
}, error = function(e) {
  message('Panel 4E fell back to placeholder: ', conditionMessage(e))
  .panel_placeholder('Panel 4E — CellChat top-15 pathways (Myocardium)')
})

## ---- Panel F: Myocardium chord, proportion-normalized (native) ----
p_4F <- tryCatch({
  stopifnot(file.exists(.cc_s14_csv), file.exists(.cc_prop_csv),
            requireNamespace('circlize',  quietly = TRUE),
            requireNamespace('ggplotify', quietly = TRUE))
  cc    <- utils::read.csv(.cc_s14_csv, stringsAsFactors = FALSE)
  myo   <- cc[cc$niche == .cc_target_niche, ]
  props <- utils::read.csv(.cc_prop_csv, stringsAsFactors = FALSE)
  cts   <- sort(unique(c(myo$source, myo$target)))
  ct_colors <- setNames(scales::hue_pal()(length(cts)), cts)  # scPalette unavailable

  build_norm_mat <- function(g) {
    Pg <- max(myo[[paste0('n_patients_', g)]], na.rm = TRUE)
    if (!is.finite(Pg) || Pg == 0) Pg <- 1
    ## avg net$weight over patients = sum over L-R of (mean_prob x n_patients) / Pg
    w   <- myo[[paste0('mean_prob_', g)]] * myo[[paste0('n_patients_', g)]] / Pg
    agg <- aggregate(w, list(source = myo$source, target = myo$target), sum)
    M   <- matrix(0, length(cts), length(cts), dimnames = list(cts, cts))
    for (k in seq_len(nrow(agg))) M[agg$source[k], agg$target[k]] <- agg$x[k]
    ## proportion normalization: divide each edge by prop_sender x prop_receiver
    pr <- setNames(props$prop[props$group == g], props$cell_type_grouped[props$group == g])
    p  <- pr[cts]; p[is.na(p)] <- 1e-6
    for (i in seq_along(cts)) for (j in seq_along(cts))
      M[i, j] <- M[i, j] / max(p[i] * p[j], 1e-10)
    M
  }
  mats <- setNames(lapply(.cc_group_order, build_norm_mat), .cc_group_order)

  draw_one <- function(M, g, min_frac = 0.02) {
    keep <- rownames(M)[rowSums(M) > 0 | colSums(M) > 0]
    if (length(keep) < 2) {
      plot.new(); title(main = sprintf('%s\n(no signal)', g), col.main = .cc_group_colors[g])
      return(invisible())
    }
    Ms <- M[keep, keep, drop = FALSE]
    Ms[Ms < max(Ms) * min_frac] <- 0
    if (sum(Ms) == 0) {
      plot.new(); title(main = sprintf('%s\n(below threshold)', g), col.main = .cc_group_colors[g])
      return(invisible())
    }
    gcol <- ct_colors[keep]
    colm <- matrix(NA_character_, length(keep), length(keep), dimnames = list(keep, keep))
    for (i in keep) colm[i, ] <- grDevices::adjustcolor(gcol[i], 0.6)  # colour links by sender
    circlize::circos.clear()
    circlize::chordDiagram(Ms, grid.col = gcol, col = colm, transparency = 0.3,
                           annotationTrack = c('grid', 'name'),
                           preAllocateTracks = list(track.height = 0.08),
                           directional = 1, direction.type = 'arrows',
                           link.arr.type = 'big.arrow', self.link = 1, scale = FALSE)
    circlize::circos.clear()
    title(main = g, col.main = .cc_group_colors[g], cex.main = 1.2, line = -1)
  }

  draw_all <- function() {
    op <- graphics::par(mfrow = c(1, 3), mar = c(2, 2, 3, 2)); on.exit(graphics::par(op))
    for (g in .cc_group_order) draw_one(mats[[g]], g)
  }
  pF <- ggplotify::as.ggplot(draw_all)
  save_figure(pF, 'Figure_4_panel_F.pdf', width = 18, height = 6)
  pF
}, error = function(e) {
  message('Panel 4F fell back to placeholder: ', conditionMessage(e))
  .panel_placeholder('Panel 4F — CellChat myocardium chord (normalized)')
})

###############################################################################
## Assemble Figure 4 (patchwork, 6 panels A-F after renumbering)
## A = schematic                                       (was A)
## B = Phase1/Phase2 key-gene heatmaps (both stacked)  (was B + C)
## C = WGCNA modules × platform × phase                (was G)
## D = DEG burden (phase × platform)                   (was D)
## E = CellChat top-15 pathways, Myocardium           (native; was CellNEST H-J)
## F = CellChat myocardium chord (normalized)         (native; was CellNEST K)
###############################################################################
row_A_B   <- p_4A | p_4BC       # Panel A + combined Panel B heatmaps
row_C     <- p_4C               # WGCNA modules
row_D     <- p_4D               # DEG burden
row_E     <- p_4E               # CellChat top-15 pathways (native, from S14 table)
row_F     <- p_4F               # CellChat chord (native, from S14 table + proportions)

fig4 <- row_A_B /
        row_C /
        row_D /
        row_E /
        row_F +
  plot_annotation(
    tag_levels = list(c('A','B','C','D','E','F')),
    theme = theme(plot.tag = element_text(
      family = FONT_FAMILY,
      size   = 16 * COMP_W / FINAL_WIDTH_IN,
      face   = 'bold'
    ))
  ) +
  plot_layout(heights = c(3, 4, 2.5, 2, 2))

save_figure(fig4, 'Figure_4.pdf', width = 14, height = 22)
message('Figure 4 complete.')

