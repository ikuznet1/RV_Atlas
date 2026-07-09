###############################################################################
## Figure 6 (v54 draft) — Subclustering and analysis of myeloid populations in RV failure
##
## Panels (from RV_snRNASeq_v54_draft.md figure legends):
##   (A) Myeloid subclustering UMAP with 6 subtypes (rMac CCR2-, iMac, Mono/Mono-derived,
##       TREM2+, DC, prolif)
##   (B) Dot plot of canonical marker genes per subtype (mean expr color, % size)
##   (C) Cell-type frequency across NF/pRV/RVF (* = KW p<0.05)
##   (D) Module dot plot of 4 myeloid-relevant bulk WGCNA modules (M1, M3, M4, M8)
##       by myeloid cell type
##   (E) ChEA TF enrichment for module genes UPregulated in RVF vs NF
##   (F) ChEA TF enrichment for module genes DOWNregulated in RVF vs NF
##   (G) Pseudobulk module-score violins for the 4 programs (GR-homeostatic + HIF/vascular
##       drift in CCR2- rMac, NF-kB MHCII/inflammasome in iMac, IFNg antigen presentation
##       pooled) across NF/pRV/RVF
##   (H) Companion dot plots of canonical genes from each program, restricted to home
##       cluster, across NF/pRV/RVF (shared expr/% legends)
##   (I) Spatial localization of myeloid subtypes in the full Xenium dataset
##       (n = 625,305 cells); non-myeloid grey; CCR2- rMacs and Macrophage C1q combined.
##
## Output: ./output/v52_figures/Figure_6_panel_<LETTER>.pdf (per-panel)
###############################################################################

source('./helper_scripts/_shared_helpers.R')
## Working-repo adaptation: per-figure output dir + standardized names.
V52_FIG_DIR <- './output/Figure_6'
dir.create(V52_FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create('./output/Xenium', showWarnings = FALSE, recursive = TRUE)

## Bootstrap the AUTHORITATIVE module x cell-type ChEA enrichment cache
## from dependencies/ into ./output (mirrors the F1 cache bootstrap).
## This cache was generated with the original enrichR; rebuilding it with
## the current enrichR diverges from the published 6E/6F, so on a clean
## ./output it MUST be restored from the committed dependency, not recomputed.
local({
  .src <- './dependencies/Figure_6/fig6_module_celltype_enrich.rds'
  .dst <- './output/fig6_module_celltype_enrich.rds'
  if (file.exists(.src) && !file.exists(.dst)) file.copy(.src, .dst)
})

## Composite figure dimensions (inches) — used by theme_v52() scaling
COMP_W <- 12
COMP_H <- 10

## Publication scales (geom widths, point sizes, text sizes) scaled to COMP_W.
## Use PS$geom_lw, PS$umap_pt, PS$scatter_pt, PS$text_mm, PS$dot_range, PS$legend_key.
PS <- pub_scales(COMP_W)

## ── Panel A cartoon (heart schematic / LV re-analysis) ────────────────────
## To regenerate cropped cartoon from the published Figure_5.pdf, run in R:
##   library(magick)
##   img <- image_read_pdf('~/Downloads/hdWGCNA_TOM/Manuscripts/Figure_5.pdf',
##                         pages = 1, density = 300)
##   ## Approx. geometry for top-left heart cartoon (17 cm composite ≈ 2008 px wide):
##   ## crop "WIDTHxHEIGHT+XOFF+YOFF" (px): heart cartoon spans ~ upper-left 18% width x 10% height
##   cart <- image_crop(img, '360x200+50+30')
##   image_write(cart, './new_scripts/assets/Figure_5_panel_A_heart_cartoon.png',
##               format = 'png', density = 300)
## Then compose with `insert_asset('Figure_5_panel_A_heart_cartoon.png')` in your layout.

library(Seurat)
## Robust AddModuleScore: Seurat's default nbin=24 (and any fixed nbin)
## throws "Insufficient data values to produce N bins" when the assay's
## avg-expression has too few distinct values (ties/sparse). Retry with a
## descending nbin ladder until it succeeds.
AddModuleScore <- function(object, ...) {
  .a <- list(...); .a[['nbin']] <- NULL
  .fn <- get('AddModuleScore', envir = asNamespace('Seurat'))
  for (.n in c(24L, 15L, 10L, 6L, 4L, 3L, 2L)) {
    .r <- tryCatch(do.call(.fn, c(list(object), .a, list(nbin = .n))),
                   error = function(e) NULL)
    if (!is.null(.r)) return(.r)
  }
  do.call(.fn, c(list(object), .a, list(nbin = 2L)))
}
library(hdWGCNA)
library(ggeasy)
library(harmony)
library(dplyr)



source('./helper_scripts/spatial_functions.R')


#######################################
###  FRESH MYELOID REANALYSIS (v52) ###
#######################################

M1 <- readRDS(file = "./dependencies/shared/myeloid_subclust.rds")

# Re-run SCTransform on RNA assay (original had multi-model SCT)
DefaultAssay(M1) <- "RNA"
.cache_sct_clustered <- './output/fig6_sct_clustered_cache.rds'
if (file.exists(.cache_sct_clustered)) {
  message('Loading cached SCTransform/Harmony/UMAP/clustering...')
  M1 <- readRDS(.cache_sct_clustered)
} else {
  M1 <- SCTransform(M1, verbose=FALSE)
  M1 <- RunPCA(M1, npcs=50, verbose=FALSE)
  M1 <- RunHarmony(M1,'patient',theta=4)
  M1 <- RunUMAP(M1,reduction = "harmony", dims = 1:30, verbose = F)
  M1 <- FindNeighbors(M1, reduction = "harmony", dims=1:30)
  M1 <- FindClusters(M1, resolution=1, verbose=FALSE)
  saveRDS(M1, .cache_sct_clustered)
}

print(table(Idents(M1), M1$patient))  

cat("Cluster counts:\n")
print(table(Idents(M1)))

# FindAllMarkers
.cache_fresh_markers <- './output/fig6_fresh_markers_cache.rds'
if (file.exists(.cache_fresh_markers)) {
  message('Loading cached FindAllMarkers...')
  markers <- readRDS(.cache_fresh_markers)
} else {
  markers <- FindAllMarkers(M1, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.5, verbose=FALSE)
  saveRDS(markers, .cache_fresh_markers)
}
write.csv(markers, "/tmp/myeloid_fresh_markers.csv", row.names=FALSE)

#######################################
###  CLUSTER IDENTITY ASSIGNMENTS   ###
#######################################

# Assign identities based on clean marker analysis
# Assign identities based on harmony-integrated marker analysis (14 clusters)
# Contamination: C1 (FB doublets DCN85%/PDGFRA61%), C6 (CM), C7 (weak), C8 (SMC), C9 (EC)
clean_ids <- c(
  "0"  = "CCR2- Resident Mac",          # HPGDS+ subtype — merged into single rMac
  "1"  = "Contamination",               # ACSM1, MEG3, DCLK1 — FB doublets (DCN 85%, PDGFRA 61%)
  "2"  = "CCR2- Resident Mac",          # broad rMac — F13A1, CD163L1, CD209, LYVE1, VSIG4
  "3"  = "CCR2- Resident Mac",          # RBPJ+ subtype — merged into single rMac
  "4"  = "Mac_Inflammatory",            # NR4A3, NR4A2, CD83, NR4A1, IL1B, CCL3
  "5"  = "Dendritic Cell",              # cDC2 — FLT3, HLA-DQA1, CIITA, CLEC10A
  "6"  = "Contamination",               # TRIM54, SCN5A, MYOM3 — cardiomyocytes
  "7"  = "Contamination",               # HTR4, NRK, SERPINE2 — weak myeloid identity
  "8"  = "Contamination",               # PTH1R, GUCY1A2, NOTCH3 — SMC/pericyte
  "9"  = "Contamination",               # MCF2L, NOTCH4, VWF — endothelial
  "10" = "TREM2+ Mac",                  # PLA2G7, GPNMB, LPL, MITF
  "11" = "Proliferating",               # CDK1, TOP2A, MKI67
  "12" = "Monocyte / Mac_Mono_Derived", # FCN1, VCAN, CD300E + IL1R2, EREG, AREG
  "13" = "Dendritic Cell"               # cDC1 — BTLA, CLEC9A, IDO1 (merged with cDC2 per Koenig)
)

## WARNING (2026-06-14): clean_ids is pinned to the ORIGINAL 14-cluster FindClusters
## result (→ canonical 4245-cell object). FindClusters has NO set.seed and is
## non-deterministic: a fresh run can yield a DIFFERENT cluster count/numbering
## (observed: 15 clusters), which makes this map misfire and corrupts the object
## (wrong cell count + scrambled labels). Treat dependencies/shared/myeloid_subclust_new.rds
## (4245 cells, from SSD canonical) as FROZEN — do NOT delete the caches and re-derive.
## See helper_scripts/load_frozen_myeloid.R proposal before ever regenerating.
M1$cluster_id <- unname(clean_ids[as.character(Idents(M1))])
M1$cluster_id[is.na(M1$cluster_id)] <- "Contamination"
Idents(M1) <- "cluster_id"

cat("\n=== CLEAN CLUSTER COUNTS ===\n")
print(sort(table(Idents(M1)), decreasing=TRUE))

# Remove remaining contamination
M1_clean <- subset(M1, idents="Contamination", invert=TRUE)
cat(sprintf("\nAfter removing contamination: %d cells\n", ncol(M1_clean)))

.cache_clean_umap <- './output/fig5_clean_umap_cache.rds'
if (file.exists(.cache_clean_umap)) {
  message('Loading cached M1_clean UMAP...')
  M1_clean <- readRDS(.cache_clean_umap)
} else {
  M1_clean <- RunUMAP(M1_clean,reduction = "harmony", dims = 1:30, verbose = F)
  saveRDS(M1_clean, .cache_clean_umap)
}

                                                                                                             
  lineage_genes <- c(                                                                                             
    # Resident mac lineage (CCR2-, self-renewing, tissue-origin)                                                  
    "LYVE1","CD163L1","VSIG4","MARCO","F13A1","FOLR2","MRC1","TIMD4","HPGDS","RBPJ",                              
    # Monocyte lineage (CCR2+, bone marrow-derived)                                                               
    "FCN1","VCAN","CD300E","S100A8","S100A9","CCR2","SELL", 
    # Pan-myeloid                                                                                                 
    "CD163","CD68","CSF1R"                                                                                        
  )                                                                                                               
  lineage_genes <- lineage_genes[lineage_genes %in% rownames(M1_clean)]                                           
                                                                                                                  
  pdf("./output/Xenium/myeloid_mac_inflam_lineage.pdf", width=11, height=4)                                       
  DotPlot(M1_clean, features=lineage_genes, assay="RNA", col.min=0, col.max=2) +
    RotatedAxis() +                                                                                               
    ggtitle("Mac_Inflammatory lineage identity: resident vs monocyte")                                            
  dev.off()    


## FREEZE the canonical myeloid object. FindClusters above is non-deterministic
## (no set.seed) and can drift (e.g. 15 vs 14 clusters), which corrupts cell count
## and labels. The SSD-canonical 4245-cell object is the source of truth: if present,
## LOAD it for all downstream panels and never overwrite it. Only save if absent.
.myeloid_canon <- "./dependencies/shared/myeloid_subclust_new.rds"
if (file.exists(.myeloid_canon)) {
  message("Using FROZEN canonical myeloid object (4245 cells); not overwriting.")
  M1_clean <- readRDS(.myeloid_canon); Idents(M1_clean) <- "cluster_id"
} else {
  saveRDS(M1_clean, .myeloid_canon)
}

# Reorder identities for plotting
ident_order <- c("CCR2- Resident Mac",
                 "Mac_Inflammatory",
                 "Monocyte / Mac_Mono_Derived",
                 "TREM2+ Mac",
                 "Dendritic Cell",
                 "Proliferating")
Idents(M1_clean) <- factor(Idents(M1_clean), levels=ident_order)

# DotPlot 1: Resident mac markers
res_genes <- c("LYVE1","FOLR2","VSIG4","CD163","F13A1","MRC1","CD163L1",
               "MARCO","HPGDS","RBPJ","TIMD4","CSF1R","CD68","CD14","C1QA","C1QB")
res_genes <- unique(res_genes[res_genes %in% rownames(M1_clean)])


pdf("./output//myeloid_audit_dotplot_resident.pdf", width=14, height=6)
DotPlot(M1_clean, features=res_genes, assay="RNA") +
  RotatedAxis() +
  ggtitle("Resident macrophage markers")
dev.off()
cat("Saved /tmp/myeloid_audit_dotplot_resident.pdf\n")

# DotPlot 2: Inflammatory / monocyte / TREM2 markers
inflam_genes <- c("CCR2","IL1B","CCL3","CCL4","CXCL3","CXCL8",
                   "FCN1","VCAN","S100A8","S100A9","CD300E",
                   "TREM2","GPNMB","SPP1","MITF","PLA2G7","LPL",
                   "TIMP1","HMOX1","FOSL1",
                   "NR4A1","NR4A2","CD83")
inflam_genes <- unique(inflam_genes[inflam_genes %in% rownames(M1_clean)])

pdf("./output/myeloid_audit_dotplot_inflam.pdf", width=16, height=6)
DotPlot(M1_clean, features=inflam_genes, assay="RNA") +
  RotatedAxis() +
  ggtitle("Inflammatory / monocyte / TREM2 markers")
dev.off()
cat("Saved /tmp/myeloid_audit_dotplot_inflam.pdf\n")

# DotPlot 3: DC markers
dc_genes <- c("FLT3","CIITA","HLA-DRA","HLA-DQA1","HLA-DQB1","CLEC10A","ITGAX",
              "CD1C","CLEC9A","BTLA","IDO1","IDO2","SLAMF7","FCER1A")
dc_genes <- unique(dc_genes[dc_genes %in% rownames(M1_clean)])

pdf("./output/myeloid_audit_dotplot_dc.pdf", width=14, height=6)
DotPlot(M1_clean, features=dc_genes, assay="RNA") +
  RotatedAxis() +
  ggtitle("DC / MHCII markers")
dev.off()
cat("Saved /tmp/myeloid_audit_dotplot_dc.pdf\n")

# DotPlot 4: All key markers combined (comprehensive audit)
all_audit <- c("LYVE1","VSIG4","CD163","F13A1","MRC1","CD163L1","MARCO","HPGDS","RBPJ",
               "CCR2","TIMP1","HMOX1",
               "FCN1","VCAN","CD300E",
               "TREM2","GPNMB","MITF","PLA2G7",
               "IL1B","CCL3","CCL4","CXCL8",
               "FLT3","CLEC10A","CLEC9A","BTLA","IDO1","CD83","NR4A1",
               "TOP2A","CDK1")
all_audit <- unique(all_audit[all_audit %in% rownames(M1_clean)])

pdf("./output/myeloid_audit_dotplot_all.pdf", width=20, height=6)
DotPlot(M1_clean, features=all_audit, assay="RNA") +
  RotatedAxis() +
  ggtitle("Comprehensive myeloid marker audit")
dev.off()
cat("Saved /tmp/myeloid_audit_dotplot_all.pdf\n")

# UMAP of clean object only
pdf("./output/myeloid_audit_umap_clean.pdf", width=10, height=8)
DimPlot(M1_clean, reduction="umap", label=TRUE, repel=TRUE, pt.size=PS$umap_pt) +
  ggtitle("Myeloid populations (contamination removed)")
dev.off()
cat("Saved /tmp/myeloid_audit_umap_clean.pdf\n")

cat("\n=== CLEAN CLUSTER COUNTS ===\n")
print(sort(table(Idents(M1_clean)), decreasing=TRUE))

# Audit confirmed: one CCR2- Resident Mac population. Proceeding with figure generation.

# ============================================================================
# Figure 6 — comprehensive panel rewrite using M1_clean (new 6-cluster scheme)
# Cluster identities (set in M1_clean@meta.data$cluster_id):
#   CCR2- Resident Mac, Mac_Inflammatory, Monocyte / Mac_Mono_Derived,
#   TREM2+ Mac, Dendritic Cell, Proliferating
# All panels in this section USE M1_clean directly — NO reload of the old
# myeloid_subclust.rds or the 21.7 GB scWGCNA_bulk2sn_projection.rds. The
# WGCNA panel (E) reads only the cached MEs object.
# ============================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(ggpubr)
  library(scales); library(viridis); library(tibble); library(cowplot)
})

.myeloid_levels <- c('CCR2- Resident Mac','Mac_Inflammatory',
                      'Monocyte / Mac_Mono_Derived','TREM2+ Mac',
                      'Dendritic Cell','Proliferating')

.myeloid_pal <- c(
  'CCR2- Resident Mac'         = '#1F77B4',
  'Mac_Inflammatory'           = '#FF7F0E',
  'Monocyte / Mac_Mono_Derived'= '#2CA02C',
  'TREM2+ Mac'                 = '#D62728',
  'Dendritic Cell'             = '#9467BD',
  'Proliferating'              = '#8C564B'
)

.disease_pal <- c('NF'='#4DAF4A','pRV'='#FF7F00','RVF'='#E41A1C')

# Ensure cluster_id is a factor in canonical order; set Idents
M1_clean$cluster_id <- factor(M1_clean$cluster_id, levels = .myeloid_levels)
Idents(M1_clean) <- 'cluster_id'


#######################################
#############  FIGURE 6A  #############
####### UMAP — 6 myeloid clusters #####
#######################################

pdf('./output/Myeloid_snUMAP.pdf', width = 5, height = 5)
print(PlotEmbedding(M1_clean, group.by = 'cluster_id',
                    point_size = PS$umap_pt, plot_under = TRUE,
                    plot_theme = umap_theme() + NoLegend(),
                    raster_dpi = 400, raster_scale = 0.5))
dev.off()
cat('Wrote: Myeloid_snUMAP.pdf\n')

# Standalone legend
.legend_p <- ggplot(data.frame(x = 1, y = seq_along(.myeloid_levels),
                                cluster = factor(.myeloid_levels, .myeloid_levels)),
                    aes(x, y, color = cluster)) +
  geom_point(size = 4) +
  scale_color_manual(values = .myeloid_pal, name = 'Subtype') +
  theme_void() + theme(legend.position = 'right')
pdf('./output/Myeloid_snUMAP_legend.pdf', width = 3.5, height = 4)
print(cowplot::ggdraw(cowplot::get_legend(.legend_p)))
dev.off()
cat('Wrote: Myeloid_snUMAP_legend.pdf\n')


#######################################
#############  FIGURE 6B  #############
####### Marker-gene DotPlot ###########
#######################################

.myeloid_markers <- c(
  # Resident / TLF (5)
  'TIMD4','LYVE1','FOLR2','F13A1','VSIG4',
  # Inflammatory (3)
  'NR4A1','IL1B','CCL3',
  # Monocyte (4) — CCR2 included as the negative-result marker that defines
  # 'CCR2- Resident Mac' and contrasts with the Monocyte cluster
  'FCN1','VCAN','CD300E','CCR2',
  # TREM2+ (3)
  'TREM2','GPNMB','SPP1',
  # DC (3)
  'FLT3','CIITA','CLEC10A',
  # Proliferating (2)
  'TOP2A','MKI67'
)
.myeloid_markers <- intersect(.myeloid_markers, rownames(M1_clean))

pdf('./output/Myeloid_dot.pdf', width = 8.4, height = 4)
print(DotPlot(M1_clean, features = .myeloid_markers, assay = 'RNA',
              col.min = 0, col.max = 2,
              cols = c('lightgrey', 'blue')) +
      xlab('') + ylab('') +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                        face = 'italic')))
dev.off()
cat('Wrote: Myeloid_dot.pdf\n')


## Note: prior Panel C (Myeloid_lineage_dotplot.pdf) was dropped because it
## restated information already in Panel B; the only unique marker (CCR2)
## was folded into Panel B's monocyte block instead. The legacy PDF is left
## on disk but no longer regenerated — see git history for the previous
## version if needed.


#######################################
#############  FIGURE 6C  #############
####### Cluster proportions by group ##
####### (v54: KW p<0.05 *)            #
#######################################

## S2C-format proportion plots: per-patient cluster proportions, with two
## complementary views (single-row dodged boxplot + facet-wrapped boxplot
## with Kruskal–Wallis p-values), matching Supplementary Figure 2 panel C.

# Per-patient cluster proportions (zero-padded for missing combos)
.pat_prop <- M1_clean@meta.data %>%
  dplyr::count(patient, group, cluster_id, name = 'n') %>%
  dplyr::group_by(patient, group) %>%
  dplyr::mutate(prop = n / sum(n)) %>%
  dplyr::ungroup()
.all_combos <- expand.grid(patient = unique(.pat_prop$patient),
                            cluster_id = .myeloid_levels,
                            stringsAsFactors = FALSE)
.all_combos <- merge(.all_combos,
                     unique(.pat_prop[, c('patient','group')]),
                     by = 'patient')
.pat_prop_full <- merge(.all_combos, .pat_prop,
                         by = c('patient','group','cluster_id'),
                         all.x = TRUE)
.pat_prop_full$prop[is.na(.pat_prop_full$prop)] <- 0
.pat_prop_full$cluster_id <- factor(.pat_prop_full$cluster_id,
                                     levels = .myeloid_levels)
.pat_prop_full$group      <- factor(.pat_prop_full$group,
                                     levels = c('NF','pRV','RVF'))

# View 1 (S2C-style compact): single-row boxplot, one box per cluster,
# dodged and coloured by disease group.
pdf('./output/Myeloid_clust_counts.pdf', width = 6.5, height = 3)
print(ggplot(.pat_prop_full,
             aes(x = cluster_id, y = prop, color = group)) +
      geom_boxplot() +
      scale_color_manual(values = .disease_pal) +
      ylab('Proportion of myeloid cells') + xlab('') +
      theme_classic())
dev.off()
cat('Wrote: Myeloid_clust_counts.pdf\n')

# View 2 (S2C-style facet-wrap): one facet per cluster, KW p-values.
pdf('./output/Myeloid_clust_freq_stats.pdf', width = 12.5, height = 8)
p <- ggboxplot(.pat_prop_full, x = 'group', y = 'prop',
               fill = 'group', group = 'group') +
  theme_classic() +
  theme(axis.text.x  = element_text(size = 14),
        axis.text.y  = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text  = element_text(size = 14),
        text         = element_text(color = 'black'),
        axis.text    = element_text(color = 'black')) +
  scale_fill_manual(values = .disease_pal) +
  facet_wrap(~ cluster_id, ncol = 3, scales = 'free_y') +
  stat_compare_means(aes(group = group), method = 'kruskal.test')
print(p)
dev.off()
cat('Wrote: Myeloid_clust_freq_stats.pdf\n')

# Stacked bar (kept as a complementary summary plot)
.prop_df <- M1_clean@meta.data %>%
  dplyr::count(group, cluster_id, name = 'n') %>%
  dplyr::group_by(group) %>%
  dplyr::mutate(prop = n / sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(group      = factor(group, levels = c('NF','pRV','RVF')),
                cluster_id = factor(cluster_id, levels = .myeloid_levels))
pdf('./output/Myeloid_prev_stacked.pdf', width = 5, height = 4.5)
print(ggplot(.prop_df, aes(x = group, y = prop, fill = cluster_id)) +
      geom_bar(stat = 'identity') +
      scale_fill_manual(values = .myeloid_pal, name = '') +
      ylab('Proportion of myeloid cells') + xlab('') +
      theme_classic(base_size = 11) +
      theme(legend.position = 'right'))
dev.off()
cat('Wrote: Myeloid_prev_stacked.pdf\n')

# Console summary (Kruskal-Wallis per cluster)
.kw_results <- .pat_prop_full %>%
  dplyr::group_by(cluster_id) %>%
  dplyr::summarise(kw_p = tryCatch(kruskal.test(prop ~ group)$p.value,
                                   error = function(e) NA_real_),
                   .groups = 'drop')
cat('\n=== Kruskal-Wallis cluster-proportion p-values (per disease group) ===\n')
print(.kw_results)


#######################################
#############  FIGURE 6D  #############
####### Bulk WGCNA module dot plot ####
####### M1/M3/M4/M8 by myeloid type ###
####### Canonical normalization:      #
####### AddModuleScore + Seurat       #
####### DotPlot (per-feature z-scaled,#
####### % expressing) — matches v57   #
####### legend + legacy Figure_6.R    #
#######################################
##
## v57 legend D: "Module expression dot plot of the four myeloid-relevant
## bulk WGCNA modules (M1, M3, M4, M8) by myeloid cell type."  The ORIGINAL
## (legacy ~/Downloads/hdWGCNA_TOM/Manuscripts/Figure_6.R) builds this by
## AddModuleScore() of the bulk-module gene sets onto the myeloid object
## then Seurat::DotPlot (col.min/col.max clamp the z-scored average score;
## dot size = % expressing).  The new_scripts variant instead z-scored
## module EIGENGENES per cluster — a different quantity AND normalization;
## that is the "normalizing incorrectly" issue, replaced here.  Needs only
## dependencies/shared/bulk_heart_modules.csv (no fig5_MEs_cache).

.bulk_mod_csv <- './dependencies/shared/bulk_heart_modules.csv'
if (!file.exists(.bulk_mod_csv)) {
  cat('Panel D: bulk_heart_modules.csv missing — skipping.\n')
} else {
  .cm <- read.csv(.bulk_mod_csv, stringsAsFactors = FALSE)[, 1:3]   # gene_name,module,color
  .cm <- subset(.cm, gene_name %in% rownames(M1_clean))
  .cm <- .cm[match(unique(.cm$gene_name), .cm$gene_name), ]
  .map      <- WGCNA::labels2colors(1:100)            # M1=turquoise, M2=blue, ...
  .want_M   <- c('M1','M3','M4','M8')
  .want_col <- .map[as.integer(sub('^M', '', .want_M))]
  .gene_lists <- lapply(.want_col, function(cc) .cm$gene_name[.cm$module == cc])
  names(.gene_lists) <- .want_M
  .ok         <- vapply(.gene_lists, function(g) length(g) >= 3, logical(1))
  .gene_lists <- .gene_lists[.ok]
  .want_M     <- names(.gene_lists)
  cat('Panel D: module gene-set sizes —',
      paste(.want_M, lengths(.gene_lists), sep = ':', collapse = '  '), '\n')

  if (length(.gene_lists) > 0) {
    M1_clean <- AddModuleScore(M1_clean, features = .gene_lists, name = 'bulkmod_')
    .src <- paste0('bulkmod_', seq_along(.gene_lists))
    .dst <- paste0('module_', .want_M)
    names(M1_clean@meta.data)[match(.src, names(M1_clean@meta.data))] <- .dst

    ## 6D pixel-match: 5 collapsed myeloid types (drop Proliferating),
    ## reference display names, order CCR2- rMac (top) ... DC (bottom),
    ## blue-white-red diverging clamped to [-1, 1].
    .lab6d <- c('CCR2- Resident Mac'          = 'CCR2⁻ rMac',
                'Mac_Inflammatory'            = 'iMac',
                'Monocyte / Mac_Mono_Derived' = 'Mono /\nMono-\nDerived',
                'TREM2+ Mac'                  = 'TREM2⁺ Mac',
                'Dendritic Cell'              = 'DC')
    .top2bot <- c('CCR2⁻ rMac','iMac','Mono /\nMono-\nDerived',
                  'TREM2⁺ Mac','DC')
    .m6d <- subset(M1_clean, subset = cluster_id %in% names(.lab6d))
    .m6d$.d6 <- factor(unname(.lab6d[as.character(.m6d$cluster_id)]),
                       levels = rev(.top2bot))   # rev: ggplot puts L1 at bottom
    Idents(.m6d) <- '.d6'
    p_6D <- DotPlot(.m6d, features = .dst,
                    dot.min = 0, col.min = -1, col.max = 1,
                    scale.min = 50, scale.max = 100) +
      scale_x_discrete(labels = .want_M) +
      scale_color_gradient2(low = 'blue', mid = 'white', high = 'red',
                            midpoint = 0, limits = c(-1, 1),
                            name = 'Avg. module\nscore (z)') +
      xlab('Bulk WGCNA module') + ylab('') +
      theme_classic() +
      theme(panel.border = element_rect(fill = NA, color = 'black'),
            axis.line    = element_blank(),
            axis.text.y  = element_text(lineheight = 0.8))
    pdf('./output/Myeloid_module_dotplot.pdf', width = 6, height = 4)
    print(p_6D)
    dev.off()
    cat('Wrote: Myeloid_module_dotplot.pdf (AddModuleScore + Seurat DotPlot)\n')

    # Companion: same module scores by cluster x disease group (Phase trends)
    .cg <- M1_clean@meta.data %>%
      tibble::rownames_to_column('cell') %>%
      dplyr::select(cell, cluster_id, group, dplyr::all_of(.dst)) %>%
      tidyr::pivot_longer(dplyr::all_of(.dst),
                          names_to = 'module', values_to = 'score') %>%
      dplyr::mutate(module = sub('^module_', '', module)) %>%
      dplyr::group_by(cluster_id, group, module) %>%
      dplyr::summarise(mean_score = mean(score, na.rm = TRUE), .groups = 'drop') %>%
      dplyr::mutate(module     = factor(module, levels = .want_M),
                    cluster_id = factor(cluster_id, levels = .myeloid_levels),
                    group      = factor(group, levels = c('NF','pRV','RVF')))
    pdf('./output/Myeloid_module_by_group.pdf', width = 8, height = 4)
    print(ggplot(.cg, aes(x = group, y = cluster_id, fill = mean_score)) +
          geom_tile() +
          facet_wrap(~ module, nrow = 1) +
          scale_fill_gradient2(low = 'blue', mid = 'grey95', high = 'red',
                               midpoint = 0, name = 'mean module\nscore') +
          xlab('') + ylab('') +
          theme_classic(base_size = 10) +
          theme(strip.text = element_text(face = 'bold'),
                axis.text.x = element_text(angle = 30, hjust = 1)))
    dev.off()
    cat('Wrote: Myeloid_module_by_group.pdf\n')
  }
}


#######################################
###########  FIGURE 6E + 6F  ###########
####### ChEA TF enrichment for ########
####### bulk-WGCNA module genes  ######
####### filtered by per-cell-type DE: #
####### 6E = UP in RVF vs NF          #
####### 6F = DOWN in RVF vs NF        #
#######################################
## Strategy mirrors the original ./Figure_6.R: for each myeloid bulk WGCNA
## module (M1, M3, M4, M8) × cell_type × disease comparison, find DEGs
## within the module's gene set, then run enrichR on UP/DOWN genes
## separately. Plot top enrichment terms per (module × cell_type) on a
## single dot plot per database/contrast/direction.

suppressPackageStartupMessages({ library(enrichR); library(WGCNA); library(stringr); library(viridis) })

.bulk_modules_csv <- './dependencies/shared/bulk_heart_modules.csv'
if (!file.exists(.bulk_modules_csv)) {
  stop('bulk_heart_modules.csv not found — needed for module-anchored enrichment')
}
bulk_modules <- read.csv(.bulk_modules_csv, stringsAsFactors = FALSE)
# bulk_heart_modules.csv has module column with WGCNA color names; map to M-numbers
.color_to_M <- setNames(paste0('M', seq_along(WGCNA::labels2colors(1:100))),
                         WGCNA::labels2colors(1:100))
bulk_modules$M <- .color_to_M[bulk_modules$module]

# Myeloid-assigned bulk modules
.myeloid_mods <- c('M1','M3','M4','M8')

# Composite identity: <cluster>_<group> + a group-only identity for module-only views
M1_clean$cluster_group <- paste0(as.character(M1_clean$cluster_id),
                                  '_', M1_clean$group)

.dbs_to_query <- c('Reactome_2022','ChEA_2022')
.contrasts    <- list(c('RVF','NF'), c('RVF','pRV'))
.cell_types   <- levels(M1_clean$cluster_id)

wrapText <- function(x, len)
  vapply(x, function(y) paste(strwrap(y, len), collapse = '\n'), character(1),
         USE.NAMES = FALSE)

## ── Run enrichments (module × cell_type × contrast) ─────────────────────
## Authoritative original cache (matches new_scripts/Figure_6.R exactly):
## present in ./output (and on the SSD), so this loads it rather than
## rebuilding — the rebuilt/modern-enrichR ordering diverged from the
## published 6E/6F, so do NOT bump this name.
.cache_module_celltype_enrich <-
  './output/fig6_module_celltype_enrich.rds'

if (file.exists(.cache_module_celltype_enrich)) {
  cat('Loading cached module × cell_type enrichment...\n')
  combined_output <- readRDS(.cache_module_celltype_enrich)
} else {
  cat('Running module × cell_type enrichment (this hits enrichR online)...\n')
  Idents(M1_clean) <- 'cluster_group'
  combined_output <- data.frame()

  for (i_M in .myeloid_mods) {
    key_genes <- bulk_modules$gene_name[bulk_modules$M == i_M]
    key_genes <- key_genes[key_genes %in% rownames(M1_clean)]
    if (length(key_genes) < 5) next

    for (j in .cell_types) {
      for (k in .contrasts) {
        id1 <- paste0(j, '_', k[1]); id2 <- paste0(j, '_', k[2])
        gene_set <- tryCatch(
          FindMarkers(M1_clean, ident.1 = id1, ident.2 = id2,
                       features = key_genes, verbose = FALSE),
          error = function(e) NULL)
        if (is.null(gene_set) || nrow(gene_set) == 0) next
        gene_set <- subset(gene_set, p_val_adj < 0.05)
        if (nrow(gene_set) == 0) next

        for (dir in c('up','down')) {
          gset <- if (dir == 'up') subset(gene_set, avg_log2FC > 0)
                  else            subset(gene_set, avg_log2FC < 0)
          if (nrow(gset) < 3) next
          enriched <- tryCatch(enrichR::enrichr(rownames(gset), .dbs_to_query),
                               error = function(e) NULL)
          Sys.sleep(2)
          if (is.null(enriched)) next
          for (db in names(enriched)) {
            cur <- enriched[[db]]
            if (is.null(cur) || nrow(cur) < 1) next
            cur$db         <- db
            cur$module     <- i_M
            cur$celltype   <- j
            cur$comparison <- paste0(k[1], '_', k[2])
            cur$direction  <- dir
            combined_output <- rbind(combined_output, cur)
          }
        }
      }
    }
  }
  saveRDS(combined_output, .cache_module_celltype_enrich)
  cat('  enrichments cached: ', nrow(combined_output), 'rows\n')
}

## ── Helper: build a single dot plot for one (db, comparison, direction) ──
## Layout (matches reference Figure 6F/G style):
##   y-axis  = Term (TF / Reactome pathway)
##   x-axis  = module_celltype (e.g., M1_CCR2- Resident Mac)
##   colour  = -log(P)             (magma palette, capped at max_p)
##   size    = log(Combined.Score)
##   bottom  = thin coloured strip grouping x-axis by WGCNA module colour
##             (M1=turquoise, M3=brown, M4=yellow, M8=pink)
.module_celltype_dotplot <- function(df, db, comparison, direction,
                                      out_path, width = 8, height = 6,
                                      title = NULL, max_p = 10) {
  sel <- df[df$db == db &
            df$comparison == comparison &
            df$direction == direction, , drop = FALSE]
  if (nrow(sel) == 0) {
    cat('  skip', basename(out_path), '(no rows)\n'); return(invisible())
  }
  sel <- subset(sel, P.value < 0.05)
  if (nrow(sel) == 0) {
    cat('  skip', basename(out_path), '(none passing P<0.05)\n'); return(invisible())
  }
  sel$module_celltype <- paste0(sel$module, '_', sel$celltype)

  # take the top hit per (module × celltype)
  idx_top <- match(unique(sel$module_celltype), sel$module_celltype)
  sel <- sel[idx_top, ]

  sel$logp <- -log(sel$P.value)
  sel$logp <- pmin(sel$logp, max_p)
  sel$Term <- str_replace(sel$Term, '\\ R-HSA.*', '')
  sel$Term <- str_replace(sel$Term, '\\ \\(.*$',  '')
  sel$wrap <- wrapText(sel$Term, 35)

  sel <- sel[order(sel$module, sel$celltype), ]
  mc_levels <- unique(sel$module_celltype)
  sel$module_celltype <- factor(sel$module_celltype, levels = mc_levels)
  sel$wrap <- factor(sel$wrap, levels = rev(unique(sel$wrap)))

  # Module colour bar data (one tile per x-axis position)
  .mod_color_map <- c('M1' = 'turquoise', 'M3' = 'brown',
                       'M4' = 'yellow',    'M8' = 'pink')
  cb_df <- data.frame(
    module_celltype = mc_levels,
    module = sub('_.*$', '', mc_levels),
    stringsAsFactors = FALSE)
  cb_df$module_color <- .mod_color_map[cb_df$module]
  cb_df$module_celltype <- factor(cb_df$module_celltype, levels = mc_levels)

  p_dot <- ggplot(sel, aes(x = module_celltype, y = wrap,
                            color = logp, size = log(Combined.Score))) +
    geom_point() +
    scale_color_stepsn(colors = rev(magma(256)), name = 'logp') +
    scale_size_continuous(name = 'log(Score)', range = c(2, 8)) +
    xlab('') + ylab('') +
    ggtitle(title %||% sprintf('%s — %s — %s',
                                  db, comparison, direction)) +
    theme_bw() +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          plot.title   = element_text(size = 11, face = 'bold'),
          panel.grid.minor = element_blank(),
          plot.margin  = margin(5, 5, 0, 5))

  p_cb <- ggplot(cb_df, aes(x = module_celltype, y = 1,
                             fill = module_color)) +
    geom_tile() +
    scale_fill_identity() +
    xlab('') + ylab(NULL) +
    theme_void() +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                       size = 8),
          axis.ticks.x = element_blank(),
          plot.margin  = margin(0, 5, 5, 5))
  # Show the actual module_celltype labels under the colour bar
  p_cb <- p_cb +
    scale_x_discrete(labels = mc_levels) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                      size = 7))

  if (!requireNamespace('patchwork', quietly = TRUE)) {
    pdf(out_path, width = width, height = height)
    print(p_dot)
    dev.off()
    cat('  wrote', out_path, '(no patchwork; dotplot only)\n')
    return(invisible())
  }
  combined <- patchwork::wrap_plots(p_dot, p_cb, ncol = 1,
                                     heights = c(20, 2))
  pdf(out_path, width = width, height = height)
  print(combined)
  dev.off()
  cat('  wrote', out_path, '\n')
}

# Define %||% if not present
`%||%` <- function(a, b) if (!is.null(a) && !is.na(a) && nzchar(a)) a else b

#######################################
####  Panel F — Reactome (4 plots) ####
#######################################
.module_celltype_dotplot(combined_output, 'Reactome_2022', 'RVF_NF', 'up',
                          './output/Myeloid_reactome_terms_cell_type_up_RVF_vs_NF.pdf',
                          width = 8, height = 7)
.module_celltype_dotplot(combined_output, 'Reactome_2022', 'RVF_NF', 'down',
                          './output/Myeloid_reactome_terms_cell_type_down_RVF_vs_NF.pdf',
                          width = 6.6, height = 5)
.module_celltype_dotplot(combined_output, 'Reactome_2022', 'RVF_pRV', 'up',
                          './output/Myeloid_reactome_terms_cell_type_up_RVF_vs_pRV.pdf',
                          width = 8, height = 7)
.module_celltype_dotplot(combined_output, 'Reactome_2022', 'RVF_pRV', 'down',
                          './output/Myeloid_reactome_terms_cell_type_down_RVF_vs_pRV.pdf',
                          width = 6.6, height = 5)

#######################################
####  Panel G — ChEA (4 plots) #########
#######################################
.module_celltype_dotplot(combined_output, 'ChEA_2022', 'RVF_NF', 'up',
                          './output/Myeloid_chea_terms_cell_type_up_RVF_vs_NF.pdf',
                          width = 10, height = 5)
.module_celltype_dotplot(combined_output, 'ChEA_2022', 'RVF_NF', 'down',
                          './output/Myeloid_chea_terms_cell_type_down_RVF_vs_NF.pdf',
                          width = 6, height = 4)
.module_celltype_dotplot(combined_output, 'ChEA_2022', 'RVF_pRV', 'up',
                          './output/Myeloid_chea_terms_cell_type_up_RVF_vs_pRV.pdf',
                          width = 10, height = 5)
.module_celltype_dotplot(combined_output, 'ChEA_2022', 'RVF_pRV', 'down',
                          './output/Myeloid_chea_terms_cell_type_down_RVF_vs_pRV.pdf',
                          width = 6, height = 4)


#######################################
####### Module-only enrichment ########
####### (Panel H/I/J — module × group ##
####### without cell-type stratification)
#######################################
.cache_module_only_enrich <- './output/fig6_module_only_enrich.rds'
if (file.exists(.cache_module_only_enrich)) {
  cat('Loading cached module-only enrichment...\n')
  module_only_output <- readRDS(.cache_module_only_enrich)
} else {
  cat('Running module-only enrichment (group only, no cell type)...\n')
  Idents(M1_clean) <- 'group'
  module_only_output <- data.frame()
  for (i_M in .myeloid_mods) {
    key_genes <- bulk_modules$gene_name[bulk_modules$M == i_M]
    key_genes <- key_genes[key_genes %in% rownames(M1_clean)]
    if (length(key_genes) < 5) next
    for (k in .contrasts) {
      gene_set <- tryCatch(
        FindMarkers(M1_clean, ident.1 = k[1], ident.2 = k[2],
                     features = key_genes, verbose = FALSE),
        error = function(e) NULL)
      if (is.null(gene_set) || nrow(gene_set) == 0) next
      gene_set <- subset(gene_set, p_val_adj < 0.05)
      if (nrow(gene_set) == 0) next
      for (dir in c('up','down')) {
        gset <- if (dir == 'up') subset(gene_set, avg_log2FC > 0)
                else            subset(gene_set, avg_log2FC < 0)
        if (nrow(gset) < 3) next
        enriched <- tryCatch(enrichR::enrichr(rownames(gset), .dbs_to_query),
                              error = function(e) NULL)
        Sys.sleep(2)
        if (is.null(enriched)) next
        for (db in names(enriched)) {
          cur <- enriched[[db]]
          if (is.null(cur) || nrow(cur) < 1) next
          cur$db         <- db
          cur$module     <- i_M
          cur$comparison <- paste0(k[1], '_', k[2])
          cur$direction  <- dir
          module_only_output <- rbind(module_only_output, cur)
        }
      }
    }
  }
  saveRDS(module_only_output, .cache_module_only_enrich)
  cat('  module-only enrichments cached: ', nrow(module_only_output), 'rows\n')
}

# Plot module-only enrichment as a per-module bar plot, top 5 terms per
# (module × comparison × direction).
if (nrow(module_only_output) > 0) {
  .outdir <- './output/scMyeloid_subclust_enrichr_plot'
  dir.create(.outdir, showWarnings = FALSE, recursive = TRUE)

  for (i_M in unique(module_only_output$module)) {
    for (cmp in unique(module_only_output$comparison)) {
      for (dir in unique(module_only_output$direction)) {
        cur <- subset(module_only_output,
                      module == i_M & comparison == cmp & direction == dir)
        if (nrow(cur) == 0) next
        cur$wrap <- wrapText(cur$Term, 45)
        for (cur_db in unique(cur$db)) {
          plot_df <- subset(cur, db == cur_db)
          plot_df <- plot_df[order(-plot_df$Combined.Score), ][seq_len(min(5, nrow(plot_df))), ]
          plot_df$Combined.Score <- log(plot_df$Combined.Score)

          p_one <- ggplot(plot_df, aes(x = Combined.Score,
                                        y = reorder(wrap, Combined.Score))) +
            geom_bar(stat = 'identity', position = 'identity',
                     color = 'white', fill = 'green') +
            ggtitle(sprintf('%s | %s | %s | %s', i_M, cur_db, cmp, dir)) +
            xlab('Enrichment log(combined score)') + ylab('') +
            theme_classic(base_size = 9) +
            theme(plot.title = element_text(size = 9, face = 'bold'))

          fname <- sprintf('%s/%s_%s_%s_%s.pdf', .outdir,
                            i_M, gsub('_2022$','',cur_db), cmp, dir)
          pdf(fname, width = 6, height = 3)
          print(p_one)
          dev.off()
        }
      }
    }
  }
  cat('  wrote per-module enrichment bars in', .outdir, '\n')
}

# Reset Idents to cluster_id
Idents(M1_clean) <- 'cluster_id'


#######################################
######  FIGURE 6 — LEGACY EXTRA  ######
####### MHCII gene expression #########
####### across cluster × group ########
####### (NOT in v54 legend; kept   ####
####### as a companion supplemental) ##
#######################################

# Composite identity: cluster × group, ordered NF→pRV→RVF within each cluster
M1_clean$cluster_group <- paste0(as.character(M1_clean$cluster_id),
                                  '_', M1_clean$group)
.cg_levels <- as.vector(outer(.myeloid_levels, c('NF','pRV','RVF'),
                               function(a,b) paste0(a,'_',b)))
M1_clean$cluster_group <- factor(M1_clean$cluster_group, levels = .cg_levels)

.mhc_genes <- c('CIITA','HLA-DRA','HLA-DRB1','HLA-DPA1','HLA-DPB1',
                'HLA-DQA1','HLA-DQB1','CD74')
.mhc_genes <- intersect(.mhc_genes, rownames(M1_clean))

pdf('./output/Myeloid_MHCII_dotplot.pdf', width = 12, height = 5)
print(DotPlot(M1_clean, features = .mhc_genes, assay = 'RNA',
              group.by = 'cluster_group',
              col.min = 0, col.max = 2,
              cols = c('grey90', 'red')) +
      RotatedAxis() + xlab('') + ylab('') +
      ggtitle('MHCII gene expression by cluster × group') +
      theme(axis.text.y = element_text(size = 9),
            axis.text.x = element_text(face = 'italic')))
dev.off()
cat('Wrote: Myeloid_MHCII_dotplot.pdf\n')


#######################################
######  FIGURE 6 — LEGACY EXTRA  ######
####### NR3C1 target genes ############
####### (glucocorticoid axis) #########
####### (NOT in v54 legend; kept   ####
####### as a companion supplemental) ##
#######################################

.nr3c1_targets <- c('FKBP5','TSC22D3','DUSP1','SGK1','PER1','KLF9','ZBTB16')
.nr3c1_targets <- intersect(.nr3c1_targets, rownames(M1_clean))

pdf('./output/Myeloid_NR3C1_targets.pdf', width = 12, height = 5)
print(DotPlot(M1_clean, features = .nr3c1_targets, assay = 'RNA',
              group.by = 'cluster_group',
              col.min = 0, col.max = 2,
              cols = c('grey90', 'red')) +
      RotatedAxis() + xlab('') + ylab('') +
      ggtitle('NR3C1 target gene expression by cluster × group') +
      theme(axis.text.y = element_text(size = 9),
            axis.text.x = element_text(face = 'italic')))
dev.off()
cat('Wrote: Myeloid_NR3C1_targets.pdf\n')

# Also small NR3C1 itself per cluster × group
pdf('./output/myeloid_nr3c1.pdf', width = 4, height = 3)
print(VlnPlot(M1_clean, features = 'NR3C1', group.by = 'cluster_group',
              pt.size = 0) +
      NoLegend() + xlab('') +
      theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 7)))
dev.off()
cat('Wrote: myeloid_nr3c1.pdf\n')

# Reset Idents to cluster_id
Idents(M1_clean) <- 'cluster_id'


#####################################################################
########### FIGURE 6G + 6H — Myeloid transcriptional programs #######
###########  GR targets (down) / MHC-II (up) /                #######
###########  Inflammasome & Type-I IFN (pertinent negatives)  #######
###########  6G = pooled-myeloid pseudobulk box plots;        #######
###########  6H = gene-level companion dot plots              #######
#####################################################################
# Two coherent, pan-myeloid Phase-1 shifts survive per-gene/per-patient
# scrutiny: (1) collapse of the canonical GR negative-feedback program
# (FKBP5, TSC22D3/GILZ, ZBTB16, KLF9, ANGPTL4 — all confirmed down per-patient;
# divergent/confounded members SGK1 [MR-shared, rises], PER1/DUSP1 [MAPK/IEG],
# DDIT4 [HIF], KLF13/GADD45B/TFCP2L1 [non-GR] excluded), and (2) CIITA-driven
# MHC-II / antigen-presentation induction. Two innate axes expected in failure
# that do NOT engage are shown as pertinent negatives: the NLRP3 inflammasome
# and the type-I interferon (antiviral ISG) response — both flat NF->disease
# in pooled myeloid. ("disease" = pRV+RVF pooled throughout.)

suppressPackageStartupMessages({library(ggpubr)})

## --- gene sets (intersected to object) ---------------------------------
## GR-transactivation negative-feedback reporters, per-gene confirmed DOWN.
GR_targets   <- c('FKBP5','TSC22D3','ZBTB16','KLF9','ANGPTL4')
## MHC-II / antigen presentation (CIITA-driven), confirmed UP.
MHCII        <- c('CIITA','CD74','HLA-DRA','HLA-DRB1','HLA-DRB5','HLA-DPA1',
                  'HLA-DPB1','HLA-DQA1','HLA-DQB1','HLA-DMA','HLA-DMB','HLA-DOA')
## Pertinent negative #1: NLRP3 inflammasome machinery (flat).
Inflammasome <- c('NLRP3','NLRP1','AIM2','CASP1','PYCARD','IL1B','IL18',
                  'GSDMD','CASP4','NAIP')
## Pertinent negative #2: type-I interferon / antiviral ISGs (flat).
TypeI_IFN    <- c('ISG15','MX1','MX2','OAS1','OAS2','OAS3','IFIT1','IFIT3',
                  'IRF7','RSAD2','USP18','IFI6','IFI44')

.in_obj <- function(g) intersect(g, rownames(M1_clean))
prog_list <- list(GR_targets   = .in_obj(GR_targets),
                  MHCII        = .in_obj(MHCII),
                  Inflammasome = .in_obj(Inflammasome),
                  TypeI_IFN    = .in_obj(TypeI_IFN))

# Supplementary table -------------------------------------------------------
supp_tbl <- do.call(rbind, lapply(names(prog_list), function(p) {
  data.frame(program = p, gene = prog_list[[p]],
    expected_direction = switch(p,
      GR_targets   = 'down (disease<NF), pan-myeloid',
      MHCII        = 'up (disease>NF), pan-myeloid',
      Inflammasome = 'unchanged (pertinent negative)',
      TypeI_IFN    = 'unchanged (pertinent negative)'),
    source = switch(p,
      GR_targets   = 'Canonical GR negative-feedback transactivation reporters (Lavin 2014; Uhlenhaut 2013); each gene confirmed down per-patient pseudobulk',
      MHCII        = 'CIITA-driven MHC-II / antigen presentation (Reith 2005)',
      Inflammasome = 'NLRP3 inflammasome machinery (Latz 2013); flat NF->disease',
      TypeI_IFN    = 'Type-I IFN / antiviral ISG signature; flat NF->disease'),
    stringsAsFactors = FALSE)
}))
write.csv(supp_tbl, './output/SuppTable_Myeloid_Programs.csv', row.names = FALSE)
cat('Wrote: SuppTable_Myeloid_Programs.csv (', nrow(supp_tbl), 'rows )\n')

# Module scores -------------------------------------------------------------
set.seed(42)
## nbin default (24) fails ("Insufficient data values to produce 24
## bins") when the active assay has few binnable avg-expression values;
## scale nbin to the feature count.
.nb <- max(2L, min(24L, floor(nrow(M1_clean) / 50)))
M1_clean <- AddModuleScore(M1_clean, features = prog_list,
                           name = paste0('Score_', names(prog_list)),
                           ctrl = 50, nbin = .nb, search = FALSE)
score_cols <- setNames(paste0('Score_', names(prog_list), seq_along(prog_list)),
                       names(prog_list))

# Heatmap of mean module score per cluster x condition ----------------------
md <- M1_clean@meta.data
hm_df <- md[!is.na(md$cluster_id) & !is.na(md$group), ]
hm_long <- do.call(rbind, lapply(names(score_cols), function(p) {
  x <- aggregate(hm_df[[score_cols[p]]],
                 by = list(cluster_id = hm_df$cluster_id, group = hm_df$group),
                 FUN = mean, na.rm = TRUE)
  x$program <- p; names(x)[3] <- 'mean_score'; x
}))
hm_long <- do.call(rbind, lapply(split(hm_long, hm_long$program), function(d) {
  d$z <- as.numeric(scale(d$mean_score)); d
}))
hm_long$program <- factor(hm_long$program,
                           levels = c('GR_targets','MHCII',
                                      'Inflammasome','TypeI_IFN'),
                           labels = c('GR targets (down)','MHC-II (up)',
                                      'Inflammasome (neg)','Type-I IFN (neg)'))
hm_long$cluster_id <- factor(hm_long$cluster_id, levels = .myeloid_levels)
hm_long$group      <- factor(hm_long$group, levels = c('NF','pRV','RVF'))

pdf('./output/Myeloid_programs_heatmap.pdf', width = 11, height = 3.5)
print(ggplot(hm_long, aes(x = group, y = cluster_id, fill = z)) +
      geom_tile(color = 'white') +
      scale_fill_gradient2(low = '#2166AC', mid = 'grey95', high = '#B2182B',
                           midpoint = 0, name = 'z(mean score)') +
      facet_wrap(~ program, nrow = 1) +
      xlab('') + ylab('') +
      theme_grey() +
      theme(strip.text = element_text(size = 9, face = 'bold')))
dev.off()
cat('Wrote: Myeloid_programs_heatmap.pdf\n')

# 6G: pooled-myeloid pseudobulk box plots -----------------------------------
# Each point = one patient's mean module score, pooled across all myeloid.
# Program name on the y-axis; title = NF vs disease (pRV+RVF) Wilcoxon p.
md$patient <- as.character(md$patient)
md$grp <- factor(md$group, levels = c('NF','pRV','RVF'))
.grp_cols <- c(NF = '#888888', pRV = '#3B82F6', RVF = '#DC2626')

.pb_pooled <- function(scorecol) {
  d <- md[!is.na(md$group), ]
  a <- aggregate(d[[scorecol]], list(patient = d$patient, grp = d$grp), mean)
  setNames(a, c('patient','grp','score'))
}
.pval_nf_dis <- function(df) {
  nf <- df$score[df$grp == 'NF']; dis <- df$score[df$grp %in% c('pRV','RVF')]
  if (length(nf) >= 2 && length(dis) >= 2)
    suppressWarnings(wilcox.test(nf, dis)$p.value) else NA
}
.g6_defs <- list(list(s = score_cols[['GR_targets']],   y = 'GR targets'),
                 list(s = score_cols[['MHCII']],        y = 'MHC-II'),
                 list(s = score_cols[['Inflammasome']], y = 'Inflammasome'),
                 list(s = score_cols[['TypeI_IFN']],    y = 'Type-I IFN'))
.make_6g <- function(d) {
  df <- .pb_pooled(d$s); pp <- .pval_nf_dis(df)
  ggplot(df, aes(x = grp, y = score, fill = grp)) +
    geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0.5) +
    geom_jitter(width = 0.12, size = 2.4, shape = 21, color = 'black', stroke = 0.4) +
    scale_fill_manual(values = .grp_cols) +
    scale_x_discrete(expand = expansion(add = 0.45)) +
    xlab('') + ylab(d$y) + ggtitle(sprintf('NF vs disease  p=%.3f', pp)) +
    theme_classic(base_size = 13.75) +
    theme(legend.position = 'none',
          plot.title   = element_text(size = 12.5, face = 'bold', color = 'black'),
          axis.text.x  = element_text(color = 'black'),
          axis.text.y  = element_text(color = 'black'),
          axis.title.y = element_text(color = 'black', face = 'bold'),
          plot.margin  = margin(1.5, 4, 1.5, 4))
}
g6_panels <- lapply(.g6_defs, .make_6g)
.gh <- 2.36 * 0.75 * length(g6_panels)
pdf('./output/Myeloid_programs_violins_home.pdf', width = 2.92, height = .gh)
print(patchwork::wrap_plots(g6_panels, ncol = 1)); dev.off()
cat('Wrote: Myeloid_programs_violins_home.pdf (6G, 4 pooled-myeloid box panels)\n')

# 6H: gene-level companion dot plots (fixed/non-z-scored color scale) --------
# Color = absolute log1p mean expression (scale = FALSE) with shared limits, so a
# truly flat gene keeps the same color across NF/pRV/RVF (negatives read flat).
.s2f_theme <- function() {
  theme_classic(base_size = 13.75) +        # fonts +25%
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border     = element_rect(linewidth = 1, fill = NA, color = 'black'),
          axis.line.x      = element_blank(),
          axis.line.y      = element_blank(),
          axis.text.x      = element_text(angle = 45, hjust = 1, vjust = 1, size = 10.6),
          axis.text.y      = element_text(size = 11.25),
          plot.title       = element_text(size = 12.5, face = 'bold'),
          plot.margin      = margin(2, 4, 2, 4))
}
Idents(M1_clean) <- factor(as.character(M1_clean$group), levels = c('NF','pRV','RVF'))
.h6_blocks <- list('GR-transactivation targets (down)'     = prog_list[['GR_targets']],
                   'MHC-II / antigen presentation (up)'     = prog_list[['MHCII']],
                   'Inflammasome (pertinent negative)'      = prog_list[['Inflammasome']],
                   'Type-I interferon (pertinent negative)' = prog_list[['TypeI_IFN']])
.make_6h <- function(nm, feats)
  DotPlot(M1_clean, features = feats, dot.scale = 5, scale = FALSE) +
    scale_color_gradient(low = 'white', high = '#B2182B', limits = c(0, 3),
                         oob = scales::squish, name = 'avg expr (log1p)') +
    scale_size(limits = c(0, 100), range = c(0.5, 5), name = '% expressed') +
    ggtitle(nm) + xlab('') + ylab('') + .s2f_theme()
h6 <- mapply(.make_6h, names(.h6_blocks), .h6_blocks, SIMPLIFY = FALSE)
composite_dot <- patchwork::wrap_plots(h6, ncol = 1) +
  patchwork::plot_layout(guides = 'collect') &
  guides(color = guide_colorbar(direction = 'horizontal', title.position = 'top',
                                barwidth = unit(2.5, 'cm'), barheight = unit(0.35, 'cm')),
         size  = guide_legend(direction = 'horizontal', title.position = 'top', nrow = 1)) &
  theme(legend.position = 'bottom', legend.box = 'horizontal',
        legend.title = element_text(size = 8), legend.text = element_text(size = 7))
pdf('./output/Myeloid_programs_dotplots_composite.pdf', width = 3.6, height = 8.6)
print(composite_dot); dev.off()
cat('Wrote: Myeloid_programs_dotplots_composite.pdf (6H, 4 program blocks)\n')

saveRDS(M1_clean@meta.data, './output/fig6_program_scores_meta.rds')
cat('Saved score-augmented metadata: fig6_program_scores_meta.rds\n')
Idents(M1_clean) <- 'cluster_id'


#######################################
#############  FIGURE 6I  #############
####### Xenium spatial myeloid ########
####### (full dataset, 625,305 cells) #
####### CCR2- rMacs + C1q combined;   #
####### non-myeloid in light grey     #
#######################################
## v54 Panel I: spatial localization of myeloid subtypes in the full Xenium
## object. Loads the clean Xenium object, combines `Macrophage_Resident`
## and `Macrophage_C1q` into a single CCR2- rMac category, leaves the other
## myeloid subtypes intact, and renders all non-myeloid cells in light grey.
##
## Inputs:
##   - ~/Documents/XeniumWorkflow/functions/output/Xenium/xenium_obj_subclustered.rds
##     (or a downstream object with `myeloid_subtype` metadata)
##
## Output:
##   - Figure_6_panel_I_xenium_spatial_myeloid.pdf
##
## NOTE: The heavy myeloid subclustering and metadata curation lives in
## additional_scripts/XeniumReanalysis.r ("Macrophage_C1q" / "Macrophage_Resident"
## sections). This block only reads the curated object and produces the panel.
## Wrapped in file.exists() so the script parses and runs even when the Xenium
## object is unavailable (placeholder PDF emitted instead).

.xen_path <- './dependencies/shared/Xenium/xenium_obj_subclustered.rds'
if (file.exists(.xen_path)) {
  message('Figure 6I: loading Xenium object from ', .xen_path)
  xen <- readRDS(.xen_path)

  ## Build the combined myeloid label (CCR2- rMac merges Resident + C1q).
  ## Robustly detect the myeloid-subtype column (schema varies by object
  ## build); fall back across known candidates and report what was found.
  meta <- xen@meta.data
  ## Fine myeloid subtypes live in cell_types_subclustering (built by
  ## additional_scripts/XeniumReanalysis.r).  cell_type_rctd_doublet is
  ## BROAD (CM/EC/FB/Myeloid/...) — using it collapses all myeloid into one
  ## blob, so prefer the subclustering column.
  .lbl_candidates <- c('cell_types_subclustering', 'cell_types_manual',
                        'myeloid_subtype', 'cell_type_seurat', 'Subsubnames',
                        'Subnames', 'cluster_id', 'cell_subtype',
                        'cell_type_rctd_doublet', 'cell_type')
  .lbl_col <- .lbl_candidates[.lbl_candidates %in% colnames(meta)][1]
  if (is.na(.lbl_col)) {
    meta$.mye_raw <- NA_character_
    message('Figure 6I: no known myeloid-subtype column found; columns = ',
            paste(head(colnames(meta), 40), collapse = ', '))
  } else {
    meta$.mye_raw <- as.character(meta[[.lbl_col]])
    message('Figure 6I: using myeloid label column "', .lbl_col, '"')
  }
  ## Normalise the spellings XeniumReanalysis.r uses; Resident(+C1q) -> CCR2-
  ## rMac, Macrophage_TREM2 -> TREM2_Mac.  Anything NOT in the myeloid
  ## allowlist becomes Non-myeloid (light grey) so CM/EC/FB/PC/VSMC/... grey.
  .norm_mye <- function(x) {
    x <- as.character(x)
    out <- rep('Non-myeloid', length(x))
    out[x %in% c('Macrophage_Resident', 'Macrophage_C1q',
                 'Macrophage_Resident_LYVE1')] <- 'CCR2- rMac (Resident + C1q)'
    out[x %in% c('Macrophage_Inflammatory', 'iMac')] <- 'Macrophage_Inflammatory'
    out[x %in% c('Monocyte', 'Mono')]                 <- 'Monocyte'
    out[x %in% c('Macrophage_TREM2', 'TREM2_Mac', 'TREM2 Mac')] <- 'TREM2_Mac'
    out[x %in% c('Dendritic_Cell', 'DCs')]            <- 'Dendritic_Cell'
    out[x %in% c('Myeloid_Proliferating')]            <- 'Myeloid_Proliferating'
    out
  }
  myeloid_lbl <- .norm_mye(meta$.mye_raw)
  message('Figure 6I: myeloid_panel_I composition:')
  print(table(myeloid_lbl, useNA = 'ifany'))
  xen$myeloid_panel_I <- factor(
    myeloid_lbl,
    levels = c('CCR2- rMac (Resident + C1q)', 'Macrophage_Inflammatory',
               'Monocyte', 'TREM2_Mac', 'Dendritic_Cell',
               'Myeloid_Proliferating', 'Non-myeloid')
  )

  myeloid_cols <- c(
    'CCR2- rMac (Resident + C1q)' = '#1f77b4',
    'Macrophage_Inflammatory'     = '#d62728',
    'Monocyte'                    = '#9467bd',
    'TREM2_Mac'                   = '#2ca02c',
    'Dendritic_Cell'              = '#ff7f0e',
    'Myeloid_Proliferating'       = '#17becf',
    'Non-myeloid'                 = 'grey85'
  )

  ## Spatial scatter: iterate over FOVs (multi-tissue Xenium) and align coords
  ## with metadata via cell name (not positional). GetTissueCoordinates without
  ## `image=` returns only the active FOV (97,493 cells), so we loop.
  .fov_names <- names(xen@images)
  message('Figure 6I: building coords from ', length(.fov_names), ' FOVs')
  .coords_list <- lapply(.fov_names, function(.fov) {
    cc <- GetTissueCoordinates(xen, image = .fov)
    cell_col <- intersect(colnames(cc), c('cell', 'cells'))[1]
    if (is.na(cell_col)) cell_col <- colnames(cc)[ncol(cc)]
    data.frame(x = cc[, 1], y = cc[, 2], cell = cc[[cell_col]],
               fov = .fov, stringsAsFactors = FALSE)
  })
  coords <- do.call(rbind, .coords_list)
  meta_all <- xen@meta.data
  meta_all$cell <- rownames(meta_all)
  df_xen <- merge(coords, meta_all[, c('cell','myeloid_panel_I','orig.ident')],
                  by = 'cell', all.x = TRUE)
  df_xen$label  <- as.character(df_xen$myeloid_panel_I)
  df_xen$label[is.na(df_xen$label)] <- 'Non-myeloid'
  df_xen$sample <- as.character(df_xen$orig.ident)

  ## Original offset technique (cf. F3 offset_coords): normalise each tissue
  ## (sample) to its own origin, then translate it into a fixed grid cell so
  ## tiles NEVER overlap.  ONE coord_equal panel, NO facet (facet + free
  ## scales is incompatible with coord_equal and left tiles stacked).
  .samps <- sort(unique(df_xen$sample))
  .ncol  <- 3
  .gap   <- 500                                   # µm gap between tiles
  .bb <- lapply(.samps, function(s) {
    ix <- df_xen$sample == s
    list(w = diff(range(df_xen$x[ix])), h = diff(range(df_xen$y[ix])))
  })
  .nrowg <- ceiling(length(.samps) / .ncol)
  .colw  <- rep(0, .ncol); .rowh <- rep(0, .nrowg)
  for (i in seq_along(.samps)) {
    r <- ceiling(i / .ncol); cc <- ((i - 1) %% .ncol) + 1
    .colw[cc] <- max(.colw[cc], .bb[[i]]$w)
    .rowh[r]  <- max(.rowh[r],  .bb[[i]]$h)
  }
  .xstart <- c(0, cumsum(.colw[-.ncol]  + .gap))
  .ystart <- c(0, cumsum(.rowh[-.nrowg] + .gap))
  .lab_df <- data.frame()
  for (i in seq_along(.samps)) {
    s  <- .samps[i]; ix <- df_xen$sample == s
    r  <- ceiling(i / .ncol); cc <- ((i - 1) %% .ncol) + 1
    df_xen$x[ix] <- df_xen$x[ix] - min(df_xen$x[ix]) + .xstart[cc]
    df_xen$y[ix] <- (max(df_xen$y[ix]) - df_xen$y[ix]) +
                    (max(.ystart) - .ystart[r])      # flip → grid top-to-bottom
    .lab_df <- rbind(.lab_df, data.frame(
      sample = s,
      x = .xstart[cc] + .colw[cc] / 2,
      y = (max(.ystart) - .ystart[r]) + .rowh[r] + .gap * 0.4))
  }
  df_xen <- df_xen[order(df_xen$label != 'Non-myeloid'), ]   # myeloid drawn on top

  p_6I <- ggplot(df_xen, aes(x = x, y = y, color = label)) +
    geom_point(size = 0.05, alpha = 0.8) +
    geom_text(data = .lab_df, aes(x = x, y = y, label = sample),
              inherit.aes = FALSE, size = 2.6, fontface = 'bold') +
    scale_color_manual(values = myeloid_cols, name = 'Myeloid subtype') +
    coord_equal() +
    theme_void(base_family = FONT_FAMILY) +
    theme(legend.position  = 'bottom',
          legend.text      = element_text(size = 7),
          plot.background  = element_rect(fill = 'white', color = NA)) +
    guides(color = guide_legend(override.aes = list(size = 2)))

  save_figure(p_6I, 'Figure_6_panel_I_xenium_spatial_myeloid.pdf',
              width = 12, height = 9)
  message('Figure 6I: wrote Figure_6_panel_I_xenium_spatial_myeloid.pdf')
  rm(xen, df_xen, coords, meta, myeloid_lbl); gc()
} else {
  message('Figure 6I: Xenium object not found at ', .xen_path,
          ' — emitting placeholder PDF.')
  p_6I_placeholder <- ggplot() +
    annotate('text', x = 0.5, y = 0.5,
             label = paste0('Figure 6I placeholder\n',
                            'Xenium object not available at:\n', .xen_path),
             size = 4, family = FONT_FAMILY) +
    theme_void()
  save_figure(p_6I_placeholder,
              'Figure_6_panel_I_xenium_spatial_myeloid.pdf',
              width = 8, height = 6)
}

cat('\n=== Figure 6 panels A-I rebuilt (v54 nine-panel scheme) ===\n')

###############################################################################
## v57 standardized per-panel emission — keyed to the v57 F6 legend
## (.figure_run_logs/v57_figure_legends.md): A myeloid UMAP | B marker dot |
## C freq by disease (Kruskal) | D module dot M1/M3/M4/M8 | E ChEA UP
## RVF-vs-NF | F ChEA DOWN RVF-vs-NF | G program violins | H companion
## gene dots | I spatial myeloid. Copy legacy ./output/<name>.pdf to
## V52_FIG_DIR/Figure_6_panel_<L>.pdf.
###############################################################################
.f6_v57 <- list(
  A = 'Myeloid_snUMAP.pdf',
  B = 'Myeloid_dot.pdf',
  C = 'Myeloid_clust_freq_stats.pdf',
  D = 'Myeloid_module_dotplot.pdf',
  E = 'Myeloid_chea_terms_cell_type_up_RVF_vs_NF.pdf',
  F = 'Myeloid_chea_terms_cell_type_down_RVF_vs_NF.pdf',
  G = 'Myeloid_programs_violins_home.pdf',
  H = 'Myeloid_programs_dotplots_composite.pdf',
  I = 'Figure_6_panel_I_xenium_spatial_myeloid.pdf')
.f6_dir <- if (exists('V52_FIG_DIR')) V52_FIG_DIR else './output/Figure_6'
dir.create(.f6_dir, showWarnings = FALSE, recursive = TRUE)
message('Figure 6 v57-keyed standardized panels:')
for (.L in names(.f6_v57)) {
  .src <- NA_character_
  for (.d in c('./output', .f6_dir, './output/Figure_6')) {
    if (file.exists(file.path(.d, .f6_v57[[.L]]))) {
      .src <- file.path(.d, .f6_v57[[.L]]); break }
  }
  if (!is.na(.src)) {
    file.copy(.src, file.path(.f6_dir,
              sprintf('Figure_6_panel_%s.pdf', .L)), overwrite = TRUE)
    message('  ', .L, ': OK (', basename(.src), ')')
  } else message('  ', .L, ': MISSING (', .f6_v57[[.L]], ')')
}
