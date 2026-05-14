###############################################################################
## Supplementary Figure 4 (v54 draft) — Left ventricular failure (DCM) comparison
##
## Panels (from RV_snRNASeq_v54_draft.md figure legends):
##   (A) UMAP of integrated human RV FB lineages (top), reference-mapped LV (bottom)
##   (B) Fibroblast subtype concordance between ventricles
##   (C) Conservation of FB Phase 1 + Phase 2 programs in LV snRNA-seq
##       (Koenig 2022; per-patient pseudobulk; Donor vs DCM Wilcoxon)
##   (D) Scatter log2FC human LV-FB(DCM/NF) vs human RV-FB(RVF/NF)
##   (E) UMAP of integrated human RV myeloid lineages
##   (F) Reference mapping of LV myeloid: snRNA (i), scRNA (ii); original LV
##       annotations (iii) snRNA, (iv) scRNA
##   (G) [TODO] Cross-chamber concordance scatter for myeloid bulk-WGCNA modules
##       (M1, M3, M4, M8): per-gene log2FC RV vs LV myeloid (Koenig 2022)
##   (H) [TODO] Bulk WGCNA myeloid modules (M1/M3/M4/M8) projected onto LV
##       (Koenig 2022) myeloid atlas via Seurat::AddModuleScore — dot plot of
##       mean per-cell module score by LV myeloid subtype × Donor/DCM
##   (I) Per-patient pseudobulk module-score box plots in LV myeloid for the
##       4 RV-derived programs (GR-homeostatic, HIF/vascular drift, NF-kB
##       MHCII/inflammasome, IFNg antigen presentation); Donor vs DCM Wilcoxon
##
## Output: ./output/Supplementary_Figure_4/v52_figures/Supplementary_Figure_4_panel_<LETTER>.pdf (per-panel)
###############################################################################

source('./helper_scripts/_shared_helpers.R')

## Per-figure output directory (introduced for consistent output paths)
V52_FIG_DIR <- './output/Supplementary_Figure_4'
dir.create(V52_FIG_DIR, showWarnings = FALSE, recursive = TRUE)


## Suppress R's default Rplots.pdf in cwd when Rscript hits a plot call
## that's outside an explicit pdf() ... dev.off() envelope.
pdf(NULL)
# SCTransform / FindIntegrationAnchors export ~1 GB of globals to future workers;
# default cap is 500 MiB. Set high enough for the largest pipeline step.
options(future.globals.maxSize = 16 * 1024^3)

COMP_W <- 12
COMP_H <- 16

## Publication geom scales (linewidths, point sizes, dot ranges, bracket widths)
PS <- pub_scales(COMP_W)

library(Seurat)
library(hdWGCNA)
library(ggeasy)
library(harmony)



source('./helper_scripts/spatial_functions.R')

## ── Cartoon / schematic asset insertion ───────────────────────────────────────
## Supplementary_Figure_3.pdf (v51 composite) is a DCM LV cross-comparison
## figure assembled from UMAPs, DotPlots, and violin plots only; the current
## composite contains no explicit cartoon/illustration panel. If a reference-
## mapping schematic (sn/sc LV  RV myeloid projection) is later added, extract
## it from the composite with the magick block below.
##
## # library(magick)
## # pdf_path <- '~/Downloads/hdWGCNA_TOM/Manuscripts/Supplementary_Figure_3.pdf'
## # pg <- image_read_pdf(pdf_path, pages = 1, density = 300)
## # # 17 cm = ~2008 px at 300 DPI; ~2008x800 band for a schematic row.
## # cartoon <- image_crop(pg, geometry = '2008x800+0+0')
## # image_write(cartoon, file.path(ASSET_DIR, 'SuppFig3_lv_mapping.png'),
## #             format = 'png')
## s3_cartoon <- insert_asset('SuppFig3_lv_mapping.png',
##                            label = 'LV -> RV reference mapping schematic')




#######################################
######  FIGURE S4 PANEL E (+ F)  ######
#######################################
# Myeloid pipeline (v54 panels E and F i-iv):
#   E: RV_Myeloid_Integrated_snUMAP.pdf
#   F i: RV_LV_myeloid_ref_mapped.pdf            (LV snRNA RefMapped)
#   F ii: RV_LV_sc_myeloid_ref_mapped.pdf        (LV scRNA RefMapped)
#   F iii/iv: LV_sn_sc_myeloid_ref_mapped.pdf    (LV original annotations)
#   G: RV_vs_LV_myeloid__dot.pdf                 (LV vs RV failure concordance)
#   H: LV_vln_modules.pdf, LV_dot_modules.pdf    (LV WGCNA module DCM vs NF)
# Memory-gated: this block plus the FB SCTransform pipeline together exceed
# ~44 GB on this 48 GB box. Set RUN_LEGACY_MYELOID=TRUE in the environment
# to enable; the master driver (run_S4_stages.sh) runs this block in a
# dedicated Rscript invocation so the heap is freed before FB/box stages.
if (isTRUE(as.logical(Sys.getenv('RUN_LEGACY_MYELOID', 'FALSE')))) {

# Use the cleaned, post-Contamination, Harmony-integrated object built in
# Figure_6.R (line 157). It has the canonical v53 cluster labels as Idents:
#   Mac_Inflammatory, CCR2- Resident Mac, Dendritic Cell,
#   Monocyte / Mac_Mono_Derived, TREM2+ Mac, Proliferating
# Reductions: pca, harmony, umap (already finalized). Skipping the legacy
# FindClusters / RenameIdents / sub-clustering / re-Harmony pipeline.
M1 <- readRDS(file = "./output/Figure_6/myeloid_subclust_new.rds")
M1$Subnames    <- Idents(M1)
M1$Subsubnames <- Idents(M1)
M1 <- SetIdent(M1, value = 'Subsubnames')
M2 <- M1                       # alias AFTER Subsubnames is set, so M2$Subsubnames exists

M1 <- AddModuleScore(M1, features=list(c("TREM2","GPNMB","MITF","SPP1")),assay="SCT",name="TREM2_Mac_Score")
M1 <- AddModuleScore(M1, features=list(c("CLEC9A","ZBTB46","CD1C","CD226")),assay="SCT",name="DC_Score")
M1 <- AddModuleScore(M1, features=list(c("FCN1","LILRB2","ITGAL","CSF3R")),assay="SCT",name="Mono_Score")
M1 <- AddModuleScore(M1, features=list(c("CCR2","CX3CR1","ITGAX")),assay="SCT",name="CCR2+_rMac_Score")
M1 <- AddModuleScore(M1, features=list(c("LYVE1","FOLR2","SIGLEC1","F13A1")),assay="SCT",name="CCR2-_rMac1_Score")
M1 <- AddModuleScore(M1, features=list(c("RBMS3","PLA2G5","EBF1")),assay="SCT",name="CCR2-_rMac2_Score")
M1 <- AddModuleScore(M1, features=list(c("IL1B","CCL3","CCL4","CXCL3","CXCL8")),assay="SCT",name="iMac_Score")


M1$SubNames_group <- paste0(M1$Subsubnames,'_',M1$group)
M1 <- SetIdent(M1, value = "group")
M1$group <- factor(M1$group,levels=c("NF","pRV","RVF"))

consensus_modules <- read.csv("./dependencies/shared/bulk_heart_modules.csv")
consensus_modules <- consensus_modules[,1:3]
consensus_modules <- subset(consensus_modules, gene_name %in% rownames(M1))
# remove duplicate gene names
consensus_modules <- consensus_modules[match(unique(consensus_modules$gene_name), consensus_modules$gene_name),]


library(dplyr)
mapping <- labels2colors(1:100)

score_calc <- consensus_modules %>% group_by(module) %>% group_split()
module_colors <- unique(unlist(lapply(score_calc,'[[','module')))
module_colors <- paste0('M',match(module_colors,mapping))


M1 <- AddModuleScore(M1,lapply(score_calc,'[[','gene_name'),name="module_score")
invisible(gc(verbose = FALSE))

cols_current <- colnames(M1@meta.data)
cols_current[startsWith(colnames(M1@meta.data),'module_score')] <- paste0('module_',module_colors)
colnames(M1@meta.data) <- cols_current

pdf(paste0('./output/Supplementary_Figure_4/', 'RV_vln_modules.pdf'), width=6, height=3)
VlnPlot(M1,c('module_M1','module_M3','module_M4','module_M8'),pt.size=0,ncol=4) &
  scale_fill_disease() & theme_v52(COMP_W)
dev.off()

pdf(paste0('./output/Supplementary_Figure_4/', 'RV_dot_modules.pdf'), width=4.5, height=2.2)
p <- DotPlot(M1,paste0('module_',
  c('M1','M3','M4','M8')),dot.min=0,col.min=0,col.max=2,
  dot.scale = PS$dot_range[2]) +
  RotatedAxis() + ylab('')+ xlab('')+
  scale_color_gradient2(high='red', mid='grey95', low='blue') +
  scale_size_continuous(range = PS$dot_range) +
  theme_v52(COMP_W) +
  theme(
    panel.border = element_rect(linewidth = PS$linewidth_mm, fill=NA, color='black'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()
)
p
dev.off()



# Panel E: harmony-based UMAP saved in myeloid_subclust_new.rds (matches
# Figure_6.R reduction='harmony', dims=1:30). Color palette pinned to .rv_pal
# below so Panel E matches F.i / F.ii left side.
.rv_levels <- levels(droplevels(factor(M2$Subsubnames)))
.rv_pal    <- setNames(scales::hue_pal()(length(.rv_levels)), .rv_levels)
M2$Subsubnames <- factor(M2$Subsubnames, levels = .rv_levels)
pdf(paste0('./output/Supplementary_Figure_4/', 'RV_Myeloid_Integrated_snUMAP.pdf'), width=5, height=5)
print(
  PlotEmbedding(M2, group.by='Subsubnames', point_size=1, plot_under=TRUE,
                plot_theme=umap_theme()+NoLegend(),
                raster_dpi=400, raster_scale=0.5) +
    scale_color_manual(values = .rv_pal, drop = FALSE)
)
dev.off()

# Refresh M2 UMAP with return.model=TRUE so MapQuery's ProjectUMAP works.
# Trained on harmony reduction (matches Fig 6 layout).
M2 <- RunUMAP(M2, reduction = "harmony", dims = 1:30, return.model = TRUE)

#######################################
############  FIGURE S4 PANEL F  ######
####### LV myeloid ref mapping ########
####### (i) snRNA + (ii) scRNA;       #
####### (iii)/(iv) original LV labels #
#######################################
# Process snLV and scLV one at a time so peak memory holds at most one
# Koenig SCT'd object at once. The joint sn+sc plot reuses ggplot objects
# captured before each Seurat object is freed.

# Koenig modules (recreate Mac1-5 / Mono / DC / Prolif via marker-argmax on
# de-novo Louvain clusters — robust to Seurat/igraph version drift).

options(future.globals.maxSize = 50 * 1024^3)

.koenig_modules <- list(
  Nonclassical_Mono  = c('FCGR3A','LILRA5','LST1'),
  Classical_Mono     = c('CD14','S100A8','S100A9','S100A12','FCN1'),
  Intermediate_Mono  = c('FCN1','OLR1','PLAUR','TRAF1'),
  Prolif             = c('MKI67','STMN1','BIRC5','TOP2A'),
  Mac1               = c('TREM2','SPP1','FABP5','LGALS3'),
  Mac2               = c('FOLR2','LYVE1','MRC1','SIGLEC1','CD163'),
  Mac3               = c('LYVE1','HSPH1','HSPA1A','HSPA1B'),
  Mac4               = c('CCL3','CCL4','PHLDA1','PMAIP1'),
  Mac5               = c('KLF2','KLF4','EGR1','RHOB'),
  DC                 = c('CD1C','CCR7','FCER1A'))
.module_cols   <- paste0(names(.koenig_modules), '1')   # AddModuleScore appends '1'
.module_labels <- names(.koenig_modules)

# .rv_pal / .rv_levels were defined for Panel E above. .lv_pal is used by the
# Koenig-labeled F.iii / F.iv sub-panels.
.lv_pal <- setNames(scales::hue_pal()(length(.module_labels)), .module_labels)

# Reference panel of M2 (used in both F.i and F.ii); render once, reuse.
p_M2 <- DimPlot(M2, reduction = "umap", group.by = "Subsubnames",
                label = TRUE, label.size = 3, repel = TRUE,
                raster = TRUE, pt.size = 1.5) +
  scale_color_manual(values = .rv_pal, drop = FALSE) +
  NoLegend() + ggtitle("Reference annotations")

# Helper: process one Koenig query (tech == 'SN' or 'SC'), plot its ref-mapped
# pair (M2 + query), then return a DimPlot of LV_names for the joint sn+sc
# panel. The Seurat object is rm()'d on exit.
.process_LV <- function(tech_id, ref_pdf) {
  q <- readRDS('./dependencies/shared/Kory_myeloid_fibroblasts.rds')
  q <- subset(q, Names == 'Myeloid')   # ~10x smaller — drop other cell types first
  invisible(gc(verbose = FALSE))
  q <- subset(q, tech == tech_id)

  DefaultAssay(q) <- "RNA"
  q[['SCT']]      <- NULL
  q <- SCTransform(q, vst.flavor = "v2", assay = 'RNA')
  q <- RunPCA(q, npcs = 50, verbose = FALSE)

  anchors <- FindTransferAnchors(
    reference = M2, query = q,
    normalization.method = "SCT", recompute.residuals = FALSE,
    reference.reduction = "pca", dims = 1:50
  )
  q <- MapQuery(anchorset = anchors, reference = M2, query = q,
                refdata = list(celltype = "Subsubnames"),
                reference.reduction = "pca", reduction.model = "umap")
  q$map_score <- MappingScore(anchors)
  rm(anchors); invisible(gc(verbose = FALSE))

  # Align the query's predicted-celltype levels with M2's Subsubnames so the
  # color palette is shared across F.i left and right (and across F.i and F.ii).
  q$predicted.celltype <- factor(q$predicted.celltype, levels = .rv_levels)

  # F.i / F.ii: M2 reference + query transferred labels
  p_query <- DimPlot(q, reduction = "ref.umap", group.by = "predicted.celltype",
                     label = TRUE, label.size = 3, repel = TRUE,
                     raster = TRUE, pt.size = 1.5) +
    scale_color_manual(values = .rv_pal, drop = FALSE) +
    NoLegend() + ggtitle("Query transferred labels")
  pdf(paste0('./output/Supplementary_Figure_4/', ref_pdf), width = 10, height = 5)
  print(p_M2 + p_query)
  dev.off()

  # Marker-module scores for argmax-labelling
  for (m in names(.koenig_modules)) {
    q <- AddModuleScore(q, list(.koenig_modules[[m]]), name = m)
  }

  # De-novo cluster + marker-argmax label (robust to Seurat/igraph drift)
  q <- RunPCA(q, npcs = 50, verbose = FALSE)
  q <- FindNeighbors(q, dims = 1:50, verbose = FALSE)
  q <- FindClusters(q,
                    resolution = if (tech_id == 'SN') 0.6 else 1,
                    verbose = FALSE)
  q <- RunUMAP(q, dims = 1:50, verbose = FALSE)

  cm <- aggregate(q@meta.data[, .module_cols],
                  by = list(cluster = as.character(Idents(q))),
                  FUN = mean)
  ids <- setNames(.module_labels[apply(cm[, .module_cols], 1, which.max)],
                  cm$cluster)
  q <- RenameIdents(q, ids)
  q$LV_names <- factor(q@active.ident, levels = .module_labels)
  Idents(q)  <- 'LV_names'

  # F.iii / F.iv panel ggplot — captured before q is freed; share .lv_pal
  # with the other panel so sn and sc are in the same color scheme.
  p_LV <- DimPlot(q, reduction = "ref.umap", label = TRUE, label.size = 3,
                  repel = TRUE, raster = TRUE, pt.size = 1.5) +
    scale_color_manual(values = .lv_pal, drop = FALSE) +
    NoLegend() +
    ggtitle(if (tech_id == 'SN') "LV snRNA-seq (Koenig labels)"
            else                  "LV scRNA-seq (Koenig labels)")
  rm(q); invisible(gc(verbose = FALSE))
  p_LV
}

p_sn <- .process_LV('SN', 'RV_LV_myeloid_ref_mapped.pdf')      # Panel F.i
p_sc <- .process_LV('SC', 'RV_LV_sc_myeloid_ref_mapped.pdf')   # Panel F.ii

# Panels F.iii + F.iv: side-by-side LV original (Koenig) annotations.
pdf(paste0('./output/Supplementary_Figure_4/', 'LV_sn_sc_myeloid_ref_mapped.pdf'), width=10, height=5)
print(p_sn + p_sc)
dev.off()
rm(p_sn, p_sc, p_M2); invisible(gc(verbose = FALSE))

}  # end of legacy S3D + S3E skip block

#######################################
########## FIGURE S4 PANEL I  #########
####### LV myeloid program boxes ######
####### (4 RV-derived programs;       #
####### Donor vs DCM Wilcoxon)        #
#######################################
# Patient-level pseudobulk module-score box plots in LV (Koenig snRNA-seq)
# myeloid cells, DCM vs Donor (NF). Mirrors Fig 6 panels G/H gene sets.
# Single horizontal panel combining all 4 program scores:
#   GR-homeostatic (down in DCM) | HIF/vascular drift (up) |
#   NF-kB MHCII / inflammasome (up) | IFNg antigen presentation (up)

# Self-loading: light NormalizeData on RNA assay (no SCTransform) keeps
# the LV myeloid object well under the 40 GB R cap on this 48 GB box.
snLV <- readRDS('./dependencies/shared/Kory_myeloid_fibroblasts.rds')
snLV <- subset(snLV, Names == 'Myeloid')   # drop ~10x cells before tech split
invisible(gc(verbose = FALSE))
snLV <- subset(snLV, tech == 'SN')
DefaultAssay(snLV) <- 'RNA'
snLV <- NormalizeData(snLV, verbose = FALSE)

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

.in_obj <- function(g, obj) intersect(g, rownames(obj))
snLV <- AddModuleScore(snLV, list(.in_obj(GR_homeostatic, snLV)), name = 'GR_hom')
snLV <- AddModuleScore(snLV, list(.in_obj(HIF_vascular,   snLV)), name = 'HIF_vasc')
snLV <- AddModuleScore(snLV, list(.in_obj(NFkB_MHCII,     snLV)), name = 'NFkB_MHCII')
snLV <- AddModuleScore(snLV, list(.in_obj(IFNg_AP,        snLV)), name = 'IFNg_AP')

snLV$condition <- factor(snLV$condition, levels = c('Donor','DCM'))

# Per-patient mean module score (pseudobulk)
.patient_scores <- aggregate(
  cbind(GR_hom     = snLV$GR_hom1,
        HIF_vasc   = snLV$HIF_vasc1,
        NFkB_MHCII = snLV$NFkB_MHCII1,
        IFNg_AP    = snLV$IFNg_AP1),
  by  = list(patient = as.character(snLV$orig.ident),
             condition = snLV$condition),
  FUN = mean
)
write.csv(.patient_scores,
  './output/Supplementary_Figure_4/fig_s4_panel_FGH_LV_myeloid_module_scores.csv', row.names = FALSE)

.box_plot <- function(df, col, title) {
  pv <- tryCatch(wilcox.test(df[[col]] ~ df$condition)$p.value,
                 error = function(e) NA_real_)
  ggplot(df, aes(x = condition, y = .data[[col]], fill = condition)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = PS$linewidth_mm) +
    geom_jitter(width = 0.12, size = PS$scatter_pt, alpha = 0.85) +
    scale_fill_disease() +
    labs(title = title, x = NULL,
         y = 'Per-patient mean module score',
         subtitle = paste0('Wilcoxon p = ', signif(pv, 2))) +
    theme_v52(COMP_W) +
    theme(legend.position = 'none')
}

# Build all 4 box plots
p_GR    <- .box_plot(.patient_scores, 'GR_hom',     'GR homeostatic (down in DCM)')
p_HIF   <- .box_plot(.patient_scores, 'HIF_vasc',   'HIF / vascular drift (up in DCM)')
p_MHCII <- .box_plot(.patient_scores, 'NFkB_MHCII', 'NF-kB MHCII / inflammasome (up in DCM)')
p_IFNg  <- .box_plot(.patient_scores, 'IFNg_AP',    'IFNg antigen presentation (up in DCM)')

# Concatenate horizontally — share the y-axis label, drop redundant per-plot y titles
.strip_y <- theme(axis.title.y = element_blank())
panel_I <- (p_GR + (p_HIF + .strip_y) +
            (p_MHCII + .strip_y) + (p_IFNg + .strip_y)) +
  plot_layout(nrow = 1)

pdf(paste0('./output/Supplementary_Figure_4/', 'fig_s4_panel_I_LV_myeloid_module_scores_box.pdf'),
  width = 5, height = 3)
print(panel_I)
dev.off()

# Free the LV myeloid pipeline objects before reloading for the FB sections.
rm(list = intersect(ls(),
  c('snLV','scLV','M1','M2','anchors','predictions','score',
    'p1','p2','p_GR','p_HIF','p_MHCII','p_IFNg','.patient_scores',
    'a','b','dataset','dataset_fb','shared','shared_fb','labs','labs_fb')))
invisible(gc(verbose = FALSE))

#######################################
############  FIGURE S4A  #############
####### RV+LV FB UMAP ref-mapped ######
#######################################
snLV <- readRDS('./dependencies/shared/Kory_myeloid_fibroblasts.rds')
snLV <- subset(snLV, Names == 'Fibroblasts')   # drop other cell types first
invisible(gc(verbose = FALSE))
snLV <- subset(snLV, tech == 'SN')

snLV <- SetIdent(snLV, value = "condition")



M1 <- readRDS(file = "./dependencies/shared/fb_subclust.rds")

new.cluster.ids <- c("Fb1","Fb2","Fb3","Fb4","Fb5","Fb6","Fb7")
names(new.cluster.ids) <- levels(M1)
M1 <- RenameIdents(M1, new.cluster.ids)

M1$Subnames <- M1@active.ident
M1$SubNames_Groups <- paste(M1$Subnames,M1$group,sep='_')



M1 <- SetIdent(M1, value = "group")
M1$group <- factor(M1$group,levels=c("NF","pRV","RVF"))


#Reference mao LV to RV

M2 <- SplitObject(M1, split.by = "patient")
M2<-PrepSCTIntegration(M2)
features<-SelectIntegrationFeatures(M2)
M2.anchors<-FindIntegrationAnchors(M2,normalization.method = 'SCT',anchor.features = features, reduction = "rpca")
M2 <- IntegrateData(anchorset = M2.anchors,normalization.method='SCT')

DefaultAssay(M2) <- "integrated"

M2 <- RunPCA(M2, npcs = 50, verbose = FALSE)
M2 <- RunUMAP(M2, reduction = "pca", dims = 1:30)



DefaultAssay(snLV) <- "RNA"
snLV[['SCT']] <- NULL

snLV <- SCTransform(snLV, vst.flavor = "v2",assay='RNA')
snLV <- RunPCA(snLV, npcs = 50, verbose = FALSE)

anchors <- FindTransferAnchors(
  reference = M2,
  query = snLV,
  normalization.method = "SCT",
  recompute.residuals=FALSE,
  reference.reduction = "pca",
  dims = 1:50
)

M2 <- RunUMAP(M2, dims = 1:50, return.model = TRUE)
snLV <- MapQuery(anchorset = anchors, reference = M2, query = snLV,
  refdata = list(celltype = "Subnames"), reference.reduction = "pca", reduction.model = "umap")

score <- MappingScore(anchors)

snLV$map_score <- score

# Pin FB palette so reference (Subnames) and query (predicted.celltype) share
# the same color per Fb cluster, with drop=FALSE for missing levels.
.fb_levels <- levels(droplevels(factor(M2$Subnames)))
.fb_pal    <- setNames(scales::hue_pal()(length(.fb_levels)), .fb_levels)
M2$Subnames          <- factor(M2$Subnames,          levels = .fb_levels)
snLV$predicted.celltype <- factor(snLV$predicted.celltype, levels = .fb_levels)

p1 <- DimPlot(M2, reduction = "umap", group.by = "Subnames", label = TRUE, label.size = 3, repel = TRUE,raster=TRUE,pt.size=1.5) +
  scale_color_manual(values = .fb_pal, drop = FALSE) +
  NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(snLV, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE, label.size = 3, pt.size=1.5,repel = TRUE,raster=TRUE) +
  scale_color_manual(values = .fb_pal, drop = FALSE) +
  NoLegend() + ggtitle("Query transferred labels")
pdf(paste0('./output/Supplementary_Figure_4/', 'RV_LV_fb_ref_mapped.pdf'), width=10, height=5)
print(p1 + p2)
dev.off()




snLV <- AddModuleScore(snLV,list(c('ACSM3','SCN7A','ABCA10','NEGR1','ABCA9')),name='Fb1')
snLV <- AddModuleScore(snLV,list(c('PCOLCE2','IGFBP6','MFAP5','S100A10','FGFBP2')),name='Fb2')
snLV <- AddModuleScore(snLV,list(c('GPX3','APOD','C3','HSPA1A','GLUL')),name='Fb3')
snLV <- AddModuleScore(snLV,list(c('PLA2G2A','RARRES1','IGFBP4','FGF7')),name='Fb4')
snLV <- AddModuleScore(snLV,list(c('ELN','GPC6','FGF14','ITGA1')),name='Fb5')
snLV <- AddModuleScore(snLV,list(c('TNC','FN1','MEOX1')),name='Fb6')
snLV <- AddModuleScore(snLV,list(c('CCL2','THBS1','CYR61','NR4A1')),name='Fb7')
snLV <- AddModuleScore(snLV,list(c('THBS4','AEBP1','POSTN','CLU','COMP')),name='Fb8')
snLV <- AddModuleScore(snLV,list(c('SERPINE1','CYR61','NFATC2','LRRFIP1')),name='Fb9')


snLV <- SetIdent(snLV, value = 'predicted.celltype')


#######################################
############  FIGURE S4B  #############
####### FB subtype concordance ########
#######################################
# RV vs LV FB annotation concordance dot.

pdf(paste0('./output/Supplementary_Figure_4/', 'RV_LV_fb_ref_mapped_dot.pdf'), width=6, height=4)
DotPlot(snLV,c('Fb31','Fb61','Fb91','Fb41','Fb71','Fb21',
  'Fb51','Fb81','Fb11'),group.by = "predicted.celltype",col.min=0,col.max=2,
  dot.scale = PS$dot_range[2]) +
  scale_size_continuous(range = PS$dot_range) +
  theme_v52(COMP_W)
dev.off()


#######################################
############  FIGURE S4C  #############
####### LV FB Phase 1 + Phase 2 #######
####### per-patient conservation ######
#######################################
# LV (Koenig SN FB) Phase-1 + Phase-2 conservation, mirroring the RV
# pseudobulk module-score panels in Supplementary_Figure_3.R (D and F).
#   Phase 1 (FB identity, Koenig 2022 donor signature):
#     GPX3, PID1, TGFBR3, ACSM3, APOD
#     Expected: erodes NF→failing  →  DOWN in DCM vs Donor
#   Phase 2 (matrifibrocyte, Fu 2018 PMID 29664017; main-text + IHC):
#     COMP, CHAD, CILP2  (Fig 10G; Fig 11C human IHC for Comp/Chad)
#     Expected: commits as failure progresses  →  UP in DCM vs Donor
#   (Prior 6-gene panel including CILP/FMOD/CTHRC1 was unsupported by Fu's
#    main text; CTHRC1 belongs to a distinct activated-profibrotic state per
#    Ruiz-Villalba et al. Circulation 2020, PMID 32972203.)

# Drop any cells with non-Donor/DCM condition (no third category exists in
# Koenig but defensive in case of NAs) before aggregating, so every column of
# the pseudobulk matrix has a known group.
.lv_fb_obj <- subset(snLV, subset = condition %in% c('Donor','DCM'))

# Pseudobulk per patient, CPM + log1p
.lv_fb_agg <- AggregateExpression(.lv_fb_obj, assays = 'RNA',
                                  group.by = 'orig.ident',
                                  slot = 'counts',
                                  return.seurat = FALSE)$RNA
.lv_lib    <- colSums(.lv_fb_agg)
.lv_logcpm <- log1p(sweep(.lv_fb_agg, 2, .lv_lib, '/') * 1e6)

# Normalize both join keys: AggregateExpression coerces orig.ident through
# make.names() (so 'TWCM-11-78' becomes 'TWCM.11.78'). Mirror that on the
# meta-data side or the join misses every patient with a dash in their ID.
.lv_pat_grp <- .lv_fb_obj@meta.data %>%
  dplyr::distinct(orig.ident, condition) %>%
  dplyr::mutate(orig.ident_key = make.names(as.character(orig.ident)),
                condition      = factor(condition, levels = c('Donor','DCM')))

cat('LV FB pseudobulk patients per condition:\n')
print(table(.lv_pat_grp$condition))

.lv_fb_pb_module <- function(gene_set, label_short) {
  have <- intersect(gene_set, rownames(.lv_logcpm))
  cat(sprintf('%s — detected: %s\n', label_short, paste(have, collapse = ', ')))
  mat <- .lv_logcpm[have, , drop = FALSE]
  z   <- t(scale(t(mat)))
  d <- data.frame(orig.ident_key = make.names(colnames(z)),
                  score          = colMeans(z, na.rm = TRUE),
                  stringsAsFactors = FALSE)
  out <- dplyr::inner_join(d, .lv_pat_grp, by = 'orig.ident_key') %>%
    dplyr::filter(!is.na(condition)) %>%
    dplyr::arrange(condition, orig.ident)
  if (nrow(out) == 0) {
    stop(sprintf('No matching pseudobulk samples for %s. ',
                 label_short),
         'colnames(z) head: ', paste(head(colnames(z), 3), collapse = ', '),
         '; pat key head: ', paste(head(.lv_pat_grp$orig.ident_key, 3),
                                   collapse = ', '))
  }
  out
}

fb_identity_genes <- c('GPX3','PID1','TGFBR3','ACSM3','APOD')
fb_matrifib_genes <- c('COMP','CHAD','CILP2')
lv_fb_phase1 <- .lv_fb_pb_module(fb_identity_genes, 'LV FB identity (Koenig 2022, Phase 1)')
lv_fb_phase2 <- .lv_fb_pb_module(fb_matrifib_genes, 'LV FB matrifibrocyte (Fu 2018, Phase 2)')

.lv_fb_box <- function(d, ylab_text, title_text, csv_tag) {
  pv <- tryCatch(
    suppressWarnings(wilcox.test(score ~ condition, data = d, exact = FALSE)$p.value),
    error = function(e) NA_real_)
  write.csv(d,
    sprintf('./output/Supplementary_Figure_4/fig_s4_panel_C_%s_score.csv', csv_tag),
    row.names = FALSE)

  y_max  <- max(d$score, na.rm = TRUE)
  y_step <- (y_max - min(d$score, na.rm = TRUE)) * 0.18
  lab    <- dplyr::case_when(
    is.na(pv)  ~ 'NA',
    pv < 0.001 ~ '***',
    pv < 0.01  ~ '**',
    pv < 0.05  ~ '*',
    TRUE       ~ sprintf('p=%s', signif(pv, 2)))

  ggplot(d, aes(x = condition, y = score, fill = condition)) +
    geom_boxplot(outlier.shape = NA, width = 0.55,
                 linewidth = PS$geom_lw, alpha = 0.85) +
    geom_jitter(width = 0.12, size = 1.4, shape = 21,
                color = 'black', stroke = 0.3) +
    scale_fill_manual(values = c(Donor = '#4477AA', DCM = '#CC3311'),
                      guide = 'none') +
    annotate('segment', x = 1, xend = 2,
             y = y_max + y_step * 1, yend = y_max + y_step * 1,
             linewidth = PS$geom_lw) +
    annotate('text', x = 1.5, y = y_max + y_step * 1.2, label = lab,
             size = PS$text_mm, family = FONT_FAMILY) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.30))) +
    coord_flip() +
    theme_v52(COMP_W) +
    theme(plot.title = element_text(face = 'bold')) +
    xlab(NULL) + ylab(ylab_text) +
    labs(title = title_text)
}

p_C_phase1 <- .lv_fb_box(lv_fb_phase1,
  ylab_text  = 'FB-identity score (z, pseudobulk)',
  title_text = 'Phase 1 conservation\n(Koenig 2022 donor signature)',
  csv_tag    = 'fb_identity_phase1')
p_C_phase2 <- .lv_fb_box(lv_fb_phase2,
  ylab_text  = 'Matrifibrocyte score (z, pseudobulk)',
  title_text = 'Phase 2 conservation\n(Fu 2018 matrifibrocyte)',
  csv_tag    = 'fb_matrifib_phase2')

pdf(paste0('./output/Supplementary_Figure_4/', 'fig_s4_panel_C_LV_FB_phase_scores_box.pdf'),
  width = 3.5, height = 4.5)
print(p_C_phase1 / p_C_phase2)   # stacked since panels are now horizontal
dev.off()


#######################################
############  FIGURE S4D  #############
####### LV vs RV FB log2FC scatter ####
#######################################
# Transcriptomic concordance in LV vs RV fibroblasts.
# Computes DCM-vs-Donor (LV) and RVF-vs-NF (RV) FB DEG sets, then plots two
# scatters: legacy unlabelled (RV_vs_LV_fb_dot.pdf) and the v53 polished
# version with axis labels and significant-gene highlighting.

snLV <- readRDS('./dependencies/shared/Kory_myeloid_fibroblasts.rds')
snLV <- subset(snLV, Names == 'Fibroblasts')   # drop other cell types first
invisible(gc(verbose = FALSE))
snLV <- subset(snLV, tech == 'SN')


snLV <- SetIdent(snLV,value='condition')

a<-FindMarkers(snLV,ident.1='DCM',ident.2='Donor',recorrect_umi=F)

M1 <- SetIdent(M1,value='group')


b<-FindMarkers(M1,ident.1='RVF',ident.2='NF')

shared <- intersect(rownames(a),rownames(b))
dataset <- data.frame(RV=b[shared,]$avg_log2FC,LV=a[shared,]$avg_log2FC)
rownames(dataset) <- shared
labs <- rownames(dataset)
#labs[abs(dataset$PAB - dataset$RV)<1] <- NA


pdf(paste0('./output/Supplementary_Figure_4/', 'RV_vs_LV_fb_dot.pdf'), width=6, height=8)
ggplot(dataset, aes(x = RV, y=LV)) + geom_point(size = PS$scatter_pt) +
  geom_text_repel(label=labs,max.overlaps=10,
    size = PS$text_mm, family = FONT_FAMILY, fontface = "italic") + theme_v52(COMP_W)
dev.off()

# --- v53 polished version (axis labels + sig highlighting) ---
# Reuses `a` (LV DCM vs Donor, snLV) and `b` (RV RVF vs NF, M1) above.

shared_fb <- intersect(rownames(a), rownames(b))
dataset_fb <- data.frame(
  RV       = b[shared_fb, ]$avg_log2FC,
  LV       = a[shared_fb, ]$avg_log2FC,
  RV_padj  = b[shared_fb, ]$p_val_adj,
  LV_padj  = a[shared_fb, ]$p_val_adj,
  row.names = shared_fb
)
dataset_fb$sig_both <- dataset_fb$RV_padj < 0.05 & dataset_fb$LV_padj < 0.05
labs_fb <- ifelse(dataset_fb$sig_both, rownames(dataset_fb), NA_character_)

pdf(paste0('./output/Supplementary_Figure_4/', 'fig_s4_panel_D_fb_LV_RV_scatter.pdf'), width = 4, height = 4)
print(
  ggplot(dataset_fb, aes(x = RV, y = LV)) +
    geom_hline(yintercept = 0, linewidth = PS$linewidth_mm, colour = 'grey80') +
    geom_vline(xintercept = 0, linewidth = PS$linewidth_mm, colour = 'grey80') +
    geom_point(aes(colour = sig_both), size = PS$scatter_pt) +
    scale_colour_manual(values = c(`TRUE` = 'firebrick', `FALSE` = 'grey60'), guide = 'none') +
    geom_text_repel(aes(label = labs_fb), max.overlaps = 20,
      size = PS$text_mm, family = FONT_FAMILY, fontface = 'italic') +
    labs(x = expression(RV~(RVF~vs~NF~Log[2]~FC)),
         y = expression(LV~(DCM~vs~NF~Log[2]~FC))) +
    theme_v52(COMP_W)
)
dev.off()

write.csv(dataset_fb,
  './output/Supplementary_Figure_4/fig_s4_panel_D_fb_LV_RV_scatter.csv', row.names = TRUE)


#######################################
