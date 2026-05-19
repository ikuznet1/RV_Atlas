###############################################################################
## Supplementary Figure 6 (v54 draft) — PAB cardiomyocyte comparison
##
## Panels (assembled order in Supplementary_Figure_6.png/.ai):
##   (A) Co-embedding UMAPs: human RV CM clusters and PAB transferred labels
##   (B) Mouse-mapped-cluster × human-RV-CM signature dot plot
##   (C) Mapping-score violin per mouse mapped cluster
##   (D) Cross-species FC scatter: human RVF/NF vs mouse PAB Sev/Nor
##   (E) WGCNA bulk-derived CM module dot plot (mouse PAB by group)
##   (F) Per-subject MitoCarta module score, mouse PAB + human RV (stacked)
##   (G) Pseudobulk DESeq2 MitoCarta volcano (RVF vs NF, n=4 vs 3, apeglm)
##   (H) Within-mouse FC trajectories:
##         (i)  mouse Sev/Nor vs Mod/Nor — internal-consistency control
##         (ii) mouse Sev/Mod vs Mod/Nor — "mild as protected state"
##   (I) Cardiac fibrosis (Sirius Red % area) by PAB stage
##   (J) Cross-species CM concordance strip (mouse + human pairs):
##         FAO program | Failure / fetal program | NR3C1 (GR axis)
##
## Source: ported from v51 Supplementary_Figure_6.R; updated for v53 with
## new CM reference (cm_subclust_new_new.rds), per-patient pseudobulk for
## the MitoCarta box and volcano, apeglm shrinkage, curated FC-scatter
## labels (CM-expression filtered, |log2FC|>1), and panels (I)/(J).
##
## Output: /SupplementaryFigure_6.pdf (composite) and
## individual panel PDFs in ./output/.
###############################################################################

source('./helper_scripts/_shared_helpers.R')
## Working-repo adaptation: per-figure output dir + standardized names.
V52_FIG_DIR <- './output/Supplementary_Figure_6'
dir.create(V52_FIG_DIR, showWarnings = FALSE, recursive = TRUE)

COMP_W <- 6
COMP_H <- 11

## Publication scaling constants (geom linewidths, point sizes, text mm, etc.)
PS <- pub_scales(COMP_W)

library(Seurat)
library(hdWGCNA)
library(ggeasy)
library(harmony)
library(tximportData)
library(tximport)
library(viridis)

options(future.globals.maxSize = 16 * 1024^3)

source('./helper_scripts/spatial_functions.R')

###############################################################################
## Cartoon / external assets
##
## Supplementary Figure 6 contains no external cartoon; all panels (A-H) are
## computed ggplot/Seurat plots. insert_asset() call retained for parity with
## other figures and will render a dashed "missing asset" box if an override
## cartoon PNG is ever added later.
##
## (If an asset needs to be produced, crop from original PDF:)
##   cd /Users/ikuz/Documents/RV_Atlas/new_scripts/assets
##   magick -density 300 ~/Downloads/hdWGCNA_TOM/Manuscripts/Supplementary_Figure_6.pdf \
##       -crop WIDTHxHEIGHT+XOFF+YOFF +repage Supplementary_Figure_6_cartoon.png
##
# p_S6_cartoon <- insert_asset('Supplementary_Figure_6_cartoon.png')


#######################################
#############  FIGURE S6A  ############
#######################################

M1 <- readRDS('./dependencies/shared/PAB_data_clean.rds')

M1 <- SetIdent(M1, value = "Names")
M1$group <- M1$orig.ident


M2 <- subset(M1, Names %in% c("CM"))

M2 <- RunPCA(M2)
M2 <- RunHarmony(M2,'patient')
M2 <- FindNeighbors(M2, dims = 1:50,reduction = "harmony")
M2 <- FindClusters(M2, resolution = 0.5,reduction = "harmony")
M2 <- RunUMAP(M2, dims = 1:50,reduction = "harmony")

M2$Names <- M2@active.ident
markers<-FindAllMarkers(M2,recorrect_umi=F)

pdf(paste0('./output/', 'PAB_CM_snUMAP.pdf'), width=5, height=5)
PlotEmbedding(M2,group.by='Names',point_size=1,plot_under=TRUE,plot_theme=umap_theme()+NoLegend(),raster_dpi=400,raster_scale=0.5)
dev.off()


M3 <- readRDS(file = "./dependencies/shared/cm_subclust_new_new.rds")

human2mouse <- read.csv('./dependencies/shared/human2mouse.csv',header=F)
idx <- match(unique(human2mouse[,2]),human2mouse[,2])
human2mouse<-human2mouse[idx,]
colnames(human2mouse) <-c('human_name', 'mouse_name')

RNA <- M2@assays$RNA['counts']
newnames <- human2mouse$human_name[match(rownames(RNA),human2mouse$mouse_name)]
newnames[is.na(newnames)] <- rownames(RNA)[is.na(newnames)]
rownames(RNA) <- newnames


M4 <- CreateAssayObject(RNA)
M2[['humanized']] <- M4
DefaultAssay(M2) <- "humanized"
M2[['SCT']] <- NULL

M2 <- SCTransform(M2, vst.flavor = "v2",assay='humanized')
M2 <- RunPCA(M2, npcs = 50, verbose = FALSE)

anchors <- FindTransferAnchors(
  reference = M3,
  query = M2,
  normalization.method = "SCT",
  recompute.residuals=FALSE,
  reference.reduction = "harmony",
  dims = 1:30
)

predictions <- TransferData(anchorset = anchors, refdata = M3$Subnames, dims = 1:30)


M3 <- RunUMAP(M3, dims = 1:30, reduction = 'harmony', return.model = TRUE)
M2 <- MapQuery(anchorset = anchors, reference = M3, query = M2,
  refdata = list(celltype = "Subnames"), reference.reduction = "pca", reduction.model = "umap")

score <- MappingScore(anchors,ndim = 30)

M2$map_score <- score

p1 <- DimPlot(M3, reduction = "umap", group.by = "Subnames", label = TRUE, label.size = 3, repel = TRUE,raster=TRUE,pt.size=1.5) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(M2, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE, label.size = 3, pt.size=1.5,repel = TRUE,raster=TRUE) + NoLegend() + ggtitle("Query transferred labels")
pdf(paste0('./output/', 'RV_PAB_CM_ref_mapped.pdf'), width=10, height=5)
p1 + p2
dev.off()


pdf(paste0('./output/', 'PAB_CM_ref_mapped.pdf'), width=5, height=5)
PlotEmbedding(M2,group.by='predicted.celltype',reduction = "ref.umap",point_size=0.2,plot_under=TRUE,plot_theme=umap_theme()+NoLegend(),raster_dpi=400,raster_scale=0.5)
dev.off()


pdf(paste0('./output/', 'RV_CM_ref_mapped.pdf'), width=5, height=5)
PlotEmbedding(M3,group.by='Subnames',point_size=0.2,plot_under=TRUE,plot_theme=umap_theme()+NoLegend(),raster_dpi=400,raster_scale=0.5)
dev.off()


FeaturePlot(M2,'map_score',reduction = "ref.umap")


#######################################
#############  FIGURE S6B  ############
#######################################

RV_marks <- FindAllMarkers(M3)

RV_marks_sig <- subset(RV_marks, p_val_adj<0.05 & avg_log2FC>0)

idx<-match(RV_marks_sig$cluster,unique(RV_marks_sig$cluster))

marks <- split(RV_marks_sig$gene,idx)

## Human CM cluster labels in the same order as RV_marks1..RV_marksN
cm_sig_levels <- as.character(unique(RV_marks_sig$cluster))

M2 <- AddModuleScore(M2,marks,name='RV_marks',ctrl=25)
DefaultAssay(M3) <- "SCT"
M3 <- AddModuleScore(M3,marks,name='RV_marks')

M2$predicted.celltype <- factor(M2$predicted.celltype,levels=levels(M3))

pdf(paste0('./output/', 'PAB_CM_marker_scores.pdf'), width=6, height=4)

DotPlot(M2,c('RV_marks1','RV_marks2','RV_marks3','RV_marks4','RV_marks5'
  ,'RV_marks6','RV_marks7','RV_marks8','RV_marks9','RV_marks10'),
  col.min = 0, scale.min = 0, scale.max = 100,group.by='predicted.celltype') +
scale_x_discrete(labels=cm_sig_levels) +
scale_size(range = PS$dot_range) + theme_v52(COMP_W) +
labs(x = 'Human RV CM signature', y = 'Murine mapped cluster') +
theme(legend.key.size = PS$legend_key,
      axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


DotPlot(M3,c('RV_marks1','RV_marks2','RV_marks3','RV_marks4','RV_marks5'
  ,'RV_marks6','RV_marks7','RV_marks8','RV_marks9','RV_marks10'),
  col.min = 0, scale.min = 0, scale.max = 100) +
scale_x_discrete(labels=cm_sig_levels) +
scale_size(range = PS$dot_range) + theme_v52(COMP_W) +
theme(legend.key.size = PS$legend_key,
      axis.text.x = element_text(angle = 45, hjust = 1))


#######################################
#############  FIGURE S6C #############
#######################################

pdf(paste0('./output/', 'PAB_CM_scores.pdf'), width=3.75, height=4.25)
VlnPlot(M2,'map_score',group.by='predicted.celltype',pt.size=0) +
  theme_v52(COMP_W) +
  theme(legend.key.size = PS$legend_key) +
  NoLegend()
dev.off()


#######################################
#############  FIGURE S6D  ############
#######################################


M3 <- SetIdent(M3,value = 'group')
M2 <- SetIdent(M2,value = 'group')


## Curated cardiac/heart-failure gene set used to label FC scatter plots.
## Drawn from sarcomere/cytoskeleton, mitochondrial/FAO/OXPHOS, glycolysis,
## Ca-handling, ECM, stress/apoptosis, hypertrophy/lineage TFs, oxidative
## stress, ion channels, plus the 10 CM_* subtype markers used in S2.
## Both legacy and modern HUGO aliases are listed (e.g. ATP5A1 + ATP5F1A,
## CTGF + CCN2); the runtime CM-expression filter below keeps whichever
## name is actually present and meaningfully expressed in CMs.
.hf_label_genes_raw <- unique(c(
  "TTN","MYH6","MYH7","MYBPC3","ACTC1","ACTA1","ACTA2","TNNT2","TNNI3",
  "TNNC1","TPM1","MYL2","MYL3","MYL7","NEBL","NEB","XIRP1","XIRP2",
  "ANKRD1","MYOM1","MYOM2","NRAP","BAG3","FLNC","OBSCN","JPH2","CNN1",
  "DES","SYNPO2L","CSRP3","MYOZ2","TCAP",
  "NPPA","NPPB","CCDC141",
  "PPARGC1A","PPARGC1B","PPARA","PPARG","NRF1","TFAM","MFN1","MFN2","OPA1",
  "ESRRA","ESRRG","ATP5A1","ATP5F1A","ATP5B","ATP5F1B","COX5B","COX6A2",
  "CKM","CKMT2","CPT1A","CPT1B","CPT2","ACADM","ACADL","ACADVL",
  "HADHA","HADHB","ECH1","ECHS1","ACOX1","HMGCS2","ETFDH","SLC25A20",
  "PDK4","PDHA1","SLC2A1","SLC2A4","HK2","LDHA","LDHB","ALDOA","ENO1",
  "PFKM","PKM",
  "RYR2","ATP2A2","PLN","CASQ2","CALR","CACNA1C","CACNB2","CAMK2D",
  "S100A1","SRL","CALM1",
  "COL1A1","COL3A1","COL4A1","FN1","POSTN","BGN","DCN","FBN1","TGFB1",
  "TGFBR1","TGFBR2","CTGF","CCN2","CYR61","CCN1","SMAD3",
  "HSPA1A","HSPA8","HSPB1","HSPB7","CRYAB","FKBP5","NR3C1","FOXO3",
  "BAX","BCL2","MYC","TRIM63","FBXO32","KLF15","KLF6",
  "GATA4","GATA6","MEF2A","MEF2C","MEF2D","NKX2-5","TBX5","HAND1","HAND2",
  "STAT3","TEAD1","YAP1","CREB1",
  "SOD1","SOD2","CAT","GPX1","GPX4","NFE2L2","NOX4","TXN1","TXN","TXNIP",
  "SLC8A1","KCNQ1","KCNH2","KCNJ2","SCN5A","HCN4","TLN1",
  "CORIN","PCSK6","NCALD","RORA","SPOCK1","EDNRA","TENM3","HTR4",
  "CACNA1E","ARNTL","SOX5","KCNJ3","HS6ST3","CAMK1D","NRXN3",
  "AGT","AGTR1","ACE","NPR3","NRG1","ERBB2","ERBB4",
  "LPL","FABP3","FABP4","CD36",
  "GDF15","TIMP1","MMP2","MMP9","TNF","IL6"
))

## Filter the curated pool to genes actually expressed in human RV CMs.
## Threshold: ≥ 5% of cells with non-zero counts (any cluster). Removes
## fibrosis/inflammation/RAS genes that aren't really CM-transcripts and
## prevents labeling them in the FC scatters.
.cm_pct_for_pool <- {
  .raw_in <- intersect(.hf_label_genes_raw, rownames(M3))
  .cm_counts <- GetAssayData(M3, assay = 'RNA', layer = 'counts')[.raw_in, , drop = FALSE]
  rowMeans(.cm_counts > 0) * 100
}
hf_label_genes <- names(.cm_pct_for_pool)[.cm_pct_for_pool >= 5]
cat(sprintf('FC-scatter label pool: %d genes after CM-expression filter (>=5%%).\n',
            length(hf_label_genes)))

## Pick up to `n` genes from the curated pool that ALSO have
## max(|FC_x|, |FC_y|) >= min_abs_fc, ranked by that score.
.pick_labels <- function(df, label_pool = hf_label_genes,
                          n = 100, min_abs_fc = 1) {
  in_pool <- rownames(df) %in% label_pool
  scores  <- pmax(abs(df$RV), abs(df$PAB))
  passes  <- in_pool & scores >= min_abs_fc
  if (!any(passes)) return(rep(NA_character_, nrow(df)))
  pool_df <- df[passes, , drop = FALSE]
  scores2 <- pmax(abs(pool_df$RV), abs(pool_df$PAB))
  ord     <- order(scores2, decreasing = TRUE)
  top_genes <- rownames(pool_df)[ord[seq_len(min(n, length(ord)))]]
  out <- rep(NA_character_, nrow(df))
  out[match(top_genes, rownames(df))] <- top_genes
  out
}

.fc_scatter <- function(dataset, xlab, ylab,
                         label_pool = hf_label_genes, n = 100,
                         min_abs_fc = 1) {
  labs <- .pick_labels(dataset, label_pool, n, min_abs_fc)
  ggplot(dataset, aes(x = RV, y = PAB)) +
    geom_hline(yintercept = 0, linewidth = PS$linewidth_mm,
               linetype = "dashed", colour = "grey80") +
    geom_vline(xintercept = 0, linewidth = PS$linewidth_mm,
               linetype = "dashed", colour = "grey80") +
    geom_point(size = PS$scatter_pt, alpha = 0.6) +
    geom_text_repel(label = labs, max.overlaps = Inf,
                    size = 2 * PS$text_mm, family = FONT_FAMILY,
                    fontface = "italic", min.segment.length = 0,
                    segment.size = PS$linewidth_mm,
                    segment.color = "grey60") +
    labs(x = xlab, y = ylab) +
    theme_v52(COMP_W) +
    theme(legend.key.size = PS$legend_key,
          aspect.ratio    = 1,
          axis.title      = element_text(size = 2 * PS$base_pt),
          axis.text       = element_text(size = 2 * PS$base_pt))
}


## Panel D: mouse Sev/Nor vs human RVF/NF (cross-species end-stage).
a <- FindMarkers(M3, ident.1 = 'RVF', ident.2 = 'NF')
b <- FindMarkers(M2, ident.1 = 'Sev', ident.2 = 'Nor')
shared  <- intersect(rownames(a), rownames(b))
dataset <- data.frame(RV  = a[shared, ]$avg_log2FC,
                      PAB = b[shared, ]$avg_log2FC,
                      row.names = shared)
cor(dataset$RV, dataset$PAB)^2
p_d <- .fc_scatter(dataset,
  xlab = 'Human RVF vs NF (avg log2FC)',
  ylab = 'Mouse PAB Sev vs Nor (avg log2FC)')
pdf(paste0('./output/', 'PAB_vs_RV_CM__dot.pdf'), width = 5, height = 5)
print(p_d); dev.off()


## Cross-species early phase: mouse Mod/Nor vs human RVF/NF.
a <- FindMarkers(M3, ident.1 = 'RVF', ident.2 = 'NF')
b <- FindMarkers(M2, ident.1 = 'Mod', ident.2 = 'Nor')
shared  <- intersect(rownames(a), rownames(b))
dataset <- data.frame(RV  = a[shared, ]$avg_log2FC,
                      PAB = b[shared, ]$avg_log2FC,
                      row.names = shared)
cor(dataset$RV, dataset$PAB)^2
p_d_mod <- .fc_scatter(dataset,
  xlab = 'Human RVF vs NF (avg log2FC)',
  ylab = 'Mouse PAB Mod vs Nor (avg log2FC)')
pdf(paste0('./output/', 'PAB_vs_RV_CM_Mod_Nor_dot.pdf'), width = 5, height = 5)
print(p_d_mod); dev.off()


## Late-phase progression matched cross-species:
## human pRV → RVF (compensated → failure) compared with mouse Mod → Sev.
a <- FindMarkers(M3, ident.1 = 'RVF', ident.2 = 'pRV')
b <- FindMarkers(M2, ident.1 = 'Sev', ident.2 = 'Mod')
shared  <- intersect(rownames(a), rownames(b))
dataset <- data.frame(RV  = a[shared, ]$avg_log2FC,
                      PAB = b[shared, ]$avg_log2FC,
                      row.names = shared)
cor(dataset$RV, dataset$PAB)^2
p_late <- .fc_scatter(dataset,
  xlab = 'Human RVF vs pRV (avg log2FC)',
  ylab = 'Mouse PAB Sev vs Mod (avg log2FC)')
pdf(paste0('./output/', 'PAB_vs_RV_CM_Sev_Mod_RVF_pRV_dot.pdf'),
    width = 5, height = 5)
print(p_late); dev.off()


## Panel H, top: mouse Sev/Nor vs Mod/Nor (within-mouse, late vs early).
a <- FindMarkers(M2, ident.1 = 'Sev', ident.2 = 'Nor')
b <- FindMarkers(M2, ident.1 = 'Mod', ident.2 = 'Nor')
shared  <- intersect(rownames(a), rownames(b))
dataset <- data.frame(RV  = a[shared, ]$avg_log2FC,
                      PAB = b[shared, ]$avg_log2FC,
                      row.names = shared)
cor(dataset$RV, dataset$PAB)^2
p_h1 <- .fc_scatter(dataset,
  xlab = 'Mouse PAB Sev vs Nor (avg log2FC)',
  ylab = 'Mouse PAB Mod vs Nor (avg log2FC)')
pdf(paste0('./output/', 'PAB_CM_Sev_Nor_Mod_Nor_dot.pdf'),
    width = 5, height = 5)
print(p_h1); dev.off()


## Panel H, bottom: mouse Sev/Mod vs Mod/Nor ("mild as protected state").
a <- FindMarkers(M2, ident.1 = 'Sev', ident.2 = 'Mod')
b <- FindMarkers(M2, ident.1 = 'Mod', ident.2 = 'Nor')
shared  <- intersect(rownames(a), rownames(b))
dataset <- data.frame(RV  = a[shared, ]$avg_log2FC,
                      PAB = b[shared, ]$avg_log2FC,
                      row.names = shared)
cor(dataset$RV, dataset$PAB)^2
p_h2 <- .fc_scatter(dataset,
  xlab = 'Mouse PAB Sev vs Mod (avg log2FC)',
  ylab = 'Mouse PAB Mod vs Nor (avg log2FC)')
pdf(paste0('./output/', 'PAB_CM_Sev_Mod_Mod_Nor_dot.pdf'),
    width = 5, height = 5)
print(p_h2); dev.off()


## Combined vertical stack: D / H1 / H2 — same horizontal dimension (5") so
## panels share x-extent in the figure layout. Each sub-panel keeps its own
## axis labels via .fc_scatter().
p_combo <- p_d / p_h1 / p_h2
ggsave('./output/PAB_CM_FC_combined.pdf',
       p_combo, width = 5, height = 15)
ggsave(file.path(V52_FIG_DIR, 'PAB_CM_FC_combined.pdf'),
       p_combo, width = 5, height = 15)


#######################################
#############  FIGURE S6E  ############
#######################################

consensus_modules <- read.csv("./dependencies/shared/bulk_heart_modules.csv")
consensus_modules <- consensus_modules[,1:3]
consensus_modules <- subset(consensus_modules, gene_name %in% rownames(M2))
consensus_modules <- consensus_modules[match(unique(consensus_modules$gene_name), consensus_modules$gene_name),]
library(dplyr)
mapping <- labels2colors(1:100)

score_calc <- consensus_modules %>% group_by(module) %>% group_split()
module_colors <- unique(unlist(lapply(score_calc,'[[','module')))
module_colors <- paste0('M',match(module_colors,mapping))


M2 <- AddModuleScore(M2,lapply(score_calc,'[[','gene_name'),name="module_score",ctrl = 50)

cols_current <- colnames(M2@meta.data)
cols_current[startsWith(colnames(M2@meta.data),'module_score')] <- paste0('module_',module_colors)
colnames(M2@meta.data) <- cols_current


M2 <- SetIdent(M2, value = "group")
M2$group <- factor(M2$group,levels=c('Nor','Mod','Sev'))


pdf(paste0('./output/', 'PAB_seurat_dot_CM.pdf'), width=5, height=2.5)

p <- DotPlot(M2,paste0('module_',
  c('M2','M10','M12','M25','M26','M28')),dot.min=0,col.min=0,col.max=2,
  scale.min=0, scale.max=100, group.by='group') +
  RotatedAxis() + ylab('')+ xlab('')+
  scale_color_gradient2(high='red', mid='grey95', low='blue') +
  scale_size(range = PS$dot_range) + theme_v52(COMP_W) +
  theme(
    panel.border = element_rect(linewidth = PS$linewidth_mm, fill=NA, color='black'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    legend.key.size = PS$legend_key
)
p

dev.off()

#######################################
############  FIGURE S6F/G  ###########
#######################################

library("readxl")
library('stringr')

Human.Mito <- read_excel("./dependencies/shared/Human.MitoCarta3.0.xls", sheet = "A Human MitoCarta3.0")
Mouse.Mito <- read_excel("./dependencies/shared/Mouse.MitoCarta3.0.xls", sheet = "A Mouse MitoCarta3.0")


M3 <- AddModuleScore(M3,list(Human.Mito$Symbol),name='mito')
M2 <- AddModuleScore(M2,list(union(Human.Mito$Symbol,Mouse.Mito$Symbol)),name='mito',ctrl = 50)

## Disease palette applied: mouse Nor/Mod/Sev and human NF/pRV/RVF both remap
## to the canonical disease palette.
mhc_pal_mouse <- setNames(disease_pal[c('NF','pRV','RVF')], c('Nor','Mod','Sev'))
mhc_pal_human <- disease_pal[c('NF','pRV','RVF')]

## ── Pseudobulk per-subject mean MitoCarta module score ────────────────────
## Per-cell mito1 (above) is averaged within (subject × group) so that the
## unit of statistical analysis is the patient/animal, not the cell.
## This avoids pseudoreplication when comparing disease groups.
library(ggpubr)

## NB: M2$orig.ident is the group label in this object; per-animal ID is M2$patient.
.mouse_mito_pat <- aggregate(
  list(mito = M2$mito1),
  by  = list(patient = as.character(M2$patient),
             group   = as.character(M2$group)),
  FUN = mean
)
.mouse_mito_pat$group <- factor(.mouse_mito_pat$group,
                                levels = c('Nor','Mod','Sev'))
write.csv(.mouse_mito_pat,
          './output/PAB_mouse_CM_mitocarto_per_animal.csv',
          row.names = FALSE)

.human_mito_pat <- aggregate(
  list(mito = M3$mito1),
  by  = list(patient = as.character(M3$patient),
             group   = as.character(M3$group)),
  FUN = mean
)
.human_mito_pat$group <- factor(.human_mito_pat$group,
                                levels = c('NF','pRV','RVF'))
write.csv(.human_mito_pat,
          './output/RV_human_CM_mitocarto_per_patient.csv',
          row.names = FALSE)

p1 <- ggplot(.mouse_mito_pat, aes(x = group, y = mito, fill = group)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = PS$linewidth_mm) +
  geom_jitter(width = 0.12, size = PS$scatter_pt, alpha = 0.85) +
  scale_fill_manual(values = mhc_pal_mouse) +
  ggpubr::stat_compare_means(
    comparisons = list(c('Nor','Mod'), c('Nor','Sev'), c('Mod','Sev')),
    method      = 'wilcox.test',
    label       = 'p.format',
    size        = PS$text_mm
  ) +
  labs(title = 'Mouse PAB — CM MitoCarta',
       x = NULL, y = 'Per-animal mean score') +
  theme_v52(COMP_W) +
  theme(legend.position = 'none')

p2 <- ggplot(.human_mito_pat, aes(x = group, y = mito, fill = group)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = PS$linewidth_mm) +
  geom_jitter(width = 0.12, size = PS$scatter_pt, alpha = 0.85) +
  scale_fill_manual(values = mhc_pal_human) +
  ggpubr::stat_compare_means(
    comparisons = list(c('NF','pRV'), c('NF','RVF'), c('pRV','RVF')),
    method      = 'wilcox.test',
    label       = 'p.format',
    size        = PS$text_mm
  ) +
  labs(title = 'Human RV — CM MitoCarta',
       x = NULL, y = 'Per-patient mean score') +
  theme_v52(COMP_W) +
  theme(legend.position = 'none')

pdf(paste0('./output/', 'PAB_RV_CM_mitocarto.pdf'), width = 3.2, height = 4.675)
print(p1 / p2)
dev.off()
ggsave(file.path(V52_FIG_DIR, 'PAB_RV_CM_mitocarto.pdf'),
       p1 / p2, width = 3.2, height = 4.675)

anno <- trimws(unlist(lapply(lapply(lapply(Human.Mito$MitoCarta3.0_MitoPathways,str_split,'>'),'[[',1),'[[',1)))

anno[anno == "Small molecule transport | Signaling"] = "Small molecule transport"

## Map MitoCarta3.0 top-level categories to display labels for the volcano
## legend. Avoid building a named vector with an empty key (R's c("" = ...)
## errors as 'attempt to use zero-length variable name'); use explicit
## boolean assignment.
.anno_remap <- function(x) {
  x   <- as.character(x)
  out <- character(length(x))
  out[is.na(x) | x == ""]                                           <- "No annotation"
  out[x == "Metabolism"]                                            <- "Metabolism"
  out[x == "OXPHOS"]                                                <- "Oxidative phos"
  out[x == "Protein homeostasis"]                                   <- "Protein import/sorting"
  out[x == "Mitochondrial central dogma"]                           <- "Mito central dogma"
  out[x == "Mitochondrial dynamics and surveillance"]               <- "Mito dynamics"
  out[x == "Signaling"]                                             <- "Signaling"
  out[x == "Small molecule transport"]                              <- "Transport"
  out[out == ""]                                                    <- "No annotation"
  out
}
anno <- .anno_remap(anno)


## Per-patient pseudobulk DESeq2 for the MitoCarta volcano: aggregate human RV
## CMs by patient, then test RVF vs NF with DESeq2 (n = 4 NF + 3 RVF).
## Replaces the prior single-cell Wilcoxon FindMarkers volcano (which inflated
## p-values via cell-level pseudoreplication).
library(DESeq2)

pb <- AggregateExpression(M3, group.by = 'patient',
                          assays   = 'RNA',
                          slot     = 'counts',
                          return.seurat = FALSE)
pb_mat <- round(as.matrix(pb$RNA))
## Seurat prepends 'g' to group names that start with a digit; strip so the
## column names match the raw patient IDs in M3@meta.data$patient.
colnames(pb_mat) <- sub('^g', '', colnames(pb_mat))

pat_meta <- unique(M3@meta.data[, c('patient', 'group')])
pat_meta$patient <- as.character(pat_meta$patient)
pat_meta <- pat_meta[match(colnames(pb_mat), pat_meta$patient), ]
rownames(pat_meta) <- pat_meta$patient
pat_meta$group <- factor(pat_meta$group, levels = c('NF', 'pRV', 'RVF'))

dds <- DESeqDataSetFromMatrix(countData = pb_mat,
                              colData   = pat_meta,
                              design    = ~ group)
dds <- DESeq(dds)
## Apply apeglm log2FC shrinkage (recommended for small-n pseudobulk):
## shrinks noisy high-FC estimates from low-count genes toward 0. Preserve
## padj from the unshrunken results() so the y-axis still reflects the
## original Wald test FDR.
suppressPackageStartupMessages(library(apeglm))
res_raw    <- as.data.frame(results(dds, contrast = c('group', 'RVF', 'NF')))
res_shrunk <- as.data.frame(lfcShrink(dds, coef = 'group_RVF_vs_NF',
                                      type = 'apeglm'))
res <- res_shrunk
res$padj <- res_raw[rownames(res), 'padj']

mito_in_res <- intersect(Human.Mito$Symbol, rownames(res))
a <- data.frame(avg_log2FC = res[mito_in_res, 'log2FoldChange'],
                p_val_adj  = res[mito_in_res, 'padj'],
                row.names  = mito_in_res)
a <- a[complete.cases(a), ]
a$p_val_adj[a$p_val_adj < 1e-50] <- 1e-50
write.csv(cbind(gene = rownames(a), a),
          './output/PAB_RV_CM_mitocarto_volcano_pseudobulk.csv',
          row.names = FALSE)

## Map each gene to its MitoCarta category, then pick one color per category
## from a *categorical* palette (Dark2, 8 perceptually-distinct colors).
## Previous code interpolated colormap::rainbow_soft, which is a continuous
## gradient and made adjacent categories look the same.
keyvals <- anno[match(rownames(a), Human.Mito$Symbol)]
library(RColorBrewer)
.cat_levels  <- unique(keyvals)
.cat_palette <- RColorBrewer::brewer.pal(
  n    = max(3, min(length(.cat_levels), 8)),
  name = 'Dark2'
)[seq_along(.cat_levels)]
names(.cat_palette) <- .cat_levels
colors <- .cat_palette[keyvals]
names(colors) <- keyvals

library(EnhancedVolcano)
p_volcano <- EnhancedVolcano(a, lab = rownames(a),
  x = 'avg_log2FC', y = 'p_val_adj',
  xlim      = c(-6, 6),
  ylim      = c(0, 6),
  FCcutoff  = 0.5, pCutoff = 0.05,
  colCustom = colors,
  labFace   = 'italic',
  title     = NULL,
  subtitle  = NULL,
  caption   = NULL,
  legendPosition = 'top'
)
p_volcano <- p_volcano +
  guides(colour = guide_legend(nrow = 4, ncol = 2, byrow = FALSE,
                               override.aes = list(size = 3))) +
  theme(legend.position = 'top',
        legend.title    = element_blank())

ggsave('./output/PAB_RV_CM_mitocarto_volcano.pdf',
       p_volcano, width = 6, height = 9.38)
ggsave(file.path(V52_FIG_DIR, 'PAB_RV_CM_mitocarto_volcano.pdf'),
       p_volcano, width = 6, height = 9.38)


###############################################################################
##############  FIGURE S6J — Cross-species CM concordance panels  #############
###############################################################################
## v54 legend J: cross-species concordance of 3 CM programs
## (FAO, failure/fetal, GR/NR3C1).
## (i)   CM-subtype composition shift, mouse Nor/Mod/Sev vs human NF/pRV/RVF
## (ii)  FAO program score per-subject (mouse + human)
## (iii) Failure/fetal program score per-subject (NPPA/NPPB/ANKRD1/MYH7/XIRP1/ACTA1/BAG3)
## (iv)  GR axis: FKBP5 + NR3C1 per-subject

mouse_cmp <- list(c('Nor','Mod'), c('Nor','Sev'), c('Mod','Sev'))
human_cmp <- list(c('NF','pRV'), c('NF','RVF'), c('pRV','RVF'))

.subject_mean <- function(obj, meta_col, group_levels) {
  df <- aggregate(
    list(value = obj@meta.data[[meta_col]]),
    by  = list(patient = as.character(obj$patient),
               group   = as.character(obj$group)),
    FUN = mean
  )
  df$group <- factor(df$group, levels = group_levels)
  df
}

.subject_box <- function(df, palette, title, ylab, comparisons) {
  ggplot(df, aes(x = group, y = value, fill = group)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = PS$linewidth_mm) +
    geom_jitter(width = 0.12, size = PS$scatter_pt, alpha = 0.85) +
    scale_fill_manual(values = palette) +
    ggpubr::stat_compare_means(comparisons = comparisons,
                                method = 'wilcox.test', label = 'p.format',
                                size = 2.2 * PS$text_mm) +
    labs(title = title, x = NULL, y = ylab) +
    theme_v52(COMP_W) +
    theme(legend.position = 'none',
          plot.title  = element_text(size = 2.2 * PS$base_pt),
          axis.title  = element_text(size = 2.2 * PS$base_pt),
          axis.text   = element_text(size = 2.2 * PS$base_pt))
}

## (i) CM-subtype composition shift
.mouse_comp <- as.data.frame(prop.table(
  table(group = M2$group, subtype = M2$predicted.celltype), 1))
.human_comp <- as.data.frame(prop.table(
  table(group = M3$group, subtype = M3$Subnames), 1))
.mouse_comp$group <- factor(.mouse_comp$group, levels = c('Nor','Mod','Sev'))
.human_comp$group <- factor(.human_comp$group, levels = c('NF','pRV','RVF'))
write.csv(.mouse_comp, './output/PAB_mouse_CM_subtype_composition.csv', row.names=FALSE)
write.csv(.human_comp, './output/RV_human_CM_subtype_composition.csv', row.names=FALSE)

p_comp_mouse <- ggplot(.mouse_comp, aes(x = group, y = Freq, fill = subtype)) +
  geom_col(width = 0.7, colour = 'black', linewidth = PS$linewidth_mm) +
  labs(title = 'Mouse PAB — CM composition', x = NULL, y = 'Fraction of CMs') +
  theme_v52(COMP_W) +
  theme(legend.key.size = PS$legend_key, legend.title = element_blank())
p_comp_human <- ggplot(.human_comp, aes(x = group, y = Freq, fill = subtype)) +
  geom_col(width = 0.7, colour = 'black', linewidth = PS$linewidth_mm) +
  labs(title = 'Human RV — CM composition', x = NULL, y = 'Fraction of CMs') +
  theme_v52(COMP_W) +
  theme(legend.key.size = PS$legend_key, legend.title = element_blank())
ggsave('./output/S6J_PAB_RV_CM_composition.pdf',
       p_comp_mouse / p_comp_human, width = 4, height = 5.5)
ggsave(file.path(V52_FIG_DIR, 'S6J_PAB_RV_CM_composition.pdf'),
       p_comp_mouse / p_comp_human, width = 4, height = 5.5)

## (ii) FAO program score per-subject
fao_genes <- c('HMGCS2','HADHA','HADHB','ACADM','ACADVL','ACADL',
               'CPT1A','CPT1B','CPT2','ACOX1','ECH1','ECHS1','ETFDH','SLC25A20')
M2 <- AddModuleScore(M2, list(intersect(fao_genes, rownames(M2))),
                     name = 'FAO', ctrl = 50)
M3 <- AddModuleScore(M3, list(intersect(fao_genes, rownames(M3))),
                     name = 'FAO')
.mouse_fao <- .subject_mean(M2, 'FAO1', c('Nor','Mod','Sev'))
.human_fao <- .subject_mean(M3, 'FAO1', c('NF','pRV','RVF'))
write.csv(.mouse_fao, './output/PAB_mouse_CM_FAO_per_animal.csv', row.names=FALSE)
write.csv(.human_fao, './output/RV_human_CM_FAO_per_patient.csv', row.names=FALSE)

p_fao_mouse <- .subject_box(.mouse_fao, mhc_pal_mouse,
  'Mouse PAB — FAO program', 'Per-animal mean score', mouse_cmp)
p_fao_human <- .subject_box(.human_fao, mhc_pal_human,
  'Human RV — FAO program', 'Per-patient mean score', human_cmp)
p_fao <- p_fao_mouse | p_fao_human
pdf('./output/S6J_PAB_RV_CM_FAO.pdf', width = 5, height = 3)
print(p_fao); dev.off()
ggsave(file.path(V52_FIG_DIR, 'S6J_PAB_RV_CM_FAO.pdf'),
       p_fao, width = 5, height = 3)

## (iii) Failure / fetal-program score per-subject
fail_genes <- c('NPPA','NPPB','ANKRD1','MYH7','XIRP1','ACTA1','BAG3')
M2 <- AddModuleScore(M2, list(intersect(fail_genes, rownames(M2))),
                     name = 'Failure', ctrl = 50)
M3 <- AddModuleScore(M3, list(intersect(fail_genes, rownames(M3))),
                     name = 'Failure')
.mouse_fail <- .subject_mean(M2, 'Failure1', c('Nor','Mod','Sev'))
.human_fail <- .subject_mean(M3, 'Failure1', c('NF','pRV','RVF'))
write.csv(.mouse_fail, './output/PAB_mouse_CM_Failure_per_animal.csv', row.names=FALSE)
write.csv(.human_fail, './output/RV_human_CM_Failure_per_patient.csv', row.names=FALSE)

p_fail_mouse <- .subject_box(.mouse_fail, mhc_pal_mouse,
  'Mouse PAB — Failure program', 'Per-animal mean score', mouse_cmp)
p_fail_human <- .subject_box(.human_fail, mhc_pal_human,
  'Human RV — Failure program', 'Per-patient mean score', human_cmp)
p_fail <- p_fail_mouse | p_fail_human
pdf('./output/S6J_PAB_RV_CM_Failure_program.pdf', width = 5, height = 3)
print(p_fail); dev.off()
ggsave(file.path(V52_FIG_DIR, 'S6J_PAB_RV_CM_Failure_program.pdf'),
       p_fail, width = 5, height = 3)

## (iv) GR axis: FKBP5 + NR3C1 per-subject (single-gene module scores)
M2 <- AddModuleScore(M2, list('FKBP5'), name = 'FKBP5score', ctrl = 50)
M3 <- AddModuleScore(M3, list('FKBP5'), name = 'FKBP5score')
M2 <- AddModuleScore(M2, list('NR3C1'), name = 'NR3C1score', ctrl = 50)
M3 <- AddModuleScore(M3, list('NR3C1'), name = 'NR3C1score')
.mouse_fkbp5 <- .subject_mean(M2, 'FKBP5score1', c('Nor','Mod','Sev'))
.human_fkbp5 <- .subject_mean(M3, 'FKBP5score1', c('NF','pRV','RVF'))
.mouse_nr3c1 <- .subject_mean(M2, 'NR3C1score1', c('Nor','Mod','Sev'))
.human_nr3c1 <- .subject_mean(M3, 'NR3C1score1', c('NF','pRV','RVF'))
.mouse_gr <- rbind(cbind(.mouse_fkbp5, gene = 'FKBP5'),
                   cbind(.mouse_nr3c1, gene = 'NR3C1'))
.human_gr <- rbind(cbind(.human_fkbp5, gene = 'FKBP5'),
                   cbind(.human_nr3c1, gene = 'NR3C1'))
write.csv(.mouse_gr, './output/PAB_mouse_CM_GR_axis_per_animal.csv', row.names=FALSE)
write.csv(.human_gr, './output/RV_human_CM_GR_axis_per_patient.csv', row.names=FALSE)

.gr_box <- function(df, palette, title, comparisons) {
  ggplot(df, aes(x = group, y = value, fill = group)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = PS$linewidth_mm) +
    geom_jitter(width = 0.12, size = PS$scatter_pt, alpha = 0.85) +
    scale_fill_manual(values = palette) +
    ggpubr::stat_compare_means(comparisons = comparisons,
                                method = 'wilcox.test', label = 'p.format',
                                size = 2.2 * PS$text_mm) +
    facet_wrap(~ gene, scales = 'free_y') +
    labs(title = title, x = NULL, y = 'Module score') +
    theme_v52(COMP_W) +
    theme(legend.position = 'none',
          plot.title    = element_text(size = 2.2 * PS$base_pt),
          axis.title    = element_text(size = 2.2 * PS$base_pt),
          axis.text     = element_text(size = 2.2 * PS$base_pt),
          strip.text    = element_text(size = 2.2 * PS$base_pt))
}
p_gr_mouse <- .gr_box(.mouse_gr, mhc_pal_mouse, 'Mouse PAB — GR axis', mouse_cmp)
p_gr_human <- .gr_box(.human_gr, mhc_pal_human, 'Human RV — GR axis', human_cmp)
p_gr <- p_gr_mouse | p_gr_human
pdf('./output/S6J_PAB_RV_CM_GR_axis.pdf', width = 7, height = 3)
print(p_gr); dev.off()
ggsave(file.path(V52_FIG_DIR, 'S6J_PAB_RV_CM_GR_axis.pdf'),
       p_gr, width = 7, height = 3)


## Combined small-panels strip: FAO | Failure | NR3C1 only (drop FKBP5).
## Build a NR3C1-only pair so all three sections have the same column count
## and patchwork gives them equal widths in the strip.
.mouse_nr3c1_only <- subset(.mouse_gr, gene == 'NR3C1')
.human_nr3c1_only <- subset(.human_gr, gene == 'NR3C1')
.mouse_nr3c1_only$gene <- NULL
.human_nr3c1_only$gene <- NULL
p_nr3c1_mouse <- .subject_box(.mouse_nr3c1_only, mhc_pal_mouse,
  'Mouse PAB — NR3C1', 'Module score', mouse_cmp)
p_nr3c1_human <- .subject_box(.human_nr3c1_only, mhc_pal_human,
  'Human RV — NR3C1', 'Module score', human_cmp)
p_nr3c1 <- p_nr3c1_mouse | p_nr3c1_human

## Flatten to 6 plots in one row with explicit equal widths so each box plot
## gets exactly 1/6 of the figure regardless of axis-label length.
p_small_combined <- patchwork::wrap_plots(
  p_fao_mouse,   p_fao_human,
  p_fail_mouse,  p_fail_human,
  p_nr3c1_mouse, p_nr3c1_human,
  nrow = 1
) + patchwork::plot_layout(widths = rep(1, 6))

ggsave('./output/S6J_PAB_RV_CM_small_panels_combined.pdf',
       p_small_combined, width = 9, height = 3)
ggsave(file.path(V52_FIG_DIR, 'S6J_PAB_RV_CM_small_panels_combined.pdf'),
       p_small_combined, width = 9, height = 3)

#######################################
##########  FIGURE S6I — PAB Sirius Red  ##########
#######################################
## v54 legend I: Cardiac fibrosis Sirius Red % area per animal
## (Sham n=6, Mod n=7, Sev n=5; Wilcoxon).
## Fibrosis quantified as Sirius Red % area by PAB stage.
## Source: ./dependencies/shared/PAB_Sirius.xlsx (Sham / PAB Mod / PAB Sev,
## 7 animals per group laid out wide; some reps NA).
library(readxl)
library(ggpubr)

.sirius_raw <- read_excel('./dependencies/shared/PAB_Sirius.xlsx', sheet = 1)
sirius_long <- data.frame(
  group  = factor(rep(c('Sham','Mod','Sev'), each = 7),
                  levels = c('Sham','Mod','Sev')),
  sirius = as.numeric(unlist(.sirius_raw[1, ]))
)
sirius_long <- sirius_long[!is.na(sirius_long$sirius), ]
write.csv(sirius_long, './output/PAB_Sirius_long.csv', row.names = FALSE)

sirius_pal <- setNames(disease_pal[c('NF','pRV','RVF')], c('Sham','Mod','Sev'))

p_sirius <- ggplot(sirius_long, aes(x = group, y = sirius, fill = group)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = PS$linewidth_mm) +
  geom_jitter(width = 0.12, size = PS$scatter_pt, alpha = 0.85) +
  scale_fill_manual(values = sirius_pal) +
  ggpubr::stat_compare_means(
    comparisons = list(c('Sham','Mod'), c('Sham','Sev'), c('Mod','Sev')),
    method      = 'wilcox.test',
    label       = 'p.format',
    size        = PS$text_mm
  ) +
  labs(title = 'Cardiac fibrosis (Sirius Red)',
       x = NULL, y = 'Sirius Red (% area)') +
  theme_v52(COMP_W) +
  theme(legend.position = 'none')

save_figure(p_sirius, 'PAB_Sirius_box.pdf', width = 2.4, height = 2.1)
ggsave('./output/S6I_PAB_Sirius_box.pdf', p_sirius, width = 2.4, height = 2.1)

## S6G (v57): MitoCarta volcano (+ Sirius companion) composite. Was saved
## under the leftover wrong name 'SupplementaryFigure_7.pdf' (port artifact);
## content is correct — give it an S6-correct name the v57 emission map picks.
p_final <- (p_volcano | p_sirius) + patchwork::plot_layout(widths = c(3, 1))
save_figure(p_final, 'S6G_mitocarto_volcano.pdf', width = 8, height = 11)

# TODO v52: add other PAB echo/histology panels (new data)



#######################################
##  FIGURE S6H — handled inside S6D (within-mouse FC trajectories)  ##
#######################################
## v54 legend H lives inside the S6D Severe-vs-Sham/Mod-vs-Sham/
## Severe-vs-Mod trajectory scatter (r^2 reporting). No separate plot.












# new.cluster.ids <- c("Cm1","Cm2","Cm3","Cm4","Cm5","Cm6","Cm7","Cm8","Cm9","Cm10")
# names(new.cluster.ids) <- levels(M3)
# M3 <- RenameIdents(M3, new.cluster.ids)

# M3$Subnames <- M3@active.ident
# M3$SubNames_Groups <- paste(M1$Subnames,M3$group,sep='_')


# M3 <- AddModuleScore(M3, features=list(c('MALAT1')),assay="SCT",name="Clust0Score")
# M3 <- AddModuleScore(M3, features=list(c('FGF12','SH3RF2','KCNMB2','PRELID2')),assay="SCT",name="Clust1Score")
# M3 <- AddModuleScore(M3, features=list(c('TNNT2','TTN','MYBPC3','MYH7')),assay="SCT",name="Clust2Score")
# M3 <- AddModuleScore(M3, features=list(c('PALLD','MYO18B','MYPN','ANKRD1')),assay="SCT",name="Clust5Score")
# M3 <- AddModuleScore(M3, features=list(c('PDE3A','CDH2','PDLIM5')),assay="SCT",name="Clust4Score")
# M3 <- AddModuleScore(M3, features=list(c('AKAP13','OBSCN','LARGE1','THSD4')),assay="SCT",name="Clust3Score")
# M3 <- AddModuleScore(M3, features=list(c('PALLD','SORBS2','CAMK2D','CCSER1','PDLIM5')),assay="SCT",name="Clust7Score")
# M3 <- AddModuleScore(M3, features=list(c('AC020637.1','LINC02388')),assay="SCT",name="Clust6Score")
# M3 <- AddModuleScore(M3, features=list(c('MIR646HG')),assay="SCT",name="Clust8Score")
# M3 <- AddModuleScore(M3, features=list(c('GPC5','HS6ST3')),assay="SCT",name="Clust9Score")

# M3 <- SplitObject(M3, split.by = "patient")
# M3<-PrepSCTIntegration(M3)
# features<-SelectIntegrationFeatures(M3)
# M3.anchors<-FindIntegrationAnchors(M3,normalization.method = 'SCT',anchor.features = features, reduction = "rpca")
# M3 <- IntegrateData(anchorset = M3.anchors,normalization.method='SCT')

# DefaultAssay(M3) <- "integrated"

# M3 <- RunPCA(M3, npcs = 50, verbose = FALSE)
# M3 <- RunUMAP(M3, reduction = "pca", dims = 1:30)



# RNA <- M2@assays$RNA['counts']
# newnames <- human2mouse$human_name[match(rownames(RNA),human2mouse$mouse_name)]
# newnames[is.na(newnames)] <- rownames(RNA)[is.na(newnames)]
# rownames(RNA) <- newnames


# M4 <- CreateAssayObject(RNA)
# M2[['humanized']] <- M4
# DefaultAssay(M2) <- "humanized"
# M2[['SCT']] <- NULL

# M2 <- SCTransform(M2, vst.flavor = "v2",assay='humanized')
# M2 <- RunPCA(M2, npcs = 50, verbose = FALSE)

# anchors <- FindTransferAnchors(
#   reference = M3,
#   query = M2,
#   normalization.method = "SCT",
#   recompute.residuals=FALSE,
#   reference.reduction = "pca",
#   dims = 1:50
# )

# predictions <- TransferData(anchorset = anchors, refdata = M3$Subnames, dims = 1:50)


# M3 <- RunUMAP(M3, dims = 1:50, return.model = TRUE)
# M2 <- MapQuery(anchorset = anchors, reference = M3, query = M2,
#   refdata = list(celltype = "Subnames"), reference.reduction = "pca", reduction.model = "umap")

# score <- MappingScore(anchors)

# M2$map_score <- score

# p1 <- DimPlot(M3, reduction = "umap", group.by = "Subnames", label = TRUE, label.size = 3, repel = TRUE,raster=TRUE,pt.size=1.5) + NoLegend() + ggtitle("Reference annotations")
# p2 <- DimPlot(M2, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE, label.size = 3, pt.size=1.5,repel = TRUE,raster=TRUE) + NoLegend() + ggtitle("Query transferred labels")
# pdf(paste0('./output/', 'RV_PAB_CM_ref_mapped.pdf'), width=10, height=5)
# p1 + p2
# dev.off()


# pdf(paste0('./output/', 'PAB_CM_ref_mapped.pdf'), width=5, height=5)
# PlotEmbedding(M2,group.by='predicted.celltype',reduction = "ref.umap",point_size=0.2,plot_under=TRUE,plot_theme=umap_theme()+NoLegend(),raster_dpi=400,raster_scale=0.5)
# dev.off()


# pdf(paste0('./output/', 'RV_CM_ref_mapped.pdf'), width=5, height=5)
# PlotEmbedding(M3,group.by='Subnames',point_size=0.2,plot_under=TRUE,plot_theme=umap_theme()+NoLegend(),raster_dpi=400,raster_scale=0.5)
# dev.off()


# FeaturePlot(M2,'map_score',reduction = "ref.umap")




















# #######################################
# ######### Consensus hdWGCNA #########
# #######################################

# M1 <- readRDS(file = "./dependencies/shared/myeloid_subclust.rds")

# M1 <- FindNeighbors(M1)
# M1 <- FindClusters(M1,resolution=1)

# new.cluster.ids <- c('CCR2- rMac2','CCR2- rMac1','CCR2- rMac2','CCR2- rMac2',
# 	'CCR2+ rMac','CCR2- rMac1','CCR2- rMac1',
# 	'CCR2- rMac1','DCs','iMac','TREM2 Mac','CCR2- rMac2')
# names(new.cluster.ids) <- levels(M1)
# M1 <- RenameIdents(M1, new.cluster.ids)

# M2 <- subset(M1,idents=c('CCR2+ rMac','DCs'))
# M2<-RunPCA(M2)
# M2 <- M2%>% 
#   RunUMAP(reduction = "pca", dims = 1:50, verbose = F) %>% 
#   FindNeighbors(reduction = "pca", dims = 1:50) %>% 
#   FindClusters(resolution=1) %>% 
#   identity()

# #0 is CD1c- CD163+ MERTK+
# #1 is CCR2+ CD11c+
# #2 is DC
# #3 is mix of classical (FCN1+), non-classical (LILRB2, ITGAL, FCGR3A), and neutrophil-like (CSF3R)

# M1$Subnames<-M1@active.ident
# M1$Subsubnames <- M1$Subnames
# cells_DC = colnames(M2)[M2$seurat_clusters==2]
# cells_mono = colnames(M2)[M2$seurat_clusters==3]
# cells_rMac = union(colnames(M2)[M2$seurat_clusters==0],colnames(M2)[M2$seurat_clusters==1])

# levels(M1$Subsubnames) <- c(levels(M1$Subsubnames), 'Mono')
# M1$Subsubnames[colnames(M1) %in% cells_DC] = 'DCs'
# M1$Subsubnames[colnames(M1) %in% cells_mono] = 'Mono'
# M1$Subsubnames[colnames(M1) %in% cells_rMac] = 'CCR2+ rMac'


# M2 <- M1
# M2<-RunPCA(M2)
# M2 <- RunHarmony(M1,'patient')
# M2 <- M2%>% 
#   RunUMAP(reduction = "harmony", dims = 1:50, verbose = F) %>% 
#   FindNeighbors(reduction = "harmony", dims = 1:50) %>% 
#   FindClusters(resolution=1) %>% 
#   identity()

# M1 <- M2
# M1 <- SetIdent(M1,value='Subsubnames')

# M1 <- AddModuleScore(M1, features=list(c("TREM2","GPNMB","MITF","SPP1")),assay="SCT",name="TREM2_Mac_Score")
# M1 <- AddModuleScore(M1, features=list(c("CLEC9A","ZBTB46","CD1C","CD226")),assay="SCT",name="DC_Score")
# M1 <- AddModuleScore(M1, features=list(c("FCN1","LILRB2","ITGAL","CSF3R")),assay="SCT",name="Mono_Score")
# M1 <- AddModuleScore(M1, features=list(c("CCR2","CX3CR1","ITGAX")),assay="SCT",name="CCR2+_rMac_Score")
# M1 <- AddModuleScore(M1, features=list(c("LYVE1","FOLR2","SIGLEC1","F13A1")),assay="SCT",name="CCR2-_rMac1_Score")
# M1 <- AddModuleScore(M1, features=list(c("RBMS3","PLA2G5","EBF1")),assay="SCT",name="CCR2-_rMac2_Score")
# M1 <- AddModuleScore(M1, features=list(c("IL1B","CCL3","CCL4","CXCL3","CXCL8")),assay="SCT",name="iMac_Score")





# M2 <- readRDS('./output/PAB_data_clean.rds')

# M2 <- SetIdent(M2, value = "Names")
# M2$group <- M2$orig.ident


# M3 <- subset(M2, Names %in% c("Myeloid"))

# M3 <- RunPCA(M3)
# M3 <- FindNeighbors(M3, dims = 1:50)
# M3 <- FindClusters(M3, resolution = 0.5)
# M3 <- RunUMAP(M3, dims = 1:10)

# #0 is CCR2- rMac
# #1 looks like cardiac contam
# #2 is
# #3 is DC and CCR, HLA high
# #4 B cells
# #5 is something non immune
# #6 EC

# M3$Subnames <- M3@active.ident
# M3 <- subset(M3, Subnames %in% c('0','1','2','3') )
# M3 <- RunUMAP(M3, dims = 1:50)

# #0 rMac
# #1
# #2 
# #3 HLA
# M3 <- subset(M3, Subnames %in% c('0','3') )
# M3 <- RunPCA(M3)
# M3 <- RunUMAP(M3, dims = 1:10)

# labels <- c('rMac','HLA')

# names(labels) <- levels(M3)
# M3 <- RenameIdents(M3, labels)
# M3$Subnames <- M3@active.ident

# human2mouse <- read.csv('./dependencies/shared/human2mouse.csv',header=F)
# idx <- match(unique(human2mouse[,2]),human2mouse[,2])
# human2mouse<-human2mouse[idx,]
# colnames(human2mouse) <-c('human_name', 'mouse_name')

# RNA <- M3@assays$RNA['counts']
# newnames <- human2mouse$human_name[match(rownames(RNA),human2mouse$mouse_name)]
# #newnames[is.na(newnames)] <- rownames(RNA)[is.na(newnames)]
# rownames(RNA) <- newnames


# M2 <- CreateAssayObject(RNA[!is.na(newnames),])
# M3[['humanized']] <- M2
# DefaultAssay(M3) <- "humanized"
# M3[['SCT']] <- NULL

# M3 <- SCTransform(M3, vst.flavor = "v2",assay='humanized',variable.features.n = 2000)
# M3 <- RunPCA(M3, npcs = 50, verbose = FALSE)

# M3$Species <- 'mouse'
# M1$Species <- 'human'

# DefaultAssay(M3) <- "humanized"
# DefaultAssay(M1) <- "RNA"
# M3[['SCT']] <- NULL
# M1[['SCT']] <- NULL
# M3$Subsubnames <- M3$Subnames


# shared <- intersect(rownames(M1),rownames(M3))


# M3[['RNA']] <- NULL
# M3[['RNA']] <- M3[['humanized']]
# DefaultAssay(M3) <- "RNA"
# M3[['humanized']] <- NULL

# M1$decontXcounts <- NULL

# M2 <- merge(M1[shared,],M3[shared,])

# M2 <- NormalizeData(M2)
# M2<-FindVariableFeatures(M2,nfeatures=3000)
# M2<-ScaleData(M2)
# M2 <- RunPCA(M2)
# M2 <- RunHarmony(M2, group.by.vars = c('Species','patient'), dims=1:50,lambda=c(1,0.01))

# M2 <- RunUMAP(M2, reduction='harmony',dims = 1:30)

# DimPlot(M2, group.by='patient', split.by='Species', raster=FALSE, label=TRUE) + umap_theme()
# DimPlot(M2, group.by='patient', split.by='Subsubnames', raster=FALSE, label=TRUE) + umap_theme()

# M2$Names <- M2$Subsubnames
# M2$Names[M2$Names == 'CCR2- rMac2'] <- 'CCR2- rMac'
# M2$Names[M2$Names == 'CCR2- rMac1'] <- 'CCR2- rMac'
# M2$Names[M2$Names == 'rMac'] <- 'CCR2- rMac'
# M2$Names[M2$Names == 'DCs'] <- 'HLA'
# M2$Names[M2$Names == 'Mono'] <- 'HLA'
# M2$Names[M2$Names == 'CCR2+ rMac'] <- 'HLA'
# M2$Names[M2$Names == 'iMac'] <- 'HLA'




# M2 <- SetupForWGCNA(
#   M2,
#   gene_select = "fraction",
#   fraction = 0.05,
#   wgcna_name = 'Myeloid_consensus'
# )

# # construct metacells:
# M2 <- MetacellsByGroups(
#   M2,
#   group.by = c("Names", 'Species'),
#   k = 25,
#   max_shared = 12,
#   min_cells = 50,
#   target_metacells = 250,
#   reduction = 'harmony',
#   ident.group = 'Names'
# )

# M2 <- NormalizeMetacells(M2)

# # setup expression matrices for each species in astrocytes
# M2 <- SetMultiExpr(
#   M2,
#   group_name = "CCR2- rMac",
#   group.by = "Names",
#   multi.group.by ="Species",
#   multi_groups = NULL
# )


# # identify soft power thresholds
# M2 <- TestSoftPowersConsensus(M2)

# # plot soft power results
# plot_list <-  PlotSoftPowers(M2)
# consensus_groups <- unique(M2$Species)
# p_list <- lapply(1:length(consensus_groups), function(i){
#   cur_group <- consensus_groups[[i]]
#   plot_list[[i]][[1]] + ggtitle(paste0(cur_group)) + theme(plot.title=element_text(hjust=0.5))
# })
# library(patchwork)
# wrap_plots(p_list, ncol=2)

# # consensus network analysis
# M2 <- ConstructNetwork(
#   M2,soft_power=c(7,7),
#   consensus=TRUE,
#   tom_name = "Species_Consensus",
#   overwrite_tom = TRUE
# )

# PlotDendrogram(M2, main='Resident macrophage cross species dendrogram')


# M2 <- ModuleEigengenes(M2, group.by.vars=c("Species","patient"))
# M2 <- ModuleConnectivity(M2, group_name ='CCR2- rMac', group.by='Names')

# # re-name modules
# M2 <- ResetModuleNames(M2, new_name = "rMac-CM")

# # visualize network with UMAP
# M2 <- RunModuleUMAP(
#   M2,
#   n_hubs = 5,
#   n_neighbors=10,
#   min_dist=0.3,
#   spread=2,
#   target_weight=0.1,
#   supervised=TRUE
# )

# # get the hub gene UMAP table from the seurat object
# umap_df <- GetModuleUMAP(M2)

# # plot with ggplot
# plot_df <- umap_df

# # compute coordinates for cluster labels
# centroid_df <- data.frame()
# for(cur_cluster in unique(plot_df[['module']])){
#     cur_meta <- plot_df[plot_df[['module']] == cur_cluster,]
#     df <- data.frame(
#       cluster = cur_cluster,
#       UMAP1 = mean(cur_meta$UMAP1),
#       UMAP2 = mean(cur_meta$UMAP2)
#     )
#   centroid_df <- rbind(centroid_df, df)
# }

# # plot with ggplot
# p <- ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
#     ggrastr::rasterise(geom_point(
#      color=umap_df$color,
#      size=umap_df$kME*2
#    ), dpi=500, scale=0.5) +
#     umap_theme() +
#     theme(
#       plot.margin = margin(0,0,0,0),
#       plot.title = element_text(hjust=0.5)
#     ) + ggtitle('HLA Consensus WGCNA') +
#     ggrepel::geom_text_repel(data = centroid_df, label=centroid_df$cluster, color='black', max.overlaps=Inf, size=3)



# hub_genes <- GetHubGenes(M2, 3)

# # add annotation
# anno_genes <- hub_genes$gene_name
# plot_df$anno <- ifelse(plot_df$gene %in% anno_genes, umap_df$gene, '')

# plot_df_anno <- subset(plot_df, anno != '')
# p <-  plot_df %>%
#     ggplot(aes(x=UMAP1, y=UMAP2, color=module)) +
#     ggrastr::rasterise(
#       geom_point(
#         inherit.aes=FALSE,
#         data=plot_df,
#         aes(x=UMAP1, y=UMAP2, color=module),
#         color=plot_df$color,
#         size=plot_df$kME*2,
#       ), dpi=500, dpi_scale=0.5) +
#     geom_point(
#       inherit.aes = FALSE,
#       data = plot_df_anno,
#       shape=21, color='black',
#       fill=plot_df_anno$color,
#       size=plot_df_anno$kME*2,
#       aes(x=UMAP1, y=UMAP2, fill=module)
#     ) +
#     # add labels
#     ggrepel::geom_text_repel(data = centroid_df, label=centroid_df$cluster, color='black', max.overlaps=Inf, size=3, fontface='bold') +
#     geom_text_repel(label=plot_df$anno, max.overlaps=Inf, color='black', fontface='italic', size=3) +
#     umap_theme() + NoLegend() +
#     coord_equal() +
#     theme(
#       plot.margin = margin(0,0,0,0)
#     )

# pdf(paste0('./output/rMac_hubgene_umap_ggplot.pdf'), width=8, height=8)
# print(p)
# dev.off()

# dbs <-c('GO_Biological_Process_2023','GO_Cellular_Component_2023','GO_Molecular_Function_2023','Reactome_2022', 'ChEA_2022',"LINCS_L1000_Chem_Pert_up")


# M2 <- RunEnrichr(
#   M2,
#   dbs=dbs, # character vector of enrichr databases to test
#   max_genes = Inf # number of genes per module to test. use max_genes = Inf to choose all genes!
# )

# enrich_df <- GetEnrichrTable(M2)

# EnrichrBarPlot(
#   M2,
#   outdir = "./output/Myeloid_rMac_consensus_modules", # name of output directory
#   n_terms = 10, # number of enriched terms to show (sometimes more show if there are ties!!!)
#   plot_size = c(5,7), # width, height of the output .pdfs
#   logscale=TRUE # do you want to show the enrichment as a log scale?
# )

# modules <- GetModules(M2) %>% subset(module != 'grey')

# cat(rownames(subset(modules,module=='rMac-CM1')),sep='\n')
# cat(rownames(subset(modules,module=='rMac-CM2')),sep='\n')
# cat(rownames(subset(modules,module=='rMac-CM6')),sep='\n')
# cat(rownames(subset(modules,module=='rMac-CM8')),sep='\n')
# cat(rownames(subset(modules,module=='rMac-CM9')),sep='\n')
# cat(rownames(subset(modules,module=='rMac-CM11')),sep='\n')
# cat(rownames(subset(modules,module=='rMac-CM12')),sep='\n')

# write.csv(modules,"./output/consensusWGCNA_rMac_modules.csv")

# modules <- GetModules(M2)
# color_df <- modules %>% subset(module!='grey') %>%
#   select(c(module, color)) %>% distinct %>%
#   rename(c(group=module, colour=color))
# mods <- levels(modules$module)
# mods <- mods[mods!='grey']

# # helper function to wrap text
# wrapText <- function(x, len) {
#     sapply(x, function(y) paste(strwrap(y, len), collapse = "\n"), USE.NAMES = FALSE)
# }

# combined_output <- GetEnrichrTable(M2)
# selected_terms <- subset(combined_output,db=="GO_Biological_Process_2023")

# # subset selected terms
# selected_terms <- subset(selected_terms, Term %in% selected_terms$Term & P.value < 0.05)
# idx_top_1 <- match(unique(selected_terms$module),selected_terms$module)
# idx_top_3 <- sort(c(idx_top_1,idx_top_1+1,idx_top_1+2))

# selected_terms<-selected_terms[idx_top_1,]
# #key_terms <- read.csv('./output/bulkRNA_GOterms_ofinterest.csv')
# #selected_terms <- subset(selected_terms,Term %in% key_terms[[1]])
# #selected_terms <- subset(combined_output, Term %in% selected_terms$Term & P.value < 0.05)


# selected_terms$group <- factor(
#   as.character(selected_terms$module),
#   levels = mods
# )


# # set max pval
# quantile(-log(selected_terms$P.value), 0.95)
# max_p <- 10

# selected_terms$logp <- -log(selected_terms$P.value)
# selected_terms$logp <- ifelse(selected_terms$logp > max_p, max_p, selected_terms$logp)

# # remove GO Term ID
# library(stringr)
# selected_terms$Term <- str_replace(selected_terms$Term, " \\s*\\([^\\)]+\\)", "")

# selected_terms <- selected_terms %>%
#   arrange(group)


# selected_terms$wrap <- wrapText(selected_terms$Term, 35)

# selected_terms$Term <- factor(
#   as.character(selected_terms$Term),
#   levels = rev(unique(as.character(selected_terms$Term)))
# )

# library(viridis)

# # GO Term dot plot
# p <- selected_terms %>%
#   ggplot(aes(x = group, y = Term, color =logp, size=log(Combined.Score))) +
#   geom_point() +
#   scale_color_stepsn(colors=rev(magma(256))) +
#   RotatedAxis() + xlab('') + ylab('') +
#   theme(
#     axis.title.x = element_blank(),
#     axis.title.y = element_blank(),
#     panel.border = element_rect(size=1, color='black', fill=NA),
#     axis.line.x = element_blank(),
#     axis.line.y = element_blank(),
#     axis.text.x = element_blank(),
#     axis.ticks.x = element_blank(),
#     plot.margin = margin(0,0,0,0),
#     panel.grid = element_line(size=0.25, color='lightgrey')
#   )


# mapping <- labels2colors(1:100)
# #color_df$group<-paste0('M',match(color_df$group,mapping))
# lvls <- stringr::str_sort(unique(color_df$group), numeric = TRUE)
# color_df$group <- factor(color_df$group, levels = lvls)
# # make the colorbar as its own heatmap
# color_df$var <- 1
# colorbar <- color_df %>%
#   ggplot(aes(x=group, y=var, fill=group)) +
#   geom_tile() +
#   scale_fill_manual(values=color_df$colour) +
#   coord_equal() +
#   NoLegend() + RotatedAxis() +
#   theme(
#     plot.title=element_blank(),
#     axis.line=element_blank(),
#     axis.ticks.y =element_blank(),
#     axis.text.y = element_blank(),
#     axis.title = element_blank(),
#     plot.margin=margin(0,0,0,0),
#   )


# pdf(paste0('./output/Consensus_WGNCA_rMac_selected_GO_terms.pdf'), width=13, height=8)
# p / colorbar #+ plot_layout(heights=c(20,1))
# dev.off()



# ### HLA

# M2 <- SetMultiExpr(
#   M2,
#   group_name = "HLA",
#   group.by = "Names",
#   multi.group.by ="Species",
#   multi_groups = NULL
# )


# # consensus network analysis
# M2 <- ConstructNetwork(
#   M2,soft_power=c(7,7),
#   consensus=TRUE,
#   tom_name = "Species_HLA_Consensus",
#   overwrite_tom = TRUE
# )

# PlotDendrogram(M2, main='HLA macrophage cross species dendrogram')

# M2 <- ModuleEigengenes(M2, group.by.vars=c("Species","patient"))
# M2 <- ModuleConnectivity(M2, group_name ='HLA', group.by='Names')

# # re-name modules
# M2 <- ResetModuleNames(M2, new_name = "HLA-CM")

# # visualize network with UMAP
# M2 <- RunModuleUMAP(
#   M2,
#   n_hubs = 5,
#   n_neighbors=10,
#   min_dist=0.3,
#   spread=2,
#   target_weight=0.1,
#   supervised=TRUE
# )

# # get the hub gene UMAP table from the seurat object
# umap_df <- GetModuleUMAP(M2)

# # plot with ggplot
# plot_df <- umap_df

# # compute coordinates for cluster labels
# centroid_df <- data.frame()
# for(cur_cluster in unique(plot_df[['module']])){
#     cur_meta <- plot_df[plot_df[['module']] == cur_cluster,]
#     df <- data.frame(
#       cluster = cur_cluster,
#       UMAP1 = mean(cur_meta$UMAP1),
#       UMAP2 = mean(cur_meta$UMAP2)
#     )
#   centroid_df <- rbind(centroid_df, df)
# }

# # plot with ggplot
# p <- ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
#     ggrastr::rasterise(geom_point(
#      color=umap_df$color,
#      size=umap_df$kME*2
#    ), dpi=500, scale=0.5) +
#     umap_theme() +
#     theme(
#       plot.margin = margin(0,0,0,0),
#       plot.title = element_text(hjust=0.5)
#     ) + ggtitle('HLA Consensus WGCNA') +
#     ggrepel::geom_text_repel(data = centroid_df, label=centroid_df$cluster, color='black', max.overlaps=Inf, size=3)



# hub_genes <- GetHubGenes(M2, 3)

# # add annotation
# anno_genes <- hub_genes$gene_name
# plot_df$anno <- ifelse(plot_df$gene %in% anno_genes, umap_df$gene, '')

# plot_df_anno <- subset(plot_df, anno != '')
# p <-  plot_df %>%
#     ggplot(aes(x=UMAP1, y=UMAP2, color=module)) +
#     ggrastr::rasterise(
#       geom_point(
#         inherit.aes=FALSE,
#         data=plot_df,
#         aes(x=UMAP1, y=UMAP2, color=module),
#         color=plot_df$color,
#         size=plot_df$kME*2,
#       ), dpi=500, dpi_scale=0.5) +
#     geom_point(
#       inherit.aes = FALSE,
#       data = plot_df_anno,
#       shape=21, color='black',
#       fill=plot_df_anno$color,
#       size=plot_df_anno$kME*2,
#       aes(x=UMAP1, y=UMAP2, fill=module)
#     ) +
#     # add labels
#     ggrepel::geom_text_repel(data = centroid_df, label=centroid_df$cluster, color='black', max.overlaps=Inf, size=3, fontface='bold') +
#     geom_text_repel(label=plot_df$anno, max.overlaps=Inf, color='black', fontface='italic', size=3) +
#     umap_theme() + NoLegend() +
#     coord_equal() +
#     theme(
#       plot.margin = margin(0,0,0,0)
#     )

# pdf(paste0('./output/HLA_hubgene_umap_ggplot.pdf'), width=8, height=8)
# print(p)
# dev.off()



# dbs <-c('GO_Biological_Process_2023','GO_Cellular_Component_2023','GO_Molecular_Function_2023','Reactome_2022', 'ChEA_2022',"LINCS_L1000_Chem_Pert_up")


# M2 <- RunEnrichr(
#   M2,
#   dbs=dbs, # character vector of enrichr databases to test
#   max_genes = Inf # number of genes per module to test. use max_genes = Inf to choose all genes!
# )

# enrich_df <- GetEnrichrTable(M2)

# EnrichrBarPlot(
#   M2,
#   outdir = "./output/Myeloid_HLA_consensus_modules", # name of output directory
#   n_terms = 10, # number of enriched terms to show (sometimes more show if there are ties!!!)
#   plot_size = c(5,7), # width, height of the output .pdfs
#   logscale=TRUE # do you want to show the enrichment as a log scale?
# )

# modules <- GetModules(M2) %>% subset(module != 'grey')

# write.csv(modules,"./output/consensusWGCNA_HLA_modules.csv")

# cat(rownames(subset(modules,module=='HLA-CM10')),sep='\n')
# cat(rownames(subset(modules,module=='HLA-CM20')),sep='\n')


# modules <- GetModules(M2)
# color_df <- modules %>% subset(module!='grey') %>%
#   select(c(module, color)) %>% distinct %>%
#   rename(c(group=module, colour=color))
# mods <- levels(modules$module)
# mods <- mods[mods!='grey']

# # helper function to wrap text
# wrapText <- function(x, len) {
#     sapply(x, function(y) paste(strwrap(y, len), collapse = "\n"), USE.NAMES = FALSE)
# }

# combined_output <- GetEnrichrTable(M2)
# selected_terms <- subset(combined_output,db=="GO_Biological_Process_2023")

# # subset selected terms
# selected_terms <- subset(selected_terms, Term %in% selected_terms$Term & P.value < 0.05)
# idx_top_1 <- match(unique(selected_terms$module),selected_terms$module)
# idx_top_3 <- sort(c(idx_top_1,idx_top_1+1,idx_top_1+2))

# selected_terms<-selected_terms[idx_top_1,]
# #key_terms <- read.csv('./output/bulkRNA_GOterms_ofinterest.csv')
# #selected_terms <- subset(selected_terms,Term %in% key_terms[[1]])
# #selected_terms <- subset(combined_output, Term %in% selected_terms$Term & P.value < 0.05)


# selected_terms$group <- factor(
#   as.character(selected_terms$module),
#   levels = mods
# )


# # set max pval
# quantile(-log(selected_terms$P.value), 0.95)
# max_p <- 10

# selected_terms$logp <- -log(selected_terms$P.value)
# selected_terms$logp <- ifelse(selected_terms$logp > max_p, max_p, selected_terms$logp)

# # remove GO Term ID
# library(stringr)
# selected_terms$Term <- str_replace(selected_terms$Term, " \\s*\\([^\\)]+\\)", "")

# selected_terms <- selected_terms %>%
#   arrange(group)


# selected_terms$wrap <- wrapText(selected_terms$Term, 35)

# selected_terms$Term <- factor(
#   as.character(selected_terms$Term),
#   levels = rev(unique(as.character(selected_terms$Term)))
# )

# library(viridis)

# # GO Term dot plot
# p <- selected_terms %>%
#   ggplot(aes(x = group, y = Term, color =logp, size=log(Combined.Score))) +
#   geom_point() +
#   scale_color_stepsn(colors=rev(magma(256))) +
#   RotatedAxis() + xlab('') + ylab('') +
#   theme(
#     axis.title.x = element_blank(),
#     axis.title.y = element_blank(),
#     panel.border = element_rect(size=1, color='black', fill=NA),
#     axis.line.x = element_blank(),
#     axis.line.y = element_blank(),
#     axis.text.x = element_blank(),
#     axis.ticks.x = element_blank(),
#     plot.margin = margin(0,0,0,0),
#     panel.grid = element_line(size=0.25, color='lightgrey')
#   )


# mapping <- labels2colors(1:100)
# #color_df$group<-paste0('M',match(color_df$group,mapping))
# lvls <- stringr::str_sort(unique(color_df$group), numeric = TRUE)
# color_df$group <- factor(color_df$group, levels = lvls)
# # make the colorbar as its own heatmap
# color_df$var <- 1
# colorbar <- color_df %>%
#   ggplot(aes(x=group, y=var, fill=group)) +
#   geom_tile() +
#   scale_fill_manual(values=color_df$colour) +
#   coord_equal() +
#   NoLegend() + RotatedAxis() +
#   theme(
#     plot.title=element_blank(),
#     axis.line=element_blank(),
#     axis.ticks.y =element_blank(),
#     axis.text.y = element_blank(),
#     axis.title = element_blank(),
#     plot.margin=margin(0,0,0,0),
#   )


# pdf(paste0('./output/Consensus_WGNCA_HLA_selected_GO_terms.pdf'), width=13, height=8)
# p / colorbar #+ plot_layout(heights=c(20,1))
# dev.off()


# #######################################
# #########  MYELOID PAB PILOT  #########
# #######################################


# #1343 outlier
# #1632,1467 RVF
# #1681 1691 1697 NF

# M1 <- readRDS('./output/PAB_data_clean.rds')

# M1 <- SetIdent(M1, value = "Names")
# M1$group <- M1$orig.ident


# M2 <- subset(M1, Names %in% c("Myeloid"))

# M2 <- RunPCA(M2)
# M2 <- FindNeighbors(M2, dims = 1:50)
# M2 <- FindClusters(M2, resolution = 0.5)
# M2 <- RunUMAP(M2, dims = 1:10)

# #0 is CCR2- rMac
# #1 looks like cardiac contam
# #2 is
# #3 is DC and CCR, HLA high
# #4 B cells
# #5 is something non immune
# #6 EC

# M2$Subnames <- M2@active.ident
# M2 <- subset(M2, Subnames %in% c('0','1','2','3') )
# M2 <- RunUMAP(M2, dims = 1:50)

# #0 rMac
# #1
# #2 
# #3 HLA
# M2 <- subset(M2, Subnames %in% c('0','3') )
# M2 <- RunPCA(M2)
# M2 <- RunUMAP(M2, dims = 1:10)

# labels <- c('rMac','HLA')

# names(labels) <- levels(M2)
# M2 <- RenameIdents(M2, labels)
# M2$Subnames <- M2@active.ident




# #### LOAD BULK DATA AND EMBED INTO ALL PAB
# dir <- '/Users/ikuz/Documents/Mouse_PAB_myeloid/output'

# samples <- read.csv(file.path(dir, "meta.csv"), header = TRUE)
# samples$Batch <- factor(samples$Batch)
# samples$Pressure.Loading <- factor(samples$Pressure.Loading)


# files <- file.path(dir,'nascent', samples$ID, "abundance.h5")
# names(files) <- samples$ID
# tx2gene <- read.table('/Users/ikuz/Documents/Mouse_PAB_myeloid/index/t2g.txt',fill=T)
# tx2gene = data.frame(TXNAME=tx2gene$V1,GENEID=tx2gene$V3)

# txi.kallisto <- tximport(files, type = "kallisto", tx2gene = tx2gene)


# rownames(samples) <- samples$ID
# bulk <- CreateSeuratObject(counts = txi.kallisto$counts, meta.data = data.frame(samples))

# bulk <- SCTransform(bulk, vst.flavor = "v2")
# bulk <- RunPCA(bulk, npcs = 10, verbose = FALSE)


# M1 <- SplitObject(M1, split.by = "patient")
# M1<-PrepSCTIntegration(M1)
# features<-SelectIntegrationFeatures(M1)
# M1.anchors<-FindIntegrationAnchors(M1,normalization.method = 'SCT',anchor.features = features, reference = c(1,2,3), reduction = "rpca")
# M1 <- IntegrateData(anchorset = M1.anchors,normalization.method='SCT')

# DefaultAssay(M1) <- "integrated"

# M1 <- RunPCA(M1, npcs = 50, verbose = FALSE)
# M1 <- RunUMAP(M1, reduction = "pca", dims = 1:30)


# anchors <- FindTransferAnchors(
#   reference = M1,
#   query = bulk,
#   normalization.method = "SCT",
#   recompute.residuals=FALSE,
#   reference.reduction = "pca",
#   dims = 1:50,
#   k.score =15
# )


# predictions <- TransferData(anchorset = anchors, refdata = M1$Names, dims = 1:50,k.weight=15)

# M1 <- RunUMAP(M1, dims = 1:50, return.model = TRUE)

# bulk <- MapQuery(anchorset = anchors, reference = M1, query = bulk,
# 	refdata = list(celltype = "Names"), reference.reduction = "pca", 
# 	reduction.model = "umap",transferdata.args = list(k.weight=15))

# #score <- MappingScore(anchors)

# #bulk$map_score <- score

# p1 <- DimPlot(M1, reduction = "umap", group.by = "Names", label = TRUE, label.size = 3, repel = TRUE,raster=TRUE,pt.size=1.5) + NoLegend() + ggtitle("Reference annotations")
# p2 <- DimPlot(bulk, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE, label.size = 3, pt.size=1.5,repel = TRUE,raster=TRUE) + NoLegend() + ggtitle("Query transferred labels")
# pdf(paste0('./output/', 'RV_PAB_bulkMyeloid_ref_mapped.pdf'), width=10, height=5)
# p1 + p2
# dev.off()



# #### LOAD BULK DATA AND EMBED INTO MYELOID PAB

# bulk <- CreateSeuratObject(counts = txi.kallisto$counts, meta.data = data.frame(samples))

# bulk <- SCTransform(bulk, vst.flavor = "v2")
# bulk <- RunPCA(bulk, npcs = 10, verbose = FALSE)


# #M2 <- SplitObject(M2, split.by = "patient")
# #M2<-PrepSCTIntegration(M2)
# #features<-SelectIntegrationFeatures(M2)
# #M2.anchors<-FindIntegrationAnchors(M2,normalization.method = 'SCT',anchor.features = features, reduction = "rpca")
# #M2 <- IntegrateData(anchorset = M2.anchors,normalization.method='SCT')

# #DefaultAssay(M2) <- "integrated"

# #M1 <- RunPCA(M2, npcs = 50, verbose = FALSE)
# #M1 <- RunUMAP(M2, reduction = "pca", dims = 1:30)
# M3 <- subset(M1, Names %in% c("Myeloid"))

# true_myeloid <- colnames(M2)
# M3 <- subset(M3,cells=true_myeloid)
# M3$Subnames = M2$Subnames

# anchors <- FindTransferAnchors(
#   reference = M3,
#   query = bulk,
#   normalization.method = "SCT",
#   recompute.residuals=FALSE,
#   reference.reduction = "pca",
#   dims = 1:50,
#   k.score =15
# )


# predictions <- TransferData(anchorset = anchors, refdata = M3$Subnames, dims = 1:50,k.weight=10)

# M3 <- RunUMAP(M3, dims = 1:50, return.model = TRUE)

# bulk <- MapQuery(anchorset = anchors, reference = M3, query = bulk,
# 	refdata = list(celltype = "Subnames"), reference.reduction = "pca", 
# 	reduction.model = "umap",transferdata.args = list(k.weight=10))

# #score <- MappingScore(anchors)

# #bulk$map_score <- score


# bulk$group <- paste0(bulk$Origin,'_',bulk$Type)
# p1 <- DimPlot(M3, reduction = "umap", group.by = "Subnames", label = TRUE, label.size = 3, repel = TRUE,raster=TRUE,pt.size=1.5) + NoLegend() + ggtitle("Reference annotations")
# p2 <- DimPlot(bulk, reduction = "ref.umap", group.by = "group", label = TRUE, label.size = 3, pt.size=1.5,repel = TRUE,raster=TRUE) + ggtitle("Query transferred labels")
# pdf(paste0('./output/', 'RV_PAB_Myeloid_Only_bulkMyeloid_ref_mapped.pdf'), width=10, height=5)
# p1 + p2
# dev.off()

# bulk_RV = subset(bulk, Origin=='RV')

# bulk_RV<-SetIdent(bulk_RV,value='Type')

# M3<-SetIdent(M3,value='group')


# DefaultAssay(M3) <- 'SCT'

# a <- FindMarkers(bulk_RV,ident.1='PAB',ident.2='Sham',min.pct=0.1,logfc.threshold=0)

# M3 <- PrepSCTFindMarkers(M3)
# b <-  FindMarkers(M3,ident.1='Sev',ident.2='Nor',min.pct=0.25,recorrect_umi=F,logfc.threshold=0)
# shared <- intersect(rownames(a),rownames(b))
# dataset <- data.frame(bulk=a[shared,]$avg_log2FC,sn=b[shared,]$avg_log2FC)
# rownames(dataset) <- shared
# labs <- rownames(dataset)

# cor(dataset[,1],dataset[,2])

# pdf(paste0('./output/', 'RV_PAB_Myeloid_Only_bulkMyeloid_scatter.pdf'), width=10, height=10)
# ggplot(dataset, aes(x = sn, y=bulk)) + geom_point() + 
#   geom_text_repel(label=labs,max.overlaps=25) + theme_classic()
# dev.off()

# cat(rownames(subset(dataset, bulk > 0.1 & sn > 0.1)),sep='\n')
# cat(rownames(subset(dataset, bulk < 0.1 & sn < 0.1)),sep='\n')



# ###Compare to human RV


# bulk_RV <- AddModuleScore(bulk_RV,list(c('Ciita','Cd74','H2-Ab1','H2-Aa',
# 	'H2-Eb1','H2-Eb2','H2-Ob','H2-DMb1','H2-DMb2','H2-DMa','H2-Oa')),name='MHCII')

# VlnPlot(bulk_RV,'MHCII1')

# human2mouse <- read.csv('./dependencies/shared/human2mouse.csv',header=F)
# idx <- match(unique(human2mouse[,2]),human2mouse[,2])
# human2mouse<-human2mouse[idx,]
# colnames(human2mouse) <-c('human_name', 'mouse_name')
# gluc_response <- "VIT;VKORC1L1;ERRFI1;AHCYL1;STEAP4;TRAF3IP2;GALNT15;SERPINE1;JADE1;SLA;CBLB;MT1X;EPS8;CCND3;BMPER;RASSF4;RPS6KA2;ANPEP;C1RL;MAP3K6;IL6R;PDGFRA;MLIP;SCARA5;IL1R1;EBF1;TTC7A;CRISPLD2;SPARCL1;FKBP5;NNMT;LPAR1;SLC1A3;PLA2G5;NID1;ACACB;ZFP36L2;PIK3R5;C3;SCFD2;LPXN;HACL1;SRGAP2;SLC38A2;SLC19A2;S100A10;KLHL29;GADD45B;ZBTB16;ELL2;CORO2B;IGF2R;NFATC4;DERA;SULT1B1;MAFB;BCL6;TMEM236;TBXAS1;NDUFAF2;RGL3;SERPINA3;MCFD2;PTPRS;ELN;PTEN;FMN1;HIF3A;TFCP2L1;PTH1R;SYNE3;CTSS;PTPRG;RNF157;ADAMTS2;C1QTNF1;IMPA2;SH3PXD2B;FLVCR2;EFHD1;AOX1;CERS6;ZHX3;KLF13;ANXA2;IFNGR1;GPX3;NCOA3;SLC39A11;NGF;OSMR;SLC39A14;TGFBR2;TGFBR3;PSMA6;ARHGAP10;MMP14;TBC1D2;SLC7A7;SLC7A8;GFOD1;DPYD;PICK1;FAM20C;COL6A3;PLIN2;ITGA5;MOCS1;ERGIC1;TMEM45A;KANK1;C1S;ADCY3;TFPI;FSTL1;TMEM165;HDAC7;KIAA0513;MTHFD1L;CLMN;PTK2B;PTPN18;GALNT6;GSN;NEGR1;TPK1;CCDC57;TXNRD1;GSR;SUSD1;LHFPL2;MERTK;KLF9;IL18R1"
# gluc_response <- str_split(gluc_response[1],';')[[1]]
# gluc_response_mouse <- human2mouse$mouse_name[match(gluc_response,human2mouse$human_name)]

# bulk_RV <- AddModuleScore(bulk_RV,list(gluc_response_mouse),name='nr3c1')

# VlnPlot(bulk_RV,'nr3c11')














# human2mouse <- read.csv('./dependencies/shared/human2mouse.csv',header=F)
# idx <- match(unique(human2mouse[,2]),human2mouse[,2])
# human2mouse<-human2mouse[idx,]
# colnames(human2mouse) <-c('human_name', 'mouse_name')

# consensus_modules <- read.csv("./dependencies/shared/bulk_heart_modules.csv")
# consensus_modules <- consensus_modules[,1:3]

# idx_match<- match(consensus_modules$gene_name,human2mouse$human_name)

# consensus_modules$gene_name <- human2mouse$mouse_name[idx_match]
# consensus_modules <- consensus_modules[!is.na(consensus_modules$gene_name),]
# rownames(consensus_modules) <- consensus_modules$gene_name

# consensus_modules <- subset(consensus_modules, gene_name %in% rownames(M2))
# # remove duplicate gene names
# consensus_modules <- consensus_modules[match(unique(consensus_modules$gene_name), consensus_modules$gene_name),]
# library(dplyr)
# score_calc <- consensus_modules %>% group_by(module) %>% group_split()
# module_colors <- unique(unlist(lapply(score_calc,'[[','module')))

# mapping <- labels2colors(1:100)

# module_colors <- paste0('M',match(module_colors,mapping))
# M2 <- AddModuleScore(M2,lapply(score_calc,'[[','gene_name'),name="module_score")
# cols_current <- colnames(M2@meta.data)
# cols_current[startsWith(colnames(M2@meta.data),'module_score')] <- paste0('module_',module_colors)
# colnames(M2@meta.data) <- cols_current

# M2 <- SetIdent(M2,value='Subnames')

# modules_int <- c('M1','M3','M4',"M8")

# pdf(paste0('./output/', 'PAB_myeloid_dot_subclust.pdf'), width=5, height=2)

# p <- DotPlot(M2,paste0('module_',modules_int),group.by='Subnames',dot.min=0,col.min=-1,col.max=1,scale.min=50,scale.max=100) +
#   RotatedAxis() + ylab('')+ xlab('')+
#   scale_color_gradient2(high='red', mid='grey95', low='blue') +
#   theme(
#     panel.border = element_rect(size=1,fill=NA, color='black'),
#     axis.line.x = element_blank(),
#     axis.line.y = element_blank()
# ) 
# p
# dev.off()

# M2$group <- factor(M2$group,levels = c('Nor','Mod','Sev'))



# pdf(paste0('./output/', 'PAB_myeloid_dot_disease.pdf'), width=5, height=2.5)

# p <- DotPlot(M2,paste0('module_',modules_int),group.by='group',dot.min=0,col.min=-1,col.max=1,scale.min=50,scale.max=100) +
#   RotatedAxis() + ylab('')+ xlab('')+
#   scale_color_gradient2(high='red', mid='grey95', low='blue') +
#   theme(
#     panel.border = element_rect(size=1,fill=NA, color='black'),
#     axis.line.x = element_blank(),
#     axis.line.y = element_blank()
# ) 
# p
# dev.off()


# #MHCII

# M2 <- AddModuleScore(M2,list(c('Ciita','Cd74','H2-Ab1','H2-Aa',
# 	'H2-Eb1','H2-Eb2','H2-Ob','H2-DMb1','H2-DMb2','H2-DMa','H2-Oa')),name='MHCII')

# pdf(paste0('./output/', 'PAB_myeloid_MHC.pdf'), width=3, height=3)

# VlnPlot(subset(M2,Subnames=='HLA'),'MHCII1',group.by='group',pt.size=0)
# dev.off()



# gluc_response <- "VIT;VKORC1L1;ERRFI1;AHCYL1;STEAP4;TRAF3IP2;GALNT15;SERPINE1;JADE1;SLA;CBLB;MT1X;EPS8;CCND3;BMPER;RASSF4;RPS6KA2;ANPEP;C1RL;MAP3K6;IL6R;PDGFRA;MLIP;SCARA5;IL1R1;EBF1;TTC7A;CRISPLD2;SPARCL1;FKBP5;NNMT;LPAR1;SLC1A3;PLA2G5;NID1;ACACB;ZFP36L2;PIK3R5;C3;SCFD2;LPXN;HACL1;SRGAP2;SLC38A2;SLC19A2;S100A10;KLHL29;GADD45B;ZBTB16;ELL2;CORO2B;IGF2R;NFATC4;DERA;SULT1B1;MAFB;BCL6;TMEM236;TBXAS1;NDUFAF2;RGL3;SERPINA3;MCFD2;PTPRS;ELN;PTEN;FMN1;HIF3A;TFCP2L1;PTH1R;SYNE3;CTSS;PTPRG;RNF157;ADAMTS2;C1QTNF1;IMPA2;SH3PXD2B;FLVCR2;EFHD1;AOX1;CERS6;ZHX3;KLF13;ANXA2;IFNGR1;GPX3;NCOA3;SLC39A11;NGF;OSMR;SLC39A14;TGFBR2;TGFBR3;PSMA6;ARHGAP10;MMP14;TBC1D2;SLC7A7;SLC7A8;GFOD1;DPYD;PICK1;FAM20C;COL6A3;PLIN2;ITGA5;MOCS1;ERGIC1;TMEM45A;KANK1;C1S;ADCY3;TFPI;FSTL1;TMEM165;HDAC7;KIAA0513;MTHFD1L;CLMN;PTK2B;PTPN18;GALNT6;GSN;NEGR1;TPK1;CCDC57;TXNRD1;GSR;SUSD1;LHFPL2;MERTK;KLF9;IL18R1"
# gluc_response <- str_split(gluc_response[1],';')[[1]]


# gluc_response_mouse <- human2mouse$mouse_name[match(gluc_response,human2mouse$human_name)]


# M2 <- AddModuleScore(M2,list(gluc_response_mouse),name='nr3c1')


# pdf(paste0('./output/', 'PAB_myeloid_nr3c1.pdf'), width=3, height=3)

# VlnPlot(M2,'nr3c11',group.by='group',pt.size=0)
# dev.off()

# M1$group <- M1$orig.ident
# pdf(paste0('./output/', 'Wnt1.pdf'), width=3, height=3)
# DimPlot(M1)
# dev.off()


# #######################################
# #############  FIGURE S6H  ############
# #######################################
# library(enrichR)


# bulk_modules <- consensus_modules
# bulk_modules$module <- match(bulk_modules$module,mapping)
# dbs <- c('GO_Biological_Process_2023','GO_Cellular_Component_2023','GO_Molecular_Function_2023','ChEA_2022','Reactome_2022')

# #Run enrichment by cell type
# M2$Names_group <- paste0(M2$Subnames,'_',M2$group)
# Idents(M2) <- "Names_group"
# combined_set <- data.frame()
# combined_output <- data.frame()

# mods_idx <- c(1,3,4,8)
# cell_types <- unique(M2$Subnames)
# comparison <- list(c("Sev","Mod"),c("Mod","Nor"),c("Sev","Nor"))
# for (i in mods_idx){
# 	for (j in cell_types){
# 		for (k in comparison){
# 			key_genes <- subset(bulk_modules,module %in% c(i))$gene_name
# 			key_genes <- key_genes[key_genes %in% rownames(M2)]

# 			gene_set <- FindMarkers(M2, ident.1 = paste0(j,"_",k[1]), ident.2 = paste0(j,"_",k[2]),features=key_genes,recorrect_umi=F)
			
# 			gene_set<-subset(gene_set,p_val_adj<0.05)
# 			if (length(rownames(gene_set))==0){next}
# 			gene_set$module <- paste0('M',i)
# 			gene_set$color <- mapping[i]
# 			gene_set$comparison <- paste0(k[1],'_',k[2])
# 			gene_set$celltype <- j

# 			if (length(combined_set) == 0){
# 				combined_set <- gene_set
# 			}
# 			else {
# 				combined_set <- rbind(combined_set,gene_set)
# 			}

# 			gene_enrich <- subset(gene_set,avg_log2FC<0)
# 			enriched <- enrichR::enrichr(rownames(gene_enrich), dbs)
# 			Sys.sleep(5)
# 			for(db in names(enriched)){
# 		  		cur_df <- enriched[[db]]
# 		  		if (nrow(cur_df) > 1){
# 		    		cur_df$db <- db
# 		    		cur_df$module <- paste0('M',i)
# 		    		cur_df$celltype <- j
# 		    		cur_df$comparison <- paste0(k[1],'_',k[2])
# 		    		cur_df$color <- mapping[i]
# 		    		cur_df$direction <- 'down'
# 		    		combined_output <- rbind(combined_output, cur_df)
# 		  		}
# 			}

# 			gene_enrich <- subset(gene_set,avg_log2FC>0)
# 			enriched <- enrichR::enrichr(rownames(gene_enrich), dbs)
# 			Sys.sleep(5)
# 			for(db in names(enriched)){
# 		  		cur_df <- enriched[[db]]
# 		  		if (nrow(cur_df) > 1){
# 		    		cur_df$db <- db
# 		    		cur_df$module <- paste0('M',i)
# 		    		cur_df$celltype <- j
# 		    		cur_df$comparison <- paste0(k[1],'_',k[2])
# 		    		cur_df$color <- mapping[i]
# 		    		cur_df$direction <- 'up'
# 		    		combined_output <- rbind(combined_output, cur_df)
# 		  		}
# 			}
# 		}
# 	}
# }


# #selected_terms <- subset(combined_output,db=="ChEA_2022")
# selected_terms <- subset(combined_output,direction=="down")
# selected_terms <- subset(selected_terms,comparison=="Sev_Nor")
# selected_terms <- subset(selected_terms,Adjusted.P.value<0.05)
# #selected_terms <- subset(selected_terms,module %in% c('M1'))
# selected_terms



# Idents(M2) <- "group"

# library(EnhancedVolcano)

# gene_set <- FindMarkers(M2, ident.1 = 'Sev', ident.2 = 'Nor',recorrect_umi=F,features=subset(bulk_modules,module==1)$gene_name)

# pdf(paste0('./output/', 'PAB_M1_Myeloid_module_volcano_all_RVF_vs_NF.pdf'), width=8, height=6)

# EnhancedVolcano(gene_set,lab=rownames(gene_set),
# 	x='avg_log2FC',y='p_val_adj',
# 	FCcutoff = 0.1,pCutoff=0.05) + coord_flip()
# dev.off()

# gene_set <- FindMarkers(M2, ident.1 = 'Sev', ident.2 = 'Nor',recorrect_umi=F,features=subset(bulk_modules,module==8)$gene_name)

# pdf(paste0('./output/', 'PAB_M8_Myeloid_module_volcano_all_RVF_vs_NF.pdf'), width=8, height=6)

# EnhancedVolcano(gene_set,lab=rownames(gene_set),
# 	x='avg_log2FC',y='p_val_adj',
# 	FCcutoff = 0.1,pCutoff=0.05) + coord_flip()
# dev.off()

# #######################################
# #############  FIGURE S6I  ############
# #######################################

# M2 <- subset(M1, Names %in% c("CM"))

# M2 <- RunPCA(M2)
# M2 <- RunHarmony(M2,'patient')
# M2 <- FindNeighbors(M2, dims = 1:50,reduction = "harmony")
# M2 <- FindClusters(M2, resolution = 0.5,reduction = "harmony")
# M2 <- RunUMAP(M2, dims = 1:50,reduction = "harmony")

# M2$Names <- M2@active.ident
# markers<-FindAllMarkers(M2,recorrect_umi=F)

# pdf(paste0('./output/', 'PAB_CM_snUMAP.pdf'), width=5, height=5)
# PlotEmbedding(M2,group.by='Names',point_size=1,plot_under=TRUE,plot_theme=umap_theme()+NoLegend(),raster_dpi=400,raster_scale=0.5)
# dev.off()



# M3 <- readRDS(file = "./dependencies/shared/cm_subclust.rds")

# new.cluster.ids <- c("Cm1","Cm2","Cm3","Cm4","Cm5","Cm6","Cm7","Cm8","Cm9","Cm10")
# names(new.cluster.ids) <- levels(M3)
# M3 <- RenameIdents(M3, new.cluster.ids)

# M3$Subnames <- M3@active.ident
# M3$SubNames_Groups <- paste(M1$Subnames,M3$group,sep='_')


# M3 <- AddModuleScore(M3, features=list(c('MALAT1')),assay="SCT",name="Clust0Score")
# M3 <- AddModuleScore(M3, features=list(c('FGF12','SH3RF2','KCNMB2','PRELID2')),assay="SCT",name="Clust1Score")
# M3 <- AddModuleScore(M3, features=list(c('TNNT2','TTN','MYBPC3','MYH7')),assay="SCT",name="Clust2Score")
# M3 <- AddModuleScore(M3, features=list(c('PALLD','MYO18B','MYPN','ANKRD1')),assay="SCT",name="Clust5Score")
# M3 <- AddModuleScore(M3, features=list(c('PDE3A','CDH2','PDLIM5')),assay="SCT",name="Clust4Score")
# M3 <- AddModuleScore(M3, features=list(c('AKAP13','OBSCN','LARGE1','THSD4')),assay="SCT",name="Clust3Score")
# M3 <- AddModuleScore(M3, features=list(c('PALLD','SORBS2','CAMK2D','CCSER1','PDLIM5')),assay="SCT",name="Clust7Score")
# M3 <- AddModuleScore(M3, features=list(c('AC020637.1','LINC02388')),assay="SCT",name="Clust6Score")
# M3 <- AddModuleScore(M3, features=list(c('MIR646HG')),assay="SCT",name="Clust8Score")
# M3 <- AddModuleScore(M3, features=list(c('GPC5','HS6ST3')),assay="SCT",name="Clust9Score")

# #DefaultAssay(M3) <- "RNA"

# #M3[["RNA"]] <- split(M3[["RNA"]], f = M3$patient)
# #M3[['SCT']] <- NULL
# #M3[['decontXcounts']] <- NULL
# #M3 <- SCTransform(M3, vst.flavor = "v2")
# #M3 <- RunPCA(M3, npcs = 50, verbose = FALSE)
# M3 <- SplitObject(M3, split.by = "patient")
# M3<-PrepSCTIntegration(M3)
# features<-SelectIntegrationFeatures(M3)
# M3.anchors<-FindIntegrationAnchors(M3,normalization.method = 'SCT',anchor.features = features, reduction = "rpca")
# M3 <- IntegrateData(anchorset = M3.anchors,normalization.method='SCT')

# DefaultAssay(M3) <- "integrated"

# M3 <- RunPCA(M3, npcs = 50, verbose = FALSE)
# M3 <- RunUMAP(M3, reduction = "pca", dims = 1:30)

# # RenameGenesSeurat <- function(obj = ls.Seurat[[i]],
# #                               newnames = HGNC.updated[[i]]$Suggested.Symbol,
# #                               assay = "RNA",
# #                               slots = c("data", "counts", "meta.features")) {
# #   #
# #   message("RenameGenesSeurat, assay: ", assay)
# #   warning("Run this before integration and downstream processing. It only attempts to change
# #           @counts, @data, and @meta.features in obj@assays$YOUR_ASSAY.", immediate. = TRUE)

# #   stopifnot(
# #     "Unequal gene name sets: nrow(assayobj) != nrow(newnames):" =
# #       length(Features(obj, assay = assay)) == length(newnames)
# #   )

# #   if (obj@version < 5) warning("obj@version < 5. Old versions are not supported. Update the obj!", immediate. = TRUE)

# #   if ("scale.data" %in% slots) {
# #     n_genes_sc_dta <- nrow(obj@assays[[assay]]$"scale.data")
# #     stopifnot(
# #       "scale.data does has different number of genes than newnames!" =
# #         n_genes_sc_dta == length(newnames)
# #     )
# #   }

# #   LayersFound <- SeuratObject::Layers(obj@assays[[assay]])
# #   #print("Present: ", sort(LayersFound))

# #   slots <- sort(intersect(slots, LayersFound))
# #   #print("Replaced: ", slots)

# #     for (slotX in slots) {
# #     print(slotX)
# #     if (slotX == "scale.data") browser()
# #     nrO <- nrow(SeuratObject::GetAssayData(object = obj, assay = assay, layer = slotX))
# #     obj <- .check_and_rename(obj, assay, newnames = newnames, layer.name = slotX)
# #     nrN <- nrow(SeuratObject::GetAssayData(object = obj, assay = assay, layer = slotX))
# #     stopifnot(nrN == nrO)
# #   }
# #   return(obj)
# # }


# # .check_and_rename <- function(obj, assay, newnames, layer.name) {
# #   cat(layer.name, fill = TRUE)

# #   length_newnames <- length(newnames)
# #   length_orig_names <- length(Features(obj, assay = assay))

# #   stopifnot(
# #     is(obj, "Seurat"),
# #     is.character(assay),
# #     is.character(layer.name),
# #     is.character(newnames),
# #     length_orig_names == length_newnames
# #   )

# #   assayobj <- obj@assays[[assay]]
# #   feature.list <- rownames(assayobj@features@.Data)

# #   if (length(feature.list) == length(newnames)) {
# #     rownames(assayobj@features@.Data) <- newnames
# #     nrX <- length(rownames(assayobj@features@.Data))
# #   } else {
# #     paste0("length feature.list", length(feature.list), "length newnames", length(newnames))
# #     stop()
# #   }

# #   if (layer.name %in% SeuratObject::Layers(assayobj)) {
# #     matrix_n <- SeuratObject::LayerData(assayobj, layer = layer.name)
# #     nr1 <- nrow(matrix_n)

# #     if (all(dim(matrix_n)) > 0) {
# #       # browser()
# #       stopifnot(nrow(matrix_n) == length(newnames))

# #       if ("dgCMatrix" %in% class(matrix_n)) {
# #         message(assay, "@", layer.name, " is of type dgeCMatrix!")
# #         matrix_n@Dimnames[[1]] <- newnames
# #       } else if ("matrix" %in% class(matrix_n)) {
# #         message(assay, "@", layer.name, " is of type Matrix!")
# #         rownames(matrix_n) <- newnames
# #       } else if ("data.frame" %in% class(matrix_n)) {
# #         message(assay, "@", layer.name, " is of type data.frame!")
# #         rownames(matrix_n) <- newnames
# #       } else {
# #         warning(">>> No renaming: ", assay, "@", layer.name,
# #           " not of type dgeCMatrix / Matrix / data.frame.",
# #           immediate. = TRUE
# #         )
# #       }
# #       stopifnot(nr1 == nrow(matrix_n))

# #       SeuratObject::LayerData(assayobj, layer = layer.name) <- matrix_n
# #       nr3 <- nrow(SeuratObject::LayerData(assayobj, layer = layer.name))
# #       stopifnot(nr3 == nrX)
# #     }
# #   } else {
# #     warning(paste(">>>", assay, "@", layer.name, "does not exist!"), immediate. = TRUE)
# #   }
# #   # obj <- SetAssayData(obj, layer = layer.name, new.data = matrix_n)
# #   obj@assays[[assay]] <- assayobj
# #   return(obj)
# # }

# # DefaultAssay(M2) <- "RNA"


# # newnames <- human2mouse$human_name[match(rownames(M2),human2mouse$mouse_name)]
# # newnames[is.na(newnames)] <- rownames(M2)[is.na(newnames)]

# # RenameGenesSeurat(M2,newnames = newnames,assay='RNA')



# RNA <- M2@assays$RNA['counts']
# newnames <- human2mouse$human_name[match(rownames(RNA),human2mouse$mouse_name)]
# newnames[is.na(newnames)] <- rownames(RNA)[is.na(newnames)]
# rownames(RNA) <- newnames


# M4 <- CreateAssayObject(RNA)
# M2[['humanized']] <- M4
# DefaultAssay(M2) <- "humanized"
# M2[['SCT']] <- NULL

# M2 <- SCTransform(M2, vst.flavor = "v2",assay='humanized')
# M2 <- RunPCA(M2, npcs = 50, verbose = FALSE)

# anchors <- FindTransferAnchors(
#   reference = M3,
#   query = M2,
#   normalization.method = "SCT",
#   recompute.residuals=FALSE,
#   reference.reduction = "pca",
#   dims = 1:50
# )

# predictions <- TransferData(anchorset = anchors, refdata = M3$Subnames, dims = 1:50)


# M3 <- RunUMAP(M3, dims = 1:50, return.model = TRUE)
# M2 <- MapQuery(anchorset = anchors, reference = M3, query = M2,
# 	refdata = list(celltype = "Subnames"), reference.reduction = "pca", reduction.model = "umap")

# score <- MappingScore(anchors)

# M2$map_score <- score

# p1 <- DimPlot(M3, reduction = "umap", group.by = "Subnames", label = TRUE, label.size = 3, repel = TRUE,raster=TRUE,pt.size=1.5) + NoLegend() + ggtitle("Reference annotations")
# p2 <- DimPlot(M2, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE, label.size = 3, pt.size=1.5,repel = TRUE,raster=TRUE) + NoLegend() + ggtitle("Query transferred labels")
# pdf(paste0('./output/', 'RV_PAB_CM_ref_mapped.pdf'), width=10, height=5)
# p1 + p2
# dev.off()


# pdf(paste0('./output/', 'PAB_CM_ref_mapped.pdf'), width=5, height=5)
# PlotEmbedding(M2,group.by='predicted.celltype',reduction = "ref.umap",point_size=0.2,plot_under=TRUE,plot_theme=umap_theme()+NoLegend(),raster_dpi=400,raster_scale=0.5)
# dev.off()


# pdf(paste0('./output/', 'RV_CM_ref_mapped.pdf'), width=5, height=5)
# PlotEmbedding(M3,group.by='Subnames',point_size=0.2,plot_under=TRUE,plot_theme=umap_theme()+NoLegend(),raster_dpi=400,raster_scale=0.5)
# dev.off()


# FeaturePlot(M2,'map_score',reduction = "ref.umap")

# #######################################
# #############  FIGURE S6J  ############
# #######################################


# RV_marks <- FindAllMarkers(M)

# RV_marks_sig <- subset(RV_marks, p_val_adj<0.05 & avg_log2FC>0)

# idx<-match(RV_marks_sig$cluster,unique(RV_marks_sig$cluster))

# marks <- split(RV_marks_sig$gene,idx)

# M2 <- AddModuleScore(M2,marks,name='RV_marks')
# DefaultAssay(M3) <- "SCT"
# M3 <- AddModuleScore(M3,marks,name='RV_marks')

# M2$predicted.celltype <- factor(M2$predicted.celltype,levels=levels(M3))

# pdf(paste0('./output/', 'PAB_CM_marker_scores.pdf'), width=6, height=4)

# DotPlot(M2,c('RV_marks1','RV_marks2','RV_marks3','RV_marks4','RV_marks5'
# 	,'RV_marks6','RV_marks7','RV_marks8','RV_marks9','RV_marks10'),
# 	col.min = 0,group.by='predicted.celltype') +
# scale_x_discrete(labels=c('Cm1','Cm2','Cm3','Cm4','Cm5','Cm6','Cm7','Cm8','Cm9','Cm10'))
# dev.off()


# DotPlot(M3,c('RV_marks1','RV_marks2','RV_marks3','RV_marks4','RV_marks5'
# 	,'RV_marks6','RV_marks7','RV_marks8','RV_marks9','RV_marks10'),col.min = 0) +
# scale_x_discrete(labels=c('Cm1','Cm2','Cm3','Cm4','Cm5','Cm6','Cm7','Cm8','Cm9','Cm10'))

# pdf(paste0('./output/', 'PAB_CM_scores.pdf'), width=3.75, height=4.25)
# VlnPlot(M2,'map_score',group.by='predicted.celltype',pt.size=0)
# dev.off()

# #######################################
# #############  FIGURE S6L  ############
# #######################################

# M3 <- SetIdent(M3,value = 'group')
# M2 <- SetIdent(M2,value = 'group')


# a <- FindMarkers(M3,ident.1='RVF',ident.2='NF')
# b <-  FindMarkers(M2,ident.1='Sev',ident.2='Nor')
# shared <- intersect(rownames(a),rownames(b))
# dataset <- data.frame(RV=a[shared,]$avg_log2FC,PAB=b[shared,]$avg_log2FC)
# rownames(dataset) <- shared
# labs <- rownames(dataset)
# #labs[abs(dataset$PAB - dataset$RV)<1] <- NA


# pdf(paste0('./output/', 'PAB_vs_RV_CM__dot.pdf'), width=6, height=8)
# ggplot(dataset, aes(x = RV, y=PAB)) + geom_point() + 
#   geom_text_repel(label=labs,max.overlaps=50) + theme_classic()
# dev.off()


# a <- FindMarkers(M3,ident.1='RVF',ident.2='NF')
# b <-  FindMarkers(M2,ident.1='Mod',ident.2='Nor')
# shared <- intersect(rownames(a),rownames(b))
# dataset <- data.frame(RV=a[shared,]$avg_log2FC,PAB=b[shared,]$avg_log2FC)
# rownames(dataset) <- shared
# labs <- rownames(dataset)
# #labs[abs(dataset$PAB - dataset$RV)<1] <- NA


# pdf(paste0('./output/', 'PAB_vs_RV_CM_Mod_Nor_dot.pdf'), width=6, height=8)
# ggplot(dataset, aes(x = RV, y=PAB)) + geom_point() + 
#   geom_text_repel(label=labs,max.overlaps=50) + theme_classic()
# dev.off()


# a <- FindMarkers(M2,ident.1='Sev',ident.2='Nor')
# b <-  FindMarkers(M2,ident.1='Mod',ident.2='Nor')
# shared <- intersect(rownames(a),rownames(b))
# dataset <- data.frame(RV=a[shared,]$avg_log2FC,PAB=b[shared,]$avg_log2FC)
# rownames(dataset) <- shared
# labs <- rownames(dataset)
# #labs[abs(dataset$PAB - dataset$RV)<1] <- NA


# pdf(paste0('./output/', 'PAB_CM_Sev_Nor_Mod_Nor_dot.pdf'), width=6, height=8)
# ggplot(dataset, aes(x = RV, y=PAB)) + geom_point() + 
#   geom_text_repel(label=labs,max.overlaps=50) + theme_classic()
# dev.off()


# a <- FindMarkers(M2,ident.1='Sev',ident.2='Mod')
# b <-  FindMarkers(M2,ident.1='Mod',ident.2='Nor')
# shared <- intersect(rownames(a),rownames(b))
# dataset <- data.frame(RV=a[shared,]$avg_log2FC,PAB=b[shared,]$avg_log2FC)
# rownames(dataset) <- shared
# labs <- rownames(dataset)
# labs[abs(dataset$PAB - dataset$RV)<1] <- NA


# pdf(paste0('./output/', 'PAB_CM_Sev_Mod_Mod_Nor_dot.pdf'), width=6, height=8)
# ggplot(dataset, aes(x = RV, y=PAB)) + geom_point() + 
#   geom_text_repel(label=labs,max.overlaps=10) + theme_classic()
# dev.off()
# #######################################
# #############  FIGURE S6M  ############
# #######################################

# consensus_modules <- read.csv("./dependencies/shared/bulk_heart_modules.csv")
# consensus_modules <- consensus_modules[,1:3]
# consensus_modules <- subset(consensus_modules, gene_name %in% rownames(M2))
# consensus_modules <- consensus_modules[match(unique(consensus_modules$gene_name), consensus_modules$gene_name),]
# library(dplyr)
# mapping <- labels2colors(1:100)

# score_calc <- consensus_modules %>% group_by(module) %>% group_split()
# module_colors <- unique(unlist(lapply(score_calc,'[[','module')))
# module_colors <- paste0('M',match(module_colors,mapping))


# M2 <- AddModuleScore(M2,lapply(score_calc,'[[','gene_name'),name="module_score",ctrl = 50)

# cols_current <- colnames(M2@meta.data)
# cols_current[startsWith(colnames(M2@meta.data),'module_score')] <- paste0('module_',module_colors)
# colnames(M2@meta.data) <- cols_current


# M2 <- SetIdent(M2, value = "group")
# M2$group <- factor(M2$group,levels=c('Nor','Mod','Sev'))


# pdf(paste0('./output/', 'PAB_seurat_dot_CM.pdf'), width=5, height=2.5)

# p <- DotPlot(M2,paste0('module_',
#   c('M2','M10','M12','M25','M26','M28')),dot.min=0,col.min=0,col.max=2,group.by='group') +
#   RotatedAxis() + ylab('')+ xlab('')+
#   scale_color_gradient2(high='red', mid='grey95', low='blue') +
#   theme(
#     panel.border = element_rect(size=1,fill=NA, color='black'),
#     axis.line.x = element_blank(),
#     axis.line.y = element_blank()
# ) 
# p

# dev.off()

# #######################################
# #############  FIGURE S6N  ############
# #######################################
# library("readxl")

# Human.Mito <- read_excel("./dependencies/shared/Human.MitoCarta3.0.xls", sheet = "A Human MitoCarta3.0")
# Mouse.Mito <- read_excel("./dependencies/shared/Mouse.MitoCarta3.0.xls", sheet = "A Mouse MitoCarta3.0")


# M3 <- AddModuleScore(M3,list(Human.Mito$Symbol),name='mito')
# M2 <- AddModuleScore(M2,list(union(Human.Mito$Symbol,Mouse.Mito$Symbol)),name='mito',ctrl = 50)

# p1<-VlnPlot(M2,'mito1',group.by='group',pt.size=0)
# p2<-VlnPlot(M3,'mito1',group.by='group',pt.size=0)

# pdf(paste0('./output/', 'PAB_RV_CM_mitocarto.pdf'), width=3, height=4)

# p1 / p2
# dev.off()

# anno <- trimws(unlist(lapply(lapply(lapply(Human.Mito$MitoCarta3.0_MitoPathways,str_split,'>'),'[[',1),'[[',1)))

# anno[anno == "Small molecule transport | Signaling"] = "Small molecule transport" 


# a<-FindMarkers(M3,ident.1='RVF',ident.2='NF',features = intersect(Human.Mito$Symbol,rownames(M3)))
# a$p_val_adj[a$p_val_adj < 1e-50] = 1e-50

# keyvals <- anno[match(rownames(a),Human.Mito$Symbol)]
# library(colormap)

# colors <- colormap(colormap=colormaps$rainbow_soft, nshades=length(unique(keyvals)))

# colors <- colors[match(keyvals,unique(keyvals))]
# names(colors) <- keyvals


# library(EnhancedVolcano)
# pdf(paste0('./output/', 'CM_RV_Mito.pdf'), width=6, height=11)

# EnhancedVolcano(a,lab=rownames(a),
#   x='avg_log2FC',y='p_val_adj',
#   FCcutoff = 0.1,pCutoff=0.05,colCustom = colors,xlim=c(-6.25,3),ylim=c(0,52))
# dev.off()

###############################################################################
## v57 standardized per-panel emission — keyed to the v57 S6 legend
## (.figure_run_logs/v57_figure_legends.md), NOT the reference's internal
## letters. Reference writes legacy ./output/<name>.pdf; copy each to
## V52_FIG_DIR/Figure_S6_panel_<L>.pdf. v57 S6:
##   A mouse→human CM ref-map | B human-CM-marker score in mouse |
##   C per-subcluster mapping score | D cross-species FC scatter |
##   E mouse-CM WGCNA module dot by disease | F per-subject MitoCarta |
##   G pseudobulk DESeq2 MitoCarta volcano (apeglm) |
##   H within-mouse FC traj (i)+(ii) | I Sirius-Red fibrosis |
##   J cross-species 3-program concordance (FAO/Failure/GR).
###############################################################################
.s6_v57 <- list(
  A   = c('RV_PAB_CM_ref_mapped.pdf','PAB_CM_ref_mapped.pdf'),
  B   = 'PAB_CM_marker_scores.pdf',
  C   = 'PAB_CM_scores.pdf',
  D   = 'PAB_vs_RV_CM__dot.pdf',
  E   = 'PAB_seurat_dot_CM.pdf',
  F   = 'PAB_RV_CM_mitocarto.pdf',
  G   = c('S6G_mitocarto_volcano.pdf', 'PAB_RV_CM_mitocarto_volcano.pdf'),
  Hi  = 'PAB_CM_Sev_Nor_Mod_Nor_dot.pdf',
  Hii = 'PAB_CM_Sev_Mod_Mod_Nor_dot.pdf',
  I   = 'PAB_Sirius_box.pdf',
  J   = 'S6J_PAB_RV_CM_FAO.pdf',
  Jii = 'S6J_PAB_RV_CM_Failure_program.pdf',
  Jiii= 'S6J_PAB_RV_CM_GR_axis.pdf')
.s6_emit_dir <- if (exists('V52_FIG_DIR')) V52_FIG_DIR else
                './output/Supplementary_Figure_6'
dir.create(.s6_emit_dir, showWarnings = FALSE, recursive = TRUE)
message('Supplementary Figure 6 v57-keyed standardized panels:')
for (.L in names(.s6_v57)) {
  .cands <- .s6_v57[[.L]]
  .src <- NA_character_
  for (.c in .cands) {
    for (.d in c('./output', .s6_emit_dir,
                 './output/Supplementary_Figure_6')) {
      if (file.exists(file.path(.d, .c))) { .src <- file.path(.d, .c); break }
    }
    if (!is.na(.src)) break
  }
  if (!is.na(.src)) {
    file.copy(.src, file.path(.s6_emit_dir,
              sprintf('Figure_S6_panel_%s.pdf', .L)), overwrite = TRUE)
    message('  ', .L, ': OK (', basename(.src), ')')
  } else message('  ', .L, ': MISSING (', .cands[1], ')')
}
