###############################################################################
## Figure 1 (v53 final) — Study overview and bulk RNA-seq landscape of RV failure
##
## Panels (final, derived from new_scripts/Figure_1.png):
##   (A) Cohort schematic: NF/pRV/RVF across 3 platforms
##       Bulk n=29/78/35 ; snRNA-seq n=4/4/3 ; Xenium n=3/3/3
##   (B) PCA of bulk RNA-seq (PC1 19%, PC2 5%) coloured by group
##   (C) snRNA-seq UMAP with cluster annotations + per-cluster GO/Reactome terms
##   (D) Per-module pathway dotplot — rows = GO/Reactome terms,
##       columns = WGCNA modules (M1..MN); colour = -log10 score, size = padj
##   (E) WGCNA module-score violins by NF/pRV/RVF for modules with a significant
##       eigengene shift (either direction) in NF->pRV or pRV->RVF (computed below
##       in mods_use; not hardcoded)
##   (F) Cross-cohort concordance dotplots (RV vs PAH/Senum):
##       RVF↓ — Respiratory ETC / APOBEC3G / Wnt
##       RVF↑ — IFN-α/β, OAS antiviral, TRAF3
##
## Outputs (./output/v52_figures/):
##   Figure_1_panel_B_pca.pdf, Figure_1_panel_C_umap.pdf,
##   Figure_1_panel_C_umap_simple.pdf, Figure_1_panel_C_dendro.pdf,
##   Figure_1_panel_D_GO.pdf, Figure_1_panel_E_violins.pdf,
##   Figure_1_panel_F_enrichr.pdf
##   (Panel A is an Illustrator asset, no PDF written.)
##   Supplementary PCAs: Figure_1_supp_pca_{age,disease,race,pRV_RVF}.pdf
###############################################################################

source('./helper_scripts/_shared_helpers.R')
## Working-repo adaptation: per-figure output dir + standardized names.
V52_FIG_DIR <- './output/Figure_1'
dir.create(V52_FIG_DIR, showWarnings = FALSE, recursive = TRUE)

## ── F1 dependency bootstrap ───────────────────────────────────────────────
## Stage committed F1 dependencies from ./dependencies/Figure_1/ into the
## working ./output locations the reference pipeline expects:
##   fig1_*_cache.rds  -> ./output/   (bypasses blocked biomaRt network +
##                                     non-deterministic WGCNA recompute)
##   bulkRV_TOM.rda    -> ./output/hdWGCNA_TOM/ and ./TOM/  (GetTOM)
## Only copies when the dest is absent, so `rm -rf output` + rerun stays
## reproducible and self-contained (no /Volumes, no network).
.f1_dep <- './dependencies/Figure_1'
if (dir.exists(.f1_dep)) {
  for (.c in list.files(.f1_dep, pattern = '^fig1_.*cache\\.rds$',
                        full.names = TRUE)) {
    .dst <- file.path('./output', basename(.c))
    if (!file.exists(.dst)) file.copy(.c, .dst)
  }
  if (file.exists(file.path(.f1_dep, 'bulkRV_TOM.rda'))) {
    for (.td in c('./output/hdWGCNA_TOM', './TOM')) {
      dir.create(.td, showWarnings = FALSE, recursive = TRUE)
      .dst <- file.path(.td, 'bulkRV_TOM.rda')
      if (!file.exists(.dst))
        file.copy(file.path(.f1_dep, 'bulkRV_TOM.rda'), .dst)
    }
  }
}

parse_ratio <- function(x) {
  parts <- strsplit(as.character(x), "/")
  sapply(parts, function(p) as.numeric(p[1]) / as.numeric(p[2]))
}

COMP_W <- 7
COMP_H <- 5.5
PS <- pub_scales(COMP_W)

## -- Cartoon asset extraction for Figure 1 (run-once, commented) ------------
## Panel A of the published Figure_1.pdf contains a hand-drawn heart-schematic
## cartoon (NF hypertrophic RV vs pRV/RVF dilated RV, plus cohort counts).
## The PDF is 17 cm wide at 300 DPI ≈ 2008 px wide. Panel A occupies the
## top-left corner. Estimated crop (adjust by eye if needed):
##   widthxheight+x_offset+y_offset = 600x500+0+0  (tune empirically)
## Run-once to extract to ./new_scripts/assets/:
# library(magick)
# fig1 <- image_read_pdf("~/Downloads/hdWGCNA_TOM/Manuscripts/Figure_1.pdf", density=300)
# heart_cartoon <- image_crop(fig1, "600x500+0+0")
# image_write(heart_cartoon, "./new_scripts/assets/Figure_1_panel_A_heart.png")

library(reticulate)
library(ggfortify)
library(edgeR)
library(RColorBrewer)
library(EnhancedVolcano)
library(DESeq2)
library(tximport)
library(biomaRt)
library(sva)
library(Seurat)
library(hdWGCNA)
library(forcats)
library(scales)
library(viridis)


#use_python('/Users/ikuz/anaconda3/envs/velocity/bin/python')

bulk <- read.csv('./dependencies/shared/BulkRNA/counts.csv')
meta <- read.csv('./dependencies/shared/BulkRNA/metadata.csv')
toDel <- seq(1,dim(meta)[1],2)
meta <- meta[-toDel,]

.cache_tx2gene <- './output/fig1_tx2gene_cache.rds'
if (file.exists(.cache_tx2gene)) {
  message('Loading cached tx2gene...')
  tx2gene <- readRDS(.cache_tx2gene)
} else {
  mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl",host = "useast.ensembl.org")
  res <- getBM(attributes = c('ensembl_transcript_id_version',
  'ensembl_gene_id',
  'external_transcript_name',
  'external_gene_name'),
  filters = 'ensembl_transcript_id_version',
  values = bulk[,1],
  mart = mart)
  tx2gene <- res[,c(1,4)]
  tx2gene <- tx2gene[tx2gene$external_gene_name != '',]
  saveRDS(tx2gene, .cache_tx2gene)
}



path = './dependencies/shared/BulkRNA/30-238740824'
files1 <- list.files(path,pattern = "\\.h5$",recursive=TRUE)
files1<-paste0(path,'/',files1)
path = './dependencies/shared/BulkRNA/30-196345105'
files2 <- list.files(path,pattern = "\\.h5$",recursive=TRUE)
files2<-paste0(path,'/',files2)
files <- c(files1,files2)

.cache_txi <- './output/fig1_tximport_cache.rds'
if (file.exists(.cache_txi)) {
  message('Loading cached tximport...')
  txi.kallisto <- readRDS(.cache_txi)
} else {
  txi.kallisto <- tximport(files, type = "kallisto", txOut = FALSE,tx2gene=tx2gene)
  saveRDS(txi.kallisto, .cache_txi)
}
#txi.kallisto <- tximport(files, type = "kallisto", txOut = FALSE,tx2gene=tx2gene,countsFromAbundance="lengthScaledTPM")


bulk <- txi.kallisto$counts
colnames(bulk) <- sapply(strsplit(sapply(strsplit(files,'P0'),'[[',2),'_R1'),'[[',1)

subjects <- colnames(bulk)

subject_map = meta[,5]
category_map = meta[,7]
disease_map = meta[,8]
category_disease_map = meta[,9]
subdisease_map = meta[,10]
sex_map = meta[,11]
age_map = meta[,12]
race_map = meta[,13]
prep_batch_map = meta[,28]
RAP_map = meta[,34]
PAS_map = meta[,35]
PAD_map = meta[,36]
PAM_map = meta[,37]
PCWP_map = meta[,38]
CI_map = meta[,39]
RAP_PCWP_map = meta[,40]
pvri_map = meta[,41]
wt_map = meta[,42]
ht_map = meta[,43]
BSA_map = meta[,44]
BMI_map = meta[,45]
thyroid_map = meta[,46]
pacer_map = meta[,47]
RIN_map = meta[,25]
sc_map = meta[,18]



category <- category_map[match(subjects,subject_map)]
category_disease <- category_disease_map[match(subjects,subject_map)]
disease <- disease_map[match(subjects,subject_map)]
disease <- disease_map[match(subjects,subject_map)]
subdisease <- subdisease_map[match(subjects,subject_map)]
sex <- sex_map[match(subjects,subject_map)]
age <- age_map[match(subjects,subject_map)]
race <- race_map[match(subjects,subject_map)]
prep_batch <- prep_batch_map[match(subjects,subject_map)]
BSA <- BSA_map[match(subjects,subject_map)]
BMI <- BMI_map[match(subjects,subject_map)]
thyroid <- thyroid_map[match(subjects,subject_map)]
pacer <- pacer_map[match(subjects,subject_map)]
RAP <- RAP_map[match(subjects,subject_map)]
PAS <- PAS_map[match(subjects,subject_map)]
PAD <- PAD_map[match(subjects,subject_map)]
PAM <- PAM_map[match(subjects,subject_map)]
PCWP <- PCWP_map[match(subjects,subject_map)]
CI <- CI_map[match(subjects,subject_map)]
RAP_PCWP <- RAP_PCWP_map[match(subjects,subject_map)]
WT <- wt_map[match(subjects,subject_map)]
HT <- ht_map[match(subjects,subject_map)]
RIN <- RIN_map[match(subjects,subject_map)]
sc <- sc_map[match(subjects,subject_map)]
pvri <- pvri_map[match(subjects,subject_map)]



genes <- rownames(bulk)
bulk <- t(bulk)
batch = rep(0,142)
batch[c(which(subjects == 1196),which(subjects == 1482),which(subjects == 1608),which(subjects == 1684),which(subjects == 1690))]=1

bulk.meta <- data.frame(subject=factor(subjects),category=factor(category),category_disease=factor(category_disease),disease=factor(disease),subdisease=factor(subdisease),sex=factor(sex),age=age,race=factor(race),batch=factor(batch),prep_batch=factor(prep_batch),BSA=BSA,BMI=BMI,thyroid=factor(thyroid),pacer=factor(pacer),RAP=RAP,PAS=PAS,PAD=PAD,PAM=PAM,PCWP=PCWP,CI=CI,RAP_PCWP=RAP_PCWP,WT=WT,HT=HT,RIN=RIN,sc=factor(sc),pvri=pvri)
bulk.meta$thyroid[is.na(bulk.meta$thyroid)] <- 'N'
bulk.meta$pacer[is.na(bulk.meta$pacer)] <- 'N'


#ddsSE <- DESeqDataSetFromTximport(txi.kallisto,bulk.meta,design=~category+sex+age+race+batch+prep_batch+BSA+thyroid+pacer+WT+HT+sc+BMI)

ddsSE <- DESeqDataSetFromTximport(txi.kallisto,bulk.meta,design=~category+batch+prep_batch)


ddsSE <- estimateSizeFactors(ddsSE)
idx <- rowSums( counts(ddsSE, normalized=TRUE) >= 5 ) >= 3
ddsSE <- ddsSE[idx,]


normalized_counts <- counts(ddsSE, normalized=TRUE)

mod  <- model.matrix(~category+sex+age+race+batch+prep_batch+BSA+thyroid+pacer+WT+HT+BMI, colData(ddsSE))
mod0 <- model.matrix(~ sex+age+race+batch+prep_batch+BSA+thyroid+pacer+WT+HT+BMI, colData(ddsSE))
.cache_svseq <- './output/fig1_svaseq_cache.rds'
if (file.exists(.cache_svseq)) {
  message('Loading cached svaseq...')
  svseq <- readRDS(.cache_svseq)
} else {
  svseq <- svaseq(normalized_counts, mod, mod0)
  saveRDS(svseq, .cache_svseq)
}


#mod  <- model.matrix(~category+batch+prep_batch, colData(ddsSE))
#mod0 <- model.matrix(~ category, colData(ddsSE))
#svseq <- svaseq(normalized_counts, mod, mod0)



bulk.meta$SV1 <- svseq$sv[,1]
bulk.meta$SV2 <- svseq$sv[,2]
bulk.meta$SV3 <- svseq$sv[,3]
bulk.meta$SV4 <- svseq$sv[,4]
bulk.meta$SV5 <- svseq$sv[,5]
bulk.meta$SV6 <- svseq$sv[,6]
bulk.meta$SV7 <- svseq$sv[,7]
bulk.meta$SV8 <- svseq$sv[,8]
bulk.meta$SV9 <- svseq$sv[,9]
bulk.meta$SV10 <- svseq$sv[,10]
bulk.meta$SV11 <- svseq$sv[,11]
bulk.meta$SV12 <- svseq$sv[,12]
bulk.meta$SV13 <- svseq$sv[,13]
bulk.meta$SV14 <- svseq$sv[,14]
bulk.meta$SV15 <- svseq$sv[,15]
bulk.meta$SV16 <- svseq$sv[,16]
bulk.meta$SV17 <- svseq$sv[,17]
bulk.meta$SV18 <- svseq$sv[,18]
bulk.meta$SV19 <- svseq$sv[,19]
bulk.meta$SV20 <- svseq$sv[,20]
bulk.meta$SV21 <- svseq$sv[,21]


ddsSE <- DESeqDataSetFromTximport(txi.kallisto,bulk.meta,design=~category+SV1+SV2+SV3+SV4+SV5+SV6+SV7+SV8+SV9+SV10+SV11+SV12+SV13+SV14+SV15+SV16+SV17+SV18+SV19+SV20+SV21)


ddsSE <- estimateSizeFactors(ddsSE)
idx <- rowSums( counts(ddsSE, normalized=TRUE) >= 5 ) >= 3
ddsSE <- ddsSE[idx,]

normalized_counts <- counts(ddsSE, normalized=TRUE)



.cache_deseq <- './output/fig1_deseq_fitted_cache.rds'
if (file.exists(.cache_deseq)) {
  message('Loading cached DESeq fitted object...')
  ddsSE <- readRDS(.cache_deseq)
} else {
  ddsSE <- DESeq(ddsSE)
  saveRDS(ddsSE, .cache_deseq)
}

.cache_vst <- './output/fig1_vst_batchcorrected_cache.rds'
if (file.exists(.cache_vst)) {
  message('Loading cached VST batch-corrected object...')
  vstSE <- readRDS(.cache_vst)
} else {
  vstSE <- vst(ddsSE,blind = FALSE)
  mat <- assay(vstSE)
  mm <- model.matrix(~category, colData(vstSE))
  mat <- limma::removeBatchEffect(mat, covariates=colData(vstSE)[,27:47], design=mm)
  assay(vstSE) <- mat
  saveRDS(vstSE, .cache_vst)
}

.cache_lfc_nf_prv <- './output/fig1_lfcshrink_nf_vs_prv_cache.rds'
if (file.exists(.cache_lfc_nf_prv)) {
  message('Loading cached lfcShrink NF vs pRV...')
  nf.vs.prv <- readRDS(.cache_lfc_nf_prv)
} else {
  nf.vs.prv <- lfcShrink(ddsSE,contrast=c('category','NF','pRV'), type="ashr")
  saveRDS(nf.vs.prv, .cache_lfc_nf_prv)
}

.cache_lfc_nf_rvf <- './output/fig1_lfcshrink_nf_vs_rvf_cache.rds'
if (file.exists(.cache_lfc_nf_rvf)) {
  message('Loading cached lfcShrink NF vs RVF...')
  nf.vs.rvf <- readRDS(.cache_lfc_nf_rvf)
} else {
  nf.vs.rvf <- lfcShrink(ddsSE,contrast=c('category','NF','RVF'), type="ashr")
  saveRDS(nf.vs.rvf, .cache_lfc_nf_rvf)
}

.cache_lfc_prv_rvf <- './output/fig1_lfcshrink_prv_vs_rvf_cache.rds'
if (file.exists(.cache_lfc_prv_rvf)) {
  message('Loading cached lfcShrink pRV vs RVF...')
  prv.vs.rvf <- readRDS(.cache_lfc_prv_rvf)
} else {
  prv.vs.rvf <- lfcShrink(ddsSE,contrast=c('category','pRV','RVF'), type="ashr")
  saveRDS(prv.vs.rvf, .cache_lfc_prv_rvf)
}

#saveRDS(vstSE,'./output/bulkRNAseq_vst.rds')

#Cast to Seurat to use hdWGCNA backend for consistency

mat <- as(assay(vstSE),'dgCMatrix')
colnames(mat)<-subjects
bulk <- CreateAssayObject(counts=mat, meta.data=bulk.meta, assay = "RNA")

bulk <-  CreateSeuratObject(bulk)
seurat_obj <- bulk
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- AddMetaData(seurat_obj, bulk.meta)
#saveRDS(seurat_obj,'./output/RV_bulkRNASeq_seurat.rds')

#######################################
#############  FIGURE 1B  #############
#######################################

p_1B <- plotPCA(vstSE,intgroup=c("category"),ntop=19355) + theme_classic() + theme(axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.title.x=element_text(size=24),axis.title.y=element_text(size=24),legend.title=element_text(size=24),legend.text=element_text(size=24),legend.key.size = PS$legend_key) + labs(color='Disease',shape='Sex') + geom_point(size=PS$scatter_pt) + scale_colour_disease()
save_figure(p_1B, 'Figure_1_panel_B_pca.pdf', width=10, height=3)

## Supplementary PCA: pRV vs RVF only (used in Supp Fig 1, not main Fig 1B)
p_1B2 <- plotPCA(vstSE[,category %in% c('pRV','RVF')],intgroup=c("category"),ntop=19355) + theme_classic() + theme(axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.title.x=element_text(size=24),axis.title.y=element_text(size=24),legend.title=element_text(size=24),legend.text=element_text(size=24),legend.key.size = PS$legend_key) + labs(color='Disease',shape='Sex') + geom_point(size=PS$scatter_pt) + scale_colour_disease()
save_figure(p_1B2, 'Figure_1_supp_pca_pRV_RVF.pdf', width=10, height=3)

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = viridis::mako(10))
save_figure(plotPCA(vstSE,intgroup=c("age"),ntop=19355) + theme_classic() + theme(axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.title.x=element_text(size=24),axis.title.y=element_text(size=24),legend.title=element_text(size=24),legend.text=element_text(size=24)) + labs(color='Age',shape='Sex')+ sc,
  'Figure_1_supp_pca_age.pdf', width=10, height=6)

save_figure(plotPCA(vstSE,intgroup=c("disease"),ntop=19355) + theme_classic() + theme(axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.title.x=element_text(size=24),axis.title.y=element_text(size=24),legend.title=element_text(size=24),legend.text=element_text(size=24)) + labs(color='Disease',shape='Etiology'),
  'Figure_1_supp_pca_disease.pdf', width=10, height=6)

save_figure(plotPCA(vstSE,intgroup=c("race"),ntop=19355)  + theme_classic() + theme(axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.title.x=element_text(size=24),axis.title.y=element_text(size=24),legend.title=element_text(size=24),legend.text=element_text(size=24)) + labs(color='Disease',shape='Race'),
  'Figure_1_supp_pca_race.pdf', width=10, height=6)


#######################################
##############  FIGURE 1C  ############
#######################################

seurat_obj <- readRDS('./dependencies/shared/RV_bulkRNASeq_seurat.rds')

.cache_wgcna <- './output/fig1_wgcna_network_cache.rds'
if (file.exists(.cache_wgcna)) {
  message('Loading cached WGCNA network object...')
  seurat_obj <- readRDS(.cache_wgcna)
} else {
  seurat_obj <- SetupForWGCNA(
        seurat_obj,
        gene_select = "fraction",
        fraction = 0.05,
        wgcna_name = 'bulkRV'
      )

  seurat_obj <- SetDatExpr(
    seurat_obj,
    mat = t(seurat_obj[['RNA']]$counts),
  )

  seurat_obj <- TestSoftPowers(seurat_obj)

  seurat_obj <- ConstructNetwork(
      seurat_obj,
      tom_dir ='./output/hdWGCNA_TOM',
      tom_name='bulkRV',
      overwrite_tom=TRUE,
      mergeCutHeight=0.15,
      soft_power=9,
  )

  # compute the MEs and kMEs
  seurat_obj <- ModuleEigengenes(seurat_obj)
  seurat_obj <- ModuleConnectivity(seurat_obj)
  saveRDS(seurat_obj, .cache_wgcna)
}


# get MEs from seurat object
.cache_MEs <- './output/fig1_module_eigengenes_cache.rds'
if (file.exists(.cache_MEs)) {
  message('Loading cached module eigengenes...')
  MEs <- readRDS(.cache_MEs)
} else {
  MEs <- GetMEs(seurat_obj)
  saveRDS(MEs, .cache_MEs)
}
mods <- colnames(MEs); mods <- mods[mods != 'grey']

# add MEs to Seurat meta-data for plotting:
meta <- seurat_obj@meta.data
seurat_obj@meta.data <- cbind(meta, MEs)
#saveRDS(seurat_obj,'./output/RV_bulkRNASeq_seurat.rds')


# plot with Seurat's DotPlot function
p <- DotPlot(seurat_obj, features=mods, group.by = 'category')


#Write module assignments to file
data_dir <- "./output/"
fig_dir <- './output/hdWGCNA/'

modules <- GetModules(seurat_obj) %>% subset(module != 'grey')

write.csv(modules, row.names=FALSE, quote=FALSE, file=paste0(data_dir, 'bulk_heart_modules.csv'))

#Enrichments
library(enrichR)

dbs <-c('GO_Biological_Process_2023','GO_Cellular_Component_2023','GO_Molecular_Function_2023','LINCS_L1000_Chem_Pert_up','LINCS_L1000_Chem_Pert_down', 'WikiPathway_2021_Human', 'KEGG_2021_Human')

# compute GO terms (fall back to cached TSV if Enrichr unreachable):
enrich_list <- list()
.cache_enrichr1 <- './output/fig1_enrichr_all_dbs_cache.rds'
if (file.exists(.cache_enrichr1)) {
  message('Loading cached RunEnrichr (all dbs)...')
  seurat_obj <- readRDS(.cache_enrichr1)
} else {
  seurat_obj <- tryCatch(RunEnrichr(seurat_obj, dbs=dbs), error = function(e) {
    message('Enrichr API unreachable, loading cached results: ', e$message)
    seurat_obj
  })
  saveRDS(seurat_obj, .cache_enrichr1)
}
enrichr_df <- tryCatch(GetEnrichrTable(seurat_obj) %>% subset(P.value < 0.05), error = function(e) {
  message('Loading cached enrichr from disk')
  read.delim(paste0(data_dir, 'bulk_heart_enrichr.tsv'), stringsAsFactors = FALSE)
})

write.table(enrichr_df, quote=FALSE, sep='\t', row.names=FALSE, file=paste0(data_dir, 'bulk_heart_enrichr.tsv'))


#Network visualizations and module eigengene plots
library(reshape2)
library(igraph)


# plot the dendrograms
pdf(file.path(V52_FIG_DIR, 'Figure_1_panel_C_dendro.pdf'), height=3, width=6)
PlotDendrogram(seurat_obj, main=' Dendrogram')
dev.off()

# module network plot
ModuleNetworkPlot(
  seurat_obj,
  mods = "all",
  outdir = paste0(data_dir, '_hubNetworks/'),
)

# module umap plot


# compute the module UMAPs
## Panel C UMAP must MATCH the published figure.  The canonical
## RV_bulkRNASeq_seurat.rds dependency already carries the EXACT frozen
## ModuleUMAP that produced the published panel.  RunModuleUMAP is a
## STOCHASTIC UMAP, so recomputing it (what the old ./output cache did)
## yields a rotated/flipped/distorted embedding that does NOT match the
## original.  Plot panel C from the canonical object's frozen embedding;
## only recompute (seeded) if that object truly lacks one.
.objC <- readRDS('./dependencies/shared/RV_bulkRNASeq_seurat.rds')
.muC  <- tryCatch(GetModuleUMAP(.objC), error = function(e) NULL)
if (is.null(.muC) || !is.data.frame(.muC) || nrow(.muC) == 0) {
  message('Panel C: canonical object lacks ModuleUMAP — recomputing (seeded)')
  set.seed(12345)
  .objC <- RunModuleUMAP(
    .objC,
    n_hubs = 5,
    n_neighbors=10,
    min_dist=0.3,
    spread=2,
    target_weight=0.1,
    supervised=TRUE
  )
} else {
  message('Panel C: using canonical frozen ModuleUMAP (',
          nrow(.muC), ' hub-gene rows) from RV_bulkRNASeq_seurat.rds')
}



# get the hub gene UMAP table from the canonical (frozen-embedding) object
umap_df <- GetModuleUMAP(.objC)

# plot with ggplot
plot_df <- umap_df

# compute coordinates for cluster labels
centroid_df <- data.frame()
for(cur_cluster in unique(plot_df[['module']])){
    cur_meta <- plot_df[plot_df[['module']] == cur_cluster,]
    df <- data.frame(
      cluster = cur_cluster,
      UMAP1 = mean(cur_meta$UMAP1),
      UMAP2 = mean(cur_meta$UMAP2)
    )
  centroid_df <- rbind(centroid_df, df)
}

# plot with ggplot
p <- ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
    ggrastr::rasterise(geom_point(
     color=umap_df$color,
     size=umap_df$kME*2
   ), dpi=500, scale=0.5) +
    umap_theme() +
    theme(
      plot.margin = margin(0,0,0,0),
      plot.title = element_text(hjust=0.5)
    ) + ggtitle('Bulk RNASeq') +
    ggrepel::geom_text_repel(data = centroid_df, label=centroid_df$cluster, color='black', max.overlaps=Inf, size=PS$text_mm, family=FONT_FAMILY)


  save_figure(p, 'Figure_1_panel_C_umap_simple.pdf', width=5, height=5)

  hub_genes <- GetHubGenes(.objC, 3)

  # add annotation
  anno_genes <- hub_genes$gene_name
  plot_df$anno <- ifelse(plot_df$gene %in% anno_genes, umap_df$gene, '')

  plot_df_anno <- subset(plot_df, anno != '')
  p <-  plot_df %>%
    ggplot(aes(x=UMAP1, y=UMAP2, color=module)) +
    ggrastr::rasterise(
      geom_point(
        inherit.aes=FALSE,
        data=plot_df,
        aes(x=UMAP1, y=UMAP2, color=module),
        color=plot_df$color,
        size=plot_df$kME*2,
      ), dpi=500, dpi_scale=0.5) +
    geom_point(
      inherit.aes = FALSE,
      data = plot_df_anno,
      shape=21, color='black',
      fill=plot_df_anno$color,
      size=plot_df_anno$kME*2,
      aes(x=UMAP1, y=UMAP2, fill=module)
    ) +
    # add labels
    ggrepel::geom_text_repel(data = centroid_df, label=centroid_df$cluster, color='black', max.overlaps=Inf, size=PS$text_mm, family=FONT_FAMILY, fontface='bold') +
    geom_text_repel(label=plot_df$anno, max.overlaps=Inf, color='black', fontface='italic', size=PS$text_mm, family=FONT_FAMILY) +
    umap_theme() + NoLegend() +
    coord_equal() +
    theme(
      plot.margin = margin(0,0,0,0)
    )

  p_1C <- p
  save_figure(p_1C, 'Figure_1_panel_C_umap.pdf', width=8, height=8)

  # plot with igraph (skip if globals too large for future framework)
  tryCatch({
    pdf(paste0(data_dir,'_hubgene_umap_igraph.pdf'), width=12, height=12)
    p <- ModuleUMAPPlot(
      .objC,
      edge.alpha=0.5,
      sample_edges=TRUE,
      keep_grey_edges=FALSE,
      edge_prop=0.075,
      label_hubs=3,
      return_graph = TRUE,
    )
    dev.off()
  }, error = function(e) {
    try(dev.off(), silent = TRUE)
    message('ModuleUMAPPlot skipped: ', e$message)
  })

#######################################
#########  FIGURE 1C (cont) ###########
#######################################


dbs <- c('GO_Biological_Process_2023','GO_Cellular_Component_2023','GO_Molecular_Function_2023')

# perform enrichment tests (fall back gracefully if Enrichr unreachable)
.cache_enrichr2 <- './output/fig1_enrichr_go_only_cache.rds'
if (file.exists(.cache_enrichr2)) {
  message('Loading cached RunEnrichr (GO only)...')
  seurat_obj <- readRDS(.cache_enrichr2)
} else {
  seurat_obj <- tryCatch(RunEnrichr(
    seurat_obj,
    dbs=dbs,
    max_genes = 100
  ), error = function(e) {
    message('Enrichr API unreachable (2nd call), continuing with prior results: ', e$message)
    seurat_obj
  })
  saveRDS(seurat_obj, .cache_enrichr2)
}

# retrieve the output table
enrich_df <- tryCatch(GetEnrichrTable(seurat_obj), error = function(e) {
  message('GetEnrichrTable failed, loading cached enrichr from disk')
  read.delim('./output/bulk_heart_enrichr.tsv', stringsAsFactors = FALSE)
})

tryCatch(EnrichrBarPlot(
  seurat_obj,
  outdir = "./output/bulk_enrichr_plot_all_genes",
  n_terms = 5,
  plot_size = c(5,4),
  logscale=TRUE
), error = function(e) message('EnrichrBarPlot skipped: ', e$message))


#######################################
#############  FIGURE 1D  #############
#######################################

## CANONICAL 29-module decomposition.  GetModules(seurat_obj) returns the
## DRIFTED 32-module network baked into the committed enrichr cache, which
## scrambled the M1..M29 staircase (Panel D "broken").  Authoritative source
## = the legacy-written dependencies/shared/bulk_heart_modules.csv (29 mods).
.bhm <- read.csv('./dependencies/shared/bulk_heart_modules.csv',
                 stringsAsFactors = FALSE)
.bhm <- subset(.bhm, module != 'grey')
modules <- data.frame(
  gene_name = .bhm$gene_name,
  module    = factor(.bhm$module, levels = unique(.bhm$module)),
  color     = .bhm$color,
  stringsAsFactors = FALSE)
color_df <- modules %>% dplyr::select(module, color) %>% dplyr::distinct() %>%
  dplyr::rename(group = module, colour = color)
mods <- as.character(unique(modules$module))
mods <- mods[mods != 'grey']

# helper function to wrap text
wrapText <- function(x, len) {
    sapply(x, function(y) paste(strwrap(y, len), collapse = "\n"), USE.NAMES = FALSE)
}

## Compute GO_BP enrichment ON THE CANONICAL 29-MODULE GENE SETS (top 100
## genes/module by that module's kME == legacy RunEnrichr max_genes=100),
## via live enrichr — NOT the drifted 32-module GetEnrichrTable cache.
## New cache name so the stale 32-module table is never reused.
.cache_pD_enr <- './output/fig1_panelD_canonical_GO_enrichr.rds'
if (file.exists(.cache_pD_enr)) {
  message('Loading cached canonical Panel-D GO enrichr...')
  combined_output <- readRDS(.cache_pD_enr)
} else {
  message('Computing canonical Panel-D GO_BP enrichr over 29 modules...')
  .mod_levels <- as.character(unique(modules$module))
  .co_list <- lapply(.mod_levels, function(.cc) {
    .kcol <- paste0('kME_', .cc)
    .sub  <- .bhm[.bhm$module == .cc, , drop = FALSE]
    if (.kcol %in% colnames(.sub))
      .sub <- .sub[order(.sub[[.kcol]], decreasing = TRUE), ]
    .genes <- head(unique(.sub$gene_name), 100)
    if (length(.genes) < 3) return(NULL)
    .e <- tryCatch(enrichR::enrichr(.genes, 'GO_Biological_Process_2023'),
                   error = function(e) NULL)
    Sys.sleep(1)
    if (is.null(.e) || is.null(.e[[1]]) || nrow(.e[[1]]) == 0) return(NULL)
    .t <- .e[[1]]
    .t <- .t[order(.t$P.value), ]
    .t$db     <- 'GO_Biological_Process_2023'
    .t$module <- .cc
    .t
  })
  combined_output <- do.call(rbind, .co_list)
  if (is.null(combined_output) || nrow(combined_output) == 0)
    stop('Panel D: canonical GO enrichr produced no rows')
  saveRDS(combined_output, .cache_pD_enr)
}
selected_terms <- subset(combined_output,db=="GO_Biological_Process_2023")

# subset selected terms
selected_terms <- subset(selected_terms, P.value < 0.05)
idx_top_1 <- match(unique(selected_terms$module),selected_terms$module)
idx_top_3 <- sort(c(idx_top_1,idx_top_1+1,idx_top_1+2))

selected_terms<-selected_terms[idx_top_3,]

# remove GO Term ID
library(stringr)
selected_terms$Term <- str_replace(selected_terms$Term, " \\s*\\([^\\)]+\\)", "")

key_terms <- read.csv('./dependencies/Figure_1/bulkRNA_GOterms_ofinterest.csv')
selected_terms <- subset(selected_terms,Term %in% key_terms[[1]])



mapping <- labels2colors(1:100)

## Full ordered module set so EVERY module M1..M29 gets an x-axis column
## (user: one entry per module even when it has no curated/enriched term).
.mod_color_lv <- mods[order(match(mods, mapping))]            # color names, M-order
.mod_M_lv     <- paste0('M', match(.mod_color_lv, mapping))   # M1..M29, same order

selected_terms$group <- factor(
  as.character(selected_terms$module),
  levels = .mod_color_lv
)


# set max pval
quantile(-log(selected_terms$Adjusted.P.value), 0.95)
max_p <- 10

selected_terms$logp <- -log(selected_terms$Adjusted.P.value)
selected_terms$logp <- ifelse(selected_terms$logp > max_p, max_p, selected_terms$logp)



selected_terms <- selected_terms %>%
  arrange(group)

selected_terms$wrap <- wrapText(selected_terms$Term, 45)

selected_terms$Term <- factor(
  as.character(selected_terms$Term),
  levels = rev(unique(as.character(selected_terms$Term)))
)

selected_terms$wrap <- factor(
  as.character(selected_terms$wrap),
  levels = rev(unique(as.character(selected_terms$wrap)))
)


library(viridis)

## Colorbar spans ALL 29 non-grey modules (one tile each), not only the
## modules that had a curated term, so every module M1..M29 shows on x.
color_df <- modules %>% subset(module != 'grey') %>%
  dplyr::select(module, color) %>% dplyr::distinct()
color_df$group  <- factor(paste0('M', match(as.character(color_df$module), mapping)),
                          levels = .mod_M_lv)
color_df$colour <- color_df$color
color_df <- color_df[order(color_df$group), ]
# make the colorbar as its own heatmap
color_df$var <- 1

c_vect <- color_df$colour
names(c_vect) <- as.character(color_df$group)


selected_terms$group_num <- factor(
  paste0('M', match(as.character(selected_terms$module), mapping)),
  levels = .mod_M_lv
)



# GO Term dot plot
p <- selected_terms %>%
  ggplot(aes(x = group, y = wrap, color =logp, size=log(Combined.Score))) +
  geom_point() +
  scale_x_discrete(drop = FALSE) +
  scale_color_stepsn(colors=rev(magma(256))) +
  scale_size(range = PS$dot_range) +
  RotatedAxis() + xlab('') + ylab('') +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_rect(size=1, color='black', fill=NA),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(0,0,0,0),
    panel.grid = element_line(size=0.25, color='lightgrey'),
    legend.key.size = PS$legend_key
  )



colorbar <- color_df %>%
  ggplot(aes(x=group, y=var, fill=group)) +
  geom_tile() +
  scale_x_discrete(drop = FALSE) +
  scale_fill_manual(values=c_vect) +
  coord_equal() +
  NoLegend() + RotatedAxis() +
  theme(
    plot.title=element_blank(),
    axis.line=element_blank(),
    axis.ticks.y =element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    plot.margin=margin(0,0,0,0),
  )


p_1D <- p / colorbar
save_figure(p_1D, 'Figure_1_panel_D_GO.pdf', width=13, height=11.7)



#######################################
#############  FIGURE 1E  #############
#######################################

source('./helper_scripts/spatial_functions.R')



group1 <- seurat_obj@meta.data %>% subset(category == 'NF') %>% rownames
group2 <- seurat_obj@meta.data %>% subset(category == 'pRV') %>% rownames
group3 <- seurat_obj@meta.data %>% subset(category == 'RVF') %>% rownames

.cache_dme_prv_rvf <- './output/fig1_dme_prv_vs_rvf_cache.rds'
if (file.exists(.cache_dme_prv_rvf)) {
  message('Loading cached FindDMEs pRV vs RVF...')
  DMEs_prv_vs_rvf <- readRDS(.cache_dme_prv_rvf)
} else {
  DMEs_prv_vs_rvf <- FindDMEs(
    seurat_obj,
    barcodes1 = group2,
    barcodes2 = group3,
    test.use='wilcox',
    pseudocount.use=0,
  )
  saveRDS(DMEs_prv_vs_rvf, .cache_dme_prv_rvf)
}

.cache_dme_nf_rvf <- './output/fig1_dme_nf_vs_rvf_cache.rds'
if (file.exists(.cache_dme_nf_rvf)) {
  message('Loading cached FindDMEs NF vs RVF...')
  DMEs_nf_vs_rvf <- readRDS(.cache_dme_nf_rvf)
} else {
  DMEs_nf_vs_rvf <- FindDMEs(
    seurat_obj,
    barcodes1 = group1,
    barcodes2 = group3,
    test.use='wilcox',
    pseudocount.use=0,
  )
  saveRDS(DMEs_nf_vs_rvf, .cache_dme_nf_rvf)
}

.cache_dme_nf_prv <- './output/fig1_dme_nf_vs_prv_cache.rds'
if (file.exists(.cache_dme_nf_prv)) {
  message('Loading cached FindDMEs NF vs pRV...')
  DMEs_nf_vs_prv <- readRDS(.cache_dme_nf_prv)
} else {
  DMEs_nf_vs_prv <- FindDMEs(
    seurat_obj,
    barcodes1 = group1,
    barcodes2 = group2,
    test.use='wilcox',
    pseudocount.use=0,
  )
  saveRDS(DMEs_nf_vs_prv, .cache_dme_nf_prv)
}

## v59 panel-E selection: keep modules with a significant module-eigengene shift
## (either direction; FindDMEs Wilcoxon, Bonferroni p_val_adj < 0.05) in at least
## one ADJACENT transition -- NF->pRV or pRV->RVF.  (Previously unioned pRV-vs-RVF
## with NF-vs-RVF listed twice, which omitted the NF->pRV transition.)
mods_use <- unique(c(
  rownames(subset(DMEs_nf_vs_prv,  p_val_adj < 0.05)),
  rownames(subset(DMEs_prv_vs_rvf, p_val_adj < 0.05))))


# seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)
.cache_MEs2 <- './output/fig1_module_eigengenes_renamed_cache.rds'
if (file.exists(.cache_MEs2)) {
  message('Loading cached module eigengenes (Fig 1E)...')
  MEs <- readRDS(.cache_MEs2)
} else {
  MEs <- GetMEs(seurat_obj)
  saveRDS(MEs, .cache_MEs2)
}
mods <- colnames(MEs); mods <- mods[mods != 'grey']


MEs_rename <- paste0('M',match(colnames(MEs),mapping))
colnames(MEs)<-MEs_rename

seurat_obj@meta.data <- seurat_obj@meta.data[,1:80]
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)
temp = data.frame(dummy=MEs[,1])
temp[1:142,]='1'
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, temp)

mods_use <-  paste0('M',match(mods_use,mapping))
## Order the stacked violins numerically M1 -> M2 -> M3 ... (top to bottom);
## custom_vln stacks `features` in the given order from top down.
mods_use <- mods_use[order(as.integer(sub('^M', '', mods_use)))]


library(ggpubr)
library(patchwork)

p <- custom_vln(
    seurat_obj,
    features = mods_use,
    group.by = 'dummy',
    groups = c('1'),
    add_boxplot=FALSE,
    split.by = 'category',
    selected_split = c('NF','pRV','RVF'),
    split_colors=c('darkorchid', 'grey','royalblue'),
    add_colorbar=FALSE,
    plot_ymin = NA,
    pval_y_adjust=0.7,
  )


p_1E <- p
save_figure(p_1E, 'Figure_1_panel_E_violins.pdf', width=5, height=10)


#######################################
#############  FIGURE 1F  #############
#######################################
## Cross-cohort PAH concordance dotplot (RV vs Senum 2023 PAH cohort).
## Two stacked Reactome dotplots:
##   (top)    Concordant DOWN — genes down in RVF vs NF in BOTH RV and PAH;
##            enriched for respiratory ETC / mitochondrial pathways.
##   (bottom) Concordant UP   — genes up in RVF vs NF in BOTH RV and PAH;
##            enriched for IFN-α/β, OAS antiviral, TRAF3 signaling.
##
## PAH DESeq2 (Control vs Decompensated, ashr-shrunken log2FC) is computed
## from GSE240921 (Senum et al. 2023) the first time the script runs and
## cached as `fig1_pah_deg_cache.rds`. Subsequent runs skip the heavy
## DESeq2 + svaseq step.
##
## Enrichment mirrors the original PAH code (additional_scripts/AnalysisPAH.R
## L466-511): enrichr over the SAME 4-db vector, the GO_Biological_Process_2025
## table (= enriched[[4]]), Adjusted.P.value < 0.05, ranked by Combined.Score,
## top 4 terms, GO-id suffix stripped.

dbs_pah <- c("ChEA_2022", "WikiPathways_2024_Human",
             "Reactome_Pathways_2024", "GO_Biological_Process_2025")

## ── Step 1. PAH DEG: cached pah.nf.vs.rvf log2FC table ──────────────────────
.cache_pah_deg <- './output/fig1_pah_deg_cache.rds'
if (file.exists(.cache_pah_deg)) {
  message('Loading cached PAH DEGs (GSE240921 Control vs Decompensated)...')
  pah.nf.vs.rvf <- readRDS(.cache_pah_deg)
} else {
  message('Computing PAH DEGs from GSE240921 (first-time only)...')
  .pah_raw_path <- './dependencies/shared/GSE240921_processed-data-human.csv'
  .pah_ann_path <- './dependencies/shared/Human.GRCh38.p13.annot.tsv'
  if (file.exists(.pah_raw_path) && file.exists(.pah_ann_path)) {
    pah_raw <- read.csv(.pah_raw_path, sep = ',')
    gene.annot <- read.csv(.pah_ann_path, sep = '\t')
    pah_raw$names <- gene.annot$Symbol[match(pah_raw$id, gene.annot$EnsemblGeneID)]
    pah_raw$names[is.na(pah_raw$names)] <- pah_raw$id[is.na(pah_raw$names)]
    pah_proc <- pah_raw[!duplicated(pah_raw$names), 2:41]
    rownames(pah_proc) <- pah_raw$names[!duplicated(pah_raw$names)]
    subj.names <- colnames(pah_raw)[2:41]
    subj.group <- stringr::str_split_i(subj.names, '_', 1)
    subj.group[subj.group == 'RV.Normal']  <- 'Control'
    subj.group[subj.group == 'RV.Compen']  <- 'Compensated'
    subj.group[subj.group == 'RV.Failing'] <- 'Decompensated'
    coldata_pah <- data.frame(group = factor(subj.group,
                                             levels = c('Control','Compensated','Decompensated')))
    rownames(coldata_pah) <- subj.names
    dds_pah <- DESeqDataSetFromMatrix(countData = pah_proc,
                                      colData   = coldata_pah,
                                      design    = ~ group)
    dds_pah <- estimateSizeFactors(dds_pah)
    keep_pah <- rowSums(counts(dds_pah, normalized = TRUE) >= 5) >= 3
    dds_pah <- dds_pah[keep_pah, ]
    mod_pah  <- model.matrix(~ group, colData(dds_pah))
    mod0_pah <- model.matrix(~ 1, colData(dds_pah))
    sv_pah <- svaseq(counts(dds_pah, normalized = TRUE), mod_pah, mod0_pah)
    coldata_pah$SV1 <- sv_pah$sv[, 1]; coldata_pah$SV2 <- sv_pah$sv[, 2]
    coldata_pah$SV3 <- sv_pah$sv[, 3]; coldata_pah$SV4 <- sv_pah$sv[, 4]
    coldata_pah$SV5 <- sv_pah$sv[, 5]
    dds_pah <- DESeqDataSetFromMatrix(countData = pah_proc[keep_pah, ],
                                      colData   = coldata_pah,
                                      design    = ~ group + SV1 + SV2 + SV3 + SV4 + SV5)
    dds_pah <- DESeq(dds_pah)
    pah.nf.vs.rvf <- lfcShrink(dds_pah,
                               contrast = c('group','Control','Decompensated'),
                               type = 'ashr')
    saveRDS(pah.nf.vs.rvf, .cache_pah_deg)
    rm(dds_pah, pah_raw, pah_proc, gene.annot, sv_pah, mod_pah, mod0_pah,
       keep_pah, coldata_pah, subj.names, subj.group); gc(verbose = FALSE)
  } else {
    message('PAH raw data missing — Panel F will be skipped: ', .pah_raw_path)
    pah.nf.vs.rvf <- NULL
  }
}

## ── Step 2. Concordant gene sets (RV ∩ PAH) ─────────────────────────────────
if (!is.null(pah.nf.vs.rvf)) {
  shared_genes <- intersect(rownames(pah.nf.vs.rvf), rownames(nf.vs.rvf))
  pah.up   <- intersect(rownames(subset(pah.nf.vs.rvf, padj < 0.1 & log2FoldChange >  0.1)), shared_genes)
  pah.down <- intersect(rownames(subset(pah.nf.vs.rvf, padj < 0.1 & log2FoldChange < -0.1)), shared_genes)
  rv.up    <- intersect(rownames(subset(nf.vs.rvf,     padj < 0.1 & log2FoldChange >  0.1)), shared_genes)
  rv.down  <- intersect(rownames(subset(nf.vs.rvf,     padj < 0.1 & log2FoldChange < -0.1)), shared_genes)
  ## Convention: positive log2FC = up in NF (= down in RVF).
  ## Concordant DOWN in RVF = up in NF in both cohorts  = pah.up ∩ rv.up.
  ## Concordant UP   in RVF = down in NF in both        = pah.down ∩ rv.down.
  concordant_down_in_RVF <- intersect(pah.up,   rv.up)
  concordant_up_in_RVF   <- intersect(pah.down, rv.down)
  message(sprintf('Panel F: %d concordant DOWN, %d concordant UP genes (RVF vs NF)',
                  length(concordant_down_in_RVF), length(concordant_up_in_RVF)))

  ## ── Step 3. GO_BP_2025 enrichr on each concordant set (cached) ───────────
  ## Use enriched[[4]] = GO_Biological_Process_2025, exactly as AnalysisPAH.R.
  ## New cache name so the stale Reactome-only cache is NOT reused.
  .cache_enrichr_pah <- './output/fig1_enrichr_pah_concordance_go_cache.rds'
  if (file.exists(.cache_enrichr_pah)) {
    message('Loading cached enrichr (PAH concordance, GO_BP_2025)...')
    enr_pah <- readRDS(.cache_enrichr_pah)
  } else {
    enr_pah <- tryCatch({
      list(down = enrichr(concordant_down_in_RVF, dbs_pah)[[4]],
           up   = enrichr(concordant_up_in_RVF,   dbs_pah)[[4]])
    }, error = function(e) {
      message('Enrichr API unreachable for Panel F: ', e$message); NULL
    })
    if (!is.null(enr_pah)) saveRDS(enr_pah, .cache_enrichr_pah)
  }

  ## ── Step 4. Build the two stacked dotplots ───────────────────────────────
  if (!is.null(enr_pah)) {
    .build_pah_dot <- function(enr_tbl, panel_title, topN = 4) {
      ## Match AnalysisPAH.R: padj<0.05, rank by Combined.Score, top 4,
      ## strip the " (GO:#######)" suffix from GO_BP_2025 term names.
      enr_tbl <- subset(enr_tbl, Adjusted.P.value < 0.05)
      enr_tbl <- enr_tbl[order(enr_tbl$Combined.Score, decreasing = TRUE), ]
      enr_tbl <- enr_tbl[rev(seq_len(min(topN, nrow(enr_tbl)))), ]
      ggplot(enr_tbl,
             aes(x = Combined.Score, y = fct_inorder(Term),
                 color = as.numeric(Adjusted.P.value),
                 size  = parse_ratio(Overlap))) +
        geom_point() +
        scale_size(range = PS$dot_range) +
        xlab('Combined Score') + ylab('Term') +
        labs(color = 'P value', size = 'Overlap') +
        ggtitle(panel_title) +
        scale_y_discrete(labels = fct_inorder(
          wrapText(sapply(strsplit(enr_tbl$Term, ' \\(GO'), `[`, 1), 35))) +
        theme_classic() +
        theme(axis.text = element_text(colour = 'black'),
              legend.key.size = PS$legend_key) +
        scale_color_stepsn(colors = rev(magma(256)))
    }
    p_1F_down <- .build_pah_dot(enr_pah$down, 'Concordant DOWN')
    p_1F_up   <- .build_pah_dot(enr_pah$up,   'Concordant UP')

    p_1F <- (p_1F_down / p_1F_up) +
      patchwork::plot_layout(guides = 'collect') &
      theme(legend.position = 'right',
            plot.title = element_text(size = 10, face = 'bold'))
    save_figure(p_1F, 'Figure_1_panel_F_enrichr.pdf', width = 6, height = 5)
  } else {
    message('Panel F skipped — enrichr cache missing and API unavailable.')
  }
} else {
  message('Panel F skipped — PAH DEG data unavailable.')
}


###############################################################################
##  Panel A — Illustrator-drawn cohort/heart schematic (asset)
##  Falls back to a dashed placeholder if the PNG asset is missing.
###############################################################################
p_1A <- insert_asset('Figure_1_panel_A_heart.png',
                     label = 'NF (n=29)  pRV (n=78)  RVF (n=35)')

message('Figure 1 (v53) per-panel PDFs written to ', V52_FIG_DIR)
