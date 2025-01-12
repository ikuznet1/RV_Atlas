library(reticulate)
library(ggfortify)
library(edgeR)
library(RColorBrewer)
library(EnhancedVolcano)
library(DESeq2)
library(tximport)
library(biomaRt)
library("sva")


use_python('/Users/ikuz/anaconda3/envs/velocity/bin/python')

bulk <- read.csv('/Volumes/RV_RNAseq/RV_snRNAseq/Final_Analysis/BulkRNA/counts.csv')
meta <- read.csv('/Volumes/RV_RNAseq/RV_snRNAseq/Final_Analysis/BulkRNA/metadata.csv')
toDel <- seq(1,dim(meta)[1],2)
meta <- meta[-toDel,]

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl",host = "useast.ensembl.org")

res <- getBM(attributes = c('ensembl_transcript_id_version',                              
'ensembl_gene_id',                              
'external_transcript_name',                           
'external_gene_name'),              
filters = 'ensembl_transcript_id_version',               
values = bulk[,1],              
mart = mart)
tx2gene <- res[,c(1,4)]


path = '/Volumes/RV_RNAseq/RV_RNAseq/30-238740824'
files1 <- list.files(path,pattern = "\\.h5$",recursive=TRUE)
files1<-paste0(path,'/',files1)
path = '/Volumes/RV_RNAseq/RV_RNAseq/30-196345105/trimmed'
files2 <- list.files(path,pattern = "\\.h5$",recursive=TRUE)
files2<-paste0(path,'/',files2)
files <- c(files1,files2)

txi.kallisto <- tximport(files, type = "kallisto", txOut = FALSE,tx2gene=tx2gene)


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


ddsSE <- DESeqDataSetFromTximport(txi.kallisto,bulk.meta,design=~category+sex+age+race+batch+prep_batch+BSA+thyroid+pacer+WT+HT+sc+BMI)

ddsSE <- DESeqDataSetFromTximport(txi.kallisto,bulk.meta,design=~category+batch+prep_batch)


ddsSE <- estimateSizeFactors(ddsSE)
idx <- rowSums( counts(ddsSE, normalized=TRUE) >= 5 ) >= 3
ddsSE <- ddsSE[idx,]

normalized_counts <- counts(ddsSE, normalized=TRUE)

mod  <- model.matrix(~category+sex+age+race+batch+prep_batch+BSA+thyroid+pacer+WT+HT+BMI, colData(ddsSE))
mod0 <- model.matrix(~ sex+age+race+batch+prep_batch+BSA+thyroid+pacer+WT+HT+BMI, colData(ddsSE))
svseq <- svaseq(normalized_counts, mod, mod0)


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



ddsSE <- DESeq(ddsSE)

vstSE <- vst(ddsSE,blind = FALSE)

mat <- assay(vstSE)
mm <- model.matrix(~category, colData(vstSE))
mat <- limma::removeBatchEffect(mat, covariates=colData(vstSE)[,27:47], design=mm)
assay(vstSE) <- mat

nf.vs.prv <- lfcShrink(ddsSE,contrast=c('category','NF','pRV'), type="ashr")
nf.vs.rvf <- lfcShrink(ddsSE,contrast=c('category','NF','RVF'), type="ashr")
prv.vs.rvf <- lfcShrink(ddsSE,contrast=c('category','pRV','RVF'), type="ashr")

saveRDS(vstSE,'~/Documents/bulkRNAseq_vst.rds')

#######################################
#############  FIGURE 1B  #############
#######################################

pdf('~/Downloads/bulk_RNAseq_pca.pdf',width=10,height=3)
plotPCA(vstSE,intgroup=c("category"),ntop=19355) + theme_classic() + theme(axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.title.x=element_text(size=24),axis.title.y=element_text(size=24),legend.title=element_text(size=24),legend.text=element_text(size=24)) + labs(color='Disease',shape='Sex') + geom_point(size=3.5)
dev.off()

pdf('~/Downloads/bulk_RNAseq_pca_pRV_RVF.pdf',width=10,height=3)
plotPCA(vstSE[,category %in% c('pRV','RVF')],intgroup=c("category"),ntop=19355) + theme_classic() + theme(axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.title.x=element_text(size=24),axis.title.y=element_text(size=24),legend.title=element_text(size=24),legend.text=element_text(size=24)) + labs(color='Disease',shape='Sex') + geom_point(size=3.5)
dev.off()

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = viridis::mako(10))
pdf('~/Downloads/bulk_pca_age.pdf',width=10,height=6)

plotPCA(vstSE,intgroup=c("age"),ntop=19355) + theme_classic() + theme(axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.title.x=element_text(size=24),axis.title.y=element_text(size=24),legend.title=element_text(size=24),legend.text=element_text(size=24)) + labs(color='Age',shape='Sex')+ sc
dev.off()

pdf('~/Downloads/bulk_pca_disease.pdf',width=10,height=6)

plotPCA(vstSE,intgroup=c("disease"),ntop=19355) + theme_classic() + theme(axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.title.x=element_text(size=24),axis.title.y=element_text(size=24),legend.title=element_text(size=24),legend.text=element_text(size=24)) + labs(color='Disease',shape='Etiology')
dev.off()

pdf('~/Downloads/bulk_pca_race.pdf',width=10,height=6)

plotPCA(vstSE,intgroup=c("race"),ntop=19355)  + theme_classic() + theme(axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.title.x=element_text(size=24),axis.title.y=element_text(size=24),legend.title=element_text(size=24),legend.text=element_text(size=24)) + labs(color='Disease',shape='Race')
dev.off()


#######################################
##############  FIGURE 1C  ############
#######################################

#Cast to Seurat to use hdWGCNA backend for consistency
library(Seurat)
library(hdWGCNA)
mat <- as(normalized_counts,'dgCMatrix')
colnames(mat)<-subjects
bulk <- CreateAssayObject(counts=normalized_counts, meta.data=bulk.meta, assay = "RNA")

bulk <-  CreateSeuratObject(bulk)
seurat_obj <- bulk
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- AddMetaData(seurat_obj, bulk.meta)
#saveRDS(seurat_obj,'~/Downloads/RV_bulkRNASeq_seurat.rds')


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
    tom_dir ='~/Downloads/hdWGCNA_TOM',
    tom_name='bulkRV', 
    overwrite_tom=TRUE,
    mergeCutHeight=0.15
)

# compute the MEs and kMEs
seurat_obj <- ModuleEigengenes(seurat_obj)
seurat_obj <- ModuleConnectivity(seurat_obj)


# get MEs from seurat object
MEs <- GetMEs(seurat_obj)
mods <- colnames(MEs); mods <- mods[mods != 'grey']

# add MEs to Seurat meta-data for plotting:
meta <- seurat_obj@meta.data
seurat_obj@meta.data <- cbind(meta, MEs)
saveRDS(seurat_obj,'~/Downloads/RV_bulkRNASeq_seurat.rds')


# plot with Seurat's DotPlot function
p <- DotPlot(seurat_obj, features=mods, group.by = 'category')


#Write module assignments to file
data_dir <- "~/Downloads/hdWGCNA_TOM/"
fig_dir <- '~/Downloads/hdWGCNA/'

modules <- GetModules(seurat_obj) %>% subset(module != 'grey')

write.csv(modules, row.names=FALSE, quote=FALSE, file=paste0(data_dir, 'bulk_heart_modules.csv'))

#Enrichments
library(enrichR)

dbs <-c('GO_Biological_Process_2023','GO_Cellular_Component_2023','GO_Molecular_Function_2023','LINCS_L1000_Chem_Pert_up','LINCS_L1000_Chem_Pert_down', 'WikiPathway_2021_Human', 'KEGG_2021_Human')

# compute GO terms:
enrich_list <- list()
seurat_obj <- RunEnrichr(seurat_obj, dbs=dbs)
enrichr_df <- GetEnrichrTable(seurat_obj) %>% subset(P.value < 0.05)

write.table(enrichr_df, quote=FALSE, sep='\t', row.names=FALSE, file=paste0(data_dir, 'bulk_heart_enrichr.tsv'))


#Network visualizations and module eigengene plots
library(reshape2)
library(igraph)


# plot the dendrograms
pdf(paste0(fig_dir, "dendro.pdf"),height=3, width=6)
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
seurat_obj <- RunModuleUMAP(
  seurat_obj,
  n_hubs = 5,
  n_neighbors=10,
  min_dist=0.3,
  spread=2,
  target_weight=0.1,
  supervised=TRUE
)



# get the hub gene UMAP table from the seurat object
umap_df <- GetModuleUMAP(seurat_obj)

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
    ggrepel::geom_text_repel(data = centroid_df, label=centroid_df$cluster, color='black', max.overlaps=Inf, size=3)


  pdf(paste0(data_dir, '_hubgene_umap_ggplot.pdf'), width=5, height=5)
  print(p)
  dev.off()

  hub_genes <- GetHubGenes(seurat_obj, 3)

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
    ggrepel::geom_text_repel(data = centroid_df, label=centroid_df$cluster, color='black', max.overlaps=Inf, size=3, fontface='bold') +
    geom_text_repel(label=plot_df$anno, max.overlaps=Inf, color='black', fontface='italic', size=3) +
    umap_theme() + NoLegend() +
    coord_equal() +
    theme(
      plot.margin = margin(0,0,0,0)
    )

  pdf(paste0('~/Downloads/', 'full', '_hubgene_umap_ggplot.pdf'), width=8, height=8)
  print(p)
  dev.off()

  # plot with igraph
  pdf(paste0(data_dir,'_hubgene_umap_igraph.pdf'), width=12, height=12)
  p <- ModuleUMAPPlot(
    seurat_obj,
    edge.alpha=0.5,
    sample_edges=TRUE,
    keep_grey_edges=FALSE,
    edge_prop=0.075, # taking the top 20% strongest edges in each module
    #label_genes = label_genes,
    label_hubs=3, # how many hub genes to plot per module?
    return_graph = TRUE,
  )
  dev.off()

}

#######################################
#########  FIGURE 1C (cont) ###########
#######################################


dbs <- c('GO_Biological_Process_2023','GO_Cellular_Component_2023','GO_Molecular_Function_2023')

# perform enrichment tests
seurat_obj <- RunEnrichr(
  seurat_obj,
  dbs=dbs, # character vector of enrichr databases to test
  max_genes = 100 # number of genes per module to test. use max_genes = Inf to choose all genes!
)

# retrieve the output table
enrich_df <- GetEnrichrTable(seurat_obj)

EnrichrBarPlot(
  seurat_obj,
  outdir = "~/Downloads/hdWGCNA_TOM/bulk_enrichr_plot_all_genes", 
  n_terms = 5,
  plot_size = c(5,4), # width, height of the output .pdfs
  logscale=TRUE # do you want to show the enrichment as a log scale?
)


#######################################
#############  FIGURE 1D  #############
#######################################

modules <- GetModules(seurat_obj)
color_df <- modules %>% subset(module!='grey') %>%
  select(c(module, color)) %>% distinct %>%
  rename(c(group=module, colour=color))
mods <- levels(modules$module)
mods <- mods[mods!='grey']

# helper function to wrap text
wrapText <- function(x, len) {
    sapply(x, function(y) paste(strwrap(y, len), collapse = "\n"), USE.NAMES = FALSE)
}

combined_output <- GetEnrichrTable(seurat_obj)
selected_terms <- subset(combined_output,db=="GO_Biological_Process_2023")

# subset selected terms
selected_terms <- subset(selected_terms, P.value < 0.05)
idx_top_1 <- match(unique(selected_terms$module),selected_terms$module)
idx_top_3 <- sort(c(idx_top_1,idx_top_1+1,idx_top_1+2))

selected_terms<-selected_terms[idx_top_3,]

# remove GO Term ID
library(stringr)
selected_terms$Term <- str_replace(selected_terms$Term, " \\s*\\([^\\)]+\\)", "")

key_terms <- read.csv('~/Downloads/hdWGCNA_TOM/bulkRNA_GOterms_ofinterest.csv')
selected_terms <- subset(selected_terms,Term %in% key_terms[[1]])



mapping <- labels2colors(1:100)


selected_terms$group <- factor(
  as.character(selected_terms$module),
  levels = mods[order(match(mods,mapping))]
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

color_df <- data.frame(group = selected_terms$group,colour = selected_terms$module)

color_df$group<-paste0('M',match(color_df$group,mapping))
lvls <- stringr::str_sort(unique(color_df$group), numeric = TRUE)
color_df$group <- factor(color_df$group, levels = lvls)
# make the colorbar as its own heatmap
color_df$var <- 1

c_vect <- color_df$colour
names(c_vect) <- color_df$group


selected_terms$group_num <- color_df$group



# GO Term dot plot
p <- selected_terms %>%
  ggplot(aes(x = group, y = wrap, color =logp, size=log(Combined.Score))) +
  geom_point() +
  scale_color_stepsn(colors=rev(magma(256))) +
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
    panel.grid = element_line(size=0.25, color='lightgrey')
  )



colorbar <- color_df %>%
  ggplot(aes(x=group, y=var, fill=group)) +
  geom_tile() +
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


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'selected_GO_terms.pdf'), width=13, height=11.7)
p / colorbar #+ plot_layout(heights=c(20,1))
dev.off()



#######################################
#############  FIGURE 1E  #############
#######################################

source('~/Downloads/hdWGCNA_TOM/spatial_functions.R')



group1 <- seurat_obj@meta.data %>% subset(category == 'NF') %>% rownames
group2 <- seurat_obj@meta.data %>% subset(category == 'pRV') %>% rownames
group3 <- seurat_obj@meta.data %>% subset(category == 'RVF') %>% rownames

DMEs_prv_vs_rvf <- FindDMEs(
  seurat_obj,
  barcodes1 = group2,
  barcodes2 = group3,
  test.use='wilcox',
  pseudocount.use=0,
)


DMEs_nf_vs_rvf <- FindDMEs(
  seurat_obj,
  barcodes1 = group1,
  barcodes2 = group3,
  test.use='wilcox',
  pseudocount.use=0,
)


DMEs_nf_vs_prv <- FindDMEs(
  seurat_obj,
  barcodes1 = group1,
  barcodes2 = group2,
  test.use='wilcox',
  pseudocount.use=0,
)

mods_use <- unique(c(rownames(
  subset(DMEs_prv_vs_rvf,p_val_adj<0.05)),
rownames(subset(DMEs_nf_vs_rvf,p_val_adj<0.05)),
rownames(subset(DMEs_nf_vs_rvf,p_val_adj<0.05))))


# seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)
MEs <- GetMEs(seurat_obj)
mods <- colnames(MEs); mods <- mods[mods != 'grey']


MEs_rename <- paste0('M',match(colnames(MEs),mapping))
colnames(MEs)<-MEs_rename

seurat_obj@meta.data <- seurat_obj@meta.data[,1:80]
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)
temp = data.frame(dummy=MEs[,1])
temp[1:142,]='1'
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, temp)

mods_use <-  paste0('M',match(mods_use,mapping))


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


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'bulk_hME_vln_stack.pdf'), width=5, height=10)
p
dev.off()


#######################################
#############  FIGURE 1F  #############
#######################################
library(JASPAR2020)
library(motifmatchr)
library(TFBSTools)
library(EnsDb.Hsapiens.v86)
library(GenomicRanges)
library(cowplot)
library(GeneOverlap)

theme_set(theme_cowplot())
set.seed(12345)

# get the pfm from JASPAR2020 using TFBSTools
pfm_core <- TFBSTools::getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

seurat_obj <- MotifScan(
  seurat_obj,
  species_genome = 'hg38',
  pfm = pfm_core,
  EnsDb = EnsDb.Hsapiens.v86
)
dim(GetMotifMatrix(seurat_obj))


# TF target genes
target_genes <- GetMotifTargets(seurat_obj)

# overlap between modules & TF target genes:
seurat_obj<- OverlapModulesMotifs(seurat_obj)

# look at the overlap data
head(GetMotifOverlap(seurat_obj))


# plot the top TFs overlapping with


n_tfs = 10
plot_size = c(5, 6)
outdir = '~/Downloads/hdWGCNA_TOM/bulk_motifs/MotifOverlaps/'
motif_font = "helvetica_regular"

if (!dir.exists(outdir)) {
         dir.create(outdir)
     }

wgcna_name <- seurat_obj@misc$active_wgcna
modules <- GetModules(seurat_obj)

mods <- levels(modules$module)
     mods <- mods[mods != "grey"]
     if (is.null(module_names)) {
         module_names <- mods
     }
     overlap_df <- GetMotifOverlap(seurat_obj, wgcna_name)
     motif_df <- GetMotifs(seurat_obj)
     pfm <- GetPFMList(seurat_obj)
     overlap_df$motif_ID <- motif_df$motif_ID[match(overlap_df$tf, 
         motif_df$motif_name)]
     overlap_df <- overlap_df %>% subset(module %in% module_names)
     for (cur_mod in module_names) {
         print(cur_mod)
         plot_df <- overlap_df %>% subset(module == cur_mod) %>% 
             top_n(n_tfs, wt = odds_ratio) %>% arrange(desc(odds_ratio))
         p1 <- plot_df %>% ggplot(aes(y = reorder(tf, odds_ratio), 
             fill = odds_ratio, x = odds_ratio)) + geom_bar(stat = "identity", 
             width = 0.7) + NoLegend() + scale_fill_gradient(high = unique(plot_df$color), 
             low = "grey90") + ylab("") + theme(axis.line.y = element_blank(), 
             axis.text.y = element_blank(), plot.margin = margin(t = 0, 
                 r = 0, b = 0, l = 0))
         plot_list <- list()
         for (i in 1:nrow(plot_df)) {
             cur_id <- plot_df[i, "motif_ID"]
             cur_name <- plot_df[i, "tf"]
             plot_list[[cur_id]] <- ggplot() + ggseqlogo::geom_logo(as.matrix(pfm[[cur_id]]), 
                 font = motif_font) + ggseqlogo::theme_logo() + 
                 xlab("") + ylab(cur_name) + theme(axis.text.x = element_blank(), 
                 axis.text.y = element_blank(), axis.title.y = element_text(angle = 0), 
                 plot.margin = margin(t = 0, r = 0, b = 0, l = 0))
         }
         patch1 <- wrap_plots(plot_list, ncol = 1)
         outplot <- (patch1 | p1) + plot_layout(ncol = 2, widths = c(1, 
             2)) + plot_annotation(title = paste0("Motif overlaps with ", 
             cur_mod), theme = theme(plot.title = element_text(hjust = 0.5)))
         pdf(paste0(outdir, "/", cur_mod, "_motif_overlaps.pdf"), 
             width = plot_size[1], height = plot_size[2], useDingbats = FALSE)
         print(outplot)
         dev.off()
     }



TF_of_interest <- c('SMAD3','SMAD4','OVOL1','OVOL2','FOXA2',
  'FOXA3','E2F1','TFCP2','ETV3','SOX12','MITF','IRF4','IRF8',
  'YY2','KLF2','KLF4','KLF5','SP1','SP3','SP4','NR2F1','Arid3a')
#store_motif_targets <- seurat_obj@misc$motifs$motif_targets

subset_motif_targets <- store_motif_targets[names(store_motif_targets) %in% TF_of_interest]
seurat_obj@misc$motifs$motif_targets <- subset_motif_targets


seurat_obj <- MotifTargetScore(
  seurat_obj,
  method='UCell',
  maxRank = 15000
)



df <- GetMotifOverlap(seurat_obj)

for (i in TF_of_interest){
  print(i)
tf_name = i
cur_df <- df %>% subset(tf == tf_name)
if (length(rownames(cur_df)) == 0) {next}
cur_df$module <- paste0('M',match(cur_df$module,mapping))
cur_df <- cur_df %>% subset(module %in% mods_use)

plot_var <- 'odds_ratio'
p <- cur_df %>%
  ggplot(aes(y=reorder(module, odds_ratio), x=odds_ratio)) +
  geom_bar(stat='identity', fill=cur_df$color) +
  geom_vline(xintercept = 1, linetype='dashed', color='gray') +
  geom_text(aes(label=Significance), color='black', size=10, hjust='center') +
  ylab('') +
  xlab("Odds Ratio") +
  ggtitle(i) +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size=30),
    axis.text.x = element_text(size=30),
    axis.text.y = element_text(size=30)
  )
pdf(paste0('~/Downloads/hdWGCNA_TOM/', tf_name,'_TF_odds.pdf'), width=5, height=10)

print(p)

dev.off()
}


dbs <- c("Human_Gene_Atlas","WikiPathway_2023_Human","Reactome_2022","GO_Biological_Process_2023")

enriched <- enrichr(NF_2_pRV_down_2_RVF_up, dbs)
enriched[[1]] <- subset(enriched[[1]],Adjusted.P.value<0.2)
p1<- ggplot(enriched[[1]][order(enriched[[1]]$Combined.Score,decreasing=T),][rev(1:3),], 
  (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), 
  size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + 
  ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + 
  ggtitle('GO Biological Process Up') + 
  scale_y_discrete(labels= fct_inorder(
    wrapText(sapply(
      strsplit(enriched[[1]][order(enriched[[1]]$Combined.Score,decreasing=T),][rev(1:3),]$Term," \\(GO"),
         `[`, 1),35))) + 
  theme(axis.text=element_text(colour="black"))+
  scale_color_stepsn(colors=rev(magma(256)))
p1

enriched[[3]] <- subset(enriched[[3]],Adjusted.P.value<0.2)
#pdf('~/Downloads/hdWGCNA_TOM/Bulk_down_down_enrichr.pdf',width=6,height=2.5)
p2<- ggplot(enriched[[3]][order(enriched[[3]]$Combined.Score,decreasing=T),][rev(1:3),], 
  (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), 
  size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + 
  ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + 
  ggtitle('GO Biological Process Up') + 
  scale_y_discrete(labels= fct_inorder(
    wrapText(sapply(
      strsplit(enriched[[3]][order(enriched[[3]]$Combined.Score,decreasing=T),][rev(1:3),]$Term," \\(GO"),
         `[`, 1),35))) + 
  theme(axis.text=element_text(colour="black"))+
  scale_color_stepsn(colors=rev(magma(256)))
p2

pdf('~/Downloads/hdWGCNA_TOM/Bulk_down_up_enrichr.pdf',width=5,height=2.5)
p2
dev.off()

#########SCRATCH PAPER OLD STUFF




####
ddsSE <- estimateSizeFactors(ddsSE)

alt.mat <- mat
norm.mat <- (alt.mat - rowMeans(alt.mat))/rowSds(alt.mat)
disease <- colData(vstSE)$category
#bulk.seurat <- CreateSeuratObject(counts=mat, meta.data=colData(vstSE), assay = "RNA")

#in.c <- rownames(subset(g5,p_val_adj<0.05))
in.c <- g2
in.bulk<-in.c[in.c %in% rownames(vstSE)]

out <- colMeans(norm.mat[in.bulk,])
out <- data.frame(disease=disease,score=out)

ggplot(out,aes(factor(disease),score,fill=factor(disease))) + geom_violin() +
  geom_point(position = position_jitter(seed = 1, width = 0.2))



NF_2_pRV_up_2_RVF_up <- read.csv(file="~/Downloads/NF_2_pRV_up_2_RVF_up.csv")[,2]
NF_2_pRV_down_2_RVF_down <- read.csv(file="~/Downloads/NF_2_pRV_down_2_RVF_down.csv")[,2]
NF_2_pRV_flat_2_RVF_up <- read.csv(file="~/Downloads/NF_2_pRV_flat_2_RVF_up.csv")[,2]
NF_2_pRV_flat_2_RVF_down <- read.csv(file="~/Downloads/NF_2_pRV_flat_2_RVF_down.csv")[,2]
NF_2_pRV_down_2_RVF_up <- read.csv(file="~/Downloads/NF_2_pRV_down_2_RVF_up.csv")[,2]
NF_2_pRV_up_2_RVF_down <- read.csv(file="~/Downloads/NF_2_pRV_up_2_RVF_down.csv")[,2]



pool = rownames(alt.mat)
nbin = 24
ctrl = 100
k = FALSE
name = "RV"
seed = 1
features = list()
features[[1]] = g1[g1 %in% rownames(vstSE)]
features[[2]] = g2[g2 %in% rownames(vstSE)]
features[[3]] = g3[g3 %in% rownames(vstSE)]
features[[4]] = g4[g4 %in% rownames(vstSE)]
features[[5]] = g5[g5 %in% rownames(vstSE)]
features[[6]] = g6[g6 %in% rownames(vstSE)]
top.nf <- rownames(head(c1,250))
top.rvf <- rownames(head(c2,250))
top.prv <- rownames(head(c3,250))
features[[7]] = top.nf[top.nf %in% rownames(vstSE)]
features[[8]] = top.rvf[top.rvf %in% rownames(vstSE)]
features[[9]] = top.prv[top.prv %in% rownames(vstSE)]

features[[10]] = intersect(a1,d1)[intersect(a1,d1) %in% rownames(vstSE)]
features[[11]] = intersect(a2,d2)[intersect(a2,d2) %in% rownames(vstSE)]
features[[12]] = intersect(a3,d3)[intersect(a3,d3) %in% rownames(vstSE)]
features[[13]] = intersect(a4,d4)[intersect(a4,d4) %in% rownames(vstSE)]
features[[14]] = intersect(a5,d5)[intersect(a5,d5) %in% rownames(vstSE)]
features[[15]] = intersect(a6,d6)[intersect(a6,d6) %in% rownames(vstSE)]
features[[16]] = NF_2_pRV_up_2_RVF_up
features[[17]] = NF_2_pRV_down_2_RVF_down
features[[18]] = NF_2_pRV_up_2_RVF_down
features[[19]] = NF_2_pRV_down_2_RVF_up
features[[20]] = NF_2_pRV_flat_2_RVF_up
features[[21]] =

intersect.all <- intersect(c(g1,g2,g3,g4,g5,g6),c(NF_2_pRV_up_2_RVF_up, NF_2_pRV_down_2_RVF_down, NF_2_pRV_up_2_RVF_down, NF_2_pRV_down_2_RVF_up, NF_2_pRV_flat_2_RVF_up, NF_2_pRV_flat_2_RVF_down))



cluster.length <- length(x = features)
assay.data <- t(as.data.frame(alt.mat))

pool = rownames(assay.data)


data.avg <- colMeans(x = assay.data[pool, ],na.rm=T)
data.avg <- data.avg[order(data.avg)]
data.cut <- ggplot2::cut_number(x = data.avg + rnorm(n = length(data.avg))/1e+30,
                                n = nbin,
                                labels = FALSE,
                                right = FALSE)
names(x = data.cut) <- names(x = data.avg)
ctrl.use <- vector(mode = "list", length = cluster.length)

# For each of the input gene lists:
for (i in 1:cluster.length) {
  features.use <- features[[i]]
  for (j in 1:length(x = features.use)) {
        ctrl.use[[i]] <- c(ctrl.use[[i]],names(x = sample(x = data.cut[which(x = data.cut == data.cut[features.use[j]])],size = ctrl,replace = FALSE)))
    }
}

ctrl.use <- lapply(X = ctrl.use, FUN = unique)
ctrl.scores <- matrix(data = numeric(length = 1L), nrow = length(x = ctrl.use), ncol = nrow(x = assay.data))

for (i in 1:length(ctrl.use)) {
  # Get control gene names as a vector  
  features.use <- setdiff(ctrl.use[[i]],"")
  # For each cell, calculate the mean expression of *all* of the control genes 
  ctrl.scores[i, ] <- rowMeans(x = assay.data[,features.use])
}

# Similar to the above, create an empty matrix
features.scores <- matrix(data = numeric(length = 1L),nrow = cluster.length,ncol = nrow(x = assay.data))

for (i in 1:cluster.length) {
    features.use <- setdiff(features[[i]],"")
    data.use <- assay.data[, features.use, drop = FALSE]
    features.scores[i, ] <- rowMeans(x = data.use)
}

features.scores.use <- features.scores - ctrl.scores
rownames(x = features.scores.use) <- paste0(name, 1:cluster.length)
features.scores.use <- as.data.frame(x = t(x = features.scores.use))

rownames(x = features.scores.use) <- rownames(x = assay.data)
features.scores.use$disease <- disease

names = c("1343", "1392", "1467", "1561", "1567", "1618", "1632", "1681", "1691", "1692", "1697")

sn <- colData(vstSE)$subject %in% names
features.scores.use$sn <- sn

ggplot(features.scores.use,aes(factor(disease),RV5,fill=factor(disease))) + geom_violin() +
  geom_point(data= features.scores.use[sn,],position = position_jitter(seed = 1, width = 0.2))




#############DEGs

sum(prv.vs.rvf$padj < 0.1, na.rm=TRUE)
resOrdered1 <- prv.vs.rvf[order(prv.vs.rvf$pvalue),]
resOrdered1$padj[is.na(resOrdered1$padj)]=1
signif1 <- rownames(subset(resOrdered1, resOrdered1$padj<0.1))
write.csv(as.data.frame(resOrdered1), 
          file="~/Downloads/pRV_vs_RVF.csv")

sum(nf.vs.rvf$padj < 0.1, na.rm=TRUE)
resOrdered2 <- nf.vs.rvf[order(nf.vs.rvf$pvalue),]
resOrdered2$padj[is.na(resOrdered2$padj)]=1
signif2 <- rownames(subset(resOrdered2, resOrdered2$padj<0.1))

write.csv(as.data.frame(resOrdered2), 
          file="~/Downloads/NF_vs_RVF.csv")

sum(nf.vs.prv$padj < 0.1, na.rm=TRUE)
resOrdered3 <- nf.vs.prv[order(nf.vs.prv$pvalue),]
resOrdered3$padj[is.na(resOrdered3$padj)]=1
signif3 <- rownames(subset(resOrdered3, resOrdered3$padj<0.1))

write.csv(as.data.frame(resOrdered3), 
          file="~/Downloads/NF_vs_pRV.csv")



#####Venn Diagrams
library(VennDiagram)
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")


x <- list(signif3, signif2,signif1)
pdf('~/Downloads/bulk_venn.pdf',width=6,height=6)

grid.newpage()
venn_object <- venn.diagram(x=x, category.names=c('pRV vs NF','RVF vs NF','pRV vs RVF'),filename = NULL,lwd = 2,lty = 'blank',fill = myCol,fontfamily = "sans",cat.fontfamily = "sans",cat.default.pos = "outer")
grid.draw(venn_object)
dev.off()


x <- list(nf.vs.prv.up, nf.vs.prv.down, nf.vs.rvf.up, nf.vs.rvf.down)
myCol <- brewer.pal(4, "Pastel2")
pdf('~/Downloads/bulk_venn_up_down_concordance.pdf',width=6,height=6)

grid.newpage()

venn_object <- venn.diagram(x=x, category.names=c('pRV vs NF up','pRV vs NF down','RVF vs NF up','RVF vs NF down'),filename = NULL,lwd = 2,lty = 'blank',fill = myCol,fontfamily = "sans",cat.fontfamily = "sans",cat.default.pos = "outer")
grid.draw(venn_object)
dev.off()

x <- list(nf.vs.prv.up, nf.vs.prv.down, nf.vs.rvf.up, nf.vs.rvf.down, prv.vs.rvf.up)
myCol <- brewer.pal(5, "Pastel2")
pdf('~/Downloads/bulk_venn_prv_rvf_up_concordance.pdf',width=6,height=6)

grid.newpage()

venn_object <- venn.diagram(x=x, category.names=c('pRV vs NF up','pRV vs NF down','RVF vs NF up','RVF vs NF down','pRV vs RVF up'),filename = NULL,lwd = 2,lty = 'blank',fill = myCol,fontfamily = "sans",cat.fontfamily = "sans",cat.default.pos = "outer")
grid.draw(venn_object)
dev.off()

#
x <- list(nf.vs.prv.up, nf.vs.prv.down, nf.vs.rvf.up, nf.vs.rvf.down, prv.vs.rvf.down)
myCol <- brewer.pal(5, "Pastel2")
pdf('~/Downloads/bulk_venn_prv_rvf_down_concordance.pdf',width=6,height=6)

grid.newpage()

venn_object <- venn.diagram(x=x, category.names=c('pRV vs NF up','pRV vs NF down','RVF vs NF up','RVF vs NF down','pRV vs RVF down'),filename = NULL,lwd = 2,lty = 'blank',fill = myCol,fontfamily = "sans",cat.fontfamily = "sans",cat.default.pos = "outer")
grid.draw(venn_object)
dev.off()










#Enrichments
library(DOSE)
library(enrichR)
library(tidyverse)

websiteLive <- getOption("enrichR.live")
if (websiteLive) {
    listEnrichrSites()
    setEnrichrSite("Enrichr") # Human genes   
}



dbs <- c("ChEA_2022","WikiPathway_2023_Human","Reactome_2016","GO_Biological_Process_2023")



####GO BP

#pRV vs RV

enriched <- enrichr(rownames(subset(prv.vs.rvf,padj<0.05 & log2FoldChange > 0.1)), dbs)
enriched[[4]] <- enriched[[4]][enriched[[4]]$Adjusted.P.value < 0.1,]
pdf('~/Downloads/bulkRNAseq_DEG_GO_BP_pRV_vs_RVF_Up.pdf',width=10,height=5.5)
p4<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:20),], (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + ggtitle('GO Biological Process Up') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:20),]$Term," \\(GO"), `[`, 1))) + theme(axis.text=element_text(colour="black"))
p4
dev.off()

enriched <- enrichr(rownames(subset(prv.vs.rvf,padj<0.05 & log2FoldChange < -0.1)), dbs)
enriched[[4]] <- enriched[[4]][enriched[[4]]$Adjusted.P.value < 0.1,]
pdf('~/Downloads/bulkRNAseq_DEG_GO_BP_pRV_vs_RVF_Down.pdf',width=10,height=3)
p4<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:8),], (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + ggtitle('GO Biological Process Down') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:20),]$Term," \\(GO"), `[`, 1))) + theme(axis.text=element_text(colour="black"))
p4
dev.off()

#NF vs pRV

enriched <- enrichr(rownames(subset(nf.vs.prv,padj<0.05 & log2FoldChange > 0.1)), dbs)
enriched[[4]] <- enriched[[4]][enriched[[4]]$Adjusted.P.value < 0.1,]
pdf('~/Downloads/bulkRNAseq_DEG_GO_BP_NF_vs_pRV_Up.pdf',width=10,height=5.5)
p4<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:20),], (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + ggtitle('GO Biological Process Up') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:20),]$Term," \\(GO"), `[`, 1))) + theme(axis.text=element_text(colour="black"))
p4
dev.off()

enriched <- enrichr(rownames(subset(nf.vs.prv,padj<0.05 & log2FoldChange < -0.1)), dbs)
enriched[[4]] <- enriched[[4]][enriched[[4]]$Adjusted.P.value < 0.1,]
pdf('~/Downloads/bulkRNAseq_DEG_GO_BP_NF_vs_pRV_Down.pdf',width=10,height=5.5)
p4<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:20),], (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + ggtitle('GO Biological Process Down') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:20),]$Term," \\(GO"), `[`, 1))) + theme(axis.text=element_text(colour="black"))
p4
dev.off()

#NF vs RVF

enriched <- enrichr(rownames(subset(nf.vs.rvf,padj<0.05 & log2FoldChange > 0.1)), dbs)
enriched[[4]] <- enriched[[4]][enriched[[4]]$Adjusted.P.value < 0.1,]
pdf('~/Downloads/bulkRNAseq_DEG_GO_BP_NF_vs_RVF_Up.pdf',width=10,height=5.5)
p4<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:20),], (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + ggtitle('GO Biological Process Up') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:20),]$Term," \\(GO"), `[`, 1))) + theme(axis.text=element_text(colour="black"))
p4
dev.off()

enriched <- enrichr(rownames(subset(nf.vs.rvf,padj<0.05 & log2FoldChange < -0.1)), dbs)
enriched[[4]] <- enriched[[4]][enriched[[4]]$Adjusted.P.value < 0.1,]
pdf('~/Downloads/bulkRNAseq_DEG_GO_BP_NF_vs_RVF_Down.pdf',width=10,height=5.5)
p4<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:20),], (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + ggtitle('GO Biological Process Down') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:20),]$Term," \\(GO"), `[`, 1))) + theme(axis.text=element_text(colour="black"))
p4
dev.off()


###CHeA

#pRV vs RV

enriched <- enrichr(rownames(subset(prv.vs.rvf,padj<0.05 & log2FoldChange > 0.1)), dbs)
enriched[[1]] <- enriched[[1]][enriched[[1]]$Adjusted.P.value < 0.1,]
pdf('~/Downloads/bulkRNAseq_DEG_TF_pRV_vs_RVF_Up.pdf',width=5,height=5.5)
p4<- ggplot(enriched[[1]][order(enriched[[1]]$Combined.Score,decreasing=T),][rev(1:20),], (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + ggtitle('ChEA Up') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[1]][order(enriched[[1]]$Combined.Score,decreasing=T),][rev(1:20),]$Term," "), `[`, 1))) + theme(axis.text=element_text(colour="black"))
p4
dev.off()

enriched <- enrichr(rownames(subset(prv.vs.rvf,padj<0.05 & log2FoldChange < -0.1)), dbs)
enriched[[1]] <- enriched[[1]][enriched[[1]]$Adjusted.P.value < 0.1,]
pdf('~/Downloads/bulkRNAseq_DEG_TF_pRV_vs_RVF_Down.pdf',width=5,height=5.5)
p4<- ggplot(enriched[[1]][order(enriched[[1]]$Combined.Score,decreasing=T),][rev(1:20),], (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + ggtitle('ChEA Down') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[1]][order(enriched[[1]]$Combined.Score,decreasing=T),][rev(1:20),]$Term," "), `[`, 1))) + theme(axis.text=element_text(colour="black"))
p4
dev.off()

#NF vs pRV

enriched <- enrichr(rownames(subset(nf.vs.prv,padj<0.05 & log2FoldChange > 0.1)), dbs)
enriched[[1]] <- enriched[[1]][enriched[[1]]$Adjusted.P.value < 0.1,]
pdf('~/Downloads/bulkRNAseq_DEG_TF_NF_vs_pRV_Up.pdf',width=5,height=5.5)
p4<- ggplot(enriched[[1]][order(enriched[[1]]$Combined.Score,decreasing=T),][rev(1:20),], (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + ggtitle('ChEA Up') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[1]][order(enriched[[1]]$Combined.Score,decreasing=T),][rev(1:20),]$Term," "), `[`, 1))) + theme(axis.text=element_text(colour="black"))
p4
dev.off()

enriched <- enrichr(rownames(subset(nf.vs.prv,padj<0.05 & log2FoldChange < -0.1)), dbs)
enriched[[1]] <- enriched[[1]][enriched[[1]]$Adjusted.P.value < 0.1,]
pdf('~/Downloads/bulkRNAseq_DEG_TF_NF_vs_pRV_Down.pdf',width=5,height=5.5)
p4<- ggplot(enriched[[1]][order(enriched[[1]]$Combined.Score,decreasing=T),][rev(1:20),], (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + ggtitle('ChEA Down') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[1]][order(enriched[[1]]$Combined.Score,decreasing=T),][rev(1:20),]$Term," "), `[`, 1))) + theme(axis.text=element_text(colour="black"))
p4
dev.off()

#NF vs RVF

enriched <- enrichr(rownames(subset(nf.vs.rvf,padj<0.05 & log2FoldChange > 0.1)), dbs)
enriched[[1]] <- enriched[[1]][enriched[[1]]$Adjusted.P.value < 0.1,]
pdf('~/Downloads/bulkRNAseq_DEG_TF_NF_vs_RVF_Up.pdf',width=5,height=5.5)
p4<- ggplot(enriched[[1]][order(enriched[[1]]$Combined.Score,decreasing=T),][rev(1:20),], (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + ggtitle('ChEA Up') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[1]][order(enriched[[1]]$Combined.Score,decreasing=T),][rev(1:20),]$Term," "), `[`, 1))) + theme(axis.text=element_text(colour="black"))
p4
dev.off()

enriched <- enrichr(rownames(subset(nf.vs.rvf,padj<0.05 & log2FoldChange < -0.1)), dbs)
enriched[[1]] <- enriched[[1]][enriched[[1]]$Adjusted.P.value < 0.1,]
pdf('~/Downloads/bulkRNAseq_DEG_TF_NF_vs_RVF_Down.pdf',width=5,height=5.5)
p4<- ggplot(enriched[[1]][order(enriched[[1]]$Combined.Score,decreasing=T),][rev(1:20),], (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + ggtitle('ChEA Down') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[1]][order(enriched[[1]]$Combined.Score,decreasing=T),][rev(1:20),]$Term," "), `[`, 1))) + theme(axis.text=element_text(colour="black"))
p4
dev.off()





####GENE ENRICHMENT FOR DIRECTIONAL ANALYSIS
####
#NF_2_pRV_up_2_RVF_up
#NF_2_pRV_down_2_RVF_down 
#NF_2_pRV_up_2_RVF_down 
#NF_2_pRV_down_2_RVF_up
#NF_2_pRV_flat_2_RVF_up
#NF_2_pRV_flat_2_RVF_down

enriched <- enrichr(NF_2_pRV_up_2_RVF_up, dbs)
enriched[[4]] <- enriched[[4]][enriched[[4]]$Adjusted.P.value < 0.1,]
pdf('~/Downloads/bulkRNAseq_DEG_GO_BP_NF_2_pRV_up_2_RVF_up.pdf',width=10,height=5.5)
p4<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:20),], (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + ggtitle('GO Biological Process') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:20),]$Term," \\(GO"), `[`, 1))) + theme(axis.text=element_text(colour="black"))
p4
dev.off()

enriched <- enrichr(NF_2_pRV_down_2_RVF_down, dbs)
enriched[[4]] <- enriched[[4]][enriched[[4]]$Adjusted.P.value < 1,]
pdf('~/Downloads/bulkRNAseq_DEG_GO_BP_NF_2_pRV_down_2_RVF_down.pdf',width=10,height=5.5)
p4<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:20),], (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + ggtitle('GO Biological Process') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:20),]$Term," \\(GO"), `[`, 1))) + theme(axis.text=element_text(colour="black"))
p4
dev.off()

enriched <- enrichr(NF_2_pRV_up_2_RVF_down, dbs)
enriched[[4]] <- enriched[[4]][enriched[[4]]$Adjusted.P.value < 1,]
pdf('~/Downloads/bulkRNAseq_DEG_GO_BP_NF_2_pRV_up_2_RVF_down.pdf',width=10,height=5.5)
p4<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:20),], (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + ggtitle('GO Biological Process') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:20),]$Term," \\(GO"), `[`, 1))) + theme(axis.text=element_text(colour="black"))
p4
dev.off()

enriched <- enrichr(NF_2_pRV_down_2_RVF_up, dbs)
enriched[[4]] <- enriched[[4]][enriched[[4]]$Adjusted.P.value < 1,]
pdf('~/Downloads/bulkRNAseq_DEG_GO_BP_NF_2_pRV_down_2_RVF_up.pdf',width=10,height=5.5)
p4<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:20),], (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + ggtitle('GO Biological Process') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:20),]$Term," \\(GO"), `[`, 1))) + theme(axis.text=element_text(colour="black"))
p4
dev.off()


enriched <- enrichr(NF_2_pRV_flat_2_RVF_up, dbs)
enriched[[4]] <- enriched[[4]][enriched[[4]]$Adjusted.P.value < 1,]
pdf('~/Downloads/bulkRNAseq_DEG_GO_BP_NF_2_pRV_flat_2_RVF_up.pdf',width=10,height=5.5)
p4<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:20),], (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + ggtitle('GO Biological Process') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:20),]$Term," \\(GO"), `[`, 1))) + theme(axis.text=element_text(colour="black"))
p4
dev.off()

enriched <- enrichr(NF_2_pRV_flat_2_RVF_down, dbs)
enriched[[4]] <- enriched[[4]][enriched[[4]]$Adjusted.P.value < 1,]
pdf('~/Downloads/bulkRNAseq_DEG_GO_BP_NF_2_pRV_flat_2_RVF_down.pdf',width=10,height=5.5)
p4<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:20),], (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + ggtitle('GO Biological Process') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:20),]$Term," \\(GO"), `[`, 1))) + theme(axis.text=element_text(colour="black"))
p4
dev.off()


##TF
enriched <- enrichr(NF_2_pRV_up_2_RVF_up, dbs)
enriched[[1]] <- enriched[[1]][enriched[[1]]$Adjusted.P.value < 0.1,]
pdf('~/Downloads/bulkRNAseq_DEG_TF_NF_2_pRV_up_2_RVF_up.pdf',width=10,height=5.5)
p4<- ggplot(enriched[[1]][order(enriched[[1]]$Combined.Score,decreasing=T),][rev(1:20),], (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + ggtitle('ChEA') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[1]][order(enriched[[1]]$Combined.Score,decreasing=T),][rev(1:20),]$Term," "), `[`, 1))) + theme(axis.text=element_text(colour="black"))
p4
dev.off()

enriched <- enrichr(NF_2_pRV_down_2_RVF_down, dbs)
enriched[[1]] <- enriched[[1]][enriched[[1]]$Adjusted.P.value < 0.1,]
pdf('~/Downloads/bulkRNAseq_DEG_TF_NF_2_pRV_down_2_RVF_down.pdf',width=10,height=5.5)
p4<- ggplot(enriched[[1]][order(enriched[[1]]$Combined.Score,decreasing=T),][rev(1:20),], (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + ggtitle('ChEA') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[1]][order(enriched[[1]]$Combined.Score,decreasing=T),][rev(1:20),]$Term," "), `[`, 1))) + theme(axis.text=element_text(colour="black"))
p4
dev.off()

enriched <- enrichr(NF_2_pRV_up_2_RVF_down, dbs)
enriched[[4]] <- enriched[[1]][enriched[[1]]$Adjusted.P.value < 0.1,]
pdf('~/Downloads/bulkRNAseq_DEG_TF_NF_2_pRV_up_2_RVF_down.pdf',width=10,height=5.5)
p4<- ggplot(enriched[[1]][order(enriched[[1]]$Combined.Score,decreasing=T),][rev(1:20),], (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + ggtitle('ChEA') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[1]][order(enriched[[1]]$Combined.Score,decreasing=T),][rev(1:20),]$Term," "), `[`, 1))) + theme(axis.text=element_text(colour="black"))
p4
dev.off()

enriched <- enrichr(NF_2_pRV_down_2_RVF_up, dbs)
enriched[[4]] <- enriched[[1]][enriched[[1]]$Adjusted.P.value < 0.1,]
pdf('~/Downloads/bulkRNAseq_DEG_TF_NF_2_pRV_down_2_RVF_up.pdf',width=10,height=5.5)
p4<- ggplot(enriched[[1]][order(enriched[[1]]$Combined.Score,decreasing=T),][rev(1:20),], (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + ggtitle('ChEA') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[1]][order(enriched[[1]]$Combined.Score,decreasing=T),][rev(1:20),]$Term," "), `[`, 1))) + theme(axis.text=element_text(colour="black"))
p4
dev.off()


enriched <- enrichr(NF_2_pRV_flat_2_RVF_up, dbs)
enriched[[4]] <- enriched[[1]][enriched[[1]]$Adjusted.P.value < 0.1,]
pdf('~/Downloads/bulkRNAseq_DEG_TF_NF_2_pRV_flat_2_RVF_up.pdf',width=10,height=5.5)
p4<- ggplot(enriched[[1]][order(enriched[[1]]$Combined.Score,decreasing=T),][rev(1:20),], (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + ggtitle('ChEA') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[1]][order(enriched[[1]]$Combined.Score,decreasing=T),][rev(1:20),]$Term," "), `[`, 1))) + theme(axis.text=element_text(colour="black"))
p4
dev.off()

enriched <- enrichr(NF_2_pRV_flat_2_RVF_down, dbs)
enriched[[4]] <- enriched[[1]][enriched[[1]]$Adjusted.P.value < 0.1,]
pdf('~/Downloads/bulkRNAseq_DEG_TF_NF_2_pRV_flat_2_RVF_down.pdf',width=10,height=5.5)
p4<- ggplot(enriched[[1]][order(enriched[[1]]$Combined.Score,decreasing=T),][rev(1:20),], (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + ggtitle('ChEA') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[1]][order(enriched[[1]]$Combined.Score,decreasing=T),][rev(1:20),]$Term," "), `[`, 1))) + theme(axis.text=element_text(colour="black"))
p4
dev.off()


#Save gene lists
#NF_2_pRV_up_2_RVF_up
#NF_2_pRV_down_2_RVF_down 
#NF_2_pRV_up_2_RVF_down 
#NF_2_pRV_down_2_RVF_up
#NF_2_pRV_flat_2_RVF_up
#NF_2_pRV_flat_2_RVF_down
write.csv(as.data.frame(NF_2_pRV_up_2_RVF_up),file="~/Downloads/NF_2_pRV_up_2_RVF_up.csv")
write.csv(as.data.frame(NF_2_pRV_down_2_RVF_down),file="~/Downloads/NF_2_pRV_down_2_RVF_down.csv")
write.csv(as.data.frame(NF_2_pRV_up_2_RVF_down),file="~/Downloads/NF_2_pRV_up_2_RVF_down.csv")
write.csv(as.data.frame(NF_2_pRV_down_2_RVF_up),file="~/Downloads/NF_2_pRV_down_2_RVF_up.csv")
write.csv(as.data.frame(NF_2_pRV_flat_2_RVF_up),file="~/Downloads/NF_2_pRV_flat_2_RVF_up.csv")
write.csv(as.data.frame(NF_2_pRV_flat_2_RVF_down),file="~/Downloads/NF_2_pRV_flat_2_RVF_down.csv")







#Compare gene lists to sc data
#Pooled data first

M1<-readRDS('/Volumes/Extreme SSD/Final_Analysis/CellTypes/Post_R3_FINAL_with_counts.rds')
M1$Names_group <- paste0(M1$Names,'_',M1$group)
M2<-AggregateExpression(M1, group.by = c("Names","group","patient"), assays = "RNA", return.seurat = TRUE)
M2$Names_group <- paste0(M2$Names,'_',M2$group)
Idents(M2) <- "group"
Idents(M1) <- "group"
M1<-PrepSCTFindMarkers(M1)


a1 <- FindMarkers(M2,ident.1=c('pRV'),ident.2=c('RVF'),logfc.threshold=0,test.use="DESeq2")
a2 <- FindMarkers(M2,ident.1=c('NF'),ident.2=c('pRV'),logfc.threshold=0,test.use="DESeq2")
a3 <- FindMarkers(M2,ident.1=c('NF'),ident.2=c('RVF'),logfc.threshold=0,test.use="DESeq2")

b1 <- FindMarkers(M1,ident.1=c('pRV'),ident.2=c('RVF'),logfc.threshold=0)
b2 <- FindMarkers(M1,ident.1=c('NF'),ident.2=c('pRV'),logfc.threshold=0)
b3 <- FindMarkers(M1,ident.1=c('NF'),ident.2=c('RVF'),logfc.threshold=0)


sc.prv.vs.rvf.up <- rownames(subset(b1, p_val_adj<0.05 & avg_log2FC>0.1))
sc.prv.vs.rvf.down <- rownames(subset(b1, p_val_adj<0.05 & avg_log2FC < -0.1))

sc.nf.vs.prv.up <- rownames(subset(b2, p_val_adj<0.05 & avg_log2FC>0.1))
sc.nf.vs.prv.down <- rownames(subset(b2, p_val_adj<0.05 & avg_log2FC < -0.1))

sc.nf.vs.rvf.up <- rownames(subset(b3, p_val_adj<0.05 & avg_log2FC>0.1))
sc.nf.vs.rvf.down <- rownames(subset(b3, p_val_adj<0.05 & avg_log2FC < -0.1))

#Gene expression gradients
sc.NF_2_pRV_up_2_RVF_up <- intersect(intersect(sc.nf.vs.prv.down, sc.nf.vs.rvf.down), sc.prv.vs.rvf.down)
sc.NF_2_pRV_down_2_RVF_down <- intersect(intersect(sc.nf.vs.prv.up, sc.nf.vs.rvf.up), sc.prv.vs.rvf.up)
sc.NF_2_pRV_up_2_RVF_down <- setdiff(intersect(sc.nf.vs.prv.down, sc.prv.vs.rvf.up), sc.nf.vs.rvf.down)
sc.NF_2_pRV_down_2_RVF_up <- setdiff(intersect(sc.nf.vs.prv.up, sc.prv.vs.rvf.down), sc.nf.vs.rvf.up)
sc.NF_2_pRV_flat_2_RVF_up <- setdiff(intersect(sc.nf.vs.rvf.up, sc.prv.vs.rvf.up), sc.nf.vs.prv.up)
sc.NF_2_pRV_flat_2_RVF_down <- setdiff(intersect(sc.nf.vs.rvf.down, sc.prv.vs.rvf.down), sc.nf.vs.prv.down)

write.csv(as.data.frame(sc.NF_2_pRV_up_2_RVF_up),file="~/Downloads/sc_NF_2_pRV_up_2_RVF_up.csv")
write.csv(as.data.frame(sc.NF_2_pRV_down_2_RVF_down),file="~/Downloads/sc_NF_2_pRV_down_2_RVF_down.csv")
write.csv(as.data.frame(sc.NF_2_pRV_up_2_RVF_down),file="~/Downloads/sc_NF_2_pRV_up_2_RVF_down.csv")
write.csv(as.data.frame(sc.NF_2_pRV_down_2_RVF_up),file="~/Downloads/sc_NF_2_pRV_down_2_RVF_up.csv")
write.csv(as.data.frame(sc.NF_2_pRV_flat_2_RVF_up),file="~/Downloads/sc_NF_2_pRV_flat_2_RVF_up.csv")
write.csv(as.data.frame(sc.NF_2_pRV_flat_2_RVF_down),file="~/Downloads/sc_NF_2_pRV_flat_2_RVF_down.csv")


x <- list(prv.vs.rvf.up, prv.vs.rvf.down, sc.prv.vs.rvf.up, sc.prv.vs.rvf.down)
myCol <- brewer.pal(4, "Pastel2")
pdf('~/Downloads/sc_vs_bulk_venn_up_down_concordance.pdf',width=6,height=6)

grid.newpage()

venn_object <- venn.diagram(x=x, category.names=c('pRV vs RVF up','pRV vs RVF down','sc pRV vs RVF up','sc pRV vs RVF down'),filename = NULL,lwd = 2,lty = 'blank',fill = myCol,fontfamily = "sans",cat.fontfamily = "sans",cat.default.pos = "outer")
grid.draw(venn_object)
dev.off()

#CMs only
Idents(M1) <- "Names_group"

b1 <- FindMarkers(M1,ident.1=c('CM_pRV'),ident.2=c('CM_RVF'),logfc.threshold=0)
b2 <- FindMarkers(M1,ident.1=c('CM_NF'),ident.2=c('CM_pRV'),logfc.threshold=0)
b3 <- FindMarkers(M1,ident.1=c('CM_NF'),ident.2=c('CM_RVF'),logfc.threshold=0)


sc.prv.vs.rvf.up <- rownames(subset(b1, p_val_adj<0.05 & avg_log2FC>0.1))
sc.prv.vs.rvf.down <- rownames(subset(b1, p_val_adj<0.05 & avg_log2FC < -0.1))

sc.nf.vs.prv.up <- rownames(subset(b2, p_val_adj<0.05 & avg_log2FC>0.1))
sc.nf.vs.prv.down <- rownames(subset(b2, p_val_adj<0.05 & avg_log2FC < -0.1))

sc.nf.vs.rvf.up <- rownames(subset(b3, p_val_adj<0.05 & avg_log2FC>0.1))
sc.nf.vs.rvf.down <- rownames(subset(b3, p_val_adj<0.05 & avg_log2FC < -0.1))

#Gene expression gradients
sc.NF_2_pRV_up_2_RVF_up <- intersect(intersect(sc.nf.vs.prv.down, sc.nf.vs.rvf.down), sc.prv.vs.rvf.down)
sc.NF_2_pRV_down_2_RVF_down <- intersect(intersect(sc.nf.vs.prv.up, sc.nf.vs.rvf.up), sc.prv.vs.rvf.up)
sc.NF_2_pRV_up_2_RVF_down <- setdiff(intersect(sc.nf.vs.prv.down, sc.prv.vs.rvf.up), sc.nf.vs.rvf.down)
sc.NF_2_pRV_down_2_RVF_up <- setdiff(intersect(sc.nf.vs.prv.up, sc.prv.vs.rvf.down), sc.nf.vs.rvf.up)
sc.NF_2_pRV_flat_2_RVF_up <- setdiff(intersect(sc.nf.vs.rvf.up, sc.prv.vs.rvf.up), sc.nf.vs.prv.up)
sc.NF_2_pRV_flat_2_RVF_down <- setdiff(intersect(sc.nf.vs.rvf.down, sc.prv.vs.rvf.down), sc.nf.vs.prv.down)

write.csv(as.data.frame(sc.NF_2_pRV_up_2_RVF_up),file="~/Downloads/sc_cm_NF_2_pRV_up_2_RVF_up.csv")
write.csv(as.data.frame(sc.NF_2_pRV_down_2_RVF_down),file="~/Downloads/sc_cm_NF_2_pRV_down_2_RVF_down.csv")
write.csv(as.data.frame(sc.NF_2_pRV_up_2_RVF_down),file="~/Downloads/sc_cm_NF_2_pRV_up_2_RVF_down.csv")
write.csv(as.data.frame(sc.NF_2_pRV_down_2_RVF_up),file="~/Downloads/sc_cm_NF_2_pRV_down_2_RVF_up.csv")
write.csv(as.data.frame(sc.NF_2_pRV_flat_2_RVF_up),file="~/Downloads/sc_cm_NF_2_pRV_flat_2_RVF_up.csv")
write.csv(as.data.frame(sc.NF_2_pRV_flat_2_RVF_down),file="~/Downloads/sc_cm_NF_2_pRV_flat_2_RVF_down.csv")



#FBs only
Idents(M1) <- "Names_group"

b1 <- FindMarkers(M1,ident.1=c('FB_pRV'),ident.2=c('FB_RVF'),logfc.threshold=0)
b2 <- FindMarkers(M1,ident.1=c('FB_NF'),ident.2=c('FB_pRV'),logfc.threshold=0)
b3 <- FindMarkers(M1,ident.1=c('FB_NF'),ident.2=c('FB_RVF'),logfc.threshold=0)


sc.prv.vs.rvf.up <- rownames(subset(b1, p_val_adj<0.05 & avg_log2FC>0.1))
sc.prv.vs.rvf.down <- rownames(subset(b1, p_val_adj<0.05 & avg_log2FC < -0.1))

sc.nf.vs.prv.up <- rownames(subset(b2, p_val_adj<0.05 & avg_log2FC>0.1))
sc.nf.vs.prv.down <- rownames(subset(b2, p_val_adj<0.05 & avg_log2FC < -0.1))

sc.nf.vs.rvf.up <- rownames(subset(b3, p_val_adj<0.05 & avg_log2FC>0.1))
sc.nf.vs.rvf.down <- rownames(subset(b3, p_val_adj<0.05 & avg_log2FC < -0.1))

#Gene expression gradients
sc.NF_2_pRV_up_2_RVF_up <- intersect(intersect(sc.nf.vs.prv.down, sc.nf.vs.rvf.down), sc.prv.vs.rvf.down)
sc.NF_2_pRV_down_2_RVF_down <- intersect(intersect(sc.nf.vs.prv.up, sc.nf.vs.rvf.up), sc.prv.vs.rvf.up)
sc.NF_2_pRV_up_2_RVF_down <- setdiff(intersect(sc.nf.vs.prv.down, sc.prv.vs.rvf.up), sc.nf.vs.rvf.down)
sc.NF_2_pRV_down_2_RVF_up <- setdiff(intersect(sc.nf.vs.prv.up, sc.prv.vs.rvf.down), sc.nf.vs.rvf.up)
sc.NF_2_pRV_flat_2_RVF_up <- setdiff(intersect(sc.nf.vs.rvf.up, sc.prv.vs.rvf.up), sc.nf.vs.prv.up)
sc.NF_2_pRV_flat_2_RVF_down <- setdiff(intersect(sc.nf.vs.rvf.down, sc.prv.vs.rvf.down), sc.nf.vs.prv.down)

write.csv(as.data.frame(sc.NF_2_pRV_up_2_RVF_up),file="~/Downloads/sc_fb_NF_2_pRV_up_2_RVF_up.csv")
write.csv(as.data.frame(sc.NF_2_pRV_down_2_RVF_down),file="~/Downloads/sc_fb_NF_2_pRV_down_2_RVF_down.csv")
write.csv(as.data.frame(sc.NF_2_pRV_up_2_RVF_down),file="~/Downloads/sc_fb_NF_2_pRV_up_2_RVF_down.csv")
write.csv(as.data.frame(sc.NF_2_pRV_down_2_RVF_up),file="~/Downloads/sc_fb_NF_2_pRV_down_2_RVF_up.csv")
write.csv(as.data.frame(sc.NF_2_pRV_flat_2_RVF_up),file="~/Downloads/sc_fb_NF_2_pRV_flat_2_RVF_up.csv")
write.csv(as.data.frame(sc.NF_2_pRV_flat_2_RVF_down),file="~/Downloads/sc_fb_NF_2_pRV_flat_2_RVF_down.csv")

#EC only
Idents(M1) <- "Names_group"

b1 <- FindMarkers(M1,ident.1=c('EC_pRV'),ident.2=c('EC_RVF'),logfc.threshold=0)
b2 <- FindMarkers(M1,ident.1=c('EC_NF'),ident.2=c('EC_pRV'),logfc.threshold=0)
b3 <- FindMarkers(M1,ident.1=c('EC_NF'),ident.2=c('EC_RVF'),logfc.threshold=0)


sc.prv.vs.rvf.up <- rownames(subset(b1, p_val_adj<0.05 & avg_log2FC>0.1))
sc.prv.vs.rvf.down <- rownames(subset(b1, p_val_adj<0.05 & avg_log2FC < -0.1))

sc.nf.vs.prv.up <- rownames(subset(b2, p_val_adj<0.05 & avg_log2FC>0.1))
sc.nf.vs.prv.down <- rownames(subset(b2, p_val_adj<0.05 & avg_log2FC < -0.1))

sc.nf.vs.rvf.up <- rownames(subset(b3, p_val_adj<0.05 & avg_log2FC>0.1))
sc.nf.vs.rvf.down <- rownames(subset(b3, p_val_adj<0.05 & avg_log2FC < -0.1))

#Gene expression gradients
sc.NF_2_pRV_up_2_RVF_up <- intersect(intersect(sc.nf.vs.prv.down, sc.nf.vs.rvf.down), sc.prv.vs.rvf.down)
sc.NF_2_pRV_down_2_RVF_down <- intersect(intersect(sc.nf.vs.prv.up, sc.nf.vs.rvf.up), sc.prv.vs.rvf.up)
sc.NF_2_pRV_up_2_RVF_down <- setdiff(intersect(sc.nf.vs.prv.down, sc.prv.vs.rvf.up), sc.nf.vs.rvf.down)
sc.NF_2_pRV_down_2_RVF_up <- setdiff(intersect(sc.nf.vs.prv.up, sc.prv.vs.rvf.down), sc.nf.vs.rvf.up)
sc.NF_2_pRV_flat_2_RVF_up <- setdiff(intersect(sc.nf.vs.rvf.up, sc.prv.vs.rvf.up), sc.nf.vs.prv.up)
sc.NF_2_pRV_flat_2_RVF_down <- setdiff(intersect(sc.nf.vs.rvf.down, sc.prv.vs.rvf.down), sc.nf.vs.prv.down)

write.csv(as.data.frame(sc.NF_2_pRV_up_2_RVF_up),file="~/Downloads/sc_ec_NF_2_pRV_up_2_RVF_up.csv")
write.csv(as.data.frame(sc.NF_2_pRV_down_2_RVF_down),file="~/Downloads/sc_ec_NF_2_pRV_down_2_RVF_down.csv")
write.csv(as.data.frame(sc.NF_2_pRV_up_2_RVF_down),file="~/Downloads/sc_ec_NF_2_pRV_up_2_RVF_down.csv")
write.csv(as.data.frame(sc.NF_2_pRV_down_2_RVF_up),file="~/Downloads/sc_ec_NF_2_pRV_down_2_RVF_up.csv")
write.csv(as.data.frame(sc.NF_2_pRV_flat_2_RVF_up),file="~/Downloads/sc_ec_NF_2_pRV_flat_2_RVF_up.csv")
write.csv(as.data.frame(sc.NF_2_pRV_flat_2_RVF_down),file="~/Downloads/sc_ec_NF_2_pRV_flat_2_RVF_down.csv")

#Myeloid only
Idents(M1) <- "Names_group"

b1 <- FindMarkers(M1,ident.1=c('Myeloid_pRV'),ident.2=c('Myeloid_RVF'),logfc.threshold=0)
b2 <- FindMarkers(M1,ident.1=c('Myeloid_NF'),ident.2=c('Myeloid_pRV'),logfc.threshold=0)
b3 <- FindMarkers(M1,ident.1=c('Myeloid_NF'),ident.2=c('Myeloid_RVF'),logfc.threshold=0)


sc.prv.vs.rvf.up <- rownames(subset(b1, p_val_adj<0.05 & avg_log2FC>0.1))
sc.prv.vs.rvf.down <- rownames(subset(b1, p_val_adj<0.05 & avg_log2FC < -0.1))

sc.nf.vs.prv.up <- rownames(subset(b2, p_val_adj<0.05 & avg_log2FC>0.1))
sc.nf.vs.prv.down <- rownames(subset(b2, p_val_adj<0.05 & avg_log2FC < -0.1))

sc.nf.vs.rvf.up <- rownames(subset(b3, p_val_adj<0.05 & avg_log2FC>0.1))
sc.nf.vs.rvf.down <- rownames(subset(b3, p_val_adj<0.05 & avg_log2FC < -0.1))

#Gene expression gradients
sc.NF_2_pRV_up_2_RVF_up <- intersect(intersect(sc.nf.vs.prv.down, sc.nf.vs.rvf.down), sc.prv.vs.rvf.down)
sc.NF_2_pRV_down_2_RVF_down <- intersect(intersect(sc.nf.vs.prv.up, sc.nf.vs.rvf.up), sc.prv.vs.rvf.up)
sc.NF_2_pRV_up_2_RVF_down <- setdiff(intersect(sc.nf.vs.prv.down, sc.prv.vs.rvf.up), sc.nf.vs.rvf.down)
sc.NF_2_pRV_down_2_RVF_up <- setdiff(intersect(sc.nf.vs.prv.up, sc.prv.vs.rvf.down), sc.nf.vs.rvf.up)
sc.NF_2_pRV_flat_2_RVF_up <- setdiff(intersect(sc.nf.vs.rvf.up, sc.prv.vs.rvf.up), sc.nf.vs.prv.up)
sc.NF_2_pRV_flat_2_RVF_down <- setdiff(intersect(sc.nf.vs.rvf.down, sc.prv.vs.rvf.down), sc.nf.vs.prv.down)

write.csv(as.data.frame(sc.NF_2_pRV_up_2_RVF_up),file="~/Downloads/sc_myeloid_NF_2_pRV_up_2_RVF_up.csv")
write.csv(as.data.frame(sc.NF_2_pRV_down_2_RVF_down),file="~/Downloads/sc_myeloid_NF_2_pRV_down_2_RVF_down.csv")
write.csv(as.data.frame(sc.NF_2_pRV_up_2_RVF_down),file="~/Downloads/sc_myeloid_NF_2_pRV_up_2_RVF_down.csv")
write.csv(as.data.frame(sc.NF_2_pRV_down_2_RVF_up),file="~/Downloads/sc_myeloid_NF_2_pRV_down_2_RVF_up.csv")
write.csv(as.data.frame(sc.NF_2_pRV_flat_2_RVF_up),file="~/Downloads/sc_myeloid_NF_2_pRV_flat_2_RVF_up.csv")
write.csv(as.data.frame(sc.NF_2_pRV_flat_2_RVF_down),file="~/Downloads/sc_myeloid_NF_2_pRV_flat_2_RVF_down.csv")

#Mural only
Idents(M1) <- "Names_group"

b1 <- FindMarkers(M1,ident.1=c('PC_pRV','SM_pRV'),ident.2=c('PC_RVF','SM_RVF'),logfc.threshold=0)
b2 <- FindMarkers(M1,ident.1=c('PC_NF','SM_NF'),ident.2=c('PC_pRV','PC_pRV'),logfc.threshold=0)
b3 <- FindMarkers(M1,ident.1=c('PC_NF','SM_NF'),ident.2=c('PC_RVF','SM_RVF'),logfc.threshold=0)


sc.prv.vs.rvf.up <- rownames(subset(b1, p_val_adj<0.05 & avg_log2FC>0.1))
sc.prv.vs.rvf.down <- rownames(subset(b1, p_val_adj<0.05 & avg_log2FC < -0.1))

sc.nf.vs.prv.up <- rownames(subset(b2, p_val_adj<0.05 & avg_log2FC>0.1))
sc.nf.vs.prv.down <- rownames(subset(b2, p_val_adj<0.05 & avg_log2FC < -0.1))

sc.nf.vs.rvf.up <- rownames(subset(b3, p_val_adj<0.05 & avg_log2FC>0.1))
sc.nf.vs.rvf.down <- rownames(subset(b3, p_val_adj<0.05 & avg_log2FC < -0.1))

#Gene expression gradients
sc.NF_2_pRV_up_2_RVF_up <- intersect(intersect(sc.nf.vs.prv.down, sc.nf.vs.rvf.down), sc.prv.vs.rvf.down)
sc.NF_2_pRV_down_2_RVF_down <- intersect(intersect(sc.nf.vs.prv.up, sc.nf.vs.rvf.up), sc.prv.vs.rvf.up)
sc.NF_2_pRV_up_2_RVF_down <- setdiff(intersect(sc.nf.vs.prv.down, sc.prv.vs.rvf.up), sc.nf.vs.rvf.down)
sc.NF_2_pRV_down_2_RVF_up <- setdiff(intersect(sc.nf.vs.prv.up, sc.prv.vs.rvf.down), sc.nf.vs.rvf.up)
sc.NF_2_pRV_flat_2_RVF_up <- setdiff(intersect(sc.nf.vs.rvf.up, sc.prv.vs.rvf.up), sc.nf.vs.prv.up)
sc.NF_2_pRV_flat_2_RVF_down <- setdiff(intersect(sc.nf.vs.rvf.down, sc.prv.vs.rvf.down), sc.nf.vs.prv.down)

write.csv(as.data.frame(sc.NF_2_pRV_up_2_RVF_up),file="~/Downloads/sc_mural_NF_2_pRV_up_2_RVF_up.csv")
write.csv(as.data.frame(sc.NF_2_pRV_down_2_RVF_down),file="~/Downloads/sc_mural_NF_2_pRV_down_2_RVF_down.csv")
write.csv(as.data.frame(sc.NF_2_pRV_up_2_RVF_down),file="~/Downloads/sc_mural_NF_2_pRV_up_2_RVF_down.csv")
write.csv(as.data.frame(sc.NF_2_pRV_down_2_RVF_up),file="~/Downloads/sc_mural_NF_2_pRV_down_2_RVF_up.csv")
write.csv(as.data.frame(sc.NF_2_pRV_flat_2_RVF_up),file="~/Downloads/sc_mural_NF_2_pRV_flat_2_RVF_up.csv")
write.csv(as.data.frame(sc.NF_2_pRV_flat_2_RVF_down),file="~/Downloads/sc_mural_NF_2_pRV_flat_2_RVF_down.csv")

##############OLD


rvf.vs.rv <- topTreat(tfit, coef=2, n=Inf)
rvf.vs.rv.topgenes <- rownames(rvf.vs.rv)[1:100]

a = order(category)
i <- which(colnames(bulk) %in% rvf.vs.rv.topgenes[order(rvf.vs.rv[1:100,]$logFC)])
i <- i[match(rvf.vs.rv.topgenes[order(rvf.vs.rv[1:100,]$logFC)],colnames(bulk)[i])]
mycol <- colorpanel(1000,"blue","white","red")


pdf('~/Downloads/bulk_RVFvsNF_heatmap.pdf',width=10,height=10)
heatmap.2(as.matrix(t(bulk[a,i])), scale="row",
   labRow=colnames(bulk)[i], labCol=category[a], 
   col=mycol, margin=c(6,6),trace="none", density.info="none", lhei=c(1,10,3),lwid=c(1,10), dendrogram='none',
   Rowv = FALSE, Colv=FALSE,srtCol=90,lmat = rbind(c(0,3),c(2,1),c(0,4)))
dev.off()



f.vs.p <- topTreat(tfit, coef=3, n=Inf)
f.vs.p.topgenes <- rownames(f.vs.p)[1:100]

a = order(category)
i <- which(colnames(bulk) %in% f.vs.p.topgenes[order(f.vs.p[1:100,]$logFC)])
i <- i[match(f.vs.p.topgenes[order(f.vs.p[1:100,]$logFC)],colnames(bulk)[i])]
mycol <- colorpanel(1000,"blue","white","red")


pdf('~/Downloads/bulk_RVFvspRV_heatmap.pdf',width=10,height=10)
heatmap.2(as.matrix(t(bulk[a,i])), scale="row",
   labRow=colnames(bulk)[i], labCol=category[a], 
   col=mycol, margin=c(6,6),trace="none", density.info="none", lhei=c(1,10,3),lwid=c(1,10), dendrogram='none',
   Rowv = FALSE, Colv=FALSE,srtCol=90,lmat = rbind(c(0,3),c(2,1),c(0,4)))
dev.off()


##############OLD


#Get var feats
#M1<-readRDS('/Volumes/Extreme SSD/Final_Analysis/CellTypes/Post_R3_FINAL_with_counts.rds')
#M1 <- SetIdent(M1,value = 'Names')
#DefaultAssay(M1) <- 'SCT'
#var_feats <- VariableFeatures(M1)


#var_keep =  colnames(bulk) %in% var_feats
#bulk_sub <- bulk[,var_keep]

#pcDat <- prcomp(bulk_sub)
#pdf('~/Downloads/bulk_pca_varfeats.pdf',width=10,height=7.5)
#autoplot(pcDat,data = bulk.meta, colour="category", shape="sex", size=4) + theme_classic() + #theme(axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_blank(),ax#is.text.y=element_blank(),axis.title.x=element_text(size=16),axis.title.y=element_text(size=16#),legend.title=element_text(size=16),legend.text=element_text(size=16)) + #labs(color='Disease',shape='Sex')
#dev.off()

#####Limma

design <- model.matrix(~0+category)

contr.matrix <- makeContrasts(
   pRVvsNF =  categorypRV- categoryNF, 
   RVFvsNF = categoryRVF- categoryNF, 
   RVFvspRV = categoryRVF- categorypRV, 
   levels = colnames(design))
contr.matrix

fit <- lmFit(t(bulk), design)
fit2 <- contrasts.fit(fit, contr.matrix)
efit <- eBayes(fit2, trend=TRUE)
plotSA(efit, main="Final model: Mean-variance trend")

summary(decideTests(efit))
dtt <- decideTests(efit)

tfit <- treat(fit2, lfc=0.25)
dt <- decideTests(tfit)
summary(dt)

de.common <- which(dt[,1]!=0 & dt[,2]!=0)
length(de.common)

pdf('~/Downloads/bulk_overlap.pdf',width=10,height=10)

vennDiagram(dtt[,1:2], circle.col=c("turquoise", "salmon"))
dev.off()

pRV.vs.NF <- topTreat(tfit, coef=1, n=Inf)
head(pRV.vs.NF)

plotMD(tfit, column=1, status=dt[,1], main=colnames(tfit)[1], 
       xlim=c(-8,13))

library(gplots)

prv.vs.rv <- topTreat(tfit, coef=1, n=Inf)
prv.vs.rv.topgenes <- rownames(prv.vs.rv)[1:100]

a = order(category)
i <- which(colnames(bulk) %in% prv.vs.rv.topgenes[order(prv.vs.rv[1:100,]$logFC)])
i <- i[match(prv.vs.rv.topgenes[order(prv.vs.rv[1:100,]$logFC)],colnames(bulk)[i])]
mycol <- colorpanel(1000,"blue","white","red")


pdf('~/Downloads/bulk_pRVvsNF_heatmap.pdf',width=10,height=10)
heatmap.2(as.matrix(t(bulk[a,i])), scale="row",
   labRow=colnames(bulk)[i], labCol=category[a], 
   col=mycol, margin=c(6,6),trace="none", density.info="none", lhei=c(1,10,3),lwid=c(1,10), dendrogram='none',
   Rowv = FALSE, Colv=FALSE,srtCol=90,lmat = rbind(c(0,3),c(2,1),c(0,4)))
dev.off()


rvf.vs.rv <- topTreat(tfit, coef=2, n=Inf)
rvf.vs.rv.topgenes <- rownames(rvf.vs.rv)[1:100]

a = order(category)
i <- which(colnames(bulk) %in% rvf.vs.rv.topgenes[order(rvf.vs.rv[1:100,]$logFC)])
i <- i[match(rvf.vs.rv.topgenes[order(rvf.vs.rv[1:100,]$logFC)],colnames(bulk)[i])]
mycol <- colorpanel(1000,"blue","white","red")


pdf('~/Downloads/bulk_RVFvsNF_heatmap.pdf',width=10,height=10)
heatmap.2(as.matrix(t(bulk[a,i])), scale="row",
   labRow=colnames(bulk)[i], labCol=category[a], 
   col=mycol, margin=c(6,6),trace="none", density.info="none", lhei=c(1,10,3),lwid=c(1,10), dendrogram='none',
   Rowv = FALSE, Colv=FALSE,srtCol=90,lmat = rbind(c(0,3),c(2,1),c(0,4)))
dev.off()



f.vs.p <- topTreat(tfit, coef=3, n=Inf)
f.vs.p.topgenes <- rownames(f.vs.p)[1:100]

a = order(category)
i <- which(colnames(bulk) %in% f.vs.p.topgenes[order(f.vs.p[1:100,]$logFC)])
i <- i[match(f.vs.p.topgenes[order(f.vs.p[1:100,]$logFC)],colnames(bulk)[i])]
mycol <- colorpanel(1000,"blue","white","red")


pdf('~/Downloads/bulk_RVFvspRV_heatmap.pdf',width=10,height=10)
heatmap.2(as.matrix(t(bulk[a,i])), scale="row",
   labRow=colnames(bulk)[i], labCol=category[a], 
   col=mycol, margin=c(6,6),trace="none", density.info="none", lhei=c(1,10,3),lwid=c(1,10), dendrogram='none',
   Rowv = FALSE, Colv=FALSE,srtCol=90,lmat = rbind(c(0,3),c(2,1),c(0,4)))
dev.off()








write.csv(topTreat(tfit, coef=1, n=Inf),file="~/Downloads/bulk_pRV_vs_NF.csv")
write.csv(topTreat(tfit, coef=2, n=Inf),file="~/Downloads/bulk_RVF_vs_NF.csv")
write.csv(topTreat(tfit, coef=3, n=Inf),file="~/Downloads/bulk_RVF_vs_pRV.csv")

write.csv(topTreat(efit, coef=3, n=Inf),file="~/Downloads/bulk_RVF_vs_pRV_nofcthresh.csv")

pdf('~/Downloads/bulk_pRV_vs_NF.pdf',width=6,height=10)
EnhancedVolcano(topTreat(tfit, coef=1, n=Inf),lab = rownames(topTreat(tfit, coef=1, n=Inf)),x = 'logFC',y = 'P.Value',pCutoff=0.05/15959,FCcutoff=0.25,title = "", borderColour = 'black') +  ggplot2::coord_cartesian(xlim=c(-4, 4)) 
dev.off()

pdf('~/Downloads/bulk_RVF_vs_NF.pdf',width=6,height=10)
EnhancedVolcano(topTreat(tfit, coef=2, n=Inf),lab = rownames(topTreat(tfit, coef=2, n=Inf)),x = 'logFC',y = 'P.Value',pCutoff=0.05/15959,FCcutoff=0.25,title = "", borderColour = 'black') +  ggplot2::coord_cartesian(xlim=c(-4, 4)) 
dev.off()

pdf('~/Downloads/bulk_RVF_vs_pRV.pdf',width=6,height=10)
EnhancedVolcano(topTreat(efit, coef=3, n=Inf),lab = rownames(topTreat(efit, coef=3, n=Inf)),x = 'logFC',y = 'P.Value',pCutoff=0.05/15959,FCcutoff=0.1,title = "", borderColour = 'black') +  ggplot2::coord_cartesian(xlim=c(-4, 4)) 
dev.off()



#########Enrich

library(DOSE)
library(enrichR)
library(tidyverse)

websiteLive <- getOption("enrichR.live")
if (websiteLive) {
    listEnrichrSites()
    setEnrichrSite("Enrichr") # Human genes   
}


RVFvspRV = topTreat(efit, coef=3, n=Inf)

dbs <- c("ChEA_2022","WikiPathway_2023_Human","Reactome_2016","GO_Biological_Process_2023")
if (websiteLive) {
    enriched <- enrichr(rownames(subset(RVFvspRV,logFC > .1)), dbs)
	p1<- ggplot(enriched[[1]][order(enriched[[1]]$Combined.Score,decreasing=T),][rev(1:20),], (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + ggtitle('ChEA Up') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[1]][order(enriched[[1]]$Combined.Score,decreasing=T),][rev(1:20),]$Term," "), `[`, 1)))
	p2<- ggplot(enriched[[2]][order(enriched[[2]]$Combined.Score,decreasing=T),][rev(1:20),], (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + ggtitle('WikiPathways Up') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[2]][order(enriched[[2]]$Combined.Score,decreasing=T),][rev(1:20),]$Term," WP"), `[`, 1)))
	p3<- ggplot(enriched[[3]][order(enriched[[3]]$Combined.Score,decreasing=T),][rev(1:20),], (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + ggtitle('Reactome Up') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[3]][order(enriched[[3]]$Combined.Score,decreasing=T),][rev(1:20),]$Term," Homo sapiens"), `[`, 1)))
	p4<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:20),], (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + ggtitle('GO Biological Process Up') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:20),]$Term," \\(GO"), `[`, 1)))

    enriched <- enrichr(rownames(subset(RVFvspRV,logFC < -.1)), dbs)
	p5<- ggplot(enriched[[1]][order(enriched[[1]]$Combined.Score,decreasing=T),][rev(1:20),], (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + ggtitle('ChEA Down') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[1]][order(enriched[[1]]$Combined.Score,decreasing=T),][rev(1:20),]$Term," "), `[`, 1)))
	p6<- ggplot(enriched[[2]][order(enriched[[2]]$Combined.Score,decreasing=T),][rev(1:20),], (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + ggtitle('WikiPathways Down') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[2]][order(enriched[[2]]$Combined.Score,decreasing=T),][rev(1:20),]$Term," WP"), `[`, 1)))
	p7<- ggplot(enriched[[3]][order(enriched[[3]]$Combined.Score,decreasing=T),][rev(1:20),], (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + ggtitle('Reactome Down') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[3]][order(enriched[[3]]$Combined.Score,decreasing=T),][rev(1:20),]$Term," Homo sapiens"), `[`, 1)))
	p8<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:20),], (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + ggtitle('GO Biological Process Down') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:20),]$Term," \\(GO"), `[`, 1)))

library(patchwork)
library(cowplot)
pdf('~/Downloads/bulk_RVFvspRV_enrich.pdf',width=40,height=12)
plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,align="v",ncol=4)
dev.off()





RVFvsNF = topTreat(tfit, coef=2, n=Inf)

dbs <- c("ChEA_2022","WikiPathway_2023_Human","Reactome_2016","GO_Biological_Process_2023")
if (websiteLive) {
    enriched <- enrichr(rownames(subset(RVFvsNF,logFC > .25)), dbs)
	p1<- ggplot(enriched[[1]][order(enriched[[1]]$Combined.Score,decreasing=T),][rev(1:20),], (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + ggtitle('ChEA Up') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[1]][order(enriched[[1]]$Combined.Score,decreasing=T),][rev(1:20),]$Term," "), `[`, 1)))
	p2<- ggplot(enriched[[2]][order(enriched[[2]]$Combined.Score,decreasing=T),][rev(1:20),], (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + ggtitle('WikiPathways Up') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[2]][order(enriched[[2]]$Combined.Score,decreasing=T),][rev(1:20),]$Term," WP"), `[`, 1)))
	p3<- ggplot(enriched[[3]][order(enriched[[3]]$Combined.Score,decreasing=T),][rev(1:20),], (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + ggtitle('Reactome Up') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[3]][order(enriched[[3]]$Combined.Score,decreasing=T),][rev(1:20),]$Term," Homo sapiens"), `[`, 1)))
	p4<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:20),], (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + ggtitle('GO Biological Process Up') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:20),]$Term," \\(GO"), `[`, 1)))

    enriched <- enrichr(rownames(subset(RVFvsNF,logFC < -.25)), dbs)
	p5<- ggplot(enriched[[1]][order(enriched[[1]]$Combined.Score,decreasing=T),][rev(1:20),], (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + ggtitle('ChEA Down') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[1]][order(enriched[[1]]$Combined.Score,decreasing=T),][rev(1:20),]$Term," "), `[`, 1)))
	p6<- ggplot(enriched[[2]][order(enriched[[2]]$Combined.Score,decreasing=T),][rev(1:20),], (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + ggtitle('WikiPathways Down') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[2]][order(enriched[[2]]$Combined.Score,decreasing=T),][rev(1:20),]$Term," WP"), `[`, 1)))
	p7<- ggplot(enriched[[3]][order(enriched[[3]]$Combined.Score,decreasing=T),][rev(1:20),], (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + ggtitle('Reactome Down') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[3]][order(enriched[[3]]$Combined.Score,decreasing=T),][rev(1:20),]$Term," Homo sapiens"), `[`, 1)))
	p8<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:20),], (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + ggtitle('GO Biological Process Down') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:20),]$Term," \\(GO"), `[`, 1)))

library(patchwork)
library(cowplot)
pdf('~/Downloads/bulk_RVFvsNF_enrich.pdf',width=40,height=12)
plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,align="v",ncol=4)
dev.off()




pRVvsNF = topTreat(tfit, coef=1, n=Inf)

dbs <- c("ChEA_2022","WikiPathway_2023_Human","Reactome_2016","GO_Biological_Process_2023")
if (websiteLive) {
    enriched <- enrichr(rownames(subset(pRVvsNF,logFC > .25)), dbs)
	p1<- ggplot(enriched[[1]][order(enriched[[1]]$Combined.Score,decreasing=T),][rev(1:20),], (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + ggtitle('ChEA Up') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[1]][order(enriched[[1]]$Combined.Score,decreasing=T),][rev(1:20),]$Term," "), `[`, 1)))
	p2<- ggplot(enriched[[2]][order(enriched[[2]]$Combined.Score,decreasing=T),][rev(1:20),], (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + ggtitle('WikiPathways Up') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[2]][order(enriched[[2]]$Combined.Score,decreasing=T),][rev(1:20),]$Term," WP"), `[`, 1)))
	p3<- ggplot(enriched[[3]][order(enriched[[3]]$Combined.Score,decreasing=T),][rev(1:20),], (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + ggtitle('Reactome Up') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[3]][order(enriched[[3]]$Combined.Score,decreasing=T),][rev(1:20),]$Term," Homo sapiens"), `[`, 1)))
	p4<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:20),], (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + ggtitle('GO Biological Process Up') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:20),]$Term," \\(GO"), `[`, 1)))

    enriched <- enrichr(rownames(subset(pRVvsNF,logFC < -.25)), dbs)
	p5<- ggplot(enriched[[1]][order(enriched[[1]]$Combined.Score,decreasing=T),][rev(1:20),], (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + ggtitle('ChEA Down') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[1]][order(enriched[[1]]$Combined.Score,decreasing=T),][rev(1:20),]$Term," "), `[`, 1)))
	p6<- ggplot(enriched[[2]][order(enriched[[2]]$Combined.Score,decreasing=T),][rev(1:20),], (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + ggtitle('WikiPathways Down') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[2]][order(enriched[[2]]$Combined.Score,decreasing=T),][rev(1:20),]$Term," WP"), `[`, 1)))
	p7<- ggplot(enriched[[3]][order(enriched[[3]]$Combined.Score,decreasing=T),][rev(1:20),], (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + ggtitle('Reactome Down') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[3]][order(enriched[[3]]$Combined.Score,decreasing=T),][rev(1:20),]$Term," Homo sapiens"), `[`, 1)))
	p8<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:20),], (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + ggtitle('GO Biological Process Down') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:20),]$Term," \\(GO"), `[`, 1)))

library(patchwork)
library(cowplot)
pdf('~/Downloads/bulk_pRVvsNF_enrich.pdf',width=40,height=12)
plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,align="v",ncol=4)
dev.off()



#####Limma with etiology, age, sex

design <- model.matrix(~0+category+sex+disease+age)

contr.matrix <- makeContrasts(
   pRVvsNF =  categorypRV- categoryNF, 
   RVFvsNF = categoryRVF- categoryNF, 
   RVFvspRV = categoryRVF- categorypRV, 
   levels = colnames(design))
contr.matrix

fit <- lmFit(t(bulk), design)
fit2 <- contrasts.fit(fit, contr.matrix)
efit <- eBayes(fit2, trend=TRUE)
plotSA(efit, main="Final model: Mean-variance trend")

summary(decideTests(efit))
dtt <- decideTests(efit)

tfit <- treat(fit2, lfc=0.25)
dt <- decideTests(tfit)
summary(dt)

de.common <- which(dt[,1]!=0 & dt[,2]!=0)
length(de.common)

pdf('~/Downloads/bulk_overlap.pdf',width=10,height=10)

vennDiagram(dtt[,1:2], circle.col=c("turquoise", "salmon"))
dev.off()

pRV.vs.NF <- topTreat(tfit, coef=1, n=Inf)
head(pRV.vs.NF)

plotMD(tfit, column=1, status=dt[,1], main=colnames(tfit)[1], 
       xlim=c(-8,13))

library(gplots)

prv.vs.rv <- topTreat(tfit, coef=1, n=Inf)
prv.vs.rv.topgenes <- rownames(prv.vs.rv)[1:1000]

a = order(category)
i <- which(colnames(bulk) %in% prv.vs.rv.topgenes)
mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(as.matrix(bulk[a,i]), scale="row",
   labRow=colnames(bulk)[i], labCol=category[a], 
   col=mycol, trace="none", density.info="none", 
   margin=c(8,6), lhei=c(2,10), dendrogram="column")



write.csv(topTreat(tfit, coef=1, n=Inf),file="~/Downloads/bulk_pRV_vs_NF.csv")
write.csv(topTreat(tfit, coef=2, n=Inf),file="~/Downloads/bulk_RVF_vs_NF.csv")
write.csv(topTreat(tfit, coef=3, n=Inf),file="~/Downloads/bulk_RVF_vs_pRV.csv")

write.csv(topTreat(fit2, coef=3, n=Inf),file="~/Downloads/bulk_RVF_vs_pRV_nofcthresh.csv")

