library(Seurat)
library(hdWGCNA)
library(ggeasy)
library(harmony)



source('~/Downloads/hdWGCNA_TOM/spatial_functions.R')

#######################################
#############  FIGURE 4A  #############
#######################################
snLV <- readRDS('~/Downloads/hdWGCNA_TOM/Kory_reprocessed_all.rds')
snLV <- subset(snLV,tech=='SN')
all_signif <- c('M1','M2','M3','M4','M5','M8','M10','M11','M12','M14','M20','M25','M26','M28')

consensus_modules <- read.csv("~/Downloads/hdWGCNA_TOM/bulk_heart_modules.csv")
consensus_modules <- consensus_modules[,1:3]
consensus_modules <- subset(consensus_modules, gene_name %in% rownames(snLV))
# remove duplicate gene names
consensus_modules <- consensus_modules[match(unique(consensus_modules$gene_name), consensus_modules$gene_name),]


library(dplyr)
mapping <- labels2colors(1:100)

score_calc <- consensus_modules %>% group_by(module) %>% group_split()
module_colors <- unique(unlist(lapply(score_calc,'[[','module')))
module_colors <- paste0('M',match(module_colors,mapping))


snLV <- AddModuleScore(snLV,lapply(score_calc,'[[','gene_name'),name="module_score")

cols_current <- colnames(snLV@meta.data)
cols_current[startsWith(colnames(snLV@meta.data),'module_score')] <- paste0('module_',module_colors)
colnames(snLV@meta.data) <- cols_current


snLV <- SetIdent(snLV, value = "Names")


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'sc_seurat_dot_bulk2lv.pdf'), width=7.5, height=4)

p <- DotPlot(snLV,paste0('module_',
  c('M20','M5','M1','M3','M4','M8','M2','M12','M25','M26','M10','M28','M14','M11')),
group.by='Names',dot.min=0,col.min=0,col.max=2,idents=c("Cardiomyocytes","Endothelium","Fibroblasts","Myeloid","Pericytes","Smooth_Muscle")) +
  RotatedAxis() + ylab('')+ xlab('')+
  scale_color_gradient2(high='red', mid='grey95', low='blue') +
  theme(
    panel.border = element_rect(size=1,fill=NA, color='black'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()
) 
p
dev.off()


snLV <- SetIdent(snLV, value = "condition")

my_levels <- c("Donor","DCM")
Idents(snLV) <- factor(Idents(snLV), levels= my_levels)


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'sc_seurat_LV_trend.pdf'), width=6.5, height=2)

p <- DotPlot(snLV,paste0('module_',
  c('M20','M5','M1','M3','M4','M8','M2','M12','M25','M26','M10','M28','M14','M11')),dot.min=0,col.min=0,col.max=2) +
  RotatedAxis() + ylab('')+ xlab('')+
  scale_color_gradient2(high='red', mid='grey95', low='blue') +
  theme(
    panel.border = element_rect(size=1,fill=NA, color='black'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()
) 
p
dev.off()


#######################################
#############  FIGURE 4B  #############
#######################################
snLV <- readRDS('~/Downloads/hdWGCNA_TOM/Kory_reprocessed_all.rds')
snLV <- subset(snLV,tech=='SN')
snLV<-subset(snLV,Names=="Cardiomyocytes")

all_signif <- c('M1','M2','M3','M4','M5','M8','M10','M11','M12','M14','M20','M25','M26','M28')

consensus_modules <- read.csv("~/Downloads/hdWGCNA_TOM/bulk_heart_modules.csv")
consensus_modules <- consensus_modules[,1:3]
consensus_modules <- subset(consensus_modules, gene_name %in% rownames(snLV))
# remove duplicate gene names
consensus_modules <- consensus_modules[match(unique(consensus_modules$gene_name), consensus_modules$gene_name),]


library(dplyr)
mapping <- labels2colors(1:100)

score_calc <- consensus_modules %>% group_by(module) %>% group_split()
module_colors <- unique(unlist(lapply(score_calc,'[[','module')))
module_colors <- paste0('M',match(module_colors,mapping))


snLV <- AddModuleScore(snLV,lapply(score_calc,'[[','gene_name'),name="module_score")

cols_current <- colnames(snLV@meta.data)
cols_current[startsWith(colnames(snLV@meta.data),'module_score')] <- paste0('module_',module_colors)
colnames(snLV@meta.data) <- cols_current

which(module_colors %in% c('M10','M25','M26','M28'))


mito_genes <- unlist(lapply(score_calc[which(module_colors %in% c('M10','M25','M26','M28'))],'[[','gene_name'))
snLV <- AddModuleScore(snLV,list(mito_genes),name="module_mito")


snLV <- SetIdent(snLV, value = "condition")
pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'LV_vln_modules.pdf'), width=4.5, height=3)
VlnPlot(snLV,c('module_M2','module_mito1'),pt.size=0)
dev.off()


p <- DotPlot(snLV,paste0('module_',
  c('M2','M12','M25','M26','M10','M28')),dot.min=0,col.min=0,col.max=2) +
  RotatedAxis() + ylab('')+ xlab('')+
  scale_color_gradient2(high='red', mid='grey95', low='blue') +
  theme(
    panel.border = element_rect(size=1,fill=NA, color='black'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()
) 
p
dev.off()


#######################################
#############  FIGURE 4C  #############
#######################################
snLV <- readRDS('~/Downloads/hdWGCNA_TOM/Kory_reprocessed_all.rds')
snLV <- subset(snLV,tech=='SN')
snLV<-subset(snLV,Names=="Cardiomyocytes")
snLV <- SetIdent(snLV, value = "condition")
snLV <- PrepSCTFindMarkers(snLV)


bulk_modules <- read.csv("~/Downloads/hdWGCNA_TOM/bulk_heart_modules.csv")
bulk_modules$module <- match(bulk_modules$module,mapping)
combined_set <- data.frame()
mods_idx <- c(2,12,28,10,25,26)
for (i in mods_idx){
  key_genes <- subset(bulk_modules,module %in% c(i))$gene_name
  key_genes <- key_genes[key_genes %in% rownames(snLV)]
  gene_set <- FindMarkers(snLV, ident.1 = "DCM", ident.2 = "Donor",features=key_genes)
  #gene_set<-subset(gene_set,p_val_adj<0.05)
  gene_set$module <- paste0('M',i)
  gene_set$color <- mapping[i]
  if (length(combined_set) == 0){
    combined_set <- gene_set
  }
  else {
    combined_set <- rbind(combined_set,gene_set)
  }
}


M1 <- readRDS(file = "/Volumes/Extreme\ SSD/Final_Analysis/CellTypes/cm_subclust.rds")
new.cluster.ids <- c("Cm1","Cm2","Cm3","Cm4","Cm5","Cm6","Cm7","Cm8","Cm9","Cm10")
names(new.cluster.ids) <- levels(M1)
M1 <- RenameIdents(M1, new.cluster.ids)
M1$Subnames <- M1@active.ident
M1$SubNames_Groups <- paste(M1$Subnames,M1$group,sep='_')
M1 <- SetIdent(M1, value = "group")

bulk_modules <- read.csv("~/Downloads/hdWGCNA_TOM/bulk_heart_modules.csv")
bulk_modules$module <- match(bulk_modules$module,mapping)
combined_set_RV <- data.frame()
mods_idx <- c(2,12,28,10,25,26)
for (i in mods_idx){
  key_genes <- subset(bulk_modules,module %in% c(i))$gene_name
  key_genes <- key_genes[key_genes %in% rownames(M1)]
  gene_set <- FindMarkers(M1, ident.1 = "RVF", ident.2 = "NF",features=key_genes)
  #gene_set<-subset(gene_set,p_val_adj<0.05)
  gene_set$module <- paste0('M',i)
  gene_set$color <- mapping[i]
  if (length(combined_set_RV) == 0){
    combined_set_RV <- gene_set
  }
  else {
    combined_set_RV <- rbind(combined_set_RV,gene_set)
  }
}



overlap_genes = intersect(rownames(M1),rownames(snLV))
gene_set <- FindMarkers(M1, ident.1 = "RVF", ident.2 = "NF",features=overlap_genes)
gene_set$module <- paste0('All')
gene_set$color <- 'grey'
combined_set_RV <- rbind(combined_set_RV,gene_set)


gene_set <- FindMarkers(snLV, ident.1 = "DCM", ident.2 = "Donor",features=overlap_genes)
gene_set$module <- paste0('All')
gene_set$color <- 'grey'
combined_set <- rbind(combined_set,gene_set)

a<-subset(combined_set,module=='M2')
b <- subset(combined_set_RV,module=='M2')
shared <- intersect(rownames(a),rownames(b))
dataset <- data.frame(LV=a[shared,]$avg_log2FC,RV=b[shared,]$avg_log2FC)
rownames(dataset) <- shared

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'LV_vs_RV_CM_module_M2_dot.pdf'), width=6, height=8)
ggplot(dataset, aes(x = RV, y=LV)) + geom_point() + 
  geom_text_repel(label=rownames(dataset)) + theme_classic()
dev.off()



a<-subset(combined_set,module=='M12')
b <- subset(combined_set_RV,module=='M12')
shared <- intersect(rownames(a),rownames(b))
dataset <- data.frame(LV=a[shared,]$avg_log2FC,RV=b[shared,]$avg_log2FC)
rownames(dataset) <- shared

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'LV_vs_RV_CM_module_M12_dot.pdf'), width=8, height=6)
ggplot(dataset, aes(x = RV, y=LV)) + geom_point() + 
  geom_text_repel(label=rownames(dataset)) + theme_classic()
dev.off()


a<-subset(combined_set,module=='M25')
b <- subset(combined_set_RV,module=='M25')
shared <- intersect(rownames(a),rownames(b))
dataset <- data.frame(LV=a[shared,]$avg_log2FC,RV=b[shared,]$avg_log2FC)
rownames(dataset) <- shared

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'LV_vs_RV_CM_module_M25_dot.pdf'), width=8, height=6)
ggplot(dataset, aes(x = RV, y=LV)) + geom_point() + 
  geom_text_repel(label=rownames(dataset)) + theme_classic()
dev.off()


a<-subset(combined_set,module=='M10')
b <- subset(combined_set_RV,module=='M10')
shared <- intersect(rownames(a),rownames(b))
dataset <- data.frame(LV=a[shared,]$avg_log2FC,RV=b[shared,]$avg_log2FC)
rownames(dataset) <- shared

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'LV_vs_RV_CM_module_M10_dot.pdf'), width=8, height=6)
ggplot(dataset, aes(x = RV, y=LV)) + geom_point() + 
  geom_text_repel(label=rownames(dataset)) + theme_classic()
dev.off()

a<-subset(combined_set,module=='M26')
b <- subset(combined_set_RV,module=='M26')
shared <- intersect(rownames(a),rownames(b))
dataset <- data.frame(LV=a[shared,]$avg_log2FC,RV=b[shared,]$avg_log2FC)
rownames(dataset) <- shared

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'LV_vs_RV_CM_module_M26_dot.pdf'), width=8, height=6)
ggplot(dataset, aes(x = RV, y=LV)) + geom_point() + 
  geom_text_repel(label=rownames(dataset)) + theme_classic()
dev.off()


a<-subset(combined_set,module=='M28')
b <- subset(combined_set_RV,module=='M28')
shared <- intersect(rownames(a),rownames(b))
dataset <- data.frame(LV=a[shared,]$avg_log2FC,RV=b[shared,]$avg_log2FC)
rownames(dataset) <- shared

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'LV_vs_RV_CM_module_M28_dot.pdf'), width=8, height=6)
ggplot(dataset, aes(x = RV, y=LV)) + geom_point() + 
  geom_text_repel(label=rownames(dataset)) + theme_classic()
dev.off()


a<-subset(combined_set,module %in% c('M10','M25','M26','M28'))
b <- subset(combined_set_RV,module %in% c('M10','M25','M26','M28'))
shared <- intersect(rownames(a),rownames(b))
dataset <- data.frame(LV=a[shared,]$avg_log2FC,RV=b[shared,]$avg_log2FC)
rownames(dataset) <- shared
#cat(rownames(subset(dataset,LV<0 & RV<0)),sep='\n')

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'LV_vs_RV_CM_module_mito_dot.pdf'), width=6, height=8)
ggplot(dataset, aes(x = RV, y=LV)) + geom_point() + 
  geom_text_repel(label=rownames(dataset)) + theme_classic()
dev.off()

a<-subset(combined_set,module=='All')
b <- subset(combined_set_RV,module=='All')
shared <- intersect(rownames(a),rownames(b))
dataset <- data.frame(LV=a[shared,]$avg_log2FC,RV=b[shared,]$avg_log2FC)
rownames(dataset) <- shared
labs <- rownames(dataset)
labs[abs(dataset$LV - dataset$RV)<0.1] <- NA


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'LV_vs_RV_CM_module_All_dot.pdf'), width=6, height=8)
ggplot(dataset, aes(x = RV, y=LV)) + geom_point() + 
  geom_text_repel(label=labs,max.overlaps=10) + theme_classic()
dev.off()



#genes_diff <- abs(dataset$LV - dataset$RV) > 1
#genes_diff <- shared[genes_diff]

#genes_diff_sig <- a[genes_diff,]$p_val_adj < 0.05 & b[genes_diff,]$p_val_adj < 0.05
#genes_diff<-genes_diff[genes_diff_sig]

#LV_diff_up <- genes_diff[a[genes_diff,]$avg_log2FC > 0]
#RV_diff_up <- genes_diff[b[genes_diff,]$avg_log2FC > 0]


#######################################
#############  FIGURE 4D  #############
#######################################
library(viridis)
a<-subset(combined_set,module %in% c('M10','M25','M26','M28'))
b <- subset(combined_set_RV,module %in% c('M10','M25','M26','M28'))
shared <- intersect(rownames(a),rownames(b))
dataset <- data.frame(LV=a[shared,]$avg_log2FC,RV=b[shared,]$avg_log2FC)
rownames(dataset) <- shared
mito_genes <- rownames(subset(dataset,LV<0 & RV<0))

dbs <- c("ChEA_2022","WikiPathway_2023_Human","Reactome_2016","GO_Biological_Process_2023")

library(enrichR)
library(forcats)

parse_ratio <- function(ratio) {
    ratio <- sub("^\\s*", "", as.character(ratio))
    ratio <- sub("\\s*$", "", ratio)
    numerator <- as.numeric(sub("/\\d+$", "", ratio))
    denominator <- as.numeric(sub("^\\d+/", "", ratio))
    return(numerator/denominator)
}

wrapText <- function(x, len) {
    sapply(x, function(y) paste(strwrap(y, len), collapse = "\n"), USE.NAMES = FALSE)
}

enriched <- enrichr(mito_genes, dbs)
pdf('~/Downloads/hdWGCNA_TOM/CM_RV_LV_mito_enrichr_up.pdf',width=5,height=2.5)
p1<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:5),], 
  (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), 
  size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + 
  ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + 
  ggtitle('GO Biological Process Up') + 
  scale_y_discrete(labels= fct_inorder(
    wrapText(sapply(
      strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:5),]$Term," \\(GO"),
         `[`, 1),35))) + 
  theme(axis.text=element_text(colour="black"))+
  scale_color_stepsn(colors=rev(magma(256)))
p1
dev.off()

#######################################
#############  FIGURE 4E  #############
#######################################

a<-subset(combined_set,module %in% c('M2'))
b <- subset(combined_set_RV,module %in% c('M2'))
shared <- intersect(rownames(a),rownames(b))
dataset <- data.frame(LV=a[shared,]$avg_log2FC,RV=b[shared,]$avg_log2FC)
rownames(dataset) <- shared
M2_genes <- rownames(subset(dataset,LV>0 & RV>0))

enriched <- enrichr(M2_genes, dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('~/Downloads/hdWGCNA_TOM/CM_RV_LV_M2_enrichr_up.pdf',width=5,height=2.5)
p2<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:5),], 
  (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), 
  size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + 
  ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + 
  ggtitle('GO Biological Process Up') + 
  scale_y_discrete(labels= fct_inorder(
    wrapText(sapply(
      strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:5),]$Term," \\(GO"),
         `[`, 1),35))) + 
  theme(axis.text=element_text(colour="black"))+
  scale_color_stepsn(colors=rev(magma(256)))
p2
dev.off()

pdf('~/Downloads/hdWGCNA_TOM/CM_RV_LV_M&mito2_enrichr_up.pdf',width=6,height=5)
p1/p2
dev.off()

#######################################
#############  FIGURE 4F  #############
#######################################


M1 <- readRDS(file = "/Volumes/Extreme\ SSD/Final_Analysis/CellTypes/cm_subclust.rds")

new.cluster.ids <- c("Cm1","Cm2","Cm3","Cm4","Cm5","Cm6","Cm7","Cm8","Cm9","Cm10")
names(new.cluster.ids) <- levels(M1)
M1 <- RenameIdents(M1, new.cluster.ids)

M1$Subnames <- M1@active.ident
M1$SubNames_Groups <- paste(M1$Subnames,M1$group,sep='_')

M1 <- SetIdent(M1, value = "group")


aggregate_snLV <- AggregateExpression(snLV, group.by = c("orig.ident",'condition'), return.seurat = TRUE,assays='RNA')

aggregate_RV <- AggregateExpression(M1, group.by = c("patient",'group'), return.seurat = TRUE,assays='RNA')

LV_counts <- aggregate_snLV$RNA['counts']
RV_counts <- aggregate_RV$RNA['counts']
LV_meta <- aggregate_snLV@meta.data
RV_meta <- aggregate_RV@meta.data

colnames(LV_meta) <- c('patient','group')
colnames(RV_meta) <- c('temp','patient','group')
RV_meta$temp <- NULL
LV_meta$origin <- 'LV'
RV_meta$origin <- 'RV'


meta <- data.frame(rbind(LV_meta,RV_meta))

genes_shared <- intersect(rownames(LV_counts),rownames(RV_counts))

countmat<-cbind(LV_counts[genes_shared,],RV_counts[genes_shared,])

meta_merge <- meta
meta_merge$group[meta_merge$group == 'Donor'] = 'NF'


library("DESeq2")
ddsSE <- DESeqDataSetFromMatrix(countData = countmat,
                              colData = meta_merge,
                              design = ~origin+group)

ddsSE <- estimateSizeFactors(ddsSE)
idx <- rowSums( counts(ddsSE, normalized=TRUE) >= 5 ) >= 3
ddsSE <- ddsSE[idx,]

normalized_counts <- counts(ddsSE, normalized=TRUE)

library("sva")

mod  <- model.matrix(~origin+group, colData(ddsSE))
mod0 <- model.matrix(~group, colData(ddsSE))
svseq <- svaseq(normalized_counts, mod, mod0)


meta_merge$SV1 <- svseq$sv[,1]
meta_merge$SV2 <- svseq$sv[,2]
meta_merge$SV3 <- svseq$sv[,3]
meta_merge$SV4 <- svseq$sv[,4]
meta_merge$SV5 <- svseq$sv[,5]
meta_merge$SV6 <- svseq$sv[,6]
meta_merge$SV7 <- svseq$sv[,7]
meta_merge$SV8 <- svseq$sv[,8]
meta_merge$SV9 <- svseq$sv[,9]
meta_merge$SV10 <- svseq$sv[,10]
meta_merge$SV11 <- svseq$sv[,11]
meta_merge$SV12 <- svseq$sv[,12]
meta_merge$SV13 <- svseq$sv[,13]
meta_merge$SV14 <- svseq$sv[,14]

ddsSE <- DESeqDataSetFromMatrix(countData = countmat,
                              colData = meta_merge,
                              design = ~group+SV1+SV2+SV3+SV4+SV5+
                              SV6+SV7+SV8+SV9+SV10+SV11+SV12+SV13+SV14)



ddsSE <- estimateSizeFactors(ddsSE)
idx <- rowSums( counts(ddsSE, normalized=TRUE) >= 5 ) >= 12
ddsSE <- ddsSE[idx,]

normalized_counts <- counts(ddsSE, normalized=TRUE)



ddsSE <- DESeq(ddsSE)

vstSE <- vst(ddsSE,blind = FALSE)


mat <- assay(vstSE)
mm <- model.matrix(~group, colData(vstSE))
mat <- limma::removeBatchEffect(mat, covariates=colData(vstSE)[,4:18], design=mm)
assay(vstSE) <- mat

plotPCA(vstSE,intgroup=c("group"),ntop=100) + theme_classic() + theme(axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.title.x=element_text(size=24),axis.title.y=element_text(size=24),legend.title=element_text(size=24),legend.text=element_text(size=24)) + labs(color='Disease',shape='Sex') + geom_point(size=3.5)

pdf('~/Downloads/hdWGCNA_TOM/CM_RV_LV_vst_batchcorrect_PCA.pdf',width=6,height=2.5)
plotPCA(vstSE[,vstSE$group %in% c('pRV','RVF','DCM')],intgroup=c("group"),ntop=100) + theme_classic() + theme(axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.title.x=element_text(size=24),axis.title.y=element_text(size=24),legend.title=element_text(size=24),legend.text=element_text(size=24)) + labs(color='Disease',shape='Sex') + geom_point(size=3.5)
dev.off()






meta_merge_pooled <- meta_merge
meta_merge_pooled[meta_merge_pooled$group %in% c('pRV','RVF'),]$group = 'RVF'
ddsSE <- DESeqDataSetFromMatrix(countData = countmat,
                              colData = meta_merge_pooled,
                              design = ~group+SV1+SV2+SV3+SV4+SV5+
                              SV6+SV7+SV8+SV9+SV10+SV11+SV12+SV13+SV14)



ddsSE <- estimateSizeFactors(ddsSE)
idx <- rowSums( counts(ddsSE, normalized=TRUE) >= 5 ) >= 12
ddsSE <- ddsSE[idx,]
normalized_counts <- counts(ddsSE, normalized=TRUE)
ddsSE <- DESeq(ddsSE)
vstSE <- vst(ddsSE,blind = FALSE)
mat <- assay(vstSE)
mm <- model.matrix(~group, colData(vstSE))
mat <- limma::removeBatchEffect(mat, covariates=colData(vstSE)[,4:18], design=mm)
assay(vstSE) <- mat

#pdf('~/Downloads/hdWGCNA_TOM/CM_RV_LV_vst_batchcorrect_PCA_temp.pdf',width=6,height=2.5)

plotPCA(vstSE[,vstSE$group %in% c('RVF','DCM')],intgroup=c("group"),ntop=100) + theme_classic() + theme(axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.title.x=element_text(size=24),axis.title.y=element_text(size=24),legend.title=element_text(size=24),legend.text=element_text(size=24)) + labs(color='Disease',shape='Sex') + geom_point(size=3.5)

#dev.off()
dcm.vs.rvf <- lfcShrink(ddsSE,contrast=c('group','DCM','RVF'), type="ashr")

dcm.vs.rvf.signif <- subset(dcm.vs.rvf,padj<0.1)


lv <- subset(dcm.vs.rvf.signif,log2FoldChange>0.1)
rv <- subset(dcm.vs.rvf.signif,log2FoldChange<0.1)


lv.genes <- rownames(lv)
lv.genes <- lv.genes[!str_detect(lv.genes,'.1') & !str_detect(lv.genes,'.2')]
cat(lv.genes,sep='\n')

rv.genes <- rownames(rv)
rv.genes <- rv.genes[!str_detect(rv.genes,'.1') & !str_detect(rv.genes,'.2') & !str_detect(rv.genes,'.3') & !str_detect(rv.genes,'.4')]
cat(rv.genes,sep='\n')



library(EnhancedVolcano)
pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'CM_LV_RV_DESeq.pdf'), width=8, height=6)

EnhancedVolcano(subset(dcm.vs.rvf,baseMean>100),lab=rownames(subset(dcm.vs.rvf,baseMean>100)),
  x='log2FoldChange',y='padj',
  FCcutoff = 0.1,pCutoff=0.05,ylim=c(0,3)) + coord_flip()
dev.off()





pdf('~/Downloads/hdWGCNA_TOM/CM_RV_LV_vst_batchcorrect_PCA_temp.pdf',width=6,height=2.5)
FeaturePlot(M1,'ZSWIM6')

dev.off()


#######################################
#############  FIGURE 4H  #############
#######################################


arvm <- read.csv('~/Downloads/hdWGCNA_TOM/ARVM_RV_vs_LV.csv')


RV_genes <- unlist(lapply(subset(arvm,logFC>0)$SYMBOL,toupper))
LV_genes <- unlist(lapply(subset(arvm,logFC<0)$SYMBOL,toupper))



pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'arvm_RV_LV.pdf'), width=8, height=6)

EnhancedVolcano(arvm,lab=toupper(arvm$SYMBOL),
  x='logFC',y='adj.P.Val',
  FCcutoff = 0.1,pCutoff=0.05,ylim=c(0,3)) + coord_flip()
dev.off()

dbs <- c("ChEA_2022","WikiPathway_2023_Human","Reactome_2016","GO_Biological_Process_2023")

library(enrichR)
library(forcats)

parse_ratio <- function(ratio) {
    ratio <- sub("^\\s*", "", as.character(ratio))
    ratio <- sub("\\s*$", "", ratio)
    numerator <- as.numeric(sub("/\\d+$", "", ratio))
    denominator <- as.numeric(sub("^\\d+/", "", ratio))
    return(numerator/denominator)
}

wrapText <- function(x, len) {
    sapply(x, function(y) paste(strwrap(y, len), collapse = "\n"), USE.NAMES = FALSE)
}

enriched <- enrichr(RV_genes, dbs)
pdf('~/Downloads/hdWGCNA_TOM/ARVM_RV_enrichr_up.pdf',width=6,height=2.5)
p1<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:5),], 
  (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), 
  size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + 
  ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + 
  ggtitle('GO Biological Process Up') + 
  scale_y_discrete(labels= fct_inorder(
    wrapText(sapply(
      strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:5),]$Term," \\(GO"),
         `[`, 1),35))) + 
  theme(axis.text=element_text(colour="black"))+
  scale_color_stepsn(colors=rev(magma(256)))
p1
dev.off()


enriched <- enrichr(LV_genes, dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('~/Downloads/hdWGCNA_TOM/ARVM_LV_enrichr_up.pdf',width=5,height=2.5)
p2<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:5),], 
  (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), 
  size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + 
  ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + 
  ggtitle('GO Biological Process Up') + 
  scale_y_discrete(labels= fct_inorder(
    wrapText(sapply(
      strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:5),]$Term," \\(GO"),
         `[`, 1),35))) + 
  theme(axis.text=element_text(colour="black"))+
  scale_color_stepsn(colors=rev(magma(256)))
p2
dev.off()

pdf('~/Downloads/hdWGCNA_TOM/ARVM_RV_LV_enrichr_up.pdf',width=6,height=5)
p1/p2
dev.off()


#######################################
#############  FIGURE 4I  #############
#######################################

M2.genes <- score_calc[which(module_colors %in% c("M2"))][[1]]$gene_name

M2.arvm <- arvm[toupper(arvm$SYMBOL) %in% M2.genes,]

mito.genes <- unlist(lapply(score_calc[which(module_colors %in% c("M10","M25","M26","M28"))],'[[',1))

mito.arvm <- arvm[toupper(arvm$SYMBOL) %in% mito.genes,]

M12.genes <- score_calc[which(module_colors %in% c("M12"))][[1]]$gene_name

M12.arvm <- arvm[toupper(arvm$SYMBOL) %in% M12.genes,]

