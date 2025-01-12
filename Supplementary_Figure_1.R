
library(reticulate)
library(ggfortify)
library(edgeR)
library(RColorBrewer)
library(EnhancedVolcano)
library(DESeq2)
library(tximport)
library(biomaRt)
library("sva")
library(Seurat)

#######################################
#############  FIGURE S1A  #############
#######################################
bulk <- read.csv('/Volumes/RV_RNAseq/RV_snRNAseq/Final_Analysis/BulkRNA/counts.csv')
meta <- read.csv('/Volumes/RV_RNAseq/RV_snRNAseq/Final_Analysis/BulkRNA/metadata.csv')
toDel <- seq(1,dim(meta)[1],2)
meta <- meta[-toDel,]

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

res <- getBM(attributes = c('ensembl_transcript_id_version',                              'ensembl_gene_id',                              
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


pdf('~/Downloads/hdWGCNA_TOM/bulk_pRV_vs_RVF_volcano.pdf',width=8,height=10)
EnhancedVolcano(prv.vs.rvf,lab = rownames(prv.vs.rvf),x = 'log2FoldChange',y = 'padj',pCutoff=0.1,FCcutoff=0.25,title = "", borderColour = 'black') +  ggplot2::coord_cartesian(xlim=c(-1, 1)) 
dev.off()

pdf('~/Downloads/hdWGCNA_TOM/bulk_NF_vs_RVF_volcano.pdf',width=8,height=10)
EnhancedVolcano(nf.vs.rvf,lab = rownames(nf.vs.rvf),x = 'log2FoldChange',y = 'padj',pCutoff=0.1,FCcutoff=0.25,title = "", borderColour = 'black') +  ggplot2::coord_cartesian(xlim=c(-6, 6)) 
dev.off()

pdf('~/Downloads/hdWGCNA_TOM/bulk_NF_vs_pRV_volcano.pdf',width=8,height=10)
EnhancedVolcano(nf.vs.prv,lab = rownames(nf.vs.prv),x = 'log2FoldChange',y = 'padj',pCutoff=0.1,FCcutoff=0.25,title = "", borderColour = 'black') +  ggplot2::coord_cartesian(xlim=c(-6, 6)) 
dev.off()


#######################################
#############  FIGURE S1B  #############
#######################################

shared <- intersect(rownames(nf.vs.prv),rownames(nf.vs.rvf))
dataset <- data.frame(pRV=nf.vs.prv[shared,]$log2FoldChange,RVF=nf.vs.rvf[shared,]$log2FoldChange)
rownames(dataset) <- shared

labs <- rownames(dataset)
labs[prv.vs.rvf$padj > 0.05 | is.na(prv.vs.rvf$padj) | prv.vs.rvf$log2FoldChange<0.25] = ""

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'pRV_vs_RVF_logFC_dot.pdf'), width=4, height=5)
ggplot(dataset, aes(x = pRV, y=RVF)) + geom_point(alpha=0.05) + 
  geom_text_repel(label=labs,max.overlaps=Inf,point.size = NA,force=10) + theme_classic() + xlim(c(-5,5)) + ylim(c(-5,5))
dev.off()


temp <- data.frame(nf.vs.rvf) %>% arrange(desc(log2FoldChange))
write.csv(temp,'~/Downloads/hdWGCNA_TOM/NF_vs_RVF_deseq.csv')

temp <- data.frame(prv.vs.rvf) %>% arrange(desc(log2FoldChange))
write.csv(temp,'~/Downloads/hdWGCNA_TOM/pRV_vs_RVF_deseq.csv')

temp <- data.frame(subset(prv.vs.rvf,padj<0.05)) %>% arrange(desc(log2FoldChange))
write.csv(temp,'~/Downloads/hdWGCNA_TOM/pRV_vs_RVF_deseq_sig.csv')


#######################################
############  FIGURE S1C/D  ############
#######################################

##########Heatmaps
library(gplots)

#pRV vs RVF
#filt <- prv.vs.rvf$padj<0.1 & rownames(prv.vs.rvf) %in% pRV.only.signif
filt <- prv.vs.rvf$padj<0.1
filt[is.na(filt)]=FALSE
prv.vs.rvf.filt <- prv.vs.rvf[filt,]
prv.vs.rvf.filt <- prv.vs.rvf.filt[order(prv.vs.rvf.filt$padj),]
prv.vs.rvf.topgenes <- rownames(prv.vs.rvf.filt)[1:100]

a = order(category)
i <- which(rownames(assay(vstSE)) %in% prv.vs.rvf.topgenes[order(prv.vs.rvf.filt[1:100,]$log2FoldChange)])
i <- i[match(prv.vs.rvf.topgenes[order(prv.vs.rvf.filt[1:100,]$log2FoldChange)],rownames(assay(vstSE))[i])]
mycol <- colorpanel(1000,"blue","white","red")


pdf('~/Downloads/hdWGCNA_TOM/bulk_pRVvsRVF_heatmap_tight.pdf',width=10,height=10)
heatmap.2(as.matrix(assay(vstSE)[i,a]), scale="row",
   labRow=rownames(assay(vstSE))[i], labCol=category[a], 
   col=mycol, margin=c(6,6),trace="none", density.info="none", lhei=c(1,10,3),lwid=c(1,10), dendrogram='none',breaks = seq(-4, 4, length.out = 1001),
   Rowv = FALSE, Colv=FALSE,srtCol=90,lmat = rbind(c(0,3),c(2,1),c(0,4)))
dev.off()


#NF vs RVF

filt <- nf.vs.rvf$padj<0.1
filt[is.na(filt)]=FALSE
nf.vs.rvf.filt <- nf.vs.rvf[filt,]
nf.vs.rvf.filt <- nf.vs.rvf.filt[order(nf.vs.rvf.filt$padj),]
nf.vs.rvf.topgenes <- rownames(nf.vs.rvf.filt)[1:100]

a = order(category)
i <- which(rownames(assay(vstSE)) %in% nf.vs.rvf.topgenes[order(nf.vs.rvf.filt[1:100,]$log2FoldChange)])
i <- i[match(nf.vs.rvf.topgenes[order(nf.vs.rvf.filt[1:100,]$log2FoldChange)],rownames(assay(vstSE))[i])]
mycol <- colorpanel(1000,"blue","white","red")


pdf('~/Downloads/hdWGCNA_TOM/bulk_NFvsRVF_heatmap_tight.pdf',width=10,height=10)
heatmap.2(as.matrix(assay(vstSE)[i,a]), scale="row",
   labRow=rownames(assay(vstSE))[i], labCol=category[a], 
   col=mycol, margin=c(6,6),trace="none", density.info="none", lhei=c(1,10,3),lwid=c(1,10), dendrogram='none',breaks = seq(-4, 4, length.out = 1001),
   Rowv = FALSE, Colv=FALSE,srtCol=90,lmat = rbind(c(0,3),c(2,1),c(0,4)))
dev.off()

#NF vs pRV

filt <- nf.vs.prv$padj<0.1
filt[is.na(filt)]=FALSE
nf.vs.prv.filt <- nf.vs.prv[filt,]
nf.vs.prv.filt <- nf.vs.prv.filt[order(nf.vs.prv.filt$padj),]
nf.vs.prv.topgenes <- rownames(nf.vs.prv.filt)[1:100]

a = order(category)
i <- which(rownames(assay(vstSE)) %in% nf.vs.prv.topgenes[order(nf.vs.prv.filt[1:100,]$log2FoldChange)])
i <- i[match(nf.vs.prv.topgenes[order(nf.vs.prv.filt[1:100,]$log2FoldChange)],rownames(assay(vstSE))[i])]
mycol <- colorpanel(1000,"blue","white","red")


pdf('~/Downloads/hdWGCNA_TOM/bulk_NFvspRV_heatmap_tight.pdf',width=10,height=10)
heatmap.2(as.matrix(assay(vstSE)[i,a]), scale="row",
   labRow=rownames(assay(vstSE))[i], labCol=category[a], 
   col=mycol, margin=c(6,6),trace="none", density.info="none", lhei=c(1,10,3),lwid=c(1,10), dendrogram='none',breaks = seq(-4, 4, length.out = 1001),
   Rowv = FALSE, Colv=FALSE,srtCol=90,lmat = rbind(c(0,3),c(2,1),c(0,4)))
dev.off()


RVF.signif <- setdiff(signif2,signif3)
pRV.signif <- setdiff(signif3,signif2)
both.signif <- setdiff(setdiff(signif2,RVF.signif),pRV.signif)
pRV.only.signif <- intersect(signif1,union(RVF.signif, pRV.signif))

nf.vs.rvf.up <- rownames(subset(nf.vs.rvf,padj<0.1 & log2FoldChange>0))
nf.vs.rvf.down <- rownames(subset(nf.vs.rvf,padj<0.1 & log2FoldChange<0))
nf.vs.prv.up <- rownames(subset(nf.vs.prv,padj<0.1 & log2FoldChange>0))
nf.vs.prv.down <- rownames(subset(nf.vs.prv,padj<0.1 & log2FoldChange<0))
prv.vs.rvf.up <- rownames(subset(prv.vs.rvf,padj<0.1 & log2FoldChange>0))
prv.vs.rvf.down <- rownames(subset(prv.vs.rvf,padj<0.1 & log2FoldChange<0))

#Gene expression gradients
NF_2_pRV_up_2_RVF_up <- intersect(intersect(nf.vs.prv.down, nf.vs.rvf.down), prv.vs.rvf.down)
NF_2_pRV_down_2_RVF_down <- intersect(intersect(nf.vs.prv.up, nf.vs.rvf.up), prv.vs.rvf.up)
NF_2_pRV_up_2_RVF_down <- setdiff(intersect(nf.vs.prv.down, prv.vs.rvf.up), nf.vs.rvf.down)
NF_2_pRV_down_2_RVF_up <- setdiff(intersect(nf.vs.prv.up, prv.vs.rvf.down), nf.vs.rvf.up)
NF_2_pRV_flat_2_RVF_up <- setdiff(intersect(nf.vs.rvf.up, prv.vs.rvf.up), nf.vs.prv.up)
NF_2_pRV_flat_2_RVF_down <- setdiff(intersect(nf.vs.rvf.down, prv.vs.rvf.down), nf.vs.prv.down)



####
#NF_2_pRV_up_2_RVF_up
#NF_2_pRV_down_2_RVF_down 
#NF_2_pRV_up_2_RVF_down 
#NF_2_pRV_down_2_RVF_up

#pRV vs RVF
filt <- rownames(prv.vs.rvf) %in% NF_2_pRV_up_2_RVF_up
filt[is.na(filt)]=FALSE
prv.vs.rvf.filt <- prv.vs.rvf[filt,]
prv.vs.rvf.filt <- prv.vs.rvf.filt[order(prv.vs.rvf.filt$padj),]
prv.vs.rvf.topgenes <- rownames(prv.vs.rvf.filt)[1:100]

a = order(category)
i <- which(rownames(assay(vstSE)) %in% prv.vs.rvf.topgenes[order(prv.vs.rvf.filt[1:100,]$log2FoldChange)])
i <- i[match(prv.vs.rvf.topgenes[order(prv.vs.rvf.filt[1:100,]$log2FoldChange)],rownames(assay(vstSE))[i])]
mycol <- colorpanel(1000,"blue","white","red")


pdf('~/Downloads/hdWGCNA_TOM/bulk_pRVvsRVF_heatmap_module_NF_2_pRV_up_2_RVF_up.pdf',width=10,height=10)
heatmap.2(as.matrix(assay(vstSE)[i,a]), scale="row",
   labRow=rownames(assay(vstSE))[i], labCol=category[a], 
   col=mycol, margin=c(6,6),trace="none", density.info="none", lhei=c(1,10,3),lwid=c(1,10), dendrogram='none',breaks = seq(-4, 4, length.out = 1001),
   Rowv = FALSE, Colv=FALSE,srtCol=90,lmat = rbind(c(0,3),c(2,1),c(0,4)))
dev.off()

#
filt <- rownames(prv.vs.rvf) %in% NF_2_pRV_down_2_RVF_down
filt[is.na(filt)]=FALSE
prv.vs.rvf.filt <- prv.vs.rvf[filt,]
prv.vs.rvf.filt <- prv.vs.rvf.filt[order(prv.vs.rvf.filt$padj),]
prv.vs.rvf.topgenes <- rownames(prv.vs.rvf.filt)[1:100]

a = order(category)
i <- which(rownames(assay(vstSE)) %in% prv.vs.rvf.topgenes[order(prv.vs.rvf.filt[1:100,]$log2FoldChange)])
i <- i[match(prv.vs.rvf.topgenes[order(prv.vs.rvf.filt[1:100,]$log2FoldChange)],rownames(assay(vstSE))[i])]
mycol <- colorpanel(1000,"blue","white","red")


pdf('~/Downloads/hdWGCNA_TOM/bulk_pRVvsRVF_heatmap_module_NF_2_pRV_down_2_RVF_down.pdf',width=10,height=10)
heatmap.2(as.matrix(assay(vstSE)[i,a]), scale="row",
   labRow=rownames(assay(vstSE))[i], labCol=category[a], 
   col=mycol, margin=c(6,6),trace="none", density.info="none", lhei=c(1,10,3),lwid=c(1,10), dendrogram='none',breaks = seq(-4, 4, length.out = 1001),
   Rowv = FALSE, Colv=FALSE,srtCol=90,lmat = rbind(c(0,3),c(2,1),c(0,4)))
dev.off()
#

filt <- rownames(prv.vs.rvf) %in% NF_2_pRV_up_2_RVF_down
filt[is.na(filt)]=FALSE
prv.vs.rvf.filt <- prv.vs.rvf[filt,]
prv.vs.rvf.filt <- prv.vs.rvf.filt[order(prv.vs.rvf.filt$padj),]
prv.vs.rvf.topgenes <- rownames(prv.vs.rvf.filt)[1:100]

a = order(category)
i <- which(rownames(assay(vstSE)) %in% prv.vs.rvf.topgenes[order(prv.vs.rvf.filt[1:100,]$log2FoldChange)])
i <- i[match(prv.vs.rvf.topgenes[order(prv.vs.rvf.filt[1:100,]$log2FoldChange)],rownames(assay(vstSE))[i])]
mycol <- colorpanel(1000,"blue","white","red")


pdf('~/Downloads/hdWGCNA_TOM/bulk_pRVvsRVF_heatmap_module_NF_2_pRV_up_2_RVF_down.pdf',width=10,height=10)
heatmap.2(as.matrix(assay(vstSE)[i,a]), scale="row",
   labRow=rownames(assay(vstSE))[i], labCol=category[a], 
   col=mycol, margin=c(6,6),trace="none", density.info="none", lhei=c(1,10,3),lwid=c(1,10), dendrogram='none',breaks = seq(-4, 4, length.out = 1001),
   Rowv = FALSE, Colv=FALSE,srtCol=90,lmat = rbind(c(0,3),c(2,1),c(0,4)))
dev.off()

#
filt <- rownames(prv.vs.rvf) %in% NF_2_pRV_down_2_RVF_up
filt[is.na(filt)]=FALSE
prv.vs.rvf.filt <- prv.vs.rvf[filt,]
prv.vs.rvf.filt <- prv.vs.rvf.filt[order(prv.vs.rvf.filt$padj),]
prv.vs.rvf.topgenes <- rownames(prv.vs.rvf.filt)[1:100]

a = order(category)
i <- which(rownames(assay(vstSE)) %in% prv.vs.rvf.topgenes[order(prv.vs.rvf.filt[1:100,]$log2FoldChange)])
i <- i[match(prv.vs.rvf.topgenes[order(prv.vs.rvf.filt[1:100,]$log2FoldChange)],rownames(assay(vstSE))[i])]
mycol <- colorpanel(1000,"blue","white","red")


pdf('~/Downloads/hdWGCNA_TOM/bulk_pRVvsRVF_heatmap_module_NF_2_pRV_down_2_RVF_up.pdf',width=10,height=10)
heatmap.2(as.matrix(assay(vstSE)[i,a]), scale="row",
   labRow=rownames(assay(vstSE))[i], labCol=category[a], 
   col=mycol, margin=c(6,6),trace="none", density.info="none", lhei=c(1,10,3),lwid=c(1,10), dendrogram='none',breaks = seq(-4, 4, length.out = 1001),
   Rowv = FALSE, Colv=FALSE,srtCol=90,lmat = rbind(c(0,3),c(2,1),c(0,4)))
dev.off()

#
filt <- rownames(prv.vs.rvf) %in% NF_2_pRV_flat_2_RVF_up
filt[is.na(filt)]=FALSE
prv.vs.rvf.filt <- prv.vs.rvf[filt,]
prv.vs.rvf.filt <- prv.vs.rvf.filt[order(prv.vs.rvf.filt$padj),]
prv.vs.rvf.topgenes <- rownames(prv.vs.rvf.filt)[1:100]

a = order(category)
i <- which(rownames(assay(vstSE)) %in% prv.vs.rvf.topgenes[order(prv.vs.rvf.filt[1:100,]$log2FoldChange)])
i <- i[match(prv.vs.rvf.topgenes[order(prv.vs.rvf.filt[1:100,]$log2FoldChange)],rownames(assay(vstSE))[i])]
mycol <- colorpanel(1000,"blue","white","red")


pdf('~/Downloads/hdWGCNA_TOM/bulk_pRVvsRVF_heatmap_module_NF_2_pRV_flat_2_RVF_up.pdf',width=10,height=10)
heatmap.2(as.matrix(assay(vstSE)[i,a]), scale="row",
   labRow=rownames(assay(vstSE))[i], labCol=category[a], 
   col=mycol, margin=c(6,6),trace="none", density.info="none", lhei=c(1,10,3),lwid=c(1,10), dendrogram='none',breaks = seq(-4, 4, length.out = 1001),
   Rowv = FALSE, Colv=FALSE,srtCol=90,lmat = rbind(c(0,3),c(2,1),c(0,4)))
dev.off()

#
filt <- rownames(prv.vs.rvf) %in% NF_2_pRV_flat_2_RVF_down
filt[is.na(filt)]=FALSE
prv.vs.rvf.filt <- prv.vs.rvf[filt,]
prv.vs.rvf.filt <- prv.vs.rvf.filt[order(prv.vs.rvf.filt$padj),]
prv.vs.rvf.topgenes <- rownames(prv.vs.rvf.filt)[1:100]

a = order(category)
i <- which(rownames(assay(vstSE)) %in% prv.vs.rvf.topgenes[order(prv.vs.rvf.filt[1:100,]$log2FoldChange)])
i <- i[match(prv.vs.rvf.topgenes[order(prv.vs.rvf.filt[1:100,]$log2FoldChange)],rownames(assay(vstSE))[i])]
mycol <- colorpanel(1000,"blue","white","red")


pdf('~/Downloads/hdWGCNA_TOM/bulk_pRVvsRVF_heatmap_module_NF_2_pRV_flat_2_RVF_down.pdf',width=10,height=10)
heatmap.2(as.matrix(assay(vstSE)[i,a]), scale="row",
   labRow=rownames(assay(vstSE))[i], labCol=category[a], 
   col=mycol, margin=c(6,6),trace="none", density.info="none", lhei=c(1,10,3),lwid=c(1,10), dendrogram='none',breaks = seq(-4, 4, length.out = 1001),
   Rowv = FALSE, Colv=FALSE,srtCol=90,lmat = rbind(c(0,3),c(2,1),c(0,4)))
dev.off()

#######################################
############  FIGURE S1E  #############
#######################################
source('~/Downloads/hdWGCNA_TOM/spatial_functions.R')
library(hdWGCNA)


bulk_modules <- read.csv("~/Downloads/hdWGCNA_TOM/bulk_heart_modules.csv")

mapping <- labels2colors(1:100)
bulk_modules$module <- match(bulk_modules$module,mapping)
levels(bulk_modules$module) <- c(1:29)

idx<-match(NF_2_pRV_up_2_RVF_up,bulk_modules$gene_name)
idx<-idx[!is.na(idx)]
p1 <- table(bulk_modules[idx,]$module)/sum(table(bulk_modules[idx,]$module))*100
cat(subset(bulk_modules[idx,],module==2)$gene_name,sep='\n')
cat(subset(bulk_modules[idx,],module==3)$gene_name,sep='\n')
cat(subset(bulk_modules[idx,],module==8)$gene_name,sep='\n')
cat(subset(bulk_modules[idx,],module==15)$gene_name,sep='\n')


idx<-match(NF_2_pRV_down_2_RVF_down,bulk_modules$gene_name)
idx<-idx[!is.na(idx)]
p2 <- table(bulk_modules[idx,]$module)/sum(table(bulk_modules[idx,]$module))*100
#table(bulk_modules$module)/sum(table(bulk_modules$module)) * 100
cat(subset(bulk_modules[idx,],module==1)$gene_name,sep='\n') #<-------- interesting
cat(subset(bulk_modules[idx,],module==10)$gene_name,sep='\n')
cat(subset(bulk_modules[idx,],module==18)$gene_name,sep='\n')
cat(subset(bulk_modules[idx,],module==25)$gene_name,sep='\n')
cat(subset(bulk_modules[idx,],module==26)$gene_name,sep='\n')


idx<-match(NF_2_pRV_flat_2_RVF_up,bulk_modules$gene_name)
idx<-idx[!is.na(idx)]
p3 <- table(bulk_modules[idx,]$module)/sum(table(bulk_modules[idx,]$module))*100
cat(subset(bulk_modules[idx,],module==10)$gene_name,sep='\n')

idx<-match(NF_2_pRV_flat_2_RVF_down,bulk_modules$gene_name)
idx<-idx[!is.na(idx)]
p4 <- table(bulk_modules[idx,]$module)/sum(table(bulk_modules[idx,]$module))*100
cat(subset(bulk_modules[idx,],module==3)$gene_name,sep='\n')

idx<-match(NF_2_pRV_up_2_RVF_down,bulk_modules$gene_name)
idx<-idx[!is.na(idx)]
p5 <-table(bulk_modules[idx,]$module)/sum(table(bulk_modules[idx,]$module))*100

idx<-match(NF_2_pRV_down_2_RVF_up,bulk_modules$gene_name)
idx<-idx[!is.na(idx)]
p6 <- table(bulk_modules[idx,]$module)/sum(table(bulk_modules[idx,]$module))*100
cat(subset(bulk_modules[idx,],module==4)$gene_name,sep='\n')

library('pracma')
mat_group <- data.frame(repmat(0,29,6))
rownames(mat_group) <- as.character(c(1:29))
colnames(mat_group) <- c('UU','DD','FU','FD','UD','DU')
mat_group[match(rownames(p1),rownames(mat_group)),1] = p1
mat_group[match(rownames(p2),rownames(mat_group)),2] = p2
mat_group[match(rownames(p3),rownames(mat_group)),3] = p3
mat_group[match(rownames(p4),rownames(mat_group)),4] = p4
mat_group[match(rownames(p5),rownames(mat_group)),5] = p5
mat_group[match(rownames(p6),rownames(mat_group)),6] = p6

percent_df <- data.frame(type=rep(colnames(mat_group),each=29),
	Freq=unlist(array(mat_group)),Var1=paste0('M',rep(rownames(mat_group),6)))


percent_df$sum <- unlist(100-cumsum(mat_group)) + array(diff(as.matrix(rbind(rep(0,6),cumsum(mat_group))))/2)
percent_df$label <- percent_df$Freq
percent_df$label[percent_df$label<2.5] = NA
percent_df$Var1 <- factor(percent_df$Var1,levels=paste0('M',c(1:29)))
percent_df$color <- mapping[c(1:29)]
percent_df$color <- factor(percent_df$color,levels=mapping[c(1:29)])
percent_df$label[!is.na(percent_df$label)] = paste0(paste(percent_df$Var1[!is.na(percent_df$label)],round(percent_df$label[!is.na(percent_df$label)],0),sep=': '),'%')

pdf('~/Downloads/hdWGCNA_TOM/Bulk_gene_program.pdf',width=9,height=6)
ggplot(percent_df, aes(fill=color, y=Freq, x=type,label=label)) +  
	geom_bar(position="stack", stat="identity",width=0.6) + 
    scale_fill_manual(values = mapping[c(1:29)]) + 
	theme_classic() + xlab("Gene Program") + ylab("Frequency") + 
	labs(fill="Module",color='black') + 
	theme(text = element_text(size=20),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),legend.text=element_text(color="black")) + 
	geom_label_repel(aes(type,sum,label=label),fill='white',nudge_x=0,direction="y")
dev.off()

#######################################
############  FIGURE S1F  #############
#######################################
idx<-match(NF_2_pRV_down_2_RVF_down,bulk_modules$gene_name)
idx<-idx[!is.na(idx)]
p2 <- table(bulk_modules[idx,]$module)/sum(table(bulk_modules[idx,]$module))*100
#table(bulk_modules$module)/sum(table(bulk_modules$module)) * 100
genes_int <- subset(bulk_modules[idx,],module==1)$gene_name 

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

library(viridis)
enriched <- enrichr(genes_int, dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('~/Downloads/hdWGCNA_TOM/Bulk_down_down_enrichr.pdf',width=6,height=2.5)
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
############  FIGURE S1G  #############
#######################################
a<-nf.vs.prv[genes_int,]
translation <- strsplit(enriched[[4]]$Genes[1],";")[[1]]
keyvals <- ifelse(rownames(a) %in% translation,'blue','red')
names(keyvals)[keyvals == 'blue'] <- 'high'
names(keyvals)[keyvals == 'red'] <- 'low'

pdf('~/Downloads/hdWGCNA_TOM/MitoRibo_down_down_bulk_enrichr_nf_prv.pdf',width=5,height=6)
EnhancedVolcano(a,lab = rownames(a),x = 'log2FoldChange',y = 'padj',
	pCutoff=0.1,FCcutoff=0.25,title = "", borderColour = 'black',
	xlim=c(0,2),colCustom = keyvals,
	selectLab = rownames(a)[which(names(keyvals) %in% c('high'))])
dev.off()


a<-nf.vs.rvf[genes_int,]
translation <- strsplit(enriched[[4]]$Genes[1],";")[[1]]
keyvals <- ifelse(rownames(a) %in% translation,'blue','red')
names(keyvals)[keyvals == 'blue'] <- 'high'
names(keyvals)[keyvals == 'red'] <- 'low'

pdf('~/Downloads/hdWGCNA_TOM/MitoRibo_down_down_bulk_enrichr_nf_rvf.pdf',width=5,height=6)
EnhancedVolcano(a,lab = rownames(a),x = 'log2FoldChange',y = 'padj',
	pCutoff=0.1,FCcutoff=0.25,title = "", borderColour = 'black',
	xlim=c(0,2),colCustom = keyvals,
	selectLab = rownames(a)[which(names(keyvals) %in% c('high'))])
dev.off()


