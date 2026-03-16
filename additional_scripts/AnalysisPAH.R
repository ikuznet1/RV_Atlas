library(reticulate)
library(ggfortify)
library(edgeR)
library(RColorBrewer)
library(EnhancedVolcano)
library(DESeq2)
library(tximport)
library(biomaRt)
library("sva")
library(stringr)


#######################################
#############  PAH DATA  ##############
#######################################

pah_raw <- read.csv('./dependencies/shared/GSE240921_processed-data-human.csv',sep=',')
gene.annot <- read.csv('./dependencies/shared/Human.GRCh38.p13.annot.tsv',sep='\t')
pah_raw$names <- gene.annot$Symbol[match(pah_raw$id,gene.annot$EnsemblGeneID)]
pah_raw$names[is.na(pah_raw$names)] <- pah_raw$id[is.na(pah_raw$names)]
subj.names <- colnames(pah_raw)[2:41]
subj.group <- str_split_i(subj.names,'_',1)
subj.group[subj.group == 'RV.Normal'] <- 'Control'
subj.group[subj.group == 'RV.Compen'] <- 'Compensated'
subj.group[subj.group == 'RV.Failing'] <- 'Decompensated'

#pah_proc <- pah_raw[,2:41]
pah_proc <- subset(pah_raw,names %in% shared_genes)[,2:41]

#rownames(pah_proc) <- pah_raw$names
rownames(pah_proc) <- subset(pah_raw,names %in% shared_genes)$names

coldata <- data.frame(group = subj.group)
rownames(coldata) <- subj.names

dds <- DESeqDataSetFromMatrix(countData = pah_proc,
                              colData = coldata,
                              design = ~ group)


dds <- estimateSizeFactors(dds)
idx <- rowSums( counts(dds, normalized=TRUE) >= 5 ) >= 3
dds <- dds[idx,]

normalized_counts <- counts(dds, normalized=TRUE)

mod  <- model.matrix(~group, colData(dds))
mod0 <- model.matrix(~1, colData(dds))
svseq <- svaseq(normalized_counts, mod, mod0)


coldata$SV1 <- svseq$sv[,1]
coldata$SV2 <- svseq$sv[,2]
coldata$SV3 <- svseq$sv[,3]
coldata$SV4 <- svseq$sv[,4]
coldata$SV5 <- svseq$sv[,5]

dds <- DESeqDataSetFromMatrix(countData = pah_proc,
                              colData = coldata,
                              design = ~ group+SV1+SV2+SV3+SV4+SV5)



dds <- estimateSizeFactors(dds)
idx <- rowSums( counts(dds, normalized=TRUE) >= 5 ) >= 3
dds <- dds[idx,]

normalized_counts <- counts(dds, normalized=TRUE)


dds <- DESeq(dds)

vst <- vst(dds,blind = FALSE)

mat <- assay(vst)
mm <- model.matrix(~group, colData(vst))
mat <- limma::removeBatchEffect(mat, covariates=colData(vst)[,2:6], design=mm)
assay(vst) <- mat

pah.nf.vs.prv <- lfcShrink(dds,contrast=c('group','Control','Compensated'), type="ashr")
pah.nf.vs.rvf <- lfcShrink(dds,contrast=c('group','Control','Decompensated'), type="ashr")
pah.prv.vs.rvf <- lfcShrink(dds,contrast=c('group','Compensated','Decompensated'), type="ashr")


write.csv(pah.nf.vs.prv,'./outputs/pah.nf.vs.prv.csv')
write.csv(pah.nf.vs.rvf,'./outputs/pah.nf.vs.rvf.csv')
write.csv(pah.prv.vs.rvf,'./outputs/pah.prv.vs.rvf.csv')


pdf('./outputs/ph.bulk_RNAseq_pca.pdf',width=10,height=3)
plotPCA(vst,intgroup=c("group"),ntop=27047) + theme_classic() + theme(axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.title.x=element_text(size=24),axis.title.y=element_text(size=24),legend.title=element_text(size=24),legend.text=element_text(size=24)) + labs(color='Disease',shape='Sex') + geom_point(size=3.5)
dev.off()



pdf('./outputs/ph.bulk_pRV_vs_RVF_volcano.pdf',width=8,height=10)
EnhancedVolcano(pah.prv.vs.rvf,lab = rownames(pah.prv.vs.rvf),x = 'log2FoldChange',y = 'padj',pCutoff=0.1,FCcutoff=0.25,title = "", borderColour = 'black') +  ggplot2::coord_cartesian(xlim=c(-1, 1)) 
dev.off()

pdf('./outputs/ph.bulk_NF_vs_RVF_volcano.pdf',width=8,height=10)
EnhancedVolcano(pah.nf.vs.rvf,lab = rownames(pah.nf.vs.rvf),x = 'log2FoldChange',y = 'padj',pCutoff=0.1,FCcutoff=0.25,title = "", borderColour = 'black') +  ggplot2::coord_cartesian(xlim=c(-6, 6)) 
dev.off()

pdf('./outputs/ph.bulk_NF_vs_pRV_volcano.pdf',width=8,height=10)
EnhancedVolcano(pah.nf.vs.prv,lab = rownames(pah.nf.vs.prv),x = 'log2FoldChange',y = 'padj',pCutoff=0.1,FCcutoff=0.25,title = "", borderColour = 'black') +  ggplot2::coord_cartesian(xlim=c(-6, 6)) 
dev.off()



######################################
#############  RV DATA  ##############
######################################

bulk <- read.csv('./dependencies/shared/BulkRNA/counts.csv')
meta <- read.csv('./dependencies/shared/BulkRNA/metadata.csv')
toDel <- seq(1,dim(meta)[1],2)
meta <- meta[-toDel,]

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

res <- getBM(attributes = c('ensembl_transcript_id_version','ensembl_gene_id','external_transcript_name',                           
'external_gene_name'),              
filters = 'ensembl_transcript_id_version',               
values = bulk[,1],              
mart = mart)
tx2gene <- res[,c(1,4)]


path = './dependencies/shared/BulkRNA/30-238740824'
files1 <- list.files(path,pattern = "\\.h5$",recursive=TRUE)
files1<-paste0(path,'/',files1)
path = './dependencies/shared/BulkRNA/30-196345105'
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
#mod  <- model.matrix(~category, colData(ddsSE))
#mod0 <- model.matrix(~1, colData(ddsSE))
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

write.csv(nf.vs.prv,'./output/rv.nf.vs.prv.csv')
write.csv(nf.vs.rvf,'./output/rv.nf.vs.rvf.csv')
write.csv(prv.vs.rvf,'./output/rv.prv.vs.rvf.csv')


######################################
#############  COMBINED ##############
######################################

raw_counts_RV <- counts(ddsSE, normalized=FALSE)
raw_counts_PAH <- counts(dds, normalized=FALSE)
shared_genes <- intersect(rownames(raw_counts_RV),rownames(raw_counts_PAH))
raw_counts_RV <- raw_counts_RV[rownames(raw_counts_RV) %in% shared_genes,]
raw_counts_PAH <- raw_counts_PAH[rownames(raw_counts_PAH) %in% shared_genes,]

raw_counts_merge = cbind(raw_counts_RV,raw_counts_PAH)
rv_grps <- colData(ddsSE)$category
pah_grps <- colData(dds)$group
pah_grps <- as.character(pah_grps)
pah_grps[pah_grps == 'Control'] = 'NF'
pah_grps[pah_grps == 'Compensated'] = 'pRV'
pah_grps[pah_grps == 'Decompensated'] = 'RVF'
pah_grps <- factor(pah_grps,level=c('NF','pRV','RVF'))
merge_grps <- c(rv_grps,pah_grps)


coldata_merge <- data.frame(group = merge_grps)
rownames(coldata_merge) <- c(as.character(ddsSE$subject),subj.names)
coldata_merge$batch <- factor(c(rep(1,142),rep(2,40)), level = c('1','2'))
coldata_merge$class <- paste0(coldata_merge$group,'_',coldata_merge$batch)


colnames(raw_counts_merge) <- rownames(coldata_merge)

dds_merge <- DESeqDataSetFromMatrix(countData = raw_counts_merge,
                              colData = coldata_merge,
                              design = ~ group + batch)

dds_merge <- estimateSizeFactors(dds_merge)
idx <- rowSums( counts(dds_merge, normalized=TRUE) >= 5 ) >= 3
dds_merge <- dds_merge[idx,]


normalized_counts <- counts(dds_merge, normalized=TRUE)

mod  <- model.matrix(~group + batch, colData(dds_merge))
mod0 <- model.matrix(~1, colData(dds_merge))
svseq <- svaseq(normalized_counts, mod, mod0)


coldata_merge$SV1 <- svseq$sv[,1]
coldata_merge$SV2 <- svseq$sv[,2]
coldata_merge$SV3 <- svseq$sv[,3]
coldata_merge$SV4 <- svseq$sv[,4]
coldata_merge$SV5 <- svseq$sv[,5]
coldata_merge$SV6 <- svseq$sv[,6]
coldata_merge$SV7 <- svseq$sv[,7]
coldata_merge$SV8 <- svseq$sv[,8]
coldata_merge$SV9 <- svseq$sv[,9]
coldata_merge$SV10 <- svseq$sv[,10]
coldata_merge$SV11 <- svseq$sv[,11]
coldata_merge$SV12 <- svseq$sv[,12]
coldata_merge$SV13 <- svseq$sv[,13]
coldata_merge$SV14 <- svseq$sv[,14]
coldata_merge$SV15 <- svseq$sv[,15]
coldata_merge$SV16 <- svseq$sv[,16]
coldata_merge$SV17 <- svseq$sv[,17]
coldata_merge$SV18 <- svseq$sv[,18]
coldata_merge$SV19 <- svseq$sv[,19]
coldata_merge$SV20 <- svseq$sv[,20]
coldata_merge$SV21 <- svseq$sv[,21]
coldata_merge$SV22 <- svseq$sv[,22]
coldata_merge$SV23 <- svseq$sv[,23]
coldata_merge$SV24 <- svseq$sv[,24]
coldata_merge$SV25 <- svseq$sv[,25]
coldata_merge$SV26 <- svseq$sv[,26]


dds_merge <- DESeqDataSetFromMatrix(countData = raw_counts_merge,
                              colData = coldata_merge,
                              design = ~ group+ batch+SV1+SV2+SV3+SV4+SV5+SV6+SV7+SV8+SV9+SV10+SV11+SV12+SV13+SV14+
                              SV15+SV16+SV17+SV18+SV19+SV20+SV21+SV22+SV23+SV24+SV25+SV26)



dds_merge <- estimateSizeFactors(dds_merge)
idx <- rowSums( counts(dds_merge, normalized=TRUE) >= 5 ) >= 3
dds_merge <- dds_merge[idx,]

normalized_counts <- counts(dds_merge, normalized=TRUE)


dds_merge <- DESeq(dds_merge)

vst_merge <- vst(dds_merge,blind = FALSE)

mat <- assay(vst_merge)
mm <- model.matrix(~group, colData(vst_merge))
mat <- limma::removeBatchEffect(mat, covariates=colData(vst_merge)[,c(2,4:29)], design=mm)
assay(vst_merge) <- mat

merge.nf.vs.prv <- lfcShrink(dds_merge,contrast=c('group','NF','pRV'), type="ashr")
merge.nf.vs.rvf <- lfcShrink(dds_merge,contrast=c('group','NF','RVF'), type="ashr")
merge.prv.vs.rvf <- lfcShrink(dds_merge,contrast=c('group','pRV','RVF'), type="ashr")

write.csv(merge.nf.vs.prv,'./output/merge.nf.vs.prv.csv')
write.csv(merge.nf.vs.rvf,'./output/merge.nf.vs.rvf.csv')
write.csv(merge.prv.vs.rvf,'./output/merge.prv.vs.rvf.csv')


plotPCA(vst_merge,intgroup=c("group","batch"),ntop=1000) + theme_classic() + theme(axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.title.x=element_text(size=24),axis.title.y=element_text(size=24),legend.title=element_text(size=24),legend.text=element_text(size=24)) + labs(color='Disease',shape='Sex') + geom_point(size=3.5)


#######################################
########### Venn Diagram  #############
#######################################

pah.nf.vs.prv <- read.csv('./dependencies/shared/pah.nf.vs.prv.csv',row.names = 1)
pah.nf.vs.rvf <- read.csv('./dependencies/shared//pah.nf.vs.rvf.csv',row.names = 1)
pah.prv.vs.rvf <- read.csv('./dependencies/shared/pah.prv.vs.rvf.csv',row.names = 1)

nf.vs.prv <- read.csv('./dependencies/shared/rv.nf.vs.prv.csv',row.names = 1)
nf.vs.rvf <- read.csv('./dependencies/shared/rv.nf.vs.rvf.csv',row.names = 1)
prv.vs.rvf <- read.csv('~./dependencies/shared/rv.prv.vs.rvf.csv',row.names = 1)


library(VennDiagram)


shared_genes_postproc <- intersect(rownames(pah.nf.vs.rvf),rownames(nf.vs.rvf))

pah.nf.vs.rvf.signif.up <- intersect(rownames(subset(pah.nf.vs.rvf, padj < 0.1 & log2FoldChange > 0.1)),shared_genes_postproc)
pah.nf.vs.rvf.signif.down <- intersect(rownames(subset(pah.nf.vs.rvf, padj < 0.1 & log2FoldChange < -0.1)),shared_genes_postproc)

rv.nf.vs.rvf.signif.up <- intersect(rownames(subset(nf.vs.rvf, padj < 0.1 & log2FoldChange > 0.1)),shared_genes_postproc)
rv.nf.vs.rvf.signif.down <- intersect(rownames(subset(nf.vs.rvf, padj < 0.1 & log2FoldChange < -0.1)),shared_genes_postproc)


nf.vs.rvf.pah.up.rv.up <- intersect(pah.nf.vs.rvf.signif.up,rv.nf.vs.rvf.signif.up)
nf.vs.rvf.pah.up.rv.down <- intersect(pah.nf.vs.rvf.signif.up,rv.nf.vs.rvf.signif.down)
nf.vs.rvf.pah.down.rv.up <- intersect(pah.nf.vs.rvf.signif.down,rv.nf.vs.rvf.signif.up)
nf.vs.rvf.pah.down.rv.down <- intersect(pah.nf.vs.rvf.signif.down,rv.nf.vs.rvf.signif.down)

venn.diagram(
  x = list(
    pah.nf.vs.rvf.signif.up,
    pah.nf.vs.rvf.signif.down,
    rv.nf.vs.rvf.signif.up,
    rv.nf.vs.rvf.signif.down
  ),
  category.names = c("PAH UP","PAH DOWN","RV UP","RV DOWN"),
  filename = './output/PAH_vs_RVF.png',
  output = TRUE
)


df.nf.vs.rvf <- data.frame(RV = nf.vs.rvf$log2FoldChange[match(shared_genes_postproc,rownames(nf.vs.rvf))[!is.na(match(shared_genes_postproc,rownames(nf.vs.rvf)))]],
  PAH = pah.nf.vs.rvf$log2FoldChange[match(shared_genes_postproc,rownames(pah.nf.vs.rvf))[!is.na(match(shared_genes_postproc,rownames(pah.nf.vs.rvf)))]],
  gene = rownames(nf.vs.rvf)[match(shared_genes_postproc,rownames(nf.vs.rvf))[!is.na(match(shared_genes_postproc,rownames(nf.vs.rvf)))]]
)


ggplot(df.nf.vs.rvf,aes(x = RV,y = PAH, label = gene)) + geom_point() + xlim(-2,2) + ylim(-1,1)


df.nf.vs.prv <- data.frame(RV = nf.vs.prv$log2FoldChange[match(shared_genes_postproc,rownames(nf.vs.prv))[!is.na(match(shared_genes_postproc,rownames(nf.vs.prv)))]],
  PAH = pah.nf.vs.prv$log2FoldChange[match(shared_genes_postproc,rownames(pah.nf.vs.prv))[!is.na(match(shared_genes_postproc,rownames(pah.nf.vs.prv)))]],
  gene = rownames(nf.vs.prv)[match(shared_genes_postproc,rownames(nf.vs.prv))[!is.na(match(shared_genes_postproc,rownames(nf.vs.prv)))]]
)


ggplot(df.nf.vs.prv,aes(x = RV,y = PAH, label = gene)) + geom_point() + xlim(-5,5) + ylim(-5,5)



df.prv.vs.rvf <- data.frame(RV = prv.vs.rvf$log2FoldChange[match(shared_genes_postproc,rownames(prv.vs.rvf))[!is.na(match(shared_genes_postproc,rownames(prv.vs.rvf)))]],
  PAH = pah.prv.vs.rvf$log2FoldChange[match(shared_genes_postproc,rownames(pah.prv.vs.rvf))[!is.na(match(shared_genes_postproc,rownames(pah.prv.vs.rvf)))]],
  gene = rownames(prv.vs.rvf)[match(shared_genes_postproc,rownames(prv.vs.rvf))[!is.na(match(shared_genes_postproc,rownames(prv.vs.rvf)))]]
)


ggplot(df.prv.vs.rvf,aes(x = RV,y = PAH, label = gene)) + geom_point() + geom_label_repel() + xlim(-5,5) + ylim(-5,5) 


dbs <- c("ChEA_2022","WikiPathways_2024_Human","Reactome_Pathways_2024","GO_Biological_Process_2025")

library(enrichR)
library(forcats)
library(viridis)

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



#######################################
########### Gene ontology  #############
#######################################

enriched <- enrichr(nf.vs.rvf.pah.up.rv.up, dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('./output/nf.vs.rvf.pah.up.rv.up.pdf',width=5,height=2.5)
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

enriched <- enrichr(nf.vs.rvf.pah.up.rv.down, dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('./output/nf.vs.rvf.pah.up.rv.down.pdf',width=5,height=2.5)
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


enriched <- enrichr(nf.vs.rvf.pah.down.rv.up, dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('./output/nf.vs.rvf.pah.down.rv.up.pdf',width=5,height=2.5)
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


enriched <- enrichr(nf.vs.rvf.pah.down.rv.down, dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('./output/nf.vs.rvf.pah.down.rv.down.pdf',width=5,height=2.5)
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

### Reactome

enriched <- enrichr(nf.vs.rvf.pah.up.rv.up, dbs)
enriched[[4]] <- subset(enriched[[3]],Adjusted.P.value<0.05)
pdf('./output/reactome_nf.vs.rvf.pah.up.rv.up.pdf',width=5,height=2.5)
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

enriched <- enrichr(nf.vs.rvf.pah.up.rv.down, dbs)
enriched[[4]] <- subset(enriched[[3]],Adjusted.P.value<0.05)
pdf('./output/reactome_nf.vs.rvf.pah.up.rv.down.pdf',width=5,height=2.5)
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


enriched <- enrichr(nf.vs.rvf.pah.down.rv.up, dbs)
enriched[[4]] <- subset(enriched[[3]],Adjusted.P.value<0.05)
pdf('./output/reactome_nf.vs.rvf.pah.down.rv.up.pdf',width=5,height=2.5)
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


enriched <- enrichr(nf.vs.rvf.pah.down.rv.down, dbs)
enriched[[4]] <- subset(enriched[[3]],Adjusted.P.value<0.05)
pdf('./output/reactome_nf.vs.rvf.pah.down.rv.down.pdf',width=5,height=2.5)
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


########################
###Analysis of pRV

library(VennDiagram)


pah.nf.vs.prv <- read.csv('./dependencies/shared/pah.nf.vs.prv.csv',row.names = 1)
pah.nf.vs.rvf <- read.csv('./dependencies/shared/pah.nf.vs.rvf.csv',row.names = 1)
pah.prv.vs.rvf <- read.csv('./dependencies/shared/pah.prv.vs.rvf.csv',row.names = 1)

nf.vs.prv <- read.csv('./dependencies/shared/rv.nf.vs.prv.csv',row.names = 1)
nf.vs.rvf <- read.csv('./dependencies/shared/rv.nf.vs.rvf.csv',row.names = 1)
prv.vs.rvf <- read.csv('./dependencies/shared/rv.prv.vs.rvf.csv',row.names = 1)

shared_genes_postproc <- intersect(rownames(pah.nf.vs.rvf),rownames(nf.vs.rvf))



pah.nf.vs.rvf.signif.up <- intersect(rownames(subset(pah.nf.vs.rvf, padj < 0.1 & log2FoldChange > 0.1)),shared_genes_postproc)
pah.nf.vs.rvf.signif.down <- intersect(rownames(subset(pah.nf.vs.rvf, padj < 0.1 & log2FoldChange < -0.1)),shared_genes_postproc)

pah.prv.vs.rvf.signif.up <- intersect(rownames(subset(pah.prv.vs.rvf, padj < 0.1 & log2FoldChange > 0.1)),shared_genes_postproc)
pah.prv.vs.rvf.signif.down <- intersect(rownames(subset(pah.prv.vs.rvf, padj < 0.1 & log2FoldChange < -0.1)),shared_genes_postproc)


rv.prv.vs.rvf.signif.up <- intersect(rownames(subset(prv.vs.rvf, padj < 0.1 & log2FoldChange > 0.1)),shared_genes_postproc)
rv.prv.vs.rvf.signif.down <- intersect(rownames(subset(prv.vs.rvf, padj < 0.1 & log2FoldChange < -0.1)),shared_genes_postproc)

pressure_related_up <- intersect(rv.prv.vs.rvf.signif.down,pah.prv.vs.rvf.signif.down)
pressure_related_down <- intersect(rv.prv.vs.rvf.signif.up,pah.prv.vs.rvf.signif.up)


not_pressure_related_down <- intersect(rv.prv.vs.rvf.signif.down,pah.prv.vs.rvf.signif.up)
not_pressure_related_up <- intersect(rv.prv.vs.rvf.signif.up,pah.prv.vs.rvf.signif.down)


venn.diagram(
  x = list(
    pah.prv.vs.rvf.signif.up,
    pah.prv.vs.rvf.signif.down,
    rv.prv.vs.rvf.signif.up,
    rv.prv.vs.rvf.signif.down
  ),
  category.names = c("PAH UP","PAH DOWN","RV UP","RV DOWN"),
  filename = './output/PAH_vs_RVF_pRV.png',
  output = TRUE
)


RV_specific_up <- setdiff(rv.prv.vs.rvf.signif.up,pah.prv.vs.rvf.signif.up)
RV_specific_down <- setdiff(rv.prv.vs.rvf.signif.down,pah.prv.vs.rvf.signif.down)



dbs <- c("ChEA_2022","WikiPathways_2024_Human","Reactome_Pathways_2024","GO_Biological_Process_2025")

library(enrichR)
library(forcats)
library(viridis)

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



enriched <- enrichr(pressure_related_up, dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('./output/prv.vs.rvf.pah.up.rv.up.pdf',width=5,height=2.5)
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


enriched <- enrichr(pressure_related_down, dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('./output/prv.vs.rvf.pah.down.rv.down.pdf',width=5,height=2.5)
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


enriched <- enrichr(RV_specific_up, dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('./output/prv.vs.rvf.rv_spec.up.pdf',width=5,height=2.5)
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


enriched <- enrichr(RV_specific_down, dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('./output//prv.vs.rvf.rv_spec.down.pdf',width=5,height=2.5)
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






