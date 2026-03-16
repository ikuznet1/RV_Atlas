library(BayesPrism)
library(biomaRt)
library(tximport)
library(DESeq2)
library(Seurat)

####### LOAD SINGLE CELL DATA ######

M1<-readRDS('./dependencies/shared/Post_R3_FINAL_with_counts.rds')
sc.dat <- t(GetAssayData(object = M1, assay = "decontXcounts", slot = "counts"))
rownames(sc.dat) <- colnames(M1)
cell.type.labels <- M1$Names

rm(M1)
sc.dat <- as.matrix(sc.dat)


####### LOAD BULK DATA ######


bulk <- read.csv('./dependencies/shared/BulkRNA/counts.csv')
meta <- read.csv('./dependencies/shared/BulkRNA/metadata.csv')
toDel <- seq(1,dim(meta)[1],2)
meta <- meta[-toDel,]

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

res <- getBM(attributes = c('ensembl_transcript_id_version',
  'ensembl_gene_id','external_transcript_name','external_gene_name'),
   filters = 'ensembl_transcript_id_version',values = bulk[,1],mart = mart)
tx2gene <- res[,c(1,4)]


path = './dependencies/shared/BulkRNA/30-238740824'
files1 <- list.files(path,pattern = "\\.h5$",recursive=TRUE)
files1<-paste0(path,'/',files1)
path = './dependencies/shared/BulkRNA/30-238740824'
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

bk.dat <- t(counts(ddsSE, normalized=FALSE))
rownames(bk.dat) <- bulk.meta$subject
#saveRDS(bk.dat,'./output/bulk_data.rds')


####### Deconvolve RNA-seq ######

sc.dat <- sc.dat[,colnames(sc.dat) %in% intersect(colnames(sc.dat),colnames(bk.dat))]


plot.cor.phi(input=sc.dat,input.labels=cell.type.labels,title="cell state correlation",cexRow=0.2, cexCol= 0.2, margins=c(2,2))

sc.stat <- plot.scRNA.outlier(
	input=sc.dat, #make sure the colnames are gene symbol or ENSMEBL ID
	cell.type.labels=cell.type.labels,
	species="hs", #currently only human(hs) and mouse(mm) annotations are supported
	return.raw=TRUE #return the data used for plotting.
)

bk.stat <- plot.bulk.outlier(
	bulk.input=bk.dat,#make sure the colnames are gene symbol or ENSMEBL ID
	sc.input=sc.dat, #make sure the colnames are gene symbol or ENSMEBL ID
	cell.type.labels=cell.type.labels,
	species="hs", #currently only human(hs) and mouse(mm) annotations are supported
	return.raw=TRUE
)


sc.dat.filtered <- cleanup.genes(input=sc.dat,
	input.type="count.matrix",
	species="hs",
	gene.group=c( "Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY") ,
	exp.cells=5)

#saveRDS(sc.dat,'~/Downloads/hdWGCNA_TOM/PAH/sc.dat.rds')
rm(sc.dat)

#note this function only works for human data. For other species, you are advised to make plots by yourself.
plot.bulk.vs.sc(sc.input = sc.dat.filtered,
	bulk.input = bk.dat
)

sc.dat.filtered.sc <- select.gene.type(sc.dat.filtered,gene.type = "protein_coding")
rm(sc.dat.filtered)

# Select marker genes (Optional)
# performing pair-wise t test for cell states from different cell types
#diff.exp.stat <- get.exp.stat(sc.dat=sc.dat.filtered.sc,# filter genes to reduce memory use
#	cell.type.labels=cell.type.labels,
#	cell.state.labels=cell.type.labels,
#	pseudo.count=0.1, #a numeric value used for log2 transformation. =0.1 for 10x data, =10 for smart-cell.count.cutoff=50, # a numeric value to exclude cell state with number of cells fewer than n.cores=1 #number of threads
#)

#sc.dat.filtered.pc.sig <- select.marker(sc.dat=sc.dat.filtered.sc,stat=diff.exp.stat,pval.max=0.01,lfc.min=0.1)
#rm(sc.dat.filtered.sc)

myPrism <- new.prism(
	reference=sc.dat.filtered.sc,
	mixture=bk.dat,
	input.type="count.matrix",
	cell.type.labels = cell.type.labels,
	cell.state.labels = cell.type.labels,
	key=NULL,
	outlier.cut=0.01,
	outlier.fraction=0.1,
)

bp.res <- run.prism(prism = myPrism, n.cores=10)

#saveRDS(bp.res,'./output/RV_RVA_deconv.rds')

# extract posterior mean of cell type fraction theta
theta <- get.fraction(bp=bp.res,
	which.theta="final",
	state.or.type="type"
)

category <- meta$maineffect[match(rownames(theta),meta$title)]

df <- data.frame(category = category, CM = theta[,'CM'], FB = theta[,'FB'],
	 Myeloid = theta[,'Myeloid'], EC = theta[,'EC'], Adipo = theta[,'Adipo'],
	  PC = theta[,'PC'], SM = theta[,'SM'], Epi = theta[,'Epi'])

df.2 <- data.frame(type = gsub('[0-9]+', '', names(unlist(df[,2:9]))), 
	freq = unlist(df[,2:9]),
	group = factor(rep(df$category,8),level=c('NF','pRV','RVF')))

library(ggpubr)

p <- ggboxplot(df.2,x="group",y="freq",fill="group",group="group")+
	theme_classic() + 
	theme(axis.text.x=element_text(size=16),
	axis.text.y=element_text(size=16),
	axis.title.x=element_text(size=16),
	axis.title.y=element_text(size=16),
	legend.title=element_text(size=16),
	legend.text=element_text(size=16),
	text=element_text(color='black'),
	axis.text=element_text(color='black')) + 
	labs(color='Group',x="Disease",y='Frequency') + 
	facet_wrap(~type,ncol=8) + 
	#stat_compare_means(aes(group=group),comparisons=my_comparisons,method="t.test",ref.group="NF")+
	stat_compare_means(aes(group=group),method="anova")

pdf('./output/deconv_cell_freq.pdf',width=12,height=5)
p
dev.off()

# extract posterior mean of cell type-specific gene expression count matrix Z
Z.CM <- get.exp(bp=bp.res,
state.or.type="type",
cell.name="CM")

Z.Myeloid <- get.exp(bp=bp.res,
state.or.type="type",
cell.name="Myeloid")

Z.EC <- get.exp(bp=bp.res,
state.or.type="type",
cell.name="EC")

Z.FB <- get.exp(bp=bp.res,
state.or.type="type",
cell.name="FB")



###############
### Myeloid ###
###############

g <- rownames(subset(diff.exp.stat[[4]],min.lfc > 0.1 & pval.up.min < 0.01))

ddsSE <- DESeqDataSetFromMatrix(countData = round(t(Z.Myeloid[,g])),
                              colData = bulk.meta,
                              design = ~ category)


ddsSE <- estimateSizeFactors(ddsSE)
idx <- rowSums( counts(ddsSE, normalized=TRUE) >= 5 ) >= 3
ddsSE <- ddsSE[idx,]

normalized_counts <- counts(ddsSE, normalized=TRUE)


mod  <- model.matrix(~category, colData(ddsSE))
mod0 <- model.matrix(~1, colData(ddsSE))
svseq <- svaseq(normalized_counts, mod, mod0)

svs <- svseq$sv[,1:svseq$n.sv]
colnames(svs) <- paste0('SV',1:svseq$n.sv)

bulk.meta.temp <- cbind(bulk.meta,svs)



ddsSE <- DESeqDataSetFromMatrix(countData = round(t(Z.Myeloid[,g])),
                              colData = bulk.meta.temp,
                              design = as.formula(paste0('~category',paste0('+','SV',1:svseq$n.sv,collapse=''))))

ddsSE <- estimateSizeFactors(ddsSE)
idx <- rowSums( counts(ddsSE, normalized=TRUE) >= 5 ) >= 3
ddsSE <- ddsSE[idx,]

normalized_counts <- counts(ddsSE, normalized=TRUE)



ddsSE <- DESeq(ddsSE)

vstSE <- vst(ddsSE,blind = FALSE)

mat <- assay(vstSE)
mm <- model.matrix(~category, colData(vstSE))
mat <- limma::removeBatchEffect(mat, covariates=colData(vstSE)[,27:(26+svseq$n.sv)], design=mm)
assay(vstSE) <- mat

plotPCA(vstSE,intgroup=c("category"),ntop=1000) + theme_classic() + theme(axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.title.x=element_text(size=24),axis.title.y=element_text(size=24),legend.title=element_text(size=24),legend.text=element_text(size=24)) + labs(color='Disease',shape='Sex') + geom_point(size=3.5)


nf.vs.prv <- lfcShrink(ddsSE,contrast=c('category','NF','pRV'), type="ashr")
nf.vs.rvf <- lfcShrink(ddsSE,contrast=c('category','NF','RVF'), type="ashr")
prv.vs.rvf <- lfcShrink(ddsSE,contrast=c('category','pRV','RVF'), type="ashr")

cat(rownames(subset(prv.vs.rvf,padj < 0.05 & log2FoldChange < -1)),sep='\n')


nf.vs.rvf.down <- rownames(subset(nf.vs.rvf,padj < 0.05 & log2FoldChange < -.25))
nf.vs.rvf.up <- rownames(subset(nf.vs.rvf,padj < 0.05 & log2FoldChange > .25))


dbs <- c("ChEA_2022","WikiPathway_2023_Human","Reactome_2024","GO_Biological_Process_2025")

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


enriched <- enrichr(nf.vs.rvf.down, dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('./output/deconvolve_Myeloid_nf.vs.rvf.down.pdf',width=5,height=2.5)
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

enriched <- enrichr(nf.vs.rvf.up, dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('./output/deconvolve_Myeloid_nf.vs.rvf.up.pdf',width=5,height=2.5)
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




prv.vs.rvf.down <- rownames(subset(prv.vs.rvf,padj < 0.05 & log2FoldChange < 0))
prv.vs.rvf.up <- rownames(subset(prv.vs.rvf,padj < 0.05 & log2FoldChange > 0))


dbs <- c("ChEA_2022","WikiPathway_2023_Human","Reactome_2024","GO_Biological_Process_2025")

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


enriched <- enrichr(prv.vs.rvf.down, dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('./output/deconvolve_Myeloid_prv.vs.rvf.down.pdf',width=5,height=2.5)
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

enriched <- enrichr(prv.vs.rvf.up, dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('./output/deconvolve_Myeloid_prv.vs.rvf.up',width=5,height=2.5)
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






###############
### CM ###
###############

#g <- rownames(subset(diff.exp.stat[[1]],min.lfc > 0.1 & pval.up.min < 0.01))

ddsSE <- DESeqDataSetFromMatrix(countData = round(t(Z.CM)),
                              colData = bulk.meta,
                              design = ~ category)


ddsSE <- estimateSizeFactors(ddsSE)
idx <- rowSums( counts(ddsSE, normalized=TRUE) >= 5 ) >= 3
ddsSE <- ddsSE[idx,]

normalized_counts <- counts(ddsSE, normalized=TRUE)


mod  <- model.matrix(~category, colData(ddsSE))
mod0 <- model.matrix(~1, colData(ddsSE))
svseq <- svaseq(normalized_counts, mod, mod0)

svs <- svseq$sv[,1:svseq$n.sv]
colnames(svs) <- paste0('SV',1:svseq$n.sv)

bulk.meta.temp <- cbind(bulk.meta,svs)



ddsSE <- DESeqDataSetFromMatrix(countData = round(t(Z.CM)),
                              colData = bulk.meta.temp,
                              design = as.formula(paste0('~category',paste0('+','SV',1:svseq$n.sv,collapse=''))))

ddsSE <- estimateSizeFactors(ddsSE)
idx <- rowSums( counts(ddsSE, normalized=TRUE) >= 5 ) >= 3
ddsSE <- ddsSE[idx,]

normalized_counts <- counts(ddsSE, normalized=TRUE)



ddsSE <- DESeq(ddsSE)

vstSE <- vst(ddsSE,blind = FALSE)

mat <- assay(vstSE)
mm <- model.matrix(~category, colData(vstSE))
mat <- limma::removeBatchEffect(mat, covariates=colData(vstSE)[,27:(26+svseq$n.sv)], design=mm)
assay(vstSE) <- mat

plotPCA(vstSE,intgroup=c("category"),ntop=1000) + theme_classic() + theme(axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.title.x=element_text(size=24),axis.title.y=element_text(size=24),legend.title=element_text(size=24),legend.text=element_text(size=24)) + labs(color='Disease',shape='Sex') + geom_point(size=3.5)


nf.vs.prv <- lfcShrink(ddsSE,contrast=c('category','NF','pRV'), type="ashr")
nf.vs.rvf <- lfcShrink(ddsSE,contrast=c('category','NF','RVF'), type="ashr")
prv.vs.rvf <- lfcShrink(ddsSE,contrast=c('category','pRV','RVF'), type="ashr")

cat(rownames(subset(prv.vs.rvf,padj < 0.05 & log2FoldChange < -1)),sep='\n')


nf.vs.rvf.down <- rownames(subset(nf.vs.rvf,padj < 0.05 & log2FoldChange < -.25))
nf.vs.rvf.up <- rownames(subset(nf.vs.rvf,padj < 0.05 & log2FoldChange > .25))


dbs <- c("ChEA_2022","WikiPathway_2023_Human","Reactome_2024","GO_Biological_Process_2025")

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


enriched <- enrichr(nf.vs.rvf.down, dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('./output/deconvolve_CM_nf.vs.rvf.down.pdf',width=5,height=2.5)
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

enriched <- enrichr(nf.vs.rvf.up, dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('./output/deconvolve_CM_nf.vs.rvf.up.pdf',width=5,height=2.5)
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




prv.vs.rvf.down <- rownames(subset(prv.vs.rvf,padj < 0.05 & log2FoldChange < 0))
prv.vs.rvf.up <- rownames(subset(prv.vs.rvf,padj < 0.05 & log2FoldChange > 0))


dbs <- c("ChEA_2022","WikiPathway_2023_Human","Reactome_2024","GO_Biological_Process_2025")

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


enriched <- enrichr(prv.vs.rvf.down, dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('./output/deconvolve_CM_prv.vs.rvf.down.pdf',width=5,height=2.5)
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

enriched <- enrichr(prv.vs.rvf.up, dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('./output/deconvolve_CM_prv.vs.rvf.up.pdf',width=5,height=2.5)
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




###############
### EC ###
###############

g <- rownames(subset(diff.exp.stat[[6]],min.lfc > 0.1 & pval.up.min < 0.01))


ddsSE <- DESeqDataSetFromMatrix(countData = ceiling(t(Z.EC[,g])),
                              colData = bulk.meta,
                              design = ~ category)


ddsSE <- estimateSizeFactors(ddsSE)
idx <- rowSums( counts(ddsSE, normalized=TRUE) >= 5 ) >= 3
ddsSE <- ddsSE[idx,]

normalized_counts <- counts(ddsSE, normalized=TRUE)


mod  <- model.matrix(~category, colData(ddsSE))
mod0 <- model.matrix(~1, colData(ddsSE))
svseq <- svaseq(normalized_counts, mod, mod0)

svs <- svseq$sv[,1:svseq$n.sv]
colnames(svs) <- paste0('SV',1:svseq$n.sv)

bulk.meta.temp <- cbind(bulk.meta,svs)



ddsSE <- DESeqDataSetFromMatrix(countData = ceiling(t(Z.EC[,g])),
                              colData = bulk.meta.temp,
                              design = as.formula(paste0('~category',paste0('+','SV',1:svseq$n.sv,collapse=''))))

ddsSE <- estimateSizeFactors(ddsSE)
idx <- rowSums( counts(ddsSE, normalized=TRUE) >= 5 ) >= 3
ddsSE <- ddsSE[idx,]

normalized_counts <- counts(ddsSE, normalized=TRUE)



ddsSE <- DESeq(ddsSE)

vstSE <- vst(ddsSE,blind = FALSE)

mat <- assay(vstSE)
mm <- model.matrix(~category, colData(vstSE))
mat <- limma::removeBatchEffect(mat, covariates=colData(vstSE)[,27:(26+svseq$n.sv)], design=mm)
assay(vstSE) <- mat

plotPCA(vstSE,intgroup=c("category"),ntop=1000) + theme_classic() + theme(axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.title.x=element_text(size=24),axis.title.y=element_text(size=24),legend.title=element_text(size=24),legend.text=element_text(size=24)) + labs(color='Disease',shape='Sex') + geom_point(size=3.5)


nf.vs.prv <- lfcShrink(ddsSE,contrast=c('category','NF','pRV'), type="ashr")
nf.vs.rvf <- lfcShrink(ddsSE,contrast=c('category','NF','RVF'), type="ashr")
prv.vs.rvf <- lfcShrink(ddsSE,contrast=c('category','pRV','RVF'), type="ashr")

cat(rownames(subset(prv.vs.rvf,padj < 0.05 & log2FoldChange < -1)),sep='\n')


nf.vs.rvf.down <- rownames(subset(nf.vs.rvf,padj < 0.05 & log2FoldChange < -.25))
nf.vs.rvf.up <- rownames(subset(nf.vs.rvf,padj < 0.05 & log2FoldChange > .25))


dbs <- c("ChEA_2022","WikiPathway_2023_Human","Reactome_2024","GO_Biological_Process_2025")

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


enriched <- enrichr(nf.vs.rvf.down, dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('./output/deconvolve_EC_nf.vs.rvf.down.pdf',width=5,height=2.5)
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

enriched <- enrichr(nf.vs.rvf.up, dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('./output/deconvolve_EC_nf.vs.rvf.up.pdf',width=5,height=2.5)
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




prv.vs.rvf.down <- rownames(subset(prv.vs.rvf,padj < 0.05 & log2FoldChange < 0))
prv.vs.rvf.up <- rownames(subset(prv.vs.rvf,padj < 0.05 & log2FoldChange > 0))


dbs <- c("ChEA_2022","WikiPathway_2023_Human","Reactome_2024","GO_Biological_Process_2025")

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


enriched <- enrichr(prv.vs.rvf.down, dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('./output/deconvolve_EC_prv.vs.rvf.down.pdf',width=5,height=2.5)
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

enriched <- enrichr(prv.vs.rvf.up, dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('./output/deconvolve_EC_prv.vs.rvf.up.pdf',width=5,height=2.5)
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



##### BISQUE
library(BisqueRNA)

bulk.eset <- Biobase::ExpressionSet(assayData = t(bk.dat))

sample.ids <- rownames(sc.dat.filtered.sc)
# individual.ids and cell.types should be in the same order as in sample.ids
sc.pheno <- data.frame(check.names=F, check.rows=F,
                       stringsAsFactors=F,
                       row.names=sample.ids,
                       SubjectName=stringr::str_split_i(sample.ids,'_',1),
                       cellType=cell.type.labels)

sc.meta <- data.frame(labelDescription=c("SubjectName",
                                         "cellType"),
                      row.names=c("SubjectName",
                                  "cellType"))
sc.pdata <- new("AnnotatedDataFrame",
                data=sc.pheno,
                varMetadata=sc.meta)
sc.eset <- Biobase::ExpressionSet(assayData=t(sc.dat.filtered.sc),
                                  phenoData=sc.pdata)

saveRDS(sc.dat.filtered.sc,'~/Downloads/hdWGCNA_TOM/PAH/sc_data_filt.rds')


res <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset, markers=NULL, use.overlap=TRUE)

ref.based.estimates <- res$bulk.props
knitr::kable(ref.based.estimates, digits=2)


category <- meta$maineffect[match(colnames(ref.based.estimates),meta$title)]

df <- as.data.frame(t(ref.based.estimates))
df$category = category

df.2 <- data.frame(type = gsub('[0-9]+', '', names(unlist(df[,1:12]))), 
	freq = unlist(df[,1:12]),
	group = factor(rep(df$category,12),level=c('NF','pRV','RVF')))

df.3 <- subset(df.2, type %in% c('Adipo','CM','EC','Myeloid','FB','SM','PC'))

library(ggpubr)

p <- ggboxplot(df.3,x="group",y="freq",fill="group",group="group")+
	theme_classic() + 
	theme(axis.text.x=element_text(size=16),
	axis.text.y=element_text(size=16),
	axis.title.x=element_text(size=16),
	axis.title.y=element_text(size=16),
	legend.title=element_text(size=16),
	legend.text=element_text(size=16),
	text=element_text(color='black'),
	axis.text=element_text(color='black')) + 
	labs(color='Group',x="Disease",y='Frequency') + 
	facet_wrap(~type,ncol=7) + 
	#stat_compare_means(aes(group=group),comparisons=my_comparisons,method="t.test",ref.group="NF")+
	stat_compare_means(aes(group=group),method="anova")

pdf('./output/deconv_cell_freq_bisque.pdf',width=12,height=5)
p
dev.off()

