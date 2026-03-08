library(Seurat)
library(DESeq2)
library(ggplot2)
library(apeglm)
library(ashr)
library(vsn)
library(dplyr)
library(tibble)
library(clipr)
library(enrichR)

#{loading data}
a <- readRDS('/Volumes/Extreme SSD/Final_Analysis/final_combined_updated.rds')
a <- subset(a,ident='CM')
Nuclei_Meta <- a@meta.data

#

#{r preprocessing}
raw.data <- as.matrix(GetAssayData(a[["RNA"]], slot = "counts"))

col.names<-as.vector(colnames(raw.data))
col.names <- substr(col.names,1,nchar(col.names)-19)

colnames(raw.data)<-col.names

sum.data<-sapply(unique(colnames(raw.data)), function(x) rowSums(raw.data[,grepl(x, colnames(raw.data))]))



col.names <-as.vector(colnames(sum.data))


Nuclei_Meta<-Nuclei_Meta[match(col.names,Nuclei_Meta$patient),]

grouping<-Nuclei_Meta$patient
col.names <- grouping(col.names)
group <- col.names
col.names <- make.names(col.names, unique=TRUE, allow_ = TRUE)

colnames(sum.data) <- col.names



#{r DESeq2 Normalization}
a <- aggregate(x = sum.data, by = list(rownames(sum.data)), FUN = sum, simplify=FALSE)

#remove cardiomyocyte contaminating genes in all populations except cardiomyocytes
#cardio<-read.delim('/Volumes/Extreme SSD/Final_Analysis/cm_genes.tsv') #File with Seurat calculated gene enrichment in cardiomyocyte population
#cardio<-rownames_to_column(cardio,var="genes")
#a <- a[-which(a$Group.1 %in% cardio$X[1:100]),]

NGSnum <- Nuclei_Meta

all(colnames(a[,-1]) == Nuclei_Meta$patient)
b <- Nuclei_Meta[,-1]

rownames(b) <- Nuclei_Meta$Samples
rownames(a) <- a$Group.1
colData <- lapply(rownames(b), as.factor)

NGSnum <- apply(X = a[,-1], MARGIN = 2, FUN = as.numeric)
NGSnum <- apply(X = NGSnum, MARGIN = 2, FUN = round)
rownames(NGSnum) <- a$Group.1
all(colnames(NGSnum) %in% rownames(b))

dds <- DESeqDataSetFromMatrix(countData = NGSnum, colData = b, design = ~0+ group)
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
## Total number of raw counts per sample
#write.table(normalized_counts, file="normalized_counts_myeloid.tsv", sep="\t", quote=F, col.names=NA)


#{r DE}
dds <- DESeq(dds)
colSums(counts(dds, normalized=T)) #Total number of normalized counts per sample

log2cutoff <- 0.5 #or 0.0 to return all genes
qvaluecutoff <- 0.05 #or 0.0 to return all genes


#{r NF vs pRV vs RVF}
#Pairwise comparison: contrast= c('factorName','numeratorLevel','denominatorLevel')

NF_pRV <- results(dds, name="NF_vs_pRV", contrast = c('group','pRV','NF'), test="Wald",independentFiltering=TRUE, alpha=0.05, pAdjustMethod="BH", parallel=FALSE)
resultsNames(dds)
NF_pRV <- lfcShrink(dds, contrast=c('group','NF','pRV'), res=NF_pRV, type = "ashr")

NF_RVF <- results(dds, name="NF_vs_RVF", contrast = c('group','RVF','NF'), test="Wald",independentFiltering=TRUE, alpha=0.05, pAdjustMethod="BH", parallel=FALSE)
resultsNames(dds)
NF_RVF <- lfcShrink(dds, contrast=c('group','NF','RVF'), res=NF_RVF, type = "ashr")

pRV_RVF <- results(dds, name="pRV_vs_RVF", contrast = c('group','RVF','pRV'), test="Wald",independentFiltering=TRUE, alpha=0.05, pAdjustMethod="BH", parallel=FALSE)
resultsNames(dds)
pRV_RVF <- lfcShrink(dds, contrast=c('group','pRV','RVF'), res=pRV_RVF, type = "ashr")

png("~/Downloads/NF_pRV.png",
    width = 10*100, height = 8*100,
    res=100)
plotMA(NF_pRV, ylim = c(-3,3), main="NF_pRV")
dev.off()

png("~/Downloads/pRV_RVF.png",
    width = 10*100, height = 8*100,
    res=100)
plotMA(pRV_RVF, ylim = c(-3,3), main="pRV_RVF")
dev.off()

png("~/Downloads/NF_RVF.png",
    width = 10*100, height = 8*100,
    res=100)
plotMA(NF_RVF, ylim = c(-3,3), main="NF_RVF")
dev.off()



summary(NF_pRV)
sum(NF_pRV$padj < 0.05, na.rm=TRUE) #number of DEGs

png("~/Downloads/Volcano_NF_pRV.png",
    width = 15*100, height = 13*100,
    res=100)
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(NF_pRV, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-4,4)
                       ,ylim=c(0,20)
                       ))
# Add colored points: blue if padj<0.05, red if log2FC>1 and padj<0.05)
with(subset(NF_pRV, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(NF_pRV, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
dev.off()

png("~/Downloads/Volcano_NF_pRV_autoscale.png",
    width = 15*100, height = 13*100,
    res=100)
#reset par
par(mfrow=c(1,1))
with(NF_pRV, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot"))
# Add colored points: blue if padj<0.05, red if log2FC>1 and padj<0.05)
with(subset(NF_pRV, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(NF_pRV, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
dev.off()


NF_pRV <- NF_pRV[order(NF_pRV$padj),] #sort summary list

png("~/Downloads/Volcano_NF_pRV_labels.png",
    width = 15*100, height = 13*100,
    res=100)
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(NF_pRV, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-4,4)
                       ,ylim=c(0,20)
                       ))
# Add colored points: blue if padj<0.05, red if log2FC>1 and padj<0.05)
#with(NF_pRV, points(log2FoldChange, -log10(pvalue), pch=20, col="black"))
with(subset(NF_pRV, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(NF_pRV, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
text(-log10(pvalue)~log2FoldChange,labels=rownames(as.matrix(subset(NF_pRV,padj<.05))), data=NF_pRV, cex=0.8)
dev.off()

png("~/Downloads/Volcano_NF_pRV_autoscale_labels.png",
    width = 15*100, height = 13*100,
    res=100)
#reset par
par(mfrow=c(1,1))
with(NF_pRV, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot"))
# Add colored points: blue if padj<0.05, red if log2FC>1 and padj<0.05)
with(subset(NF_pRV, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(NF_pRV, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
text(-log10(pvalue)~log2FoldChange,labels=rownames(as.matrix(subset(NF_pRV,padj<.05))), data=NF_pRV, cex=0.8)
dev.off()

NF_pRV <- unique(subset(NF_pRV, padj<=qvaluecutoff & abs(log2FoldChange)>=log2cutoff))
write.table(NF_pRV, "~/Downloads/DE_NF_pRV.tsv", sep = '\t')

png("~/Downloads/plotCounts_NF_pRV.png",
    width = 13*100, height = 10*100,
    res=100)
par(mfrow=c(2,3))
plotCounts(dds, gene=rownames(as.matrix(NF_pRV))[1], intgroup="group")
plotCounts(dds, gene=rownames(as.matrix(NF_pRV))[2], intgroup="group")
plotCounts(dds, gene=rownames(as.matrix(NF_pRV))[3], intgroup="group")
plotCounts(dds, gene=rownames(as.matrix(NF_pRV))[4], intgroup="group")
plotCounts(dds, gene=rownames(as.matrix(NF_pRV))[5], intgroup="group")
plotCounts(dds, gene=rownames(as.matrix(NF_pRV))[6], intgroup="group")
dev.off()

summary(NF_RVF)
sum(NF_RVF$padj < 0.05, na.rm=TRUE) #number of DEGs

png("~/Downloads/Volcano_NF_RVF.png",
    width = 15*100, height = 13*100,
    res=100)
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(NF_RVF, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-4,4)
                       ,ylim=c(0,20)
                       ))
# Add colored points: blue if padj<0.05, red if log2FC>1 and padj<0.05)
with(subset(NF_RVF, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(NF_RVF, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
dev.off()

png("~/Downloads/Volcano_NF_RVF_autoscale.png",
    width = 15*100, height = 13*100,
    res=100)
#reset par
par(mfrow=c(1,1))
with(NF_RVF, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot"))
# Add colored points: blue if padj<0.05, red if log2FC>1 and padj<0.05)
with(subset(NF_RVF, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(NF_RVF, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
dev.off()


NF_RVF <- NF_RVF[order(NF_RVF$padj),] #sort summary list

png("~/Downloads/Volcano_NF_RVF_labels.png",
    width = 15*100, height = 13*100,
    res=100)
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(NF_RVF, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-4,4)
                       ,ylim=c(0,20)
                       ))
# Add colored points: blue if padj<0.05, red if log2FC>1 and padj<0.05)
#with(NF_RVF, points(log2FoldChange, -log10(pvalue), pch=20, col="black"))
with(subset(NF_RVF, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(NF_RVF, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
text(-log10(pvalue)~log2FoldChange,labels=rownames(as.matrix(subset(NF_RVF,padj<.05))), data=NF_RVF, cex=0.8)
dev.off()

png("~/Downloads/Volcano_NF_RVF_autoscale_labels.png",
    width = 15*100, height = 13*100,
    res=100)
#reset par
par(mfrow=c(1,1))
with(NF_RVF, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot"))
# Add colored points: blue if padj<0.05, red if log2FC>1 and padj<0.05)
with(subset(NF_RVF, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(NF_RVF, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
text(-log10(pvalue)~log2FoldChange,labels=rownames(as.matrix(subset(NF_RVF,padj<.05))), data=NF_RVF, cex=0.8)
dev.off()

NF_RVF <- unique(subset(NF_RVF, padj<=qvaluecutoff & abs(log2FoldChange)>=log2cutoff))
write.table(NF_RVF, "~/Downloads/DE_NF_RVF.tsv", sep = '\t')

png("~/Downloads/plotCounts_NF_RVF.png",
    width = 13*100, height = 10*100,
    res=100)
par(mfrow=c(2,3))
plotCounts(dds, gene=rownames(as.matrix(NF_RVF))[1], intgroup="group")
plotCounts(dds, gene=rownames(as.matrix(NF_RVF))[2], intgroup="group")
plotCounts(dds, gene=rownames(as.matrix(NF_RVF))[3], intgroup="group")
plotCounts(dds, gene=rownames(as.matrix(NF_RVF))[4], intgroup="group")
plotCounts(dds, gene=rownames(as.matrix(NF_RVF))[5], intgroup="group")
plotCounts(dds, gene=rownames(as.matrix(NF_RVF))[6], intgroup="group")
dev.off()

summary(pRV_RVF)
sum(pRV_RVF$padj < 0.05, na.rm=TRUE) #number of DEGs

png("~/Downloads/Volcano_pRV_RVF.png",
    width = 15*100, height = 13*100,
    res=100)
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(pRV_RVF, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-4,4)
                       ,ylim=c(0,20)
                       ))
# Add colored points: blue if padj<0.05, red if log2FC>1 and padj<0.05)
with(subset(pRV_RVF, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(pRV_RVF, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
dev.off()

png("~/Downloads/Volcano_pRV_RVF_autoscale.png",
    width = 15*100, height = 13*100,
    res=100)
#reset par
par(mfrow=c(1,1))
with(pRV_RVF, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot"))
# Add colored points: blue if padj<0.05, red if log2FC>1 and padj<0.05)
with(subset(pRV_RVF, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(pRV_RVF, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
dev.off()


pRV_RVF <- pRV_RVF[order(pRV_RVF$padj),] #sort summary list

png("~/Downloads/Volcano_pRV_RVF_labels.png",
    width = 15*100, height = 13*100,
    res=100)
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(pRV_RVF, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-4,4)
                       ,ylim=c(0,20)
                       ))
# Add colored points: blue if padj<0.05, red if log2FC>1 and padj<0.05)
#with(pRV_RVF, points(log2FoldChange, -log10(pvalue), pch=20, col="black"))
with(subset(pRV_RVF, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(pRV_RVF, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
text(-log10(pvalue)~log2FoldChange,labels=rownames(as.matrix(subset(pRV_RVF,padj<.05))), data=pRV_RVF, cex=0.8)
dev.off()

png("~/Downloads/Volcano_pRV_RVF_autoscale_labels.png",
    width = 15*100, height = 13*100,
    res=100)
#reset par
par(mfrow=c(1,1))
with(pRV_RVF, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot"))
# Add colored points: blue if padj<0.05, red if log2FC>1 and padj<0.05)
with(subset(pRV_RVF, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(pRV_RVF, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
text(-log10(pvalue)~log2FoldChange,labels=rownames(as.matrix(subset(pRV_RVF,padj<.05))), data=pRV_RVF, cex=0.8)
dev.off()

pRV_RVF <- unique(subset(pRV_RVF, padj<=qvaluecutoff & abs(log2FoldChange)>=log2cutoff))
write.table(pRV_RVF, "~/Downloads/DE_pRV_RVF.tsv", sep = '\t')

png("~/Downloads/plotCounts_pRV_RVF.png",
    width = 13*100, height = 10*100,
    res=100)
par(mfrow=c(2,3))
plotCounts(dds, gene=rownames(as.matrix(pRV_RVF))[1], intgroup="group")
plotCounts(dds, gene=rownames(as.matrix(pRV_RVF))[2], intgroup="group")
plotCounts(dds, gene=rownames(as.matrix(pRV_RVF))[3], intgroup="group")
plotCounts(dds, gene=rownames(as.matrix(pRV_RVF))[4], intgroup="group")
plotCounts(dds, gene=rownames(as.matrix(pRV_RVF))[5], intgroup="group")
plotCounts(dds, gene=rownames(as.matrix(pRV_RVF))[6], intgroup="group")
dev.off()




#Next steps in exploring these data...BLAST to database to find associated gene function



###{r PCA}
#First we need to transform the raw count data
#vst function will perform variance stabilizing transformation or rlog

vsdata <- vst(dds, blind=FALSE)
rlongdata <- rlog(dds, blind = FALSE)

png("~/Downloads/PCA_Cond.png", width = 13*100, height = 10*100, res=100)
plotPCA(rlongdata, intgroup="group") +geom_point(size=12)+ theme(panel.background = element_rect(fill="white"))
dev.off()


###enrichR
library(DOSE)
library(enrichR)
library(tidyverse)

websiteLive <- getOption("enrichR.live")
if (websiteLive) {
    listEnrichrSites()
    setEnrichrSite("Enrichr") # Human genes   
}

dbs <- c("ChEA_2022","WikiPathway_2023_Human","Reactome_2016","GO_Biological_Process_2023")
if (websiteLive) {
    enriched <- enrichr(rownames(subset(NF_pRV,log2FoldChange > 1)), dbs)
	p1<- ggplot(enriched[[1]][order(enriched[[1]]$Combined.Score,decreasing=T),][rev(1:10),], (aes(x=Combined.Score, y=fct_inorder(sapply(strsplit(Term," "), `[`, 1)), color = as.numeric(P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + ggtitle('ChEA Up')

	p2 <- ggplot(enriched[[2]][order(enriched[[2]]$Combined.Score,decreasing=T),][rev(1:10),], (aes(x=Combined.Score, y=fct_inorder(sapply(strsplit(Term," WP"), `[`, 1)), color = as.numeric(P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic() + ggtitle('WikiPathway Up')

	p3 <- ggplot(enriched[[3]][order(enriched[[3]]$Combined.Score,decreasing=T),][rev(1:10),], (aes(x=Combined.Score, y=fct_inorder(sapply(strsplit(Term," Homo sapiens"), `[`, 1)), color = as.numeric(P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic() + ggtitle('Reactome Up')

	p4 <- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:10),], (aes(x=Combined.Score, y=fct_inorder(sapply(strsplit(Term," \\(GO"), `[`, 1)), color = as.numeric(P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + ggtitle('GO Biological Process Up')


	png("~/Downloads/Enrich_NF_pRV_up.png",
    width = 30*100, height = 10*100,
    res=100)
plot_grid(p1,p2,p3,p4,align="v",ncol=2)
dev.off()

enriched <- enrichr(rownames(subset(NF_RVF,log2FoldChange > 1)), dbs)
	p1<- ggplot(enriched[[1]][order(enriched[[1]]$Combined.Score,decreasing=T),][rev(1:10),], (aes(x=Combined.Score, y=fct_inorder(sapply(strsplit(Term," "), `[`, 1)), color = as.numeric(P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic() + ggtitle('ChEA')

	p2 <- ggplot(enriched[[2]][order(enriched[[2]]$Combined.Score,decreasing=T),][rev(1:10),], (aes(x=Combined.Score, y=fct_inorder(sapply(strsplit(Term," WP"), `[`, 1)), color = as.numeric(P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic() + ggtitle('WikiPathway')

	p3 <- ggplot(enriched[[3]][order(enriched[[3]]$Combined.Score,decreasing=T),][rev(1:10),], (aes(x=Combined.Score, y=fct_inorder(sapply(strsplit(Term," Homo sapiens"), `[`, 1)), color = as.numeric(P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic() + ggtitle('Reactome')

	p4 <- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:10),], (aes(x=Combined.Score, y=fct_inorder(sapply(strsplit(Term," \\(GO"), `[`, 1)), color = as.numeric(P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + ggtitle('GO Biological Process')


	png("~/Downloads/Enrich_NF_RVF_up.png",
    width = 30*100, height = 10*100,
    res=100)
plot_grid(p1,p2,p3,p4,align="v",ncol=2)
dev.off()

enriched <- enrichr(rownames(subset(pRV_RVF,log2FoldChange > 1)), dbs)
	p1<- ggplot(enriched[[1]][order(enriched[[1]]$Combined.Score,decreasing=T),][rev(1:10),], (aes(x=Combined.Score, y=fct_inorder(sapply(strsplit(Term," "), `[`, 1)), color = as.numeric(P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic() + ggtitle('ChEA')

	p2 <- ggplot(enriched[[2]][order(enriched[[2]]$Combined.Score,decreasing=T),][rev(1:10),], (aes(x=Combined.Score, y=fct_inorder(sapply(strsplit(Term," WP"), `[`, 1)), color = as.numeric(P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic() + ggtitle('WikiPathway')

	p3 <- ggplot(enriched[[3]][order(enriched[[3]]$Combined.Score,decreasing=T),][rev(1:10),], (aes(x=Combined.Score, y=fct_inorder(sapply(strsplit(Term," Homo sapiens"), `[`, 1)), color = as.numeric(P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic() + ggtitle('Reactome')
	p3 <- ggplot(enriched[[3]][order(enriched[[3]]$Combined.Score,decreasing=T),][rev(1:10),], (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic() + ggtitle('Reactome')

	p4 <- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:10),], (aes(x=Combined.Score, y=fct_inorder(sapply(strsplit(Term," \\(GO"), `[`, 1)), color = as.numeric(P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic() + ggtitle('GO Biological Process')


	png("~/Downloads/Enrich_pRV_RVF_up.png",
    width = 30*100, height = 10*100,
    res=100)
plot_grid(p1,p2,p3,p4,align="v",ncol=2)
dev.off()

    enriched <- enrichr(rownames(subset(NF_pRV,log2FoldChange < -1)), dbs)
	p1<- ggplot(enriched[[1]][order(enriched[[1]]$Combined.Score,decreasing=T),][rev(1:10),], (aes(x=Combined.Score, y=fct_inorder(sapply(strsplit(Term," "), `[`, 1)), color = as.numeric(P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + ggtitle('ChEA Up')

	p2 <- ggplot(enriched[[2]][order(enriched[[2]]$Combined.Score,decreasing=T),][rev(1:10),], (aes(x=Combined.Score, y=fct_inorder(sapply(strsplit(Term," WP"), `[`, 1)), color = as.numeric(P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic() + ggtitle('WikiPathway Up')

	p3 <- ggplot(enriched[[3]][order(enriched[[3]]$Combined.Score,decreasing=T),][rev(1:10),], (aes(x=Combined.Score, y=fct_inorder(sapply(strsplit(Term," Homo sapiens"), `[`, 1)), color = as.numeric(P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic() + ggtitle('Reactome Up')

	p4 <- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:10),], (aes(x=Combined.Score, y=fct_inorder(sapply(strsplit(Term," \\(GO"), `[`, 1)), color = as.numeric(P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + ggtitle('GO Biological Process Up')


	png("~/Downloads/Enrich_NF_pRV_down.png",
    width = 30*100, height = 10*100,
    res=100)
plot_grid(p1,p2,p3,p4,align="v",ncol=2)
dev.off()

enriched <- enrichr(rownames(subset(NF_RVF,log2FoldChange < -1)), dbs)
	p1<- ggplot(enriched[[1]][order(enriched[[1]]$Combined.Score,decreasing=T),][rev(1:10),], (aes(x=Combined.Score, y=fct_inorder(sapply(strsplit(Term," "), `[`, 1)), color = as.numeric(P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic() + ggtitle('ChEA')

	p2 <- ggplot(enriched[[2]][order(enriched[[2]]$Combined.Score,decreasing=T),][rev(1:10),], (aes(x=Combined.Score, y=fct_inorder(sapply(strsplit(Term," WP"), `[`, 1)), color = as.numeric(P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic() + ggtitle('WikiPathway')

	p3 <- ggplot(enriched[[3]][order(enriched[[3]]$Combined.Score,decreasing=T),][rev(1:10),], (aes(x=Combined.Score, y=fct_inorder(sapply(strsplit(Term," Homo sapiens"), `[`, 1)), color = as.numeric(P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic() + ggtitle('Reactome')

	p4 <- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:10),], (aes(x=Combined.Score, y=fct_inorder(sapply(strsplit(Term," \\(GO"), `[`, 1)), color = as.numeric(P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + ggtitle('GO Biological Process')


	png("~/Downloads/Enrich_NF_RVF_down.png",
    width = 30*100, height = 10*100,
    res=100)
plot_grid(p1,p2,p3,p4,align="v",ncol=2)
dev.off()

enriched <- enrichr(rownames(subset(pRV_RVF,log2FoldChange < -1)), dbs)
	p1<- ggplot(enriched[[1]][order(enriched[[1]]$Combined.Score,decreasing=T),][rev(1:10),], (aes(x=Combined.Score, y=fct_inorder(sapply(strsplit(Term," "), `[`, 1)), color = as.numeric(P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic() + ggtitle('ChEA')

	p2 <- ggplot(enriched[[2]][order(enriched[[2]]$Combined.Score,decreasing=T),][rev(1:10),], (aes(x=Combined.Score, y=fct_inorder(sapply(strsplit(Term," WP"), `[`, 1)), color = as.numeric(P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic() + ggtitle('WikiPathway')

	p3 <- ggplot(enriched[[3]][order(enriched[[3]]$Combined.Score,decreasing=T),][rev(1:10),], (aes(x=Combined.Score, y=fct_inorder(sapply(strsplit(Term," Homo sapiens"), `[`, 1)), color = as.numeric(P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic() + ggtitle('Reactome')

	p4 <- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:10),], (aes(x=Combined.Score, y=fct_inorder(sapply(strsplit(Term," \\(GO"), `[`, 1)), color = as.numeric(P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic() + ggtitle('GO Biological Process')


	png("~/Downloads/Enrich_pRV_RVF_down.png",
    width = 30*100, height = 10*100,
    res=100)
plot_grid(p1,p2,p3,p4,align="v",ncol=2)
dev.off()

}


