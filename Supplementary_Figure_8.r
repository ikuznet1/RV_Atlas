library(Seurat)
library(hdWGCNA)
library(ggeasy)
library(dplyr)


source('~/Downloads/hdWGCNA_TOM/spatial_functions.R')



#######################################
#############  FIGURE 8A  #############
#######################################
M2 <- readRDS('~/Downloads/hdWGCNA_TOM/RV_Peds_merge.rds')
M2$condition[M2$condition == "NF"] = "pRV" 
M2$condition[M2$condition == "Donor"] = "NF" 
M2$condition[M2$condition == "SystolicHF"] = "RVF" 

M2$group[is.na(M2$group)] <- M2$condition[is.na(M2$group)]

M1 <- subset(M2,group=='NF')

rm(M2)

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'sn_Peds_adult_NF.pdf'), width=5, height=5)
PlotEmbedding(M1,group.by='origin',point_size=0.2,plot_under=TRUE,plot_theme=umap_theme()+NoLegend(),raster_dpi=400,raster_scale=0.5)
dev.off()



#######################################
#############  FIGURE 8B  #############
#######################################

t(table(M1$CombinedNames,M1$origin))/rowSums(t(table(M1$CombinedNames,M1$origin)))



cells <- table(M1$CombinedNames,M1$patient)
cells <- cells[c('CM','FB','EC','Myeloid'),]
cells <- sweep(cells,2,colSums(cells),'/')
cells <- data.frame(cells)

cells$origin<-M1$origin[match(cells$Var2,M1$patient)]

cells$origin<-as.factor(cells$origin)

library(ggpubr)
pdf('~/Downloads/hdWGCNA_TOM/Peds_RV_NF_clust_freq.pdf',width=8,height=5)
p <- ggboxplot(cells[length(cells$origin):1,],x="origin",y="Freq",fill="origin",group="origin")+
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
  facet_wrap(~Var1,ncol=7) + 
  stat_compare_means(aes(group=origin),method="anova")
p
dev.off()


#######################################
#############  FIGURE 8C  #############
#######################################


########Embed bulk in single nuc peds


######Load module
consensus_modules <- read.csv("~/Downloads/hdWGCNA_TOM/bulk_heart_modules.csv")
consensus_modules <- consensus_modules[,1:3]
consensus_modules <- subset(consensus_modules, gene_name %in% rownames(M2))


# remove duplicate gene names
consensus_modules <- consensus_modules[match(unique(consensus_modules$gene_name), consensus_modules$gene_name),]


M2 <- M1

#A bit hacky since ref here should be NULL and wgcna_name should be None but function errors then. This will be overwritten by modules=consensus_modules anyway
DefaultAssay(M2) <- 'RNA'
M2<-FindVariableFeatures(M2)
M2 <- ScaleData(M2,block.size=1000)

M2 <- SetupForWGCNA(
  M2,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "Cardiomyocyte" # the name of the hdWGCNA experiment
)


M2 <- ProjectModules(
  M2,
  modules = consensus_modules,
  group.by.vars = "patient",
  seurat_ref = M2,
  wgcna_name = "Cardiomyocyte",
  wgcna_name_proj = 'bulk2sn'
)



M2 <- SetActiveWGCNA(M2, 'bulk2sn')
mapping <- labels2colors(1:100)
MEs <- GetMEs(M2, harmonized=TRUE)
mods <- colnames(MEs); mods <- mods[mods != 'grey']
mods_num <- paste0('M',match(mods,mapping))
all_signif <- c('M1','M2','M3','M4','M5','M8','M10','M11','M12','M14','M20','M25','M26','M28')


colnames(MEs)<-paste0('M',match(colnames(MEs),mapping))
M2@meta.data <- cbind(M2@meta.data, MEs)
M2 <- SetIdent(M2, value = "CombinedNames")

library(dplyr)
score_calc <- consensus_modules %>% group_by(module) %>% group_split()
module_colors <- unique(unlist(lapply(score_calc,'[[','module')))
module_colors <- paste0('M',match(module_colors,mapping))
#saveRDS(M2, '~/Downloads/hdWGCNA_TOM/scWGCNA_RV_Peds_bulk2sn_projection.rds')
#M2<- readRDS('~/Downloads/hdWGCNA_TOM/scWGCNA_RV_Peds_bulk2sn_projection.rds')

DefaultAssay(M2) <- 'SCT'



#rm(seurat_ref)
#gc()
#seurat_ref<-readRDS('/Volumes/Extreme SSD/Final_Analysis/CellTypes/Post_R3_FINAL_with_counts.rds')
#seurat_ref <- SetIdent(seurat_ref, value = "Names")
#seurat_ref@meta.data <- cbind(seurat_ref@meta.data, MEs)



M2 <- AddModuleScore(M2,lapply(score_calc,'[[','gene_name'),name="module_score")


cols_current <- colnames(M2@meta.data)
cols_current[startsWith(colnames(M2@meta.data),'module_score')] <- paste0('module_',module_colors)
colnames(M2@meta.data) <- cols_current

M2$origin[M2$origin] = 'Peds'
M2$origin[M2$origin == FALSE] = 'Adult'

M2$CombinedNamesSplit <- paste0(M2$CombinedNames,'_',M2$origin)

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'RV_Peds_NF_Dot.pdf'), width=7, height=5)
p <- DotPlot(M2,paste0('module_',all_signif),group.by='CombinedNamesSplit',dot.min=0,col.min=0,col.max=2,idents=c("CM","EC","FB","Myeloid","PC","SM")) +
  RotatedAxis() + ylab('')+ xlab('')+
  scale_color_gradient2(high='red', mid='grey95', low='blue') +
  theme(
    panel.border = element_rect(size=1,fill=NA, color='black'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()
) 
p
dev.off()


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'RV_Peds_Dot_NF_ordered.pdf'), width=7, height=5)

p <- DotPlot(M2,paste0('module_',c('M20','M5','M1','M3','M4','M8','M2','M12','M25','M26','M10','M28','M14','M11')),group.by='CombinedNamesSplit',dot.min=0,col.min=0,col.max=2,idents=c("CM","EC","FB","Myeloid","PC","SM")) +
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
#############  FIGURE 8D  #############
#######################################
library(VennDiagram)

cm_peds <- colnames(subset(M1, origin==TRUE & CombinedNames == 'CM'))
all_other_peds <- colnames(subset(M1, origin==TRUE & CombinedNames != 'CM'))
cm_adult <- colnames(subset(M1, origin==FALSE & CombinedNames == 'CM'))
all_other_adult <- colnames(subset(M1, origin==FALSE & CombinedNames != 'CM'))


peds_cm_marks <- FindMarkers(M1, ident.1 = cm_peds, ident.2 = all_other_peds,recorrect_umi = FALSE,logfc.threshold=0)
adult_cm_marks <- FindMarkers(M1, ident.1 = cm_adult, ident.2 = all_other_adult,recorrect_umi = FALSE,logfc.threshold=0)

peds_cm_marks_s <- subset(peds_cm_marks,p_val_adj < 0.05)
adult_cm_marks_s <- subset(adult_cm_marks,p_val_adj < 0.05)


x <- list(adult = rownames(adult_cm_marks_s), peds = rownames(peds_cm_marks_s))

venn.diagram(x, filename = paste0('~/Downloads/hdWGCNA_TOM/', 'RV_Peds_CM_Venn.png'),
	imagetype = "png",fill = c("#56B4E9", "#E69F00"))


shared <- intersect(rownames(adult_cm_marks_s),rownames(peds_cm_marks_s))
adult_unique <- setdiff(rownames(adult_cm_marks_s),shared)
peds_unique <- setdiff(rownames(peds_cm_marks_s),shared)

dataset1 <- data.frame(Peds=peds_cm_marks[adult_unique,]$avg_log2FC,Adult=adult_cm_marks[adult_unique,]$avg_log2FC)
rownames(dataset1) <- adult_unique
dataset1 <- dataset1[complete.cases(dataset1), ]


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'Peds_vs_Adult_Adult_CM_Specific.pdf'), width=8, height=8)
ggplot(dataset1, aes(x = Peds, y=Adult)) + geom_point() + 
  geom_text_repel(label=rownames(dataset1),max.overlaps = 10) + theme_classic()
dev.off()

dataset2 <- data.frame(Peds=peds_cm_marks[peds_unique,]$avg_log2FC,Adult=adult_cm_marks[peds_unique,]$avg_log2FC)
rownames(dataset2) <- peds_unique
dataset2 <- dataset2[complete.cases(dataset2), ]


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'Peds_vs_Adult_Peds_CM_Specific.pdf'), width=8, height=8)
ggplot(dataset2, aes(x = Peds, y=Adult)) + geom_point() + 
  geom_text_repel(label=rownames(dataset2),max.overlaps = 10) + theme_classic()
dev.off()


dataset1$origin <- 'Adult'
dataset2$origin <- 'Peds'
dataset <- rbind(dataset1,dataset2)

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'Peds_vs_Adult_CM_Specific.pdf'), width=10, height=8)
ggplot(dataset, aes(x = Adult, y=Peds,color = origin)) + geom_point() + 
  geom_text_repel(label=rownames(dataset),max.overlaps = 15) + theme_classic() + 
  scale_colour_manual(name="",values = c("Adult" = "red","Peds" = "blue"))
dev.off()



#######################################
#############  FIGURE 8E  #############
#######################################

dbs <- c("ChEA_2022","WikiPathway_2023_Human","Reactome_2016","GO_Biological_Process_2023")

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

adult_cm_marks_s_u <- subset(adult_cm_marks_s,avg_log2FC > 0)
adult_cm_marks_s_d <- subset(adult_cm_marks_s,avg_log2FC < 0)


peds_cm_marks_s_u <- subset(peds_cm_marks_s,avg_log2FC > 0)
peds_cm_marks_s_d <- subset(peds_cm_marks_s,avg_log2FC < 0)

adult_unique_u <- rownames(adult_cm_marks_s_u[intersect(adult_unique,rownames(adult_cm_marks_s_u)),])
adult_unique_d <- rownames(adult_cm_marks_s_d[intersect(adult_unique,rownames(adult_cm_marks_s_d)),])

# Make sure genes actually exist in peds dataset
adult_unique_u <- intersect(adult_unique_u,rownames(peds_cm_marks))
adult_unique_d <- intersect(adult_unique_d,rownames(peds_cm_marks))


peds_unique_u <- rownames(peds_cm_marks_s_u[intersect(peds_unique,rownames(peds_cm_marks_s_u)),])
peds_unique_d <- rownames(peds_cm_marks_s_d[intersect(peds_unique,rownames(peds_cm_marks_s_d)),])

# Make sure genes actually exist in adult dataset
peds_unique_u <- intersect(peds_unique_u,rownames(adult_cm_marks))
peds_unique_d <- intersect(peds_unique_d,rownames(adult_cm_marks))


enriched <- enrichr(adult_unique_u, dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
enriched_au <- enriched
pdf('~/Downloads/hdWGCNA_TOM/CM_Adults_vs_Peds_Adult_unique_up.pdf',width=5,height=2.5)
p1<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),], 
  (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), 
  size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + 
  ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + 
  ggtitle('GO BP') + 
  scale_y_discrete(labels= fct_inorder(
    wrapText(sapply(
      strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),]$Term," \\(GO"),
         `[`, 1),35))) + 
  theme(axis.text=element_text(colour="black"))+
  scale_color_stepsn(colors=rev(magma(256)))
p1
dev.off()

enriched <- enrichr(adult_unique_d, dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
enriched_ad <- enriched
pdf('~/Downloads/hdWGCNA_TOM/CM_Adults_vs_Peds_Adult_unique_down.pdf',width=5,height=2.5)
p2<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),], 
  (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), 
  size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + 
  ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + 
  ggtitle('GO BP') + 
  scale_y_discrete(labels= fct_inorder(
    wrapText(sapply(
      strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),]$Term," \\(GO"),
         `[`, 1),35))) + 
  theme(axis.text=element_text(colour="black"))+
  scale_color_stepsn(colors=rev(magma(256)))
p2
dev.off()

enriched <- enrichr(peds_unique_u, dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
enriched_pu <- enriched
pdf('~/Downloads/hdWGCNA_TOM/CM_Adults_vs_Peds_Peds_unique_up.pdf',width=5,height=2.5)
p1<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),], 
  (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), 
  size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + 
  ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + 
  ggtitle('GO BP') + 
  scale_y_discrete(labels= fct_inorder(
    wrapText(sapply(
      strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),]$Term," \\(GO"),
         `[`, 1),35))) + 
  theme(axis.text=element_text(colour="black"))+
  scale_color_stepsn(colors=rev(magma(256)))
p1
dev.off()

enriched <- enrichr(peds_unique_d, dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
enriched_pd <- enriched
pdf('~/Downloads/hdWGCNA_TOM/CM_Adults_vs_Peds_Peds_unique_down.pdf',width=5,height=2.5)
p2<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),], 
  (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), 
  size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + 
  ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + 
  ggtitle('GO BP') + 
  scale_y_discrete(labels= fct_inorder(
    wrapText(sapply(
      strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),]$Term," \\(GO"),
         `[`, 1),35))) + 
  theme(axis.text=element_text(colour="black"))+
  scale_color_stepsn(colors=rev(magma(256)))
p2
dev.off()








#######################################
#############  FIGURE 8F  #############
#######################################
library(VennDiagram)

cm_peds <- colnames(subset(M1, origin==TRUE & CombinedNames == 'Myeloid'))
all_other_peds <- colnames(subset(M1, origin==TRUE & CombinedNames != 'Myeloid'))
cm_adult <- colnames(subset(M1, origin==FALSE & CombinedNames == 'Myeloid'))
all_other_adult <- colnames(subset(M1, origin==FALSE & CombinedNames != 'Myeloid'))


peds_cm_marks <- FindMarkers(M1, ident.1 = cm_peds, ident.2 = all_other_peds,recorrect_umi = FALSE,logfc.threshold=0)
adult_cm_marks <- FindMarkers(M1, ident.1 = cm_adult, ident.2 = all_other_adult,recorrect_umi = FALSE,logfc.threshold=0)

peds_cm_marks_s <- subset(peds_cm_marks,p_val_adj < 0.05)
adult_cm_marks_s <- subset(adult_cm_marks,p_val_adj < 0.05)


x <- list(adult = rownames(adult_cm_marks_s), peds = rownames(peds_cm_marks_s))

venn.diagram(x, filename = paste0('~/Downloads/hdWGCNA_TOM/', 'RV_Peds_Myeloid_Venn.png'),
	imagetype = "png",fill = c("#56B4E9", "#E69F00"))


shared <- intersect(rownames(adult_cm_marks_s),rownames(peds_cm_marks_s))
adult_unique <- setdiff(rownames(adult_cm_marks_s),shared)
peds_unique <- setdiff(rownames(peds_cm_marks_s),shared)

dataset1 <- data.frame(Peds=peds_cm_marks[adult_unique,]$avg_log2FC,Adult=adult_cm_marks[adult_unique,]$avg_log2FC)
rownames(dataset1) <- adult_unique
dataset1 <- dataset1[complete.cases(dataset1), ]


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'Peds_vs_Adult_Adult_Myeloid_Specific.pdf'), width=8, height=8)
ggplot(dataset1, aes(x = Peds, y=Adult)) + geom_point() + 
  geom_text_repel(label=rownames(dataset1),max.overlaps = 10) + theme_classic()
dev.off()

dataset2 <- data.frame(Peds=peds_cm_marks[peds_unique,]$avg_log2FC,Adult=adult_cm_marks[peds_unique,]$avg_log2FC)
rownames(dataset2) <- peds_unique
dataset2 <- dataset2[complete.cases(dataset2), ]


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'Peds_vs_Adult_Peds_Myeloid_Specific.pdf'), width=8, height=8)
ggplot(dataset2, aes(x = Peds, y=Adult)) + geom_point() + 
  geom_text_repel(label=rownames(dataset2),max.overlaps = 10) + theme_classic()
dev.off()


dataset1$origin <- 'Adult'
dataset2$origin <- 'Peds'
dataset <- rbind(dataset1,dataset2)

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'Peds_vs_Adult_Myeloid_Specific.pdf'), width=10, height=8)
ggplot(dataset, aes(x = Adult, y=Peds,color = origin)) + geom_point() + 
  geom_text_repel(label=rownames(dataset),max.overlaps = 15) + theme_classic() + 
  scale_colour_manual(name="",values = c("Adult" = "red","Peds" = "blue"))
dev.off()



dbs <- c("ChEA_2022","WikiPathway_2023_Human","Reactome_2016","GO_Biological_Process_2023")

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

adult_cm_marks_s_u <- subset(adult_cm_marks_s,avg_log2FC > 0)
adult_cm_marks_s_d <- subset(adult_cm_marks_s,avg_log2FC < 0)


peds_cm_marks_s_u <- subset(peds_cm_marks_s,avg_log2FC > 0)
peds_cm_marks_s_d <- subset(peds_cm_marks_s,avg_log2FC < 0)

adult_unique_u <- rownames(adult_cm_marks_s_u[intersect(adult_unique,rownames(adult_cm_marks_s_u)),])
adult_unique_d <- rownames(adult_cm_marks_s_d[intersect(adult_unique,rownames(adult_cm_marks_s_d)),])

# Make sure genes actually exist in peds dataset
adult_unique_u <- intersect(adult_unique_u,rownames(peds_cm_marks))
adult_unique_d <- intersect(adult_unique_d,rownames(peds_cm_marks))


peds_unique_u <- rownames(peds_cm_marks_s_u[intersect(peds_unique,rownames(peds_cm_marks_s_u)),])
peds_unique_d <- rownames(peds_cm_marks_s_d[intersect(peds_unique,rownames(peds_cm_marks_s_d)),])

# Make sure genes actually exist in adult dataset
peds_unique_u <- intersect(peds_unique_u,rownames(adult_cm_marks))
peds_unique_d <- intersect(peds_unique_d,rownames(adult_cm_marks))


enriched <- enrichr(adult_unique_u, dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
enriched_au <- enriched
pdf('~/Downloads/hdWGCNA_TOM/Myeloid_Adults_vs_Peds_Adult_unique_up.pdf',width=5,height=2.5)
p1<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),], 
  (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), 
  size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + 
  ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + 
  ggtitle('GO BP') + 
  scale_y_discrete(labels= fct_inorder(
    wrapText(sapply(
      strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),]$Term," \\(GO"),
         `[`, 1),35))) + 
  theme(axis.text=element_text(colour="black"))+
  scale_color_stepsn(colors=rev(magma(256)))
p1
dev.off()

enriched <- enrichr(adult_unique_d, dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
enriched_ad <- enriched
pdf('~/Downloads/hdWGCNA_TOM/Myeloid_Adults_vs_Peds_Adult_unique_down.pdf',width=5,height=2.5)
p2<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),], 
  (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), 
  size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + 
  ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + 
  ggtitle('GO BP') + 
  scale_y_discrete(labels= fct_inorder(
    wrapText(sapply(
      strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),]$Term," \\(GO"),
         `[`, 1),35))) + 
  theme(axis.text=element_text(colour="black"))+
  scale_color_stepsn(colors=rev(magma(256)))
p2
dev.off()

enriched <- enrichr(peds_unique_u, dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
enriched_pu <- enriched
pdf('~/Downloads/hdWGCNA_TOM/Myeloid_Adults_vs_Peds_Peds_unique_up.pdf',width=5,height=2.5)
p1<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),], 
  (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), 
  size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + 
  ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + 
  ggtitle('GO BP') + 
  scale_y_discrete(labels= fct_inorder(
    wrapText(sapply(
      strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),]$Term," \\(GO"),
         `[`, 1),35))) + 
  theme(axis.text=element_text(colour="black"))+
  scale_color_stepsn(colors=rev(magma(256)))
p1
dev.off()

enriched <- enrichr(peds_unique_d, dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
enriched_pd <- enriched
pdf('~/Downloads/hdWGCNA_TOM/Myeloid_Adults_vs_Peds_Peds_unique_down.pdf',width=5,height=2.5)
p2<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),], 
  (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), 
  size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + 
  ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + 
  ggtitle('GO BP') + 
  scale_y_discrete(labels= fct_inorder(
    wrapText(sapply(
      strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),]$Term," \\(GO"),
         `[`, 1),35))) + 
  theme(axis.text=element_text(colour="black"))+
  scale_color_stepsn(colors=rev(magma(256)))
p2
dev.off()




#######################################
#############  FIGURE 8G  #############
#######################################
library(VennDiagram)

cm_peds <- colnames(subset(M1, origin==TRUE & CombinedNames == 'EC'))
all_other_peds <- colnames(subset(M1, origin==TRUE & CombinedNames != 'EC'))
cm_adult <- colnames(subset(M1, origin==FALSE & CombinedNames == 'EC'))
all_other_adult <- colnames(subset(M1, origin==FALSE & CombinedNames != 'EC'))


peds_cm_marks <- FindMarkers(M1, ident.1 = cm_peds, ident.2 = all_other_peds,recorrect_umi = FALSE,logfc.threshold=0)
adult_cm_marks <- FindMarkers(M1, ident.1 = cm_adult, ident.2 = all_other_adult,recorrect_umi = FALSE,logfc.threshold=0)

peds_cm_marks_s <- subset(peds_cm_marks,p_val_adj < 0.05)
adult_cm_marks_s <- subset(adult_cm_marks,p_val_adj < 0.05)


x <- list(adult = rownames(adult_cm_marks_s), peds = rownames(peds_cm_marks_s))

venn.diagram(x, filename = paste0('~/Downloads/hdWGCNA_TOM/', 'RV_Peds_EC_Venn.png'),
	imagetype = "png",fill = c("#56B4E9", "#E69F00"))


shared <- intersect(rownames(adult_cm_marks_s),rownames(peds_cm_marks_s))
adult_unique <- setdiff(rownames(adult_cm_marks_s),shared)
peds_unique <- setdiff(rownames(peds_cm_marks_s),shared)

dataset1 <- data.frame(Peds=peds_cm_marks[adult_unique,]$avg_log2FC,Adult=adult_cm_marks[adult_unique,]$avg_log2FC)
rownames(dataset1) <- adult_unique
dataset1 <- dataset1[complete.cases(dataset1), ]


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'Peds_vs_Adult_Adult_EC_Specific.pdf'), width=8, height=8)
ggplot(dataset1, aes(x = Peds, y=Adult)) + geom_point() + 
  geom_text_repel(label=rownames(dataset1),max.overlaps = 10) + theme_classic()
dev.off()

dataset2 <- data.frame(Peds=peds_cm_marks[peds_unique,]$avg_log2FC,Adult=adult_cm_marks[peds_unique,]$avg_log2FC)
rownames(dataset2) <- peds_unique
dataset2 <- dataset2[complete.cases(dataset2), ]


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'Peds_vs_Adult_Peds_EC_Specific.pdf'), width=8, height=8)
ggplot(dataset2, aes(x = Peds, y=Adult)) + geom_point() + 
  geom_text_repel(label=rownames(dataset2),max.overlaps = 10) + theme_classic()
dev.off()


dataset1$origin <- 'Adult'
dataset2$origin <- 'Peds'
dataset <- rbind(dataset1,dataset2)

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'Peds_vs_Adult_EC_Specific.pdf'), width=10, height=8)
ggplot(dataset, aes(x = Adult, y=Peds,color = origin)) + geom_point() + 
  geom_text_repel(label=rownames(dataset),max.overlaps = 15) + theme_classic() + 
  scale_colour_manual(name="",values = c("Adult" = "red","Peds" = "blue"))
dev.off()



dbs <- c("ChEA_2022","WikiPathway_2023_Human","Reactome_2016","GO_Biological_Process_2023")

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

adult_cm_marks_s_u <- subset(adult_cm_marks_s,avg_log2FC > 0)
adult_cm_marks_s_d <- subset(adult_cm_marks_s,avg_log2FC < 0)


peds_cm_marks_s_u <- subset(peds_cm_marks_s,avg_log2FC > 0)
peds_cm_marks_s_d <- subset(peds_cm_marks_s,avg_log2FC < 0)

adult_unique_u <- rownames(adult_cm_marks_s_u[intersect(adult_unique,rownames(adult_cm_marks_s_u)),])
adult_unique_d <- rownames(adult_cm_marks_s_d[intersect(adult_unique,rownames(adult_cm_marks_s_d)),])

# Make sure genes actually exist in peds dataset
adult_unique_u <- intersect(adult_unique_u,rownames(peds_cm_marks))
adult_unique_d <- intersect(adult_unique_d,rownames(peds_cm_marks))


peds_unique_u <- rownames(peds_cm_marks_s_u[intersect(peds_unique,rownames(peds_cm_marks_s_u)),])
peds_unique_d <- rownames(peds_cm_marks_s_d[intersect(peds_unique,rownames(peds_cm_marks_s_d)),])

# Make sure genes actually exist in adult dataset
peds_unique_u <- intersect(peds_unique_u,rownames(adult_cm_marks))
peds_unique_d <- intersect(peds_unique_d,rownames(adult_cm_marks))


enriched <- enrichr(adult_unique_u, dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
enriched_au <- enriched
pdf('~/Downloads/hdWGCNA_TOM/EC_Adults_vs_Peds_Adult_unique_up.pdf',width=5,height=2.5)
p1<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),], 
  (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), 
  size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + 
  ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + 
  ggtitle('GO BP') + 
  scale_y_discrete(labels= fct_inorder(
    wrapText(sapply(
      strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),]$Term," \\(GO"),
         `[`, 1),35))) + 
  theme(axis.text=element_text(colour="black"))+
  scale_color_stepsn(colors=rev(magma(256)))
p1
dev.off()

enriched <- enrichr(adult_unique_d, dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
enriched_ad <- enriched
pdf('~/Downloads/hdWGCNA_TOM/EC_Adults_vs_Peds_Adult_unique_down.pdf',width=5,height=2.5)
p2<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),], 
  (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), 
  size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + 
  ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + 
  ggtitle('GO BP') + 
  scale_y_discrete(labels= fct_inorder(
    wrapText(sapply(
      strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),]$Term," \\(GO"),
         `[`, 1),35))) + 
  theme(axis.text=element_text(colour="black"))+
  scale_color_stepsn(colors=rev(magma(256)))
p2
dev.off()

enriched <- enrichr(peds_unique_u, dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
enriched_pu <- enriched
pdf('~/Downloads/hdWGCNA_TOM/EC_Adults_vs_Peds_Peds_unique_up.pdf',width=5,height=2.5)
p1<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),], 
  (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), 
  size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + 
  ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + 
  ggtitle('GO BP') + 
  scale_y_discrete(labels= fct_inorder(
    wrapText(sapply(
      strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),]$Term," \\(GO"),
         `[`, 1),35))) + 
  theme(axis.text=element_text(colour="black"))+
  scale_color_stepsn(colors=rev(magma(256)))
p1
dev.off()

enriched <- enrichr(peds_unique_d, dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
enriched_pd <- enriched
pdf('~/Downloads/hdWGCNA_TOM/EC_Adults_vs_Peds_Peds_unique_down.pdf',width=5,height=2.5)
p2<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),], 
  (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), 
  size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + 
  ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + 
  ggtitle('GO BP') + 
  scale_y_discrete(labels= fct_inorder(
    wrapText(sapply(
      strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),]$Term," \\(GO"),
         `[`, 1),35))) + 
  theme(axis.text=element_text(colour="black"))+
  scale_color_stepsn(colors=rev(magma(256)))
p2
dev.off()

