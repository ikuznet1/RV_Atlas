library(Seurat)
library(hdWGCNA)
library(ggeasy)
library(dplyr)


source('~/Downloads/hdWGCNA_TOM/spatial_functions.R')



#######################################
#############  FIGURE 7A  #############
#######################################
M1 <- readRDS('~/Downloads/hdWGCNA_TOM/Kory_Peds_Hearts/objects/all_data.rds')
M1$Names <- M1$cell.type
M1$NewNames <- M1$Names
M1 <- SetIdent(M1, value = "NewNames")

cluster.ids.new <- c("EC", "FB", "CM", "Myeloid","PC","Endo","NKT","SM","LEC","Adipo","Neuron","B","Mast")
names(cluster.ids.new) <- levels(M1)
M1 <- RenameIdents(M1, cluster.ids.new)
reorder_levels <- c("Adipo","CM","EC","Endo","FB","LEC","Myeloid","Neuron","NKT","PC","SM","Mast","B")
levels(M1) <- reorder_levels
M1$NewNames <- M1@active.ident


#Do some extra doublet removal


feats <- c("PLIN1","RYR2","VWF","LEPR","DCN","CCL21","CSF1R","NRXN1","CD2","PDGFRB","MYH11","KIT","MS4A1")
p<-VlnPlot(M1, features = feats,ncol=1,pt.size=F,group.by="NewNames")



for(i in 1:13) {  
   p[[i]] <- p[[i]] + NoLegend() + easy_remove_axes(which="y",what = c("ticks", "text","line")) + ggtitle("") + ylab(feats[i])
   if(i<13){p[[i]]<-p[[i]]+easy_remove_axes(which="x")}
   
}


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'snPeds_Vln.pdf'), width=4, height=20)
p
dev.off()



M2 <- readRDS('/Volumes/Extreme SSD/Final_Analysis/CellTypes/Post_R3_FINAL_with_counts.rds')


M2 <- merge(M1,M2)
VariableFeatures(M2[["SCT"]]) <- rownames(M2[["SCT"]]@scale.data)
M2$origin <- M2$orig.ident == "SeuratProject"
M2$patient[is.na(M2$patient)] = M2$sample[is.na(M2$patient)]

M2 <- RunPCA(M2)

M2 <- RunHarmony(M2,'patient')

M2 <- M2%>% 
  RunUMAP(reduction = "harmony", dims = 1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:50) %>% 
  FindClusters(resolution=0.5) %>% 
  identity()

M2$CombinedNames <- M2$NewNames
M2$CombinedNames[is.na(M2$NewNames)] <- M2$Names[is.na(M2$NewNames)]
#saveRDS(M2,'~/Downloads/hdWGCNA_TOM/RV_Peds_merge.rds')
#M2 <- readRDS('~/Downloads/hdWGCNA_TOM/RV_Peds_merge.rds')
#M1 <- subset(M2,origin==TRUE)
#saveRDS(M1,'~/Downloads/hdWGCNA_TOM/Peds_clean.rds')

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'sn_RV_Peds_UMAP.pdf'), width=5, height=5)
PlotEmbedding(M2,group.by='CombinedNames',point_size=0.2,plot_under=TRUE,plot_theme=umap_theme()+NoLegend(),raster_dpi=400,raster_scale=0.5)
dev.off()

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'sn_Peds_UMAP_reprojected.pdf'), width=5, height=5)
PlotEmbedding(M1,group.by='NewNames',point_size=0.2,plot_under=TRUE,plot_theme=umap_theme()+NoLegend(),raster_dpi=400,raster_scale=0.5)
dev.off()

M2$condition[M2$condition == "NF"] = "pRV" 
M2$condition[M2$condition == "Donor"] = "NF" 
M2$condition[M2$condition == "SystolicHF"] = "RVF" 



M2$group[is.na(M2$group)] <- M2$condition[is.na(M2$group)]


#######################################
#############  FIGURE 7B  #############
#######################################

t(table(M2$CombinedNames,M2$group))/rowSums(t(table(M2$CombinedNames,M2$group)))


NF_percent_cell <- cbind(as.data.frame(table(subset(M2,group=="NF")@active.ident)/length(subset(M2,group=="NF")@active.ident)*100),type = "NF")
NF_percent_cell$sum <- (rev(cumsum(rev(NF_percent_cell$Freq))) - NF_percent_cell$Freq/2)/100
NF_percent_cell$Freq <- NF_percent_cell$Freq/100


pRV_percent_cell <- cbind(as.data.frame(table(subset(M2,group=="pRV")@active.ident)/length(subset(M2,group=="pRV")@active.ident)*100),type = "pRV")
pRV_percent_cell$sum <- (rev(cumsum(rev(pRV_percent_cell$Freq))) - pRV_percent_cell$Freq/2)/100
pRV_percent_cell$Freq <- pRV_percent_cell$Freq/100


RVF_percent_cell <- cbind(as.data.frame(table(subset(M2,group=="RVF")@active.ident)/length(subset(M2,group=="RVF")@active.ident)*100),type = "RVF")
RVF_percent_cell$sum <- (rev(cumsum(rev(RVF_percent_cell$Freq))) - RVF_percent_cell$Freq/2)/100
RVF_percent_cell$Freq <- RVF_percent_cell$Freq/100

percent_cell_df <- rbind(NF_percent_cell,pRV_percent_cell,RVF_percent_cell)
pdf('~/Downloads/hdWGCNA_TOM/RV_Peds_prev_stacked.pdf',width=6,height=2.5)
ggplot(percent_cell_df, aes(fill=Var1, y=Freq, x=type,label=round(sum,1))) +  
geom_bar(position="stack", stat="identity",width=0.6) + theme_classic() + coord_flip()+
xlab("Disease State") + ylab("Frequency") + labs(fill="Cell type",color='black') + theme(text = element_text(size=20),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),legend.text=element_text(color="black")) + scale_y_continuous(expand=c(0,0)) + geom_label_repel(aes(type,sum,label=scales::percent(round(Freq,2))),fill=NA,nudge_x=0.5,direction="y")
dev.off()

cells <- table(M2$CombinedNames,M2$patient)
cells <- cells[c('CM','FB','EC','Myeloid'),]
cells <- sweep(cells,2,colSums(cells),'/')
cells <- data.frame(cells)

cells$group<-M2$group[match(cells$Var2,M2$patient)]

cells$group<-as.factor(cells$group)

library(ggpubr)
pdf('~/Downloads/hdWGCNA_TOM/Peds_RV_clust_freq.pdf',width=8,height=5)
p <- ggboxplot(cells[length(cells$group):1,],x="group",y="Freq",fill="group",group="group")+
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
  stat_compare_means(aes(group=group),method="anova")
p
dev.off()

cells <- table(M2$NewNames,M2$patient)
cells <- cells[,12:25]
cells <- cells[c('CM','FB','EC','Myeloid'),]
cells <- sweep(cells,2,colSums(cells),'/')
cells <- data.frame(cells)

cells$group<-M2$group[match(cells$Var2,M2$patient)]

cells$group<-as.factor(cells$group)

library(ggpubr)
pdf('~/Downloads/hdWGCNA_TOM/Peds_clust_freq.pdf',width=8,height=5)
p <- ggboxplot(cells[length(cells$group):1,],x="group",y="Freq",fill="group",group="group")+
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
  stat_compare_means(aes(group=group),method="anova")
p
dev.off()

#######################################
#############  FIGURE 7C  #############
#######################################

########Embed bulk in single nuc peds


######Load module
consensus_modules <- read.csv("~/Downloads/hdWGCNA_TOM/bulk_heart_modules.csv")
consensus_modules <- consensus_modules[,1:3]
consensus_modules <- subset(consensus_modules, gene_name %in% rownames(M2))


# remove duplicate gene names
consensus_modules <- consensus_modules[match(unique(consensus_modules$gene_name), consensus_modules$gene_name),]



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
prv_vs_rv_signif <- c('M3','M4','M5','M10','M11','M12','M14')
all_signif <- c('M1','M2','M3','M4','M5','M8','M10','M11','M12','M14','M20','M25','M26','M28')


colnames(MEs)<-paste0('M',match(colnames(MEs),mapping))
M2@meta.data <- cbind(M2@meta.data, MEs)
M2 <- SetIdent(M2, value = "CombinedNames")


#consensus_modules <- read.csv("~/Downloads/hdWGCNA_TOM/bulk_heart_modules.csv")
#consensus_modules <- consensus_modules[,1:3]
#consensus_modules <- subset(consensus_modules, gene_name %in% rownames(M2))
# remove duplicate gene names
#consensus_modules <- consensus_modules[match(unique(consensus_modules$gene_name), consensus_modules$gene_name),]


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
M2$origin[M2$origin == FALSE] = 'RV'

M2$CombinedNamesSplit <- paste0(M2$CombinedNames,'_',M2$origin)

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'RV_Peds_Dot.pdf'), width=7, height=5)
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


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'RV_Peds_Dot_ordered.pdf'), width=7, height=5)

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

M2$groupSplit <- paste0(M2$group,'_',M2$origin)


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'RV_Peds_Dot_disease_ordered.pdf'), width=7, height=4)

p <- DotPlot(M2,paste0('module_',c('M20','M5','M1','M3','M4','M8','M2','M12','M25','M26','M10','M28','M14','M11')),group.by='groupSplit',dot.min=0,col.min=0,col.max=2,idents=c("CM","EC","FB","Myeloid","PC","SM")) +
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
#############  FIGURE 7D  #############
#######################################

M2$group_split <- paste0(M2$group,'_',M2$origin)
M2$group_split <- factor(M2$group_split,levels=c('NF_RV','pRV_RV','RVF_RV','NF_Peds','pRV_Peds','RVF_Peds'))

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'RV_Peds_Dot_CM_Peds.pdf'), width=7, height=3)
p <- DotPlot(subset(M2,CombinedNames=='CM' & origin=='Peds'),paste0('module_',c('M2','M12','M28','M10','M25','M26')),group.by='group_split',dot.min=0,col.min=-2,col.max=2,idents=c("CM","EC","FB","Myeloid","PC","SM")) +
  RotatedAxis() + ylab('')+ xlab('')+
  scale_color_gradient2(high='red', mid='grey95', low='blue') +
  theme(
    panel.border = element_rect(size=1,fill=NA, color='black'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()
) 
p
dev.off()

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'RV_Peds_Dot_CM_RV.pdf'), width=7, height=3)
p <- DotPlot(subset(M2,CombinedNames=='CM' & origin=='RV'),paste0('module_',c('M2','M12','M28','M10','M25','M26')),group.by='group_split',dot.min=0,col.min=-2,col.max=2,idents=c("CM","EC","FB","Myeloid","PC","SM")) +
  RotatedAxis() + ylab('')+ xlab('')+
  scale_color_gradient2(high='red', mid='grey95', low='blue') +
  theme(
    panel.border = element_rect(size=1,fill=NA, color='black'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()
) 
p
dev.off()

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'RV_Peds_Dot_Myeloid_Peds.pdf'), width=7, height=3)
p <- DotPlot(subset(M2,CombinedNames=='Myeloid' & origin=='Peds'),paste0('module_',c('M1','M3','M4','M8')),group.by='group_split',dot.min=0,col.min=-2,col.max=2,idents=c("CM","EC","FB","Myeloid","PC","SM")) +
  RotatedAxis() + ylab('')+ xlab('')+
  scale_color_gradient2(high='red', mid='grey95', low='blue') +
  theme(
    panel.border = element_rect(size=1,fill=NA, color='black'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()
) 
p
dev.off()

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'RV_Peds_Dot_Myeloid_RV.pdf'), width=7, height=3)
p <- DotPlot(subset(M2,CombinedNames=='Myeloid' & origin=='RV'),paste0('module_',c('M1','M3','M4','M8')),group.by='group_split',dot.min=0,col.min=-2,col.max=2,idents=c("CM","EC","FB","Myeloid","PC","SM")) +
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
#############  FIGURE 7E  #############
#######################################

M2 <- SetIdent(M2, value = "groupSplit")
slot(M2$SCT@SCTModel.list[[1]], 'median_umi') = median(M2$SCT@SCTModel.list[[1]]@cell.attributes$umi)

gene_set_RV <- FindMarkers(M2, ident.1 = "RVF_RV", ident.2 = "NF_RV",recorrect_umi=F)
gene_set_Peds <- FindMarkers(M2, ident.1 = "RVF_Peds", ident.2 = "NF_Peds",recorrect_umi=F)




shared <- intersect(rownames(gene_set_RV),rownames(gene_set_Peds))
dataset <- data.frame(Peds=gene_set_Peds[shared,]$avg_log2FC,RV=gene_set_RV[shared,]$avg_log2FC)
rownames(dataset) <- shared


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'Peds_vs_RV.pdf'), width=8, height=8)
ggplot(dataset, aes(x = RV, y=Peds)) + geom_point() + 
  geom_text_repel(label=rownames(dataset),max.overlaps = 50) + theme_classic()
dev.off()

shared <- intersect(rownames(subset(gene_set_RV,pct.1>0.05 & pct.2>0.05)),
  rownames(subset(gene_set_Peds,pct.1>0.05 & pct.2>0.05)))

dataset <- data.frame(Peds=gene_set_Peds[shared,]$avg_log2FC,RV=gene_set_RV[shared,]$avg_log2FC)
rownames(dataset) <- shared


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'Peds_vs_RV_5percent.pdf'), width=8, height=8)
ggplot(dataset, aes(x = RV, y=Peds)) + geom_point() + 
  geom_text_repel(label=rownames(dataset),max.overlaps = 20) + theme_classic()
dev.off()


#######################################
#############  FIGURE 7F  #############
#######################################

M1 <- readRDS('~/Downloads/hdWGCNA_TOM/Kory_Peds_Hearts/objects/cardiomyocyte annotated.rds')
M1$Names <- M1$cell.type
M1$NewNames <- M1$Names
M1$Subnames <- M1$sub.type
M1$NewSubnames <- M1$Subnames
M1 <- SetIdent(M1, value = "NewSubnames")
DefaultAssay(M1) <- 'SCT'


consensus_modules <- read.csv("~/Downloads/hdWGCNA_TOM/bulk_heart_modules.csv")
consensus_modules <- consensus_modules[,1:3]
consensus_modules <- subset(consensus_modules, gene_name %in% rownames(M1))
# remove duplicate gene names
consensus_modules <- consensus_modules[match(unique(consensus_modules$gene_name), consensus_modules$gene_name),]
library(dplyr)
score_calc <- consensus_modules %>% group_by(module) %>% group_split()
module_colors <- unique(unlist(lapply(score_calc,'[[','module')))
module_colors <- paste0('M',match(module_colors,mapping))
M1 <- AddModuleScore(M1,lapply(score_calc,'[[','gene_name'),name="module_score")
cols_current <- colnames(M1@meta.data)
cols_current[startsWith(colnames(M1@meta.data),'module_score')] <- paste0('module_',module_colors)
colnames(M1@meta.data) <- cols_current

Idents(M1) <- factor(x = Idents(M1), levels = sort(levels(M1)))


#Dot Plot of enrichment cell type of CM enriched modules

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'peds_mod_trend_subcluster_CM.pdf'), width=4.5, height=3)

p <- DotPlot(M1,paste0('module_',
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


M1 <- SetIdent(M1, value = "condition")
Idents(M1) <- factor(x = Idents(M1), levels = c('Donor','NF','SystolicHF'))


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'peds_mod_trend_condition_CM.pdf'), width=5, height=2.5)

p <- DotPlot(M1,paste0('module_',
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
#############  FIGURE 7G  #############
#######################################


M1 <- SetIdent(M1, value = "condition")
Idents(M1) <- factor(x = Idents(M1), levels = c('Donor','NF','SystolicHF'))

dbs <-c('GO_Biological_Process_2023','GO_Cellular_Component_2023','GO_Molecular_Function_2023','Reactome_2022', 'ChEA_2022',"LINCS_L1000_Chem_Pert_Consensus_Sigs")

library(enrichR)
#Run enrichment by cell type
combined_set <- data.frame()
combined_output <- data.frame()
bulk_modules <- consensus_modules
bulk_modules$module <- match(consensus_modules$module,mapping)
mods_idx <- list(2,12,10,25,26,28)
cell_types <- unique(M1$NewNames)
comparison <- list(c("SystolicHF","Donor"),c("SystolicHF","NF"))
for (i in mods_idx){
  for (j in cell_types){
    for (k in comparison){
      key_genes <- subset(bulk_modules,module %in% c(i))$gene_name
      key_genes <- key_genes[key_genes %in% rownames(M1)]

      gene_set <- FindMarkers(M1, ident.1 = paste0(k[1]), ident.2 = paste0(k[2]),features=key_genes)
      
      gene_set<-subset(gene_set,p_val_adj<0.05)
      if (length(rownames(gene_set))==0){next}
      gene_set$module <- paste0('M',i)
      gene_set$color <- mapping[i]
      gene_set$comparison <- paste0(k[1],'_',k[2])
      gene_set$celltype <- j

      if (length(combined_set) == 0){
        combined_set <- gene_set
      }
      else {
        combined_set <- rbind(combined_set,gene_set)
      }

      gene_enrich <- subset(gene_set,avg_log2FC<0)
      enriched <- enrichR::enrichr(rownames(gene_enrich), dbs)
      Sys.sleep(5)
      for(db in names(enriched)){
          cur_df <- enriched[[db]]
          if (nrow(cur_df) > 1){
            cur_df$db <- db
            cur_df$module <- paste0('M',i)
            cur_df$celltype <- j
            cur_df$comparison <- paste0(k[1],'_',k[2])
            cur_df$color <- mapping[i]
            cur_df$direction <- 'down'
            combined_output <- rbind(combined_output, cur_df)
          }
      }

      gene_enrich <- subset(gene_set,avg_log2FC>0)
      enriched <- enrichR::enrichr(rownames(gene_enrich), dbs)
      Sys.sleep(5)
      for(db in names(enriched)){
          cur_df <- enriched[[db]]
          if (nrow(cur_df) > 1){
            cur_df$db <- db
            cur_df$module <- paste0('M',i)
            cur_df$celltype <- j
            cur_df$comparison <- paste0(k[1],'_',k[2])
            cur_df$color <- mapping[i]
            cur_df$direction <- 'up'
            combined_output <- rbind(combined_output, cur_df)
          }
      }
    }
  }
}

wrapText <- function(x, len) {
    sapply(x, function(y) paste(strwrap(y, len), collapse = "\n"), USE.NAMES = FALSE)
}


###Systolic HF vs Donor

#Up
selected_terms <- subset(combined_output,db=="GO_Biological_Process_2023")
selected_terms <- subset(selected_terms,direction=="up")
selected_terms <- subset(selected_terms,comparison=="SystolicHF_Donor")
selected_terms <- subset(selected_terms,color %in% mapping[c(2,12)])


# subset selected terms
selected_terms <- subset(selected_terms, P.value < 0.05)
selected_terms$module_celltype <- paste0(selected_terms$module,'_',selected_terms$celltype)
#idx_top_1 <- match(unique(selected_terms$module_celltype),selected_terms$module_celltype)
idx_top_5 <- sort(c(idx_top_1,idx_top_1+1,idx_top_1+2,idx_top_1+3,idx_top_1+4))

selected_terms<-selected_terms[idx_top_5,]


selected_terms$group <- factor(
  as.character(selected_terms$module_celltype),
  levels = unique(selected_terms$module_celltype)
)


# set max pval
quantile(-log(selected_terms$P.value), 0.95)
max_p <- 10

selected_terms$logp <- -log(selected_terms$P.value)
selected_terms$logp <- ifelse(selected_terms$logp > max_p, max_p, selected_terms$logp)

# remove Reactome Term ID
library(stringr)
selected_terms$Term <- str_replace(selected_terms$Term, "\\(GO.*", "")

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

# Reactome Term dot plot
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


mapping <- labels2colors(1:100)


color_df <- data.frame(group = selected_terms$group,colour = selected_terms$color)



color_df$group <- factor(
  as.character(color_df$group),
  levels = unique(color_df$group)
)


c_vect <- color_df$colour
names(c_vect) <- color_df$group



# make the colorbar as its own heatmap
color_df$var <- 1
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


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'peds_CM_by_cluster_terms_cell_type_up_SystolicHF_vs_Donor.pdf'), width=6, height=4)
p / colorbar 
dev.off()



###RVF vs NF

#Down
selected_terms <- subset(combined_output,db=="GO_Biological_Process_2023")
selected_terms <- subset(selected_terms,direction=="down")
selected_terms <- subset(selected_terms,comparison=="SystolicHF_Donor")
selected_terms <- subset(selected_terms,color %in% mapping[c(10,25,26,28)])


# subset selected terms
selected_terms <- subset(selected_terms, P.value < 0.05)
selected_terms$module_celltype <- paste0(selected_terms$module,'_',selected_terms$celltype)
#idx_top_1 <- match(unique(selected_terms$module_celltype),selected_terms$module_celltype)
idx_top_5 <- sort(c(idx_top_1,idx_top_1+1,idx_top_1+2,idx_top_1+3,idx_top_1+4))

selected_terms<-selected_terms[idx_top_5,]


selected_terms$group <- factor(
  as.character(selected_terms$module_celltype),
  levels = unique(selected_terms$module_celltype)
)


# set max pval
quantile(-log(selected_terms$P.value), 0.95)
max_p <- 10

selected_terms$logp <- -log(selected_terms$P.value)
selected_terms$logp <- ifelse(selected_terms$logp > max_p, max_p, selected_terms$logp)

# remove Reactome Term ID
library(stringr)
selected_terms$Term <- str_replace(selected_terms$Term, "\\(GO.*", "")

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

# Reactome Term dot plot
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


mapping <- labels2colors(1:100)


color_df <- data.frame(group = selected_terms$group,colour = selected_terms$color)



color_df$group <- factor(
  as.character(color_df$group),
  levels = unique(color_df$group)
)


c_vect <- color_df$colour
names(c_vect) <- color_df$group



# make the colorbar as its own heatmap
color_df$var <- 1
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


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'peds_CM_by_cluster_terms_cell_type_down_SystolicHF_vs_Donor.pdf'), width=6, height=4)
p / colorbar 
dev.off()



#Both
selected_terms <- subset(combined_output,db=="GO_Biological_Process_2023")
#selected_terms <- subset(selected_terms,direction=="up")
selected_terms <- subset(selected_terms,comparison=="SystolicHF_Donor")
selected_terms <- subset(selected_terms,color %in% mapping[c(2,12,10,25,26,28)])


# subset selected terms
selected_terms <- subset(selected_terms, P.value < 0.05)
selected_terms$module_celltype <- paste0(selected_terms$module,'_',selected_terms$celltype)
idx_top_1 <- match(unique(selected_terms$module_celltype),selected_terms$module_celltype)
idx_top_5 <- sort(c(idx_top_1,idx_top_1+1,idx_top_1+2,idx_top_1+3,idx_top_1+4))

selected_terms<-selected_terms[idx_top_5,]


selected_terms$group <- factor(
  as.character(selected_terms$module_celltype),
  levels = unique(selected_terms$module_celltype)
)


# set max pval
quantile(-log(selected_terms$P.value), 0.95)
max_p <- 10

selected_terms$logp <- -log(selected_terms$P.value)
selected_terms$logp <- ifelse(selected_terms$logp > max_p, max_p, selected_terms$logp)

# remove Reactome Term ID
library(stringr)
selected_terms$Term <- str_replace(selected_terms$Term, "\\(GO.*", "")

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

# Reactome Term dot plot
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


mapping <- labels2colors(1:100)


color_df <- data.frame(group = selected_terms$group,colour = selected_terms$color)



color_df$group <- factor(
  as.character(color_df$group),
  levels = unique(color_df$group)
)


c_vect <- color_df$colour
names(c_vect) <- color_df$group



# make the colorbar as its own heatmap
color_df$var <- 1
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


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'peds_CM_by_cluster_terms_cell_type_both_SystolicHF_vs_Donor.pdf'), width=6, height=9)
p / colorbar 
dev.off()

#######################################
#############  FIGURE 7H  #############
#######################################
#Deep dive M2

bulk_modules <- read.csv("~/Downloads/hdWGCNA_TOM/bulk_heart_modules.csv")
bulk_modules$module <- match(bulk_modules$module,mapping)


combined_set <- data.frame()
mods_idx <- c(2,12,28,10,25,26)
for (i in mods_idx){
  key_genes <- subset(bulk_modules,module %in% c(i))$gene_name
  key_genes <- key_genes[key_genes %in% rownames(M1)]
  gene_set <- FindMarkers(M1, ident.1 = "SystolicHF", ident.2 = "Donor",features=key_genes)
  gene_set<-subset(gene_set,p_val_adj<0.05)
  gene_set$module <- paste0('M',i)
  gene_set$color <- mapping[i]
  if (length(combined_set) == 0){
    combined_set <- gene_set
  }
  else {
    combined_set <- rbind(combined_set,gene_set)
  }
}

#cat(rownames(subset(combined_set,module=='M2' & avg_log2FC>0)),sep='\n')

M2_genes_up <- rownames(subset(combined_set,module=="M2" & avg_log2FC>0))
M10_genes_up <- rownames(subset(combined_set,module=="M10" & avg_log2FC>0))
M26_genes_up <- rownames(subset(combined_set,module=="M26" & avg_log2FC>0))
M25_genes_up <- rownames(subset(combined_set,module=="M25" & avg_log2FC>0))
M28_genes_up <- rownames(subset(combined_set,module=="M28" & avg_log2FC>0))


M2_genes_down <- rownames(subset(combined_set,module=="M2" & avg_log2FC<0))
M10_genes_down <- rownames(subset(combined_set,module=="M10" & avg_log2FC<0))
M26_genes_down <- rownames(subset(combined_set,module=="M26" & avg_log2FC<0))


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

enriched <- enrichr(M2_genes_up, dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('~/Downloads/hdWGCNA_TOM/CM_Peds_M2_enrichr_up.pdf',width=5,height=2.5)
p1<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),], 
  (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), 
  size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + 
  ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + 
  ggtitle('GO Biological Process Up') + 
  scale_y_discrete(labels= fct_inorder(
    wrapText(sapply(
      strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),]$Term," \\(GO"),
         `[`, 1),35))) + 
  theme(axis.text=element_text(colour="black"))+
  scale_color_stepsn(colors=rev(magma(256)))
p1
dev.off()

enriched <- enrichr(M2_genes_down, dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('~/Downloads/hdWGCNA_TOM/CM_Peds_M2_enrichr_down.pdf',width=5,height=2.5)
p2<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),], 
  (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), 
  size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + 
  ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + 
  ggtitle('GO Biological Process Up') + 
  scale_y_discrete(labels= fct_inorder(
    wrapText(sapply(
      strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),]$Term," \\(GO"),
         `[`, 1),35))) + 
  theme(axis.text=element_text(colour="black"))+
  scale_color_stepsn(colors=rev(magma(256)))
p2
dev.off()

pdf('~/Downloads/hdWGCNA_TOM/CM_Peds_M2_enrichr_up_down.pdf',width=6,height=4)
p1/p2
dev.off()



enriched <- enrichr(c(M10_genes_up,M25_genes_up,M26_genes_up,M28_genes_up), dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('~/Downloads/hdWGCNA_TOM/CM_Peds_mito_enrichr_up.pdf',width=6,height=3)
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




###PRV


bulk_modules <- read.csv("~/Downloads/hdWGCNA_TOM/bulk_heart_modules.csv")
bulk_modules$module <- match(bulk_modules$module,mapping)


combined_set <- data.frame()
mods_idx <- c(2,12,28,10,25,26)
for (i in mods_idx){
  key_genes <- subset(bulk_modules,module %in% c(i))$gene_name
  key_genes <- key_genes[key_genes %in% rownames(M1)]
  gene_set <- FindMarkers(M1, ident.1 = "SystolicHF", ident.2 = "NF",features=key_genes)
  gene_set<-subset(gene_set,p_val_adj<0.05)
  gene_set$module <- paste0('M',i)
  gene_set$color <- mapping[i]
  if (length(combined_set) == 0){
    combined_set <- gene_set
  }
  else {
    combined_set <- rbind(combined_set,gene_set)
  }
}

#cat(rownames(subset(combined_set,module=='M2' & avg_log2FC>0)),sep='\n')

M2_genes_up <- rownames(subset(combined_set,module=="M2" & avg_log2FC>0))
M10_genes_up <- rownames(subset(combined_set,module=="M10" & avg_log2FC>0))
M26_genes_up <- rownames(subset(combined_set,module=="M26" & avg_log2FC>0))
M25_genes_up <- rownames(subset(combined_set,module=="M25" & avg_log2FC>0))
M28_genes_up <- rownames(subset(combined_set,module=="M28" & avg_log2FC>0))


M2_genes_down <- rownames(subset(combined_set,module=="M2" & avg_log2FC<0))
M10_genes_down <- rownames(subset(combined_set,module=="M10" & avg_log2FC<0))
M26_genes_down <- rownames(subset(combined_set,module=="M26" & avg_log2FC<0))


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

enriched <- enrichr(M2_genes_up, dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('~/Downloads/hdWGCNA_TOM/CM_Peds_M2_enrichr_up_RVF_vs_pRV.pdf.pdf',width=5,height=2.5)
p1<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),], 
  (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), 
  size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + 
  ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + 
  ggtitle('GO Biological Process Up') + 
  scale_y_discrete(labels= fct_inorder(
    wrapText(sapply(
      strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),]$Term," \\(GO"),
         `[`, 1),35))) + 
  theme(axis.text=element_text(colour="black"))+
  scale_color_stepsn(colors=rev(magma(256)))
p1
dev.off()

enriched <- enrichr(M2_genes_down, dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('~/Downloads/hdWGCNA_TOM/CM_Peds_M2_enrichr_down_RVF_vs_pRV.pdf.pdf',width=5,height=2.5)
p2<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),], 
  (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), 
  size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + 
  ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + 
  ggtitle('GO Biological Process Up') + 
  scale_y_discrete(labels= fct_inorder(
    wrapText(sapply(
      strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),]$Term," \\(GO"),
         `[`, 1),35))) + 
  theme(axis.text=element_text(colour="black"))+
  scale_color_stepsn(colors=rev(magma(256)))
p2
dev.off()

pdf('~/Downloads/hdWGCNA_TOM/CM_Peds_M2_enrichr_up_down_RVF_vs_pRV.pdf',width=6,height=4)
p1/p2
dev.off()



enriched <- enrichr(c(M10_genes_up,M25_genes_up,M26_genes_up,M28_genes_up), dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('~/Downloads/hdWGCNA_TOM/CM_Peds_mito_enrichr_up_RVF_vs_pRV.pdf.pdf',width=6,height=3)
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



###SV vs NF


bulk_modules <- read.csv("~/Downloads/hdWGCNA_TOM/bulk_heart_modules.csv")
bulk_modules$module <- match(bulk_modules$module,mapping)


combined_set <- data.frame()
mods_idx <- c(2,12,28,10,25,26)
for (i in mods_idx){
  key_genes <- subset(bulk_modules,module %in% c(i))$gene_name
  key_genes <- key_genes[key_genes %in% rownames(M1)]
  gene_set <- FindMarkers(M1, ident.1 = "NF", ident.2 = "Donor",features=key_genes)
  gene_set<-subset(gene_set,p_val_adj<0.05)
  gene_set$module <- paste0('M',i)
  gene_set$color <- mapping[i]
  if (length(combined_set) == 0){
    combined_set <- gene_set
  }
  else {
    combined_set <- rbind(combined_set,gene_set)
  }
}

#cat(rownames(subset(combined_set,module=='M2' & avg_log2FC>0)),sep='\n')

M2_genes_up <- rownames(subset(combined_set,module=="M2" & avg_log2FC>0))
M10_genes_up <- rownames(subset(combined_set,module=="M10" & avg_log2FC>0))
M26_genes_up <- rownames(subset(combined_set,module=="M26" & avg_log2FC>0))
M25_genes_up <- rownames(subset(combined_set,module=="M25" & avg_log2FC>0))
M28_genes_up <- rownames(subset(combined_set,module=="M28" & avg_log2FC>0))


M2_genes_down <- rownames(subset(combined_set,module=="M2" & avg_log2FC<0))
M10_genes_down <- rownames(subset(combined_set,module=="M10" & avg_log2FC<0))
M26_genes_down <- rownames(subset(combined_set,module=="M26" & avg_log2FC<0))


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

enriched <- enrichr(M2_genes_up, dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('~/Downloads/hdWGCNA_TOM/CM_Peds_M2_enrichr_up_pRV_vs_NF.pdf.pdf',width=5,height=2.5)
p1<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),], 
  (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), 
  size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + 
  ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + 
  ggtitle('GO Biological Process Up') + 
  scale_y_discrete(labels= fct_inorder(
    wrapText(sapply(
      strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),]$Term," \\(GO"),
         `[`, 1),35))) + 
  theme(axis.text=element_text(colour="black"))+
  scale_color_stepsn(colors=rev(magma(256)))
p1
dev.off()

enriched <- enrichr(M2_genes_down, dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('~/Downloads/hdWGCNA_TOM/CM_Peds_M2_enrichr_down_pRV_vs_NF.pdf.pdf',width=5,height=2.5)
p2<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),], 
  (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), 
  size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + 
  ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + 
  ggtitle('GO Biological Process Up') + 
  scale_y_discrete(labels= fct_inorder(
    wrapText(sapply(
      strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),]$Term," \\(GO"),
         `[`, 1),35))) + 
  theme(axis.text=element_text(colour="black"))+
  scale_color_stepsn(colors=rev(magma(256)))
p2
dev.off()

pdf('~/Downloads/hdWGCNA_TOM/CM_Peds_M2_enrichr_up_down_pRV_vs_NF.pdf',width=6,height=4)
p1/p2
dev.off()



enriched <- enrichr(c(M10_genes_up,M25_genes_up,M26_genes_up,M28_genes_up), dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('~/Downloads/hdWGCNA_TOM/CM_Peds_mito_enrichr_up_pRV_vs_NF.pdf.pdf',width=6,height=3)
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
#############  FIGURE 7I  #############
#######################################

M1 <- readRDS('~/Downloads/hdWGCNA_TOM/Kory_Peds_Hearts/objects/myeloid annotated.rds')
M1$Names <- M1$cell.type
M1$NewNames <- M1$Names
M1$Subnames <- M1$sub.type
M1$NewSubnames <- M1$Subnames
M1 <- SetIdent(M1, value = "NewSubnames")
DefaultAssay(M1) <- 'SCT'


M1 <- AddModuleScore(M1, features=list(c("TREM2","GPNMB","MITF","SPP1")),assay="SCT",name="TREM2_Mac_Score")
M1 <- AddModuleScore(M1, features=list(c("CLEC9A","ZBTB46","CD1C","CD226")),assay="SCT",name="DC_Score")
M1 <- AddModuleScore(M1, features=list(c("FCN1","LILRB2","ITGAL","CSF3R")),assay="SCT",name="Mono_Score")
M1 <- AddModuleScore(M1, features=list(c("CCR2","CX3CR1","ITGAX")),assay="SCT",name="CCR2+_rMac_Score")
M1 <- AddModuleScore(M1, features=list(c("LYVE1","FOLR2","SIGLEC1","F13A1")),assay="SCT",name="CCR2-_rMac1_Score")
M1 <- AddModuleScore(M1, features=list(c("RBMS3","PLA2G5","EBF1")),assay="SCT",name="CCR2-_rMac2_Score")
M1 <- AddModuleScore(M1, features=list(c("IL1B","CCL3","CCL4","CXCL3","CXCL8")),assay="SCT",name="iMac_Score")

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'Peds_myeloid_dot.pdf'), width=7, height=6)
DotPlot(M1, c("Mono_Score1","CCR2-_rMac2_Score1","iMac_Score1","TREM2_Mac_Score1",
  "CCR2+_rMac_Score1","CCR2-_rMac1_Score1","DC_Score1"))+xlab('Marker Score')+
  scale_x_discrete(labels=c('Mono','CCR2- 2','iMac','TREM2','CCR2+','CCR2- 1','DC'))
dev.off()

M1 <- SetIdent(M1, value = "condition")
DotPlot(M1, c("Mono_Score1","CCR2-_rMac2_Score1","iMac_Score1","TREM2_Mac_Score1",
  "CCR2+_rMac_Score1","CCR2-_rMac1_Score1","DC_Score1"))+xlab('Marker Score')+
  scale_x_discrete(labels=c('Mono','CCR2- 2','iMac','TREM2','CCR2+','CCR2- 1','DC'))

a <- FindMarkers(M1,ident.1='SystolicHF',ident.2='NF')

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'Peds_myeloid_volcano_RVF_pRV.pdf'), width=6, height=5)
p1<-EnhancedVolcano(a,lab=rownames(a),
  x='avg_log2FC',y='p_val_adj',
  FCcutoff = 0.1) + coord_flip()
dev.off()

a <- FindMarkers(M1,ident.1='NF',ident.2='Donor')

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'Peds_myeloid_volcano_pRV_NF.pdf'), width=6, height=5)
p2<-EnhancedVolcano(a,lab=rownames(a),
  x='avg_log2FC',y='p_val_adj',
  FCcutoff = 0.1) + coord_flip()
dev.off()

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'Peds_myeloid_volcano.pdf'), width=12, height=5)
p2+p1
dev.off()

#######################################
#############  FIGURE 7J  #############
#######################################

consensus_modules <- read.csv("~/Downloads/hdWGCNA_TOM/bulk_heart_modules.csv")
consensus_modules <- consensus_modules[,1:3]
consensus_modules <- subset(consensus_modules, gene_name %in% rownames(M1))
# remove duplicate gene names
consensus_modules <- consensus_modules[match(unique(consensus_modules$gene_name), consensus_modules$gene_name),]
library(dplyr)
score_calc <- consensus_modules %>% group_by(module) %>% group_split()
module_colors <- unique(unlist(lapply(score_calc,'[[','module')))
module_colors <- paste0('M',match(module_colors,mapping))
M1 <- AddModuleScore(M1,lapply(score_calc,'[[','gene_name'),name="module_score")
cols_current <- colnames(M1@meta.data)
cols_current[startsWith(colnames(M1@meta.data),'module_score')] <- paste0('module_',module_colors)
colnames(M1@meta.data) <- cols_current


Idents(M1) <- factor(x = Idents(M1), levels = sort(levels(M1)))


#Dot Plot of enrichment cell type of CM enriched modules

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'peds_mod_trend_subcluster_myeloid.pdf'), width=5, height=3)

p <- DotPlot(M1,paste0('module_',
  c('M1','M3','M4','M8')),dot.min=0,col.min=0,col.max=2) +
  RotatedAxis() + ylab('')+ xlab('')+
  scale_color_gradient2(high='red', mid='grey95', low='blue') +
  theme(
    panel.border = element_rect(size=1,fill=NA, color='black'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()
) 
p
dev.off()


M1 <- SetIdent(M1, value = "condition")
Idents(M1) <- factor(x = Idents(M1), levels = c('Donor','NF','SystolicHF'))


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'peds_mod_trend_condition_myeloid.pdf'), width=5, height=2.5)

p <- DotPlot(M1,paste0('module_',
  c('M1','M3','M4','M8')),dot.min=0,col.min=0,col.max=2) +
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
#############  FIGURE 7K  #############
#######################################


bulk_modules <- read.csv("~/Downloads/hdWGCNA_TOM/bulk_heart_modules.csv")
bulk_modules$module <- match(bulk_modules$module,mapping)


combined_set <- data.frame()
mods_idx <- c(1,3,4,8)
for (i in mods_idx){
  key_genes <- subset(bulk_modules,module %in% c(i))$gene_name
  key_genes <- key_genes[key_genes %in% rownames(M1)]
  gene_set <- FindMarkers(M1, ident.1 = "SystolicHF", ident.2 = "Donor",features=key_genes)
  gene_set<-subset(gene_set,p_val_adj<0.05)
  gene_set$module <- paste0('M',i)
  gene_set$color <- mapping[i]
  if (length(combined_set) == 0){
    combined_set <- gene_set
  }
  else {
    combined_set <- rbind(combined_set,gene_set)
  }
}

#cat(rownames(subset(combined_set,module=='M2' & avg_log2FC>0)),sep='\n')

M1_genes_up <- rownames(subset(combined_set,module=="M1" & avg_log2FC>0))
M8_genes_up <- rownames(subset(combined_set,module=="M8" & avg_log2FC>0))



M1_genes_down <- rownames(subset(combined_set,module=="M1" & avg_log2FC<0))
M8_genes_down <- rownames(subset(combined_set,module=="M8" & avg_log2FC<0))


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

enriched <- enrichr(M1_genes_up, dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('~/Downloads/hdWGCNA_TOM/myeloid_Peds_M1_enrichr_up.pdf',width=5,height=2.5)
p1<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),], 
  (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), 
  size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + 
  ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + 
  ggtitle('GO Biological Process Up') + 
  scale_y_discrete(labels= fct_inorder(
    wrapText(sapply(
      strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),]$Term," \\(GO"),
         `[`, 1),35))) + 
  theme(axis.text=element_text(colour="black"))+
  scale_color_stepsn(colors=rev(magma(256)))
p1
dev.off()

enriched <- enrichr(M1_genes_down, dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('~/Downloads/hdWGCNA_TOM/myeloid_Peds_M1_enrichr_down.pdf',width=5,height=2.5)
p2<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),], 
  (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), 
  size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + 
  ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + 
  ggtitle('GO Biological Process Up') + 
  scale_y_discrete(labels= fct_inorder(
    wrapText(sapply(
      strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),]$Term," \\(GO"),
         `[`, 1),35))) + 
  theme(axis.text=element_text(colour="black"))+
  scale_color_stepsn(colors=rev(magma(256)))
p2
dev.off()

pdf('~/Downloads/hdWGCNA_TOM/myeloid_Peds_M1_enrichr_up_down.pdf',width=6,height=4)
p1/p2
dev.off()



enriched <- enrichr(M8_genes_up, dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('~/Downloads/hdWGCNA_TOM/myeloid_Peds_M8_enrichr_up.pdf',width=6,height=3)
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

enriched <- enrichr(M8_genes_down, dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('~/Downloads/hdWGCNA_TOM/myeloid_Peds_M8_enrichr_down.pdf',width=6,height=3)
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

pdf('~/Downloads/hdWGCNA_TOM/myeloid_Peds_M8_enrichr_up_down.pdf',width=6,height=4)
p1/p2
dev.off()









###PRV


bulk_modules <- read.csv("~/Downloads/hdWGCNA_TOM/bulk_heart_modules.csv")
bulk_modules$module <- match(bulk_modules$module,mapping)


combined_set <- data.frame()
mods_idx <- c(1,3,4,8)
for (i in mods_idx){
  key_genes <- subset(bulk_modules,module %in% c(i))$gene_name
  key_genes <- key_genes[key_genes %in% rownames(M1)]
  gene_set <- FindMarkers(M1, ident.1 = "SystolicHF", ident.2 = "NF",features=key_genes)
  gene_set<-subset(gene_set,p_val_adj<0.05)
  gene_set$module <- paste0('M',i)
  gene_set$color <- mapping[i]
  if (length(combined_set) == 0){
    combined_set <- gene_set
  }
  else {
    combined_set <- rbind(combined_set,gene_set)
  }
}

#cat(rownames(subset(combined_set,module=='M2' & avg_log2FC>0)),sep='\n')

M1_genes_up <- rownames(subset(combined_set,module=="M1" & avg_log2FC>0))
M8_genes_up <- rownames(subset(combined_set,module=="M8" & avg_log2FC>0))



M1_genes_down <- rownames(subset(combined_set,module=="M1" & avg_log2FC<0))
M8_genes_down <- rownames(subset(combined_set,module=="M8" & avg_log2FC<0))


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

enriched <- enrichr(M1_genes_up, dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('~/Downloads/hdWGCNA_TOM/myeloid_Peds_M1_enrichr_up_RVF_vs_pRV.pdf.pdf',width=5,height=2.5)
p1<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),], 
  (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), 
  size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + 
  ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + 
  ggtitle('GO Biological Process Up') + 
  scale_y_discrete(labels= fct_inorder(
    wrapText(sapply(
      strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),]$Term," \\(GO"),
         `[`, 1),35))) + 
  theme(axis.text=element_text(colour="black"))+
  scale_color_stepsn(colors=rev(magma(256)))
p1
dev.off()

enriched <- enrichr(M1_genes_down, dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('~/Downloads/hdWGCNA_TOM/myeloid_Peds_M1_enrichr_down_RVF_vs_pRV.pdf.pdf',width=5,height=2.5)
p2<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),], 
  (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), 
  size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + 
  ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + 
  ggtitle('GO Biological Process Up') + 
  scale_y_discrete(labels= fct_inorder(
    wrapText(sapply(
      strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),]$Term," \\(GO"),
         `[`, 1),35))) + 
  theme(axis.text=element_text(colour="black"))+
  scale_color_stepsn(colors=rev(magma(256)))
p2
dev.off()

pdf('~/Downloads/hdWGCNA_TOM/myeloid_Peds_M1_enrichr_up_down_RVF_vs_pRV.pdf',width=6,height=4)
p1/p2
dev.off()


enriched <- enrichr(M8_genes_up, dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('~/Downloads/hdWGCNA_TOM/myeloid_Peds_M8_enrichr_up_RVF_vs_pRV.pdf.pdf',width=5,height=2.5)
p1<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),], 
  (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), 
  size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + 
  ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + 
  ggtitle('GO Biological Process Up') + 
  scale_y_discrete(labels= fct_inorder(
    wrapText(sapply(
      strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),]$Term," \\(GO"),
         `[`, 1),35))) + 
  theme(axis.text=element_text(colour="black"))+
  scale_color_stepsn(colors=rev(magma(256)))
p1
dev.off()

enriched <- enrichr(M8_genes_down, dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('~/Downloads/hdWGCNA_TOM/myeloid_Peds_M8_enrichr_down_RVF_vs_pRV.pdf.pdf',width=5,height=2.5)
p2<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),], 
  (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), 
  size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + 
  ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + 
  ggtitle('GO Biological Process Up') + 
  scale_y_discrete(labels= fct_inorder(
    wrapText(sapply(
      strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),]$Term," \\(GO"),
         `[`, 1),35))) + 
  theme(axis.text=element_text(colour="black"))+
  scale_color_stepsn(colors=rev(magma(256)))
p2
dev.off()

pdf('~/Downloads/hdWGCNA_TOM/myeloid_Peds_M8_enrichr_up_down_RVF_vs_pRV.pdf',width=6,height=4)
p1/p2
dev.off()




###SV vs NF

bulk_modules <- read.csv("~/Downloads/hdWGCNA_TOM/bulk_heart_modules.csv")
bulk_modules$module <- match(bulk_modules$module,mapping)


combined_set <- data.frame()
mods_idx <- c(1,3,4,8)
for (i in mods_idx){
  key_genes <- subset(bulk_modules,module %in% c(i))$gene_name
  key_genes <- key_genes[key_genes %in% rownames(M1)]
  gene_set <- FindMarkers(M1, ident.1 = "NF", ident.2 = "Donor",features=key_genes)
  gene_set<-subset(gene_set,p_val_adj<0.05)
  gene_set$module <- paste0('M',i)
  gene_set$color <- mapping[i]
  if (length(combined_set) == 0){
    combined_set <- gene_set
  }
  else {
    combined_set <- rbind(combined_set,gene_set)
  }
}

#cat(rownames(subset(combined_set,module=='M2' & avg_log2FC>0)),sep='\n')

M1_genes_up <- rownames(subset(combined_set,module=="M1" & avg_log2FC>0))
M8_genes_up <- rownames(subset(combined_set,module=="M8" & avg_log2FC>0))



M1_genes_down <- rownames(subset(combined_set,module=="M1" & avg_log2FC<0))
M8_genes_down <- rownames(subset(combined_set,module=="M8" & avg_log2FC<0))


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

enriched <- enrichr(M1_genes_up, dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('~/Downloads/hdWGCNA_TOM/myeloid_Peds_M1_enrichr_up_pRV_vs_NF.pdf.pdf',width=5,height=2.5)
p1<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),], 
  (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), 
  size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + 
  ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + 
  ggtitle('GO Biological Process Up') + 
  scale_y_discrete(labels= fct_inorder(
    wrapText(sapply(
      strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),]$Term," \\(GO"),
         `[`, 1),35))) + 
  theme(axis.text=element_text(colour="black"))+
  scale_color_stepsn(colors=rev(magma(256)))
p1
dev.off()

enriched <- enrichr(M1_genes_down, dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('~/Downloads/hdWGCNA_TOM/myeloid_Peds_M1_enrichr_down_pRV_vs_NFpdf.pdf',width=5,height=2.5)
p2<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),], 
  (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), 
  size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + 
  ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + 
  ggtitle('GO Biological Process Up') + 
  scale_y_discrete(labels= fct_inorder(
    wrapText(sapply(
      strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),]$Term," \\(GO"),
         `[`, 1),35))) + 
  theme(axis.text=element_text(colour="black"))+
  scale_color_stepsn(colors=rev(magma(256)))
p2
dev.off()

pdf('~/Downloads/hdWGCNA_TOM/myeloid_Peds_M1_enrichr_up_down_pRV_vs_NF.pdf',width=6,height=4)
p1/p2
dev.off()


enriched <- enrichr(M8_genes_up, dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('~/Downloads/hdWGCNA_TOM/myeloid_Peds_M8_enrichr_up_pRV_vs_NF.pdf.pdf',width=5,height=2.5)
p1<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),], 
  (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), 
  size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + 
  ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + 
  ggtitle('GO Biological Process Up') + 
  scale_y_discrete(labels= fct_inorder(
    wrapText(sapply(
      strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),]$Term," \\(GO"),
         `[`, 1),35))) + 
  theme(axis.text=element_text(colour="black"))+
  scale_color_stepsn(colors=rev(magma(256)))
p1
dev.off()

enriched <- enrichr(M8_genes_down, dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('~/Downloads/hdWGCNA_TOM/myeloid_Peds_M8_enrichr_down_pRV_vs_NF.pdf.pdf',width=5,height=2.5)
p2<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),], 
  (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), 
  size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + 
  ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + 
  ggtitle('GO Biological Process Up') + 
  scale_y_discrete(labels= fct_inorder(
    wrapText(sapply(
      strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),]$Term," \\(GO"),
         `[`, 1),35))) + 
  theme(axis.text=element_text(colour="black"))+
  scale_color_stepsn(colors=rev(magma(256)))
p2
dev.off()

pdf('~/Downloads/hdWGCNA_TOM/myeloid_Peds_M8_enrichr_up_down_pRV_vs_NF.pdf',width=6,height=4)
p1/p2
dev.off()



#######################################
#############  FIGURE 7L  #############
#######################################

M1$Names_group <- paste0(M1$sub.type,'_',M1$condition)
#Run enrichment by cell type
Idents(M1) <- "Names_group"
combined_set <- data.frame()
combined_output <- data.frame()

mods_idx <- c(1,3,4,8)
cell_types <- unique(M1$sub.type)
comparison <- list(c("SystolicHF","Donor"),c("SystolicHF","NF"),c("NF","Donor"))
for (i in mods_idx){
  for (j in cell_types){
    for (k in comparison){
      key_genes <- subset(bulk_modules,module %in% c(i))$gene_name
      key_genes <- key_genes[key_genes %in% rownames(M1)]

      gene_set <- FindMarkers(M1, ident.1 = paste0(j,"_",k[1]), ident.2 = paste0(j,"_",k[2]),features=key_genes)
      
      gene_set<-subset(gene_set,p_val_adj<0.05)
      if (length(rownames(gene_set))==0){next}
      gene_set$module <- paste0('M',i)
      gene_set$color <- mapping[i]
      gene_set$comparison <- paste0(k[1],'_',k[2])
      gene_set$celltype <- j

      if (length(combined_set) == 0){
        combined_set <- gene_set
      }
      else {
        combined_set <- rbind(combined_set,gene_set)
      }

      gene_enrich <- subset(gene_set,avg_log2FC<0)
      enriched <- enrichR::enrichr(rownames(gene_enrich), dbs)
      Sys.sleep(5)
      for(db in names(enriched)){
          cur_df <- enriched[[db]]
          if (nrow(cur_df) > 1){
            cur_df$db <- db
            cur_df$module <- paste0('M',i)
            cur_df$celltype <- j
            cur_df$comparison <- paste0(k[1],'_',k[2])
            cur_df$color <- mapping[i]
            cur_df$direction <- 'down'
            combined_output <- rbind(combined_output, cur_df)
          }
      }

      gene_enrich <- subset(gene_set,avg_log2FC>0)
      enriched <- enrichR::enrichr(rownames(gene_enrich), dbs)
      Sys.sleep(5)
      for(db in names(enriched)){
          cur_df <- enriched[[db]]
          if (nrow(cur_df) > 1){
            cur_df$db <- db
            cur_df$module <- paste0('M',i)
            cur_df$celltype <- j
            cur_df$comparison <- paste0(k[1],'_',k[2])
            cur_df$color <- mapping[i]
            cur_df$direction <- 'up'
            combined_output <- rbind(combined_output, cur_df)
          }
      }
    }
  }
}



###RVF vs NF

#Up
selected_terms <- subset(combined_output,db=="Reactome_2016")
selected_terms <- subset(selected_terms,direction=="up")
selected_terms <- subset(selected_terms,comparison=="SystolicHF_Donor")
selected_terms <- subset(selected_terms,color %in% mapping[c(1,3,4,8)])


# subset selected terms
selected_terms <- subset(selected_terms, P.value < 0.05)
selected_terms$module_celltype <- paste0(selected_terms$module,'_',selected_terms$celltype)
idx_top_1 <- match(unique(selected_terms$module_celltype),selected_terms$module_celltype)
#idx_top_5 <- sort(c(idx_top_1,idx_top_1+1,idx_top_1+2,idx_top_1+3,idx_top_1+4))

selected_terms<-selected_terms[idx_top_1,]


selected_terms$group <- factor(
  as.character(selected_terms$module_celltype),
  levels = unique(selected_terms$module_celltype)
)


# set max pval
quantile(-log(selected_terms$P.value), 0.95)
max_p <- 10

selected_terms$logp <- -log(selected_terms$P.value)
selected_terms$logp <- ifelse(selected_terms$logp > max_p, max_p, selected_terms$logp)

# remove Reactome Term ID
library(stringr)
selected_terms$Term <- str_replace(selected_terms$Term, "\\ R-HSA.*", "")

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

# Reactome Term dot plot
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


mapping <- labels2colors(1:100)


color_df <- data.frame(group = selected_terms$group,colour = selected_terms$color)



color_df$group <- factor(
  as.character(color_df$group),
  levels = unique(color_df$group)
)


c_vect <- color_df$colour
names(c_vect) <- color_df$group



# make the colorbar as its own heatmap
color_df$var <- 1
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


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'peds_Myeloid_reactome_terms_cell_type_up_RVF_vs_NF.pdf'), width=8, height=7)
p / colorbar 
dev.off()

#Down
selected_terms <- subset(combined_output,db=="Reactome_2016")
selected_terms <- subset(selected_terms,direction=="down")
selected_terms <- subset(selected_terms,comparison=="SystolicHF_Donor")

selected_terms <- subset(selected_terms,color %in% mapping[c(1,3,4,8)])


# subset selected terms
selected_terms <- subset(selected_terms, P.value < 0.05)
selected_terms$module_celltype <- paste0(selected_terms$module,'_',selected_terms$celltype)
idx_top_1 <- match(unique(selected_terms$module_celltype),selected_terms$module_celltype)
#idx_top_5 <- sort(c(idx_top_1,idx_top_1+1,idx_top_1+2,idx_top_1+3,idx_top_1+4))

selected_terms<-selected_terms[idx_top_1,]


selected_terms$group <- factor(
  as.character(selected_terms$module_celltype),
  levels = unique(selected_terms$module_celltype)
)


# set max pval
quantile(-log(selected_terms$P.value), 0.95)
max_p <- 10

selected_terms$logp <- -log(selected_terms$P.value)
selected_terms$logp <- ifelse(selected_terms$logp > max_p, max_p, selected_terms$logp)

# remove Reactome Term ID
library(stringr)
selected_terms$Term <- str_replace(selected_terms$Term, "\\ R-HSA.*", "")

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

# Reactome Term dot plot
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


mapping <- labels2colors(1:100)


color_df <- data.frame(group = selected_terms$group,colour = selected_terms$color)



color_df$group <- factor(
  as.character(color_df$group),
  levels = unique(color_df$group)
)


c_vect <- color_df$colour
names(c_vect) <- color_df$group



# make the colorbar as its own heatmap
color_df$var <- 1
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


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'peds_Myeloid_reactome_terms_cell_type_down_RVF_vs_NF.pdf'), width=6.6, height=5)
p / colorbar 
dev.off()


###RVF vs pRV

#Up
selected_terms <- subset(combined_output,db=="Reactome_2016")
selected_terms <- subset(selected_terms,direction=="up")
selected_terms <- subset(selected_terms,comparison=="SystolicHF_NF")
selected_terms <- subset(selected_terms,color %in% mapping[c(1,3,4,8)])


# subset selected terms
selected_terms <- subset(selected_terms, P.value < 0.05)
selected_terms$module_celltype <- paste0(selected_terms$module,'_',selected_terms$celltype)
idx_top_1 <- match(unique(selected_terms$module_celltype),selected_terms$module_celltype)
#idx_top_5 <- sort(c(idx_top_1,idx_top_1+1,idx_top_1+2,idx_top_1+3,idx_top_1+4))

selected_terms<-selected_terms[idx_top_1,]


selected_terms$group <- factor(
  as.character(selected_terms$module_celltype),
  levels = unique(selected_terms$module_celltype)
)


# set max pval
quantile(-log(selected_terms$P.value), 0.95)
max_p <- 10

selected_terms$logp <- -log(selected_terms$P.value)
selected_terms$logp <- ifelse(selected_terms$logp > max_p, max_p, selected_terms$logp)

# remove Reactome Term ID
library(stringr)
selected_terms$Term <- str_replace(selected_terms$Term, "\\ R-HSA.*", "")

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

# Reactome Term dot plot
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


mapping <- labels2colors(1:100)


color_df <- data.frame(group = selected_terms$group,colour = selected_terms$color)



color_df$group <- factor(
  as.character(color_df$group),
  levels = unique(color_df$group)
)


c_vect <- color_df$colour
names(c_vect) <- color_df$group



# make the colorbar as its own heatmap
color_df$var <- 1
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


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'peds_Myeloid_reactome_terms_cell_type_up_RVF_vs_pRV.pdf'), width=8, height=7)
p / colorbar 
dev.off()

#Down
selected_terms <- subset(combined_output,db=="Reactome_2016")
selected_terms <- subset(selected_terms,direction=="down")
selected_terms <- subset(selected_terms,comparison=="SystolicHF_NF")

selected_terms <- subset(selected_terms,color %in% mapping[c(1,3,4,8)])


# subset selected terms
selected_terms <- subset(selected_terms, P.value < 0.05)
selected_terms$module_celltype <- paste0(selected_terms$module,'_',selected_terms$celltype)
idx_top_1 <- match(unique(selected_terms$module_celltype),selected_terms$module_celltype)
#idx_top_5 <- sort(c(idx_top_1,idx_top_1+1,idx_top_1+2,idx_top_1+3,idx_top_1+4))

selected_terms<-selected_terms[idx_top_1,]


selected_terms$group <- factor(
  as.character(selected_terms$module_celltype),
  levels = unique(selected_terms$module_celltype)
)


# set max pval
quantile(-log(selected_terms$P.value), 0.95)
max_p <- 10

selected_terms$logp <- -log(selected_terms$P.value)
selected_terms$logp <- ifelse(selected_terms$logp > max_p, max_p, selected_terms$logp)

# remove Reactome Term ID
library(stringr)
selected_terms$Term <- str_replace(selected_terms$Term, "\\ R-HSA.*", "")

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

# Reactome Term dot plot
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


mapping <- labels2colors(1:100)


color_df <- data.frame(group = selected_terms$group,colour = selected_terms$color)



color_df$group <- factor(
  as.character(color_df$group),
  levels = unique(color_df$group)
)


c_vect <- color_df$colour
names(c_vect) <- color_df$group



# make the colorbar as its own heatmap
color_df$var <- 1
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


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'peds_Myeloid_reactome_terms_cell_type_down_RVF_vs_pRV.pdf'), width=6.6, height=5)
p / colorbar 
dev.off()


###RVF vs pRV

#Up
selected_terms <- subset(combined_output,db=="Reactome_2016")
selected_terms <- subset(selected_terms,direction=="up")
selected_terms <- subset(selected_terms,comparison=="NF_Donor")
selected_terms <- subset(selected_terms,color %in% mapping[c(1,3,4,8)])


# subset selected terms
selected_terms <- subset(selected_terms, P.value < 0.05)
selected_terms$module_celltype <- paste0(selected_terms$module,'_',selected_terms$celltype)
idx_top_1 <- match(unique(selected_terms$module_celltype),selected_terms$module_celltype)
#idx_top_5 <- sort(c(idx_top_1,idx_top_1+1,idx_top_1+2,idx_top_1+3,idx_top_1+4))

selected_terms<-selected_terms[idx_top_1,]


selected_terms$group <- factor(
  as.character(selected_terms$module_celltype),
  levels = unique(selected_terms$module_celltype)
)


# set max pval
quantile(-log(selected_terms$P.value), 0.95)
max_p <- 10

selected_terms$logp <- -log(selected_terms$P.value)
selected_terms$logp <- ifelse(selected_terms$logp > max_p, max_p, selected_terms$logp)

# remove Reactome Term ID
library(stringr)
selected_terms$Term <- str_replace(selected_terms$Term, "\\ R-HSA.*", "")

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

# Reactome Term dot plot
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


mapping <- labels2colors(1:100)


color_df <- data.frame(group = selected_terms$group,colour = selected_terms$color)



color_df$group <- factor(
  as.character(color_df$group),
  levels = unique(color_df$group)
)


c_vect <- color_df$colour
names(c_vect) <- color_df$group



# make the colorbar as its own heatmap
color_df$var <- 1
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


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'peds_Myeloid_reactome_terms_cell_type_up_pRV_vs_NF.pdf'), width=8, height=7)
p / colorbar 
dev.off()

#Down
selected_terms <- subset(combined_output,db=="Reactome_2016")
selected_terms <- subset(selected_terms,direction=="down")
selected_terms <- subset(selected_terms,comparison=="NF_Donor")

selected_terms <- subset(selected_terms,color %in% mapping[c(1,3,4,8)])


# subset selected terms
selected_terms <- subset(selected_terms, P.value < 0.05)
selected_terms$module_celltype <- paste0(selected_terms$module,'_',selected_terms$celltype)
idx_top_1 <- match(unique(selected_terms$module_celltype),selected_terms$module_celltype)
#idx_top_5 <- sort(c(idx_top_1,idx_top_1+1,idx_top_1+2,idx_top_1+3,idx_top_1+4))

selected_terms<-selected_terms[idx_top_1,]


selected_terms$group <- factor(
  as.character(selected_terms$module_celltype),
  levels = unique(selected_terms$module_celltype)
)


# set max pval
quantile(-log(selected_terms$P.value), 0.95)
max_p <- 10

selected_terms$logp <- -log(selected_terms$P.value)
selected_terms$logp <- ifelse(selected_terms$logp > max_p, max_p, selected_terms$logp)

# remove Reactome Term ID
library(stringr)
selected_terms$Term <- str_replace(selected_terms$Term, "\\ R-HSA.*", "")

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

# Reactome Term dot plot
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


mapping <- labels2colors(1:100)


color_df <- data.frame(group = selected_terms$group,colour = selected_terms$color)



color_df$group <- factor(
  as.character(color_df$group),
  levels = unique(color_df$group)
)


c_vect <- color_df$colour
names(c_vect) <- color_df$group



# make the colorbar as its own heatmap
color_df$var <- 1
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


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'peds_Myeloid_reactome_terms_cell_type_down_pRV_vs_NF.pdf'), width=6.6, height=5)
p / colorbar 
dev.off()


#######################################
#############  FIGURE 7M  #############
#######################################


###RVF vs NF

###Chea
#Up
selected_terms <- subset(combined_output,db=="ChEA_2022")
selected_terms <- subset(selected_terms,direction=="up")
selected_terms <- subset(selected_terms,comparison=="SystolicHF_Donor")
selected_terms <- subset(selected_terms,color %in% mapping[c(1,3,4,8)])


# subset selected terms
selected_terms <- subset(selected_terms, P.value < 0.05)
selected_terms$module_celltype <- paste0(selected_terms$module,'_',selected_terms$celltype)
idx_top_1 <- match(unique(selected_terms$module_celltype),selected_terms$module_celltype)
#idx_top_5 <- sort(c(idx_top_1,idx_top_1+1,idx_top_1+2,idx_top_1+3,idx_top_1+4))

selected_terms<-selected_terms[idx_top_1,]


selected_terms$group <- factor(
  as.character(selected_terms$module_celltype),
  levels = unique(selected_terms$module_celltype)
)


# set max pval
quantile(-log(selected_terms$P.value), 0.95)
max_p <- 10

selected_terms$logp <- -log(selected_terms$P.value)
selected_terms$logp <- ifelse(selected_terms$logp > max_p, max_p, selected_terms$logp)

# remove Reactome Term ID
library(stringr)
selected_terms$Term <- str_replace(selected_terms$Term, "\\ .*", "")

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

# Reactome Term dot plot
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


mapping <- labels2colors(1:100)


color_df <- data.frame(group = selected_terms$group,colour = selected_terms$color)



color_df$group <- factor(
  as.character(color_df$group),
  levels = unique(color_df$group)
)


c_vect <- color_df$colour
names(c_vect) <- color_df$group



# make the colorbar as its own heatmap
color_df$var <- 1
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


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'peds_Myeloid_chea_terms_cell_type_up_RVF_vs_NF.pdf'), width=10, height=5)
p / colorbar 
dev.off()

#Down
selected_terms <- subset(combined_output,db=="ChEA_2022")
selected_terms <- subset(selected_terms,direction=="down")
selected_terms <- subset(selected_terms,comparison=="SystolicHF_Donor")

selected_terms <- subset(selected_terms,color %in% mapping[c(1,3,4,8)])


# subset selected terms
selected_terms <- subset(selected_terms, P.value < 0.05)
selected_terms$module_celltype <- paste0(selected_terms$module,'_',selected_terms$celltype)
idx_top_1 <- match(unique(selected_terms$module_celltype),selected_terms$module_celltype)
#idx_top_5 <- sort(c(idx_top_1,idx_top_1+1,idx_top_1+2,idx_top_1+3,idx_top_1+4))

selected_terms<-selected_terms[idx_top_1,]


selected_terms$group <- factor(
  as.character(selected_terms$module_celltype),
  levels = unique(selected_terms$module_celltype)
)


# set max pval
quantile(-log(selected_terms$P.value), 0.95)
max_p <- 10

selected_terms$logp <- -log(selected_terms$P.value)
selected_terms$logp <- ifelse(selected_terms$logp > max_p, max_p, selected_terms$logp)

# remove Reactome Term ID
library(stringr)
selected_terms$Term <- str_replace(selected_terms$Term, "\\ .*", "")

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

# Reactome Term dot plot
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


mapping <- labels2colors(1:100)


color_df <- data.frame(group = selected_terms$group,colour = selected_terms$color)



color_df$group <- factor(
  as.character(color_df$group),
  levels = unique(color_df$group)
)


c_vect <- color_df$colour
names(c_vect) <- color_df$group



# make the colorbar as its own heatmap
color_df$var <- 1
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


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'peds_Myeloid_chea_terms_cell_type_down_RVF_vs_NF.pdf'), width=6, height=4)
p / colorbar 
dev.off()


###RVF vs pRV
###Chea
#Up
selected_terms <- subset(combined_output,db=="ChEA_2022")
selected_terms <- subset(selected_terms,direction=="up")
selected_terms <- subset(selected_terms,comparison=="SystolicHF_NF")
selected_terms <- subset(selected_terms,color %in% mapping[c(1,3,4,8)])


# subset selected terms
selected_terms <- subset(selected_terms, P.value < 0.05)
selected_terms$module_celltype <- paste0(selected_terms$module,'_',selected_terms$celltype)
idx_top_1 <- match(unique(selected_terms$module_celltype),selected_terms$module_celltype)
#idx_top_5 <- sort(c(idx_top_1,idx_top_1+1,idx_top_1+2,idx_top_1+3,idx_top_1+4))

selected_terms<-selected_terms[idx_top_1,]


selected_terms$group <- factor(
  as.character(selected_terms$module_celltype),
  levels = unique(selected_terms$module_celltype)
)


# set max pval
quantile(-log(selected_terms$P.value), 0.95)
max_p <- 10

selected_terms$logp <- -log(selected_terms$P.value)
selected_terms$logp <- ifelse(selected_terms$logp > max_p, max_p, selected_terms$logp)

# remove Reactome Term ID
library(stringr)
selected_terms$Term <- str_replace(selected_terms$Term, "\\ .*", "")

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

# Reactome Term dot plot
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


mapping <- labels2colors(1:100)


color_df <- data.frame(group = selected_terms$group,colour = selected_terms$color)



color_df$group <- factor(
  as.character(color_df$group),
  levels = unique(color_df$group)
)


c_vect <- color_df$colour
names(c_vect) <- color_df$group



# make the colorbar as its own heatmap
color_df$var <- 1
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


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'peds_Myeloid_chea_terms_cell_type_up_RVF_vs_pRV.pdf'), width=10, height=5)
p / colorbar 
dev.off()

#Down
selected_terms <- subset(combined_output,db=="ChEA_2022")
selected_terms <- subset(selected_terms,direction=="down")
selected_terms <- subset(selected_terms,comparison=="SystolicHF_NF")

selected_terms <- subset(selected_terms,color %in% mapping[c(1,3,4,8)])


# subset selected terms
selected_terms <- subset(selected_terms, P.value < 0.05)
selected_terms$module_celltype <- paste0(selected_terms$module,'_',selected_terms$celltype)
idx_top_1 <- match(unique(selected_terms$module_celltype),selected_terms$module_celltype)
#idx_top_5 <- sort(c(idx_top_1,idx_top_1+1,idx_top_1+2,idx_top_1+3,idx_top_1+4))

selected_terms<-selected_terms[idx_top_1,]


selected_terms$group <- factor(
  as.character(selected_terms$module_celltype),
  levels = unique(selected_terms$module_celltype)
)


# set max pval
quantile(-log(selected_terms$P.value), 0.95)
max_p <- 10

selected_terms$logp <- -log(selected_terms$P.value)
selected_terms$logp <- ifelse(selected_terms$logp > max_p, max_p, selected_terms$logp)

# remove Reactome Term ID
library(stringr)
selected_terms$Term <- str_replace(selected_terms$Term, "\\ .*", "")

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

# Reactome Term dot plot
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


mapping <- labels2colors(1:100)


color_df <- data.frame(group = selected_terms$group,colour = selected_terms$color)



color_df$group <- factor(
  as.character(color_df$group),
  levels = unique(color_df$group)
)


c_vect <- color_df$colour
names(c_vect) <- color_df$group



# make the colorbar as its own heatmap
color_df$var <- 1
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


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'peds_Myeloid_chea_terms_cell_type_down_RVF_vs_pRV.pdf'), width=6, height=4)
p / colorbar 
dev.off()





###pRV vs NF
###Chea
#Up
selected_terms <- subset(combined_output,db=="ChEA_2022")
selected_terms <- subset(selected_terms,direction=="up")
selected_terms <- subset(selected_terms,comparison=="NF_Donor")
selected_terms <- subset(selected_terms,color %in% mapping[c(1,3,4,8)])


# subset selected terms
selected_terms <- subset(selected_terms, P.value < 0.05)
selected_terms$module_celltype <- paste0(selected_terms$module,'_',selected_terms$celltype)
idx_top_1 <- match(unique(selected_terms$module_celltype),selected_terms$module_celltype)
#idx_top_5 <- sort(c(idx_top_1,idx_top_1+1,idx_top_1+2,idx_top_1+3,idx_top_1+4))

selected_terms<-selected_terms[idx_top_1,]


selected_terms$group <- factor(
  as.character(selected_terms$module_celltype),
  levels = unique(selected_terms$module_celltype)
)


# set max pval
quantile(-log(selected_terms$P.value), 0.95)
max_p <- 10

selected_terms$logp <- -log(selected_terms$P.value)
selected_terms$logp <- ifelse(selected_terms$logp > max_p, max_p, selected_terms$logp)

# remove Reactome Term ID
library(stringr)
selected_terms$Term <- str_replace(selected_terms$Term, "\\ .*", "")

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

# Reactome Term dot plot
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


mapping <- labels2colors(1:100)


color_df <- data.frame(group = selected_terms$group,colour = selected_terms$color)



color_df$group <- factor(
  as.character(color_df$group),
  levels = unique(color_df$group)
)


c_vect <- color_df$colour
names(c_vect) <- color_df$group



# make the colorbar as its own heatmap
color_df$var <- 1
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


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'peds_Myeloid_chea_terms_cell_type_up_pRV_vs_NF.pdf'), width=10, height=5)
p / colorbar 
dev.off()

#Down
selected_terms <- subset(combined_output,db=="ChEA_2022")
selected_terms <- subset(selected_terms,direction=="down")
selected_terms <- subset(selected_terms,comparison=="NF_Donor")

selected_terms <- subset(selected_terms,color %in% mapping[c(1,3,4,8)])


# subset selected terms
selected_terms <- subset(selected_terms, P.value < 0.05)
selected_terms$module_celltype <- paste0(selected_terms$module,'_',selected_terms$celltype)
idx_top_1 <- match(unique(selected_terms$module_celltype),selected_terms$module_celltype)
#idx_top_5 <- sort(c(idx_top_1,idx_top_1+1,idx_top_1+2,idx_top_1+3,idx_top_1+4))

selected_terms<-selected_terms[idx_top_1,]


selected_terms$group <- factor(
  as.character(selected_terms$module_celltype),
  levels = unique(selected_terms$module_celltype)
)


# set max pval
quantile(-log(selected_terms$P.value), 0.95)
max_p <- 10

selected_terms$logp <- -log(selected_terms$P.value)
selected_terms$logp <- ifelse(selected_terms$logp > max_p, max_p, selected_terms$logp)

# remove Reactome Term ID
library(stringr)
selected_terms$Term <- str_replace(selected_terms$Term, "\\ .*", "")

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

# Reactome Term dot plot
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


mapping <- labels2colors(1:100)


color_df <- data.frame(group = selected_terms$group,colour = selected_terms$color)



color_df$group <- factor(
  as.character(color_df$group),
  levels = unique(color_df$group)
)


c_vect <- color_df$colour
names(c_vect) <- color_df$group



# make the colorbar as its own heatmap
color_df$var <- 1
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


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'peds_Myeloid_chea_terms_cell_type_down_NF_vs_Donor.pdf'), width=6, height=4)
p / colorbar 
dev.off()


#######################################
#############  FIGURE 7N  #############
#######################################


M2 <- readRDS('~/Downloads/hdWGCNA_TOM/EC_hdWGCNA_by_celltype.rds')

M1 <- readRDS('~/Downloads/hdWGCNA_TOM/Kory_Peds_Hearts/objects/endothelium annotated.rds')

consensus_modules <- GetModules(M2) %>% subset(module != 'grey')
consensus_modules <- consensus_modules[,1:3]
consensus_modules <- subset(consensus_modules, gene_name %in% rownames(M1))
# remove duplicate gene names
consensus_modules <- consensus_modules[match(unique(consensus_modules$gene_name), consensus_modules$gene_name),]
library(dplyr)
score_calc <- consensus_modules %>% group_by(module) %>% group_split()
module_colors <- unique(unlist(lapply(score_calc,'[[','module')))
module_colors <- paste0('M',match(module_colors,mapping))
M1 <- AddModuleScore(M1,lapply(score_calc,'[[','gene_name'),name="module_score")
cols_current <- colnames(M1@meta.data)
cols_current[startsWith(colnames(M1@meta.data),'module_score')] <- paste0('module_',module_colors)
colnames(M1@meta.data) <- cols_current


M1<-SetIdent(M1,value='sub.type')


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'peds_mod_trend_subcluster_EC.pdf'), width=4.5, height=3)

p <- DotPlot(M1,paste0('module_',
  c('M1','M2','M3','M4','M5','M6','M7')),dot.min=0,col.min=0,col.max=2) +
  RotatedAxis() + ylab('')+ xlab('')+
  scale_color_gradient2(high='red', mid='grey95', low='blue') +
  theme(
    panel.border = element_rect(size=1,fill=NA, color='black'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()
) 
p
dev.off()


M1 <- SetIdent(M1, value = "condition")
Idents(M1) <- factor(x = Idents(M1), levels = c('Donor','NF','SystolicHF'))


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'peds_mod_trend_condition_EC.pdf'), width=5, height=2.5)

p <- DotPlot(M1,paste0('module_',
  c('M1','M2','M3','M4','M5','M6','M7')),dot.min=0,col.min=0,col.max=2) +
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
#############  FIGURE 7O  #############
#######################################




M1 <- SetIdent(M1, value = "sub.type")



all_percent_cell <- as.data.frame(table(M1@active.ident)/length(M1@active.ident)*100)

NF_percent_cell <- cbind(as.data.frame(table(subset(M1,condition=="Donor")@active.ident)/length(subset(M1,condition=="Donor")@active.ident)*100),type = "Donor")
NF_percent_cell$sum <- (rev(cumsum(rev(NF_percent_cell$Freq))) - NF_percent_cell$Freq/2)/100
NF_percent_cell$Freq <- NF_percent_cell$Freq/100


pRV_percent_cell <- cbind(as.data.frame(table(subset(M1,condition=="NF")@active.ident)/length(subset(M1,condition=="NF")@active.ident)*100),type = "NF")
pRV_percent_cell$sum <- (rev(cumsum(rev(pRV_percent_cell$Freq))) - pRV_percent_cell$Freq/2)/100
pRV_percent_cell$Freq <- pRV_percent_cell$Freq/100


RVF_percent_cell <- cbind(as.data.frame(table(subset(M1,condition=="SystolicHF")@active.ident)/length(subset(M1,condition=="SystolicHF")@active.ident)*100),type = "SystolicHF")
RVF_percent_cell$sum <- (rev(cumsum(rev(RVF_percent_cell$Freq))) - RVF_percent_cell$Freq/2)/100
RVF_percent_cell$Freq <- RVF_percent_cell$Freq/100


percent_cell_df <- rbind(NF_percent_cell,pRV_percent_cell,RVF_percent_cell)

pdf('~/Downloads/hdWGCNA_TOM/peds_EC_subclust_prev_stacked.pdf',width=5,height=5)
ggplot(percent_cell_df, aes(fill=Var1, y=Freq, x=type,label=round(sum,1))) +  geom_bar(position="stack", stat="identity",width=0.6) + theme_classic() + xlab("Disease State") + ylab("Frequency") + labs(fill="Cell type",color='black') + theme(text = element_text(size=20),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),legend.text=element_text(color="black")) + scale_y_continuous(expand=c(0,0)) + geom_label_repel(aes(type,sum,label=scales::percent(round(Freq,2))),fill=NA,nudge_x=0.5,direction="y")
dev.off()


#######################################
#############  FIGURE 7M  #############
#######################################

M1<-SetIdent(M1,value='condition')

#Get all terms for + regulation of angiogenesis from GO
library(biomaRt)
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl") #uses human ensembl annotations
#gets gene symbol, transcript_id and go_id for all genes annotated with GO:0007507
gene.data <- getBM(attributes=c('hgnc_symbol', 'ensembl_transcript_id', 'go_id'),
                   filters = 'go', values = 'GO:0045765', mart = ensembl)
angio_genes <- unique(gene.data$hgnc_symbol)
angio_genes <- intersect(angio_genes,rownames(M1))
M1 <- AddModuleScore(M1,list(angio_genes),name="AngioGenes")

gene.data <- getBM(attributes=c('hgnc_symbol', 'ensembl_transcript_id', 'go_id'),
                   filters = 'go', values = 'GO:0045766', mart = ensembl)
angio_genes <- unique(gene.data$hgnc_symbol)
angio_genes <- intersect(angio_genes,rownames(M1))
M1 <- AddModuleScore(M1,list(angio_genes),name="PositiveAngioGenes")

gene.data <- getBM(attributes=c('hgnc_symbol', 'ensembl_transcript_id', 'go_id'),
                   filters = 'go', values = 'GO:0016525', mart = ensembl)
angio_genes <- unique(gene.data$hgnc_symbol)
angio_genes <- intersect(angio_genes,rownames(M1))
M1 <- AddModuleScore(M1,list(angio_genes),name="NegativeAngioGenes")

gene.data <- getBM(attributes=c('hgnc_symbol', 'ensembl_transcript_id', 'go_id'),
                   filters = 'go', values = 'GO:0060978', mart = ensembl)
angio_genes <- unique(gene.data$hgnc_symbol)
angio_genes <- intersect(angio_genes,rownames(M1))
M1 <- AddModuleScore(M1,list(angio_genes),name="CoronaryAngioGenes")

gene.data <- getBM(attributes=c('hgnc_symbol', 'ensembl_transcript_id', 'go_id'),
                   filters = 'go', values = 'GO:0001525', mart = ensembl)
angio_genes <- unique(gene.data$hgnc_symbol)
angio_genes <- intersect(angio_genes,rownames(M1))
M1 <- AddModuleScore(M1,list(angio_genes),name="AllAngioGenes")


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'peds_EC_all_angio_genes.pdf'), width=8, height=6)
VlnPlot(M1,'AllAngioGenes1',pt.size=0)
dev.off()

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'peds_EC_coronary_angio_genes.pdf'), width=8, height=6)
VlnPlot(M1,'CoronaryAngioGenes1',pt.size=0)
dev.off()

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'peds_EC_negative_angio_genes.pdf'), width=8, height=6)
VlnPlot(M1,'NegativeAngioGenes1',pt.size=0)
dev.off()

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'peds_EC_positive_angio_genes.pdf'), width=8, height=6)
VlnPlot(M1,'PositiveAngioGenes1',pt.size=0)
dev.off()

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'peds_EC_regulatory_angio_genes.pdf'), width=8, height=6)
VlnPlot(M1,'AngioGenes1',pt.size=0)
dev.off()


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'peds_EC_combined_angio_genes.pdf'), width=10, height=4)
VlnPlot(M1,c('AllAngioGenes1','CoronaryAngioGenes1','PositiveAngioGenes1','NegativeAngioGenes1'),group.by='condition',ncol=4,pt.size=0)
dev.off()

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'peds_EC_MECOM.pdf'), width=3, height=3)
VlnPlot(subset(M1,sub.type=='Artery1'),'MECOM',group.by='condition',pt.size=0)
dev.off()

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'peds_EC_MECOM_feat.pdf'), width=4, height=3)
DotPlot(M1,'MECOM',group.by='sub.type',col.min=0)
dev.off()

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'peds_EC_SMAD1.pdf'), width=3, height=3)
VlnPlot(subset(M1,sub.type=='Capillary'),'SMAD1',group.by='condition',pt.size=0)
dev.off()

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'peds_EC_SMAD1_feat.pdf'), width=4, height=3)
DotPlot(M1,'SMAD1',group.by='sub.type',col.min=0)
dev.off()



#######################################
#############  FIGURE 7N  #############
#######################################


M1$Names_group <- paste0(M1$sub.type,'_',M1$condition)
#Run enrichment by cell type
Idents(M1) <- "Names_group"
combined_set <- data.frame()
combined_output <- data.frame()


bulk_modules<-consensus_modules
bulk_modules$module <- match(bulk_modules$module,mapping)

mods_idx <- c(1,4,7)
cell_types <- unique(M1$sub.type)
comparison <- list(c("SystolicHF","Donor"),c("SystolicHF","NF"),c("NF","Donor"))
for (i in mods_idx){
  for (j in cell_types){
    for (k in comparison){
      key_genes <- subset(bulk_modules,module %in% c(i))$gene_name
      key_genes <- key_genes[key_genes %in% rownames(M1)]

      gene_set <- FindMarkers(M1, ident.1 = paste0(j,"_",k[1]), ident.2 = paste0(j,"_",k[2]),features=key_genes)
      
      gene_set<-subset(gene_set,p_val_adj<0.05)
      if (length(rownames(gene_set))==0){next}
      gene_set$module <- paste0('M',i)
      gene_set$color <- mapping[i]
      gene_set$comparison <- paste0(k[1],'_',k[2])
      gene_set$celltype <- j

      if (length(combined_set) == 0){
        combined_set <- gene_set
      }
      else {
        combined_set <- rbind(combined_set,gene_set)
      }

      gene_enrich <- subset(gene_set,avg_log2FC<0)
      enriched <- enrichR::enrichr(rownames(gene_enrich), dbs)
      Sys.sleep(5)
      for(db in names(enriched)){
          cur_df <- enriched[[db]]
          if (nrow(cur_df) > 1){
            cur_df$db <- db
            cur_df$module <- paste0('M',i)
            cur_df$celltype <- j
            cur_df$comparison <- paste0(k[1],'_',k[2])
            cur_df$color <- mapping[i]
            cur_df$direction <- 'down'
            combined_output <- rbind(combined_output, cur_df)
          }
      }

      gene_enrich <- subset(gene_set,avg_log2FC>0)
      enriched <- enrichR::enrichr(rownames(gene_enrich), dbs)
      Sys.sleep(5)
      for(db in names(enriched)){
          cur_df <- enriched[[db]]
          if (nrow(cur_df) > 1){
            cur_df$db <- db
            cur_df$module <- paste0('M',i)
            cur_df$celltype <- j
            cur_df$comparison <- paste0(k[1],'_',k[2])
            cur_df$color <- mapping[i]
            cur_df$direction <- 'up'
            combined_output <- rbind(combined_output, cur_df)
          }
      }
    }
  }
}


selected_terms <- subset(combined_output,db=="GO_Biological_Process_2023")
selected_terms <- subset(selected_terms,module == 'M7')
selected_terms <- subset(selected_terms,direction == 'up')
selected_terms <- subset(selected_terms,comparison == 'SystolicHF_Donor')
selected_terms <- subset(selected_terms,Adjusted.P.value < 0.05)

selected_terms <- subset(combined_output,db=="GO_Biological_Process_2023")
selected_terms <- subset(selected_terms,module == 'M1')
selected_terms <- subset(selected_terms,direction == 'up')
selected_terms <- subset(selected_terms,comparison == 'SystolicHF_Donor')
selected_terms <- subset(selected_terms,Adjusted.P.value < 0.05)


selected_terms <- subset(combined_output,db=="GO_Biological_Process_2023")
selected_terms <- subset(selected_terms,module == 'M4')
selected_terms <- subset(selected_terms,direction == 'up')
selected_terms <- subset(selected_terms,comparison == 'SystolicHF_Donor')
selected_terms <- subset(selected_terms,Adjusted.P.value < 0.05)

selected_terms <- subset(combined_output,db=="GO_Biological_Process_2023")
selected_terms <- subset(selected_terms,Adjusted.P.value < 0.05)
selected_terms <- subset(selected_terms,module %in% c('M1','M4','M7'))
selected_terms <- subset(selected_terms,direction == 'up')
selected_terms <- subset(selected_terms,comparison %in% c('SystolicHF_NF','NF_Donor'))



# subset selected terms
selected_terms$module_celltype_comparison_direction <- paste0(selected_terms$module,'_',
  selected_terms$celltype,'_',selected_terms$comparison,'_',selected_terms$direction)
idx_top_1 <- match(unique(selected_terms$module_celltype),selected_terms$module_celltype)
#idx_top_5 <- sort(c(idx_top_1,idx_top_1+1,idx_top_1+2,idx_top_1+3,idx_top_1+4))

selected_terms<-selected_terms[idx_top_1,]


selected_terms$group <- factor(
  as.character(selected_terms$module_celltype_comparison_direction),
  levels = unique(selected_terms$module_celltype_comparison_direction)
)


# set max pval
quantile(-log(selected_terms$P.value), 0.95)
max_p <- 10

selected_terms$logp <- -log(selected_terms$P.value)
selected_terms$logp <- ifelse(selected_terms$logp > max_p, max_p, selected_terms$logp)

# remove Reactome Term ID
library(stringr)
selected_terms$Term <- str_replace(selected_terms$Term, "\\(G.*", "")

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

# Reactome Term dot plot
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


mapping <- labels2colors(1:100)


color_df <- data.frame(group = selected_terms$group,colour = selected_terms$color)



color_df$group <- factor(
  as.character(color_df$group),
  levels = unique(color_df$group)
)


c_vect <- color_df$colour
names(c_vect) <- color_df$group



# make the colorbar as its own heatmap
color_df$var <- 1
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


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'pedsEC_GO_terms_cell_type.pdf'), width=8, height=7)
p / colorbar 
dev.off()


selected_terms <- subset(combined_output,db=="ChEA_2022")
selected_terms <- subset(selected_terms,Adjusted.P.value < 0.1)
selected_terms <- subset(selected_terms,module %in% c('M7'))
selected_terms <- subset(selected_terms,celltype %in% c('Vein'))
selected_terms <- subset(selected_terms,direction == c('up','down'))
selected_terms <- subset(selected_terms,comparison %in% c('SystolicHF_Donor','SystolicHF_NF','NF_Donor'))



# subset selected terms
selected_terms$module_celltype_comparison_direction <- paste0(selected_terms$module,'_',
  selected_terms$celltype,'_',selected_terms$comparison,'_',selected_terms$direction)
idx_top_1 <- match(unique(selected_terms$module_celltype),selected_terms$module_celltype)
idx_top_5 <- sort(c(idx_top_1,idx_top_1+1,idx_top_1+2,idx_top_1+3,idx_top_1+4))

selected_terms<-selected_terms[idx_top_5,]


selected_terms$group <- factor(
  as.character(selected_terms$module_celltype_comparison_direction),
  levels = unique(selected_terms$module_celltype_comparison_direction)
)


# set max pval
quantile(-log(selected_terms$P.value), 0.95)
max_p <- 10

selected_terms$logp <- -log(selected_terms$P.value)
selected_terms$logp <- ifelse(selected_terms$logp > max_p, max_p, selected_terms$logp)

# remove Reactome Term ID
library(stringr)
selected_terms$Term <- str_replace(selected_terms$Term, "\\(G.*", "")

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

# Reactome Term dot plot
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


mapping <- labels2colors(1:100)


color_df <- data.frame(group = selected_terms$group,colour = selected_terms$color)



color_df$group <- factor(
  as.character(color_df$group),
  levels = unique(color_df$group)
)


c_vect <- color_df$colour
names(c_vect) <- color_df$group



# make the colorbar as its own heatmap
color_df$var <- 1
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


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'pedsEC_GO_terms_M7.pdf'), width=8, height=7)
p / colorbar 
dev.off()


#######################################
#############  FIGURE 7O  #############
#######################################

DefaultAssay(M1) <- 'SCT'
DefaultAssay(M2) <- 'SCT'


#From https://ashpublications.org/blood/article/122/24/3982/31962/Unraveling-a-novel-transcription-factor-code
#Arterial ECs= Noth induces EFNB2 and blocks EPHB4
M1 <- AddModuleScore(M1,list(c('AFF3','HEY2','SOX17','MSX1','EMX2','TOX2','PRDM16')),name='ArtTF')
M2 <- AddModuleScore(M2,list(c('AFF3','HEY2','SOX17','MSX1','TOX2','PRDM16')),name='ArtTF')


p1 <- VlnPlot(M1,'ArtTF1',group.by='condition',pt.size=0)
p2 <- VlnPlot(M2,'ArtTF1',group.by='group',pt.size=0)

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'peds_RV_art_TF.pdf'), width=3, height=5)
p1/p2
dev.off()


#From https://www.nature.com/articles/s41598-022-05666-1, Fig S2
#Canonical Notch
set <- c('RND1','RASSF10','HEY2','EFNB2',
  'PRICKLE2','SLC46A3','DLL4','RAPGEF5','PLPP3','HEY1','LRRC32','SAT1','HES4','ANKRD33B','FLT1')
M1 <- AddModuleScore(M1,list(set),name='Notch')
M2 <- AddModuleScore(M2,list(set),name='Notch')

p1 <- VlnPlot(M1,'Notch1',group.by='condition',pt.size=0)
p2 <- VlnPlot(M2,'Notch1',group.by='group',pt.size=0)

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'peds_RV_Notch.pdf'), width=3, height=5)
p1/p2
dev.off()


#SMAD1 capullary
p1 <- VlnPlot(M1,'SMAD1',group.by='condition',pt.size=0)
p2 <- VlnPlot(M2,'SMAD1',group.by='group',pt.size=0)

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'peds_RV_Smad.pdf'), width=3, height=5)
p1/p2
dev.off()

#NR2F2 Vein from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7212523/
#Table1
#>=2 bindings sies for NR2F2
#H3K27AC down common
#Fold change > -1
set <- c('CDC6','RRM2','PRIM1','HJURP','RTEL1-TNFRSF6B',
  'RACGAP1','MBOAT1','SELL','GEMIN2','PLVAP','CHAF1A',
  'C1orf112','SIRPB2','SELE','RGCC','ABHD10',
  'CORO2B','SH3BP4','GPRIN3','HAUS7','KPNA2','GCH1',
  'SLC19A2','KDM4D','RGS19','TRDMT1','S1PR4','XRCC3','PARP12',
  'ITGA11','MYCN','C15orf54','OPRL1','CKB','ZNF438','HPCAL1','TMEM223',
  'NUCKS1','PCOLCE','LINC00176','GTPBP3','ELK3','AFF4',
  'NFATC2','SOX18','PPARGC1B','HTR7P1','CYB5R2',
  'ARPP19','FLT4','HIP1','ANKLE1','EXD2','MYRIP','NARS','ZNF22','FILIP1',
  'SPRY4','EMP1','GGT5','ARL15','STX11','NTN4','ACE','ARHGEF37','GAN',
  'TAF1B','TRIM14','PTGES3','SASH1','LIPG','MYO5A','PURB','RPRD1B',
  'STMN3','CUBN','SNHG1','SLC45A3','SIRPB1','MYO5C','CLN6','SRSF3')

M1 <- AddModuleScore(M1,list(set),name='Vein')
M2 <- AddModuleScore(M2,list(set),name='Vein')

p1 <- VlnPlot(subset(M1,sub.type=='Vein'),'Vein1',group.by='condition',pt.size=0)
p2 <- VlnPlot(subset(M2,Names=='Venous'),'Vein1',group.by='group',pt.size=0)

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'peds_RV_vein.pdf'), width=3, height=5)
p1/p2
dev.off()


p1<-VlnPlot(subset(M1,sub.type=='Vein'),'NR2F2',group.by='condition',pt.size=0)
p2<-VlnPlot(subset(M2,Names=='Venous'),'NR2F2',group.by='group',pt.size=0)

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'peds_RV_NR2F2.pdf'), width=3, height=5)
p1/p2
dev.off()




#RBPJ is a notch repressor
VlnPlot(subset(M1,subtype='Vein'),'RBPJ',group.by='condition')
VlnPlot(subset(M2,subtype='Vein'),'RBPJ',group.by='group')


#Arterial specifiers, canonical NOTCH up
DotPlot(M1,c('EFNB2','HEY1','HEY2'),group.by='condition')

#Venous markers
DotPlot(M1,c('NR2F2','NRP2','EPHB4'),group.by='condition')

DotPlot(M1,c('RBPJ'),group.by='condition')


DotPlot(M1,c('DLL1','DLL4','NR2F2','HEY1','EFNB2','EPHB4','MECOM','HEY2'),group.by='condition')







##### OLD


M3 <- subset(M2,CombinedNames=='CM')
DefaultAssay(M3) <- 'RNA'
M3[["RNA"]] <- split(M3[["RNA"]], f = M3$patient)


M3 <- NormalizeData(M3)
M3 <- FindVariableFeatures(M3)
M3 <- ScaleData(M3)
M3 <- RunPCA(M3)


M3 <- IntegrateLayers(
  object = M3, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = FALSE
)

M3 <- FindNeighbors(M3, reduction = "integrated.cca", dims = 1:50)
M3 <- FindClusters(M3, resolution = 2, cluster.name = "cca_clusters")
M3 <- RunUMAP(M3, reduction = "integrated.cca", dims = 1:50, reduction.name = "umap.cca")

p1 <- DimPlot(
  M3,
  reduction = "umap.cca",
  group.by = c("origin"),
  combine = FALSE, label.size = 2
)


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'sn_RV_Peds_CM_Harmonized_UMAP.pdf'), width=5, height=5)
PlotEmbedding(M3,group.by='origin',reduction="umap.cca",point_size=0.2,plot_under=TRUE,plot_theme=umap_theme()+NoLegend(),raster_dpi=400,raster_scale=0.5)
dev.off()


M3 <- subset(M2,CombinedNames=='Myeloid')
DefaultAssay(M3) <- 'RNA'
M3[["RNA"]] <- split(M3[["RNA"]], f = M3$patient)


M3 <- NormalizeData(M3)
M3 <- FindVariableFeatures(M3)
M3 <- ScaleData(M3)
M3 <- RunPCA(M3)


M3 <- IntegrateLayers(
  object = M3, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = FALSE
)

M3 <- FindNeighbors(M3, reduction = "integrated.cca", dims = 1:50)
M3 <- FindClusters(M3, resolution = 2, cluster.name = "cca_clusters")
M3 <- RunUMAP(M3, reduction = "integrated.cca", dims = 1:50, reduction.name = "umap.cca")

p1 <- DimPlot(
  M3,
  reduction = "umap.cca",
  group.by = c("origin"),
  combine = FALSE, label.size = 2
)


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'sn_RV_Peds_Myeloid_Harmonized_UMAP.pdf'), width=5, height=5)
PlotEmbedding(M3,group.by='origin',reduction="umap.cca",point_size=0.2,plot_under=TRUE,plot_theme=umap_theme()+NoLegend(),raster_dpi=400,raster_scale=0.5)
dev.off()



















DefaultAssay(M3) <- 'RNA'
M3 <- NormalizeData(M3)
M3 <- FindVariableFeatures(M3)
M3 <- ScaleData(M3)
M3 <- RunPCA(M3)
M3 <- RunHarmony(M3,'origin')
M3 <- RunUMAP(M3,reduction = "harmony", dims = 1:50, verbose = F) 

RV_and_Peds <- intersect(rownames(subset(M3,origin=='RV')),rownames(subset(M3,origin=='Peds')))



M3 <- SetIdent(M3, value = "origin")

cm_compare <- FindMarkers(M3,ident.1='RV',ident.2='Peds')


write.csv(cm_compare,'~/Downloads/hdWGCNA_TOM/cm_RV_vs_Peds_DEG.csv')



M3 <- subset(M2,CombinedNames=='Myeloid')
DefaultAssay(M3) <- 'RNA'
M3 <- SetIdent(M3, value = "origin")
myeloid_compare <- FindMarkers(M3,ident.1='RV',ident.2='Peds')









RefMerge<-readRDS('~/Downloads/hdWGCNA_TOM/Kory_with_RV_modules_projected.rds')




####CM
RefMerge<-SetActiveWGCNA(RefMerge, "CM_RV2LV")

#RefMerge <- ModuleConnectivity(RefMerge,group.by = 'Names', group_name = 'Cardiomyocytes')


modules <- GetModules(RefMerge)
MEs <- GetMEs(RefMerge, T)
genes_use <- as.character(modules$gene_name)
params <- GetWGCNAParams(RefMerge)

# get the assay
assay <- DefaultAssay(RefMerge)

cells.use <- RefMerge@meta.data %>% subset(get('Names') %in% 'Cardiomyocytes') %>% rownames
MEs <- MEs[cells.use,]

exp_mat <- SeuratObject::LayerData(RefMerge, assay=assay, layer='data')[genes_use,cells.use]

kMEs <- corSparse(X = Matrix::t(exp_mat),Y = as.matrix(MEs))
rownames(kMEs) <- genes_use
kMEs <- as.data.frame(kMEs)

modules <- modules[,1:3]
mods <- levels(modules$module)
colnames(kMEs) <- colnames(MEs)
kMEs <- kMEs[,mods]
colnames(kMEs) <- paste0("kME_", colnames(kMEs))
kMEs <- cbind(modules, kMEs)
RefMerge <- SetModules(RefMerge, kMEs, "CM_RV2LV")



RefMerge <- ModuleExprScore(RefMerge,method='Seurat')


mod_scores <-  GetModuleScores(RefMerge)
mod_scores <- mod_scores[,colnames(mod_scores) != 'grey']

# add hMEs to Seurat meta-data:
RefMerge@meta.data <- cbind(
  RefMerge@meta.data,
  mod_scores
)

# plot with Seurat's DotPlot function
p1 <- DotPlot(
    RefMerge,
    features = mixedsort(colnames(mod_scores)),
    group.by = 'Names'
)

# flip the x/y axes, rotate the axis labels, and change color scheme:
p1 <- p1 +
  RotatedAxis() +
  scale_color_gradient2(high='darkorchid4', mid='grey95', low='seagreen') +
  xlab('') + ylab('')


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'RV2LV_LV_CM_modules_dot.pdf'), width=7.5, height=5)
p1
dev.off()


seurat_ref <- readRDS('~/Downloads/hdWGCNA_TOM/scWGCNA_all_celltypes.rds')
seurat_ref<-SetActiveWGCNA(seurat_ref, "CM")
seurat_ref <- ModuleConnectivity(seurat_ref,group.by = 'Names', group_name = 'CM')

seurat_ref <- ModuleExprScore(seurat_ref,method='Seurat')


mod_scores <-  GetModuleScores(seurat_ref)
mod_scores <- mod_scores[,colnames(mod_scores) != 'grey']

# add hMEs to Seurat meta-data:
seurat_ref@meta.data <- cbind(
  seurat_ref@meta.data,
  mod_scores
)

# plot with Seurat's DotPlot function
p2 <- DotPlot(
    seurat_ref,
    features = mixedsort(colnames(mod_scores)),
    group.by = 'Names'
)

# flip the x/y axes, rotate the axis labels, and change color scheme:
p2 <- p2 +
  RotatedAxis() +
  scale_color_gradient2(high='darkorchid4', mid='grey95', low='seagreen') +
  xlab('') + ylab('')


pdf(paste0('~/Downloads/hdWGCNA_TOM/', "LV_CM_modules_dot.pdf"), width=7.5, height=5)
p2
dev.off()

pdf(paste0('~/Downloads/hdWGCNA_TOM/', "LV_CM_LV2RM_modules_dot.pdf"), width=7.5, height=8)
wrap_plots(list(p1,p2),ncol=1)
dev.off()

###EC

RefMerge<-SetActiveWGCNA(RefMerge, "EC_RV2LV")

#RefMerge <- ModuleConnectivity(RefMerge,group.by = 'Names', group_name = 'Cardiomyocytes')


modules <- GetModules(RefMerge)
MEs <- GetMEs(RefMerge, T)
genes_use <- as.character(modules$gene_name)
params <- GetWGCNAParams(RefMerge)

# get the assay
assay <- DefaultAssay(RefMerge)

cells.use <- RefMerge@meta.data %>% subset(get('Names') %in% 'Endothelium') %>% rownames
MEs <- MEs[cells.use,]

exp_mat <- SeuratObject::LayerData(RefMerge, assay=assay, layer='data')[genes_use,cells.use]

kMEs <- corSparse(X = Matrix::t(exp_mat),Y = as.matrix(MEs))
rownames(kMEs) <- genes_use
kMEs <- as.data.frame(kMEs)

modules <- modules[,1:3]
mods <- levels(modules$module)
colnames(kMEs) <- colnames(MEs)
kMEs <- kMEs[,mods]
colnames(kMEs) <- paste0("kME_", colnames(kMEs))
kMEs <- cbind(modules, kMEs)
RefMerge <- SetModules(RefMerge, kMEs, "EC_RV2LV")



RefMerge <- ModuleExprScore(RefMerge,method='Seurat')


mod_scores <-  GetModuleScores(RefMerge)
mod_scores <- mod_scores[,colnames(mod_scores) != 'grey']

# add hMEs to Seurat meta-data:
RefMerge@meta.data <- cbind(
  RefMerge@meta.data,
  mod_scores
)

# plot with Seurat's DotPlot function
p1 <- DotPlot(
    RefMerge,
    features = mixedsort(colnames(mod_scores)),
    group.by = 'Names'
)

# flip the x/y axes, rotate the axis labels, and change color scheme:
p1 <- p1 +
  RotatedAxis() +
  scale_color_gradient2(high='darkorchid4', mid='grey95', low='seagreen') +
  xlab('') + ylab('')


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'RV2LV_LV_EC_modules_dot.pdf'), width=7.5, height=5)
p1
dev.off()


seurat_ref<-SetActiveWGCNA(seurat_ref, "EC")
seurat_ref <- ModuleConnectivity(seurat_ref,group.by = 'Names', group_name = 'EC')

seurat_ref <- ModuleExprScore(seurat_ref,method='Seurat')


mod_scores <-  GetModuleScores(seurat_ref)
mod_scores <- mod_scores[,colnames(mod_scores) != 'grey']

# add hMEs to Seurat meta-data:
seurat_ref@meta.data <- cbind(
  seurat_ref@meta.data,
  mod_scores
)

# plot with Seurat's DotPlot function
p2 <- DotPlot(
    seurat_ref,
    features = mixedsort(colnames(mod_scores)),
    group.by = 'Names'
)

# flip the x/y axes, rotate the axis labels, and change color scheme:
p2 <- p2 +
  RotatedAxis() +
  scale_color_gradient2(high='darkorchid4', mid='grey95', low='seagreen') +
  xlab('') + ylab('')


pdf(paste0('~/Downloads/hdWGCNA_TOM/', "LV_EC_modules_dot.pdf"), width=7.5, height=5)
p2
dev.off()

pdf(paste0('~/Downloads/hdWGCNA_TOM/', "LV_EC_LV2RM_modules_dot.pdf"), width=7.5, height=8)
wrap_plots(list(p1,p2),ncol=1)
dev.off()



###FB

RefMerge<-SetActiveWGCNA(RefMerge, "FB_RV2LV")

#RefMerge <- ModuleConnectivity(RefMerge,group.by = 'Names', group_name = 'Cardiomyocytes')


modules <- GetModules(RefMerge)
MEs <- GetMEs(RefMerge, T)
genes_use <- as.character(modules$gene_name)
params <- GetWGCNAParams(RefMerge)

# get the assay
assay <- DefaultAssay(RefMerge)

cells.use <- RefMerge@meta.data %>% subset(get('Names') %in% 'Fibroblasts') %>% rownames
MEs <- MEs[cells.use,]

exp_mat <- SeuratObject::LayerData(RefMerge, assay=assay, layer='data')[genes_use,cells.use]

kMEs <- corSparse(X = Matrix::t(exp_mat),Y = as.matrix(MEs))
rownames(kMEs) <- genes_use
kMEs <- as.data.frame(kMEs)

modules <- modules[,1:3]
mods <- levels(modules$module)
colnames(kMEs) <- colnames(MEs)
kMEs <- kMEs[,mods]
colnames(kMEs) <- paste0("kME_", colnames(kMEs))
kMEs <- cbind(modules, kMEs)
RefMerge <- SetModules(RefMerge, kMEs, "FB_RV2LV")



RefMerge <- ModuleExprScore(RefMerge,method='Seurat')


mod_scores <-  GetModuleScores(RefMerge)
mod_scores <- mod_scores[,colnames(mod_scores) != 'grey']

# add hMEs to Seurat meta-data:
RefMerge@meta.data <- cbind(
  RefMerge@meta.data,
  mod_scores
)

# plot with Seurat's DotPlot function
p1 <- DotPlot(
    RefMerge,
    features = mixedsort(colnames(mod_scores)),
    group.by = 'Names'
)

# flip the x/y axes, rotate the axis labels, and change color scheme:
p1 <- p1 +
  RotatedAxis() +
  scale_color_gradient2(high='darkorchid4', mid='grey95', low='seagreen') +
  xlab('') + ylab('')


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'RV2LV_LV_FB_modules_dot.pdf'), width=7.5, height=5)
p1
dev.off()


seurat_ref<-SetActiveWGCNA(seurat_ref, "FB")
seurat_ref <- ModuleConnectivity(seurat_ref,group.by = 'Names', group_name = 'FB')

seurat_ref <- ModuleExprScore(seurat_ref,method='Seurat')


mod_scores <-  GetModuleScores(seurat_ref)
mod_scores <- mod_scores[,colnames(mod_scores) != 'grey']

# add hMEs to Seurat meta-data:
seurat_ref@meta.data <- cbind(
  seurat_ref@meta.data,
  mod_scores
)

# plot with Seurat's DotPlot function
p2 <- DotPlot(
    seurat_ref,
    features = mixedsort(colnames(mod_scores)),
    group.by = 'Names'
)

# flip the x/y axes, rotate the axis labels, and change color scheme:
p2 <- p2 +
  RotatedAxis() +
  scale_color_gradient2(high='darkorchid4', mid='grey95', low='seagreen') +
  xlab('') + ylab('')


pdf(paste0('~/Downloads/hdWGCNA_TOM/', "LV_FB_modules_dot.pdf"), width=7.5, height=5)
p2
dev.off()

pdf(paste0('~/Downloads/hdWGCNA_TOM/', "LV_FB_LV2RM_modules_dot.pdf"), width=7.5, height=8)
wrap_plots(list(p1,p2),ncol=1)
dev.off()



###Myeloid

RefMerge<-SetActiveWGCNA(RefMerge, "Myeloid_RV2LV")

#RefMerge <- ModuleConnectivity(RefMerge,group.by = 'Names', group_name = 'Cardiomyocytes')


modules <- GetModules(RefMerge)
MEs <- GetMEs(RefMerge, T)
genes_use <- as.character(modules$gene_name)
params <- GetWGCNAParams(RefMerge)

# get the assay
assay <- DefaultAssay(RefMerge)

cells.use <- RefMerge@meta.data %>% subset(get('Names') %in% 'Myeloid') %>% rownames
MEs <- MEs[cells.use,]

exp_mat <- SeuratObject::LayerData(RefMerge, assay=assay, layer='data')[genes_use,cells.use]

kMEs <- corSparse(X = Matrix::t(exp_mat),Y = as.matrix(MEs))
rownames(kMEs) <- genes_use
kMEs <- as.data.frame(kMEs)

modules <- modules[,1:3]
mods <- levels(modules$module)
colnames(kMEs) <- colnames(MEs)
kMEs <- kMEs[,mods]
colnames(kMEs) <- paste0("kME_", colnames(kMEs))
kMEs <- cbind(modules, kMEs)
RefMerge <- SetModules(RefMerge, kMEs, "Myeloid_RV2LV")



RefMerge <- ModuleExprScore(RefMerge,method='Seurat')


mod_scores <-  GetModuleScores(RefMerge)
mod_scores <- mod_scores[,colnames(mod_scores) != 'grey']

# add hMEs to Seurat meta-data:
RefMerge@meta.data <- cbind(
  RefMerge@meta.data,
  mod_scores
)

# plot with Seurat's DotPlot function
p1 <- DotPlot(
    RefMerge,
    features = mixedsort(colnames(mod_scores)),
    group.by = 'Names'
)

# flip the x/y axes, rotate the axis labels, and change color scheme:
p1 <- p1 +
  RotatedAxis() +
  scale_color_gradient2(high='darkorchid4', mid='grey95', low='seagreen') +
  xlab('') + ylab('')


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'RV2LV_LV_Myeloid_modules_dot.pdf'), width=7.5, height=5)
p1
dev.off()


seurat_ref<-SetActiveWGCNA(seurat_ref, "Myeloid")
seurat_ref <- ModuleConnectivity(seurat_ref,group.by = 'Names', group_name = 'Myeloid')

seurat_ref <- ModuleExprScore(seurat_ref,method='Seurat')


mod_scores <-  GetModuleScores(seurat_ref)
mod_scores <- mod_scores[,colnames(mod_scores) != 'grey']

# add hMEs to Seurat meta-data:
seurat_ref@meta.data <- cbind(
  seurat_ref@meta.data,
  mod_scores
)

# plot with Seurat's DotPlot function
p2 <- DotPlot(
    seurat_ref,
    features = mixedsort(colnames(mod_scores)),
    group.by = 'Names'
)

# flip the x/y axes, rotate the axis labels, and change color scheme:
p2 <- p2 +
  RotatedAxis() +
  scale_color_gradient2(high='darkorchid4', mid='grey95', low='seagreen') +
  xlab('') + ylab('')


pdf(paste0('~/Downloads/hdWGCNA_TOM/', "LV_Myeloid_modules_dot.pdf"), width=7.5, height=5)
p2
dev.off()

pdf(paste0('~/Downloads/hdWGCNA_TOM/', "LV_Myeloid_LV2RM_modules_dot.pdf"), width=7.5, height=8)
wrap_plots(list(p1,p2),ncol=1)
dev.off()



#######################################
#############  FIGURE 7B  #############
#######################################


RefMerge<-readRDS('~/Downloads/hdWGCNA_TOM/Kory_with_RV_modules_projected.rds')

RefMerge<-SetActiveWGCNA(RefMerge, "CM_RV2LV")


plot_list <- PlotModulePreservation(
  RefMerge,
  name="CM_pres",
  statistics = "summary"
)

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'RV2LV_CM_pres.pdf'), width=5, height=5)
wrap_plots(plot_list, ncol=2)
dev.off()

#######################################
#############  FIGURE 7C  #############
#######################################
#Library EnrichR

#CM not preserved - M1, M5
#CM preserved - M7, M9, M10, M11, M14, M14

#EC not preserved - M1, M2, M4, M5, M10, M11, M12
#EC preserved - M8

seurat_ref <- readRDS('~/Downloads/hdWGCNA_TOM/scWGCNA_all_celltypes.rds')
seurat_ref<-SetActiveWGCNA(seurat_ref, "CM")


dbs <- c('GO_Biological_Process_2023','GO_Cellular_Component_2023','GO_Molecular_Function_2023','ChEA_2022','Reactome_2022')

# perform enrichment tests
seurat_ref <- RunEnrichr(
  seurat_ref,
  dbs=dbs, # character vector of enrichr databases to test
  max_genes = 2000 # number of genes per module to test. use max_genes = Inf to choose all genes!
)

# retrieve the output table
enrich_df <- GetEnrichrTable(seurat_ref)

EnrichrBarPlot(
  seurat_ref,
  outdir = "~/Downloads/hdWGCNA_TOM/sc_enrichr_plot_CM", 
  n_terms = 5,
  plot_size = c(5,4), # width, height of the output .pdfs
  logscale=TRUE # do you want to show the enrichment as a log scale?
)


seurat_ref<-SetActiveWGCNA(seurat_ref, "EC")


dbs <- c('GO_Biological_Process_2023','GO_Cellular_Component_2023','GO_Molecular_Function_2023','ChEA_2022','Reactome_2022')

# perform enrichment tests
seurat_ref <- RunEnrichr(
  seurat_ref,
  dbs=dbs, # character vector of enrichr databases to test
  max_genes = 2000 # number of genes per module to test. use max_genes = Inf to choose all genes!
)

# retrieve the output table
enrich_df <- GetEnrichrTable(seurat_ref)

EnrichrBarPlot(
  seurat_ref,
  outdir = "~/Downloads/hdWGCNA_TOM/sc_enrichr_plot_EC", 
  n_terms = 5,
  plot_size = c(5,4), # width, height of the output .pdfs
  logscale=TRUE # do you want to show the enrichment as a log scale?
)


#######################################
#############  FIGURE 7D  #############
#######################################
library(GeneOverlap)

seurat_ref <- readRDS('~/Downloads/hdWGCNA_TOM/scWGCNA_all_celltypes.rds')
seurat_ref<-SetActiveWGCNA(seurat_ref, "CM")


bulk_rv_vs_lv <- read.csv(paste0('~/Downloads/hdWGCNA_TOM/', 'RV_LV_align_NRVM_NRMV_ARVM.csv'))
genes <- toupper(bulk_rv_vs_lv[,2])
#nrvm <- bulk_rv_vs_lv[2:348,14]
arvm <- bulk_rv_vs_lv[,20]
arvm.p <-  bulk_rv_vs_lv[,22]
arvm_rv <-genes[arvm.p < 0.05 & arvm < 0 & genes %in% rownames(seurat_ref)]
arvm_lv <- genes[arvm.p < 0.05 & arvm > 0 & genes %in% rownames(seurat_ref)]
arvm_share <-genes[arvm.p >= 1 & arvm < 0.5 & arvm > -0.5 & genes %in% rownames(seurat_ref)]


# load modules
modules <- GetModules(seurat_ref)
mods <- levels(modules$module)
genome.size <- nrow(modules)


overlap_df_rv <- do.call(rbind, lapply(mods, function(cur_mod){

  cur_genes <- modules %>% subset(module == cur_mod) %>% .$gene_name

  cur_overlap <- testGeneOverlap(newGeneOverlap(
      cur_genes,
      arvm_rv,
      genome.size=genome.size
  ))

  cur_overlap <- data.frame(
    'odds.ratio' = cur_overlap@odds.ratio,
    'pval' = cur_overlap@pval,
    'Jaccard' = cur_overlap@Jaccard,
    'size_intersection' = length(cur_overlap@intersection),
    'module' = cur_mod
  )

  cur_overlap

})) %>% as.data.frame()

overlap_df_lv <- do.call(rbind, lapply(mods, function(cur_mod){

  cur_genes <- modules %>% subset(module == cur_mod) %>% .$gene_name

  cur_overlap <- testGeneOverlap(newGeneOverlap(
      cur_genes,
      arvm_lv,
      genome.size=genome.size
  ))

  cur_overlap <- data.frame(
    'odds.ratio' = cur_overlap@odds.ratio,
    'pval' = cur_overlap@pval,
    'Jaccard' = cur_overlap@Jaccard,
    'size_intersection' = length(cur_overlap@intersection),
    'module' = cur_mod
  )

  cur_overlap

})) %>% as.data.frame()

overlap_df_share <- do.call(rbind, lapply(mods, function(cur_mod){

  cur_genes <- modules %>% subset(module == cur_mod) %>% .$gene_name

  cur_overlap <- testGeneOverlap(newGeneOverlap(
      cur_genes,
      arvm_share,
      genome.size=genome.size
  ))

  cur_overlap <- data.frame(
    'odds.ratio' = cur_overlap@odds.ratio,
    'pval' = cur_overlap@pval,
    'Jaccard' = cur_overlap@Jaccard,
    'size_intersection' = length(cur_overlap@intersection),
    'module' = cur_mod
  )

  cur_overlap

})) %>% as.data.frame()

overlap_df_rv <- overlap_df_rv %>% mutate(fdr=p.adjust(pval, method='fdr'))
overlap_df_rv <- overlap_df_rv %>% subset(module != 'grey')



overlap_df_lv <- overlap_df_lv %>% mutate(fdr=p.adjust(pval, method='fdr'))
overlap_df_lv <- overlap_df_lv %>% subset(module != 'grey')

overlap_df_share <- overlap_df_share %>% mutate(fdr=p.adjust(pval, method='fdr'))
overlap_df_share <- overlap_df_share %>% subset(module != 'grey')

# Plot as a lollipop

overlap_df_rv$shape <- ifelse(overlap_df_rv$fdr < 0.05, 21, 4)
overlap_df_rv <- overlap_df_rv %>% arrange(odds.ratio, descending=TRUE)
overlap_df_rv$module <- factor(as.character(overlap_df_rv$module), levels=as.character(overlap_df_rv$module))

mod_colors <- dplyr::select(modules, c(module, color)) %>%
  distinct
cp <- mod_colors$color; names(cp) <- mod_colors$module

p <- overlap_df_rv %>%
  ggplot(aes(y=module, x=odds.ratio, size= size_intersection, color=module)) +
  geom_segment(aes(y=module, yend=module, x=0, xend=odds.ratio), size=0.5, color='grey') +
  geom_point() +
  geom_point(shape=overlap_df_rv$shape, color='black', fill=NA) +
  scale_color_manual(values=cp, guide='none') +
  ylab('') + xlab("Odds ratio") +
  scale_x_continuous(breaks = c(0, 1, 2,3)) +
  labs(size='Size\nintersection') +
  ggtitle('Overlap with SFARI genes') +

  theme(
    panel.border = element_rect(size=1, color='black', fill=NA),
    axis.line.y = element_blank(),
    axis.line.x = element_blank(),
    plot.title = element_text(hjust=0.5, face='plain')
  )

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'RVbulk_R2LV_CM_overlap.pdf'), width=4.5, height=5)
p
dev.off()


#######################################
#############  FIGURE 7D  #############
#######################################
rm(RefMerge)
rm(seurat_ref)
seurat_ref <- readRDS('/Volumes/Extreme SSD/Final_Analysis/CellTypes/Post_R3_FINAL_with_counts.rds')
load('~/Downloads/hdWGCNA_TOM/GSE183852_DCM_Integrated.Robj')


sc_modules <- read.csv(paste0('~/Downloads/hdWGCNA_TOM/', 'sc_heart_modules.csv'))
seurat_ref <- AddModuleScore(seurat_ref,list(
  subset(sc_modules,module=="CM-M1")$gene_name,
  subset(sc_modules,module=="CM-M2")$gene_name,
  subset(sc_modules,module=="CM-M3")$gene_name,
  subset(sc_modules,module=="CM-M4")$gene_name,
  subset(sc_modules,module=="CM-M5")$gene_name,
  subset(sc_modules,module=="CM-M6")$gene_name,
  subset(sc_modules,module=="CM-M7")$gene_name
  ),name='CM_score')
RefMerge <- AddModuleScore(RefMerge,list(
  subset(sc_modules,module=="CM-M1")$gene_name,
  subset(sc_modules,module=="CM-M2")$gene_name,
  subset(sc_modules,module=="CM-M3")$gene_name,
  subset(sc_modules,module=="CM-M4")$gene_name,
  subset(sc_modules,module=="CM-M5")$gene_name,
  subset(sc_modules,module=="CM-M6")$gene_name,
  subset(sc_modules,module=="CM-M7")$gene_name
  )
  ,name='CM_score')


FeaturePlot(RefMerge,'CM_score1')
DotPlot(RefMerge,c('CM_score1','CM_score2','CM_score3','CM_score4','CM_score5','CM_score6','CM_score7'),group.by='Names',split.by='condition')
DotPlot(seurat_ref,c('CM_score1','CM_score2','CM_score3','CM_score4','CM_score5','CM_score6','CM_score7'),group.by='Names',split.by='group')


cellsofint<- colnames(RefMerge)[RefMerge@meta.data$'CM-M6'>5 & RefMerge@meta.data$'EC-M7'>5 & RefMerge@meta.data$'CM-M11'>5]
DimPlot(RefMerge, cells.highlight=cellsofint)
FindMarkers(RefMerge,ident.1=cellsofint)



#cellsofint<- colnames(seurat_ref)[seurat_ref@meta.data$'CM-M6'>0 & seurat_ref@meta.data$'EC-M7'>0 & seurat_ref@meta.data$'CM-M11'>0]
cellsofint<- colnames(seurat_ref)[seurat_ref@meta.data$'CM-M6'>10]

cellstocomp<- colnames(seurat_ref)[seurat_ref@meta.data$'CM-M6'<10 & seurat_ref@meta.data$Names == "CM"]

DimPlot(seurat_ref, cells.highlight=cellsofint,raster=T)
FindMarkers(seurat_ref,ident.1=cellsofint,ident.2=cellstocomp)




