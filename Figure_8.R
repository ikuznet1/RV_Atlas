library(Seurat)
library(hdWGCNA)
library(ggeasy)
library(dplyr)


source('./dependencies/shared/spatial_functions.R')



#######################################
#############  FIGURE 8A  #############
#######################################
M1 <- readRDS('./dependencies/shared/all_peds_data.rds')
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


pdf(paste0('./output/', 'snPeds_Vln.pdf'), width=4, height=20)
p
dev.off()



M2 <- readRDS('./dependencies/shared/Post_R3_FINAL_with_counts.rds')


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
#saveRDS(M2,'./output/RV_Peds_merge.rds')
#M2 <- readRDS('./output/RV_Peds_merge.rds')
#M1 <- subset(M2,origin==TRUE)
#saveRDS(M1,'./output/Peds_clean.rds')
#M2 <- readRDS('./dependencies/Figure_8/RV_Peds_merge.rds')


pdf(paste0('./output/', 'sn_RV_Peds_UMAP.pdf'), width=5, height=5)
PlotEmbedding(M2,group.by='CombinedNames',point_size=0.2,plot_under=TRUE,plot_theme=umap_theme()+NoLegend(),raster_dpi=400,raster_scale=0.5)
dev.off()

pdf(paste0('./output/', 'sn_RV_Peds_UMAP_origin.pdf'), width=5, height=5)
PlotEmbedding(M2,group.by='origin',point_size=0.2,plot_under=TRUE,plot_theme=umap_theme()+NoLegend(),raster_dpi=400,raster_scale=0.5)
dev.off()


pdf(paste0('./output/', 'sn_Peds_UMAP_reprojected.pdf'), width=5, height=5)
PlotEmbedding(M1,group.by='NewNames',point_size=0.2,plot_under=TRUE,plot_theme=umap_theme()+NoLegend(),raster_dpi=400,raster_scale=0.5)
dev.off()

M2$condition[M2$condition == "NF"] = "pRV" 
M2$condition[M2$condition == "Donor"] = "NF" 
M2$condition[M2$condition == "SystolicHF"] = "RVF" 



M2$group[is.na(M2$group)] <- M2$condition[is.na(M2$group)]


#######################################
#############  FIGURE 8B  #############
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
pdf('./output/RV_Peds_prev_stacked.pdf',width=6,height=2.5)
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
pdf('./output/Peds_RV_clust_freq.pdf',width=8,height=5)
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
pdf('./output/Peds_clust_freq.pdf',width=8,height=5)
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
#############  FIGURE 8C  #############
#######################################

########Embed bulk in single nuc peds


######Load module
consensus_modules <- read.csv("./dependencies/shared/bulk_heart_modules.csv")
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


#consensus_modules <- read.csv("./dependencies/shared/bulk_heart_modules.R")
#consensus_modules <- consensus_modules[,1:3]
#consensus_modules <- subset(consensus_modules, gene_name %in% rownames(M2))
# remove duplicate gene names
#consensus_modules <- consensus_modules[match(unique(consensus_modules$gene_name), consensus_modules$gene_name),]


library(dplyr)
score_calc <- consensus_modules %>% group_by(module) %>% group_split()
module_colors <- unique(unlist(lapply(score_calc,'[[','module')))
module_colors <- paste0('M',match(module_colors,mapping))
#saveRDS(M2, './output/scWGCNA_RV_Peds_bulk2sn_projection.rds')
#M2<- readRDS('./output/scWGCNA_RV_Peds_bulk2sn_projection.rds')

DefaultAssay(M2) <- 'SCT'



#rm(seurat_ref)
#gc()
#seurat_ref<-readRDS('./dependencies/shared/Post_R3_FINAL_with_counts.rds')
#seurat_ref <- SetIdent(seurat_ref, value = "Names")
#seurat_ref@meta.data <- cbind(seurat_ref@meta.data, MEs)



M2 <- AddModuleScore(M2,lapply(score_calc,'[[','gene_name'),name="module_score")


cols_current <- colnames(M2@meta.data)
cols_current[startsWith(colnames(M2@meta.data),'module_score')] <- paste0('module_',module_colors)
colnames(M2@meta.data) <- cols_current

M2$origin[M2$origin] = 'Peds'
M2$origin[M2$origin == FALSE] = 'RV'

M2$CombinedNamesSplit <- paste0(M2$CombinedNames,'_',M2$origin)

pdf(paste0('./output/', 'RV_Peds_Dot.pdf'), width=7, height=5)
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


pdf(paste0('./output/', 'RV_Peds_Dot_ordered.pdf'), width=7, height=5)

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


pdf(paste0('./output/', 'RV_Peds_Dot_disease_ordered.pdf'), width=7, height=4)

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



M2$group_split <- paste0(M2$group,'_',M2$origin)
M2$group_split <- factor(M2$group_split,levels=c('NF_RV','pRV_RV','RVF_RV','NF_Peds','pRV_Peds','RVF_Peds'))

pdf(paste0('./output/', 'RV_Peds_Dot_CM_Peds.pdf'), width=7, height=3)
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

pdf(paste0('./output/', 'RV_Peds_Dot_CM_RV.pdf'), width=7, height=3)
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

pdf(paste0('./output/', 'RV_Peds_Dot_Myeloid_Peds.pdf'), width=7, height=3)
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

pdf(paste0('./output/', 'RV_Peds_Dot_Myeloid_RV.pdf'), width=7, height=3)
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
#############  FIGURE 8D  #############
#######################################

M2 <- SetIdent(M2, value = "groupSplit")
slot(M2$SCT@SCTModel.list[[1]], 'median_umi') = median(M2$SCT@SCTModel.list[[1]]@cell.attributes$umi)

gene_set_RV <- FindMarkers(M2, ident.1 = "RVF_RV", ident.2 = "NF_RV",recorrect_umi=F)
gene_set_Peds <- FindMarkers(M2, ident.1 = "RVF_Peds", ident.2 = "NF_Peds",recorrect_umi=F)




shared <- intersect(rownames(gene_set_RV),rownames(gene_set_Peds))
dataset <- data.frame(Peds=gene_set_Peds[shared,]$avg_log2FC,RV=gene_set_RV[shared,]$avg_log2FC)
rownames(dataset) <- shared
#Cor 0.5288722

pdf(paste0('./output/', 'Peds_vs_RV.pdf'), width=8, height=8)
ggplot(dataset, aes(x = RV, y=Peds)) + geom_point() + 
  geom_text_repel(label=rownames(dataset),max.overlaps = 50) + theme_classic()
dev.off()

shared <- intersect(rownames(subset(gene_set_RV,pct.1>0.05 & pct.2>0.05)),
  rownames(subset(gene_set_Peds,pct.1>0.05 & pct.2>0.05)))

dataset <- data.frame(Peds=gene_set_Peds[shared,]$avg_log2FC,RV=gene_set_RV[shared,]$avg_log2FC)
rownames(dataset) <- shared


pdf(paste0('./output/', 'Peds_vs_RV_5percent.pdf'), width=8, height=8)
ggplot(dataset, aes(x = RV, y=Peds)) + geom_point() + 
  geom_text_repel(label=rownames(dataset),max.overlaps = 20) + theme_classic()
dev.off()
#Cor 0.4706919

#######################################
#############  FIGURE 8E  #############
#######################################

M1 <- readRDS('./dependencies/Figure_8/cardiomyocyte annotated.rds')
M1$Names <- M1$cell.type
M1$NewNames <- M1$Names
M1$Subnames <- M1$sub.type
M1$NewSubnames <- M1$Subnames
M1 <- SetIdent(M1, value = "NewSubnames")
DefaultAssay(M1) <- 'SCT'


consensus_modules <- read.csv("./dependencies/shared/bulk_heart_modules.csv")
consensus_modules <- consensus_modules[,1:3]
consensus_modules <- subset(consensus_modules, gene_name %in% rownames(M1))
# remove duplicate gene names
consensus_modules <- consensus_modules[match(unique(consensus_modules$gene_name), consensus_modules$gene_name),]
library(dplyr)
score_calc <- consensus_modules %>% group_by(module) %>% group_split()
module_colors <- unique(unlist(lapply(score_calc,'[[','module')))
mapping <- labels2colors(1:100)
module_colors <- paste0('M',match(module_colors,mapping))
M1 <- AddModuleScore(M1,lapply(score_calc,'[[','gene_name'),name="module_score")
cols_current <- colnames(M1@meta.data)
cols_current[startsWith(colnames(M1@meta.data),'module_score')] <- paste0('module_',module_colors)
colnames(M1@meta.data) <- cols_current

Idents(M1) <- factor(x = Idents(M1), levels = sort(levels(M1)))


#Dot Plot of enrichment cell type of CM enriched modules

pdf(paste0('./output/', 'peds_mod_trend_subcluster_CM.pdf'), width=4.5, height=3)

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


pdf(paste0('./output/', 'peds_mod_trend_condition_CM.pdf'), width=5, height=2.5)

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
#############  FIGURE 8F  #############
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


pdf(paste0('./output/', 'peds_CM_by_cluster_terms_cell_type_up_SystolicHF_vs_Donor.pdf'), width=6, height=4)
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


pdf(paste0('./output/', 'peds_CM_by_cluster_terms_cell_type_down_SystolicHF_vs_Donor.pdf'), width=6, height=4)
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


pdf(paste0('./output/', 'peds_CM_by_cluster_terms_cell_type_both_SystolicHF_vs_Donor.pdf'), width=6, height=9)
p / colorbar 
dev.off()



#Deep dive M2

bulk_modules <- read.csv("./dependencies/shared/bulk_heart_modules.csv")
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
pdf('./output/CM_Peds_M2_enrichr_up.pdf',width=5,height=2.5)
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
pdf('./output/CM_Peds_M2_enrichr_down.pdf',width=5,height=2.5)
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

pdf('./output/CM_Peds_M2_enrichr_up_down.pdf',width=6,height=4)
p1/p2
dev.off()



enriched <- enrichr(c(M10_genes_up,M25_genes_up,M26_genes_up,M28_genes_up), dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('./output/CM_Peds_mito_enrichr_up.pdf',width=6,height=3)
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


bulk_modules <- read.csv("./dependencies/shared/bulk_heart_modules.R")
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
pdf('./output/CM_Peds_M2_enrichr_up_RVF_vs_pRV.pdf.pdf',width=5,height=2.5)
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
pdf('./output/CM_Peds_M2_enrichr_down_RVF_vs_pRV.pdf.pdf',width=5,height=2.5)
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

pdf('./output/CM_Peds_M2_enrichr_up_down_RVF_vs_pRV.pdf',width=6,height=4)
p1/p2
dev.off()



enriched <- enrichr(c(M10_genes_up,M25_genes_up,M26_genes_up,M28_genes_up), dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('./output/CM_Peds_mito_enrichr_up_RVF_vs_pRV.pdf.pdf',width=6,height=3)
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


bulk_modules <- read.csv("./dependencies/shared/bulk_heart_modules.R")
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
pdf('./output/CM_Peds_M2_enrichr_up_pRV_vs_NF.pdf.pdf',width=5,height=2.5)
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
pdf('./output/CM_Peds_M2_enrichr_down_pRV_vs_NF.pdf.pdf',width=5,height=2.5)
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

pdf('./output/CM_Peds_M2_enrichr_up_down_pRV_vs_NF.pdf',width=6,height=4)
p1/p2
dev.off()



enriched <- enrichr(c(M10_genes_up,M25_genes_up,M26_genes_up,M28_genes_up), dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('./output/CM_Peds_mito_enrichr_up_pRV_vs_NF.pdf.pdf',width=6,height=3)
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

















# ##### OLD


# M3 <- subset(M2,CombinedNames=='CM')
# DefaultAssay(M3) <- 'RNA'
# M3[["RNA"]] <- split(M3[["RNA"]], f = M3$patient)


# M3 <- NormalizeData(M3)
# M3 <- FindVariableFeatures(M3)
# M3 <- ScaleData(M3)
# M3 <- RunPCA(M3)


# M3 <- IntegrateLayers(
#   object = M3, method = HarmonyIntegration,
#   orig.reduction = "pca", new.reduction = "integrated.cca",
#   verbose = FALSE
# )

# M3 <- FindNeighbors(M3, reduction = "integrated.cca", dims = 1:50)
# M3 <- FindClusters(M3, resolution = 2, cluster.name = "cca_clusters")
# M3 <- RunUMAP(M3, reduction = "integrated.cca", dims = 1:50, reduction.name = "umap.cca")

# p1 <- DimPlot(
#   M3,
#   reduction = "umap.cca",
#   group.by = c("origin"),
#   combine = FALSE, label.size = 2
# )


# pdf(paste0('./output/', 'sn_RV_Peds_CM_Harmonized_UMAP.pdf'), width=5, height=5)
# PlotEmbedding(M3,group.by='origin',reduction="umap.cca",point_size=0.2,plot_under=TRUE,plot_theme=umap_theme()+NoLegend(),raster_dpi=400,raster_scale=0.5)
# dev.off()


# M3 <- subset(M2,CombinedNames=='Myeloid')
# DefaultAssay(M3) <- 'RNA'
# M3[["RNA"]] <- split(M3[["RNA"]], f = M3$patient)


# M3 <- NormalizeData(M3)
# M3 <- FindVariableFeatures(M3)
# M3 <- ScaleData(M3)
# M3 <- RunPCA(M3)


# M3 <- IntegrateLayers(
#   object = M3, method = HarmonyIntegration,
#   orig.reduction = "pca", new.reduction = "integrated.cca",
#   verbose = FALSE
# )

# M3 <- FindNeighbors(M3, reduction = "integrated.cca", dims = 1:50)
# M3 <- FindClusters(M3, resolution = 2, cluster.name = "cca_clusters")
# M3 <- RunUMAP(M3, reduction = "integrated.cca", dims = 1:50, reduction.name = "umap.cca")

# p1 <- DimPlot(
#   M3,
#   reduction = "umap.cca",
#   group.by = c("origin"),
#   combine = FALSE, label.size = 2
# )


# pdf(paste0('./output/', 'sn_RV_Peds_Myeloid_Harmonized_UMAP.pdf'), width=5, height=5)
# PlotEmbedding(M3,group.by='origin',reduction="umap.cca",point_size=0.2,plot_under=TRUE,plot_theme=umap_theme()+NoLegend(),raster_dpi=400,raster_scale=0.5)
# dev.off()



















# DefaultAssay(M3) <- 'RNA'
# M3 <- NormalizeData(M3)
# M3 <- FindVariableFeatures(M3)
# M3 <- ScaleData(M3)
# M3 <- RunPCA(M3)
# M3 <- RunHarmony(M3,'origin')
# M3 <- RunUMAP(M3,reduction = "harmony", dims = 1:50, verbose = F) 

# RV_and_Peds <- intersect(rownames(subset(M3,origin=='RV')),rownames(subset(M3,origin=='Peds')))



# M3 <- SetIdent(M3, value = "origin")

# cm_compare <- FindMarkers(M3,ident.1='RV',ident.2='Peds')


# write.csv(cm_compare,'./output/cm_RV_vs_Peds_DEG.csv')



# M3 <- subset(M2,CombinedNames=='Myeloid')
# DefaultAssay(M3) <- 'RNA'
# M3 <- SetIdent(M3, value = "origin")
# myeloid_compare <- FindMarkers(M3,ident.1='RV',ident.2='Peds')









# RefMerge<-readRDS('./output/Kory_with_RV_modules_projected.rds')




# ####CM
# RefMerge<-SetActiveWGCNA(RefMerge, "CM_RV2LV")

# #RefMerge <- ModuleConnectivity(RefMerge,group.by = 'Names', group_name = 'Cardiomyocytes')


# modules <- GetModules(RefMerge)
# MEs <- GetMEs(RefMerge, T)
# genes_use <- as.character(modules$gene_name)
# params <- GetWGCNAParams(RefMerge)

# # get the assay
# assay <- DefaultAssay(RefMerge)

# cells.use <- RefMerge@meta.data %>% subset(get('Names') %in% 'Cardiomyocytes') %>% rownames
# MEs <- MEs[cells.use,]

# exp_mat <- SeuratObject::LayerData(RefMerge, assay=assay, layer='data')[genes_use,cells.use]

# kMEs <- corSparse(X = Matrix::t(exp_mat),Y = as.matrix(MEs))
# rownames(kMEs) <- genes_use
# kMEs <- as.data.frame(kMEs)

# modules <- modules[,1:3]
# mods <- levels(modules$module)
# colnames(kMEs) <- colnames(MEs)
# kMEs <- kMEs[,mods]
# colnames(kMEs) <- paste0("kME_", colnames(kMEs))
# kMEs <- cbind(modules, kMEs)
# RefMerge <- SetModules(RefMerge, kMEs, "CM_RV2LV")



# RefMerge <- ModuleExprScore(RefMerge,method='Seurat')


# mod_scores <-  GetModuleScores(RefMerge)
# mod_scores <- mod_scores[,colnames(mod_scores) != 'grey']

# # add hMEs to Seurat meta-data:
# RefMerge@meta.data <- cbind(
#   RefMerge@meta.data,
#   mod_scores
# )

# # plot with Seurat's DotPlot function
# p1 <- DotPlot(
#     RefMerge,
#     features = mixedsort(colnames(mod_scores)),
#     group.by = 'Names'
# )

# # flip the x/y axes, rotate the axis labels, and change color scheme:
# p1 <- p1 +
#   RotatedAxis() +
#   scale_color_gradient2(high='darkorchid4', mid='grey95', low='seagreen') +
#   xlab('') + ylab('')


# pdf(paste0('./output/', 'RV2LV_LV_CM_modules_dot.pdf'), width=7.5, height=5)
# p1
# dev.off()


# seurat_ref <- readRDS('./output/scWGCNA_all_celltypes.rds')
# seurat_ref<-SetActiveWGCNA(seurat_ref, "CM")
# seurat_ref <- ModuleConnectivity(seurat_ref,group.by = 'Names', group_name = 'CM')

# seurat_ref <- ModuleExprScore(seurat_ref,method='Seurat')


# mod_scores <-  GetModuleScores(seurat_ref)
# mod_scores <- mod_scores[,colnames(mod_scores) != 'grey']

# # add hMEs to Seurat meta-data:
# seurat_ref@meta.data <- cbind(
#   seurat_ref@meta.data,
#   mod_scores
# )

# # plot with Seurat's DotPlot function
# p2 <- DotPlot(
#     seurat_ref,
#     features = mixedsort(colnames(mod_scores)),
#     group.by = 'Names'
# )

# # flip the x/y axes, rotate the axis labels, and change color scheme:
# p2 <- p2 +
#   RotatedAxis() +
#   scale_color_gradient2(high='darkorchid4', mid='grey95', low='seagreen') +
#   xlab('') + ylab('')


# pdf(paste0('./output/', "LV_CM_modules_dot.pdf"), width=7.5, height=5)
# p2
# dev.off()

# pdf(paste0('./output/', "LV_CM_LV2RM_modules_dot.pdf"), width=7.5, height=8)
# wrap_plots(list(p1,p2),ncol=1)
# dev.off()

# ###EC

# RefMerge<-SetActiveWGCNA(RefMerge, "EC_RV2LV")

# #RefMerge <- ModuleConnectivity(RefMerge,group.by = 'Names', group_name = 'Cardiomyocytes')


# modules <- GetModules(RefMerge)
# MEs <- GetMEs(RefMerge, T)
# genes_use <- as.character(modules$gene_name)
# params <- GetWGCNAParams(RefMerge)

# # get the assay
# assay <- DefaultAssay(RefMerge)

# cells.use <- RefMerge@meta.data %>% subset(get('Names') %in% 'Endothelium') %>% rownames
# MEs <- MEs[cells.use,]

# exp_mat <- SeuratObject::LayerData(RefMerge, assay=assay, layer='data')[genes_use,cells.use]

# kMEs <- corSparse(X = Matrix::t(exp_mat),Y = as.matrix(MEs))
# rownames(kMEs) <- genes_use
# kMEs <- as.data.frame(kMEs)

# modules <- modules[,1:3]
# mods <- levels(modules$module)
# colnames(kMEs) <- colnames(MEs)
# kMEs <- kMEs[,mods]
# colnames(kMEs) <- paste0("kME_", colnames(kMEs))
# kMEs <- cbind(modules, kMEs)
# RefMerge <- SetModules(RefMerge, kMEs, "EC_RV2LV")



# RefMerge <- ModuleExprScore(RefMerge,method='Seurat')


# mod_scores <-  GetModuleScores(RefMerge)
# mod_scores <- mod_scores[,colnames(mod_scores) != 'grey']

# # add hMEs to Seurat meta-data:
# RefMerge@meta.data <- cbind(
#   RefMerge@meta.data,
#   mod_scores
# )

# # plot with Seurat's DotPlot function
# p1 <- DotPlot(
#     RefMerge,
#     features = mixedsort(colnames(mod_scores)),
#     group.by = 'Names'
# )

# # flip the x/y axes, rotate the axis labels, and change color scheme:
# p1 <- p1 +
#   RotatedAxis() +
#   scale_color_gradient2(high='darkorchid4', mid='grey95', low='seagreen') +
#   xlab('') + ylab('')


# pdf(paste0('./output/', 'RV2LV_LV_EC_modules_dot.pdf'), width=7.5, height=5)
# p1
# dev.off()


# seurat_ref<-SetActiveWGCNA(seurat_ref, "EC")
# seurat_ref <- ModuleConnectivity(seurat_ref,group.by = 'Names', group_name = 'EC')

# seurat_ref <- ModuleExprScore(seurat_ref,method='Seurat')


# mod_scores <-  GetModuleScores(seurat_ref)
# mod_scores <- mod_scores[,colnames(mod_scores) != 'grey']

# # add hMEs to Seurat meta-data:
# seurat_ref@meta.data <- cbind(
#   seurat_ref@meta.data,
#   mod_scores
# )

# # plot with Seurat's DotPlot function
# p2 <- DotPlot(
#     seurat_ref,
#     features = mixedsort(colnames(mod_scores)),
#     group.by = 'Names'
# )

# # flip the x/y axes, rotate the axis labels, and change color scheme:
# p2 <- p2 +
#   RotatedAxis() +
#   scale_color_gradient2(high='darkorchid4', mid='grey95', low='seagreen') +
#   xlab('') + ylab('')


# pdf(paste0('./output/', "LV_EC_modules_dot.pdf"), width=7.5, height=5)
# p2
# dev.off()

# pdf(paste0('./output/', "LV_EC_LV2RM_modules_dot.pdf"), width=7.5, height=8)
# wrap_plots(list(p1,p2),ncol=1)
# dev.off()



# ###FB

# RefMerge<-SetActiveWGCNA(RefMerge, "FB_RV2LV")

# #RefMerge <- ModuleConnectivity(RefMerge,group.by = 'Names', group_name = 'Cardiomyocytes')


# modules <- GetModules(RefMerge)
# MEs <- GetMEs(RefMerge, T)
# genes_use <- as.character(modules$gene_name)
# params <- GetWGCNAParams(RefMerge)

# # get the assay
# assay <- DefaultAssay(RefMerge)

# cells.use <- RefMerge@meta.data %>% subset(get('Names') %in% 'Fibroblasts') %>% rownames
# MEs <- MEs[cells.use,]

# exp_mat <- SeuratObject::LayerData(RefMerge, assay=assay, layer='data')[genes_use,cells.use]

# kMEs <- corSparse(X = Matrix::t(exp_mat),Y = as.matrix(MEs))
# rownames(kMEs) <- genes_use
# kMEs <- as.data.frame(kMEs)

# modules <- modules[,1:3]
# mods <- levels(modules$module)
# colnames(kMEs) <- colnames(MEs)
# kMEs <- kMEs[,mods]
# colnames(kMEs) <- paste0("kME_", colnames(kMEs))
# kMEs <- cbind(modules, kMEs)
# RefMerge <- SetModules(RefMerge, kMEs, "FB_RV2LV")



# RefMerge <- ModuleExprScore(RefMerge,method='Seurat')


# mod_scores <-  GetModuleScores(RefMerge)
# mod_scores <- mod_scores[,colnames(mod_scores) != 'grey']

# # add hMEs to Seurat meta-data:
# RefMerge@meta.data <- cbind(
#   RefMerge@meta.data,
#   mod_scores
# )

# # plot with Seurat's DotPlot function
# p1 <- DotPlot(
#     RefMerge,
#     features = mixedsort(colnames(mod_scores)),
#     group.by = 'Names'
# )

# # flip the x/y axes, rotate the axis labels, and change color scheme:
# p1 <- p1 +
#   RotatedAxis() +
#   scale_color_gradient2(high='darkorchid4', mid='grey95', low='seagreen') +
#   xlab('') + ylab('')


# pdf(paste0('./output/', 'RV2LV_LV_FB_modules_dot.pdf'), width=7.5, height=5)
# p1
# dev.off()


# seurat_ref<-SetActiveWGCNA(seurat_ref, "FB")
# seurat_ref <- ModuleConnectivity(seurat_ref,group.by = 'Names', group_name = 'FB')

# seurat_ref <- ModuleExprScore(seurat_ref,method='Seurat')


# mod_scores <-  GetModuleScores(seurat_ref)
# mod_scores <- mod_scores[,colnames(mod_scores) != 'grey']

# # add hMEs to Seurat meta-data:
# seurat_ref@meta.data <- cbind(
#   seurat_ref@meta.data,
#   mod_scores
# )

# # plot with Seurat's DotPlot function
# p2 <- DotPlot(
#     seurat_ref,
#     features = mixedsort(colnames(mod_scores)),
#     group.by = 'Names'
# )

# # flip the x/y axes, rotate the axis labels, and change color scheme:
# p2 <- p2 +
#   RotatedAxis() +
#   scale_color_gradient2(high='darkorchid4', mid='grey95', low='seagreen') +
#   xlab('') + ylab('')


# pdf(paste0('./output/', "LV_FB_modules_dot.pdf"), width=7.5, height=5)
# p2
# dev.off()

# pdf(paste0('./output/', "LV_FB_LV2RM_modules_dot.pdf"), width=7.5, height=8)
# wrap_plots(list(p1,p2),ncol=1)
# dev.off()



# ###Myeloid

# RefMerge<-SetActiveWGCNA(RefMerge, "Myeloid_RV2LV")

# #RefMerge <- ModuleConnectivity(RefMerge,group.by = 'Names', group_name = 'Cardiomyocytes')


# modules <- GetModules(RefMerge)
# MEs <- GetMEs(RefMerge, T)
# genes_use <- as.character(modules$gene_name)
# params <- GetWGCNAParams(RefMerge)

# # get the assay
# assay <- DefaultAssay(RefMerge)

# cells.use <- RefMerge@meta.data %>% subset(get('Names') %in% 'Myeloid') %>% rownames
# MEs <- MEs[cells.use,]

# exp_mat <- SeuratObject::LayerData(RefMerge, assay=assay, layer='data')[genes_use,cells.use]

# kMEs <- corSparse(X = Matrix::t(exp_mat),Y = as.matrix(MEs))
# rownames(kMEs) <- genes_use
# kMEs <- as.data.frame(kMEs)

# modules <- modules[,1:3]
# mods <- levels(modules$module)
# colnames(kMEs) <- colnames(MEs)
# kMEs <- kMEs[,mods]
# colnames(kMEs) <- paste0("kME_", colnames(kMEs))
# kMEs <- cbind(modules, kMEs)
# RefMerge <- SetModules(RefMerge, kMEs, "Myeloid_RV2LV")



# RefMerge <- ModuleExprScore(RefMerge,method='Seurat')


# mod_scores <-  GetModuleScores(RefMerge)
# mod_scores <- mod_scores[,colnames(mod_scores) != 'grey']

# # add hMEs to Seurat meta-data:
# RefMerge@meta.data <- cbind(
#   RefMerge@meta.data,
#   mod_scores
# )

# # plot with Seurat's DotPlot function
# p1 <- DotPlot(
#     RefMerge,
#     features = mixedsort(colnames(mod_scores)),
#     group.by = 'Names'
# )

# # flip the x/y axes, rotate the axis labels, and change color scheme:
# p1 <- p1 +
#   RotatedAxis() +
#   scale_color_gradient2(high='darkorchid4', mid='grey95', low='seagreen') +
#   xlab('') + ylab('')


# pdf(paste0('./output/', 'RV2LV_LV_Myeloid_modules_dot.pdf'), width=7.5, height=5)
# p1
# dev.off()


# seurat_ref<-SetActiveWGCNA(seurat_ref, "Myeloid")
# seurat_ref <- ModuleConnectivity(seurat_ref,group.by = 'Names', group_name = 'Myeloid')

# seurat_ref <- ModuleExprScore(seurat_ref,method='Seurat')


# mod_scores <-  GetModuleScores(seurat_ref)
# mod_scores <- mod_scores[,colnames(mod_scores) != 'grey']

# # add hMEs to Seurat meta-data:
# seurat_ref@meta.data <- cbind(
#   seurat_ref@meta.data,
#   mod_scores
# )

# # plot with Seurat's DotPlot function
# p2 <- DotPlot(
#     seurat_ref,
#     features = mixedsort(colnames(mod_scores)),
#     group.by = 'Names'
# )

# # flip the x/y axes, rotate the axis labels, and change color scheme:
# p2 <- p2 +
#   RotatedAxis() +
#   scale_color_gradient2(high='darkorchid4', mid='grey95', low='seagreen') +
#   xlab('') + ylab('')


# pdf(paste0('./output/', "LV_Myeloid_modules_dot.pdf"), width=7.5, height=5)
# p2
# dev.off()

# pdf(paste0('./output/', "LV_Myeloid_LV2RM_modules_dot.pdf"), width=7.5, height=8)
# wrap_plots(list(p1,p2),ncol=1)
# dev.off()



# #######################################
# #############  FIGURE 7B  #############
# #######################################


# RefMerge<-readRDS('./output/Kory_with_RV_modules_projected.rds')

# RefMerge<-SetActiveWGCNA(RefMerge, "CM_RV2LV")


# plot_list <- PlotModulePreservation(
#   RefMerge,
#   name="CM_pres",
#   statistics = "summary"
# )

# pdf(paste0('./output/', 'RV2LV_CM_pres.pdf'), width=5, height=5)
# wrap_plots(plot_list, ncol=2)
# dev.off()

#######################################
#############  FIGURE 7C  #############
#######################################
#Library EnrichR

#CM not preserved - M1, M5
#CM preserved - M7, M9, M10, M11, M14, M14

#EC not preserved - M1, M2, M4, M5, M10, M11, M12
#EC preserved - M8

# seurat_ref <- readRDS('./output/scWGCNA_all_celltypes.rds')
# seurat_ref<-SetActiveWGCNA(seurat_ref, "CM")


# dbs <- c('GO_Biological_Process_2023','GO_Cellular_Component_2023','GO_Molecular_Function_2023','ChEA_2022','Reactome_2022')

# # perform enrichment tests
# seurat_ref <- RunEnrichr(
#   seurat_ref,
#   dbs=dbs, # character vector of enrichr databases to test
#   max_genes = 2000 # number of genes per module to test. use max_genes = Inf to choose all genes!
# )

# # retrieve the output table
# enrich_df <- GetEnrichrTable(seurat_ref)

# EnrichrBarPlot(
#   seurat_ref,
#   outdir = "./output/sc_enrichr_plot_CM", 
#   n_terms = 5,
#   plot_size = c(5,4), # width, height of the output .pdfs
#   logscale=TRUE # do you want to show the enrichment as a log scale?
# )


# seurat_ref<-SetActiveWGCNA(seurat_ref, "EC")


# dbs <- c('GO_Biological_Process_2023','GO_Cellular_Component_2023','GO_Molecular_Function_2023','ChEA_2022','Reactome_2022')

# # perform enrichment tests
# seurat_ref <- RunEnrichr(
#   seurat_ref,
#   dbs=dbs, # character vector of enrichr databases to test
#   max_genes = 2000 # number of genes per module to test. use max_genes = Inf to choose all genes!
# )

# # retrieve the output table
# enrich_df <- GetEnrichrTable(seurat_ref)

# EnrichrBarPlot(
#   seurat_ref,
#   outdir = "./output/sc_enrichr_plot_EC", 
#   n_terms = 5,
#   plot_size = c(5,4), # width, height of the output .pdfs
#   logscale=TRUE # do you want to show the enrichment as a log scale?
# )


# #######################################
# #############  FIGURE 7D  #############
# #######################################
# library(GeneOverlap)

# seurat_ref <- readRDS('./output/scWGCNA_all_celltypes.rds')
# seurat_ref<-SetActiveWGCNA(seurat_ref, "CM")


# bulk_rv_vs_lv <- read.csv(paste0('./output/', 'RV_LV_align_NRVM_NRMV_ARVM.csv'))
# genes <- toupper(bulk_rv_vs_lv[,2])
# #nrvm <- bulk_rv_vs_lv[2:348,14]
# arvm <- bulk_rv_vs_lv[,20]
# arvm.p <-  bulk_rv_vs_lv[,22]
# arvm_rv <-genes[arvm.p < 0.05 & arvm < 0 & genes %in% rownames(seurat_ref)]
# arvm_lv <- genes[arvm.p < 0.05 & arvm > 0 & genes %in% rownames(seurat_ref)]
# arvm_share <-genes[arvm.p >= 1 & arvm < 0.5 & arvm > -0.5 & genes %in% rownames(seurat_ref)]


# # load modules
# modules <- GetModules(seurat_ref)
# mods <- levels(modules$module)
# genome.size <- nrow(modules)


# overlap_df_rv <- do.call(rbind, lapply(mods, function(cur_mod){

#   cur_genes <- modules %>% subset(module == cur_mod) %>% .$gene_name

#   cur_overlap <- testGeneOverlap(newGeneOverlap(
#       cur_genes,
#       arvm_rv,
#       genome.size=genome.size
#   ))

#   cur_overlap <- data.frame(
#     'odds.ratio' = cur_overlap@odds.ratio,
#     'pval' = cur_overlap@pval,
#     'Jaccard' = cur_overlap@Jaccard,
#     'size_intersection' = length(cur_overlap@intersection),
#     'module' = cur_mod
#   )

#   cur_overlap

# })) %>% as.data.frame()

# overlap_df_lv <- do.call(rbind, lapply(mods, function(cur_mod){

#   cur_genes <- modules %>% subset(module == cur_mod) %>% .$gene_name

#   cur_overlap <- testGeneOverlap(newGeneOverlap(
#       cur_genes,
#       arvm_lv,
#       genome.size=genome.size
#   ))

#   cur_overlap <- data.frame(
#     'odds.ratio' = cur_overlap@odds.ratio,
#     'pval' = cur_overlap@pval,
#     'Jaccard' = cur_overlap@Jaccard,
#     'size_intersection' = length(cur_overlap@intersection),
#     'module' = cur_mod
#   )

#   cur_overlap

# })) %>% as.data.frame()

# overlap_df_share <- do.call(rbind, lapply(mods, function(cur_mod){

#   cur_genes <- modules %>% subset(module == cur_mod) %>% .$gene_name

#   cur_overlap <- testGeneOverlap(newGeneOverlap(
#       cur_genes,
#       arvm_share,
#       genome.size=genome.size
#   ))

#   cur_overlap <- data.frame(
#     'odds.ratio' = cur_overlap@odds.ratio,
#     'pval' = cur_overlap@pval,
#     'Jaccard' = cur_overlap@Jaccard,
#     'size_intersection' = length(cur_overlap@intersection),
#     'module' = cur_mod
#   )

#   cur_overlap

# })) %>% as.data.frame()

# overlap_df_rv <- overlap_df_rv %>% mutate(fdr=p.adjust(pval, method='fdr'))
# overlap_df_rv <- overlap_df_rv %>% subset(module != 'grey')



# overlap_df_lv <- overlap_df_lv %>% mutate(fdr=p.adjust(pval, method='fdr'))
# overlap_df_lv <- overlap_df_lv %>% subset(module != 'grey')

# overlap_df_share <- overlap_df_share %>% mutate(fdr=p.adjust(pval, method='fdr'))
# overlap_df_share <- overlap_df_share %>% subset(module != 'grey')

# # Plot as a lollipop

# overlap_df_rv$shape <- ifelse(overlap_df_rv$fdr < 0.05, 21, 4)
# overlap_df_rv <- overlap_df_rv %>% arrange(odds.ratio, descending=TRUE)
# overlap_df_rv$module <- factor(as.character(overlap_df_rv$module), levels=as.character(overlap_df_rv$module))

# mod_colors <- dplyr::select(modules, c(module, color)) %>%
#   distinct
# cp <- mod_colors$color; names(cp) <- mod_colors$module

# p <- overlap_df_rv %>%
#   ggplot(aes(y=module, x=odds.ratio, size= size_intersection, color=module)) +
#   geom_segment(aes(y=module, yend=module, x=0, xend=odds.ratio), size=0.5, color='grey') +
#   geom_point() +
#   geom_point(shape=overlap_df_rv$shape, color='black', fill=NA) +
#   scale_color_manual(values=cp, guide='none') +
#   ylab('') + xlab("Odds ratio") +
#   scale_x_continuous(breaks = c(0, 1, 2,3)) +
#   labs(size='Size\nintersection') +
#   ggtitle('Overlap with SFARI genes') +

#   theme(
#     panel.border = element_rect(size=1, color='black', fill=NA),
#     axis.line.y = element_blank(),
#     axis.line.x = element_blank(),
#     plot.title = element_text(hjust=0.5, face='plain')
#   )

# pdf(paste0('./output/', 'RVbulk_R2LV_CM_overlap.pdf'), width=4.5, height=5)
# p
# dev.off()


# #######################################
# #############  FIGURE 7D  #############
# #######################################
# rm(RefMerge)
# rm(seurat_ref)
# seurat_ref <- readRDS('./dependencies/shared/Post_R3_FINAL_with_counts.rds')
# load('./output/GSE183852_DCM_Integrated.Robj')


# sc_modules <- read.csv(paste0('./output/', 'sc_heart_modules.csv'))
# seurat_ref <- AddModuleScore(seurat_ref,list(
#   subset(sc_modules,module=="CM-M1")$gene_name,
#   subset(sc_modules,module=="CM-M2")$gene_name,
#   subset(sc_modules,module=="CM-M3")$gene_name,
#   subset(sc_modules,module=="CM-M4")$gene_name,
#   subset(sc_modules,module=="CM-M5")$gene_name,
#   subset(sc_modules,module=="CM-M6")$gene_name,
#   subset(sc_modules,module=="CM-M7")$gene_name
#   ),name='CM_score')
# RefMerge <- AddModuleScore(RefMerge,list(
#   subset(sc_modules,module=="CM-M1")$gene_name,
#   subset(sc_modules,module=="CM-M2")$gene_name,
#   subset(sc_modules,module=="CM-M3")$gene_name,
#   subset(sc_modules,module=="CM-M4")$gene_name,
#   subset(sc_modules,module=="CM-M5")$gene_name,
#   subset(sc_modules,module=="CM-M6")$gene_name,
#   subset(sc_modules,module=="CM-M7")$gene_name
#   )
#   ,name='CM_score')


# FeaturePlot(RefMerge,'CM_score1')
# DotPlot(RefMerge,c('CM_score1','CM_score2','CM_score3','CM_score4','CM_score5','CM_score6','CM_score7'),group.by='Names',split.by='condition')
# DotPlot(seurat_ref,c('CM_score1','CM_score2','CM_score3','CM_score4','CM_score5','CM_score6','CM_score7'),group.by='Names',split.by='group')


# cellsofint<- colnames(RefMerge)[RefMerge@meta.data$'CM-M6'>5 & RefMerge@meta.data$'EC-M7'>5 & RefMerge@meta.data$'CM-M11'>5]
# DimPlot(RefMerge, cells.highlight=cellsofint)
# FindMarkers(RefMerge,ident.1=cellsofint)



# #cellsofint<- colnames(seurat_ref)[seurat_ref@meta.data$'CM-M6'>0 & seurat_ref@meta.data$'EC-M7'>0 & seurat_ref@meta.data$'CM-M11'>0]
# cellsofint<- colnames(seurat_ref)[seurat_ref@meta.data$'CM-M6'>10]

# cellstocomp<- colnames(seurat_ref)[seurat_ref@meta.data$'CM-M6'<10 & seurat_ref@meta.data$Names == "CM"]

# DimPlot(seurat_ref, cells.highlight=cellsofint,raster=T)
# FindMarkers(seurat_ref,ident.1=cellsofint,ident.2=cellstocomp)




