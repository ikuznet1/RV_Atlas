library(Seurat)
library(hdWGCNA)
library(ggeasy)
library(harmony)



source('~/Downloads/hdWGCNA_TOM/spatial_functions.R')



#######################################
#############  FIGURE 5A  #############
#######################################

M1 <- readRDS(file = "/Volumes/Extreme\ SSD/Final_Analysis/CellTypes/myeloid_subclust.rds")

M1 <- FindNeighbors(M1)
M1 <- FindClusters(M1,resolution=1)

new.cluster.ids <- c('CCR2- rMac2','CCR2- rMac1','CCR2- rMac2','CCR2- rMac2',
	'CCR2+ rMac','CCR2- rMac1','CCR2- rMac1',
	'CCR2- rMac1','DCs','iMac','TREM2 Mac','CCR2- rMac2')
names(new.cluster.ids) <- levels(M1)
M1 <- RenameIdents(M1, new.cluster.ids)

M2 <- subset(M1,idents=c('CCR2+ rMac','DCs'))
M2<-RunPCA(M2)
M2 <- M2%>% 
  RunUMAP(reduction = "pca", dims = 1:50, verbose = F) %>% 
  FindNeighbors(reduction = "pca", dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()

#0 is CD1c- CD163+ MERTK+
#1 is CCR2+ CD11c+
#2 is DC
#3 is mix of classical (FCN1+), non-classical (LILRB2, ITGAL, FCGR3A), and neutrophil-like (CSF3R)

M1$Subnames<-M1@active.ident
M1$Subsubnames <- M1$Subnames
cells_DC = colnames(M2)[M2$seurat_clusters==2]
cells_mono = colnames(M2)[M2$seurat_clusters==3]
cells_rMac = union(colnames(M2)[M2$seurat_clusters==0],colnames(M2)[M2$seurat_clusters==1])

levels(M1$Subsubnames) <- c(levels(M1$Subsubnames), 'Mono')
M1$Subsubnames[colnames(M1) %in% cells_DC] = 'DCs'
M1$Subsubnames[colnames(M1) %in% cells_mono] = 'Mono'
M1$Subsubnames[colnames(M1) %in% cells_rMac] = 'CCR2+ rMac'


M2 <- M1
M2<-RunPCA(M2)
M2 <- RunHarmony(M1,'patient')
M2 <- M2%>% 
  RunUMAP(reduction = "harmony", dims = 1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()

M1 <- M2
M1 <- SetIdent(M1,value='Subsubnames')

M1 <- AddModuleScore(M1, features=list(c("TREM2","GPNMB","MITF","SPP1")),assay="SCT",name="TREM2_Mac_Score")
M1 <- AddModuleScore(M1, features=list(c("CLEC9A","ZBTB46","CD1C","CD226")),assay="SCT",name="DC_Score")
M1 <- AddModuleScore(M1, features=list(c("FCN1","LILRB2","ITGAL","CSF3R")),assay="SCT",name="Mono_Score")
M1 <- AddModuleScore(M1, features=list(c("CCR2","CX3CR1","ITGAX")),assay="SCT",name="CCR2+_rMac_Score")
M1 <- AddModuleScore(M1, features=list(c("LYVE1","FOLR2","SIGLEC1","F13A1")),assay="SCT",name="CCR2-_rMac1_Score")
M1 <- AddModuleScore(M1, features=list(c("RBMS3","PLA2G5","EBF1")),assay="SCT",name="CCR2-_rMac2_Score")
M1 <- AddModuleScore(M1, features=list(c("IL1B","CCL3","CCL4","CXCL3","CXCL8")),assay="SCT",name="iMac_Score")



# 0, 2,3,10,11 mystery


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'Myeloid_snUMAP.pdf'), width=5, height=5)
PlotEmbedding(M1,group.by='Subsubnames',point_size=1,plot_under=TRUE,plot_theme=umap_theme()+NoLegend(),raster_dpi=400,raster_scale=0.5)
dev.off()

#######################################
#############  FIGURE 5B  #############
#######################################

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'Myeloid_dot.pdf'), width=6, height=5)
DotPlot(M1, c("Mono_Score1","TREM2_Mac_Score1","iMac_Score1",
	"DC_Score1","CCR2+_rMac_Score1","CCR2-_rMac1_Score1",
	"CCR2-_rMac2_Score1"))+xlab('Marker Score')+
	scale_x_discrete(labels=c('1','2','3','4','5','6','7'))
dev.off()

M1$Names_group <- paste0(M1$Subsubnames,'_',M1$group)


table(M1$group,M1$Subsubnames)/rowSums(table(M1$group,M1$Subsubnames))

NF_percent_cell <- cbind(as.data.frame(table(subset(M1,group=="NF")@active.ident)/length(subset(M1,group=="NF")@active.ident)*100),type = "NF")
NF_percent_cell$sum <- (rev(cumsum(rev(NF_percent_cell$Freq))) - NF_percent_cell$Freq/2)/100
NF_percent_cell$Freq <- NF_percent_cell$Freq/100


pRV_percent_cell <- cbind(as.data.frame(table(subset(M1,group=="pRV")@active.ident)/length(subset(M1,group=="pRV")@active.ident)*100),type = "pRV")
pRV_percent_cell$sum <- (rev(cumsum(rev(pRV_percent_cell$Freq))) - pRV_percent_cell$Freq/2)/100
pRV_percent_cell$Freq <- pRV_percent_cell$Freq/100


RVF_percent_cell <- cbind(as.data.frame(table(subset(M1,group=="RVF")@active.ident)/length(subset(M1,group=="RVF")@active.ident)*100),type = "RVF")
RVF_percent_cell$sum <- (rev(cumsum(rev(RVF_percent_cell$Freq))) - RVF_percent_cell$Freq/2)/100
RVF_percent_cell$Freq <- RVF_percent_cell$Freq/100

percent_cell_df <- rbind(NF_percent_cell,pRV_percent_cell,RVF_percent_cell)
pdf('~/Downloads/hdWGCNA_TOM/Myeloid_prev_stacked.pdf',width=6,height=2.5)
ggplot(percent_cell_df, aes(fill=Var1, y=Freq, x=type,label=round(sum,1))) +  
geom_bar(position="stack", stat="identity",width=0.6) + theme_classic() + coord_flip()+
xlab("Disease State") + ylab("Frequency") + labs(fill="Cell type",color='black') + theme(text = element_text(size=20),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),legend.text=element_text(color="black")) + scale_y_continuous(expand=c(0,0)) + geom_label_repel(aes(type,sum,label=scales::percent(round(Freq,2))),fill=NA,nudge_x=0.5,direction="y")
dev.off()

#######################################
#############  FIGURE 5C #############
#######################################

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'Myeloid_dot_disease.pdf'), width=6, height=2.5)

DotPlot(M1, c("Mono_Score1","TREM2_Mac_Score1","iMac_Score1",
	"DC_Score1","CCR2+_rMac_Score1","CCR2-_rMac1_Score1",
	"CCR2-_rMac2_Score1"),col.min=0,col.max=1,group.by='group')+xlab('Marker Score') +
	scale_x_discrete(labels=c('Mono','TREM2 Mac','iMac','DCs','CCR2+ rMac','CCR2- rMac1','CCR2- rMac2'))+
  scale_color_gradient2(high='red', mid='grey95', low='blue')+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()



#######################################
#############  FIGURE 5D  #############
#######################################
library(monocle3)
library(SeuratWrappers)
cds <- as.cell_data_set(M1)
cds <- cluster_cells(cds, resolution=1e-3)

p1 <- plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
wrap_plots(p1, p2)


cds <- learn_graph(cds, use_partition = TRUE, verbose = FALSE)
plot_cells(cds,
           color_cells_by = "cluster",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

cds <- order_cells(cds, root_cells = colnames(cds[,clusters(cds) == 6]))
plot_cells(cds,
           color_cells_by = "pseudotime",
           group_cells_by = "cluster",
           label_cell_groups = FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE,
           trajectory_graph_color = "grey60")
pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'Myeloid_rMac_pseudo.pdf'), width=5, height=5)
M1$pseudotime <- pseudotime(cds)
FeaturePlot(M1,'pseudotime')
dev.off()

mono_cells = colnames(M1)[M1$Subsubnames == 'Mono']
cds <- order_cells(cds, root_cells = mono_cells)
plot_cells(cds,
           color_cells_by = "pseudotime",
           group_cells_by = "cluster",
           label_cell_groups = FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE,
           trajectory_graph_color = "grey60")
pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'Myeloid_mono_pseudo.pdf'), width=5, height=5)
M1$pseudotime <- pseudotime(cds)
FeaturePlot(M1,'pseudotime')
dev.off()

cds <- order_cells(cds, root_cells = union(colnames(cds[,clusters(cds) == 6]),mono_cells))
plot_cells(cds,
           color_cells_by = "pseudotime",
           group_cells_by = "cluster",
           label_cell_groups = FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE,
           trajectory_graph_color = "grey60")
pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'Myeloid_pseudo.pdf'), width=5, height=5)
M1$pseudotime <- pseudotime(cds)
FeaturePlot(M1,'pseudotime')
dev.off()
#######################################
#############  FIGURE 5E  #############
#######################################

seurat_ref <- readRDS('~/Downloads/hdWGCNA_TOM/scWGCNA_bulk2sn_projection.rds')


seurat_ref <- SetActiveWGCNA(seurat_ref , 'bulk2sn')


mapping <- labels2colors(1:100)
MEs <- GetMEs(seurat_ref, harmonized=TRUE)
mods <- colnames(MEs); mods <- mods[mods != 'grey']
mods_num <- paste0('M',match(mods,mapping))




colnames(MEs)<-paste0('M',match(colnames(MEs),mapping))
seurat_ref@meta.data <- cbind(seurat_ref@meta.data, MEs)
seurat_ref <- SetIdent(seurat_ref, value = "Names")


consensus_modules <- read.csv("~/Downloads/hdWGCNA_TOM/bulk_heart_modules.csv")
consensus_modules <- consensus_modules[,1:3]
consensus_modules <- subset(consensus_modules, gene_name %in% rownames(seurat_ref))
# remove duplicate gene names
consensus_modules <- consensus_modules[match(unique(consensus_modules$gene_name), consensus_modules$gene_name),]


library(dplyr)
score_calc <- consensus_modules %>% group_by(module) %>% group_split()
module_colors <- unique(unlist(lapply(score_calc,'[[','module')))
module_colors <- paste0('M',match(module_colors,mapping))


rm(seurat_ref)
gc()
seurat_ref<-M1
#seurat_ref<-readRDS('~/Downloads/hdWGCNA_TOM/RV_data.rds')


idx <- rownames(MEs) %in% colnames(seurat_ref)
seurat_ref$Subnames <- seurat_ref@active.ident
seurat_ref$SubNames_Groups <- paste(seurat_ref$Subnames,seurat_ref$group,sep='_')

seurat_ref <- SetIdent(seurat_ref, value = "Subnames")



seurat_ref@meta.data <- cbind(seurat_ref@meta.data, MEs[idx,])



seurat_ref <- AddModuleScore(seurat_ref,lapply(score_calc,'[[','gene_name'),name="module_score")


cols_current <- colnames(seurat_ref@meta.data)
cols_current[startsWith(colnames(seurat_ref@meta.data),'module_score')] <- paste0('module_',module_colors)
colnames(seurat_ref@meta.data) <- cols_current

#Dot Plot of enrichment cell type of myeloid enriched modules
modules_up <- c('M3','M4',"M8")
modules_down <- c('M1')

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'Myeloid_dot_subclust_up.pdf'), width=5, height=5)

p <- DotPlot(seurat_ref,paste0('module_',modules_up),group.by='Subnames',dot.min=0,col.min=-1,col.max=1,scale.min=50,scale.max=100) +
  RotatedAxis() + ylab('')+ xlab('')+
  scale_color_gradient2(high='red', mid='grey95', low='blue') +
  theme(
    panel.border = element_rect(size=1,fill=NA, color='black'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()
) 
p
dev.off()

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'Myeloid_dot_subclust_down.pdf'), width=5, height=5)

p <- DotPlot(seurat_ref,paste0('module_',modules_down),group.by='Subnames',dot.min=0.5,col.min=-1,col.max=1,scale.min=50,scale.max=100) +
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
#############  FIGURE 5F  #############
#######################################


seurat_obj <- readRDS('~/Downloads/RV_bulkRNASeq_seurat.rds')

modules <- GetModules(seurat_obj)
#color_df <- modules %>% subset(module!='grey') %>%
#  select(c(module, color)) %>% distinct %>%
#  rename(c(group=module, colour=color))
#mods <- levels(modules$module)
#mods <- mods[mods!='grey']

mods <- mapping[c(1,3,4,8)]



#Enrichments
library(enrichR)

dbs <-c('GO_Biological_Process_2023','GO_Cellular_Component_2023','GO_Molecular_Function_2023','Reactome_2022', 'ChEA_2022',"LINCS_L1000_Chem_Pert_up")

# compute GO terms:
enrich_list <- list()
seurat_obj <- RunEnrichr(seurat_obj, dbs=dbs,max_genes = 10000)



# helper function to wrap text
wrapText <- function(x, len) {
    sapply(x, function(y) paste(strwrap(y, len), collapse = "\n"), USE.NAMES = FALSE)
}

combined_output <- GetEnrichrTable(seurat_obj)
selected_terms <- subset(combined_output,db=="Reactome_2022")
selected_terms <- subset(selected_terms,module %in% mapping[c(1,3,4,8)])


# subset selected terms
selected_terms <- subset(selected_terms, P.value < 0.05)
idx_top_1 <- match(unique(selected_terms$module),selected_terms$module)
idx_top_3 <- sort(c(idx_top_1,idx_top_1+1,idx_top_1+2))

selected_terms<-selected_terms[idx_top_3,]


selected_terms$group <- factor(
  as.character(selected_terms$module),
  levels = mods
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


selected_terms$wrap <- wrapText(selected_terms$Term, 35)

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
  RotatedAxis() + xlab('') + ylab('') + coord_flip()
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


color_df <- modules %>% subset(module!='grey') %>%
  select(c(module, color)) %>% distinct
colnames(color_df) = c('group','colour')
color_df <- subset(color_df,colour %in% mods)


color_df$group<-paste0('M',match(color_df$group,mapping))

color_df$group <- factor(
  as.character(color_df$group),
  levels = paste0('M',match(mods,mapping))
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
  NoLegend() + RotatedAxis() + coord_flip()+
  theme(
    plot.title=element_blank(),
    axis.line=element_blank(),
    axis.ticks.y =element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    plot.margin=margin(0,0,0,0),
  )


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'Myeloid_reactome_terms.pdf'), width=10, height=6)
colorbar + p
dev.off()


#Run enrichment by cell type
Idents(seurat_ref) <- "Names_group"
combined_set <- data.frame()
combined_output <- data.frame()

mods_idx <- c(1,3,4,8)
cell_types <- unique(seurat_ref$Subsubnames)
comparison <- list(c("RVF","NF"),c("RVF","pRV"))
for (i in mods_idx){
	for (j in cell_types){
		for (k in comparison){
			key_genes <- subset(bulk_modules,module %in% c(i))$gene_name
			key_genes <- key_genes[key_genes %in% rownames(seurat_ref)]

			gene_set <- FindMarkers(seurat_ref, ident.1 = paste0(j,"_",k[1]), ident.2 = paste0(j,"_",k[2]),features=key_genes)
			
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
selected_terms <- subset(combined_output,db=="Reactome_2022")
selected_terms <- subset(selected_terms,direction=="up")
selected_terms <- subset(selected_terms,comparison=="RVF_NF")
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


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'Myeloid_reactome_terms_cell_type_up_RVF_vs_NF.pdf'), width=8, height=7)
p / colorbar 
dev.off()

#Down
selected_terms <- subset(combined_output,db=="Reactome_2022")
selected_terms <- subset(selected_terms,direction=="down")
selected_terms <- subset(selected_terms,comparison=="RVF_NF")

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


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'Myeloid_reactome_terms_cell_type_down_RVF_vs_NF.pdf'), width=6.6, height=5)
p / colorbar 
dev.off()


###RVF vs pRV

#Up
selected_terms <- subset(combined_output,db=="Reactome_2022")
selected_terms <- subset(selected_terms,direction=="up")
selected_terms <- subset(selected_terms,comparison=="RVF_pRV")
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


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'Myeloid_reactome_terms_cell_type_up_RVF_vs_pRV.pdf'), width=8, height=7)
p / colorbar 
dev.off()

#Down
selected_terms <- subset(combined_output,db=="Reactome_2022")
selected_terms <- subset(selected_terms,direction=="down")
selected_terms <- subset(selected_terms,comparison=="RVF_pRV")

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


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'Myeloid_reactome_terms_cell_type_down_RVF_vs_pRV.pdf'), width=6.6, height=5)
p / colorbar 
dev.off()

#######################################
#############  FIGURE 5G  #############
#######################################


###RVF vs NF

###Chea
#Up
selected_terms <- subset(combined_output,db=="ChEA_2022")
selected_terms <- subset(selected_terms,direction=="up")
selected_terms <- subset(selected_terms,comparison=="RVF_NF")
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


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'Myeloid_chea_terms_cell_type_up_RVF_vs_NF.pdf'), width=10, height=5)
p / colorbar 
dev.off()

#Down
selected_terms <- subset(combined_output,db=="ChEA_2022")
selected_terms <- subset(selected_terms,direction=="down")
selected_terms <- subset(selected_terms,comparison=="RVF_NF")

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


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'Myeloid_chea_terms_cell_type_down_RVF_vs_NF.pdf'), width=6, height=4)
p / colorbar 
dev.off()


###RVF vs pRV
###Chea
#Up
selected_terms <- subset(combined_output,db=="ChEA_2022")
selected_terms <- subset(selected_terms,direction=="up")
selected_terms <- subset(selected_terms,comparison=="RVF_pRV")
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


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'Myeloid_chea_terms_cell_type_up_RVF_vs_pRV.pdf'), width=10, height=5)
p / colorbar 
dev.off()

#Down
selected_terms <- subset(combined_output,db=="ChEA_2022")
selected_terms <- subset(selected_terms,direction=="down")
selected_terms <- subset(selected_terms,comparison=="RVF_pRV")

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


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'Myeloid_chea_terms_cell_type_down_RVF_vs_pRV.pdf'), width=6, height=4)
p / colorbar 
dev.off()

#######################################
#############  FIGURE 5H  #############
#######################################

dbs <-c('GO_Biological_Process_2023','GO_Cellular_Component_2023','GO_Molecular_Function_2023','Reactome_2022', 'ChEA_2022',"LINCS_L1000_Chem_Pert_Consensus_Sigs")

#Run enrichment by module
Idents(seurat_ref) <- "group"
combined_set <- data.frame()
combined_output <- data.frame()

mods_idx <- c(1,3,4,8)
cell_types <- unique(seurat_ref$Subsubnames)
comparison <- list(c("RVF","NF"),c("RVF","pRV"))
for (i in mods_idx){
	for (k in comparison){
		key_genes <- subset(bulk_modules,module %in% c(i))$gene_name
		key_genes <- key_genes[key_genes %in% rownames(seurat_ref)]

		gene_set <- FindMarkers(seurat_ref, ident.1 = paste0(k[1]), ident.2 = paste0(k[2]),features=key_genes)
		
		gene_set<-subset(gene_set,p_val_adj<0.05)
		if (length(rownames(gene_set))==0){next}
		gene_set$module <- paste0('M',i)
		gene_set$color <- mapping[i]
		gene_set$comparison <- paste0(k[1],'_',k[2])

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
	    		cur_df$comparison <- paste0(k[1],'_',k[2])
	    		cur_df$color <- mapping[i]
	    		cur_df$direction <- 'up'
	    		combined_output <- rbind(combined_output, cur_df)
	  		}
		}
	}
}



outdir = '~/Downloads/hdWGCNA_TOM/scMyeloid_subclust_enrichr_plot'
#dir.create("~/Downloads/hdWGCNA_TOM/scMyeloid_subclust_enrichr_plot")


enrichr_df <- combined_output

for (i in 1:length(unique(enrichr_df$module))) {
	for (j in unique(enrichr_df$comparison)){
		for (k in unique(enrichr_df$direction)){
		    cur_mod <- unique(enrichr_df$module)[i]
		    cur_terms <- subset(enrichr_df, module == cur_mod & comparison == j & direction == k)
		    print(cur_mod)

		    if (nrow(cur_terms) == 0) {
		        next
		    }
		    cur_terms$wrap <- wrapText(cur_terms$Term, 45)
		    plot_list <- list()
		    for (cur_db in dbs) {
		        plot_df <- subset(cur_terms, db == cur_db) %>% top_n(5, 
		            wt = Combined.Score)
		        text_color = "black"
		        cur_color = 'green'

		        plot_df$Combined.Score <- log(plot_df$Combined.Score)
		        lab <- "Enrichment log(combined score)"
		        x <- 0.2

		        plot_list[[cur_db]] <- ggplot(plot_df, aes(x = Combined.Score, 
		            y = reorder(wrap, Combined.Score))) + geom_bar(stat = "identity", 
		            position = "identity", color = "white", fill = cur_color) + 
		            geom_text(aes(label = wrap), x = x, color = text_color, 
		              size = 3.5, hjust = "left") + ylab("Term") + 
		            xlab(lab) + ggtitle(cur_db) + theme(panel.grid.major = element_blank(), 
		            panel.grid.minor = element_blank(), legend.title = element_blank(), 
		            axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
		            plot.title = element_text(hjust = 0.5))
		    }
		    pdf(paste0(outdir, "/", cur_mod, "_",j,'_',k,".pdf"), width = 5, 
		        height = 4)
		    for (plot in plot_list) {
		        print(plot)
		    }
		    dev.off()
		}
	}
}


#LINCs only
enrichr_df <- combined_output

for (i in 1:length(unique(enrichr_df$module))) {
	for (j in unique(enrichr_df$comparison)){
		for (k in unique(enrichr_df$direction)){
		    cur_mod <- unique(enrichr_df$module)[i]
		    cur_terms <- subset(enrichr_df, module == cur_mod & comparison == j & direction == k & db == "LINCS_L1000_Chem_Pert_Consensus_Sigs")
		    print(cur_mod)

		    if (nrow(cur_terms) == 0) {
		        next
		    }
		    cur_terms$wrap <- wrapText(cur_terms$Term, 45)
		    plot_list <- list()
		    for (cur_db in c('LINCS_L1000_Chem_Pert_Consensus_Sigs')) {
		        plot_df <- subset(cur_terms, db == cur_db) %>% top_n(5, 
		            wt = Combined.Score)
		        text_color = "black"
		        cur_color = 'green'

		        plot_df$Combined.Score <- log(plot_df$Combined.Score)
		        lab <- "Enrichment log(combined score)"
		        x <- 0.2

		        plot_list[[cur_db]] <- ggplot(plot_df, aes(x = Combined.Score, 
		            y = reorder(wrap, Combined.Score))) + geom_bar(stat = "identity", 
		            position = "identity", color = "white", fill = cur_color) + 
		            geom_text(aes(label = wrap), x = x, color = text_color, 
		              size = 3.5, hjust = "left") + ylab("Term") + 
		            xlab(lab) + ggtitle(cur_db) + theme(panel.grid.major = element_blank(), 
		            panel.grid.minor = element_blank(), legend.title = element_blank(), 
		            axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
		            plot.title = element_text(hjust = 0.5))
		    }
		    pdf(paste0(outdir, "/LINCS_", cur_mod, "_",j,'_',k,".pdf"), width = 2.5, 
		        height = 3)
		    for (plot in plot_list) {
		        print(plot)
		    }
		    dev.off()
		}
	}
}



#######################################
#############  FIGURE 5I  #############
#######################################

#RVF vs NF
Idents(seurat_ref) <- "group"
bulk_modules <- consensus_modules
bulk_modules$module <- match(consensus_modules$module,mapping)

combined_set <- data.frame()
mods_idx <- c(1,3,4,8)
for (i in mods_idx){
	key_genes <- subset(bulk_modules,module %in% c(i))$gene_name
	key_genes <- key_genes[key_genes %in% rownames(seurat_ref)]
	gene_set <- FindMarkers(seurat_ref, ident.1 = "RVF", ident.2 = "NF",features=key_genes)

	gene_set<-subset(gene_set,p_val_adj<0.05)
	if (length(rownames(gene_set))==0){next}
	gene_set$module <- paste0('M',i)
	gene_set$color <- mapping[i]
	if (length(combined_set) == 0){
		combined_set <- gene_set
	}
	else {
		combined_set <- rbind(combined_set,gene_set)
	}
}

library(EnhancedVolcano)

combined_set <- combined_set[!grepl('MT-',rownames(combined_set)),]
keyvals <- combined_set$color
names(keyvals)  <- combined_set$module

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'Myeloid_module_volcano_all_RVF_vs_NF.pdf'), width=8, height=6)

EnhancedVolcano(combined_set,lab=rownames(combined_set),
	x='avg_log2FC',y='p_val_adj',
	FCcutoff = 0.1,
	colCustom=keyvals) + coord_flip()
dev.off()


sub_set <- sub_set[!grepl('MT-',rownames(sub_set)),]
sub_set <- subset(combined_set,combined_set$module != 'M1')
keyvals <- sub_set$color
names(keyvals)  <- sub_set$module

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'Myeloid_module_volcano_no_M1_all_RVF_vs_NF.pdf'), width=8, height=6)

EnhancedVolcano(sub_set,lab=rownames(sub_set),
	x='avg_log2FC',y='p_val_adj',
	FCcutoff = 0.1,
	colCustom=keyvals) + coord_flip()
dev.off()

subsub_set <- subset(sub_set,sub_set$module == 'M8' )
keyvals <- subsub_set$color
names(keyvals)  <- subsub_set$module

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'Myeloid_module_volcano_M8_only_RVF_vs_NF.pdf'), width=8, height=6)

EnhancedVolcano(subsub_set,lab=rownames(subsub_set),
	x='avg_log2FC',y='p_val_adj',
	FCcutoff = 0.1,
	colCustom=keyvals) + coord_flip(xlim=c(-5, 5)) 
dev.off()

#RVF vs pRV
combined_set <- data.frame()
mods_idx <- c(1,3,4,8)
for (i in mods_idx){
	key_genes <- subset(bulk_modules,module %in% c(i))$gene_name
	key_genes <- key_genes[key_genes %in% rownames(seurat_ref)]
	gene_set <- FindMarkers(seurat_ref, ident.1 = "RVF", ident.2 = "pRV",features=key_genes)
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

library(EnhancedVolcano)

combined_set <- combined_set[!grepl('MT-',rownames(combined_set)),]
keyvals <- combined_set$color
names(keyvals)  <- combined_set$module

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'Myeloid_module_volcano_all_RVF_vs_pRV.pdf'), width=8, height=6)

EnhancedVolcano(combined_set,lab=rownames(combined_set),
	x='avg_log2FC',y='p_val_adj',
	FCcutoff = 0.1,
	colCustom=keyvals) + coord_flip()
dev.off()


sub_set <- sub_set[!grepl('MT-',rownames(sub_set)),]
sub_set <- subset(combined_set,combined_set$module != 'M1')
keyvals <- sub_set$color
names(keyvals)  <- sub_set$module

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'Myeloid_module_volcano_no_M1_all_RVF_vs_pRV.pdf'), width=8, height=6)

EnhancedVolcano(sub_set,lab=rownames(sub_set),
	x='avg_log2FC',y='p_val_adj',
	FCcutoff = 0.1,
	colCustom=keyvals) + coord_flip()
dev.off()

subsub_set <- subset(sub_set,sub_set$module == 'M8' )
keyvals <- subsub_set$color
names(keyvals)  <- subsub_set$module

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'Myeloid_module_volcano_M8_only_RVF_vs_pRV.pdf'), width=8, height=6)

EnhancedVolcano(subsub_set,lab=rownames(subsub_set),
	x='avg_log2FC',y='p_val_adj',
	FCcutoff = 0.1,
	colCustom=keyvals) + coord_flip(xlim=c(-5, 5)) 
dev.off()




#######################################
#############  FIGURE 5I  #############
#######################################


#NR3C1 gene_set
#"VIT;VKORC1L1;ERRFI1;AHCYL1;STEAP4;TRAF3IP2;GALNT15;SERPINE1;JADE1;SLA;CBLB;MT1X;EPS8;CCND3;BMPER;RASSF4;RPS6KA2;ANPEP;C1RL;MAP3K6;IL6R;PDGFRA;MLIP;SCARA5;IL1R1;EBF1;TTC7A;CRISPLD2;SPARCL1;FKBP5;NNMT;LPAR1;SLC1A3;PLA2G5;NID1;ACACB;ZFP36L2;PIK3R5;C3;SCFD2;LPXN;HACL1;SRGAP2;SLC38A2;SLC19A2;S100A10;KLHL29;GADD45B;ZBTB16;ELL2;CORO2B;IGF2R;NFATC4;DERA;SULT1B1;MAFB;BCL6;TMEM236;TBXAS1;NDUFAF2;RGL3;SERPINA3;MCFD2;PTPRS;ELN;PTEN;FMN1;HIF3A;TFCP2L1;PTH1R;SYNE3;CTSS;PTPRG;RNF157;ADAMTS2;C1QTNF1;IMPA2;SH3PXD2B;FLVCR2;EFHD1;AOX1;CERS6;ZHX3;KLF13;ANXA2;IFNGR1;GPX3;NCOA3;SLC39A11;NGF;OSMR;SLC39A14;TGFBR2;TGFBR3;PSMA6;ARHGAP10;MMP14;TBC1D2;SLC7A7;SLC7A8;GFOD1;DPYD;PICK1;FAM20C;COL6A3;PLIN2;ITGA5;MOCS1;ERGIC1;TMEM45A;KANK1;C1S;ADCY3;TFPI;FSTL1;TMEM165;HDAC7;KIAA0513;MTHFD1L;CLMN;PTK2B;PTPN18;GALNT6;GSN;NEGR1;TPK1;CCDC57;TXNRD1;GSR;SUSD1;LHFPL2;MERTK;KLF9;IL18R1"

gluc_response <- "VIT;VKORC1L1;ERRFI1;AHCYL1;STEAP4;TRAF3IP2;GALNT15;SERPINE1;JADE1;SLA;CBLB;MT1X;EPS8;CCND3;BMPER;RASSF4;RPS6KA2;ANPEP;C1RL;MAP3K6;IL6R;PDGFRA;MLIP;SCARA5;IL1R1;EBF1;TTC7A;CRISPLD2;SPARCL1;FKBP5;NNMT;LPAR1;SLC1A3;PLA2G5;NID1;ACACB;ZFP36L2;PIK3R5;C3;SCFD2;LPXN;HACL1;SRGAP2;SLC38A2;SLC19A2;S100A10;KLHL29;GADD45B;ZBTB16;ELL2;CORO2B;IGF2R;NFATC4;DERA;SULT1B1;MAFB;BCL6;TMEM236;TBXAS1;NDUFAF2;RGL3;SERPINA3;MCFD2;PTPRS;ELN;PTEN;FMN1;HIF3A;TFCP2L1;PTH1R;SYNE3;CTSS;PTPRG;RNF157;ADAMTS2;C1QTNF1;IMPA2;SH3PXD2B;FLVCR2;EFHD1;AOX1;CERS6;ZHX3;KLF13;ANXA2;IFNGR1;GPX3;NCOA3;SLC39A11;NGF;OSMR;SLC39A14;TGFBR2;TGFBR3;PSMA6;ARHGAP10;MMP14;TBC1D2;SLC7A7;SLC7A8;GFOD1;DPYD;PICK1;FAM20C;COL6A3;PLIN2;ITGA5;MOCS1;ERGIC1;TMEM45A;KANK1;C1S;ADCY3;TFPI;FSTL1;TMEM165;HDAC7;KIAA0513;MTHFD1L;CLMN;PTK2B;PTPN18;GALNT6;GSN;NEGR1;TPK1;CCDC57;TXNRD1;GSR;SUSD1;LHFPL2;MERTK;KLF9;IL18R1"
gluc_response <- str_split(gluc_response[1],';')[[1]]

M1 <- AddModuleScore(M1,list(gluc_response),name='nr3c1')


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'myeloid_nr3c1.pdf'), width=3, height=3)

VlnPlot(M1,'nr3c11',group.by='group',pt.size=0)
dev.off()














