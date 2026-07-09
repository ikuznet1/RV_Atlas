###############################################################################
## Supplementary Figure 6 (v53 draft) — Murine PAB model: immune and vascular characterization
##
## Panels (from RV_snRNASeq_v52_draft.md figure legends):
##   (A) PAB echocardiographic characterization (sham, mild RVF, severe RVF)
##   (B) Cell population recovery mouse vs human RV
##   (C) EC expansion with RVF in PAB mice
##   (D) EC subtype subclustering in mouse
##   (E-G) Myeloid DEGs between conditions in mouse
##   (H) Pooled MHCII-associated expression decreased with mouse RVF (contrasting human)
##   (I) NR3C1 target expression unchanged in mouse RVF
##   (J) CX3CR1-tdTomato bulk RNA-seq of sorted PAB RV macrophages
##
## Source: copied from v51 Supplementary_Figure_5.R on 2026-04-10
## Status: SKELETON — v52 porting pending
##
## Output: ./output/Supplementary_Figure_5/v52_figures/SupplementaryFigure_6.pdf
###############################################################################

source('./helper_scripts/_shared_helpers.R')

## Per-figure output directory (introduced for consistent output paths)
V52_FIG_DIR <- './output/Supplementary_Figure_5'
dir.create(V52_FIG_DIR, showWarnings = FALSE, recursive = TRUE)


## Suppress R's default Rplots.pdf in cwd when Rscript hits a plot call
## that's outside an explicit pdf() ... dev.off() envelope.
pdf(NULL)
COMP_W <- 12
COMP_H <- 16

## Publication scaling constants (geom linewidths, point sizes, text mm, etc.)
PS <- pub_scales(COMP_W)

library(Seurat)
library(hdWGCNA)
library(ggeasy)
library(harmony)
library(readxl)
library(tidyr)
library(ggpubr)



source('./helper_scripts/spatial_functions.R')

###############################################################################
## Cartoon / external assets
##
## Supplementary Figure 5 Panel A contains the PAB mouse heart cartoon (top-left
## of the original PDF). To regenerate asset from the original PDF, run at a
## terminal (once -- output saved to ./new_scripts/assets/):
##
##   cd /Users/ikuz/Documents/RV_Atlas/new_scripts/assets
##   # 17 cm ≈ 2008 px at 300 DPI; source PDF rasterised at 300 DPI →
##   # the PAB cartoon occupies roughly the upper-left ~20 % × ~18 % of page.
##   # Adjust -crop WxH+Xoff+Yoff after visually confirming on the rendered png.
##   magick -density 300 ~/Downloads/hdWGCNA_TOM/Manuscripts/Supplementary_Figure_5.pdf \
##       -crop 400x360+0+0 +repage Supplementary_Figure_5_PAB_cartoon.png
##
p_S5_cartoon <- insert_asset('Supplementary_Figure_5_PAB_cartoon.png',
                             label = 'PAB model schematic')

#######################################
#############  FIGURE S5A  ############
# v53: PAB snRNA-seq UMAP across all recovered cell populations
# (sham, mild RVF, severe RVF). Echo characterization → Table S5.
#######################################

M1 <- readRDS('./dependencies/shared/PAB_data_clean.rds')

M1 <- SetIdent(M1, value = "Names")
M1$group <- M1$orig.ident


pdf(paste0('./output/Supplementary_Figure_5/', 'PAB_snUMAP.pdf'), width=5, height=5)
print(PlotEmbedding(M1,group.by='Names',point_size=1,plot_under=TRUE,plot_theme=umap_theme()+NoLegend(),raster_dpi=400,raster_scale=0.5))
dev.off()

#######################################
#############  FIGURE S5B  ############
# v53: marker gene expression defining the PAB cell-type annotation.
#######################################


#FeaturePlot(M1,reduction = "umap",'Ryr2',label=T)
#FeaturePlot(M1,reduction = "umap",'Dcn',label=T)
#FeaturePlot(M1,reduction = "umap",'Csf1r',label=T)
#FeaturePlot(M1,reduction = "umap",'Myh11',label=T)
#FeaturePlot(M1,reduction = "umap",'Pecam1',label=T) <- better EC
#FeaturePlot(M1,reduction = "umap",'Wt1',label=T)
#FeaturePlot(M1,reduction = "umap",'Pdgfrb',label=T)
#FeaturePlot(M1,reduction = "umap",'Prox1',label=T) <- LEC
#FeaturePlot(M1,reduction = "umap",'Cdh11',label=T) <- Endo


M1$Names <- factor(M1$Names, levels = c("CM","Atria","FB","Myeloid","SM","Epi","PC","LEC","Endo","EC","Neuron"))
feats <- c("Ryr2","Nppa","Dcn","Csf1r","Myh11","Wt1","Pdgfrb","Prox1","Cdh11","Pecam1","Slc35f1")

## v53 panel B: DotPlot of marker genes per cell type using the default Seurat
## blue gradient (lightgrey → blue). Dot size = % cells expressing,
## color = scaled mean expression. Vertical italic gene labels;
## legend at bottom in 2 horizontal rows.
pdf(paste0('./output/Supplementary_Figure_5/', 'PAB_sn_Dot.pdf'), width=5, height=5)
print(
  DotPlot(M1, features = feats, group.by = "Names",
          cols = c("lightgrey", "blue")) +
    xlab("") + ylab("") +
    theme_v52(COMP_W) +
    guides(
      colour = guide_colorbar(title = "Average Expression",
                              title.position = "top", order = 1),
      size   = guide_legend(   title = "Percent Expressed",
                               title.position = "top", order = 2, nrow = 1)
    ) +
    theme(
      panel.border    = element_rect(linewidth = PS$linewidth_mm, fill = NA, color = "black"),
      axis.line.x     = element_blank(),
      axis.line.y     = element_blank(),
      axis.text.x     = element_text(face = "italic", angle = 90, hjust = 1, vjust = 0.5),
      legend.position = "bottom",
      legend.box      = "horizontal",
      legend.key.size = PS$legend_key
    )
)
dev.off()

## Legacy violin alternative — kept for reversion.
# p <- VlnPlot(M1, features = feats, ncol = 1, pt.size = F, group.by = "Names")
# for(i in 1:11) {
#   p[[i]] <- p[[i]] + NoLegend() +
#     easy_remove_axes(which = "y", what = c("ticks","text","line")) +
#     ggtitle("") + ylab(feats[i])
#   if (i < 11) p[[i]] <- p[[i]] + easy_remove_axes(which = "x")
# }
# pdf(paste0('./output/Supplementary_Figure_5/', 'PAB_sn_Vln.pdf'), width = 4, height = 19)
# print(p)
# dev.off()


#######################################
#############  FIGURE S5C  ############
# v53: Cell population recovery per mouse (per-patient pseudobulk averaging,
# matching legacy ~/Downloads/hdWGCNA_TOM/Manuscripts/Supplementary_Figure_6.r).
#######################################


Sham_percent_cell <- t(table(subset(M1,group=="Nor")$patient,subset(M1,group=="Nor")$Names)) %*% diag(1/table(subset(M1,group=="Nor")$patient))
Sham_percent_cell <- data.frame(Freq = rowMeans(Sham_percent_cell), type = 'Sham', Var1 = rownames(Sham_percent_cell))
Sham_percent_cell$sum <- (rev(cumsum(rev(Sham_percent_cell$Freq))) - Sham_percent_cell$Freq/2)

Mild_percent_cell <- t(table(subset(M1,group=="Mod")$patient,subset(M1,group=="Mod")$Names)) %*% diag(1/table(subset(M1,group=="Mod")$patient))
Mild_percent_cell <- data.frame(Freq = rowMeans(Mild_percent_cell), type = 'Mild RVF', Var1 = rownames(Mild_percent_cell))
Mild_percent_cell$sum <- (rev(cumsum(rev(Mild_percent_cell$Freq))) - Mild_percent_cell$Freq/2)

Severe_percent_cell <- t(table(subset(M1,group=="Sev")$patient,subset(M1,group=="Sev")$Names)) %*% diag(1/table(subset(M1,group=="Sev")$patient))
Severe_percent_cell <- data.frame(Freq = rowMeans(Severe_percent_cell), type = 'Severe RVF', Var1 = rownames(Severe_percent_cell))
Severe_percent_cell$sum <- (rev(cumsum(rev(Severe_percent_cell$Freq))) - Severe_percent_cell$Freq/2)


percent_cell_df <- rbind(Sham_percent_cell, Mild_percent_cell, Severe_percent_cell)
percent_cell_df$Var1 <- factor(percent_cell_df$Var1, levels = rownames(Sham_percent_cell))
percent_cell_df$type <- factor(percent_cell_df$type, levels = c('Sham','Mild RVF','Severe RVF'))

percent_cell_df$label <- round(percent_cell_df$Freq, 2)
percent_cell_df$label[percent_cell_df$label < 0.03] <- NA
percent_cell_df$label <- scales::percent(percent_cell_df$label)

pdf('./output/Supplementary_Figure_5/PAB_prev_stacked.pdf',width=4.5,height=5)
print(
ggplot(percent_cell_df, aes(fill=Var1, y=Freq, x=type,label=round(sum,1))) +
geom_bar(position="stack", stat="identity", width=0.6, linewidth = PS$geom_lw) + theme_v52(COMP_W) +
xlab("Disease State") + ylab("Frequency") + labs(fill="Cell type",color='black') +
scale_y_touch() +
geom_label_repel(aes(type,sum,label=label),fill=NA,nudge_x=0.5,direction="y",
                 size = PS$text_mm, family = FONT_FAMILY) +
theme(legend.key.size = PS$legend_key)
)
dev.off()

#######################################
#############  FIGURE S5D-E  ##########
# v53: D = EC expansion with RVF in PAB; E = EC subtype subclustering
# (with human EC hdWGCNA module inset).
#######################################

M2 <- subset(M1, Names %in% c("EC","Endo","LEC"))

M2 <- RunPCA(M2)
M2 <- FindNeighbors(M2, dims = 1:50)
M2 <- FindClusters(M2, resolution = 0.5)
M2 <- RunUMAP(M2, dims = 1:10)


#7 LEC
#6 KIT+
#5 Art - hey1, jag1
#4 cap micro cd36
#2 vein - nr2f2
#3 endo - npr3
#4,8 - cap, rgcc
#0 gCaps - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10187412/ Aplnr, Cdh5, Myo10, and Ccdc85a
M2$Subnames <- M2@active.ident
labels <- c('Cap1','Cap2','Vein','Endo','Cap3','Art','Cap4','LEC','Cap5','Cap6')

names(labels) <- levels(M2)
M2 <- RenameIdents(M2, labels)
M2$Subnames <- M2@active.ident

M2$Subsubnames <- M2@active.ident

labels <- c('Cap','Cap','Vein','Endo','Cap','Art','Cap','LEC','Cap','Cap')

names(labels) <- levels(M2)
M2 <- RenameIdents(M2, labels)
M2$Subsubnames <- M2@active.ident

pdf(paste0('./output/Supplementary_Figure_5/', 'PAB_EC_snUMAP.pdf'), width=5, height=5)
print(PlotEmbedding(M2,group.by='Subsubnames',point_size=1,plot_under=TRUE,plot_theme=umap_theme()+NoLegend(),raster_dpi=400,raster_scale=0.5))
dev.off()


human2mouse <- read.csv('./dependencies/shared/human2mouse.csv',header=F)
idx <- match(unique(human2mouse[,2]),human2mouse[,2])
human2mouse<-human2mouse[idx,]
colnames(human2mouse) <-c('human_name', 'mouse_name')

M3 <- readRDS('./dependencies/shared/EC_hdWGCNA_by_celltype.rds')

consensus_modules <- GetModules(M3) %>% subset(module != 'grey')
consensus_modules <- consensus_modules[,1:3]

idx_match<- match(consensus_modules$gene_name,human2mouse$human_name)

consensus_modules$gene_name <- human2mouse$mouse_name[idx_match]
consensus_modules <- consensus_modules[!is.na(consensus_modules$gene_name),]
rownames(consensus_modules) <- consensus_modules$gene_name

consensus_modules <- subset(consensus_modules, gene_name %in% rownames(M2))
# remove duplicate gene names
consensus_modules <- consensus_modules[match(unique(consensus_modules$gene_name), consensus_modules$gene_name),]
library(dplyr)
score_calc <- consensus_modules %>% group_by(module) %>% group_split()
module_colors <- unique(unlist(lapply(score_calc,'[[','module')))

mapping <- labels2colors(1:100)

module_colors <- paste0('M',match(module_colors,mapping))
M2 <- AddModuleScore(M2,lapply(score_calc,'[[','gene_name'),name="module_score")
cols_current <- colnames(M2@meta.data)
cols_current[startsWith(colnames(M2@meta.data),'module_score')] <- paste0('module_',module_colors)
colnames(M2@meta.data) <- cols_current


M2<-SetIdent(M2,value='Subsubnames')


pdf(paste0('./output/Supplementary_Figure_5/', 'PAB_EC_trend_subcluster.pdf'), width=4.5, height=3)

p <- DotPlot(M2,paste0('module_',
  c('M1','M2','M3','M4','M5','M6','M7')),dot.min=0,col.min=0,col.max=2,
  scale.min=0, scale.max=100) +
  RotatedAxis() + ylab('')+ xlab('')+
  scale_color_gradient2(high='red', mid='grey95', low='blue') +
  scale_size(range = PS$dot_range) + theme_v52(COMP_W) +
  theme(
    panel.border = element_rect(linewidth = PS$linewidth_mm, fill=NA, color='black'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    legend.key.size = PS$legend_key
)
p
print(p)
dev.off()


M2 <- SetIdent(M2, value = "group")
Idents(M2) <- factor(x = Idents(M2), levels = c('Nor','Mod','Sev'))

pdf(paste0('./output/Supplementary_Figure_5/', 'PAB_trend_condition_EC.pdf'), width=5, height=2.5)

p <- DotPlot(M2,paste0('module_',
  c('M1','M2','M3','M4','M5','M6','M7')),dot.min=0,col.min=0,col.max=2,
  scale.min=0, scale.max=100) +
  RotatedAxis() + ylab('')+ xlab('')+
  scale_color_gradient2(high='red', mid='grey95', low='blue') +
  scale_size(range = PS$dot_range) + theme_v52(COMP_W) +
  theme(
    panel.border = element_rect(linewidth = PS$linewidth_mm, fill=NA, color='black'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    legend.key.size = PS$legend_key
)
p
print(p)
dev.off()


#######################################
#############  FIGURE S5F  ############
# v53: Per-mouse pseudobulk EC subtype proportions. Tabulates over Subsubnames
# (Cap/Vein/Endo/Art/LEC) per patient → averages across patients per group.
#######################################


Sham_percent_cell <- t(table(subset(M2,group=="Nor")$patient,subset(M2,group=="Nor")$Subsubnames)) %*% diag(1/table(subset(M2,group=="Nor")$patient))
Sham_percent_cell <- data.frame(Freq = rowMeans(Sham_percent_cell), type = 'Sham', Var1 = rownames(Sham_percent_cell))
Sham_percent_cell$sum <- (rev(cumsum(rev(Sham_percent_cell$Freq))) - Sham_percent_cell$Freq/2)

Mild_percent_cell <- t(table(subset(M2,group=="Mod")$patient,subset(M2,group=="Mod")$Subsubnames)) %*% diag(1/table(subset(M2,group=="Mod")$patient))
Mild_percent_cell <- data.frame(Freq = rowMeans(Mild_percent_cell), type = 'Mild RVF', Var1 = rownames(Mild_percent_cell))
Mild_percent_cell$sum <- (rev(cumsum(rev(Mild_percent_cell$Freq))) - Mild_percent_cell$Freq/2)

Severe_percent_cell <- t(table(subset(M2,group=="Sev")$patient,subset(M2,group=="Sev")$Subsubnames)) %*% diag(1/table(subset(M2,group=="Sev")$patient))
Severe_percent_cell <- data.frame(Freq = rowMeans(Severe_percent_cell), type = 'Severe RVF', Var1 = rownames(Severe_percent_cell))
Severe_percent_cell$sum <- (rev(cumsum(rev(Severe_percent_cell$Freq))) - Severe_percent_cell$Freq/2)

percent_cell_df <- rbind(Sham_percent_cell, Mild_percent_cell, Severe_percent_cell)
percent_cell_df$Var1 <- factor(percent_cell_df$Var1, levels = rownames(Sham_percent_cell))
percent_cell_df$type <- factor(percent_cell_df$type, levels = c('Sham','Mild RVF','Severe RVF'))

percent_cell_df$label <- round(percent_cell_df$Freq, 2)
percent_cell_df$label[percent_cell_df$label < 0.03] <- NA
percent_cell_df$label <- scales::percent(percent_cell_df$label)

pdf('./output/Supplementary_Figure_5/PAB_EC_prev_stacked.pdf',width=4.5,height=5)
print(
ggplot(percent_cell_df, aes(fill=Var1, y=Freq, x=type,label=round(sum,1))) +
geom_bar(position="stack", stat="identity", width=0.6, linewidth = PS$geom_lw) + theme_v52(COMP_W) +
xlab("Disease State") + ylab("Frequency") + labs(fill="Cell type",color='black') +
scale_y_touch() +
geom_label_repel(aes(type,sum,label=label),fill=NA,nudge_x=0.5,direction="y",
                 size = PS$text_mm, family = FONT_FAMILY) +
theme(legend.key.size = PS$legend_key)
)
dev.off()



#######################################
#############  FIGURE S5G  ############
# v53: Myeloid WGCNA module expression (M3, M8) in mouse APC and resident mac.
#######################################


M2 <- subset(M1, Names %in% c("Myeloid"))

M2 <- RunPCA(M2)
M2 <- FindNeighbors(M2, dims = 1:50)
M2 <- FindClusters(M2, resolution = 0.5)
M2 <- RunUMAP(M2, dims = 1:10)

#0 is CCR2- rMac
#1 looks like cardiac contam
#2 is
#3 is DC and CCR, HLA high
#4 B cells
#5 is something non immune
#6 EC

M2$Subnames <- M2@active.ident
M2 <- subset(M2, Subnames %in% c('0','1','2','3') )
M2 <- RunUMAP(M2, dims = 1:50)

#0 rMac
#1
#2 
#3 HLA
M2 <- subset(M2, Subnames %in% c('0','3') )
M2 <- RunPCA(M2)
M2 <- RunUMAP(M2, dims = 1:10)

labels <- c('rMac','HLA')

names(labels) <- levels(M2)
M2 <- RenameIdents(M2, labels)
M2$Subnames <- M2@active.ident


human2mouse <- read.csv('./dependencies/shared/human2mouse.csv',header=F)
idx <- match(unique(human2mouse[,2]),human2mouse[,2])
human2mouse<-human2mouse[idx,]
colnames(human2mouse) <-c('human_name', 'mouse_name')

consensus_modules <- read.csv("./dependencies/shared/bulk_heart_modules.csv")
consensus_modules <- consensus_modules[,1:3]

idx_match<- match(consensus_modules$gene_name,human2mouse$human_name)

consensus_modules$gene_name <- human2mouse$mouse_name[idx_match]
consensus_modules <- consensus_modules[!is.na(consensus_modules$gene_name),]
rownames(consensus_modules) <- consensus_modules$gene_name

consensus_modules <- subset(consensus_modules, gene_name %in% rownames(M2))
# remove duplicate gene names
consensus_modules <- consensus_modules[match(unique(consensus_modules$gene_name), consensus_modules$gene_name),]
library(dplyr)
score_calc <- consensus_modules %>% group_by(module) %>% group_split()
module_colors <- unique(unlist(lapply(score_calc,'[[','module')))

mapping <- labels2colors(1:100)

module_colors <- paste0('M',match(module_colors,mapping))
M2 <- AddModuleScore(M2,lapply(score_calc,'[[','gene_name'),name="module_score")
cols_current <- colnames(M2@meta.data)
cols_current[startsWith(colnames(M2@meta.data),'module_score')] <- paste0('module_',module_colors)
colnames(M2@meta.data) <- cols_current

M2 <- SetIdent(M2,value='Subnames')

modules_int <- c('M1','M3','M4',"M8")

pdf(paste0('./output/Supplementary_Figure_5/', 'PAB_myeloid_dot_subclust.pdf'), width=5, height=2)

p <- DotPlot(M2,paste0('module_',modules_int),group.by='Subnames',dot.min=0,col.min=-1,col.max=1,scale.min=50,scale.max=100) +
  RotatedAxis() + ylab('')+ xlab('')+
  scale_color_gradient2(high='red', mid='grey95', low='blue') +
  scale_size(range = PS$dot_range) + theme_v52(COMP_W) +
  theme(
    panel.border = element_rect(linewidth = PS$linewidth_mm, fill=NA, color='black'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    legend.key.size = PS$legend_key
)
p
print(p)
dev.off()

M2$group <- factor(M2$group,levels = c('Nor','Mod','Sev'))



pdf(paste0('./output/Supplementary_Figure_5/', 'PAB_myeloid_dot_disease.pdf'), width=5, height=2.5)

p <- DotPlot(M2,paste0('module_',modules_int),group.by='group',dot.min=0,col.min=-1,col.max=1,scale.min=50,scale.max=100) +
  RotatedAxis() + ylab('')+ xlab('')+
  scale_color_gradient2(high='red', mid='grey95', low='blue') +
  scale_size(range = PS$dot_range) + theme_v52(COMP_W) +
  theme(
    panel.border = element_rect(linewidth = PS$linewidth_mm, fill=NA, color='black'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    legend.key.size = PS$legend_key
)
p
print(p)
dev.off()

#######################################
#############  FIGURE S5H  ############
# v53: DEGs (volcano) between conditions in mouse myeloid cells.
#######################################

library(enrichR)


bulk_modules <- consensus_modules
bulk_modules$module <- match(bulk_modules$module,mapping)
dbs <- c('GO_Biological_Process_2023','GO_Cellular_Component_2023','GO_Molecular_Function_2023','ChEA_2022','Reactome_2022')

#Run enrichment by cell type
M2$Names_group <- paste0(M2$Subnames,'_',M2$group)
Idents(M2) <- "Names_group"
combined_set <- data.frame()
combined_output <- data.frame()

mods_idx <- c(1,3,4,8)
cell_types <- unique(M2$Subnames)
comparison <- list(c("Sev","Mod"),c("Mod","Nor"),c("Sev","Nor"))
for (i in mods_idx){
	for (j in cell_types){
		for (k in comparison){
			key_genes <- subset(bulk_modules,module %in% c(i))$gene_name
			key_genes <- key_genes[key_genes %in% rownames(M2)]

			gene_set <- FindMarkers(M2, ident.1 = paste0(j,"_",k[1]), ident.2 = paste0(j,"_",k[2]),features=key_genes,recorrect_umi=F)
			
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
			if (length(rownames(gene_enrich))==0){next}
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
			if (length(rownames(gene_enrich))==0){next}
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


selected_terms <- subset(combined_output,direction=="down")
selected_terms <- subset(selected_terms,comparison=="Sev_Nor")
selected_terms <- subset(selected_terms,Adjusted.P.value<0.05)
selected_terms



Idents(M2) <- "group"

library(EnhancedVolcano)

gene_set <- FindMarkers(M2, ident.1 = 'Sev', ident.2 = 'Nor',recorrect_umi=F,features=subset(bulk_modules,module==1)$gene_name)

pdf(paste0('./output/Supplementary_Figure_5/', 'PAB_M1_Myeloid_module_volcano_all_RVF_vs_NF.pdf'), width=8, height=6)

print(EnhancedVolcano(gene_set,lab=rownames(gene_set),
	x='avg_log2FC',y='p_val_adj',
	FCcutoff = 0.1,pCutoff=0.05, labFace = "italic") + coord_flip())
dev.off()

gene_set <- FindMarkers(M2, ident.1 = 'Sev', ident.2 = 'Nor',recorrect_umi=F,features=subset(bulk_modules,module==8)$gene_name)

pdf(paste0('./output/Supplementary_Figure_5/', 'PAB_M8_Myeloid_module_volcano_all_RVF_vs_NF.pdf'), width=8, height=6)

print(EnhancedVolcano(gene_set,lab=rownames(gene_set),
	x='avg_log2FC',y='p_val_adj',
	FCcutoff = 0.1,pCutoff=0.05, labFace = "italic") + coord_flip())
dev.off()

#######################################
############  FIGURE S5I/J  ###########
# v53: per-mouse pseudobulk (mean module score per patient) box+jitter with KW.
# Replaces single-cell VlnPlot to avoid pseudo-replication (n is the mouse).
#######################################


#MHCII

M2 <- AddModuleScore(M2,list(c('Ciita','Cd74','H2-Ab1','H2-Aa',
	'H2-Eb1','H2-Eb2','H2-Ob','H2-DMb1','H2-DMb2','H2-DMa','H2-Oa')),name='MHCII')

mhc_pb <- M2@meta.data %>%
  dplyr::filter(Subnames == 'HLA') %>%
  dplyr::group_by(patient, group) %>%
  dplyr::summarise(score = mean(MHCII1), .groups = 'drop') %>%
  dplyr::mutate(group = factor(dplyr::recode(group, Nor='Sham', Mod='Mild RVF', Sev='Severe RVF'),
                               levels = c('Sham','Mild RVF','Severe RVF')))
mhc_kw <- kruskal.test(score ~ group, data = mhc_pb)$p.value

pdf(paste0('./output/Supplementary_Figure_5/', 'PAB_myeloid_MHC.pdf'), width=3, height=3)
print(
  ggplot(mhc_pb, aes(x=group, y=score, fill=group)) +
    geom_boxplot(outlier.shape=NA, alpha=0.6) +
    geom_jitter(width=0.15, size=2) +
    theme_v52(COMP_W) + NoLegend() +
    ylab('MHCII module score') + xlab('') +
    ggtitle(sprintf('MHCII (HLA mac, KW p=%.2g)', mhc_kw)) +
    theme(plot.title = element_text(size=8),
          axis.text.x = element_text(angle=30, hjust=1, color='black'))
)
dev.off()



gluc_response <- "VIT;VKORC1L1;ERRFI1;AHCYL1;STEAP4;TRAF3IP2;GALNT15;SERPINE1;JADE1;SLA;CBLB;MT1X;EPS8;CCND3;BMPER;RASSF4;RPS6KA2;ANPEP;C1RL;MAP3K6;IL6R;PDGFRA;MLIP;SCARA5;IL1R1;EBF1;TTC7A;CRISPLD2;SPARCL1;FKBP5;NNMT;LPAR1;SLC1A3;PLA2G5;NID1;ACACB;ZFP36L2;PIK3R5;C3;SCFD2;LPXN;HACL1;SRGAP2;SLC38A2;SLC19A2;S100A10;KLHL29;GADD45B;ZBTB16;ELL2;CORO2B;IGF2R;NFATC4;DERA;SULT1B1;MAFB;BCL6;TMEM236;TBXAS1;NDUFAF2;RGL3;SERPINA3;MCFD2;PTPRS;ELN;PTEN;FMN1;HIF3A;TFCP2L1;PTH1R;SYNE3;CTSS;PTPRG;RNF157;ADAMTS2;C1QTNF1;IMPA2;SH3PXD2B;FLVCR2;EFHD1;AOX1;CERS6;ZHX3;KLF13;ANXA2;IFNGR1;GPX3;NCOA3;SLC39A11;NGF;OSMR;SLC39A14;TGFBR2;TGFBR3;PSMA6;ARHGAP10;MMP14;TBC1D2;SLC7A7;SLC7A8;GFOD1;DPYD;PICK1;FAM20C;COL6A3;PLIN2;ITGA5;MOCS1;ERGIC1;TMEM45A;KANK1;C1S;ADCY3;TFPI;FSTL1;TMEM165;HDAC7;KIAA0513;MTHFD1L;CLMN;PTK2B;PTPN18;GALNT6;GSN;NEGR1;TPK1;CCDC57;TXNRD1;GSR;SUSD1;LHFPL2;MERTK;KLF9;IL18R1"
gluc_response <- stringr::str_split(gluc_response[1],';')[[1]]


gluc_response_mouse <- human2mouse$mouse_name[match(gluc_response,human2mouse$human_name)]


M2 <- AddModuleScore(M2,list(gluc_response_mouse),name='nr3c1')


nr3c1_pb <- M2@meta.data %>%
  dplyr::group_by(patient, group) %>%
  dplyr::summarise(score = mean(nr3c11), .groups = 'drop') %>%
  dplyr::mutate(group = factor(dplyr::recode(group, Nor='Sham', Mod='Mild RVF', Sev='Severe RVF'),
                               levels = c('Sham','Mild RVF','Severe RVF')))
nr3c1_kw <- kruskal.test(score ~ group, data = nr3c1_pb)$p.value

pdf(paste0('./output/Supplementary_Figure_5/', 'PAB_myeloid_nr3c1.pdf'), width=3, height=3)
print(
  ggplot(nr3c1_pb, aes(x=group, y=score, fill=group)) +
    geom_boxplot(outlier.shape=NA, alpha=0.6) +
    geom_jitter(width=0.15, size=2) +
    theme_v52(COMP_W) + NoLegend() +
    ylab('NR3C1 target score') + xlab('') +
    ggtitle(sprintf('NR3C1 targets (all myeloid, KW p=%.2g)', nr3c1_kw)) +
    theme(plot.title = element_text(size=8),
          axis.text.x = element_text(angle=30, hjust=1, color='black'))
)
dev.off()

#######################################
#############  FIGURE S5I  ############
# v53: CX3CR1-tdTomato lineage tracing — bulk RNA-seq of sorted RV macs.
# (Legend I = bulk-RNAseq myeloid; banner relabel from K -> I per v54 legend.)
#######################################

### Plot J i

M1 <- readRDS('./dependencies/shared/PAB_data_clean.rds')

M1 <- SetIdent(M1, value = "Names")
M1$group <- M1$orig.ident


M2 <- subset(M1, Names %in% c("Myeloid"))

M2 <- RunPCA(M2)
M2 <- FindNeighbors(M2, dims = 1:50)
M2 <- FindClusters(M2, resolution = 0.5)
M2 <- RunUMAP(M2, dims = 1:10)

#0 is CCR2- rMac
#1 looks like cardiac contam
#2 is
#3 is DC and CCR, HLA high
#4 B cells
#5 is something non immune
#6 EC

M2$Subnames <- M2@active.ident
M2 <- subset(M2, Subnames %in% c('0','1','2','3') )
M2 <- RunUMAP(M2, dims = 1:50)

#0 rMac
#1
#2 
#3 HLA
M2 <- subset(M2, Subnames %in% c('0','3') )
M2 <- RunPCA(M2)
M2 <- RunUMAP(M2, dims = 1:10)

labels <- c('rMac','HLA')

names(labels) <- levels(M2)
M2 <- RenameIdents(M2, labels)
M2$Subnames <- M2@active.ident


library(tximport)

##### LOAD BULK DATA ##### 
dir <- './dependencies/shared/Mouse_PAB_myeloid'

samples <- read.csv(file.path(dir, "meta.csv"), header = TRUE)
samples$Batch <- factor(samples$Batch)
samples$Pressure.Loading <- factor(samples$Pressure.Loading)


files <- file.path(dir,'nascent', samples$ID, "abundance.h5")
names(files) <- samples$ID
tx2gene <- read.table('./dependencies/shared/Mouse_PAB_myeloid/t2g.txt',fill=T)
tx2gene = data.frame(TXNAME=tx2gene$V1,GENEID=tx2gene$V3)

txi.kallisto <- tximport(files, type = "kallisto", tx2gene = tx2gene)

rownames(samples) <- samples$ID

#### EMBED INTO ALL PAB #####
# bulk <- CreateSeuratObject(counts = txi.kallisto$counts, meta.data = data.frame(samples))

# bulk <- SCTransform(bulk, vst.flavor = "v2")
# bulk <- RunPCA(bulk, npcs = 10, verbose = FALSE)


# M1 <- SplitObject(M1, split.by = "patient")
# M1<-PrepSCTIntegration(M1)
# features<-SelectIntegrationFeatures(M1)
# M1.anchors<-FindIntegrationAnchors(M1,normalization.method = 'SCT',anchor.features = features, reference = c(1,2,3), reduction = "rpca")
# M1 <- IntegrateData(anchorset = M1.anchors,normalization.method='SCT')

# DefaultAssay(M1) <- "integrated"

# M1 <- RunPCA(M1, npcs = 50, verbose = FALSE)
# M1 <- RunUMAP(M1, reduction = "pca", dims = 1:30)


# anchors <- FindTransferAnchors(
#   reference = M1,
#   query = bulk,
#   normalization.method = "SCT",
#   recompute.residuals=FALSE,
#   reference.reduction = "pca",
#   dims = 1:50,
#   k.score =15
# )


# predictions <- TransferData(anchorset = anchors, refdata = M1$Names, dims = 1:50,k.weight=15)

# M1 <- RunUMAP(M1, dims = 1:50, return.model = TRUE)

# bulk <- MapQuery(anchorset = anchors, reference = M1, query = bulk,
# 	refdata = list(celltype = "Names"), reference.reduction = "pca", 
# 	reduction.model = "umap",transferdata.args = list(k.weight=15))

# #score <- MappingScore(anchors)

# #bulk$map_score <- score

# p1 <- DimPlot(M1, reduction = "umap", group.by = "Names", label = TRUE, label.size = 3, repel = TRUE,raster=TRUE,pt.size=1.5) + NoLegend() + ggtitle("Reference annotations")
# p2 <- DimPlot(bulk, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE, label.size = 3, pt.size=1.5,repel = TRUE,raster=TRUE) + NoLegend() + ggtitle("Query transferred labels")
# pdf(paste0('./output/Supplementary_Figure_5/', 'RV_PAB_bulkMyeloid_ref_mapped.pdf'), width=10, height=5)
# p1 + p2
# dev.off()



tryCatch({  # v53: skip ref-mapping if multi-batch SCT error fires
  #### LOAD BULK DATA AND EMBED INTO MYELOID PAB
  
  bulk <- CreateSeuratObject(counts = txi.kallisto$counts, meta.data = data.frame(samples))
  
  bulk <- SCTransform(bulk, vst.flavor = "v2")
  bulk <- RunPCA(bulk, npcs = 10, verbose = FALSE)
  
  
  #M2 <- SplitObject(M2, split.by = "patient")
  #M2<-PrepSCTIntegration(M2)
  #features<-SelectIntegrationFeatures(M2)
  #M2.anchors<-FindIntegrationAnchors(M2,normalization.method = 'SCT',anchor.features = features, reduction = "rpca")
  #M2 <- IntegrateData(anchorset = M2.anchors,normalization.method='SCT')
  
  #DefaultAssay(M2) <- "integrated"
  
  #M1 <- RunPCA(M2, npcs = 50, verbose = FALSE)
  #M1 <- RunUMAP(M2, reduction = "pca", dims = 1:30)
  M3 <- subset(M1, Names %in% c("Myeloid"))
  
  true_myeloid <- colnames(M2)
  M3 <- subset(M3,cells=true_myeloid)
  M3$Subnames = M2$Subnames
  
  anchors <- FindTransferAnchors(
    reference = M3,
    query = bulk,
    normalization.method = "SCT",
    recompute.residuals=FALSE,
    reference.reduction = "pca",
    dims = 1:50,
    k.score =15
  )
  
  
  predictions <- TransferData(anchorset = anchors, refdata = M3$Subnames, dims = 1:50,k.weight=10)
  
  M3 <- RunUMAP(M3, dims = 1:50, return.model = TRUE)
  
  bulk <- MapQuery(anchorset = anchors, reference = M3, query = bulk,
  	refdata = list(celltype = "Subnames"), reference.reduction = "pca", 
  	reduction.model = "umap",transferdata.args = list(k.weight=10))
  
  #score <- MappingScore(anchors)
  
  #bulk$map_score <- score
  
  
  bulk$group <- paste0(bulk$Origin,'_',bulk$Type)
  p1 <- DimPlot(M3, reduction = "umap", group.by = "Subnames", label = TRUE, label.size = 3, repel = TRUE,raster=TRUE,pt.size=1.5) + NoLegend() + ggtitle("Reference annotations")
  p2 <- DimPlot(bulk, reduction = "ref.umap", group.by = "group", label = TRUE, label.size = 3, pt.size=1.5,repel = TRUE,raster=TRUE) + ggtitle("Query transferred labels")
  pdf(paste0('./output/Supplementary_Figure_5/', 'RV_PAB_Myeloid_Only_bulkMyeloid_ref_mapped.pdf'), width=10, height=5)
  print(p1 + p2)
  dev.off()
  
  bulk_RV = subset(bulk, Origin=='RV')
  
  bulk_RV<-SetIdent(bulk_RV,value='Type')
  
  M3<-SetIdent(M3,value='group')
  
  
  DefaultAssay(M3) <- 'SCT'
  
  a <- FindMarkers(bulk_RV,ident.1='PAB',ident.2='Sham',min.pct=0.1,logfc.threshold=0)
  
  M3 <- PrepSCTFindMarkers(M3)
  b <-  FindMarkers(M3,ident.1='Sev',ident.2='Nor',min.pct=0.25,recorrect_umi=F,logfc.threshold=0)
  shared <- intersect(rownames(a),rownames(b))
  dataset <- data.frame(bulk=a[shared,]$avg_log2FC,sn=b[shared,]$avg_log2FC)
  rownames(dataset) <- shared
  labs <- rownames(dataset)
  
  cor(dataset[,1],dataset[,2])
  
  pdf(paste0('./output/Supplementary_Figure_5/', 'RV_PAB_Myeloid_Only_bulkMyeloid_scatter.pdf'), width=10, height=10)
  print(ggplot(dataset, aes(x = sn, y=bulk)) + geom_point(size = PS$scatter_pt) +
    geom_text_repel(label=labs,max.overlaps=25,
                    size = PS$text_mm, family = FONT_FAMILY, fontface = "italic") + theme_v52(COMP_W) +
    theme(legend.key.size = PS$legend_key))
  dev.off()
  
  cat(rownames(subset(dataset, bulk > 0.1 & sn > 0.1)),sep='\n')
  cat(rownames(subset(dataset, bulk < 0.1 & sn < 0.1)),sep='\n')

  ## v57 S5 Panel I sub-panels (I-i = bulk-vs-sn scatter above):
  ##  I-ii  = PCA of bulk RV macrophages — PAB vs Sham cluster separation
  ##  I-iii = bulk RV PAB-vs-Sham DEG volcano (MHCII transcripts down)
  tryCatch({
    .v  <- bulk_RV[['pca']]@stdev^2
    .v  <- round(100 * .v / sum(.v))
    .pd <- data.frame(PC1 = Embeddings(bulk_RV, 'pca')[, 1],
                      PC2 = Embeddings(bulk_RV, 'pca')[, 2],
                      grp = paste0('RV ', as.character(bulk_RV$Type)))
    p_S5Iii <- ggplot(.pd, aes(PC1, PC2, color = grp)) +
      geom_point(size = 3) +
      scale_color_manual(values = c('RV PAB' = '#E41A1C',
                                    'RV Sham' = '#4DAF4A'), name = NULL) +
      labs(x = sprintf('PC1: %d%% variance', .v[1]),
           y = sprintf('PC2: %d%% variance', .v[2]),
           title = 'Module score') +
      theme_v52(COMP_W)
    save_figure(p_S5Iii, 'Figure_S5_panel_Iii.pdf', width = 4.5, height = 3.5)
    ## I-iii is the script's own bulk_PAB_myeloid_RV_PAB_vs_RV_Sham.pdf
    ## (generated later, L~989) — copied to Figure_S5_panel_Iiii.pdf by
    ## the rename-map at the end. Do NOT build a custom volcano here.
    message('S5 panel I-ii (PCA) built')
  }, error = function(e)
     message('S5 I-ii build failed: ', conditionMessage(e)))
  
  
  
  ###Compare to human RV
  
  
  bulk_RV <- AddModuleScore(bulk_RV,list(c('Ciita','Cd74','H2-Ab1','H2-Aa',
  	'H2-Eb1','H2-Eb2','H2-Ob','H2-DMb1','H2-DMb2','H2-DMa','H2-Oa')),name='MHCII')
  
  VlnPlot(bulk_RV,'MHCII1')
  
  human2mouse <- read.csv('./dependencies/shared/human2mouse.csv',header=F)
  idx <- match(unique(human2mouse[,2]),human2mouse[,2])
  human2mouse<-human2mouse[idx,]
  colnames(human2mouse) <-c('human_name', 'mouse_name')
  gluc_response <- "VIT;VKORC1L1;ERRFI1;AHCYL1;STEAP4;TRAF3IP2;GALNT15;SERPINE1;JADE1;SLA;CBLB;MT1X;EPS8;CCND3;BMPER;RASSF4;RPS6KA2;ANPEP;C1RL;MAP3K6;IL6R;PDGFRA;MLIP;SCARA5;IL1R1;EBF1;TTC7A;CRISPLD2;SPARCL1;FKBP5;NNMT;LPAR1;SLC1A3;PLA2G5;NID1;ACACB;ZFP36L2;PIK3R5;C3;SCFD2;LPXN;HACL1;SRGAP2;SLC38A2;SLC19A2;S100A10;KLHL29;GADD45B;ZBTB16;ELL2;CORO2B;IGF2R;NFATC4;DERA;SULT1B1;MAFB;BCL6;TMEM236;TBXAS1;NDUFAF2;RGL3;SERPINA3;MCFD2;PTPRS;ELN;PTEN;FMN1;HIF3A;TFCP2L1;PTH1R;SYNE3;CTSS;PTPRG;RNF157;ADAMTS2;C1QTNF1;IMPA2;SH3PXD2B;FLVCR2;EFHD1;AOX1;CERS6;ZHX3;KLF13;ANXA2;IFNGR1;GPX3;NCOA3;SLC39A11;NGF;OSMR;SLC39A14;TGFBR2;TGFBR3;PSMA6;ARHGAP10;MMP14;TBC1D2;SLC7A7;SLC7A8;GFOD1;DPYD;PICK1;FAM20C;COL6A3;PLIN2;ITGA5;MOCS1;ERGIC1;TMEM45A;KANK1;C1S;ADCY3;TFPI;FSTL1;TMEM165;HDAC7;KIAA0513;MTHFD1L;CLMN;PTK2B;PTPN18;GALNT6;GSN;NEGR1;TPK1;CCDC57;TXNRD1;GSR;SUSD1;LHFPL2;MERTK;KLF9;IL18R1"
  gluc_response <- stringr::str_split(gluc_response[1],';')[[1]]
  gluc_response_mouse <- human2mouse$mouse_name[match(gluc_response,human2mouse$human_name)]
  
  bulk_RV <- AddModuleScore(bulk_RV,list(gluc_response_mouse),name='nr3c1')
  
  VlnPlot(bulk_RV,'nr3c11')
  
}, error = function(e) message("[S5K skip ref-mapping] ", conditionMessage(e)))

### Plot J ii and iii



library('tximportData')
library('tximport')
library(DESeq2)
library(EnhancedVolcano)
library("sva")



dir <- './dependencies/shared/Mouse_PAB_myeloid'

samples <- read.csv(file.path(dir, "meta.csv"), header = TRUE)
samples$Batch <- factor(samples$Batch)
samples$Pressure.Loading <- factor(samples$Pressure.Loading)


files <- file.path(dir,'nascent', samples$ID, "abundance.h5")
names(files) <- samples$ID
tx2gene <- read.table('./dependencies/shared/Mouse_PAB_myeloid/t2g.txt',fill=T)
tx2gene = data.frame(TXNAME=tx2gene$V1,GENEID=tx2gene$V3)

txi.kallisto <- tximport(files, type = "kallisto", tx2gene = tx2gene)


#RV VS LV

dds <- DESeqDataSetFromTximport(txi.kallisto, samples, ~Origin)


smallestGroupSize <- 8
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

dds$Origin <- relevel(dds$Origin, ref = "LV")

dds <- DESeq(dds)
res <- results(dds)

resOrdered <- res[order(res$pvalue),]
resLFC <- lfcShrink(dds, coef="Origin_RV_vs_LV", type="apeglm")

pdf('./output/Supplementary_Figure_5/bulk_PAB_myeloid_RV_vs_LV.pdf',width=10,height=10)
print(EnhancedVolcano(resLFC,
    lab = rownames(res),
    x = 'log2FoldChange',
    y = 'padj',
    pCutoff = 0.05,
    ylim=c(0,4),
    labFace = "italic"))
dev.off()

lv <- cat(rownames(subset(resOrdered,padj<0.05 & log2FoldChange<0)),sep="\n")
rv <- cat(rownames(subset(resOrdered,padj<0.05 & log2FoldChange>0)),sep="\n")


#PAB VS SHAM

dds <- DESeqDataSetFromTximport(txi.kallisto, samples, ~Type)


smallestGroupSize <- 8
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

dds$Type <- relevel(dds$Type, ref = "Sham")

dds <- DESeq(dds)
res <- results(dds)

resOrdered <- res[order(res$pvalue),]
resLFC <- lfcShrink(dds, coef="Type_PAB_vs_Sham", type="apeglm")

pdf('./output/Supplementary_Figure_5/bulk_PAB_myeloid_PAB_vs_Sham.pdf',width=10,height=10)
print(EnhancedVolcano(resLFC,
    lab = rownames(res),
    x = 'log2FoldChange',
    y = 'padj',
    pCutoff = 0.05,
    ylim=c(0,4),
    labFace = "italic"))
dev.off()

#ALL COMPARISONS

samples$group <- paste0(samples$Origin,'_',samples$Type)
dds <- DESeqDataSetFromTximport(txi.kallisto, samples, ~group+Batch)


smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

dds$group <- relevel(dds$group, ref = "RV_Sham")


dds <- DESeq(dds)
res <- results(dds,contrast=c("group","RV_PAB","RV_Sham"))

resOrdered <- res[order(res$pvalue),]
resLFC <- lfcShrink(dds, coef="group_RV_PAB_vs_RV_Sham", type="apeglm")

pdf('./output/Supplementary_Figure_5/bulk_PAB_myeloid_RV_PAB_vs_RV_Sham.pdf',width=10,height=10)
print(EnhancedVolcano(resLFC,
    lab = rownames(res),
    x = 'log2FoldChange',
    y = 'padj',
    pCutoff = 0.05,
    ylim=c(0,4),
    labFace = "italic"))
dev.off()

rv_pab_up <- cat(rownames(subset(resOrdered,padj<0.05 & log2FoldChange>0.25 & baseMean > 100)),sep="\n")
rv_pab_down <- cat(rownames(subset(resOrdered,padj<0.05 & log2FoldChange<0)),sep="\n")


cat(intersect(rownames(M1),rownames(subset(resOrdered,padj<0.05 & log2FoldChange>0.25))),sep="\n")


dds <- DESeqDataSetFromTximport(txi.kallisto, samples, ~group)


smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

dds$group <- relevel(dds$group, ref = "LV_Sham")

dds <- DESeq(dds)


res <- results(dds,contrast=c("group","LV_PAB","LV_Sham"))


resOrdered <- res[order(res$pvalue),]
resLFC <- lfcShrink(dds, coef="group_LV_PAB_vs_LV_Sham", type="apeglm")

pdf('./output/Supplementary_Figure_5/bulk_PAB_myeloid_LV_PAB_vs_LV_Sham.pdf',width=10,height=10)
print(EnhancedVolcano(res,
    lab = rownames(res),
    x = 'log2FoldChange',
    y = 'padj',
    pCutoff = 0.05,
    ylim=c(0,4),
    labFace = "italic"))
dev.off()

lv_pab_up <- cat(rownames(subset(resOrdered,padj<0.05 & log2FoldChange>0)),sep="\n")


# plotCounts(dds, gene='Wnt3a', intgroup="group")  # diagnostic; gene may be filtered out — disabled


#Regress out variability

library(sva)
library(ggplot2)
library(magrittr)
samples$group <- paste0(samples$Origin,'_',samples$Type)


dds <- DESeqDataSetFromTximport(txi.kallisto, samples, ~group)
dds <- estimateSizeFactors(dds)
smallestGroupSize <- 5
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

normalized_counts <- counts(dds, normalized=TRUE)


mod  <- model.matrix(~Pressure.Loading+Surgical.BW+Echo.BW+Cells+Age.at.Sac+Batch+Rejected.Cells, colData(dds))
mod0 <- model.matrix(~Surgical.BW+Echo.BW+Cells+Age.at.Sac+Batch+Rejected.Cells, colData(dds))
svseq <- svaseq(normalized_counts, mod, mod0)

samples$SV1 <- svseq$sv[,1]
samples$SV2 <- svseq$sv[,2]
samples$SV3 <- svseq$sv[,3]


dds <- DESeqDataSetFromTximport(txi.kallisto, samples, ~group+SV1+SV2+SV3)

dds <- estimateSizeFactors(dds)

smallestGroupSize <- 5
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

normalized_counts <- counts(dds, normalized=TRUE)

#dds$Pressure.Loading <- relevel(dds$Pressure.Loading, ref = "0")


dds <- DESeq(dds)

vstSE <- vst(dds,blind = FALSE)

mat <- assay(vstSE)
mm <- model.matrix(~group, colData(vstSE))
mat <- limma::removeBatchEffect(mat, covariates=colData(vstSE)[,24:26], design=mm)
assay(vstSE) <- mat


plotPCA(vstSE,intgroup="group",ntop=1000)

rv_pab <- lfcShrink(dds,contrast=c('group','RV_PAB','RV_Sham'), type="ashr")
lv_pab <- lfcShrink(dds,contrast=c('group','LV_PAB','LV_Sham'), type="ashr")

#rv_pab <- lfcShrink(dds,contrast=c('Pressure.Loading','1','0'), type="ashr")


EnhancedVolcano(rv_pab,
    lab = rownames(rv_pab),
    x = 'log2FoldChange',
    y = 'padj',
    pCutoff = 0.05,
    labFace = "italic")

## Diagnostic plotCounts — skip silently if a gene is missing from filtered dds.
for (.g in c('Tbxa2r','Gdf6','Ugt1a8','Slc35d3','1700123L14Rik')) {
  if (.g %in% rownames(dds)) {
    try(plotCounts(dds, gene = .g, intgroup = "group"), silent = TRUE)
  }
}

EnhancedVolcano(rv_pab,
    lab = rownames(rv_pab),
    x = 'log2FoldChange',
    y = 'padj',
    pCutoff = 0.05,
    xlim=c(-4,4),ylim=c(0,5),
    labFace = "italic")

cat(rownames(subset(rv_pab,padj<0.05 & log2FoldChange>0)),sep="\n")

cat(rownames(subset(rv_pab,padj<0.05 & log2FoldChange<0)),sep="\n")























# M2 <- subset(M1, Names %in% c("CM"))

# M2 <- RunPCA(M2)
# M2 <- RunHarmony(M2,'patient')
# M2 <- FindNeighbors(M2, dims = 1:50,reduction = "harmony")
# M2 <- FindClusters(M2, resolution = 0.5,reduction = "harmony")
# M2 <- RunUMAP(M2, dims = 1:50,reduction = "harmony")

# M2$Names <- M2@active.ident
# markers<-FindAllMarkers(M2,recorrect_umi=F)

# pdf(paste0('./output/Supplementary_Figure_5/', 'PAB_CM_snUMAP.pdf'), width=5, height=5)
# PlotEmbedding(M2,group.by='Names',point_size=1,plot_under=TRUE,plot_theme=umap_theme()+NoLegend(),raster_dpi=400,raster_scale=0.5)
# dev.off()



# M3 <- readRDS(file = "./output/Supplementary_Figure_5/cm_new_subclust.rds")

# #new.cluster.ids <- c("Cm1","Cm2","Cm3","Cm4","Cm5","Cm6","Cm7","Cm8","Cm9","Cm10")
# #names(new.cluster.ids) <- levels(M3)
# #M3 <- RenameIdents(M3, new.cluster.ids)

# #M3$Subnames <- M3@active.ident
# #M3$SubNames_Groups <- paste(M1$Subnames,M3$group,sep='_')


# #M3 <- AddModuleScore(M3, features=list(c('MALAT1')),assay="SCT",name="Clust0Score")
# #M3 <- AddModuleScore(M3, features=list(c('FGF12','SH3RF2','KCNMB2','PRELID2')),assay="SCT",name="Clust1Score")
# #M3 <- AddModuleScore(M3, features=list(c('TNNT2','TTN','MYBPC3','MYH7')),assay="SCT",name="Clust2Score")
# #M3 <- AddModuleScore(M3, features=list(c('PALLD','MYO18B','MYPN','ANKRD1')),assay="SCT",name="Clust5Score")
# #M3 <- AddModuleScore(M3, features=list(c('PDE3A','CDH2','PDLIM5')),assay="SCT",name="Clust4Score")
# #M3 <- AddModuleScore(M3, features=list(c('AKAP13','OBSCN','LARGE1','THSD4')),assay="SCT",name="Clust3Score")
# #M3 <- AddModuleScore(M3, features=list(c('PALLD','SORBS2','CAMK2D','CCSER1','PDLIM5')),assay="SCT",name="Clust7Score")
# #M3 <- AddModuleScore(M3, features=list(c('AC020637.1','LINC02388')),assay="SCT",name="Clust6Score")
# #M3 <- AddModuleScore(M3, features=list(c('MIR646HG')),assay="SCT",name="Clust8Score")
# #M3 <- AddModuleScore(M3, features=list(c('GPC5','HS6ST3')),assay="SCT",name="Clust9Score")

# #DefaultAssay(M3) <- "RNA"

# #M3[["RNA"]] <- split(M3[["RNA"]], f = M3$patient)
# #M3[['SCT']] <- NULL
# #M3[['decontXcounts']] <- NULL
# #M3 <- SCTransform(M3, vst.flavor = "v2")
# #M3 <- RunPCA(M3, npcs = 50, verbose = FALSE)
# #M3 <- SplitObject(M3, split.by = "patient")
# #M3<-PrepSCTIntegration(M3)
# #features<-SelectIntegrationFeatures(M3)
# #M3.anchors<-FindIntegrationAnchors(M3,normalization.method = 'SCT',anchor.features = features, reduction = "rpca")
# #M3 <- IntegrateData(anchorset = M3.anchors,normalization.method='SCT')

# #DefaultAssay(M3) <- "integrated"

# #M3 <- RunPCA(M3, npcs = 50, verbose = FALSE)
# #M3 <- RunUMAP(M3, reduction = "pca", dims = 1:30)

# # RenameGenesSeurat <- function(obj = ls.Seurat[[i]],
# #                               newnames = HGNC.updated[[i]]$Suggested.Symbol,
# #                               assay = "RNA",
# #                               slots = c("data", "counts", "meta.features")) {
# #   #
# #   message("RenameGenesSeurat, assay: ", assay)
# #   warning("Run this before integration and downstream processing. It only attempts to change
# #           @counts, @data, and @meta.features in obj@assays$YOUR_ASSAY.", immediate. = TRUE)

# #   stopifnot(
# #     "Unequal gene name sets: nrow(assayobj) != nrow(newnames):" =
# #       length(Features(obj, assay = assay)) == length(newnames)
# #   )

# #   if (obj@version < 5) warning("obj@version < 5. Old versions are not supported. Update the obj!", immediate. = TRUE)

# #   if ("scale.data" %in% slots) {
# #     n_genes_sc_dta <- nrow(obj@assays[[assay]]$"scale.data")
# #     stopifnot(
# #       "scale.data does has different number of genes than newnames!" =
# #         n_genes_sc_dta == length(newnames)
# #     )
# #   }

# #   LayersFound <- SeuratObject::Layers(obj@assays[[assay]])
# #   #print("Present: ", sort(LayersFound))

# #   slots <- sort(intersect(slots, LayersFound))
# #   #print("Replaced: ", slots)

# #     for (slotX in slots) {
# #     print(slotX)
# #     if (slotX == "scale.data") browser()
# #     nrO <- nrow(SeuratObject::GetAssayData(object = obj, assay = assay, layer = slotX))
# #     obj <- .check_and_rename(obj, assay, newnames = newnames, layer.name = slotX)
# #     nrN <- nrow(SeuratObject::GetAssayData(object = obj, assay = assay, layer = slotX))
# #     stopifnot(nrN == nrO)
# #   }
# #   return(obj)
# # }


# # .check_and_rename <- function(obj, assay, newnames, layer.name) {
# #   cat(layer.name, fill = TRUE)

# #   length_newnames <- length(newnames)
# #   length_orig_names <- length(Features(obj, assay = assay))

# #   stopifnot(
# #     is(obj, "Seurat"),
# #     is.character(assay),
# #     is.character(layer.name),
# #     is.character(newnames),
# #     length_orig_names == length_newnames
# #   )

# #   assayobj <- obj@assays[[assay]]
# #   feature.list <- rownames(assayobj@features@.Data)

# #   if (length(feature.list) == length(newnames)) {
# #     rownames(assayobj@features@.Data) <- newnames
# #     nrX <- length(rownames(assayobj@features@.Data))
# #   } else {
# #     paste0("length feature.list", length(feature.list), "length newnames", length(newnames))
# #     stop()
# #   }

# #   if (layer.name %in% SeuratObject::Layers(assayobj)) {
# #     matrix_n <- SeuratObject::LayerData(assayobj, layer = layer.name)
# #     nr1 <- nrow(matrix_n)

# #     if (all(dim(matrix_n)) > 0) {
# #       # browser()
# #       stopifnot(nrow(matrix_n) == length(newnames))

# #       if ("dgCMatrix" %in% class(matrix_n)) {
# #         message(assay, "@", layer.name, " is of type dgeCMatrix!")
# #         matrix_n@Dimnames[[1]] <- newnames
# #       } else if ("matrix" %in% class(matrix_n)) {
# #         message(assay, "@", layer.name, " is of type Matrix!")
# #         rownames(matrix_n) <- newnames
# #       } else if ("data.frame" %in% class(matrix_n)) {
# #         message(assay, "@", layer.name, " is of type data.frame!")
# #         rownames(matrix_n) <- newnames
# #       } else {
# #         warning(">>> No renaming: ", assay, "@", layer.name,
# #           " not of type dgeCMatrix / Matrix / data.frame.",
# #           immediate. = TRUE
# #         )
# #       }
# #       stopifnot(nr1 == nrow(matrix_n))

# #       SeuratObject::LayerData(assayobj, layer = layer.name) <- matrix_n
# #       nr3 <- nrow(SeuratObject::LayerData(assayobj, layer = layer.name))
# #       stopifnot(nr3 == nrX)
# #     }
# #   } else {
# #     warning(paste(">>>", assay, "@", layer.name, "does not exist!"), immediate. = TRUE)
# #   }
# #   # obj <- SetAssayData(obj, layer = layer.name, new.data = matrix_n)
# #   obj@assays[[assay]] <- assayobj
# #   return(obj)
# # }

# # DefaultAssay(M2) <- "RNA"


# # newnames <- human2mouse$human_name[match(rownames(M2),human2mouse$mouse_name)]
# # newnames[is.na(newnames)] <- rownames(M2)[is.na(newnames)]

# # RenameGenesSeurat(M2,newnames = newnames,assay='RNA')



# RNA <- M2@assays$RNA['counts']
# newnames <- human2mouse$human_name[match(rownames(RNA),human2mouse$mouse_name)]
# newnames[is.na(newnames)] <- rownames(RNA)[is.na(newnames)]
# rownames(RNA) <- newnames


# M4 <- CreateAssayObject(RNA)
# M2[['humanized']] <- M4
# DefaultAssay(M2) <- "humanized"
# M2[['SCT']] <- NULL

# M2 <- SCTransform(M2, vst.flavor = "v2",assay='humanized')
# M2 <- RunPCA(M2, npcs = 50, verbose = FALSE)

# anchors <- FindTransferAnchors(
#   reference = M3,
#   query = M2,
#   normalization.method = "SCT",
#   recompute.residuals=FALSE,
#   reference.reduction = "harmony",
#   dims = 1:30
# )

# predictions <- TransferData(anchorset = anchors, refdata = M3$Subnames, dims = 1:30)


# M3 <- RunUMAP(M3, dims = 1:30, reduction = 'harmony', return.model = TRUE)
# M2 <- MapQuery(anchorset = anchors, reference = M3, query = M2,
# 	refdata = list(celltype = "Subnames"), reference.reduction = "pca", reduction.model = "umap")

# score <- MappingScore(anchors,ndim = 30)

# M2$map_score <- score

# p1 <- DimPlot(M3, reduction = "umap", group.by = "Subnames", label = TRUE, label.size = 3, repel = TRUE,raster=TRUE,pt.size=1.5) + NoLegend() + ggtitle("Reference annotations")
# p2 <- DimPlot(M2, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE, label.size = 3, pt.size=1.5,repel = TRUE,raster=TRUE) + NoLegend() + ggtitle("Query transferred labels")
# pdf(paste0('./output/Supplementary_Figure_5/', 'RV_PAB_CM_ref_mapped.pdf'), width=10, height=5)
# p1 + p2
# dev.off()


# pdf(paste0('./output/Supplementary_Figure_5/', 'PAB_CM_ref_mapped.pdf'), width=5, height=5)
# PlotEmbedding(M2,group.by='predicted.celltype',reduction = "ref.umap",point_size=0.2,plot_under=TRUE,plot_theme=umap_theme()+NoLegend(),raster_dpi=400,raster_scale=0.5)
# dev.off()


# pdf(paste0('./output/Supplementary_Figure_5/', 'RV_CM_ref_mapped.pdf'), width=5, height=5)
# PlotEmbedding(M3,group.by='Subnames',point_size=0.2,plot_under=TRUE,plot_theme=umap_theme()+NoLegend(),raster_dpi=400,raster_scale=0.5)
# dev.off()


# FeaturePlot(M2,'map_score',reduction = "ref.umap")

.s5_dir <- './output/Supplementary_Figure_5/'
dir.create(.s5_dir, showWarnings = FALSE, recursive = TRUE)

###############################################################################
## FIGURE S5G — cross-species myeloid programs in PAB (gene sets v60)
## Per-mouse pseudobulk module scores for the 4 corrected Figure-6 / Data-S17
## myeloid programs (GR targets, MHC-II, Inflammasome, Type-I IFN), mouse PAB
## myeloid, box + per-mouse jitter by Nor/Mod/Sev, KW + pairwise Wilcoxon.
###############################################################################
tryCatch({
  suppressPackageStartupMessages({library(Seurat); library(ggplot2)
                                  library(dplyr); library(patchwork)})
  ## Corrected Figure-6 / Data-S17 myeloid programs (v60). MHC-II given as MOUSE
  ## orthologs (human HLA-* do not toupper-map to mouse H2-*); GR / Inflammasome /
  ## Type-I-IFN map fine from human symbols via the .mm() matcher below.
  GR_targets   <- c('FKBP5','TSC22D3','ZBTB16','KLF9','ANGPTL4')
  MHCII        <- c('Ciita','Cd74','H2-Ab1','H2-Aa','H2-Eb1','H2-Eb2','H2-Ob',
                    'H2-DMb1','H2-DMb2','H2-DMa','H2-Oa')
  Inflammasome <- c('NLRP3','NLRP1','AIM2','CASP1','PYCARD','IL1B','IL18',
                    'GSDMD','CASP4','NAIP')
  TypeI_IFN    <- c('ISG15','MX1','MX2','OAS1','OAS2','OAS3','IFIT1','IFIT3',
                    'IRF7','RSAD2','USP18','IFI6','IFI44')
  .pab <- readRDS('./dependencies/shared/PAB_data_clean.rds')
  .pab <- subset(.pab, Names == 'Myeloid')
  DefaultAssay(.pab) <- 'RNA'
  if (!'data' %in% SeuratObject::Layers(.pab[['RNA']]))
    .pab <- NormalizeData(.pab, verbose = FALSE)
  ## human symbol -> match against PAB rownames (try as-is, UPPER, Title)
  .rn <- rownames(.pab); .rnU <- toupper(.rn)
  .mm <- function(g) {
    h <- intersect(g, .rn)
    h <- union(h, .rn[match(setdiff(toupper(g), toupper(h)), .rnU)])
    h[!is.na(h)]
  }
  ## --- POOLED myeloid: all four programs scored over ALL myeloid cells ---
  ## MHC-II is computed over pooled myeloid (not the HLA subset): per-mouse the
  ## HLA-high subcluster is too sparse (Sham contributes a single mouse), so the
  ## subset boxplot is not robust. Pooled myeloid keeps all four panels uniform.
  .poolprogs <- list(GR_targets = GR_targets, MHCII = MHCII,
                     Inflammasome = Inflammasome, TypeI_IFN = TypeI_IFN)
  for (nm in names(.poolprogs)) {
    gg <- .mm(.poolprogs[[nm]])
    message('  S5G ', nm, ' (pooled myeloid): ', length(gg), '/', length(.poolprogs[[nm]]), ' genes mapped')
    .pab <- AddModuleScore(.pab, features = list(gg), name = nm, assay = 'RNA')
  }
  ## per-mouse pseudobulk: orig.ident -> Sham/Mild/Severe; patient = mouse
  .pseudob <- function(obj, scorecol) {
    m <- obj@meta.data
    m$.pt <- if ('patient' %in% colnames(m)) as.character(m$patient) else as.character(m$orig.ident)
    m$.gr <- factor(dplyr::recode(as.character(m$orig.ident), Nor='Sham', Mod='Mild RVF', Sev='Severe RVF'),
                    levels = c('Sham','Mild RVF','Severe RVF'))
    m <- m[!is.na(m$.gr) & !is.na(m$.pt), ]
    m %>% group_by(.pt, .gr) %>% summarise(score = mean(.data[[scorecol]]), .groups = 'drop')
  }
  ## panel definitions (title, per-mouse pseudobulk); all over pooled myeloid
  .panels <- list(list(t = 'GR targets',   pb = .pseudob(.pab, 'GR_targets1')),
                  list(t = 'MHC-II',       pb = .pseudob(.pab, 'MHCII1')),
                  list(t = 'Inflammasome', pb = .pseudob(.pab, 'Inflammasome1')),
                  list(t = 'Type-I IFN',   pb = .pseudob(.pab, 'TypeI_IFN1')))
  .plist <- lapply(seq_along(.panels), function(i) {
    d <- data.frame(score = .panels[[i]]$pb$score, grp = .panels[[i]]$pb$.gr)
    .nf  <- d$score[d$grp == 'Sham']                       # NF = Sham
    .dis <- d$score[d$grp %in% c('Mild RVF','Severe RVF')]  # disease = mild or severe RVF
    pp <- tryCatch(suppressWarnings(wilcox.test(.nf, .dis)$p.value), error = function(e) NA)
    ggplot(d, aes(grp, score, fill = grp)) +
      geom_boxplot(outlier.shape = NA, width = 0.78) +
      geom_jitter(width = 0.12, size = 3.6, color = 'black') +
      scale_fill_manual(values = c(Sham='#F8766D',`Mild RVF`='#00BA38',`Severe RVF`='#619CFF')) +
      coord_flip() +                                  # 90-deg rotation (match S5H)
      ggtitle(sprintf('%s (p=%.2g)', .panels[[i]]$t, pp)) +   # NF vs disease Wilcoxon; comparison defined in caption
      ylab(if (i == length(.panels)) 'Per-mouse mean score' else NULL) + xlab(NULL) +
      theme_classic(base_size = 24) +                 # all fonts uniform (24 pt)
      theme(legend.position = 'none',
            plot.title = element_text(size = 24, face = 'bold'),
            axis.text  = element_text(size = 24, color = 'black'),
            axis.title = element_text(size = 24))
  })
  p_S5G <- wrap_plots(.plist, ncol = 1)                # vertical stack (match S5H)
  save_figure(p_S5G, 'Figure_S5_panel_G.pdf', width = 6.4, height = 11)
  rm(.pab); invisible(gc(verbose = FALSE))
  message('S5 panel G (cross-species myeloid programs) built')
}, error = function(e) message('S5 G build failed: ', conditionMessage(e)))

###############################################################################
## FIGURE S5H — PAB flow-cytometry validation (RE-GATED, v60)
## Re-gated from raw .fcs (Box:Jonathan) with CORRECTED gating: residency = CCR2
## alone (CCR2- resident = quadrants Q1+Q4); MHCII is a READOUT (gMFI), not a gate
## -- the legacy 'resident = CCR2-MHCII+' gate was circular for the MHCII question.
## Upstream GatingSet build + per-sample extraction live in helper_scripts/pab_flow/
## (steps 01-10; run ONCE on the raw .fcs / external drive -> emits the cached
## per-sample table read below; the 6.3 GB GatingSets are NOT a dependency).
## n=22 RV (Set1 RV10 excluded: a Sham with 76% CCR2+ / 94% MHCII- = debris/dead-
## cell artifact, no viability dye in panel). Per-experiment Sham-MEDIAN
## normalization (display); set-adjusted log-scale trend test (MFI/fold-changes
## are multiplicative) shown in titles.
###############################################################################
tryCatch({
  suppressPackageStartupMessages({library(ggplot2); library(patchwork)})
  .pf <- read.csv('./dependencies/shared/pab_flow_6pop_RV.csv', check.names = FALSE)
  .pf$group <- factor(.pf$group, levels = c('Sham','Mild RVF','Severe RVF'))
  .pf$ord   <- as.integer(.pf$group)
  .mets <- c(Leukocytes = 'Pan-immune (CD45+)', Myeloid = 'Myeloid (CD11b+)',
             Macrophages = 'Macrophage (CD64+)', Resident = 'Resident mac (CCR2-)',
             TREM2 = 'TREM2+ mac', MHCII_gMFI_res = 'MHCII gMFI (resident)')
  ## set-adjusted log-scale trend test (p in titles)
  .ltp <- sapply(names(.mets), function(m)
    summary(lm(log(.pf[[m]]) ~ .pf$ord + .pf$set))$coef['.pf$ord', 'Pr(>|t|)'])
  ## per-experiment Sham-MEDIAN normalization (display only; robust at small n)
  for (m in names(.mets)) for (s in unique(.pf$set)) {
    .sh <- median(.pf[[m]][.pf$set == s & .pf$group == 'Sham'], na.rm = TRUE)
    .pf[[m]][.pf$set == s] <- .pf[[m]][.pf$set == s] / .sh }
  .cols <- c(Sham = '#F8766D', `Mild RVF` = '#00BA38', `Severe RVF` = '#619CFF')
  .mk <- function(m, showx = FALSE)
    ggplot(.pf, aes(group, .data[[m]], fill = group)) +
      geom_hline(yintercept = 1, linetype = 'dashed', color = 'grey40', linewidth = .3) +
      geom_boxplot(outlier.shape = NA, width = 0.78) +
      geom_jitter(width = 0.12, size = 1.8, color = 'black') +
      scale_fill_manual(values = .cols) +
      scale_x_discrete(expand = expansion(add = 0.55)) + coord_flip() +
      ggtitle(sprintf('%s (trend p=%.2g)', .mets[m], .ltp[m])) +
      ylab(if (showx) 'Fold-change vs Sham median' else NULL) + xlab(NULL) +
      theme_classic(base_size = 20) +
      theme(legend.position = 'none', plot.title = element_text(size = 17, face = 'bold'),
            axis.text = element_text(color = 'black'))
  .ml <- names(.mets)
  p_S5H <- wrap_plots(lapply(seq_along(.ml), function(i) .mk(.ml[i], i == length(.ml))), ncol = 1)
  ## replicate the pab_flow folder outputs alongside the S5 panel
  dir.create('./output/pab_flow', showWarnings = FALSE, recursive = TRUE)
  save_figure(p_S5H, 'Figure_S5_panel_H.pdf', width = 5.99, height = 12.15)
  ggsave('./output/pab_flow/PAB_flow_S5H_exact_RV.pdf', p_S5H, width = 5.99, height = 12.15)
  ggsave('./output/pab_flow/PAB_flow_S5H_exact_RV.png', p_S5H, width = 5.99, height = 12.15, dpi = 150)
  message('S5 panel H (PAB flow, re-gated) built -- MHCII gMFI trend p=',
          signif(.ltp['MHCII_gMFI_res'], 2))
}, error = function(e) message('S5 H (re-gated) build failed: ', conditionMessage(e)))

###############################################################################
## v57 standardized panel names — copy legacy PDFs to Figure_S5_panel_<L>.
## G & H are BUILT above (saved directly) — excluded here so they are not
## overwritten. S5 header is SKELETON; absent legacy PDFs flagged.
###############################################################################
## Mapping reconciled to v57 legend (.figure_run_logs/v57_figure_legends.md):
## A=PAB UMAP, B=marker dot, C=celltype-prop bar, D=EC UMAP+inset,
## E=EC-prop bar, F=mouse myeloid WGCNA module dot, I(i/ii/iii)=
## bulk-sorted-mac orthogonal validation (I-ii PCA = task TODO).
.s5_map <- c(
  A   = 'PAB_snUMAP.pdf',
  B   = 'PAB_sn_Dot.pdf',
  C   = 'PAB_prev_stacked.pdf',
  D   = 'PAB_EC_snUMAP.pdf',
  E   = 'PAB_EC_prev_stacked.pdf',
  F   = 'PAB_myeloid_dot_subclust.pdf',
  Ii   = 'RV_PAB_Myeloid_Only_bulkMyeloid_scatter.pdf',
  Iiii = 'bulk_PAB_myeloid_RV_PAB_vs_RV_Sham.pdf')
## I-ii (PCA) is BUILT directly above as Figure_S5_panel_Iii.pdf.
## I-iii is the script's own bulk_PAB_myeloid_RV_PAB_vs_RV_Sham.pdf
## (generated ~L989) -> copied to Figure_S5_panel_Iiii.pdf here.
for (.L in names(.s5_map)) {
  .src <- file.path(.s5_dir, .s5_map[[.L]])
  if (file.exists(.src)) {
    file.copy(.src, file.path(.s5_dir, sprintf('Figure_S5_panel_%s.pdf', .L)),
              overwrite = TRUE)
  } else message('S5 panel ', .L, ': legacy PDF missing (',
                 .s5_map[[.L]], ') — likely a SKELETON porting TODO')
}

#######################################
