library(Seurat)
library(hdWGCNA)
library(ggeasy)
library(harmony)



source('~/Downloads/hdWGCNA_TOM/spatial_functions.R')

#######################################
#############  FIGURE S6A  ############
#######################################

M1 <- readRDS('~/Downloads/hdWGCNA_TOM/PAB_data_clean.rds')

M1 <- SetIdent(M1, value = "Names")
M1$group <- M1$orig.ident


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'PAB_snUMAP.pdf'), width=5, height=5)
PlotEmbedding(M1,group.by='Names',point_size=1,plot_under=TRUE,plot_theme=umap_theme()+NoLegend(),raster_dpi=400,raster_scale=0.5)
dev.off()

#######################################
#############  FIGURE S6B  ############
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
p<-VlnPlot(M1, features = feats,ncol=1,pt.size=F,group.by="Names")



for(i in 1:11) {  
	 p[[i]] <- p[[i]] + NoLegend() + easy_remove_axes(which="y",what = c("ticks", "text","line")) + ggtitle("") + ylab(feats[i])
	 if(i<11){p[[i]]<-p[[i]]+easy_remove_axes(which="x")}
	 
}


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'PAB_sn_Vln.pdf'), width=4, height=19)
p
dev.off()


#######################################
#############  FIGURE S6C  ############
#######################################



NF_percent_cell <- cbind(as.data.frame(table(subset(M1,group=="Nor")@active.ident)/length(subset(M1,group=="Nor")@active.ident)*100),type = "NF")
NF_percent_cell$sum <- (rev(cumsum(rev(NF_percent_cell$Freq))) - NF_percent_cell$Freq/2)/100
NF_percent_cell$Freq <- NF_percent_cell$Freq/100


pRV_percent_cell <- cbind(as.data.frame(table(subset(M1,group=="Mod")@active.ident)/length(subset(M1,group=="Mod")@active.ident)*100),type = "pRV")
pRV_percent_cell$sum <- (rev(cumsum(rev(pRV_percent_cell$Freq))) - pRV_percent_cell$Freq/2)/100
pRV_percent_cell$Freq <- pRV_percent_cell$Freq/100


RVF_percent_cell <- cbind(as.data.frame(table(subset(M1,group=="Sev")@active.ident)/length(subset(M1,group=="Sev")@active.ident)*100),type = "RVF")
RVF_percent_cell$sum <- (rev(cumsum(rev(RVF_percent_cell$Freq))) - RVF_percent_cell$Freq/2)/100
RVF_percent_cell$Freq <- RVF_percent_cell$Freq/100

percent_cell_df <- rbind(NF_percent_cell,pRV_percent_cell,RVF_percent_cell)

percent_cell_df$label <- round(percent_cell_df$Freq,2)
percent_cell_df$label[percent_cell_df$label<0.03] = NA
percent_cell_df$label <- scales::percent(percent_cell_df$label)

pdf('~/Downloads/hdWGCNA_TOM/PAB_prev_stacked.pdf',width=4.5,height=5)
ggplot(percent_cell_df, aes(fill=Var1, y=Freq, x=type,label=round(sum,1))) +  
geom_bar(position="stack", stat="identity",width=0.6) + theme_classic() + 
xlab("Disease State") + ylab("Frequency") + labs(fill="Cell type",color='black') + 
theme(text = element_text(size=20),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),legend.text=element_text(color="black")) + 
scale_y_continuous(expand=c(0,0)) + 
geom_label_repel(aes(type,sum,label=label),fill=NA,nudge_x=0.5,direction="y")
dev.off()

#######################################
#############  FIGURE S6D  ############
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

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'PAB_EC_snUMAP.pdf'), width=5, height=5)
PlotEmbedding(M2,group.by='Subsubnames',point_size=1,plot_under=TRUE,plot_theme=umap_theme()+NoLegend(),raster_dpi=400,raster_scale=0.5)
dev.off()

#######################################
#############  FIGURE S6E  ############
#######################################


NF_percent_cell <- cbind(as.data.frame(table(subset(M2,group=="Nor")@active.ident)/length(subset(M2,group=="Nor")@active.ident)*100),type = "NF")
NF_percent_cell$sum <- (rev(cumsum(rev(NF_percent_cell$Freq))) - NF_percent_cell$Freq/2)/100
NF_percent_cell$Freq <- NF_percent_cell$Freq/100


pRV_percent_cell <- cbind(as.data.frame(table(subset(M2,group=="Mod")@active.ident)/length(subset(M2,group=="Mod")@active.ident)*100),type = "pRV")
pRV_percent_cell$sum <- (rev(cumsum(rev(pRV_percent_cell$Freq))) - pRV_percent_cell$Freq/2)/100
pRV_percent_cell$Freq <- pRV_percent_cell$Freq/100


RVF_percent_cell <- cbind(as.data.frame(table(subset(M2,group=="Sev")@active.ident)/length(subset(M2,group=="Sev")@active.ident)*100),type = "RVF")
RVF_percent_cell$sum <- (rev(cumsum(rev(RVF_percent_cell$Freq))) - RVF_percent_cell$Freq/2)/100
RVF_percent_cell$Freq <- RVF_percent_cell$Freq/100

percent_cell_df <- rbind(NF_percent_cell,pRV_percent_cell,RVF_percent_cell)

percent_cell_df$label <- round(percent_cell_df$Freq,2)
percent_cell_df$label[percent_cell_df$label<0.03] = NA
percent_cell_df$label <- scales::percent(percent_cell_df$label)

pdf('~/Downloads/hdWGCNA_TOM/PAB_EC_prev_stacked.pdf',width=4.5,height=5)
ggplot(percent_cell_df, aes(fill=Var1, y=Freq, x=type,label=round(sum,1))) +  
geom_bar(position="stack", stat="identity",width=0.6) + theme_classic() + 
xlab("Disease State") + ylab("Frequency") + labs(fill="Cell type",color='black') + 
theme(text = element_text(size=20),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),legend.text=element_text(color="black")) + 
scale_y_continuous(expand=c(0,0)) + 
geom_label_repel(aes(type,sum,label=label),fill=NA,nudge_x=0.5,direction="y")
dev.off()

#######################################
#############  FIGURE S6F  ############
#######################################
human2mouse <- read.csv('~/Downloads/hdWGCNA_TOM/human2mouse.csv',header=F)
idx <- match(unique(human2mouse[,2]),human2mouse[,2])
human2mouse<-human2mouse[idx,]
colnames(human2mouse) <-c('human_name', 'mouse_name')

M3 <- readRDS('~/Downloads/hdWGCNA_TOM/EC_hdWGCNA_by_celltype.rds')

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


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'PAB_EC_trend_subcluster.pdf'), width=4.5, height=3)

p <- DotPlot(M2,paste0('module_',
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


M2 <- SetIdent(M2, value = "group")
Idents(M2) <- factor(x = Idents(M2), levels = c('Nor','Mod','Sev'))

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'PAB_trend_condition_EC.pdf'), width=5, height=2.5)

p <- DotPlot(M2,paste0('module_',
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
#############  FIGURE S6G  ############
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


human2mouse <- read.csv('~/Downloads/hdWGCNA_TOM/human2mouse.csv',header=F)
idx <- match(unique(human2mouse[,2]),human2mouse[,2])
human2mouse<-human2mouse[idx,]
colnames(human2mouse) <-c('human_name', 'mouse_name')

consensus_modules <- read.csv("~/Downloads/hdWGCNA_TOM/bulk_heart_modules.csv")
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

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'PAB_myeloid_dot_subclust.pdf'), width=5, height=2)

p <- DotPlot(M2,paste0('module_',modules_int),group.by='Subnames',dot.min=0,col.min=-1,col.max=1,scale.min=50,scale.max=100) +
  RotatedAxis() + ylab('')+ xlab('')+
  scale_color_gradient2(high='red', mid='grey95', low='blue') +
  theme(
    panel.border = element_rect(size=1,fill=NA, color='black'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()
) 
p
dev.off()

M2$group <- factor(M2$group,levels = c('Nor','Mod','Sev'))



pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'PAB_myeloid_dot_disease.pdf'), width=5, height=2.5)

p <- DotPlot(M2,paste0('module_',modules_int),group.by='group',dot.min=0,col.min=-1,col.max=1,scale.min=50,scale.max=100) +
  RotatedAxis() + ylab('')+ xlab('')+
  scale_color_gradient2(high='red', mid='grey95', low='blue') +
  theme(
    panel.border = element_rect(size=1,fill=NA, color='black'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()
) 
p
dev.off()


#MHCII

M2 <- AddModuleScore(M2,list(c('Ciita','Cd74','H2-Ab1','H2-Aa',
	'H2-Eb1','H2-Eb2','H2-Ob','H2-DMb1','H2-DMb2','H2-DMa','H2-Oa')),name='MHCII')

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'PAB_myeloid_MHC.pdf'), width=3, height=3)

VlnPlot(subset(M2,Subnames=='HLA'),'MHCII1',group.by='group',pt.size=0)
dev.off()



gluc_response <- "VIT;VKORC1L1;ERRFI1;AHCYL1;STEAP4;TRAF3IP2;GALNT15;SERPINE1;JADE1;SLA;CBLB;MT1X;EPS8;CCND3;BMPER;RASSF4;RPS6KA2;ANPEP;C1RL;MAP3K6;IL6R;PDGFRA;MLIP;SCARA5;IL1R1;EBF1;TTC7A;CRISPLD2;SPARCL1;FKBP5;NNMT;LPAR1;SLC1A3;PLA2G5;NID1;ACACB;ZFP36L2;PIK3R5;C3;SCFD2;LPXN;HACL1;SRGAP2;SLC38A2;SLC19A2;S100A10;KLHL29;GADD45B;ZBTB16;ELL2;CORO2B;IGF2R;NFATC4;DERA;SULT1B1;MAFB;BCL6;TMEM236;TBXAS1;NDUFAF2;RGL3;SERPINA3;MCFD2;PTPRS;ELN;PTEN;FMN1;HIF3A;TFCP2L1;PTH1R;SYNE3;CTSS;PTPRG;RNF157;ADAMTS2;C1QTNF1;IMPA2;SH3PXD2B;FLVCR2;EFHD1;AOX1;CERS6;ZHX3;KLF13;ANXA2;IFNGR1;GPX3;NCOA3;SLC39A11;NGF;OSMR;SLC39A14;TGFBR2;TGFBR3;PSMA6;ARHGAP10;MMP14;TBC1D2;SLC7A7;SLC7A8;GFOD1;DPYD;PICK1;FAM20C;COL6A3;PLIN2;ITGA5;MOCS1;ERGIC1;TMEM45A;KANK1;C1S;ADCY3;TFPI;FSTL1;TMEM165;HDAC7;KIAA0513;MTHFD1L;CLMN;PTK2B;PTPN18;GALNT6;GSN;NEGR1;TPK1;CCDC57;TXNRD1;GSR;SUSD1;LHFPL2;MERTK;KLF9;IL18R1"
gluc_response <- str_split(gluc_response[1],';')[[1]]


gluc_response_mouse <- human2mouse$mouse_name[match(gluc_response,human2mouse$human_name)]


M2 <- AddModuleScore(M2,list(gluc_response_mouse),name='nr3c1')


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'PAB_myeloid_nr3c1.pdf'), width=3, height=3)

VlnPlot(M2,'nr3c11',group.by='group',pt.size=0)
dev.off()



#######################################
#############  FIGURE S6H  ############
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


#selected_terms <- subset(combined_output,db=="ChEA_2022")
selected_terms <- subset(combined_output,direction=="down")
selected_terms <- subset(selected_terms,comparison=="Sev_Nor")
selected_terms <- subset(selected_terms,Adjusted.P.value<0.05)
#selected_terms <- subset(selected_terms,module %in% c('M1'))
selected_terms



Idents(M2) <- "group"

library(EnhancedVolcano)

gene_set <- FindMarkers(M2, ident.1 = 'Sev', ident.2 = 'Nor',recorrect_umi=F,features=subset(bulk_modules,module==1)$gene_name)

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'PAB_M1_Myeloid_module_volcano_all_RVF_vs_NF.pdf'), width=8, height=6)

EnhancedVolcano(gene_set,lab=rownames(gene_set),
	x='avg_log2FC',y='p_val_adj',
	FCcutoff = 0.1,pCutoff=0.05) + coord_flip()
dev.off()

gene_set <- FindMarkers(M2, ident.1 = 'Sev', ident.2 = 'Nor',recorrect_umi=F,features=subset(bulk_modules,module==8)$gene_name)

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'PAB_M8_Myeloid_module_volcano_all_RVF_vs_NF.pdf'), width=8, height=6)

EnhancedVolcano(gene_set,lab=rownames(gene_set),
	x='avg_log2FC',y='p_val_adj',
	FCcutoff = 0.1,pCutoff=0.05) + coord_flip()
dev.off()

#######################################
#############  FIGURE S6I  ############
#######################################

M2 <- subset(M1, Names %in% c("CM"))

M2 <- RunPCA(M2)
M2 <- RunHarmony(M2,'patient')
M2 <- FindNeighbors(M2, dims = 1:50,reduction = "harmony")
M2 <- FindClusters(M2, resolution = 0.5,reduction = "harmony")
M2 <- RunUMAP(M2, dims = 1:50,reduction = "harmony")

M2$Names <- M2@active.ident
markers<-FindAllMarkers(M2,recorrect_umi=F)

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'PAB_CM_snUMAP.pdf'), width=5, height=5)
PlotEmbedding(M2,group.by='Names',point_size=1,plot_under=TRUE,plot_theme=umap_theme()+NoLegend(),raster_dpi=400,raster_scale=0.5)
dev.off()



M3 <- readRDS(file = "~/Downloads/hdWGCNA_TOM/cm_new_subclust.rds")

#new.cluster.ids <- c("Cm1","Cm2","Cm3","Cm4","Cm5","Cm6","Cm7","Cm8","Cm9","Cm10")
#names(new.cluster.ids) <- levels(M3)
#M3 <- RenameIdents(M3, new.cluster.ids)

#M3$Subnames <- M3@active.ident
#M3$SubNames_Groups <- paste(M1$Subnames,M3$group,sep='_')


#M3 <- AddModuleScore(M3, features=list(c('MALAT1')),assay="SCT",name="Clust0Score")
#M3 <- AddModuleScore(M3, features=list(c('FGF12','SH3RF2','KCNMB2','PRELID2')),assay="SCT",name="Clust1Score")
#M3 <- AddModuleScore(M3, features=list(c('TNNT2','TTN','MYBPC3','MYH7')),assay="SCT",name="Clust2Score")
#M3 <- AddModuleScore(M3, features=list(c('PALLD','MYO18B','MYPN','ANKRD1')),assay="SCT",name="Clust5Score")
#M3 <- AddModuleScore(M3, features=list(c('PDE3A','CDH2','PDLIM5')),assay="SCT",name="Clust4Score")
#M3 <- AddModuleScore(M3, features=list(c('AKAP13','OBSCN','LARGE1','THSD4')),assay="SCT",name="Clust3Score")
#M3 <- AddModuleScore(M3, features=list(c('PALLD','SORBS2','CAMK2D','CCSER1','PDLIM5')),assay="SCT",name="Clust7Score")
#M3 <- AddModuleScore(M3, features=list(c('AC020637.1','LINC02388')),assay="SCT",name="Clust6Score")
#M3 <- AddModuleScore(M3, features=list(c('MIR646HG')),assay="SCT",name="Clust8Score")
#M3 <- AddModuleScore(M3, features=list(c('GPC5','HS6ST3')),assay="SCT",name="Clust9Score")

#DefaultAssay(M3) <- "RNA"

#M3[["RNA"]] <- split(M3[["RNA"]], f = M3$patient)
#M3[['SCT']] <- NULL
#M3[['decontXcounts']] <- NULL
#M3 <- SCTransform(M3, vst.flavor = "v2")
#M3 <- RunPCA(M3, npcs = 50, verbose = FALSE)
#M3 <- SplitObject(M3, split.by = "patient")
#M3<-PrepSCTIntegration(M3)
#features<-SelectIntegrationFeatures(M3)
#M3.anchors<-FindIntegrationAnchors(M3,normalization.method = 'SCT',anchor.features = features, reduction = "rpca")
#M3 <- IntegrateData(anchorset = M3.anchors,normalization.method='SCT')

#DefaultAssay(M3) <- "integrated"

#M3 <- RunPCA(M3, npcs = 50, verbose = FALSE)
#M3 <- RunUMAP(M3, reduction = "pca", dims = 1:30)

# RenameGenesSeurat <- function(obj = ls.Seurat[[i]],
#                               newnames = HGNC.updated[[i]]$Suggested.Symbol,
#                               assay = "RNA",
#                               slots = c("data", "counts", "meta.features")) {
#   #
#   message("RenameGenesSeurat, assay: ", assay)
#   warning("Run this before integration and downstream processing. It only attempts to change
#           @counts, @data, and @meta.features in obj@assays$YOUR_ASSAY.", immediate. = TRUE)

#   stopifnot(
#     "Unequal gene name sets: nrow(assayobj) != nrow(newnames):" =
#       length(Features(obj, assay = assay)) == length(newnames)
#   )

#   if (obj@version < 5) warning("obj@version < 5. Old versions are not supported. Update the obj!", immediate. = TRUE)

#   if ("scale.data" %in% slots) {
#     n_genes_sc_dta <- nrow(obj@assays[[assay]]$"scale.data")
#     stopifnot(
#       "scale.data does has different number of genes than newnames!" =
#         n_genes_sc_dta == length(newnames)
#     )
#   }

#   LayersFound <- SeuratObject::Layers(obj@assays[[assay]])
#   #print("Present: ", sort(LayersFound))

#   slots <- sort(intersect(slots, LayersFound))
#   #print("Replaced: ", slots)

#     for (slotX in slots) {
#     print(slotX)
#     if (slotX == "scale.data") browser()
#     nrO <- nrow(SeuratObject::GetAssayData(object = obj, assay = assay, layer = slotX))
#     obj <- .check_and_rename(obj, assay, newnames = newnames, layer.name = slotX)
#     nrN <- nrow(SeuratObject::GetAssayData(object = obj, assay = assay, layer = slotX))
#     stopifnot(nrN == nrO)
#   }
#   return(obj)
# }


# .check_and_rename <- function(obj, assay, newnames, layer.name) {
#   cat(layer.name, fill = TRUE)

#   length_newnames <- length(newnames)
#   length_orig_names <- length(Features(obj, assay = assay))

#   stopifnot(
#     is(obj, "Seurat"),
#     is.character(assay),
#     is.character(layer.name),
#     is.character(newnames),
#     length_orig_names == length_newnames
#   )

#   assayobj <- obj@assays[[assay]]
#   feature.list <- rownames(assayobj@features@.Data)

#   if (length(feature.list) == length(newnames)) {
#     rownames(assayobj@features@.Data) <- newnames
#     nrX <- length(rownames(assayobj@features@.Data))
#   } else {
#     paste0("length feature.list", length(feature.list), "length newnames", length(newnames))
#     stop()
#   }

#   if (layer.name %in% SeuratObject::Layers(assayobj)) {
#     matrix_n <- SeuratObject::LayerData(assayobj, layer = layer.name)
#     nr1 <- nrow(matrix_n)

#     if (all(dim(matrix_n)) > 0) {
#       # browser()
#       stopifnot(nrow(matrix_n) == length(newnames))

#       if ("dgCMatrix" %in% class(matrix_n)) {
#         message(assay, "@", layer.name, " is of type dgeCMatrix!")
#         matrix_n@Dimnames[[1]] <- newnames
#       } else if ("matrix" %in% class(matrix_n)) {
#         message(assay, "@", layer.name, " is of type Matrix!")
#         rownames(matrix_n) <- newnames
#       } else if ("data.frame" %in% class(matrix_n)) {
#         message(assay, "@", layer.name, " is of type data.frame!")
#         rownames(matrix_n) <- newnames
#       } else {
#         warning(">>> No renaming: ", assay, "@", layer.name,
#           " not of type dgeCMatrix / Matrix / data.frame.",
#           immediate. = TRUE
#         )
#       }
#       stopifnot(nr1 == nrow(matrix_n))

#       SeuratObject::LayerData(assayobj, layer = layer.name) <- matrix_n
#       nr3 <- nrow(SeuratObject::LayerData(assayobj, layer = layer.name))
#       stopifnot(nr3 == nrX)
#     }
#   } else {
#     warning(paste(">>>", assay, "@", layer.name, "does not exist!"), immediate. = TRUE)
#   }
#   # obj <- SetAssayData(obj, layer = layer.name, new.data = matrix_n)
#   obj@assays[[assay]] <- assayobj
#   return(obj)
# }

# DefaultAssay(M2) <- "RNA"


# newnames <- human2mouse$human_name[match(rownames(M2),human2mouse$mouse_name)]
# newnames[is.na(newnames)] <- rownames(M2)[is.na(newnames)]

# RenameGenesSeurat(M2,newnames = newnames,assay='RNA')



RNA <- M2@assays$RNA['counts']
newnames <- human2mouse$human_name[match(rownames(RNA),human2mouse$mouse_name)]
newnames[is.na(newnames)] <- rownames(RNA)[is.na(newnames)]
rownames(RNA) <- newnames


M4 <- CreateAssayObject(RNA)
M2[['humanized']] <- M4
DefaultAssay(M2) <- "humanized"
M2[['SCT']] <- NULL

M2 <- SCTransform(M2, vst.flavor = "v2",assay='humanized')
M2 <- RunPCA(M2, npcs = 50, verbose = FALSE)

anchors <- FindTransferAnchors(
  reference = M3,
  query = M2,
  normalization.method = "SCT",
  recompute.residuals=FALSE,
  reference.reduction = "harmony",
  dims = 1:30
)

predictions <- TransferData(anchorset = anchors, refdata = M3$Subnames, dims = 1:30)


M3 <- RunUMAP(M3, dims = 1:30, reduction = 'harmony', return.model = TRUE)
M2 <- MapQuery(anchorset = anchors, reference = M3, query = M2,
	refdata = list(celltype = "Subnames"), reference.reduction = "pca", reduction.model = "umap")

score <- MappingScore(anchors,ndim = 30)

M2$map_score <- score

p1 <- DimPlot(M3, reduction = "umap", group.by = "Subnames", label = TRUE, label.size = 3, repel = TRUE,raster=TRUE,pt.size=1.5) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(M2, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE, label.size = 3, pt.size=1.5,repel = TRUE,raster=TRUE) + NoLegend() + ggtitle("Query transferred labels")
pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'RV_PAB_CM_ref_mapped.pdf'), width=10, height=5)
p1 + p2
dev.off()


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'PAB_CM_ref_mapped.pdf'), width=5, height=5)
PlotEmbedding(M2,group.by='predicted.celltype',reduction = "ref.umap",point_size=0.2,plot_under=TRUE,plot_theme=umap_theme()+NoLegend(),raster_dpi=400,raster_scale=0.5)
dev.off()


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'RV_CM_ref_mapped.pdf'), width=5, height=5)
PlotEmbedding(M3,group.by='Subnames',point_size=0.2,plot_under=TRUE,plot_theme=umap_theme()+NoLegend(),raster_dpi=400,raster_scale=0.5)
dev.off()


FeaturePlot(M2,'map_score',reduction = "ref.umap")

#######################################
#############  FIGURE S6J  ############
#######################################


RV_marks <- FindAllMarkers(M3)

RV_marks_sig <- subset(RV_marks, p_val_adj<0.05 & avg_log2FC>0)

idx<-match(RV_marks_sig$cluster,unique(RV_marks_sig$cluster))

marks <- split(RV_marks_sig$gene,idx)

M2 <- AddModuleScore(M2,marks,name='RV_marks',ctrl=25)
DefaultAssay(M3) <- "SCT"
M3 <- AddModuleScore(M3,marks,name='RV_marks')

M2$predicted.celltype <- factor(M2$predicted.celltype,levels=levels(M3))

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'PAB_CM_marker_scores.pdf'), width=6, height=4)

DotPlot(M2,c('RV_marks1','RV_marks2','RV_marks3','RV_marks4','RV_marks5'
	,'RV_marks6','RV_marks7','RV_marks8','RV_marks9','RV_marks10'),
	col.min = 0,group.by='predicted.celltype') +
scale_x_discrete(labels=c('Cm1','Cm2','Cm3','Cm4','Cm5','Cm6','Cm7','Cm8','Cm9','Cm10'))
dev.off()


DotPlot(M3,c('RV_marks1','RV_marks2','RV_marks3','RV_marks4','RV_marks5'
	,'RV_marks6','RV_marks7','RV_marks8','RV_marks9','RV_marks10'),col.min = 0) +
scale_x_discrete(labels=c('Cm1','Cm2','Cm3','Cm4','Cm5','Cm6','Cm7','Cm8','Cm9','Cm10'))

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'PAB_CM_scores.pdf'), width=3.75, height=4.25)
VlnPlot(M2,'map_score',group.by='predicted.celltype',pt.size=0)
dev.off()

#######################################
#############  FIGURE S6L  ############
#######################################

M3 <- SetIdent(M3,value = 'group')
M2 <- SetIdent(M2,value = 'group')


a <- FindMarkers(M3,ident.1='RVF',ident.2='NF')
b <-  FindMarkers(M2,ident.1='Sev',ident.2='Nor')
shared <- intersect(rownames(a),rownames(b))
dataset <- data.frame(RV=a[shared,]$avg_log2FC,PAB=b[shared,]$avg_log2FC)
rownames(dataset) <- shared
labs <- rownames(dataset)
#labs[abs(dataset$PAB - dataset$RV)<1] <- NA


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'PAB_vs_RV_CM__dot.pdf'), width=6, height=8)
ggplot(dataset, aes(x = RV, y=PAB)) + geom_point() + 
  geom_text_repel(label=labs,max.overlaps=50) + theme_classic()
dev.off()


a <- FindMarkers(M3,ident.1='RVF',ident.2='NF')
b <-  FindMarkers(M2,ident.1='Mod',ident.2='Nor')
shared <- intersect(rownames(a),rownames(b))
dataset <- data.frame(RV=a[shared,]$avg_log2FC,PAB=b[shared,]$avg_log2FC)
rownames(dataset) <- shared
labs <- rownames(dataset)
#labs[abs(dataset$PAB - dataset$RV)<1] <- NA


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'PAB_vs_RV_CM_Mod_Nor_dot.pdf'), width=6, height=8)
ggplot(dataset, aes(x = RV, y=PAB)) + geom_point() + 
  geom_text_repel(label=labs,max.overlaps=50) + theme_classic()
dev.off()


a <- FindMarkers(M2,ident.1='Sev',ident.2='Nor')
b <-  FindMarkers(M2,ident.1='Mod',ident.2='Nor')
shared <- intersect(rownames(a),rownames(b))
dataset <- data.frame(RV=a[shared,]$avg_log2FC,PAB=b[shared,]$avg_log2FC)
rownames(dataset) <- shared
labs <- rownames(dataset)
#labs[abs(dataset$PAB - dataset$RV)<1] <- NA


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'PAB_CM_Sev_Nor_Mod_Nor_dot.pdf'), width=6, height=8)
ggplot(dataset, aes(x = RV, y=PAB)) + geom_point() + 
  geom_text_repel(label=labs,max.overlaps=50) + theme_classic()
dev.off()


a <- FindMarkers(M2,ident.1='Sev',ident.2='Mod')
b <-  FindMarkers(M2,ident.1='Mod',ident.2='Nor')
shared <- intersect(rownames(a),rownames(b))
dataset <- data.frame(RV=a[shared,]$avg_log2FC,PAB=b[shared,]$avg_log2FC)
rownames(dataset) <- shared
labs <- rownames(dataset)
labs[abs(dataset$PAB - dataset$RV)<1] <- NA


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'PAB_CM_Sev_Mod_Mod_Nor_dot.pdf'), width=6, height=8)
ggplot(dataset, aes(x = RV, y=PAB)) + geom_point() + 
  geom_text_repel(label=labs,max.overlaps=10) + theme_classic()
dev.off()
#######################################
#############  FIGURE S6M  ############
#######################################

consensus_modules <- read.csv("~/Downloads/hdWGCNA_TOM/bulk_heart_modules.csv")
consensus_modules <- consensus_modules[,1:3]
consensus_modules <- subset(consensus_modules, gene_name %in% rownames(M2))
consensus_modules <- consensus_modules[match(unique(consensus_modules$gene_name), consensus_modules$gene_name),]
library(dplyr)
mapping <- labels2colors(1:100)

score_calc <- consensus_modules %>% group_by(module) %>% group_split()
module_colors <- unique(unlist(lapply(score_calc,'[[','module')))
module_colors <- paste0('M',match(module_colors,mapping))


M2 <- AddModuleScore(M2,lapply(score_calc,'[[','gene_name'),name="module_score",ctrl = 50)

cols_current <- colnames(M2@meta.data)
cols_current[startsWith(colnames(M2@meta.data),'module_score')] <- paste0('module_',module_colors)
colnames(M2@meta.data) <- cols_current


M2 <- SetIdent(M2, value = "group")
M2$group <- factor(M2$group,levels=c('Nor','Mod','Sev'))


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'PAB_seurat_dot_CM.pdf'), width=5, height=2.5)

p <- DotPlot(M2,paste0('module_',
  c('M2','M10','M12','M25','M26','M28')),dot.min=0,col.min=0,col.max=2,group.by='group') +
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
#############  FIGURE S6N  ############
#######################################
library("readxl")

Human.Mito <- read_excel("~/Downloads/hdWGCNA_TOM/Human.MitoCarta3.0.xls", sheet = "A Human MitoCarta3.0")
Mouse.Mito <- read_excel("~/Downloads/hdWGCNA_TOM/Mouse.MitoCarta3.0.xls", sheet = "A Mouse MitoCarta3.0")


M3 <- AddModuleScore(M3,list(Human.Mito$Symbol),name='mito')
M2 <- AddModuleScore(M2,list(union(Human.Mito$Symbol,Mouse.Mito$Symbol)),name='mito',ctrl = 50)

p1<-VlnPlot(M2,'mito1',group.by='group',pt.size=0)
p2<-VlnPlot(M3,'mito1',group.by='group',pt.size=0)

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'PAB_RV_CM_mitocarto.pdf'), width=3, height=4)

p1 / p2
dev.off()

anno <- trimws(unlist(lapply(lapply(lapply(Human.Mito$MitoCarta3.0_MitoPathways,str_split,'>'),'[[',1),'[[',1)))

anno[anno == "Small molecule transport | Signaling"] = "Small molecule transport" 


a<-FindMarkers(M3,ident.1='RVF',ident.2='NF',features = intersect(Human.Mito$Symbol,rownames(M3)))
a$p_val_adj[a$p_val_adj < 1e-50] = 1e-50

keyvals <- anno[match(rownames(a),Human.Mito$Symbol)]
library(colormap)

colors <- colormap(colormap=colormaps$rainbow_soft, nshades=length(unique(keyvals)))

colors <- colors[match(keyvals,unique(keyvals))]
names(colors) <- keyvals


library(EnhancedVolcano)
pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'CM_RV_Mito.pdf'), width=6, height=11)

EnhancedVolcano(a,lab=rownames(a),
  x='avg_log2FC',y='p_val_adj',
  FCcutoff = 0.1,pCutoff=0.05,colCustom = colors,xlim=c(-6.25,3),ylim=c(0,52))
dev.off()
