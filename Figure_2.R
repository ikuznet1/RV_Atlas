library(Seurat)
library(hdWGCNA)
library(ggeasy)


source('~/Downloads/hdWGCNA_TOM/spatial_functions.R')



#######################################
#############  FIGURE 2A  #############
#######################################

M1<-readRDS('/Volumes/Extreme SSD/Final_Analysis/CellTypes/Post_R3_FINAL_with_counts.rds')

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'snUMAP.pdf'), width=5, height=5)
PlotEmbedding(M1,group.by='Names',point_size=1,plot_under=TRUE,plot_theme=umap_theme()+NoLegend(),raster_dpi=400,raster_scale=0.5)
dev.off()


#######################################
#############  FIGURE 2B  #############
#######################################

feats <- c("PLIN1","RYR2","LEPR","VWF","HAS1","DCN","CCL21","CSF1R","NRXN1","CD2","PDGFRB","MYH11")
p<-VlnPlot(M1, features = feats,ncol=1,pt.size=F,group.by="Names")



for(i in 1:12) {  
	 p[[i]] <- p[[i]] + NoLegend() + easy_remove_axes(which="y",what = c("ticks", "text","line")) + ggtitle("") + ylab(feats[i])
	 if(i<12){p[[i]]<-p[[i]]+easy_remove_axes(which="x")}
	 
}

#for(i in 1:12) {  
#	 p[[i]] <- p[[i]] + NoLegend() + easy_remove_axes(which="y",what = c("ticks", "text","line")) + ggtitle("") + ylab(feats[i])
#	 if(i<11){p[[i]]<-p[[i]]+easy_remove_axes(which="x")}
#	 
#}

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'sn_Vln.pdf'), width=4, height=20)
p
dev.off()


#######################################
#############  FIGURE 2C  #############
#######################################
M1 <- SetIdent(M1, value = "Names")



all_percent_cell <- as.data.frame(table(M1@active.ident)/length(M1@active.ident)*100)

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

pdf('~/Downloads/CM_prev_stacked.pdf',width=5,height=5)
ggplot(percent_cell_df, aes(fill=Var1, y=Freq, x=type,label=round(sum,1))) +  geom_bar(position="stack", stat="identity",width=0.6) + theme_classic() + xlab("Disease State") + ylab("Frequency") + labs(fill="Cell type",color='black') + theme(text = element_text(size=20),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),legend.text=element_text(color="black")) + scale_y_continuous(expand=c(0,0)) + geom_label_repel(aes(type,sum,label=scales::percent(round(Freq,2))),fill=NA,nudge_x=0.5,direction="y")
dev.off()


#Prevalence comparisons
cells <- table(M1@active.ident,M1@meta.data$patient)
cells <- cells[c('CM','FB','Myeloid','EC','Adipo','PC','SM'),]
cells <- sweep(cells,2,colSums(cells),'/')
cells <- data.frame(cells)


#ggboxplot(cells[77:1,],x="Var2",y="Freq",fill="Var2",group="Var2")+theme_classic() + theme(axis.text.x=element_text(size=16),axis.text.y=element_text(size=16),axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),legend.title=element_text(size=16),legend.text=element_text(size=16),text=element_text(color='black'),axis.text=element_text(color='black')) + labs(color='Group',x="Disease",y='Frequency') + facet_wrap(~Var1,ncol=7)


cells$group = rep(c("RVF","pRV","RVF","NF","pRV","pRV","RVF","NF","NF","pRV","NF"),each=7)


#my_comparisons <- list( c("NF", "pRV"),c("pRV", "RVF"),c("NF", "RVF"))

library(ggpubr)
pdf('~/Downloads/hdWGCNA_TOM/clust_freq.pdf',width=12.5,height=5)
p <- ggboxplot(cells[77:1,],x="group",y="Freq",fill="group",group="group")+
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
	#stat_compare_means(aes(group=group),comparisons=my_comparisons,method="t.test",ref.group="NF")+
	stat_compare_means(aes(group=group),method="anova")
p
dev.off()






#######################################
#############  FIGURE 2D  #############
#######################################
seurat_ref <- readRDS('~/Downloads/hdWGCNA_TOM/scWGCNA_bulk2sn_projection.rds')


seurat_ref <- SetActiveWGCNA(seurat_ref , 'bulk2sn')

mapping <- labels2colors(1:100)
MEs <- GetMEs(seurat_ref, harmonized=TRUE)
mods <- colnames(MEs); mods <- mods[mods != 'grey']
mods_num <- paste0('M',match(mods,mapping))
prv_vs_rv_signif <- c('M3','M4','M5','M10','M11','M12','M14')
all_signif <- c('M1','M2','M3','M4','M5','M8','M10','M11','M12','M14','M20','M25','M26','M28')


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
seurat_ref<-readRDS('/Volumes/Extreme SSD/Final_Analysis/CellTypes/Post_R3_FINAL_with_counts.rds')
seurat_ref <- SetIdent(seurat_ref, value = "Names")
seurat_ref@meta.data <- cbind(seurat_ref@meta.data, MEs)



seurat_ref <- AddModuleScore(seurat_ref,lapply(score_calc,'[[','gene_name'),name="module_score")


cols_current <- colnames(seurat_ref@meta.data)
cols_current[startsWith(colnames(seurat_ref@meta.data),'module_score')] <- paste0('module_',module_colors)
colnames(seurat_ref@meta.data) <- cols_current


p <- DotPlot(seurat_ref,paste0('module_',all_signif),group.by='Names',dot.min=0,col.min=0,col.max=2,idents=c("CM","EC","FB","Myeloid","PC","SM")) +
  RotatedAxis() + ylab('')+ xlab('')+
  scale_color_gradient2(high='red', mid='grey95', low='blue') +
  theme(
    panel.border = element_rect(size=1,fill=NA, color='black'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()
) 


source('~/Downloads/hdWGCNA_TOM/spatial_functions.R')
library(matrixStats)
#mods_use <- mapping[as.integer(unlist(lapply(strsplit(all_signif,'M'),'[[',2)))]
mods_use <- paste0('module_',all_signif)
seurat_ref@meta.data$dummy <- 1


#for (i in paste0('module_',all_signif)){
#	print(i)
#	seurat_ref@meta.data[i]<-(seurat_ref@meta.data[i] - colMeans(seurat_ref@meta.data[i]))/colSds(as.matrix(seurat_ref@meta.data[i]))
#}

#SM
p <- DotPlot(seurat_ref,paste0('module_',c('M20')),group.by='Names',dot.min=0,col.min=0,col.max=2,idents=c("CM","EC","FB","Myeloid","PC","SM")) +
  RotatedAxis() + ylab('')+ xlab('')+
  scale_color_gradient2(high='red', mid='grey95', low='blue') +
  theme(
    panel.border = element_rect(size=1,fill=NA, color='black'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()
)

#PC
p <- DotPlot(seurat_ref,paste0('module_',c('M11')),group.by='Names',dot.min=0,col.min=0,col.max=2,idents=c("CM","EC","FB","Myeloid","PC","SM")) +
  RotatedAxis() + ylab('')+ xlab('')+
  scale_color_gradient2(high='red', mid='grey95', low='blue') +
  theme(
    panel.border = element_rect(size=1,fill=NA, color='black'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()
)  

#Myeloid
p <- DotPlot(seurat_ref,paste0('module_',c('M1','M3','M4','M8')),group.by='Names',dot.min=0,col.min=0,col.max=2,idents=c("CM","EC","FB","Myeloid","PC","SM")) +
  RotatedAxis() + ylab('')+ xlab('')+
  scale_color_gradient2(high='red', mid='grey95', low='blue') +
  theme(
    panel.border = element_rect(size=1,fill=NA, color='black'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()
) 

#CM specific
p <- DotPlot(seurat_ref,paste0('module_',c('M2','M5','M12','M25','M26')),group.by='Names',dot.min=0,col.min=0,col.max=2,idents=c("CM","EC","FB","Myeloid","PC","SM")) +
  RotatedAxis() + ylab('')+ xlab('')+
  scale_color_gradient2(high='red', mid='grey95', low='blue') +
  theme(
    panel.border = element_rect(size=1,fill=NA, color='black'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()
) 

#Mixed
p <- DotPlot(seurat_ref,paste0('module_',c('M10','M14','M28')),group.by='Names',dot.min=0,col.min=0,col.max=2,idents=c("CM","EC","FB","Myeloid","PC","SM")) +
  RotatedAxis() + ylab('')+ xlab('')+
  scale_color_gradient2(high='red', mid='grey95', low='blue') +
  theme(
    panel.border = element_rect(size=1,fill=NA, color='black'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()
) 

#All
pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'sc_seurat_dit_bulk2rv.pdf'), width=7, height=4)

p <- DotPlot(seurat_ref,paste0('module_',
	c('M20','M5','M1','M3','M4','M8','M2','M12','M25','M26','M10','M28','M14','M11')),
group.by='Names',dot.min=0,col.min=0,col.max=2,idents=c("CM","EC","FB","Myeloid","PC","SM")) +
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
#############  FIGURE 2E  #############
#######################################







library(ggpubr)
p <- custom_vln(
    seurat_ref,
    features = c('M1','M2'),
    group.by = 'dummy',
    groups = c('1'),
    add_boxplot=FALSE,
    split.by = 'group',
    selected_split = c('NF','pRV','RVF'),
    split_colors=c('darkorchid', 'grey','royalblue'),
    add_colorbar=FALSE,
    plot_ymin = NA,
    pval_y_adjust=0.7,
  )


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'sc_seurat_trend.pdf'), width=6.75, height=2.5)

p <- DotPlot(subset(snLV,Names=='CM'),paste0('module_',
	c('M20','M5','M1','M3','M4','M8','M2','M12','M25','M26','M10','M28','M14','M11')),
	group.by='group',dot.min=0,col.min=0,col.max=2) +
  RotatedAxis() + ylab('')+ xlab('')+
  scale_color_gradient2(high='red', mid='grey95', low='blue') +
  theme(
    panel.border = element_rect(size=1,fill=NA, color='black'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()
) 
p
dev.off()



pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'sc_seurat_vln_stack.pdf'), width=5, height=10)
p
dev.off()


VlnPlot(seurat_ref,paste0('module_',c('M1','M2','M12','M25','M26')),pt.size=0,group.by='group')






library(Nebulosa)
plot_density(seurat_ref, "module_M2")


#######################################
#############  FIGURE 2E  #############
#######################################





#######################################
#############  FIGURE 2F  #############
#######################################