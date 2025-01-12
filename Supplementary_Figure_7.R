library(Seurat)
library(hdWGCNA)
library(ggeasy)
library(harmony)



source('~/Downloads/hdWGCNA_TOM/spatial_functions.R')

#######################################
#############  FIGURE S7A  ############
#######################################

M1 <- readRDS('~/Downloads/hdWGCNA_TOM/Kory_Peds_Hearts/objects/all_data.rds')
M1$Names <- M1$cell.type


M1 <- SetIdent(M1, value = "Names")


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'Peds_snUMAP.pdf'), width=5, height=5)
PlotEmbedding(M1,group.by='Names',point_size=1,plot_under=TRUE,plot_theme=umap_theme()+NoLegend(),raster_dpi=400,raster_scale=0.5)
dev.off()

#######################################
#############  FIGURE S7B  ############
#######################################

all_percent_cell <- as.data.frame(table(M1@active.ident)/length(M1@active.ident)*100)

NF_percent_cell <- cbind(as.data.frame(table(subset(M1,condition=="NF")@active.ident)/length(subset(M1,condition=="NF")@active.ident)*100),type = "NF")
NF_percent_cell$sum <- (rev(cumsum(rev(NF_percent_cell$Freq))) - NF_percent_cell$Freq/2)/100
NF_percent_cell$Freq <- NF_percent_cell$Freq/100


Donor_percent_cell <- cbind(as.data.frame(table(subset(M1,condition=="Donor")@active.ident)/length(subset(M1,condition=="Donor")@active.ident)*100),type = "Donor")
Donor_percent_cell$sum <- (rev(cumsum(rev(Donor_percent_cell$Freq))) - Donor_percent_cell$Freq/2)/100
Donor_percent_cell$Freq <- Donor_percent_cell$Freq/100


HF_percent_cell <- cbind(as.data.frame(table(subset(M1,condition=="SystolicHF")@active.ident)/length(subset(M1,condition=="SystolicHF")@active.ident)*100),type = "SystolicHF")
HF_percent_cell$sum <- (rev(cumsum(rev(HF_percent_cell$Freq))) - HF_percent_cell$Freq/2)/100
HF_percent_cell$Freq <- HF_percent_cell$Freq/100


percent_cell_df <- rbind(NF_percent_cell,Donor_percent_cell,HF_percent_cell)

ggplot(percent_cell_df, aes(fill=Var1, y=Freq, x=type,label=round(sum,1))) +  geom_bar(position="stack", stat="identity",width=0.6) + theme_classic() + xlab("Disease State") + ylab("Frequency") + labs(fill="Cell type",color='black') + theme(text = element_text(size=20),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),legend.text=element_text(color="black")) + scale_y_continuous(expand=c(0,0)) + geom_label_repel(aes(type,sum,label=scales::percent(round(Freq,2))),fill=NA,nudge_x=0.5,direction="y")


#Prevalence comparisons
cells <- table(M1@active.ident,M1@meta.data$sample)
cells <- cells[c('Cardiomyocyte','Endothelium','Fibroblast','Myeloid','Adipocyte','Pericyte','SMC'),]
cells <- sweep(cells,2,colSums(cells),'/')
cells <- data.frame(cells)


#ggboxplot(cells[77:1,],x="Var2",y="Freq",fill="Var2",group="Var2")+theme_classic() + theme(axis.text.x=element_text(size=16),axis.text.y=element_text(size=16),axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),legend.title=element_text(size=16),legend.text=element_text(size=16),text=element_text(color='black'),axis.text=element_text(color='black')) + labs(color='Group',x="Disease",y='Frequency') + facet_wrap(~Var1,ncol=7)

cells$group = rep(M1$condition[match(unique(M1@meta.data$sample),M1@meta.data$sample)],each=7)
cells$group = factor(cells$group,levels=c('Donor','NF','SystolicHF'))

#my_comparisons <- list( c("NF", "pRV"),c("pRV", "RVF"),c("NF", "RVF"))

library(ggpubr)
pdf('~/Downloads/hdWGCNA_TOM/peds_clust_freq.pdf',width=12.5,height=5)
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
#############  FIGURE S7C  ############
#######################################
M1 <- readRDS('~/Downloads/hdWGCNA_TOM/Kory_Peds_Hearts/objects/all_data.rds')

