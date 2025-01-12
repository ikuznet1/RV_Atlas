##############################################
##############################################
#### Figure S2A
##############################################
##############################################

RV_data <- readRDS('~/Downloads/hdWGCNA_TOM//RV_data.rds')
RV_data <- subset(RV_data, Names=='CM')

RV_data <- SCTransform(RV_data,vst.flavor="v2", assay = "RNA",vars.to.regress=c('nFeature_RNA','percent.mt'))
RV_data <- RunPCA(RV_data, npcs = 30)
RV_data <- RunHarmony(RV_data,'patient')
RV_data <- RunUMAP(RV_data, reduction = "harmony", dims = 1:30)
RV_data <- FindNeighbors(RV_data, reduction = "harmony", dims = 1:30)
RV_data <- FindClusters(RV_data, resolution = 0.3)

RV_data <- subset(RV_data,idents = c('4','8'),invert= T)

RV_data <- SCTransform(RV_data,vst.flavor="v2", assay = "decontXcounts",vars.to.regress=c('nFeature_RNA','percent.mt'))
RV_data <- RunPCA(RV_data, npcs = 30,assay='SCT')
RV_data <- RunHarmony(RV_data,'patient')
RV_data <- FindNeighbors(RV_data, reduction = "harmony", dims = 1:30)
RV_data <- RunUMAP(RV_data, reduction = "harmony", dims = 1:30)
RV_data <- FindClusters(RV_data, resolution = 0.3)

RV_data <- subset(RV_data,idents = c('8'),invert= T)

feats <- setdiff(rownames(RV_data),c('XIST','TTTY14','TTTY10','UTY'))
feats <- setdiff(feats,grep('^LINC', feats, value = TRUE))
feats <- grep("-", feats, invert=TRUE, value = TRUE)
feats <- grep("\\.", feats, invert=TRUE, value = TRUE)


RV_data <- SCTransform(RV_data[feats,],vst.flavor="v2", assay = "decontXcounts",vars.to.regress=c('nFeature_RNA','percent.mt'))
RV_data <- RunPCA(RV_data, npcs = 30,assay='SCT')
RV_data <- RunHarmony(RV_data,'patient')
RV_data <- FindNeighbors(RV_data, reduction = "harmony", dims = 1:30)
RV_data <- RunUMAP(RV_data, reduction = "harmony", dims = 1:30)
RV_data <- FindClusters(RV_data, resolution = 0.3)


RV.cm.marks <- FindAllMarkers(RV_data)
head(subset(RV.cm.marks,p_val_adj < 0.05 & cluster == 0))

marker.genes <- unique(subset(RV.cm.marks, p_val_adj < 0.05)$gene)


new.cluster.ids <- c("Cm1","Cm2","Cm3","Cm4","Cm5","Cm6","Cm7","Cm8","Cm9","Cm10")
names(new.cluster.ids) <- levels(RV_data)
RV_data <- RenameIdents(RV_data, new.cluster.ids)

RV_data$Subnames <- RV_data@active.ident
RV_data$SubNames_Groups <- paste(RV_data$Subnames,RV_data$group,sep='_')

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'CM_New_snUMAP.pdf'), width=5, height=5)
PlotEmbedding(RV_data,group.by='Subnames',point_size=1,plot_under=TRUE,plot_theme=umap_theme()+NoLegend(),raster_dpi=400,raster_scale=0.5)
dev.off()

#saveRDS(RV_data,'~/Downloads/hdWGCNA_TOM/cm_new_subclust.rds')
##############################################
##############################################
#### Figure S2B
##############################################
##############################################

RV_data <- AddModuleScore(RV_data, features=list(c('ACTA1','PLCE1','TNNT1')),assay="SCT",name="Clust0Score")
RV_data <- AddModuleScore(RV_data, features=list(c('TOGARAM2','CPNE5','CORIN','FGF12','CACNB2','PDE3A')),assay="SCT",name="Clust1Score")
RV_data <- AddModuleScore(RV_data, features=list(c('GRK5','TSIX')),assay="SCT",name="Clust2Score")
RV_data <- AddModuleScore(RV_data, features=list(c('ARL17B','SOX5','CCND3')),assay="SCT",name="Clust3Score")
RV_data <- AddModuleScore(RV_data, features=list(c('PRELID2','SH3RF2','GPC5','GRXCR2')),assay="SCT",name="Clust4Score")
RV_data <- AddModuleScore(RV_data, features=list(c('CCSER1','PALLD','SORBS2','CAMK2D')),assay="SCT",name="Clust5Score")
RV_data <- AddModuleScore(RV_data, features=list(c('CNN1','MYPN','C4orf54','MYO18B')),assay="SCT",name="Clust6Score")
RV_data <- AddModuleScore(RV_data, features=list(c('NPPA','NPPB')),assay="SCT",name="Clust7Score")
RV_data <- AddModuleScore(RV_data, features=list(c('HS6ST3','GRIK2','KCNJ3')),assay="SCT",name="Clust8Score")
RV_data <- AddModuleScore(RV_data, features=list(c('COX7A1','COX6A2')),assay="SCT",name="Clust9Score")


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'CM_new_sn_Dot.pdf'), width=8, height=5)
DotPlot(RV_data,features=c("Clust0Score1",
  "Clust1Score1","Clust2Score1","Clust3Score1",
  "Clust4Score1","Clust5Score1","Clust6Score1",
  "Clust7Score1","Clust8Score1","Clust9Score1"),
col.min=0,col.max=2)+ylab('Condition')+xlab('Marker Score')+
scale_x_discrete(labels=c('Clust 1','Clust 2','Clust 3',
  'Clust 4','Clust 5','Clust 6','Clust 7','Clust 8','Clust 9','Clust 10'))
dev.off()


##############################################
##############################################
#### Figure S2C
##############################################
##############################################


cm.patient <- table(RV_data$Subnames,RV_data$patient)
cm.patient <- t(t(cm.patient) / colSums(cm.patient))

disease = c("RVF","pRV","RVF","NF","pRV","pRV","RVF","NF","NF","pRV","NF")

disease <- c(t(replicate(10,disease)))

cm.patient <- data.frame(disease = disease,cm.patient)

cm.patient$disease <- factor(cm.patient$disease, levels=c('NF','pRV','RVF'))

#cm.patient$Var1 <- factor(cm.patient$Var1, levels=rev(names(sort(niche.counts))))



pdf('~/Downloads/hdWGCNA_TOM/CM_clust_counts.pdf',width=10,height=3)
ggplot(cm.patient,aes(Var1,Freq,color = disease))+geom_boxplot() + theme_classic()
dev.off()

library(ggpubr)
pdf('~/Downloads/hdWGCNA_TOM/CM_clust_freq_stats.pdf',width=12.5,height=15)
p <- ggboxplot(cm.patient,x="disease",y="Freq",fill="disease",group="disease")+
  theme_classic() + 
  theme(axis.text.x=element_text(size=16),
  axis.text.y=element_text(size=16),
  axis.title.x=element_text(size=16),
  axis.title.y=element_text(size=16),
  legend.title=element_text(size=16),
  legend.text=element_text(size=16),
  text=element_text(color='black'),
  axis.text=element_text(color='black')) + 
  facet_wrap(~Var1,ncol=5) + 
  #stat_compare_means(aes(group=group),comparisons=my_comparisons,method="t.test",ref.group="NF")+
  stat_compare_means(aes(group=disease),method="kruskal.test")
p
dev.off()



##############################################
##############################################
#### Figure S2D
##############################################
##############################################


library(enrichR)
#CAMK2D, SORBS2, TNNT1. ACTA1, PLCE1, NPPA, NPPB

M1 <- readRDS(file = "~/Downloads/hdWGCNA_TOM/cm_new_subclust.rds")
a <- FindAllMarkers(M1)


dbs <-c('GO_Biological_Process_2023','GO_Cellular_Component_2021','GO_Molecular_Function_2021','LINCS_L1000_Chem_Pert_up','LINCS_L1000_Chem_Pert_down', 'WikiPathway_2021_Human', 'KEGG_2021_Human')

# compute GO terms:
enrich_list <- list()

combined_output <- data.frame()
for (i in 1:length(unique(M1$Subnames))) {
    cur_mod <- unique(M1$Subnames)[i]
    cur_info <- subset(a, cluster == cur_mod & avg_log2FC>0 & p_val_adj < 0.05)
    cur_info <- cur_info[, c("gene")]
    
    enriched <- enrichR::enrichr(cur_info, dbs)
    Sys.sleep(5)
    for (db in names(enriched)) {
        cur_df <- enriched[[db]]
        if (nrow(cur_df) > 1) {
            cur_df$db <- db
            cur_df$module <- cur_mod
            combined_output <- rbind(combined_output, cur_df)
        }
    }
}

M1 <- SetupForWGCNA(M1,wgcna_name='temp')
M1 <- SetEnrichrTable(M1, combined_output)

outdir = '~/Downloads/hdWGCNA_TOM/scCM_subclust_enrichr_plot'


wrapText <- function(x, len) {
        sapply(x, function(y) paste(strwrap(y, len), collapse = "\n"), 
            USE.NAMES = FALSE)
}

enrichr_df <- GetEnrichrTable(M1)

for (i in 1:length(unique(M1$Subnames))) {
    cur_mod <- unique(M1$Subnames)[i]
    cur_terms <- subset(enrichr_df, module == cur_mod)
    print(cur_mod)
    #cur_color <- modules %>% subset(module == cur_mod) %>% 
    #    .$color %>% unique %>% as.character
    #if (!is.null(plot_bar_color)) {
    #    cur_color <- plot_bar_color
    #}
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
    pdf(paste0(outdir, "/", cur_mod, ".pdf"), width = 5, 
        height = 4)
    for (plot in plot_list) {
        print(plot)
    }
    dev.off()
}


#Plot enrichments by cluster
selected_terms <- subset(combined_output,db=="GO_Biological_Process_2023")
selected_terms <- subset(selected_terms, Adjusted.P.value < 0.05)



# subset selected terms
idx_top_1 <- match(unique(selected_terms$module),selected_terms$module)
idx_top_3 <- sort(c(idx_top_1,idx_top_1+1,idx_top_1+2))

selected_terms<-selected_terms[idx_top_3,]


selected_terms$group <- factor(
  as.character(selected_terms$module),
  levels = levels(selected_terms$module)
)


# set max pval
quantile(-log(selected_terms$P.value), 0.95)
max_p <- 10

selected_terms$logp <- -log(selected_terms$Adjusted.P.value)
selected_terms$logp <- ifelse(selected_terms$logp > max_p, max_p, selected_terms$logp)

# remove GO Term ID
library(stringr)
selected_terms$Term <- str_replace(selected_terms$Term, "\\(GO.*", "")

selected_terms <- selected_terms %>%
  arrange(group)


selected_terms$wrap <- wrapText(selected_terms$Term, 30)

selected_terms$Term <- factor(
  as.character(selected_terms$Term),
  levels = rev(unique(as.character(selected_terms$Term)))
)

selected_terms$wrap <- factor(
  as.character(selected_terms$wrap),
  levels = rev(unique(as.character(selected_terms$wrap)))
)

library(viridis)

# GO Term dot plot
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
    plot.margin = margin(0,0,0,0),
    panel.grid = element_line(size=0.25, color='lightgrey')
  )


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'CM_by_cluster_New_GO.pdf'), width=5.25, height=8)
p
dev.off()


##############################################
##############################################
#### Figure S2E
##############################################
##############################################

seurat_ref <- readRDS(file = "~/Downloads/hdWGCNA_TOM/cm_new_subclust.rds")
seurat_ref <- SetIdent(seurat_ref, value = "Subnames")

modules_up <- c('M2','M12','M28')
modules_down <- c('M10','M25','M26')

mapping <- labels2colors(1:100)


consensus_modules <- read.csv("~/Downloads/hdWGCNA_TOM/bulk_heart_modules.csv")
consensus_modules <- consensus_modules[,1:3]
consensus_modules <- subset(consensus_modules, gene_name %in% rownames(seurat_ref))
# remove duplicate gene names
consensus_modules <- consensus_modules[match(unique(consensus_modules$gene_name), consensus_modules$gene_name),]


library(dplyr)
score_calc <- consensus_modules %>% group_by(module) %>% group_split()
module_colors <- unique(unlist(lapply(score_calc,'[[','module')))
module_colors <- paste0('M',match(module_colors,mapping))


#MEs <- GetMEs(seurat_ref, harmonized=TRUE)
#mods <- colnames(MEs); mods <- mods[mods != 'grey']
#mods_num <- paste0('M',match(mods,mapping))

#colnames(MEs)<-paste0('M',match(colnames(MEs),mapping))


#idx <- rownames(MEs) %in% colnames(seurat_ref)

#seurat_ref@meta.data <- cbind(seurat_ref@meta.data, MEs[idx,])



seurat_ref <- AddModuleScore(seurat_ref,lapply(score_calc,'[[','gene_name'),name="module_score")


cols_current <- colnames(seurat_ref@meta.data)
cols_current[startsWith(colnames(seurat_ref@meta.data),'module_score')] <- paste0('module_',module_colors)
colnames(seurat_ref@meta.data) <- cols_current

#Dot Plot of enrichment cell type of CM enriched modules

modules_all <- c('M2','M12',"M10","M25","M26","M28")


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'CM_cluster_dot_new_subclust_mods.pdf'), width=6, height=2.5)

p <- DotPlot(seurat_ref,paste0('module_',modules_all),group.by='Subnames',dot.min=0,col.min=0) +
  RotatedAxis() + ylab('')+ xlab('')+
  scale_color_gradient2(high='red', mid='grey95', low='blue') +coord_flip()+
  theme(
    panel.border = element_rect(size=1,fill=NA, color='black'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()
) 
p
dev.off()


##############################################
##############################################
#### Figure S2F
##############################################
##############################################


Xenium.cm <- readRDS(file = "~/Downloads/hdWGCNA_TOM/Xenium/cm.rds")


Xenium.cm <- SCTransform(Xenium.cm,vst.flavor="v2", assay = "Xenium",vars.to.regress=c('nFeature_Xenium'))
Xenium.cm <- RunPCA(Xenium.cm, npcs = 30)
Xenium.cm <- RunHarmony(Xenium.cm,'patient')
Xenium.cm <- FindNeighbors(Xenium.cm, reduction = "harmony", dims = 1:30)
Xenium.cm <- RunUMAP(Xenium.cm, reduction = "harmony", dims = 1:30)
Xenium.cm <- FindClusters(Xenium.cm, resolution = 0.3)

a <- FindAllMarkers(Xenium.cm)
head(subset(a,p_val_adj < 0.05 & cluster == 0))

#3 is contaminated by EC's
#4 is contaminated by FBs
#8 is contam by PBs
#9 FB markers
#11 Toss - mixed
#13 toss
#14 toss
#15 toss
#16 toss
#18 toss
#20 toss
#21 toss

Xenium.cm <- subset(Xenium.cm,idents = c('3','4','8','9','11','13','14','15','16','18','20','21'),invert=T)

Xenium.cm <- SCTransform(Xenium.cm,vst.flavor="v2", assay = "Xenium",vars.to.regress=c('nFeature_Xenium'))
Xenium.cm <- RunPCA(Xenium.cm, npcs = 30)
Xenium.cm <- RunHarmony(Xenium.cm,'patient')
Xenium.cm <- FindNeighbors(Xenium.cm, reduction = "harmony", dims = 1:30)
Xenium.cm <- RunUMAP(Xenium.cm, reduction = "harmony", dims = 1:30)
Xenium.cm <- FindClusters(Xenium.cm, resolution = 0.3)

a <- FindAllMarkers(Xenium.cm)
head(subset(a,p_val_adj < 0.05 & cluster == 0))

#7 toss 
#12 toss
#13 toss

Xenium.cm <- subset(Xenium.cm,idents = c('7','12','13'),invert=T)
Xenium.cm <- subset(Xenium.cm,idents = c('10','15','16'),invert=T)


new.cluster.ids <- c('TNNT1/PLCE1','Mixed','MYH6','NPPA/NPPB','PALLD/ANKRD1/MYH7','NPPA/NPPB','CNN1','TIMP3','PPP1R1A','HMGCS2','ANKRD1/ACTA1')
names(new.cluster.ids) <- levels(Xenium.cm)
Xenium.cm <- RenameIdents(Xenium.cm, new.cluster.ids)

Xenium.cm$Subnames <- Xenium.cm@active.ident

Xenium.cm <- AddModuleScore(Xenium.cm, features=list(c('TNNT1','PLCE1')),assay="SCT",name="Clust0Score",ctrl=10)
Xenium.cm <- AddModuleScore(Xenium.cm, features=list(c('MYL2','TNNI3','ACTA1','MYH7')),assay="SCT",name="Clust1Score",ctrl=10)
Xenium.cm <- AddModuleScore(Xenium.cm, features=list(c('MYH6')),assay="SCT",name="Clust2Score",ctrl=10)
Xenium.cm <- AddModuleScore(Xenium.cm, features=list(c('NPPA','NPPB')),assay="SCT",name="Clust3Score",ctrl=10)
Xenium.cm <- AddModuleScore(Xenium.cm, features=list(c('PALLD','ANKRD1','MYH7')),assay="SCT",name="Clust4Score",ctrl=10)
Xenium.cm <- AddModuleScore(Xenium.cm, features=list(c('CNN1')),assay="SCT",name="Clust5Score",ctrl=10)
Xenium.cm <- AddModuleScore(Xenium.cm, features=list(c('TIMP3')),assay="SCT",name="Clust6Score",ctrl=10)
Xenium.cm <- AddModuleScore(Xenium.cm, features=list(c('PPP1R1A')),assay="SCT",name="Clust7Score",ctrl=10)
Xenium.cm <- AddModuleScore(Xenium.cm, features=list(c('HMGCS2')),assay="SCT",name="Clust8Score",ctrl=10)
Xenium.cm <- AddModuleScore(Xenium.cm, features=list(c('ANKRD1','ACTA1')),assay="SCT",name="Clust9Score",ctrl=10)


#PTN

pdf(paste0('~/Downloads/hdWGCNA_TOM/Xenium/', 'Xenium_cm_marks.pdf'), width=7, height=3.5)
DotPlot(Xenium.cm,features=c("Clust0Score1",
  "Clust1Score1","Clust2Score1","Clust3Score1",
  "Clust4Score1","Clust5Score1","Clust6Score1",
  "Clust7Score1","Clust8Score1","Clust9Score1"),
col.min=0,col.max=2,dot.min=0)+ylab('Condition')+xlab('Marker Score')+
scale_color_gradient2(high='red', mid='grey95', low='blue') +
scale_x_discrete(labels=c('1','2','3',
  '4','5','6','7','8','9','10')) +
  
  theme(
    panel.border = element_rect(size=1,fill=NA, color='black'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()
) 

dev.off()


##############################################
##############################################
#### Figure S2G
##############################################
##############################################

meta.data <- read.csv('~/Downloads/hdWGCNA_TOM/Xenium/metadata.csv')
niche_manual <- meta.data$niche_manual
names(niche_manual) <- rownames(meta.data)

Xenium.cm$niche_manual <- niche_manual[colnames(Xenium.cm)]
cm.niche <- table(Xenium.cm$niche_manual,Xenium.cm$Subnames)
cm.niche <- t(t(cm.niche)/colSums(cm.niche))

feats <- c('Endocardial CMs','Isolated CMs','Peri-immune CMs','Peri-lymphatic/venous CMs',
  'Perivascular CMs','Stromal CMs 2','Stromal CMs/THBS4 FBs')

cm.niche <- cm.niche[feats,]

mycol <- colorpanel(1000,"blue","white","red")


#saveRDS(Xenium.cm,'~/Downloads/hdWGCNA_TOM/Xenium/cm_minimalist.rds')


pdf('~/Downloads/hdWGCNA_TOM/Xenium/xenium_heatmap_niche.pdf',width=3.5,height=4.5)

heatmap.2(as.matrix(cm.niche), scale="row",
   labRow=rownames(cm.niche), 
   col=mycol, margin=c(6,6),trace="none", density.info="none", lhei=c(1,10,3),lwid=c(1,10), dendrogram='none',breaks = seq(-4, 4, length.out = 1001),
   Rowv = FALSE, Colv=FALSE,srtCol=90,lmat = rbind(c(0,3),c(2,1),c(0,4)))

dev.off()




cm.niche <- data.frame(cm.niche)





pdf('~/Downloads/hdWGCNA_TOM/Xenium/CM_niche_counts.pdf',width=10,height=3)
ggplot(cm.niche,aes(Var1,Freq,color = Var2))+geom_boxplot() + theme_classic()
dev.off()

##############################################
##############################################
#### Figure S2H
##############################################
##############################################

Xenium.cm$NPPA <- as.character(Xenium.cm$Subnames)
Xenium.cm$NPPA[Xenium.cm$NPPA != 'NPPA/NPPB'] = 'Other'

pdf('~/Downloads/hdWGCNA_TOM/Xenium/CM_niche_spatial.pdf',width=5,height=5)

ImageDimPlot(subset(Xenium.cm,niche_manual == 'Endocardial CMs'), 
  group.by = 'NPPA',fov = "fov.5",
  axes = F, dark.background=F,size=1)
dev.off()

pdf('~/Downloads/hdWGCNA_TOM/Xenium/CM_niche_spatial_all.pdf',width=5,height=5)

ImageDimPlot(Xenium.cm, 
  group.by = 'NPPA',fov = "fov.5",
  axes = F, dark.background=F,size=1)
dev.off()



pdf('~/Downloads/hdWGCNA_TOM/Xenium/CM_niche_spatial_alt.pdf',width=5,height=5)

ImageDimPlot(subset(Xenium.cm,niche_manual == 'Endocardial CMs'), 
  group.by = 'NPPA',fov = "fov.7",
  axes = F, dark.background=F,size=1)
dev.off()

pdf('~/Downloads/hdWGCNA_TOM/Xenium/CM_niche_spatial_all_alt.pdf',width=5,height=5)

ImageDimPlot(Xenium.cm, 
  group.by = 'NPPA',fov = "fov.7",
  axes = F, dark.background=F,size=1)
dev.off()