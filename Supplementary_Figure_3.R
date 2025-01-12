library(Seurat)
library(hdWGCNA)
library(ggeasy)
library(harmony)



source('~/Downloads/hdWGCNA_TOM/spatial_functions.R')


#######################################
############  FIGURE S3A  #############
#######################################
snLV <- readRDS('~/Downloads/hdWGCNA_TOM/Kory_reprocessed_all.rds')
snLV <- subset(snLV,tech=='SN')
snLV<-subset(snLV,Names=="Myeloid")

all_signif <- c('M1','M2','M3','M4','M5','M8','M10','M11','M12','M14','M20','M25','M26','M28')

consensus_modules <- read.csv("~/Downloads/hdWGCNA_TOM/bulk_heart_modules.csv")
consensus_modules <- consensus_modules[,1:3]
consensus_modules <- subset(consensus_modules, gene_name %in% rownames(snLV))
# remove duplicate gene names
consensus_modules <- consensus_modules[match(unique(consensus_modules$gene_name), consensus_modules$gene_name),]


library(dplyr)
mapping <- labels2colors(1:100)

score_calc <- consensus_modules %>% group_by(module) %>% group_split()
module_colors <- unique(unlist(lapply(score_calc,'[[','module')))
module_colors <- paste0('M',match(module_colors,mapping))


snLV <- AddModuleScore(snLV,lapply(score_calc,'[[','gene_name'),name="module_score")

cols_current <- colnames(snLV@meta.data)
cols_current[startsWith(colnames(snLV@meta.data),'module_score')] <- paste0('module_',module_colors)
colnames(snLV@meta.data) <- cols_current

which(module_colors %in% c('M1','M3','M4','M8'))



snLV <- SetIdent(snLV, value = "condition")
pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'LV_vln_modules.pdf'), width=6, height=3)
VlnPlot(snLV,c('module_M1','module_M3','module_M4','module_M8'),pt.size=0,ncol=4)
dev.off()

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'LV_dot_modules.pdf'), width=4.5, height=2.2)
p <- DotPlot(snLV,paste0('module_',
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
############  FIGURE S3B  #############
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


M1$SubNames_group <- paste0(M1$Subsubnames,'_',M1$group)
M1 <- SetIdent(M1, value = "group")
M1$group <- factor(M1$group,levels=c("NF","pRV","RVF"))

consensus_modules <- read.csv("~/Downloads/hdWGCNA_TOM/bulk_heart_modules.csv")
consensus_modules <- consensus_modules[,1:3]
consensus_modules <- subset(consensus_modules, gene_name %in% rownames(M1))
# remove duplicate gene names
consensus_modules <- consensus_modules[match(unique(consensus_modules$gene_name), consensus_modules$gene_name),]


library(dplyr)
mapping <- labels2colors(1:100)

score_calc <- consensus_modules %>% group_by(module) %>% group_split()
module_colors <- unique(unlist(lapply(score_calc,'[[','module')))
module_colors <- paste0('M',match(module_colors,mapping))


M1 <- AddModuleScore(M1,lapply(score_calc,'[[','gene_name'),name="module_score")

cols_current <- colnames(M1@meta.data)
cols_current[startsWith(colnames(M1@meta.data),'module_score')] <- paste0('module_',module_colors)
colnames(M1@meta.data) <- cols_current

which(module_colors %in% c('M1','M3','M4','M8'))

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'RV_vln_modules.pdf'), width=6, height=3)
VlnPlot(M1,c('module_M1','module_M3','module_M4','module_M8'),pt.size=0,ncol=4)
dev.off()

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'RV_dot_modules.pdf'), width=4.5, height=2.2)
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
############  FIGURE S3C  #############
#######################################


#Reference mao LV to RV
scLV <- readRDS('~/Downloads/hdWGCNA_TOM/Kory_reprocessed_all.rds')
scLV <- subset(scLV,tech=='SC')
scLV<-subset(scLV,Names=="Myeloid")


M2 <- SplitObject(M1, split.by = "patient")
M2<-PrepSCTIntegration(M2)
features<-SelectIntegrationFeatures(M2)
M2.anchors<-FindIntegrationAnchors(M2,normalization.method = 'SCT',anchor.features = features, reduction = "rpca")
M2 <- IntegrateData(anchorset = M2.anchors,normalization.method='SCT')

DefaultAssay(M2) <- "integrated"

M2 <- RunPCA(M2, npcs = 50, verbose = FALSE)
M2 <- RunUMAP(M2, reduction = "pca", dims = 1:30)


human2mouse <- read.csv('~/Downloads/hdWGCNA_TOM/human2mouse.csv',header=F)
idx <- match(unique(human2mouse[,2]),human2mouse[,2])
human2mouse<-human2mouse[idx,]
colnames(human2mouse) <-c('human_name', 'mouse_name')

RNA <- snLV@assays$RNA['counts']
newnames <- human2mouse$human_name[match(rownames(RNA),human2mouse$mouse_name)]
newnames[is.na(newnames)] <- rownames(RNA)[is.na(newnames)]
rownames(RNA) <- newnames


#M4 <- CreateAssayObject(RNA)
#snLV[['humanized']] <- M4
DefaultAssay(snLV) <- "RNA"
snLV[['SCT']] <- NULL

snLV <- SCTransform(snLV, vst.flavor = "v2",assay='RNA')
snLV <- RunPCA(snLV, npcs = 50, verbose = FALSE)

anchors <- FindTransferAnchors(
  reference = M2,
  query = snLV,
  normalization.method = "SCT",
  recompute.residuals=FALSE,
  reference.reduction = "pca",
  dims = 1:50
)

predictions <- TransferData(anchorset = anchors, refdata = M2$Subsubnames, dims = 1:50)

M2 <- RunUMAP(M2, dims = 1:50, return.model = TRUE)
snLV <- MapQuery(anchorset = anchors, reference = M2, query = snLV,
  refdata = list(celltype = "Subsubnames"), reference.reduction = "pca", reduction.model = "umap")

score <- MappingScore(anchors)

snLV$map_score <- score

p1 <- DimPlot(M2, reduction = "umap", group.by = "Subsubnames", label = TRUE, label.size = 3, repel = TRUE,raster=TRUE,pt.size=1.5) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(snLV, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE, label.size = 3, pt.size=1.5,repel = TRUE,raster=TRUE) + NoLegend() + ggtitle("Query transferred labels")
pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'RV_LV_myeloid_ref_mapped.pdf'), width=10, height=5)
p1 + p2
dev.off()




#RNA <- scLV@assays$RNA['counts']

#M4 <- CreateAssayObject(RNA)
#scLV[['humanized']] <- M4
DefaultAssay(scLV) <- "RNA"
scLV[['SCT']] <- NULL

scLV <- SCTransform(scLV, vst.flavor = "v2",assay='RNA')
scLV <- RunPCA(scLV, npcs = 50, verbose = FALSE)

anchors <- FindTransferAnchors(
  reference = M2,
  query = scLV,
  normalization.method = "SCT",
  recompute.residuals=FALSE,
  reference.reduction = "pca",
  dims = 1:50
)

predictions <- TransferData(anchorset = anchors, refdata = M2$Subsubnames, dims = 1:50)

#M2 <- RunUMAP(M2, dims = 1:50, return.model = TRUE)
scLV <- MapQuery(anchorset = anchors, reference = M2, query = scLV,
  refdata = list(celltype = "Subsubnames"), reference.reduction = "pca", reduction.model = "umap")

score <- MappingScore(anchors)

scLV$map_score <- score

p1 <- DimPlot(M2, reduction = "umap", group.by = "Subsubnames", label = TRUE, label.size = 3, repel = TRUE,raster=TRUE,pt.size=1.5) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(scLV, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE, label.size = 3, pt.size=1.5,repel = TRUE,raster=TRUE) + NoLegend() + ggtitle("Query transferred labels")
pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'RV_LV_sc_myeloid_ref_mapped.pdf'), width=10, height=5)
p1 + p2
dev.off()


scLV <- AddModuleScore(scLV,list(c('FCGR3A','LILRA5','LST1')),name='Nonclassical_Mono')
scLV <- AddModuleScore(scLV,list(c('CD14','S100A8','S100A9','S100A12','FCN1')),name='Classical_Mono')
scLV <- AddModuleScore(scLV,list(c('FCN1','OLR1','PLAUR','TRAF1')),name='Intermediate_Mono')
scLV <- AddModuleScore(scLV,list(c('MKI67','STMN1','BIRC5','TOP2A')),name='Prolif')
scLV <- AddModuleScore(scLV,list(c('TREM2','SPP1','FABP5','LGALS3')),name='Mac1')
scLV <- AddModuleScore(scLV,list(c('FOLR2','LYVE1','MRC1','SIGLEC1','CD163')),name='Mac2')
scLV <- AddModuleScore(scLV,list(c('LYVE1','HSPH1','HSPA1A','HSPA1B')),name='Mac3')
scLV <- AddModuleScore(scLV,list(c('CCL3','CCL4','PHLDA1','PMAIP1')),name='Mac4')
scLV <- AddModuleScore(scLV,list(c('KLF2','KLF4','EGR1','RHOB')),name='Mac5')
scLV <- AddModuleScore(scLV,list(c('CD1C','CCR7','FCER1A')),name='DC')


scLV <- SetIdent(scLV,value='predicted.celltype')
FeaturePlot(scLV,reduction='ref.umap',label=TRUE,'Classical_Mono1')
FeaturePlot(scLV,reduction='ref.umap',label=TRUE,'Nonclassical_Mono1')
FeaturePlot(scLV,reduction='ref.umap',label=TRUE,'Intermediate_Mono1')
FeaturePlot(scLV,reduction='ref.umap',label=TRUE,'Prolif1')
FeaturePlot(scLV,reduction='ref.umap',label=TRUE,'Mac11')
FeaturePlot(scLV,reduction='ref.umap',label=TRUE,'Mac21')
FeaturePlot(scLV,reduction='ref.umap',label=TRUE,'Mac31')
FeaturePlot(scLV,reduction='ref.umap',label=TRUE,'Mac41')
FeaturePlot(scLV,reduction='ref.umap',label=TRUE,'Mac51')
FeaturePlot(scLV,reduction='ref.umap',label=TRUE,'DC1')


DotPlot(scLV,c('Classical_Mono1','Nonclassical_Mono1','Intermediate_Mono1',
  'Prolif1','Mac11','Mac21','Mac31','Mac41','Mac51','DC1'))


scLV <- RunPCA(scLV, npcs=50, verbose = TRUE)
scLV <- FindNeighbors(scLV, dims = 1:50)
scLV <- FindClusters(scLV, resolution = 1)
scLV <- RunUMAP(scLV, dims = 1:50)



DotPlot(scLV,c('Classical_Mono1','Nonclassical_Mono1','Intermediate_Mono1',
  'Prolif1','Mac11','Mac21','Mac31','Mac41','Mac51','DC1'))



new.cluster.ids <- c('Mac4','Mac2','Mac4','Mac1','Intermediate_Mono','Classical_Mono','Prolif',
  'DC','Nonclassical_Mono','Mac3/5','Mac1','Mac4','Intermediate_Mono')

names(new.cluster.ids) <- levels(scLV)
scLV <- RenameIdents(scLV, new.cluster.ids)
scLV$LV_names <- scLV@active.ident

snLV <- RunPCA(snLV, npcs=50, verbose = TRUE)
snLV <- FindNeighbors(snLV, dims = 1:50)
snLV <- FindClusters(snLV, resolution = 0.6)
snLV <- RunUMAP(snLV, dims = 1:50)

snLV <- AddModuleScore(snLV,list(c('FCGR3A','LILRA5','LST1')),name='Nonclassical_Mono')
snLV <- AddModuleScore(snLV,list(c('CD14','S100A8','S100A9','S100A12','FCN1')),name='Classical_Mono')
snLV <- AddModuleScore(snLV,list(c('FCN1','OLR1','PLAUR','TRAF1')),name='Intermediate_Mono')
snLV <- AddModuleScore(snLV,list(c('MKI67','STMN1','BIRC5','TOP2A')),name='Prolif')
snLV <- AddModuleScore(snLV,list(c('TREM2','SPP1','FABP5','LGALS3')),name='Mac1')
snLV <- AddModuleScore(snLV,list(c('FOLR2','LYVE1','MRC1','SIGLEC1','CD163')),name='Mac2')
snLV <- AddModuleScore(snLV,list(c('LYVE1','HSPH1','HSPA1A','HSPA1B')),name='Mac3')
snLV <- AddModuleScore(snLV,list(c('CCL3','CCL4','PHLDA1','PMAIP1')),name='Mac4')
snLV <- AddModuleScore(snLV,list(c('KLF2','KLF4','EGR1','RHOB')),name='Mac5')
snLV <- AddModuleScore(snLV,list(c('CD1C','CCR7','FCER1A')),name='DC')


DotPlot(snLV,c('Classical_Mono1','Nonclassical_Mono1','Intermediate_Mono1',
  'Prolif1','Mac11','Mac21','Mac31','Mac41','Mac51','DC1'))


new.cluster.ids <- c('Mac2/3','Mac2/3','Mac2/3','Mac2/3','Mac2/3',
  'DC','Classical_Mono','Mac2/3','Mac5','Mac2/3',
  'Mac2/3','Prolif','Intermediate_Mono','Nonclassical_Mono','Mac2',
  'Mac1','Mac2/3','Mac2/3','Mac2/3')

names(new.cluster.ids) <- levels(snLV)
snLV <- RenameIdents(snLV, new.cluster.ids)
snLV$LV_names <- snLV@active.ident



p1 <- DimPlot(snLV, reduction = "ref.umap", label = TRUE, label.size = 3, repel = TRUE,raster=TRUE,pt.size=1.5) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(scLV, reduction = "ref.umap", label = TRUE, label.size = 3, pt.size=1.5,repel = TRUE,raster=TRUE) + NoLegend() + ggtitle("Query transferred labels")
pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'LV_sn_sc_myeloid_ref_mapped.pdf'), width=10, height=5)
p1 + p2
dev.off()


#######################################
############  FIGURE S3D  #############
#######################################

snLV <- SetIdent(snLV,value='condition')

a<-FindMarkers(snLV,ident.1='DCM',ident.2='Donor')

M1 <- SetIdent(M1,value='group')


b<-FindMarkers(M1,ident.1='RVF',ident.2='NF')

shared <- intersect(rownames(a),rownames(b))
dataset <- data.frame(RV=b[shared,]$avg_log2FC,LV=a[shared,]$avg_log2FC)
rownames(dataset) <- shared
labs <- rownames(dataset)
#labs[abs(dataset$PAB - dataset$RV)<1] <- NA


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'RV_vs_LV_myeloid__dot.pdf'), width=6, height=8)
ggplot(dataset, aes(x = RV, y=LV)) + geom_point() + 
  geom_text_repel(label=labs,max.overlaps=15) + theme_classic()
dev.off()



#######################################
############  FIGURE S3D  #############
#######################################


gluc_response <- "VIT;VKORC1L1;ERRFI1;AHCYL1;STEAP4;TRAF3IP2;GALNT15;SERPINE1;JADE1;SLA;CBLB;MT1X;EPS8;CCND3;BMPER;RASSF4;RPS6KA2;ANPEP;C1RL;MAP3K6;IL6R;PDGFRA;MLIP;SCARA5;IL1R1;EBF1;TTC7A;CRISPLD2;SPARCL1;FKBP5;NNMT;LPAR1;SLC1A3;PLA2G5;NID1;ACACB;ZFP36L2;PIK3R5;C3;SCFD2;LPXN;HACL1;SRGAP2;SLC38A2;SLC19A2;S100A10;KLHL29;GADD45B;ZBTB16;ELL2;CORO2B;IGF2R;NFATC4;DERA;SULT1B1;MAFB;BCL6;TMEM236;TBXAS1;NDUFAF2;RGL3;SERPINA3;MCFD2;PTPRS;ELN;PTEN;FMN1;HIF3A;TFCP2L1;PTH1R;SYNE3;CTSS;PTPRG;RNF157;ADAMTS2;C1QTNF1;IMPA2;SH3PXD2B;FLVCR2;EFHD1;AOX1;CERS6;ZHX3;KLF13;ANXA2;IFNGR1;GPX3;NCOA3;SLC39A11;NGF;OSMR;SLC39A14;TGFBR2;TGFBR3;PSMA6;ARHGAP10;MMP14;TBC1D2;SLC7A7;SLC7A8;GFOD1;DPYD;PICK1;FAM20C;COL6A3;PLIN2;ITGA5;MOCS1;ERGIC1;TMEM45A;KANK1;C1S;ADCY3;TFPI;FSTL1;TMEM165;HDAC7;KIAA0513;MTHFD1L;CLMN;PTK2B;PTPN18;GALNT6;GSN;NEGR1;TPK1;CCDC57;TXNRD1;GSR;SUSD1;LHFPL2;MERTK;KLF9;IL18R1"
gluc_response <- strsplit(gluc_response[1],';')[[1]]




M1 <- AddModuleScore(M1,list(gluc_response),name='nr3c1')
snLV <- AddModuleScore(snLV,list(gluc_response),name='nr3c1')
scLV <- AddModuleScore(scLV,list(gluc_response),name='nr3c1')



pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'RV_myeloid_nr3c1.pdf'), width=3, height=3)

VlnPlot(M1,'nr3c11',group.by='group',pt.size=0)
dev.off()

snLV$condition <- factor(snLV$condition,levels=c('Donor','DCM'))
scLV$condition <- factor(scLV$condition,levels=c('Donor','DCM'))


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'LV_sn_myeloid_nr3c1.pdf'), width=3, height=3)
VlnPlot(snLV,'nr3c11',group.by='condition',pt.size=0)
dev.off()


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'LV_sc_myeloid_nr3c1.pdf'), width=3, height=3)
VlnPlot(scLV,'nr3c11',group.by='condition',pt.size=0)
dev.off()




mhcII <- c('CIITA','CD74','HLA-DMA','HLA-DMB','HLA-DOA',
  'HLA-DOB','HLA-DPA1','HLA-DPB1','HLA-DQA1','HLA-DQA2',
  'HLA-DQB1','HLA-DQB2','HLA-DRA','HLA-DRB1','HLA-DRB3',
  'HLA-DRB4','HLA-DRB5')


M1 <- AddModuleScore(M1,list(mhcII),name='mhcII')
snLV <- AddModuleScore(snLV,list(mhcII),name='mhcII')
scLV <- AddModuleScore(scLV,list(mhcII),name='mhcII')



pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'RV_myeloid_mhcII.pdf'), width=3, height=3)

VlnPlot(M1,'mhcII1',group.by='group',pt.size=0)
dev.off()

snLV$condition <- factor(snLV$condition,levels=c('Donor','DCM'))
scLV$condition <- factor(scLV$condition,levels=c('Donor','DCM'))


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'LV_sn_myeloid_mhcII.pdf'), width=3, height=3)
VlnPlot(snLV,'mhcII1',group.by='condition',pt.size=0)
dev.off()


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'LV_sc_myeloid_mhcII.pdf'), width=3, height=3)
VlnPlot(scLV,'mhcII1',group.by='condition',pt.size=0)
dev.off()





#######################################
############  FIGURE S3E  #############
#######################################
snLV <- readRDS('~/Downloads/hdWGCNA_TOM/Kory_reprocessed_all.rds')
snLV <- subset(snLV,tech=='SN')
snLV<-subset(snLV,Names=="Fibroblasts")

snLV <- SetIdent(snLV, value = "condition")


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'LV_FB_dot_myoFb.pdf'), width=4, height=3.5)
p <- DotPlot(snLV,
  c('ACTA2','CDH11','TAGLN','SLIT3','MINDY2','MYO1B','LIMS2','GARS'),dot.min=0,col.min=0,col.max=2) +
  RotatedAxis() + ylab('')+ xlab('')+
  scale_color_gradient2(high='red', mid='grey95', low='blue') +
  theme(
    panel.border = element_rect(size=1,fill=NA, color='black'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()
) + coord_flip()
p
dev.off()

#######################################
############  FIGURE S3F  #############
#######################################


M1 <- readRDS(file = "/Volumes/Extreme\ SSD/Final_Analysis/CellTypes/fb_subclust.rds")

new.cluster.ids <- c("Fb1","Fb2","Fb3","Fb4","Fb5","Fb6","Fb7")
names(new.cluster.ids) <- levels(M1)
M1 <- RenameIdents(M1, new.cluster.ids)

M1$Subnames <- M1@active.ident
M1$SubNames_Groups <- paste(M1$Subnames,M1$group,sep='_')



M1 <- SetIdent(M1, value = "group")
M1$group <- factor(M1$group,levels=c("NF","pRV","RVF"))



#Reference mao LV to RV

M2 <- SplitObject(M1, split.by = "patient")
M2<-PrepSCTIntegration(M2)
features<-SelectIntegrationFeatures(M2)
M2.anchors<-FindIntegrationAnchors(M2,normalization.method = 'SCT',anchor.features = features, reduction = "rpca")
M2 <- IntegrateData(anchorset = M2.anchors,normalization.method='SCT')

DefaultAssay(M2) <- "integrated"

M2 <- RunPCA(M2, npcs = 50, verbose = FALSE)
M2 <- RunUMAP(M2, reduction = "pca", dims = 1:30)



DefaultAssay(snLV) <- "RNA"
snLV[['SCT']] <- NULL

snLV <- SCTransform(snLV, vst.flavor = "v2",assay='RNA')
snLV <- RunPCA(snLV, npcs = 50, verbose = FALSE)

anchors <- FindTransferAnchors(
  reference = M2,
  query = snLV,
  normalization.method = "SCT",
  recompute.residuals=FALSE,
  reference.reduction = "pca",
  dims = 1:50
)

predictions <- TransferData(anchorset = anchors, refdata = M2$Subnames, dims = 1:50)

M2 <- RunUMAP(M2, dims = 1:50, return.model = TRUE)
snLV <- MapQuery(anchorset = anchors, reference = M2, query = snLV,
  refdata = list(celltype = "Subnames"), reference.reduction = "pca", reduction.model = "umap")

score <- MappingScore(anchors)

snLV$map_score <- score

p1 <- DimPlot(M2, reduction = "umap", group.by = "Subnames", label = TRUE, label.size = 3, repel = TRUE,raster=TRUE,pt.size=1.5) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(snLV, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE, label.size = 3, pt.size=1.5,repel = TRUE,raster=TRUE) + NoLegend() + ggtitle("Query transferred labels")
pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'RV_LV_fb_ref_mapped.pdf'), width=10, height=5)
p1 + p2
dev.off()




snLV <- AddModuleScore(snLV,list(c('ACSM3','SCN7A','ABCA10','NEGR1','ABCA9')),name='Fb1')
snLV <- AddModuleScore(snLV,list(c('PCOLCE2','IGFBP6','MFAP5','S100A10','FGFBP2')),name='Fb2')
snLV <- AddModuleScore(snLV,list(c('GPX3','APOD','C3','HSPA1A','GLUL')),name='Fb3')
snLV <- AddModuleScore(snLV,list(c('PLA2G2A','RARRES1','IGFBP4','FGF7')),name='Fb4')
snLV <- AddModuleScore(snLV,list(c('ELN','GPC6','FGF14','ITGA1')),name='Fb5')
snLV <- AddModuleScore(snLV,list(c('TNC','FN1','MEOX1')),name='Fb6')
snLV <- AddModuleScore(snLV,list(c('CCL2','THBS1','CYR61','NR4A1')),name='Fb7')
snLV <- AddModuleScore(snLV,list(c('THBS4','AEBP1','POSTN','CLU','COMP')),name='Fb8')
snLV <- AddModuleScore(snLV,list(c('SERPINE1','CYR61','NFATC2','LRRFIP1')),name='Fb9')


snLV <- SetIdent(snLV,value='predicted.celltype')
FeaturePlot(snLV,reduction='ref.umap',label=TRUE,'Fb11')
FeaturePlot(snLV,reduction='ref.umap',label=TRUE,'Fb81')
FeaturePlot(snLV,reduction='ref.umap',label=TRUE,'Fb91')

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'RV_LV_fb_ref_mapped_dot.pdf'), width=6, height=4)

DotPlot(snLV,c('Fb31','Fb61','Fb91','Fb41','Fb71','Fb21',
  'Fb51','Fb81','Fb11'),group.by = "predicted.celltype",col.min=0,col.max=2)
dev.off()


#######################################
############  FIGURE S3G  #############
#######################################



snLV <- SetIdent(snLV,value='condition')

a<-FindMarkers(snLV,ident.1='DCM',ident.2='Donor')

M1 <- SetIdent(M1,value='group')


b<-FindMarkers(M1,ident.1='RVF',ident.2='NF')

shared <- intersect(rownames(a),rownames(b))
dataset <- data.frame(RV=b[shared,]$avg_log2FC,LV=a[shared,]$avg_log2FC)
rownames(dataset) <- shared
labs <- rownames(dataset)
#labs[abs(dataset$PAB - dataset$RV)<1] <- NA


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'RV_vs_LV_fb_dot.pdf'), width=6, height=8)
ggplot(dataset, aes(x = RV, y=LV)) + geom_point() + 
  geom_text_repel(label=labs,max.overlaps=10) + theme_classic()
dev.off()












scLV <- RunPCA(scLV, npcs=50, verbose = TRUE)
scLV <- FindNeighbors(scLV, dims = 1:50)
scLV <- FindClusters(scLV, resolution = 1)
scLV <- RunUMAP(scLV, dims = 1:50)



DotPlot(scLV,c('Classical_Mono1','Nonclassical_Mono1','Intermediate_Mono1',
  'Prolif1','Mac11','Mac21','Mac31','Mac41','Mac51','DC1'))



new.cluster.ids <- c('Mac4','Mac2','Mac4','Mac1','Intermediate_Mono','Classical_Mono','Prolif',
  'DC','Nonclassical_Mono','Mac3/5','Mac1','Mac4','Intermediate_Mono')

names(new.cluster.ids) <- levels(scLV)
scLV <- RenameIdents(scLV, new.cluster.ids)
scLV$LV_names <- scLV@active.ident

snLV <- RunPCA(snLV, npcs=50, verbose = TRUE)
snLV <- FindNeighbors(snLV, dims = 1:50)
snLV <- FindClusters(snLV, resolution = 0.6)
snLV <- RunUMAP(snLV, dims = 1:50)

snLV <- AddModuleScore(snLV,list(c('FCGR3A','LILRA5','LST1')),name='Nonclassical_Mono')
snLV <- AddModuleScore(snLV,list(c('CD14','S100A8','S100A9','S100A12','FCN1')),name='Classical_Mono')
snLV <- AddModuleScore(snLV,list(c('FCN1','OLR1','PLAUR','TRAF1')),name='Intermediate_Mono')
snLV <- AddModuleScore(snLV,list(c('MKI67','STMN1','BIRC5','TOP2A')),name='Prolif')
snLV <- AddModuleScore(snLV,list(c('TREM2','SPP1','FABP5','LGALS3')),name='Mac1')
snLV <- AddModuleScore(snLV,list(c('FOLR2','LYVE1','MRC1','SIGLEC1','CD163')),name='Mac2')
snLV <- AddModuleScore(snLV,list(c('LYVE1','HSPH1','HSPA1A','HSPA1B')),name='Mac3')
snLV <- AddModuleScore(snLV,list(c('CCL3','CCL4','PHLDA1','PMAIP1')),name='Mac4')
snLV <- AddModuleScore(snLV,list(c('KLF2','KLF4','EGR1','RHOB')),name='Mac5')
snLV <- AddModuleScore(snLV,list(c('CD1C','CCR7','FCER1A')),name='DC')


DotPlot(snLV,c('Classical_Mono1','Nonclassical_Mono1','Intermediate_Mono1',
  'Prolif1','Mac11','Mac21','Mac31','Mac41','Mac51','DC1'))


new.cluster.ids <- c('Mac2/3','Mac2/3','Mac2/3','Mac2/3','Mac2/3',
  'DC','Classical_Mono','Mac2/3','Mac5','Mac2/3',
  'Mac2/3','Prolif','Intermediate_Mono','Nonclassical_Mono','Mac2',
  'Mac1','Mac2/3','Mac2/3','Mac2/3')

names(new.cluster.ids) <- levels(snLV)
snLV <- RenameIdents(snLV, new.cluster.ids)
snLV$LV_names <- snLV@active.ident



p1 <- DimPlot(snLV, reduction = "ref.umap", label = TRUE, label.size = 3, repel = TRUE,raster=TRUE,pt.size=1.5) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(scLV, reduction = "ref.umap", label = TRUE, label.size = 3, pt.size=1.5,repel = TRUE,raster=TRUE) + NoLegend() + ggtitle("Query transferred labels")
pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'LV_sn_sc_myeloid_ref_mapped.pdf'), width=10, height=5)
p1 + p2
dev.off()


































bulk_modules <- read.csv("~/Downloads/hdWGCNA_TOM/bulk_heart_modules.csv")
bulk_modules$module <- match(bulk_modules$module,mapping)
combined_set_RV <- data.frame()
mods_idx <- c(1,3,4,8)
for (i in mods_idx){
  key_genes <- subset(bulk_modules,module %in% c(i))$gene_name
  key_genes <- key_genes[key_genes %in% rownames(M1)]
  gene_set <- FindMarkers(M1, ident.1 = "RVF", ident.2 = "NF",features=key_genes)
  #gene_set<-subset(gene_set,p_val_adj<0.05)
  gene_set$module <- paste0('M',i)
  gene_set$color <- mapping[i]
  if (length(combined_set_RV) == 0){
    combined_set_RV <- gene_set
  }
  else {
    combined_set_RV <- rbind(combined_set_RV,gene_set)
  }
}

snLV <- readRDS('~/Downloads/hdWGCNA_TOM/Kory_reprocessed_all.rds')
snLV <- subset(snLV,tech=='SN')
snLV<-subset(snLV,Names=="Myeloid")
snLV <- SetIdent(snLV, value = "condition")
snLV <- PrepSCTFindMarkers(snLV)


bulk_modules <- read.csv("~/Downloads/hdWGCNA_TOM/bulk_heart_modules.csv")
bulk_modules$module <- match(bulk_modules$module,mapping)
combined_set <- data.frame()
mods_idx <- c(1,3,4,8)
for (i in mods_idx){
  key_genes <- subset(bulk_modules,module %in% c(i))$gene_name
  key_genes <- key_genes[key_genes %in% rownames(snLV)]
  gene_set <- FindMarkers(snLV, ident.1 = "DCM", ident.2 = "Donor",features=key_genes,recorrect_umi=FALSE)
  #gene_set<-subset(gene_set,p_val_adj<0.05)
  gene_set$module <- paste0('M',i)
  gene_set$color <- mapping[i]
  if (length(combined_set) == 0){
    combined_set <- gene_set
  }
  else {
    combined_set <- rbind(combined_set,gene_set)
  }
}




overlap_genes = intersect(rownames(M1),rownames(snLV))
gene_set <- FindMarkers(M1, ident.1 = "RVF", ident.2 = "NF",features=overlap_genes)
gene_set$module <- paste0('All')
gene_set$color <- 'grey'
combined_set_RV <- rbind(combined_set_RV,gene_set)


gene_set <- FindMarkers(snLV, ident.1 = "DCM", ident.2 = "Donor",features=overlap_genes,recorrect_umi=FALSE)
gene_set$module <- paste0('All')
gene_set$color <- 'grey'
combined_set <- rbind(combined_set,gene_set)

a<-subset(combined_set,module=='M1')
b <- subset(combined_set_RV,module=='M1')
shared <- intersect(rownames(a),rownames(b))
dataset <- data.frame(LV=a[shared,]$avg_log2FC,RV=b[shared,]$avg_log2FC)
rownames(dataset) <- shared

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'LV_vs_RV_Myeloid_module_M1_dot.pdf'), width=6, height=8)
ggplot(dataset, aes(x = RV, y=LV)) + geom_point() + 
  geom_text_repel(label=rownames(dataset)) + theme_classic()
dev.off()

a<-subset(combined_set,module=='M3')
b <- subset(combined_set_RV,module=='M3')
shared <- intersect(rownames(a),rownames(b))
dataset <- data.frame(LV=a[shared,]$avg_log2FC,RV=b[shared,]$avg_log2FC)
rownames(dataset) <- shared

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'LV_vs_RV_Myeloid_module_M3_dot.pdf'), width=6, height=8)
ggplot(dataset, aes(x = RV, y=LV)) + geom_point() + 
  geom_text_repel(label=rownames(dataset)) + theme_classic()
dev.off()


a<-subset(combined_set,module=='M4')
b <- subset(combined_set_RV,module=='M4')
shared <- intersect(rownames(a),rownames(b))
dataset <- data.frame(LV=a[shared,]$avg_log2FC,RV=b[shared,]$avg_log2FC)
rownames(dataset) <- shared

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'LV_vs_RV_Myeloid_module_M4_dot.pdf'), width=6, height=8)
ggplot(dataset, aes(x = RV, y=LV)) + geom_point() + 
  geom_text_repel(label=rownames(dataset)) + theme_classic()
dev.off()

a<-subset(combined_set,module=='M8')
b <- subset(combined_set_RV,module=='M8')
shared <- intersect(rownames(a),rownames(b))
dataset <- data.frame(LV=a[shared,]$avg_log2FC,RV=b[shared,]$avg_log2FC)
rownames(dataset) <- shared

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'LV_vs_RV_Myeloid_module_M8_dot.pdf'), width=6, height=8)
ggplot(dataset, aes(x = RV, y=LV)) + geom_point() + 
  geom_text_repel(label=rownames(dataset)) + theme_classic()
dev.off()


a<-subset(combined_set,module=='All')
b <- subset(combined_set_RV,module=='All')
shared <- intersect(rownames(a),rownames(b))
dataset <- data.frame(LV=a[shared,]$avg_log2FC,RV=b[shared,]$avg_log2FC)
rownames(dataset) <- shared
labs <- rownames(dataset)
labs[abs(dataset$LV - dataset$RV)<0.1] <- NA


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'LV_vs_RV_CM_module_All_dot.pdf'), width=6, height=8)
ggplot(dataset, aes(x = RV, y=LV)) + geom_point() + 
  geom_text_repel(label=labs,max.overlaps=10) + theme_classic()
dev.off()





#genes_diff <- abs(dataset$LV - dataset$RV) > 1
#genes_diff <- shared[genes_diff]

#genes_diff_sig <- a[genes_diff,]$p_val_adj < 0.05 & b[genes_diff,]$p_val_adj < 0.05
#genes_diff<-genes_diff[genes_diff_sig]

#LV_diff_up <- genes_diff[a[genes_diff,]$avg_log2FC > 0]
#RV_diff_up <- genes_diff[b[genes_diff,]$avg_log2FC > 0]


#######################################
#############  FIGURE 4D  #############
#######################################
library(viridis)
a<-subset(combined_set,module %in% c('M10','M25','M26','M28'))
b <- subset(combined_set_RV,module %in% c('M10','M25','M26','M28'))
shared <- intersect(rownames(a),rownames(b))
dataset <- data.frame(LV=a[shared,]$avg_log2FC,RV=b[shared,]$avg_log2FC)
rownames(dataset) <- shared
mito_genes <- rownames(subset(dataset,LV<0 & RV<0))

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

enriched <- enrichr(mito_genes, dbs)
pdf('~/Downloads/hdWGCNA_TOM/CM_RV_LV_mito_enrichr_up.pdf',width=5,height=2.5)
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

a<-subset(combined_set,module %in% c('M2'))
b <- subset(combined_set_RV,module %in% c('M2'))
shared <- intersect(rownames(a),rownames(b))
dataset <- data.frame(LV=a[shared,]$avg_log2FC,RV=b[shared,]$avg_log2FC)
rownames(dataset) <- shared
M2_genes <- rownames(subset(dataset,LV>0 & RV>0))

enriched <- enrichr(M2_genes, dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('~/Downloads/hdWGCNA_TOM/CM_RV_LV_M2_enrichr_up.pdf',width=5,height=2.5)
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

pdf('~/Downloads/hdWGCNA_TOM/CM_RV_LV_M&mito2_enrichr_up.pdf',width=6,height=5)
p1/p2
dev.off()

#######################################
#############  FIGURE 4F  #############
#######################################


M1 <- readRDS(file = "/Volumes/Extreme\ SSD/Final_Analysis/CellTypes/cm_subclust.rds")

new.cluster.ids <- c("Cm1","Cm2","Cm3","Cm4","Cm5","Cm6","Cm7","Cm8","Cm9","Cm10")
names(new.cluster.ids) <- levels(M1)
M1 <- RenameIdents(M1, new.cluster.ids)

M1$Subnames <- M1@active.ident
M1$SubNames_Groups <- paste(M1$Subnames,M1$group,sep='_')

M1 <- SetIdent(M1, value = "group")


aggregate_snLV <- AggregateExpression(snLV, group.by = c("orig.ident",'condition'), return.seurat = TRUE,assays='RNA')

aggregate_RV <- AggregateExpression(M1, group.by = c("patient",'group'), return.seurat = TRUE,assays='RNA')

LV_counts <- aggregate_snLV$RNA['counts']
RV_counts <- aggregate_RV$RNA['counts']
LV_meta <- aggregate_snLV@meta.data
RV_meta <- aggregate_RV@meta.data

colnames(LV_meta) <- c('patient','group')
colnames(RV_meta) <- c('temp','patient','group')
RV_meta$temp <- NULL
LV_meta$origin <- 'LV'
RV_meta$origin <- 'RV'


meta <- data.frame(rbind(LV_meta,RV_meta))

genes_shared <- intersect(rownames(LV_counts),rownames(RV_counts))

countmat<-cbind(LV_counts[genes_shared,],RV_counts[genes_shared,])

meta_merge <- meta
meta_merge$group[meta_merge$group == 'Donor'] = 'NF'


library("DESeq2")
ddsSE <- DESeqDataSetFromMatrix(countData = countmat,
                              colData = meta_merge,
                              design = ~origin+group)

ddsSE <- estimateSizeFactors(ddsSE)
idx <- rowSums( counts(ddsSE, normalized=TRUE) >= 5 ) >= 3
ddsSE <- ddsSE[idx,]

normalized_counts <- counts(ddsSE, normalized=TRUE)

library("sva")

mod  <- model.matrix(~origin+group, colData(ddsSE))
mod0 <- model.matrix(~group, colData(ddsSE))
svseq <- svaseq(normalized_counts, mod, mod0)


meta_merge$SV1 <- svseq$sv[,1]
meta_merge$SV2 <- svseq$sv[,2]
meta_merge$SV3 <- svseq$sv[,3]
meta_merge$SV4 <- svseq$sv[,4]
meta_merge$SV5 <- svseq$sv[,5]
meta_merge$SV6 <- svseq$sv[,6]
meta_merge$SV7 <- svseq$sv[,7]
meta_merge$SV8 <- svseq$sv[,8]
meta_merge$SV9 <- svseq$sv[,9]
meta_merge$SV10 <- svseq$sv[,10]
meta_merge$SV11 <- svseq$sv[,11]
meta_merge$SV12 <- svseq$sv[,12]
meta_merge$SV13 <- svseq$sv[,13]
meta_merge$SV14 <- svseq$sv[,14]

ddsSE <- DESeqDataSetFromMatrix(countData = countmat,
                              colData = meta_merge,
                              design = ~group+SV1+SV2+SV3+SV4+SV5+
                              SV6+SV7+SV8+SV9+SV10+SV11+SV12+SV13+SV14)



ddsSE <- estimateSizeFactors(ddsSE)
idx <- rowSums( counts(ddsSE, normalized=TRUE) >= 5 ) >= 12
ddsSE <- ddsSE[idx,]

normalized_counts <- counts(ddsSE, normalized=TRUE)



ddsSE <- DESeq(ddsSE)

vstSE <- vst(ddsSE,blind = FALSE)


mat <- assay(vstSE)
mm <- model.matrix(~group, colData(vstSE))
mat <- limma::removeBatchEffect(mat, covariates=colData(vstSE)[,4:18], design=mm)
assay(vstSE) <- mat

plotPCA(vstSE,intgroup=c("group"),ntop=100) + theme_classic() + theme(axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.title.x=element_text(size=24),axis.title.y=element_text(size=24),legend.title=element_text(size=24),legend.text=element_text(size=24)) + labs(color='Disease',shape='Sex') + geom_point(size=3.5)

pdf('~/Downloads/hdWGCNA_TOM/CM_RV_LV_vst_batchcorrect_PCA.pdf',width=6,height=2.5)
plotPCA(vstSE[,vstSE$group %in% c('pRV','RVF','DCM')],intgroup=c("group"),ntop=100) + theme_classic() + theme(axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.title.x=element_text(size=24),axis.title.y=element_text(size=24),legend.title=element_text(size=24),legend.text=element_text(size=24)) + labs(color='Disease',shape='Sex') + geom_point(size=3.5)
dev.off()






meta_merge_pooled <- meta_merge
meta_merge_pooled[meta_merge_pooled$group %in% c('pRV','RVF'),]$group = 'RVF'
ddsSE <- DESeqDataSetFromMatrix(countData = countmat,
                              colData = meta_merge_pooled,
                              design = ~group+SV1+SV2+SV3+SV4+SV5+
                              SV6+SV7+SV8+SV9+SV10+SV11+SV12+SV13+SV14)



ddsSE <- estimateSizeFactors(ddsSE)
idx <- rowSums( counts(ddsSE, normalized=TRUE) >= 5 ) >= 12
ddsSE <- ddsSE[idx,]
normalized_counts <- counts(ddsSE, normalized=TRUE)
ddsSE <- DESeq(ddsSE)
vstSE <- vst(ddsSE,blind = FALSE)
mat <- assay(vstSE)
mm <- model.matrix(~group, colData(vstSE))
mat <- limma::removeBatchEffect(mat, covariates=colData(vstSE)[,4:18], design=mm)
assay(vstSE) <- mat

#pdf('~/Downloads/hdWGCNA_TOM/CM_RV_LV_vst_batchcorrect_PCA_temp.pdf',width=6,height=2.5)

plotPCA(vstSE[,vstSE$group %in% c('RVF','DCM')],intgroup=c("group"),ntop=100) + theme_classic() + theme(axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.title.x=element_text(size=24),axis.title.y=element_text(size=24),legend.title=element_text(size=24),legend.text=element_text(size=24)) + labs(color='Disease',shape='Sex') + geom_point(size=3.5)

#dev.off()
dcm.vs.rvf <- lfcShrink(ddsSE,contrast=c('group','DCM','RVF'), type="ashr")

dcm.vs.rvf.signif <- subset(dcm.vs.rvf,padj<0.1)


lv <- subset(dcm.vs.rvf.signif,log2FoldChange>0.1)
rv <- subset(dcm.vs.rvf.signif,log2FoldChange<0.1)


lv.genes <- rownames(lv)
lv.genes <- lv.genes[!str_detect(lv.genes,'.1') & !str_detect(lv.genes,'.2')]
cat(lv.genes,sep='\n')

rv.genes <- rownames(rv)
rv.genes <- rv.genes[!str_detect(rv.genes,'.1') & !str_detect(rv.genes,'.2') & !str_detect(rv.genes,'.3') & !str_detect(rv.genes,'.4')]
cat(rv.genes,sep='\n')



library(EnhancedVolcano)
pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'CM_LV_RV_DESeq.pdf'), width=8, height=6)

EnhancedVolcano(subset(dcm.vs.rvf,baseMean>100),lab=rownames(subset(dcm.vs.rvf,baseMean>100)),
  x='log2FoldChange',y='padj',
  FCcutoff = 0.1,pCutoff=0.05,ylim=c(0,3)) + coord_flip()
dev.off()





pdf('~/Downloads/hdWGCNA_TOM/CM_RV_LV_vst_batchcorrect_PCA_temp.pdf',width=6,height=2.5)
FeaturePlot(M1,'ZSWIM6')

dev.off()


#######################################
#############  FIGURE 4H  #############
#######################################


arvm <- read.csv('~/Downloads/hdWGCNA_TOM/ARVM_RV_vs_LV.csv')


RV_genes <- unlist(lapply(subset(arvm,logFC>0)$SYMBOL,toupper))
LV_genes <- unlist(lapply(subset(arvm,logFC<0)$SYMBOL,toupper))



pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'arvm_RV_LV.pdf'), width=8, height=6)

EnhancedVolcano(arvm,lab=toupper(arvm$SYMBOL),
  x='logFC',y='adj.P.Val',
  FCcutoff = 0.1,pCutoff=0.05,ylim=c(0,3)) + coord_flip()
dev.off()

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

enriched <- enrichr(RV_genes, dbs)
pdf('~/Downloads/hdWGCNA_TOM/ARVM_RV_enrichr_up.pdf',width=6,height=2.5)
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


enriched <- enrichr(LV_genes, dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('~/Downloads/hdWGCNA_TOM/ARVM_LV_enrichr_up.pdf',width=5,height=2.5)
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

pdf('~/Downloads/hdWGCNA_TOM/ARVM_RV_LV_enrichr_up.pdf',width=6,height=5)
p1/p2
dev.off()


#######################################
#############  FIGURE 4I  #############
#######################################

M2.genes <- score_calc[which(module_colors %in% c("M2"))][[1]]$gene_name

M2.arvm <- arvm[toupper(arvm$SYMBOL) %in% M2.genes,]

mito.genes <- unlist(lapply(score_calc[which(module_colors %in% c("M10","M25","M26","M28"))],'[[',1))

mito.arvm <- arvm[toupper(arvm$SYMBOL) %in% mito.genes,]

M12.genes <- score_calc[which(module_colors %in% c("M12"))][[1]]$gene_name

M12.arvm <- arvm[toupper(arvm$SYMBOL) %in% M12.genes,]

