library(Seurat)
library(decoupleR)
library(dplyr)
library(tibble)
library(tidyr)
library(patchwork)
library(ggplot2)
library(pheatmap)
library(harmony)

#M1 <- readRDS(file = "/Volumes/Extreme\ SSD/Final_Analysis/CellTypes/cm_subclust.rds")
#new.cluster.ids <- c("Cm1","Cm2","Cm3","Cm4","Cm5","Cm6","Cm7","Cm8","Cm9","Cm10")
#names(new.cluster.ids) <- levels(M1)
#M1 <- RenameIdents(M1, new.cluster.ids)
#M1$Subnames <- M1@active.ident
#M1$SubNames_Groups <- paste(M1$Subnames,M1$group,sep='_')
#M1 <- SetIdent(M1, value = "group")


M1<-readRDS('./dependencies/shared/Post_R3_FINAL_with_counts.rds')
M1 <- subset(M1, Names == 'CM')



mem.maxVSize(vsize=3e12)
net <- decoupleR::get_collectri(organism = 'human', 
                                split_complexes = FALSE)

# Extract the normalized log-transformed counts
mat <- as.matrix(M1@assays$SCT@data)

# Run ulm
acts <- decoupleR::run_ulm(mat = mat, 
                           net = net, 
                           .source = 'source', 
                           .target = 'target',
                           .mor='mor', 
                           minsize = 5)
saveRDS(acts,'./dependencies/shared/TF_activity.rds')

offset <- min(acts$score)
#acts$score <- acts$score + offset

# Extract ulm and store it in tfsulm in pbmc
M1[['tfsulm']] <- acts %>%
                    tidyr::pivot_wider(id_cols = 'source', 
                                       names_from = 'condition',
                                       values_from = 'score') %>%
                    tibble::column_to_rownames('source') %>%
                    Seurat::CreateAssayObject(.)


# Change assay
DefaultAssay(object = M1) <- "tfsulm"

# Scale the data
M1 <- Seurat::ScaleData(M1)
M1@assays$tfsulm@data <- M1@assays$tfsulm@scale.data


n_tfs <- 25

# Extract activities from object as a long dataframe
df <- t(as.matrix(M1@assays$tfsulm@data)) %>%
      as.data.frame() %>%
      dplyr::mutate(cluster =factor(M1$group,levels=c('NF','pRV','RVF'))) %>%
      tidyr::pivot_longer(cols = -cluster, 
                          names_to = "source", 
                          values_to = "score") %>%
      dplyr::group_by(cluster, source) %>%
      dplyr::summarise(mean = mean(score))

# Get top tfs with more variable means across clusters
all_tfs <- df %>%
       dplyr::group_by(source) %>%
       dplyr::summarise(std = stats::sd(mean)) %>%
       dplyr::arrange(-abs(std)) 


tfs <- df %>%
       dplyr::group_by(source) %>%
       dplyr::summarise(std = stats::sd(mean)) %>%
       dplyr::arrange(-abs(std)) %>%
       head(n_tfs) %>%
       dplyr::pull(source)

# Subset long data frame to top tfs and transform to wide matrix
top_acts_mat <- df %>%
                dplyr::filter(source %in% tfs) %>%
                tidyr::pivot_wider(id_cols = 'cluster', 
                                   names_from = 'source',
                                   values_from = 'mean') %>%
                tibble::column_to_rownames('cluster') %>%
                as.matrix()

# Choose color palette
colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
colors.use <- grDevices::colorRampPalette(colors = colors)(100)

my_breaks <- c(seq(-2, 0, length.out = ceiling(100 / 2) + 1),
               seq(0.05, 2, length.out = floor(100 / 2)))

# Plot
p1 <- pheatmap::pheatmap(mat = top_acts_mat,
                   color = colors.use,
                   border_color = "white",
                   breaks = my_breaks,
                   cellwidth = 15,
                   cellheight = 15,
                   treeheight_row = 20,
                   treeheight_col = 20) 
pdf('./output/CM_TF.pdf',width=8,height=2.5)
p1
dev.off()


#AEBP1, NR3C1, PPARA, ELF5, PPARDm BATF, PHF5A, FOXG1, PLALGL1, NCOR2, BACH2, POU2AF1, NROB2, MBD2, NCOR1, DNMT3B, MESP2, SIN3A, TEAD1, TEAD4, ZNF804A, TRPS1, TEF, JARID2, PURA

# Extract activities from object as a long dataframe
df <- t(as.matrix(M1@assays$tfsulm@data)) %>%
      as.data.frame() %>%
      dplyr::mutate(cluster =factor(M1$patient,levels=c('1561','1681','1691','1697','1392','1567','1618','1692','1343','1467','1632'))) %>%
      tidyr::pivot_longer(cols = -cluster, 
                          names_to = "source", 
                          values_to = "score") %>%
      dplyr::group_by(cluster, source) %>%
      dplyr::summarise(mean = mean(score))


top_acts_mat <- df %>%
                dplyr::filter(source %in% tfs) %>%
                tidyr::pivot_wider(id_cols = 'cluster', 
                                   names_from = 'source',
                                   values_from = 'mean') %>%
                tibble::column_to_rownames('cluster') %>%
                as.matrix()


top_acts_mat <- top_acts_mat[,c('AEBP1', 'NR3C1', 'PPARA', 'ELF5', 'PPARD', 'BATF', 'PHF5A', 'FOXG1', 'PLAGL1', 'NCOR2', 'BACH2', 'POU2AF1', 'NR0B2', 'MBD2', 'NCOR1', 'DNMT3B', 'MESP2', 'SIN3A', 'TEAD1', 'TEAD4', 'ZNF804A', 'TRPS1', 'TEF', 'JARID2', 'PURA')]

# Choose color palette
colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
colors.use <- grDevices::colorRampPalette(colors = colors)(100)

my_breaks <- c(seq(-2, 0, length.out = ceiling(100 / 2) + 1),
               seq(0.05, 2, length.out = floor(100 / 2)))

# Plot
p1 <- pheatmap::pheatmap(mat = top_acts_mat,
                   color = colors.use,
                   border_color = "white",
                   breaks = my_breaks,
                   cellwidth = 15,
                   cellheight = 15,
                   treeheight_row = 20,
                   treeheight_col = 20,
                   cluster_rows = FALSE,
                   cluster_cols = FALSE) 
pdf('./output/CM_TF_patient_breakdown.pdf',width=8,height=2.5)
p1
dev.off()








n_tfs <- 25

# Extract activities from object as a long dataframe
df <- t(as.matrix(M1@assays$tfsulm@data)) %>%
      as.data.frame() %>%
      dplyr::mutate(cluster =factor(M1$patient,levels=c('1561','1681','1691','1697','1392','1567','1618','1692','1343','1467','1632'))) %>%
      tidyr::pivot_longer(cols = -cluster, 
                          names_to = "source", 
                          values_to = "score") %>%
      dplyr::group_by(cluster, source) %>%
      dplyr::summarise(mean = mean(score))

# Get top tfs with more variable means across clusters
all_tfs <- df %>%
       dplyr::group_by(source) %>%
       dplyr::summarise(std = stats::sd(mean)) %>%
       dplyr::arrange(-abs(std)) 


tfs <- df %>%
       dplyr::group_by(source) %>%
       dplyr::summarise(std = stats::sd(mean)) %>%
       dplyr::arrange(-abs(std)) %>%
       head(n_tfs) %>%
       dplyr::pull(source)

# Subset long data frame to top tfs and transform to wide matrix
top_acts_mat <- df %>%
                dplyr::filter(source %in% tfs) %>%
                tidyr::pivot_wider(id_cols = 'cluster', 
                                   names_from = 'source',
                                   values_from = 'mean') %>%
                tibble::column_to_rownames('cluster') %>%
                as.matrix()

# Choose color palette
colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
colors.use <- grDevices::colorRampPalette(colors = colors)(100)

my_breaks <- c(seq(-2, 0, length.out = ceiling(100 / 2) + 1),
               seq(0.05, 2, length.out = floor(100 / 2)))

# Plot
p1 <- pheatmap::pheatmap(mat = top_acts_mat,
                   color = colors.use,
                   border_color = "white",
                   breaks = my_breaks,
                   cellwidth = 15,
                   cellheight = 15,
                   treeheight_row = 20,
                   treeheight_col = 20,
                   cluster_rows = FALSE) 
pdf('./output/CM_TF_patient_diff.pdf',width=8,height=2.5)
p1
dev.off()


n_tfs <- 25

# Extract activities from object as a long dataframe
df <- t(as.matrix(M1@assays$tfsulm@data)) %>%
      as.data.frame() %>%
      dplyr::mutate(cluster =factor(M1$patient,levels=c('1561','1681','1691','1697','1392','1567','1618','1692','1343','1467','1632'))) %>%
      tidyr::pivot_longer(cols = -cluster, 
                          names_to = "source", 
                          values_to = "score") %>%
      dplyr::group_by(cluster, source) %>%
      dplyr::summarise(mean = mean(score))


df_t <- df %>%
      dplyr::group_by(source) %>%
      dplyr::summarise(mean = median(mean))

df_t <- subset(df_t,mean > .1)

# Get top tfs with more variable means across clusters
all_tfs <- df %>%
       dplyr::group_by(source) %>%
       dplyr::filter(source %in% df_t$source) %>%
       dplyr::summarise(std = stats::sd(mean)) %>%
       dplyr::arrange(abs(std)) 


tfs <- df %>%
       dplyr::group_by(source) %>%
       dplyr::filter(source %in% df_t$source) %>%
       dplyr::summarise(std = stats::sd(mean)) %>%
       dplyr::arrange(abs(std)) %>%
       head(n_tfs) %>%
       dplyr::pull(source)

# Subset long data frame to top tfs and transform to wide matrix
top_acts_mat <- df %>%
                dplyr::filter(source %in% tfs) %>%
                tidyr::pivot_wider(id_cols = 'cluster', 
                                   names_from = 'source',
                                   values_from = 'mean') %>%
                tibble::column_to_rownames('cluster') %>%
                as.matrix()

# Choose color palette
colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
colors.use <- grDevices::colorRampPalette(colors = colors)(100)

my_breaks <- c(seq(-2, 0, length.out = ceiling(100 / 2) + 1),
               seq(0.05, 2, length.out = floor(100 / 2)))

# Plot
p1 <- pheatmap::pheatmap(mat = top_acts_mat,
                   color = colors.use,
                   border_color = "white",
                   breaks = my_breaks,
                   cellwidth = 15,
                   cellheight = 15,
                   treeheight_row = 20,
                   treeheight_col = 20,
                   cluster_rows = FALSE) 
pdf('./output/CM_TF_patient_same.pdf',width=8,height=2.5)
p1
dev.off()








#NRF1 and NRF2 are very up
subset(df, source == 'NFE2L2')
subset(df, source == 'NRF1')


# Not uch change in PGC1a
subset(df, source == 'PPARGC1A')

# Estrogen receptors are down
subset(df, source == 'ESRRA')
subset(df, source == 'ESRRB')
subset(df, source == 'ESRRG')
subset(df, source == 'PPARA')
subset(df, source == 'PPARD')

#PPARG pretty flat
subset(df, source == 'PPARG')


p1 <- VlnPlot(M1,'ESRRA',group.by='group',pt.size=0,layer='scale.data')
p2 <- VlnPlot(M1,'ESRRB',group.by='group',pt.size=0,layer='scale.data')
p3 <- VlnPlot(M1,'ESRRG',group.by='group',pt.size=0,layer='scale.data')
p4 <- VlnPlot(M1,'PPARA',group.by='group',pt.size=0,layer='scale.data')
p5 <- VlnPlot(M1,'PPARD',group.by='group',pt.size=0,layer='scale.data')
p6 <- VlnPlot(M1,'PPARG',group.by='group',pt.size=0,layer='scale.data')
p7 <- VlnPlot(M1,'NRF1',group.by='group',pt.size=0,layer='scale.data')
p8 <- VlnPlot(M1,'NFE2L2',group.by='group',pt.size=0,layer='scale.data')
p9 <- VlnPlot(M1,'PPARGC1A',group.by='group',pt.size=0,layer='scale.data')

pdf('./output/CM_TF_VLN.pdf',width=7,height=7)
p1 + p4 + p7| p2  + p5 + p8 | p3 + p6 + p9
dev.off()

Idents(M1) <- 'group'
nf.vs.prv <- FindMarkers(M1, ident.1 = "NF", ident.2 = "pRV", verbose = FALSE,slot = 'scale.data')
rvf.vs.prv <- FindMarkers(M1, ident.1 = "RVF", ident.2 = "pRV", verbose = FALSE,slot = 'scale.data')


M2 <- AggregateExpression(M1, assays = "tfsulm", return.seurat = T, group.by = c("patient", "group"))

gene_name <- 'PPARGC1A'
a <- subset(M1,group=='NF')[['tfsulm']]@scale.data[gene_name,]
b <- subset(M1,group=='pRV')[['tfsulm']]@scale.data[gene_name,]
c <- subset(M1,group=='RVF')[['tfsulm']]@scale.data[gene_name,]
wilcox.test(a,b,alternative="two.sided")
wilcox.test(b,c,alternative="two.sided")


p1 <- VlnPlot(M2,'ESRRA',group.by='group')
p2 <- VlnPlot(M2,'ESRRB',group.by='group')
p3 <- VlnPlot(M2,'ESRRG',group.by='group')
p4 <- VlnPlot(M2,'PPARA',group.by='group')
p5 <- VlnPlot(M2,'PPARD',group.by='group')
p6 <- VlnPlot(M2,'PPARG',group.by='group')
p7 <- VlnPlot(M2,'NRF1',group.by='group')
p8 <- VlnPlot(M2,'NFE2L2',group.by='group')
p9 <- VlnPlot(M2,'PPARGC1A',group.by='group')

pdf('./outputs/CM_TF_VLN_pseudobulk.pdf',width=7,height=7)
p1 + p4 + p7| p2  + p5 + p8 | p3 + p6 + p9
dev.off()

M2 <- FindVariableFeatures(M2, selection.method = "vst", nfeatures = 50)
M2 <- RunPCA(M2,npcs=5)
M2 <- RunHarmony(M2, "patient",dims.use=1:5,nclust = 5)
M2 <- RunUMAP(M2, dims = 1:5,n.neighbors=10,reduction = "harmony")

pdf('~/Downloads/hdWGCNA_TOM/PAH/CM_TF_UMAP.pdf',width=7,height=7)
DimPlot(M2,group.by='group')
dev.off()

M1 <- FindVariableFeatures(M1, selection.method = "vst", nfeatures = 50)
M1 <- RunPCA(M1,npcs=30)
M1 <- RunHarmony(M1, "patient",dims.use=1:30)
M1 <- RunUMAP(M1, dims = 1:30,reduction = "harmony")

pdf('./outputs/CM_TF_UMAP_ss.pdf',width=7,height=7)
DimPlot(M1,group.by='group')
dev.off()










