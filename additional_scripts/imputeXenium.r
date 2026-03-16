#############################
### WRITE DATA FOR PYTHON ###
#############################


library(Seurat)
library(SeuratDisk)
library(Matrix)

M1<-readRDS('~/Downloads/hdWGCNA_TOM/PAH/Post_R3_FINAL_with_counts.rds')
sc.dat <- t(GetAssayData(object = M1, assay = "decontXcounts", slot = "counts"))
rownames(sc.dat) <- colnames(M1)
rm(M1)

write.csv(sc.dat,'~/Downloads/hdWGCNA_TOM/PAH/RV_decont_raw.csv')


M1<-readRDS('~/Downloads/hdWGCNA_TOM/Xenium/collated_xenium_data.rds')
M1.loom <- as.loom(M1, filename = "~/Downloads/hdWGCNA_TOM/PAH/Xenium.loom", verbose = FALSE)
M1.loom$close_all()

# Max x is ~12Km Max & is ~23K
a1 <- GetTissueCoordinates(M1[['fov']])
a2 <- GetTissueCoordinates(M1[['fov.2']])
a3 <- GetTissueCoordinates(M1[['fov.3']])
a4 <- GetTissueCoordinates(M1[['fov.4']])
a5 <- GetTissueCoordinates(M1[['fov.5']])
a6 <- GetTissueCoordinates(M1[['fov.6']])
a7 <- GetTissueCoordinates(M1[['fov.7']])
a8 <- GetTissueCoordinates(M1[['fov.8']])
a9 <- GetTissueCoordinates(M1[['fov.9']])

a <- rbind(cbind(a1,z = 1),
	cbind(a2,z = 2),
	cbind(a3,z = 3),
	cbind(a4,z = 4),
	cbind(a5,z = 5),
	cbind(a6,z = 6),
	cbind(a7,z = 7),
	cbind(a8,z = 8),
	cbind(a9,z = 9))



write.csv(a,'~/Downloads/hdWGCNA_TOM/PAH/Xenium_coords.csv')

sc.dat <- t(GetAssayData(object = M1, assay = "Xenium", slot = "counts"))
write.csv(sc.dat,'~/Downloads/hdWGCNA_TOM/PAH/Xenium_raw.csv')

####################

# Recombine data
library(nanoparquet)
library(pillar)
library(mltools)


mem.maxVSize(vsize=3e12)


M1 <- read.csv('~/Downloads/hdWGCNA_TOM/PAH/Xenium_raw.csv',row.names = 1)
M1 <- as.data.frame(M1)
M2 <- read_parquet_info('~/Downloads/hdWGCNA_TOM/PAH/XeniumImpute/imputed_genes_unnorm.parquet')
M2 <- read_parquet_metadata('~/Downloads/hdWGCNA_TOM/PAH/XeniumImpute/imputed_genes_unnorm.parquet')
M2 <- read_parquet('~/Downloads/hdWGCNA_TOM/PAH/XeniumImpute/imputed_genes_unnorm.parquet',col_select=c(1))
M2 <- read_parquet('~/Downloads/hdWGCNA_TOM/PAH/XeniumImpute/imputed_genes_unnorm.parquet')

M2$HLA.DQB2 <- NULL
M2 <- cbind(M2,M1)

rm(M1)
gc()

write.table(M2, file="~/Downloads/hdWGCNA_TOM/PAH/XeniumImpute/full_geneset_imputed.txt", row.names=TRUE, col.names=TRUE)

#cat <(head -1 'full_geneset_imputed.txt' | cut -d' ' -f2-) <(tail -n+2 'full_geneset_imputed.txt') > /Volumes/Extreme\ SSD/RV_Paper_Revisions/full_geneset_imputed_trim.txt




library(BPCells)
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(Azimuth)


#M1 <- Matrix::readMM('/Users/ikuz/Downloads/hdWGCNA_TOM/PAH/XeniumImpute/full_geneset_imputed_sparse_1.mtx')

genes <- array(unlist(read.table('~/Downloads/hdWGCNA_TOM/PAH/XeniumImpute/genes.txt',sep=' ')))
cells <- array(unlist(read.table('~/Downloads/hdWGCNA_TOM/PAH/XeniumImpute/cells.txt',sep=' ')))


M1 <- import_matrix_market('/Volumes/Extreme\ SSD/RV_Paper_Revisions/full_geneset_imputed_sparse_1.mtx',row_names=cells[1:100000],col_names = c('index',genes))
M1 <-  M1[,2:31062]

# Write the matrix to a directory
write_matrix_dir(
   mat = M1,
   dir = '~/Downloads/hdWGCNA_TOM/PAH/XeniumImpute/seurat_data_1')

M1 <- import_matrix_market('/Volumes/Extreme\ SSD/RV_Paper_Revisions/full_geneset_imputed_sparse_2.mtx',row_names=cells[100001:200000],col_names = c('index',genes))
M1 <-  M1[,2:31062]

write_matrix_dir(
   mat = M1,
   dir = '~/Downloads/hdWGCNA_TOM/PAH/XeniumImpute/seurat_data_2')

M1 <- import_matrix_market('/Volumes/Extreme\ SSD/RV_Paper_Revisions/full_geneset_imputed_sparse_3.mtx',row_names=cells[200001:300000],col_names = c('index',genes))
M1 <-  M1[,2:31062]

write_matrix_dir(
   mat = M1,
   dir = '~/Downloads/hdWGCNA_TOM/PAH/XeniumImpute/seurat_data_3')

M1 <- import_matrix_market('/Volumes/Extreme\ SSD/RV_Paper_Revisions/full_geneset_imputed_sparse_4.mtx',row_names=cells[300001:400000],col_names = c('index',genes))
M1 <-  M1[,2:31062]

write_matrix_dir(
   mat = M1,
   dir = '~/Downloads/hdWGCNA_TOM/PAH/XeniumImpute/seurat_data_4')

M1 <- import_matrix_market('/Volumes/Extreme\ SSD/RV_Paper_Revisions/full_geneset_imputed_sparse_5.mtx',row_names=cells[400001:500000],col_names = c('index',genes))
M1 <-  M1[,2:31062]

write_matrix_dir(
   mat = M1,
   dir = '~/Downloads/hdWGCNA_TOM/PAH/XeniumImpute/seurat_data_5')

M1 <- import_matrix_market('/Volumes/Extreme\ SSD/RV_Paper_Revisions/full_geneset_imputed_sparse_6.mtx',row_names=cells[500001:600000],col_names = c('index',genes))
M1 <-  M1[,2:31062]

write_matrix_dir(
   mat = M1,
   dir = '~/Downloads/hdWGCNA_TOM/PAH/XeniumImpute/seurat_data_6')

M1 <- import_matrix_market('/Volumes/Extreme\ SSD/RV_Paper_Revisions/full_geneset_imputed_sparse_7.mtx',row_names=cells[600001:length(cells)],col_names = c('index',genes))
M1 <-  M1[,2:31062]

write_matrix_dir(
   mat = M1,
   dir = '~/Downloads/hdWGCNA_TOM/PAH/XeniumImpute/seurat_data_7')

data.list <- c()
data.list[[1]] <- t(open_matrix_dir(dir = '~/Downloads/hdWGCNA_TOM/PAH/XeniumImpute/seurat_data_1'))
data.list[[2]] <- t(open_matrix_dir(dir = '~/Downloads/hdWGCNA_TOM/PAH/XeniumImpute/seurat_data_2'))
data.list[[3]] <- t(open_matrix_dir(dir = '~/Downloads/hdWGCNA_TOM/PAH/XeniumImpute/seurat_data_3'))
data.list[[4]] <- t(open_matrix_dir(dir = '~/Downloads/hdWGCNA_TOM/PAH/XeniumImpute/seurat_data_4'))
data.list[[5]] <- t(open_matrix_dir(dir = '~/Downloads/hdWGCNA_TOM/PAH/XeniumImpute/seurat_data_5'))
data.list[[6]] <- t(open_matrix_dir(dir = '~/Downloads/hdWGCNA_TOM/PAH/XeniumImpute/seurat_data_6'))
data.list[[7]] <- t(open_matrix_dir(dir = '~/Downloads/hdWGCNA_TOM/PAH/XeniumImpute/seurat_data_7'))

M2<-readRDS('~/Downloads/hdWGCNA_TOM/Xenium/collated_xenium_data.rds')
meta.data <-  M2@meta.data

M1 <- CreateSeuratObject(counts = data.list,meta.data = meta.data,assay = "Xenium")
M1[["BlankCodeword"]] <- M2[["BlankCodeword"]]
M1[["ControlCodeword"]] <- M2[["ControlCodeword"]]
M1[["ControlProbe"]] <- M2[["ControlProbe"]]
M1[['fov']] <- M2[['fov']]
M1[['fov.2']] <- M2[['fov.2']]
M1[['fov.3']] <- M2[['fov.3']]
M1[['fov.4']] <- M2[['fov.4']]
M1[['fov.5']] <- M2[['fov.5']]
M1[['fov.6']] <- M2[['fov.6']]
M1[['fov.7']] <- M2[['fov.7']]
M1[['fov.8']] <- M2[['fov.8']]
M1[['fov.9']] <- M2[['fov.9']]
M1[['seg1']] <- M2[['seg1']]
M1[['seg2']] <- M2[['seg2']]
M1[['seg3']] <- M2[['seg3']]
M1[['seg1.9']] <- M2[['seg1.9']]
M1[['seg2.8']] <- M2[['seg2.8']]

saveRDS(
  object = M1,
  file = "/Volumes/Extreme\ SSD/RV_Paper_Revisions/ImputedXenium.rds"
)

ImageDimPlot(M1, fov = "fov", molecules = c("NPPA"), nmols = 20000)



segmentations.data <- list(
  centroids = CreateCentroids(data$centroids),
  segmentation = CreateSegmentation(data$segmentations))
coords <- CreateFOV(
  coords = segmentations.data, 
  type = c("segmentation", "centroids"), 
  molecules = data$microns, 
  assay = "Xenium")


xenium.obj[["fov"]] <- coords









