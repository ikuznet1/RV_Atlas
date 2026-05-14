'''
Bash 



### Corrected, non_imputed

python spatialdata2seurat.py \
~/.xenium_explorer_cache/proseg_proseg-output-XETG00217__0038291__Region_1__20241212__142808_2dec91fb/spatialdata_proseg.zarr \
proseg_proseg-output-XETG00217__0038291__Region_1__20241212__142808_2dec91fb_corrected.rds \
--cache_dir ~/.xenium_explorer_cache/ \
--assay Xenium --fov 1467 --no_nucleus  --no_imputed --use_corrected 

python spatialdata2seurat.py \
~/.xenium_explorer_cache/proseg_proseg-output-XETG00217__0038291__Region_2__20241212__142808_2dec91fb/spatialdata_proseg.zarr \
proseg_proseg-output-XETG00217__0038291__Region_2__20241212__142808_2dec91fb_corrected.rds \
--cache_dir ~/.xenium_explorer_cache/ \
--assay Xenium --fov 1632 --no_nucleus  --no_imputed --use_corrected 

python spatialdata2seurat.py \
~/.xenium_explorer_cache/proseg_proseg-output-XETG00217__0038290__Region_1__20241212__142808_2dec91fb/spatialdata_proseg.zarr \
proseg_proseg-output-XETG00217__0038290__Region_1__20241212__142808_2dec91fb_corrected.rds \
--cache_dir ~/.xenium_explorer_cache/ \
--assay Xenium --fov 1561 --no_nucleus  --no_imputed --use_corrected 

python spatialdata2seurat.py \
~/.xenium_explorer_cache/proseg_proseg-output-XETG00217__0038290__Region_2__20241212__142808_2dec91fb/spatialdata_proseg.zarr \
proseg_proseg-output-XETG00217__0038290__Region_2__20241212__142808_2dec91fb_corrected.rds \
--cache_dir ~/.xenium_explorer_cache/ \
--assay Xenium --fov 1343 --no_nucleus  --no_imputed --use_corrected 

python spatialdata2seurat.py \
~/.xenium_explorer_cache/proseg_proseg-output-XETG00217__0038213__Region_1__20241206__182124_2dec91fb/spatialdata_proseg.zarr \
proseg_proseg-output-XETG00217__0038213__Region_1__20241206__182124_2dec91fb_top_corrected.rds \
--cache_dir ~/.xenium_explorer_cache/ \
--assay Xenium --roi ~/.xenium_explorer_cache/rois_d72b1263.json  \
--roi_cls Sample --roi_name Top --fov 1697 --no_nucleus  --no_imputed --use_corrected 

python spatialdata2seurat.py \
~/.xenium_explorer_cache/proseg_proseg-output-XETG00217__0038213__Region_1__20241206__182124_2dec91fb/spatialdata_proseg.zarr \
proseg_proseg-output-XETG00217__0038213__Region_1__20241206__182124_2dec91fb_middle_corrected.rds \
--cache_dir ~/.xenium_explorer_cache/ \
--assay Xenium --roi ~/.xenium_explorer_cache/rois_d72b1263.json  \
--roi_cls Sample --roi_name Middle --fov 1691 --no_nucleus  --no_imputed --use_corrected 

python spatialdata2seurat.py \
~/.xenium_explorer_cache/proseg_proseg-output-XETG00217__0038213__Region_1__20241206__182124_2dec91fb/spatialdata_proseg.zarr \
proseg_proseg-output-XETG00217__0038213__Region_1__20241206__182124_2dec91fb_bottom_corrected.rds \
--cache_dir ~/.xenium_explorer_cache/ \
--assay Xenium --roi ~/.xenium_explorer_cache/rois_d72b1263.json  \
--roi_cls Sample --roi_name Bottom --fov 1618 --no_nucleus  --no_imputed --use_corrected 

python spatialdata2seurat.py \
~/.xenium_explorer_cache/proseg_proseg-output-XETG00217__0038216__Region_1__20241206__182124_2dec91fb/spatialdata_proseg.zarr \
proseg_proseg-output-XETG00217__0038216__Region_1__20241206__182124_2dec91fb_top_corrected.rds \
--cache_dir ~/.xenium_explorer_cache/ \
--assay Xenium --roi ~/.xenium_explorer_cache/rois_16a362cb.json  \
--roi_cls Sample --roi_name Top --fov 1567 --no_nucleus  --no_imputed --use_corrected 

python spatialdata2seurat.py \
~/.xenium_explorer_cache/proseg_proseg-output-XETG00217__0038216__Region_1__20241206__182124_2dec91fb/spatialdata_proseg.zarr \
proseg_proseg-output-XETG00217__0038216__Region_1__20241206__182124_2dec91fb_bottom_corrected.rds \
--cache_dir ~/.xenium_explorer_cache/ \
--assay Xenium --roi ~/.xenium_explorer_cache/rois_16a362cb.json  \
--roi_cls Sample --roi_name Bottom --fov 1692 --no_nucleus  --no_imputed --use_corrected 


### Non-corrected, imputed

python spatialdata2seurat.py \
~/.xenium_explorer_cache/proseg_proseg-output-XETG00217__0038291__Region_1__20241212__142808_2dec91fb/spatialdata_proseg.zarr \
proseg_proseg-output-XETG00217__0038291__Region_1__20241212__142808_2dec91fb_imputed.rds \
--cache_dir ~/.xenium_explorer_cache/ \
--assay Xenium --fov 1467 --no_nucleus --merge_assays

python spatialdata2seurat.py \
~/.xenium_explorer_cache/proseg_proseg-output-XETG00217__0038291__Region_2__20241212__142808_2dec91fb/spatialdata_proseg.zarr \
proseg_proseg-output-XETG00217__0038291__Region_2__20241212__142808_2dec91fb_imputed.rds \
--cache_dir ~/.xenium_explorer_cache/ \
--assay Xenium --fov 1632 --no_nucleus --merge_assays

python spatialdata2seurat.py \
~/.xenium_explorer_cache/proseg_proseg-output-XETG00217__0038290__Region_1__20241212__142808_2dec91fb/spatialdata_proseg.zarr \
proseg_proseg-output-XETG00217__0038290__Region_1__20241212__142808_2dec91fb_imputed.rds \
--cache_dir ~/.xenium_explorer_cache/ \
--assay Xenium --fov 1561 --no_nucleus --merge_assays

python spatialdata2seurat.py \
~/.xenium_explorer_cache/proseg_proseg-output-XETG00217__0038290__Region_2__20241212__142808_2dec91fb/spatialdata_proseg.zarr \
proseg_proseg-output-XETG00217__0038290__Region_2__20241212__142808_2dec91fb_imputed.rds \
--cache_dir ~/.xenium_explorer_cache/ \
--assay Xenium --fov 1343 --no_nucleus --merge_assays

python spatialdata2seurat.py \
~/.xenium_explorer_cache/proseg_proseg-output-XETG00217__0038213__Region_1__20241206__182124_2dec91fb/spatialdata_proseg.zarr \
proseg_proseg-output-XETG00217__0038213__Region_1__20241206__182124_2dec91fb_top_imputed.rds \
--cache_dir ~/.xenium_explorer_cache/ \
--assay Xenium --roi ~/.xenium_explorer_cache/rois_d72b1263.json  \
--roi_cls Sample --roi_name Top --fov 1697 --no_nucleus --merge_assays

python spatialdata2seurat.py \
~/.xenium_explorer_cache/proseg_proseg-output-XETG00217__0038213__Region_1__20241206__182124_2dec91fb/spatialdata_proseg.zarr \
proseg_proseg-output-XETG00217__0038213__Region_1__20241206__182124_2dec91fb_middle_imputed.rds \
--cache_dir ~/.xenium_explorer_cache/ \
--assay Xenium --roi ~/.xenium_explorer_cache/rois_d72b1263.json  \
--roi_cls Sample --roi_name Middle --fov 1691 --no_nucleus --merge_assays

python spatialdata2seurat.py \
~/.xenium_explorer_cache/proseg_proseg-output-XETG00217__0038213__Region_1__20241206__182124_2dec91fb/spatialdata_proseg.zarr \
proseg_proseg-output-XETG00217__0038213__Region_1__20241206__182124_2dec91fb_bottom_imputed.rds \
--cache_dir ~/.xenium_explorer_cache/ \
--assay Xenium --roi ~/.xenium_explorer_cache/rois_d72b1263.json  \
--roi_cls Sample --roi_name Bottom --fov 1618 --no_nucleus --merge_assays

python spatialdata2seurat.py \
~/.xenium_explorer_cache/proseg_proseg-output-XETG00217__0038216__Region_1__20241206__182124_2dec91fb/spatialdata_proseg.zarr \
proseg_proseg-output-XETG00217__0038216__Region_1__20241206__182124_2dec91fb_top_imputed.rds \
--cache_dir ~/.xenium_explorer_cache/ \
--assay Xenium --roi ~/.xenium_explorer_cache/rois_16a362cb.json  \
--roi_cls Sample --roi_name Top --fov 1567 --no_nucleus --merge_assays

python spatialdata2seurat.py \
~/.xenium_explorer_cache/proseg_proseg-output-XETG00217__0038216__Region_1__20241206__182124_2dec91fb/spatialdata_proseg.zarr \
proseg_proseg-output-XETG00217__0038216__Region_1__20241206__182124_2dec91fb_bottom_imputed.rds \
--cache_dir ~/.xenium_explorer_cache/20241206_KristinaRawData/ \
--assay Xenium --roi ~/.xenium_explorer_cache/rois_16a362cb.json  \
--roi_cls Sample --roi_name Bottom --fov 1692 --no_nucleus --merge_assays


### Assign original Xenium nuclei


python transfer_seg_idents.py ~/.xenium_explorer_cache/proseg_proseg-output-XETG00217__0038213__Region_1__20241206__182124_2dec91fb/spatialdata_proseg.zarr \
    /Volumes/Extreme\ SSD/Xenium/20241206_KristinaRawData/output-XETG00217__0038213__Region_1__20241206__182124 \
    --identity_column cell_type_rctd_doublet \
    --query_boundary nucleus \
    --output_csv output-XETG00217__0038213__Region_1__20241206__182124_top_nuclei_assignments.csv \
    --roi ~/.xenium_explorer_cache/rois_d72b1263.json  --roi_cls Sample --roi_name Top

python transfer_seg_idents.py ~/.xenium_explorer_cache/proseg_proseg-output-XETG00217__0038213__Region_1__20241206__182124_2dec91fb/spatialdata_proseg.zarr \
    /Volumes/Extreme\ SSD/Xenium/20241206_KristinaRawData/output-XETG00217__0038213__Region_1__20241206__182124 \
    --identity_column cell_type_rctd_doublet \
    --query_boundary nucleus \
    --output_csv output-XETG00217__0038213__Region_1__20241206__182124_middle_nuclei_assignments.csv \
    --roi ~/.xenium_explorer_cache/rois_d72b1263.json  --roi_cls Sample --roi_name Middle

python transfer_seg_idents.py ~/.xenium_explorer_cache/proseg_proseg-output-XETG00217__0038213__Region_1__20241206__182124_2dec91fb/spatialdata_proseg.zarr \
    /Volumes/Extreme\ SSD/Xenium/20241206_KristinaRawData/output-XETG00217__0038213__Region_1__20241206__182124 \
    --identity_column cell_type_rctd_doublet \
    --query_boundary nucleus \
    --output_csv output-XETG00217__0038213__Region_1__20241206__182124_bottom_nuclei_assignments.csv \
    --roi ~/.xenium_explorer_cache/rois_d72b1263.json  --roi_cls Sample --roi_name Bottom
                                                                                                                 
python transfer_seg_idents.py ~/.xenium_explorer_cache/proseg_proseg-output-XETG00217__0038216__Region_1__20241206__182124_2dec91fb/spatialdata_proseg.zarr \
    /Volumes/Extreme\ SSD/Xenium/20241206_KristinaRawData/output-XETG00217__0038216__Region_1__20241206__182124 \
    --identity_column cell_type_rctd_doublet \
    --query_boundary nucleus \
    --output_csv output-XETG00217__0038216__Region_1__20241206__182124_top_nuclei_assignments.csv \
    --roi ~/.xenium_explorer_cache/rois_16a362cb.json  --roi_cls Sample --roi_name Top

python transfer_seg_idents.py ~/.xenium_explorer_cache/proseg_proseg-output-XETG00217__0038216__Region_1__20241206__182124_2dec91fb/spatialdata_proseg.zarr \
    /Volumes/Extreme\ SSD/Xenium/20241206_KristinaRawData/output-XETG00217__0038216__Region_1__20241206__182124 \
    --identity_column cell_type_rctd_doublet \
    --query_boundary nucleus \
    --output_csv output-XETG00217__0038216__Region_1__20241206__182124_bottom_nuclei_assignments.csv \
    --roi ~/.xenium_explorer_cache/rois_16a362cb.json  --roi_cls Sample --roi_name Bottom
                                                                                                                 
python transfer_seg_idents.py ~/.xenium_explorer_cache/proseg_proseg-output-XETG00217__0038290__Region_1__20241212__142808_2dec91fb/spatialdata_proseg.zarr \
    /Volumes/Extreme\ SSD/Xenium/20241212_KristinaRawData/output-XETG00217__0038290__Region_1__20241212__142808 \
    --identity_column cell_type_rctd_doublet \
    --query_boundary nucleus \
    --output_csv output-XETG00217__0038290__Region_1__20241212__142808_nuclei_assignments.csv
                                                                                                                 
python transfer_seg_idents.py ~/.xenium_explorer_cache/proseg_proseg-output-XETG00217__0038290__Region_2__20241212__142808_2dec91fb/spatialdata_proseg.zarr \
    /Volumes/Extreme\ SSD/Xenium/20241212_KristinaRawData/output-XETG00217__0038290__Region_2__20241212__142808 \
    --identity_column cell_type_rctd_doublet \
    --query_boundary nucleus \
    --output_csv output-XETG00217__0038290__Region_2__20241212__142808_nuclei_assignments.csv
                                                                                                                 
python transfer_seg_idents.py ~/.xenium_explorer_cache/proseg_proseg-output-XETG00217__0038291__Region_1__20241212__142808_2dec91fb/spatialdata_proseg.zarr \
    /Volumes/Extreme\ SSD/Xenium/20241212_KristinaRawData/output-XETG00217__0038291__Region_1__20241212__142808 \
    --identity_column cell_type_rctd_doublet \
    --query_boundary nucleus \
    --output_csv output-XETG00217__0038291__Region_1__20241212__142808_nuclei_assignments.csv
                                                                                                                 
python transfer_seg_idents.py ~/.xenium_explorer_cache/proseg_proseg-output-XETG00217__0038291__Region_2__20241212__142808_2dec91fb/spatialdata_proseg.zarr \
    /Volumes/Extreme\ SSD/Xenium/20241212_KristinaRawData/output-XETG00217__0038291__Region_2__20241212__142808 \
    --identity_column cell_type_rctd_doublet \
    --query_boundary nucleus \
    --output_csv output-XETG00217__0038291__Region_2__20241212__142808_nuclei_assignments.csv


python seurat2cellnest.py Xenium_resegmented_imputed_final.rds --output_dir cellnest_input/ --fov_col patient --cell_type_col cell_type_rctd_doublet

'''
options(future.globals.maxSize= 48012896000)


library(Seurat)
library(harmony)
library(ggplot2)
library(BPCells)
library(hdWGCNA)
library(dplyr)
#source('./dependencies/shared/spatial_functions.R')
source('~/Downloads/hdWGCNA_TOM/spatial_functions.R')


xenium.obj.1 <- readRDS('proseg_proseg-output-XETG00217__0038291__Region_1__20241212__142808_2dec91fb_corrected.rds')
xenium.obj.1$patient <- '1467'
xenium.obj.1$group <- 'RVF'
idx <- rownames(xenium.obj.1)[1:477]
xenium.obj.1 <- subset(xenium.obj.1,features=idx)

xenium.obj.2 <- readRDS('proseg_proseg-output-XETG00217__0038291__Region_2__20241212__142808_2dec91fb_corrected.rds')
xenium.obj.2$patient <- '1632'
xenium.obj.2$group <- 'RVF'
idx <- rownames(xenium.obj.2)[1:477]
xenium.obj.2 <- subset(xenium.obj.2,features=idx)

xenium.obj.3 <- readRDS('proseg_proseg-output-XETG00217__0038290__Region_1__20241212__142808_2dec91fb_corrected.rds')
xenium.obj.3$patient <- '1561'
xenium.obj.3$group <- 'NF'
idx <- rownames(xenium.obj.3)[1:477]
xenium.obj.3 <- subset(xenium.obj.3,features=idx)

xenium.obj.4 <- readRDS('proseg_proseg-output-XETG00217__0038290__Region_2__20241212__142808_2dec91fb_corrected.rds')
xenium.obj.4$patient <- '1343'
xenium.obj.4$group <- 'RVF'
idx <- rownames(xenium.obj.4)[1:477]
xenium.obj.4 <- subset(xenium.obj.4,features=idx)

xenium.obj.5 <- readRDS('proseg_proseg-output-XETG00217__0038213__Region_1__20241206__182124_2dec91fb_top_corrected.rds')
xenium.obj.5$patient <- '1697'
xenium.obj.5$group <- 'NF'
idx <- rownames(xenium.obj.5)[1:477]
xenium.obj.5 <- subset(xenium.obj.5,features=idx)

xenium.obj.6 <- readRDS('proseg_proseg-output-XETG00217__0038213__Region_1__20241206__182124_2dec91fb_middle_corrected.rds')
xenium.obj.6$patient <- '1691'
xenium.obj.6$group <- 'NF'
idx <- rownames(xenium.obj.6)[1:477]
xenium.obj.6 <- subset(xenium.obj.6,features=idx)

xenium.obj.7 <- readRDS('proseg_proseg-output-XETG00217__0038213__Region_1__20241206__182124_2dec91fb_bottom_corrected.rds')
xenium.obj.7$patient <- '1618'
xenium.obj.7$group <- 'pRV'
idx <- rownames(xenium.obj.7)[1:477]
xenium.obj.7 <- subset(xenium.obj.7,features=idx)

xenium.obj.8 <- readRDS('proseg_proseg-output-XETG00217__0038216__Region_1__20241206__182124_2dec91fb_top_corrected.rds')
xenium.obj.8$patient <- '1567'
xenium.obj.8$group <- 'pRV'
idx <- rownames(xenium.obj.8)[1:477]
xenium.obj.8 <- subset(xenium.obj.8,features=idx)

xenium.obj.9 <- readRDS('proseg_proseg-output-XETG00217__0038216__Region_1__20241206__182124_2dec91fb_bottom_corrected.rds')
xenium.obj.9$patient <- '1692'
xenium.obj.9$group <- 'pRV'
idx <- rownames(xenium.obj.9)[1:477]
xenium.obj.9 <- subset(xenium.obj.9,features=idx)




xenium.obj <- merge(xenium.obj.1,xenium.obj.2)
xenium.obj <- merge(xenium.obj,xenium.obj.3)
xenium.obj <- merge(xenium.obj,xenium.obj.4)
xenium.obj <- merge(xenium.obj,xenium.obj.5)
xenium.obj <- merge(xenium.obj,xenium.obj.6)
xenium.obj <- merge(xenium.obj,xenium.obj.7)
xenium.obj <- merge(xenium.obj,xenium.obj.8)
xenium.obj <- merge(xenium.obj,xenium.obj.9)

xenium.obj <- NormalizeData(xenium.obj)
xenium.obj <- FindVariableFeatures(xenium.obj)
xenium.obj <- ScaleData(xenium.obj,split.by='patient')
xenium.obj <- RunPCA(xenium.obj, npcs = 30)
xenium.obj <- JoinLayers(xenium.obj)
xenium.obj <- RunHarmony(xenium.obj, group.by.vars = 'patient',reduction.save = "harmony",kmeans_init_nstart=20, kmeans_init_iter_max=100)

xenium.obj <- FindNeighbors(xenium.obj, dims = 1:30, reduction = "harmony")
xenium.obj <- FindClusters(xenium.obj, resolution = 1)
xenium.obj <- RunUMAP(xenium.obj, reduction = "pca", dims = 1:30,reduction.name='umap_orig')
xenium.obj <- RunUMAP(xenium.obj, reduction = "harmony", dims = 1:30,reduction.name='umap')

saveRDS(xenium.obj,'./Xenium_resegmented_corrected.rds')
#n = 568651

rm(xenium.obj.1)
rm(xenium.obj.2)
rm(xenium.obj.3)
rm(xenium.obj.4)
rm(xenium.obj.5)
rm(xenium.obj.6)
rm(xenium.obj.7)
rm(xenium.obj.8)
rm(xenium.obj.9)
gc()

############ UMAP (Panel A)

pdf(paste0('./output/Xenium/', 'Xenium_snUMAP_new.pdf'), width=5, height=5)
PlotEmbedding(xenium.obj,group.by='cell_type_rctd_doublet',reduction='umap_orig',point_size=.05,plot_under=TRUE,plot_theme=umap_theme()+NoLegend(),raster_dpi=400,raster_scale=0.5)
dev.off()

pdf(paste0('./output/Xenium/', 'Xenium_snUMAP_new_pt.pdf'), width=5, height=5)
PlotEmbedding(xenium.obj,group.by='patient',reduction='umap_orig',point_size=.05,plot_under=TRUE,plot_theme=umap_theme()+NoLegend(),raster_dpi=400,raster_scale=0.5)
dev.off()

############ Cell Markers 

Idents(xenium.obj) <- 'cell_type_rctd_doublet'
cm.marks <- FindMarkers(xenium.obj,ident.1 = 'CM',only.pos=T)

fb.marks <- FindMarkers(xenium.obj,ident.1 = 'FB',only.pos=T)

sm.marks <- FindMarkers(xenium.obj,ident.1 = 'SM',only.pos=T)

pc.marks <- FindMarkers(xenium.obj,ident.1 = 'PC',only.pos=T)

ec.marks <- FindMarkers(xenium.obj,ident.1 = 'EC',only.pos=T)

endo.marks <- FindMarkers(xenium.obj,ident.1 = 'Endo',only.pos=T)

lec.marks <- FindMarkers(xenium.obj,ident.1 = 'LEC',only.pos=T)

adipo.marks <- FindMarkers(xenium.obj,ident.1 =  'Adipo',only.pos=T)

epi.marks <- FindMarkers(xenium.obj,ident.1 =  'Epi',only.pos=T)

#CD28
nkt.marks <- FindMarkers(xenium.obj,ident.1 = 'NKT',only.pos=T)

myeloid.marks <- FindMarkers(xenium.obj,ident.1 = 'Myeloid',only.pos=T)

neuron.marks <- FindMarkers(xenium.obj,ident.1 = 'Neuron',only.pos=T)



########### Cell frequency

cell_freq <- table(xenium.obj$cell_type_rctd_doublet,xenium.obj$patient)
cell_freq <- t(t(cell_freq)/colSums(cell_freq))
idx <- rep(c('NF','pRV','RVF'),9)[t(table(xenium.obj$patient,xenium.obj$group) > 0)]
cell_freq <- data.frame(t(cell_freq))
cell_freq$group <- rep(idx,12)

ggplot(cell_freq, aes(y=Freq, fill=group)) +
  geom_boxplot() +
  facet_wrap(~Var2)

cell_freq <- table(xenium.obj$cell_type_seurat,xenium.obj$patient)
cell_freq <- t(t(cell_freq)/colSums(cell_freq))
idx <- rep(c('NF','pRV','RVF'),9)[t(table(xenium.obj$patient,xenium.obj$group) > 0)]
cell_freq <- data.frame(t(cell_freq))
cell_freq$group <- rep(idx,34)

ggplot(cell_freq, aes(y=Freq, fill=group)) +
  geom_boxplot() +
  facet_wrap(~Var2)



########### Cell Density

df <- data.frame(celltype = xenium.obj$cell_type_rctd_doublet, cellsubtype = xenium.obj$cell_type_seurat,
    group = xenium.obj$group, patient = xenium.obj$patient, volume = xenium.obj$cell_area)
df$radius <- (df$volume *3/4/pi)^(1/3)
df$area <- pi*df$radius^2

# Normalize out systematic segmentation differences between samples
# Only CMs and maybe FBs should change in size
df.t <- df %>% group_by(patient,celltype) %>%
  mutate(mean_volume = mean(volume, na.rm = TRUE)) %>% 
  filter(!(celltype %in% c('CM','FB'))) %>% 
  group_by(patient) %>%
  summarize(mean_mean_volume = mean(mean_volume, na.rm = TRUE))

df$mean_mean_volume = 0
for (i in 1:dim(df.t)[1]){
  pt <- as.numeric(df.t[i,1][[1]])
  df[df$patient == (pt),]$mean_mean_volume = df.t[i,2][[1]]

}


ggplot(df, aes(y=mean_mean_volume, fill=patient)) +
  geom_boxplot() +
  facet_wrap(~celltype)


# Average cell volume
df.sum <- df %>% group_by(patient,celltype,group) %>%
  summarise(mean_volume = sd(volume), .groups = 'drop')


ggplot(df.sum, aes(y=mean_volume, fill=group)) +
  geom_boxplot() +
  facet_wrap(~celltype)

# STD cell volume
df.sum <- df %>% group_by(patient,celltype,group) %>%
  summarise(std_volume = sd(volume), .groups = 'drop')

ggplot(df.sum, aes(y=std_volume, fill=group)) +
  geom_boxplot() +
  facet_wrap(~celltype)


#Average cell volume (normalized by segmentation differences)
df.sum <- df %>% group_by(patient,celltype,group) %>%
  summarise(norm_volume = mean(volume)/mean(mean_mean_volume), .groups = 'drop')

ggplot(df.sum, aes(y=norm_volume, fill=group)) +
  geom_boxplot() +
  facet_wrap(~celltype)

# Total volume by celltype

df.sum <- df %>% group_by(patient,celltype,group) %>%
  summarise(total_volume = sum(volume), .groups = 'drop') %>%
  group_by(patient)  %>%
  mutate(volume_normalized = total_volume / sum(total_volume, na.rm = TRUE))


ggplot(df.sum, aes(y=volume_normalized, fill=group)) +
  geom_boxplot() +
  facet_wrap(~celltype)


# Total normalized (by segmentation differences and total volume) volume by celltype

df.sum <- df %>% group_by(patient,celltype,group) %>%
  summarise(total_volume = sum(volume),total_norm_volume = sum(volume)/mean(mean_mean_volume), .groups = 'drop') %>%
  group_by(patient)  %>%
  mutate(volume_normalized = total_norm_volume / sum(total_volume, na.rm = TRUE))


ggplot(df.sum, aes(y=volume_normalized, fill=group)) +
  geom_boxplot() +
  facet_wrap(~celltype)











# Cell count normalized by area

df.sum <- df %>% group_by(patient,celltype,group,mean_area) %>%
  count(celltype) %>%
  group_by(patient)  %>%
  mutate(n_norm = n / sum(n, na.rm = TRUE)) %>%
  group_by(patient,celltype)  %>%
  mutate(n_norm_norm = n_norm / mean_area)

df.sum <- df %>% group_by(patient,celltype,group,mean_area) %>%
  count(celltype) %>%
  group_by(patient,celltype) %>%
  mutate(count_corrected = n / mean_area) %>%
  group_by(patient) %>%
  mutate(n_norm = count_corrected / sum(count_corrected, na.rm = TRUE)) 

df.sum <- df %>% group_by(patient,celltype,group,mean_area) %>%
  count(celltype) %>%
  group_by(patient,celltype) %>%
  mutate(area_counts = n * mean_area) %>%
  group_by(patient) %>%
  mutate(n_norm = area_counts / sum(area_counts, na.rm = TRUE)) 


df.sum <- df %>% group_by(patient) %>%
  mutate(total_area = sum(area)) %>% 
  ungroup() %>% group_by(patient,group,total_area) %>%
  count(celltype) %>%
  mutate(n_norm = n / total_area, na.rm = TRUE)


ggplot(df.sum, aes(y=n_norm, fill=group)) +
  geom_boxplot() +
  facet_wrap(~celltype)







df.sum <- df %>% group_by(patient,cellsubtype,celltype,group) %>%
  summarise(total_area = sum(area), .groups = 'drop') %>%
  group_by(patient)  %>%
  mutate(area_normalized = total_area / sum(total_area, na.rm = TRUE)) %>%
  filter(cellsubtype %in% c('EC_Arterial','EC_Capillary','EC_Venous')) %>%
  filter(celltype %in% c('EC'))



ggplot(df.sum, aes(y=area_normalized, fill=group)) +
  geom_boxplot() +
  facet_wrap(~cellsubtype)


df.sum <- df %>% group_by(patient,cellsubtype,group) %>%
  summarise(total_area = sum(area), .groups = 'drop') %>%
  group_by(patient)  %>%
  mutate(area_normalized = total_area / sum(total_area, na.rm = TRUE)) %>%
  filter(cellsubtype %in% c('Mono','NK_T','TREM2 Mac','CCR2- rMac2','CCR2- rMac1','iMac','CCR2+ rMac'))



ggplot(df.sum, aes(y=area_normalized, fill=group)) +
  geom_boxplot() +
  facet_wrap(~cellsubtype)

# CM hypertrophy? Vs nuclei

# Should I look at nuclei counts?





#### Imputed data — on-disk merge via BPCells
# Each imputed RDS is ~20 GB; loading all 9 simultaneously requires ~180 GB RAM.
# BPCells converts sparse matrices to memory-mapped on-disk format so merge/join
# operate without holding all data in memory.

library(BPCells)

imputed_files <- c(
  'proseg_proseg-output-XETG00217__0038291__Region_1__20241212__142808_2dec91fb_imputed.rds',
  'proseg_proseg-output-XETG00217__0038291__Region_2__20241212__142808_2dec91fb_imputed.rds',
  'proseg_proseg-output-XETG00217__0038290__Region_1__20241212__142808_2dec91fb_imputed.rds',
  'proseg_proseg-output-XETG00217__0038290__Region_2__20241212__142808_2dec91fb_imputed.rds',
  'proseg_proseg-output-XETG00217__0038213__Region_1__20241206__182124_2dec91fb_top_imputed.rds',
  'proseg_proseg-output-XETG00217__0038213__Region_1__20241206__182124_2dec91fb_middle_imputed.rds',
  'proseg_proseg-output-XETG00217__0038213__Region_1__20241206__182124_2dec91fb_bottom_imputed.rds',
  'proseg_proseg-output-XETG00217__0038216__Region_1__20241206__182124_2dec91fb_top_imputed.rds',
  'proseg_proseg-output-XETG00217__0038216__Region_1__20241206__182124_2dec91fb_bottom_imputed.rds'
)
patients <- c('1467','1632','1561','1343','1697','1691','1618','1567','1692')
groups   <- c('RVF','RVF','NF','RVF','NF','NF','pRV','pRV','pRV')

bp_dir <- './bpcells_imputed'
dir.create(bp_dir, showWarnings = FALSE, recursive = TRUE)

# Pass 1: discover shared genes across all 9 files (reads only rownames, not matrix)
message("=== Pass 1: discovering shared gene set ===")
gene_sets <- list()
for (i in seq_along(imputed_files)) {
  message("  Reading gene names from sample ", i, "/", length(imputed_files))
  obj <- readRDS(imputed_files[i])
  gene_sets[[i]] <- rownames(obj)
  rm(obj); gc(verbose = FALSE)
}
shared_genes <- Reduce(intersect, gene_sets)
message("Shared genes: ", length(shared_genes))
rm(gene_sets); gc(verbose = FALSE)

# Pass 2: convert each RDS to on-disk BPCells, keep lightweight Seurat shell
message("=== Pass 2: converting to on-disk BPCells format ===")
obj_list <- vector("list", length(imputed_files))
for (i in seq_along(imputed_files)) {
  message("Processing sample ", i, "/", length(imputed_files), ": ", patients[i])
  obj <- readRDS(imputed_files[i])
  obj$patient <- patients[i]
  obj$group   <- groups[i]

  # Subset to shared genes
  obj <- subset(obj, features = shared_genes)

  # Write counts to disk and replace with on-disk pointer
  assay_name <- DefaultAssay(obj)
  counts_path <- file.path(bp_dir, paste0("sample_", i, "_counts"))
  mat <- GetAssayData(obj, assay = assay_name, layer = "counts")
  BPCells::write_matrix_dir(mat, dir = counts_path, overwrite = TRUE)
  obj[[assay_name]]$counts <- BPCells::open_matrix_dir(counts_path)

  obj_list[[i]] <- obj
  rm(obj); gc(verbose = FALSE)
  message("  Done — on-disk at ", bp_dir, "/sample_", i, "_*")
}

# Merge (lightweight: matrices are on-disk pointers, only metadata is in memory)
message("=== Merging ===")
xenium.imp <- merge(obj_list[[1]], obj_list[2:length(obj_list)])
rm(obj_list); gc(verbose = FALSE)



# Standard Seurat v5 workflow (all on-disk via BPCells)
xenium.imp <- NormalizeData(xenium.imp)
xenium.imp <- FindVariableFeatures(xenium.imp)
xenium.imp <- ScaleData(xenium.imp, split.by = 'patient')
xenium.imp <- RunPCA(xenium.imp, npcs = 30)
xenium.imp <- JoinLayers(xenium.imp)
xenium.imp <- RunHarmony(xenium.imp, group.by.vars = 'patient',
                         reduction.save = "harmony")
xenium.imp <- FindNeighbors(xenium.imp, dims = 1:30, reduction = "harmony")
xenium.imp <- FindClusters(xenium.imp, resolution = 1)
xenium.imp <- RunUMAP(xenium.imp, reduction = "pca", dims = 1:30, reduction.name = 'umap_orig')
xenium.imp <- RunUMAP(xenium.imp, reduction = "harmony", dims = 1:30, reduction.name = 'umap')

saveRDS(xenium.imp, './Xenium_resegmented_imputed.rds')
message("=== Done — saved Xenium_resegmented_imputed.rds ===")


# Remove cells with <10 transcripts
xenium.imp <- xenium.imp[,colSums(GetAssayData(xenium.imp, assay = assay_name, layer = "counts")[rownames(xenium.obj),]) > 10]
xenium.imp <- xenium.imp[,!(xenium.imp$cell_type_rctd_doublet %in% c('None','Unknown'))]
saveRDS(xenium.imp, './Xenium_resegmented_imputed_final.rds')
# n = 625305

DimPlot(xenium.imp,group.by='cell_type_rctd_doublet',label=T)

xenium.obj <- readRDS('./Xenium_resegmented_corrected.rds')


options(future.globals.maxSize = Inf)  # allow arbitrarily large globals (BPCells + large objects)
plan("sequential")
anchors <- FindTransferAnchors(reference = xenium.obj, query = xenium.imp, dims = 1:30, reference.reduction = "harmony")

xenium.obj <- RunUMAP(xenium.obj, reduction = "pca", dims = 1:30,reduction.name='umap',return.model = TRUE)


xenium.imp <- MapQuery(anchorset = anchors, reference = xenium.obj, query = xenium.imp,
    refdata = list(cell_type_rctd_doublet_pred = "cell_type_rctd_doublet"), reference.reduction = "harmony", reduction.model = "umap")



############ Cell Markers (Panel B)

xenium.obj <- readRDS('~/Documents/XeniumWorkflow/functions/Xenium_resegmented_imputed_final.rds')
xenium.base <- readRDS('./Xenium_resegmented_corrected.rds')
non_imputed_genes <- intersect(rownames(xenium.obj),rownames(xenium.base))
rm(xenium.base)
gc()

Idents(xenium.obj) <- 'cell_type_rctd_doublet'
levels(xenium.obj) <- sort(levels(xenium.obj))                                                                   


cm.marks <- FindMarkers(xenium.obj,ident.1 = 'CM',only.pos=T, features = non_imputed_genes)

fb.marks <- FindMarkers(xenium.obj,ident.1 = 'FB',only.pos=T, features = non_imputed_genes)

sm.marks <- FindMarkers(xenium.obj,ident.1 = 'SM',only.pos=T, features = non_imputed_genes)

pc.marks <- FindMarkers(xenium.obj,ident.1 = 'PC',only.pos=T, features = non_imputed_genes)

ec.marks <- FindMarkers(xenium.obj,ident.1 = 'EC',only.pos=T, features = non_imputed_genes)

endo.marks <- FindMarkers(xenium.obj,ident.1 = 'Endo',only.pos=T, features = non_imputed_genes)

lec.marks <- FindMarkers(xenium.obj,ident.1 = 'LEC',only.pos=T, features = non_imputed_genes)

adipo.marks <- FindMarkers(xenium.obj,ident.1 =  'Adipo',only.pos=T, features = non_imputed_genes)

epi.marks <- FindMarkers(xenium.obj,ident.1 =  'Epi',only.pos=T, features = non_imputed_genes)

nkt.marks <- FindMarkers(xenium.obj,ident.1 = 'NKT',only.pos=T, features = non_imputed_genes)

myeloid.marks <- FindMarkers(xenium.obj,ident.1 = 'Myeloid',only.pos=T, features = non_imputed_genes)

neuron.marks <- FindMarkers(xenium.obj,ident.1 = 'Neuron',only.pos=T, features = non_imputed_genes)

pdf(paste0('./output/Xenium/', 'Xenium_marks.pdf'), width=8, height=5)
DotPlot(xenium.obj, features = c("ADIPOQ","ACTC1","VWF",
  "PECAM1","BMP4","WT1","DCN","PROX1","CD68","F13A1","SCN7A",
  "CD3D","CD69","PDGFRB","ACTA2","MYH11"),col.min = 0, col.max = 2)
dev.off()



##############################################                                                                   
#### SPATIAL DENSITY (cells per µm²)
##############################################                                                                   
                                                                                                                 
# Compute tissue area from cell segmentation polygons (concave hull per FOV)                                      
library(concaveman)
library(sf)

# Build convex hulls + centroids for all FOVs                                                                    
hull_list <- list()                       
fov_names <- Images(xenium.obj)

for (fov_name in fov_names) {
  coords <- GetTissueCoordinates(xenium.obj[[fov_name]][["centroids"]])
  if (nrow(coords) < 3) next                                                                                     
   
  hull_pts <- concaveman(as.matrix(coords[, c("x", "y")]), concavity = 0.7)                                        
  poly <- st_polygon(list(hull_pts))
  hull_list[[fov_name]] <- st_sf(fov = fov_name, geometry = st_sfc(poly))                                        
}                                                                                                                
   
hulls_sf <- do.call(rbind, hull_list)    
      
fov_name <- fov_names[5]  # pick your FOV

coords <- GetTissueCoordinates(xenium.obj[[fov_name]][["centroids"]])                                            
hull_single <- hulls_sf[hulls_sf$fov == fov_name, ]
                                                                                                                   
ggplot() +      
    geom_point(data = coords, aes(x = x, y = -y), size = 0.1, alpha = 0.2) +                                       
    geom_sf(data = hull_single, fill = NA, color = "red", linewidth = 0.8) +                                       
    theme_minimal() +
    labs(title = fov_name)




                                                                                                                 
compute_tissue_area <- function(obj, fov) {                                                                      
    coords <- GetTissueCoordinates(obj[[fov]][["centroids"]])
    if (nrow(coords) < 3) return(0)                                                                                
    hull_pts <- concaveman(as.matrix(coords[, c("x", "y")]), concavity = 0.7)                                      
    poly <- st_polygon(list(hull_pts))
    return(as.numeric(st_area(poly)))  # µm²                                                                       
} 

fov_names <- Images(xenium.obj)
tissue_area_um2 <- list()
for (fov_name in fov_names) {
  tissue_area_um2[fov_name] <- compute_tissue_area(xenium.obj,fov_name)                                                               
}

tissue_area_mm2 <- unlist(tissue_area_um2) / 1e6
tissue_area_um2 <- unlist(tissue_area_um2) 


# Cell type density
ct_counts <- table(xenium.obj$patient,xenium.obj$cell_type_rctd_doublet)                                                            
ct_density_mm2 <- as.data.frame(ct_counts / tissue_area_mm2)  # cells per mm²
ct_density_mm2['group'] <- (unique(xenium.obj@meta.data[, c("patient", "group")]) %>% arrange(patient))$group



ggplot(ct_density_mm2, aes(x = group, y = Freq, fill = group)) +                                              
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 1.5) +                                                                         
    facet_wrap(~Var2, scales = "free_y") +
    labs(y = expression("Cell density (cells/mm"^2*")"), x = NULL) +                                               
    theme_bw() +                                                                                                   
    theme(legend.position = "none",                                                                                
          strip.text = element_text(face = "bold")) +                                                              
    scale_x_discrete(limits = c("NF", "pRV", "RVF"))                                                               
   



# Neighborhood enrichment

library(RANN)  # for nn2() fast nearest neighbors                                                                
library(ggplot2)                                                                                                 
library(reshape2)                                                                                                
                                                                                                                 
compute_spatial_proximity <- function(obj, label_col = "cell_type_rctd_doublet",                                 
                                       n_neighbors = 20) {
  # Get coordinates and labels                                                                                   
  coords <- cbind(obj$x_centroid, obj$y_centroid)                                                                
  labels <- obj@meta.data[[label_col]]
                                                                                                                 
  # Drop NAs                                                                                                     
  valid <- !is.na(labels)
  coords <- coords[valid, ]                                                                                      
  labels <- labels[valid]

  cell_types <- sort(unique(labels))                                                                             
  n_types <- length(cell_types)
  label_idx <- match(labels, cell_types)  # integer labels 1..n_types                                            
                                                                                                                 
  # Build kNN graph
  message(sprintf("Computing %d nearest neighbors for %d cells...", n_neighbors, nrow(coords)))                  
  nn <- nn2(coords, k = n_neighbors + 1)  # +1 because first neighbor is self                                    
  neighbor_indices <- nn$nn.idx[, -1]      # drop self                                                           
                                                                                                                 
  # Build raw count matrix (how often type A neighbors type B)                                                   
  message("Building neighbor type matrix...")                                                                    
  raw_matrix <- matrix(0, nrow = n_types, ncol = n_types)                                                        
  src <- rep(label_idx, each = n_neighbors)                                                                      
  tgt <- label_idx[as.vector(t(neighbor_indices))]
  for (i in seq_along(src)) {                                                                                    
    raw_matrix[src[i], tgt[i]] <- raw_matrix[src[i], tgt[i]] + 1                                                 
  }                                                                                                              
  rownames(raw_matrix) <- colnames(raw_matrix) <- cell_types                                                     
                                                                                                                 
  # Cell type frequencies
  counts <- table(factor(labels, levels = cell_types))
  freqs <- as.numeric(counts) / sum(counts)                                                                      

  # Co-occurrence: log2(observed / expected)                                                                     
  total_pairs <- sum(raw_matrix)
  expected <- outer(freqs, freqs) * total_pairs * 2                                                              
  cooccurrence <- log2(raw_matrix / (expected + 1e-9) + 1e-9)
                                                                                                                 
  # Neighborhood enrichment: normalize rows by source count, cols by target freq                                 
  enrichment <- raw_matrix                                                                                       
  safe_counts <- ifelse(as.numeric(counts) > 0, as.numeric(counts), 1)                                           
  enrichment <- enrichment / safe_counts  # normalize by source cell count                                       
  safe_freqs <- ifelse(freqs > 0, freqs, 1e-9)                                                                   
  enrichment <- log2(t(t(enrichment) / safe_freqs) + 1e-9)                                                       
                                                                                                                 
  rownames(cooccurrence) <- colnames(cooccurrence) <- cell_types
  rownames(enrichment) <- colnames(enrichment) <- cell_types                                                     
                                                                                                                 
  return(list(
    cooccurrence = cooccurrence,                                                                                 
    enrichment = enrichment,
    raw_counts = raw_matrix,
    cell_type_counts = counts
  ))                                                                                                             
}
                                                                                                                 
# Run it        
result <- compute_spatial_proximity(xenium.obj, label_col = 'cell_type_seurat', n_neighbors = 30)
  
pdf(paste0('./output/Xenium/', 'Xenium_spatial_prox.pdf'), width=8, height=7)                                                                                                               
# Plot co-occurrence heatmap
cooc_df <- melt(result$cooccurrence)                                                                             
ggplot(cooc_df, aes(x = Var2, y = Var1, fill = value)) +                                                         
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,                                  
                       name = "log2(O/E)") +                                                                     
  theme_minimal() +                                                                                              
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +                                                     
  labs(title = "Co-occurrence (log2 observed/expected)", x = NULL, y = NULL) + 
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,                                    
                       limits = c(-2, 2), oob = scales::squish,                                                    
                       name = "log2(O/E)") 
dev.off()                                     
                                                                                                                 
# Plot neighborhood enrichment heatmap
enrich_df <- melt(result$enrichment)                                                                             
ggplot(enrich_df, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +                                                                                                  
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                       name = "log2(O/E)") +                                                                     
  theme_minimal() +                                                                                              
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Neighborhood Enrichment (log2 observed/expected)", x = NULL, y = NULL)  




##############################################
#### 1. CHORD DIAGRAM
##############################################
library(circlize)

# Use enrichment (log2 O/E) — normalizes for cell type frequency
# Symmetrize and shift to non-negative (chordDiagram requires non-negative values)
enrich_sym <- (result$enrichment + t(result$enrichment)) / 2
enrich_pos <- log2(enrich_sym)
enrich_pos[is.na(enrich_pos)] = 0

# Threshold weak interactions
enrich_pos[enrich_pos < quantile(enrich_pos, 0.9)] <- 0

pdf('./output/Xenium/spatial_chord.pdf', width = 8, height = 8)
chordDiagram(enrich_pos, symmetric = TRUE,
             annotationTrack = "grid",
             preAllocateTracks = list(track.height = 0.05))
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 1.4)
}, bg.border = NA)
dev.off()


##############################################
#### 4. NICHE COMPOSITION (microenvironment profile)
##############################################

# Normalize raw neighbor fractions by global cell type frequency (O/E)
# so common cell types don't dominate just by abundance
global_freq <- as.numeric(table(factor(xenium.obj$cell_type_rctd_doublet,
                                        levels = rownames(result$raw_counts)))) /
               sum(!is.na(xenium.obj$cell_type_rctd_doublet))

niche_profile_raw <- result$raw_counts / rowSums(result$raw_counts)
niche_profile_oe  <- t(t(niche_profile_raw) / global_freq)  # divide each col by global freq
niche_profile_oe[is.nan(niche_profile_oe) | is.infinite(niche_profile_oe)] <- 0
# Renormalize rows to sum to 1 for stacked bar readability
niche_profile_oe  <- niche_profile_oe / rowSums(niche_profile_oe)

niche_df <- melt(niche_profile_oe)
colnames(niche_df) <- c("source", "neighbor", "fraction_oe")

pdf('./output/Xenium/niche_composition.pdf', width = 25, height = 10)
ggplot(niche_df, aes(x = neighbor, y = fraction_oe, fill = neighbor)) +
  geom_col() +
  facet_wrap(~source) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        strip.text = element_text(face = "bold")) +
  labs(x = "Neighbor cell type", y = "Frequency-normalized neighbor fraction (O/E)",
       title = "Niche composition per cell type (normalized for abundance)")
dev.off()


##############################################
#### 5. PER-PATIENT INTERACTION SCORES (boxplots by group)
##############################################

# Run compute_spatial_proximity per patient, extract all pairwise enrichment scores
patient_meta <- unique(xenium.obj@meta.data[, c("patient", "group")]) %>% arrange(patient)

interaction_list <- list()
for (pat in patient_meta$patient) {
  message(sprintf("Computing spatial proximity for patient %s...", pat))
  obj_sub <- subset(xenium.obj, subset = patient == pat)
  res_pat <- compute_spatial_proximity(obj_sub, n_neighbors = 20)
  enrich_mat <- res_pat$enrichment
  df <- melt(enrich_mat)
  colnames(df) <- c("source", "target", "enrichment")
  df$patient <- pat
  df$group <- patient_meta$group[patient_meta$patient == pat]
  df$pair <- paste0(df$source, " -> ", df$target)
  interaction_list[[pat]] <- df
}

interaction_df <- do.call(rbind, interaction_list)
interaction_df$group <- factor(interaction_df$group, levels = c("NF", "pRV", "RVF"))

pdf('./output/Xenium/interaction_boxplots.pdf', width = 18, height = 14)
ggplot(interaction_df, aes(x = group, y = enrichment, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1) +
  facet_wrap(~pair, scales = "free_y") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  theme_bw() +
  theme(legend.position = "none",
        strip.text = element_text(size = 6),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Neighborhood enrichment (log2 O/E)", x = NULL,
       title = "Per-patient spatial interaction scores by group")
dev.off()


##############################################
#### 6. NICHE CLUSTERING (unsupervised, per-cell neighbor composition)
##############################################
library(RANN)
library(uwot)

# Build per-cell neighbor type composition vector
coords  <- cbind(xenium.obj$x_centroid, xenium.obj$y_centroid)
labels  <- xenium.obj$cell_type_rctd_doublet
valid   <- !is.na(labels)
coords  <- coords[valid, ]
labels  <- labels[valid]
cell_ids <- colnames(xenium.obj)[valid]

cell_types <- sort(unique(labels))
n_types    <- length(cell_types)
label_idx  <- match(labels, cell_types)

message("Computing kNN for niche clustering...")
nn <- nn2(coords, k = 21)
neighbor_indices <- nn$nn.idx[, -1]  # drop self

# For each cell, fraction of neighbors that are each type
message("Building per-cell niche vectors...")
niche_mat <- matrix(0, nrow = length(labels), ncol = n_types)
colnames(niche_mat) <- cell_types
for (i in seq_len(nrow(neighbor_indices))) {
  neighbor_types <- label_idx[neighbor_indices[i, ]]
  tab <- tabulate(neighbor_types, nbins = n_types)
  niche_mat[i, ] <- tab / sum(tab)
}

# Normalize each column by global cell type frequency (O/E)
# so abundant types don't dominate clustering
global_freq_cells <- as.numeric(table(factor(labels, levels = cell_types))) / length(labels)
niche_mat_oe <- t(t(niche_mat) / global_freq_cells)
# Renormalize rows to sum to 1
niche_mat_oe <- niche_mat_oe / rowSums(niche_mat_oe)

# Cluster niche vectors with kmeans
set.seed(42)
n_niches <- 10
message(sprintf("Clustering into %d niches...", n_niches))
km <- kmeans(niche_mat_oe, centers = n_niches, nstart = 25, iter.max = 100)
niche_labels <- km$cluster

# UMAP of niche vectors
message("Running UMAP on niche vectors...")
niche_umap <- umap(niche_mat_oe, n_neighbors = 15, min_dist = 0.3, n_threads = 4)

niche_df_plot <- data.frame(
  UMAP1      = niche_umap[, 1],
  UMAP2      = niche_umap[, 2],
  niche      = factor(niche_labels),
  cell_type  = labels,
  patient    = xenium.obj$patient[valid],
  group      = xenium.obj$group[valid]
)

pdf('./output/Xenium/niche_umap.pdf', width = 12, height = 5)

# UMAP colored by niche cluster
p1 <- ggplot(niche_df_plot, aes(x = UMAP1, y = UMAP2, color = niche)) +
  geom_point(size = 0.3, alpha = 0.4) +
  theme_void() +
  labs(title = "Niche clusters") +
  guides(color = guide_legend(override.aes = list(size = 3)))

# UMAP colored by cell type
p2 <- ggplot(niche_df_plot, aes(x = UMAP1, y = UMAP2, color = cell_type)) +
  geom_point(size = 0.3, alpha = 0.4) +
  theme_void() +
  labs(title = "Cell type") +
  guides(color = guide_legend(override.aes = list(size = 3)))

print(p1 + p2)
dev.off()

# Niche composition per cluster (what cell types populate each niche)
niche_comp <- aggregate(niche_mat_oe, by = list(niche = niche_labels), FUN = mean)
niche_comp_df <- melt(niche_comp, id.vars = "niche")
colnames(niche_comp_df) <- c("niche", "cell_type", "fraction")
niche_comp_df$niche <- factor(niche_comp_df$niche)

pdf('./output/Xenium/niche_composition_clusters.pdf', width = 12, height = 6)
ggplot(niche_comp_df, aes(x = niche, y = fraction, fill = cell_type)) +
  geom_col(position = "stack") +
  theme_bw() +
  labs(x = "Niche cluster", y = "Mean neighbor fraction",
       title = "Cell type composition per niche cluster")
dev.off()

# Store niche labels back onto the Seurat object
xenium.obj$spatial_niche <- NA
xenium.obj$spatial_niche[valid] <- niche_labels


##############################################
#### 7. Area normalization
##############################################

# Compute 2D polygon area per cell from segmentation boundaries
# (Proseg reports volumes not areas, so we recompute from 2D polygons)
library(sf)

get_cell_areas_2d <- function(obj, fov_name) {
  print(fov_name)
  seg <- obj[[fov_name]][["segmentation"]]
  coords_list <- GetTissueCoordinates(seg, which = "segmentation")
  areas <- tapply(seq_len(nrow(coords_list)), coords_list$cell, function(idx) {
    pts <- coords_list[idx, c("x", "y")]
    if (nrow(pts) < 3) return(NA)
    poly <- st_polygon(list(as.matrix(rbind(pts, pts[1, ]))))
    st_area(poly)
  })
  data.frame(cell = names(areas), cell_area_2d = as.numeric(areas))
}

area_list <- lapply(Images(xenium.obj), function(fov) get_cell_areas_2d(xenium.obj, fov))
cell_areas_2d <- do.call(rbind, area_list)

write.csv(cell_areas_2d,'./output/Xenium/cell_areas_2d.csv')

# Add 2D areas to metadata
xenium.obj$cell_area_2d <- cell_areas_2d$cell_area_2d[match(colnames(xenium.obj), cell_areas_2d$cell)]

# Compute area fraction: sum of 2D cell area per type / total 2D cell area per patient
area_by_type <- aggregate(cell_area_2d ~ patient + cell_type_rctd_doublet,
                          data = xenium.obj@meta.data, FUN = sum, na.rm = TRUE)

total_cell_area <- aggregate(cell_area_2d ~ patient,
                             data = xenium.obj@meta.data, FUN = sum, na.rm = TRUE)
names(total_cell_area)[2] <- "total_cell_area_2d"

area_by_type <- merge(area_by_type, total_cell_area, by = "patient")
area_by_type$area_fraction <- area_by_type$cell_area_2d / area_by_type$total_cell_area_2d

# Add group info
patient_meta <- unique(xenium.obj@meta.data[, c("patient", "group")])
area_by_type <- merge(area_by_type, patient_meta, by = "patient")
area_by_type$group <- factor(area_by_type$group, levels = c("NF", "pRV", "RVF"))

pdf('./output/Xenium/cell_area_norm.pdf', width = 6, height = 6)
ggplot(area_by_type, aes(x = group, y = area_fraction, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.5) +
  facet_wrap(~cell_type_rctd_doublet, scales = "free_y") +
  theme_bw() +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Fraction of total segmented area (2D)", x = NULL,
       title = "Cell type area coverage by group")
dev.off()

##############################################
#### 8. Myeloid subclust
##############################################

ref <- readRDS('../snRV_ref.rds')                          
myeloid.ref <- subset(ref,subset = Names %in% c('Myeloid'))



myeloid.obj <- subset(xenium.obj, subset = cell_type_rctd_doublet == "Myeloid")
myeloid.obj <- JoinLayers(myeloid.obj)

# Purge stale UMAP coordinate columns inherited from the parent — DimPlot's
# FetchData checks meta.data before @reductions, so leftover umap_1/umap_2
# shadow any freshly computed subcluster UMAP.
umap_cols <- grep("^umap_\\d+$|^UMAP_\\d+$|^umaporig_\\d+$",
                  colnames(myeloid.obj@meta.data), value = TRUE)
if (length(umap_cols)) myeloid.obj@meta.data[umap_cols] <- NULL

myeloid.obj[["pca"]]     <- NULL
myeloid.obj[["harmony"]]     <- NULL
myeloid.obj[["umap"]]     <- NULL
myeloid.obj[["umap_orig"]]     <- NULL

myeloid.obj <- FindVariableFeatures(myeloid.obj)
myeloid.obj <- ScaleData(myeloid.obj)
myeloid.obj <- RunPCA(myeloid.obj, npcs = 30)
myeloid.obj <- RunHarmony(myeloid.obj, "patient")

myeloid.obj <- RunUMAP(myeloid.obj, reduction = "harmony", dims = 1:30,reduction.name='umap')
myeloid.obj <- RunUMAP(myeloid.obj, reduction = "pca", dims = 1:30,reduction.name='umap_orig')

myeloid.obj <- FindNeighbors(myeloid.obj, reduction = "harmony", dims = 1:30)
myeloid.obj <- FindClusters(myeloid.obj, resolution = 0.5)


DimPlot(myeloid.obj,reduction='umap_orig',label=T)

m.marks <- FindAllMarkers(myeloid.obj)
write.csv(m.marks,'./output/Xenium/m.marks.csv')
write.csv(m.marks,'./output/Xenium/m.marks.21.csv')


myeloid_labels <- c(
  '0'  = 'Macrophage_Resident',
  '1'  = 'Macrophage_C1q',
  '2'  = 'Fibroblast',
  '3'  = 'Cardiomyocyte',
  '4'  = 'Unknown',
  '5'  = 'Monocyte',
  '6'  = 'Fibroblast_Doublet',
  '7'  = 'Neutrophil',
  '8'  = 'Pericyte',
  '9'  = 'Pericyte',
  '10' = 'Fibroblast',
  '11' = 'Cardiomyocyte',
  '12' = 'Adipocyte',
  '13' = 'TNK_Cell',
  '14' = 'Basophil',
  '15' = 'Macrophage_Inflammatory',
  '16' = 'Myeloid_Proliferating',
  '17' = 'Unknown',
  '18' = 'Schwann_Cell',
  '19' = 'LymphEndothelial',
  '20' = 'Neutrophil',
  '21' = 'Unknown'
)

  # ┌─────────┬─────────────────────────┬──────────────────────────────────────────────┐                
  # │ Cluster │        Cell Type        │                 Key Markers                  │                
  # ├─────────┼─────────────────────────┼──────────────────────────────────────────────┤                
  # │ 0       │ Macrophage_Resident     │ PLA2G15, HPS1, WASHC2A (lysosomal/endosomal) │                
  # ├─────────┼─────────────────────────┼──────────────────────────────────────────────┤                
  # │ 1       │ Macrophage_C1q          │ C1QA, CD33                                   │                
  # ├─────────┼─────────────────────────┼──────────────────────────────────────────────┤                
  # │ 2       │ Fibroblast              │ CILP, FNDC1, SERPINA3                        │                
  # ├─────────┼─────────────────────────┼──────────────────────────────────────────────┤                
  # │ 3       │ Cardiomyocyte           │ HSPB3, CAVIN4                                │                
  # ├─────────┼─────────────────────────┼──────────────────────────────────────────────┤                
  # │ 4       │ Unknown                 │ SEMA3G, SEMA3F, ADGRG1                       │                
  # ├─────────┼─────────────────────────┼──────────────────────────────────────────────┤                
  # │ 5       │ Monocyte                │ SIRPD, CD300LF, IFI30, LTA4H                 │                
  # ├─────────┼─────────────────────────┼──────────────────────────────────────────────┤                
  # │ 6       │ Fibroblast_Doublet      │ SFRP2, RPL10 (ribosomal)                     │                
  # ├─────────┼─────────────────────────┼──────────────────────────────────────────────┤                
  # │ 7       │ Neutrophil              │ HVCN1, SLC22A4, SLC22A5                      │                
  # ├─────────┼─────────────────────────┼──────────────────────────────────────────────┤                
  # │ 8       │ Pericyte                │ RGS5, ENPEP, NOTCH3                          │                
  # ├─────────┼─────────────────────────┼──────────────────────────────────────────────┤                
  # │ 9       │ Pericyte                │ KCNJ8, KCNE4                                 │                
  # ├─────────┼─────────────────────────┼──────────────────────────────────────────────┤                
  # │ 10      │ Fibroblast              │ WFDC1, CSRP1                                 │                
  # ├─────────┼─────────────────────────┼──────────────────────────────────────────────┤                
  # │ 11      │ Cardiomyocyte           │ SCN5A, CKMT2, RGS5                           │
  # ├─────────┼─────────────────────────┼──────────────────────────────────────────────┤
  # │ 12      │ Adipocyte               │ GPD1, CIDEA, PFKFB1                          │
  # ├─────────┼─────────────────────────┼──────────────────────────────────────────────┤
  # │ 13      │ TNK_Cell                │ RUNX3, CD6, KLRK1                            │
  # ├─────────┼─────────────────────────┼──────────────────────────────────────────────┤
  # │ 14      │ Basophil                │ HDC, HPGD, RAB27B                            │
  # ├─────────┼─────────────────────────┼──────────────────────────────────────────────┤
  # │ 15      │ Macrophage_Inflammatory │ CCL3, CCL4, CCL8, EGR2, IER3                 │
  # ├─────────┼─────────────────────────┼──────────────────────────────────────────────┤
  # │ 16      │ Myeloid_Proliferating   │ BUB1B, TPX2, ANLN, KIF11, MELK               │
  # ├─────────┼─────────────────────────┼──────────────────────────────────────────────┤
  # │ 17      │ Unknown                 │ PTGS1, SLC6A4, COLEC11                       │
  # ├─────────┼─────────────────────────┼──────────────────────────────────────────────┤
  # │ 18      │ Schwann_Cell            │ SOX10, GRIK3                                 │
  # ├─────────┼─────────────────────────┼──────────────────────────────────────────────┤
  # │ 19      │ LymphEndothelial        │ STAB2, CCL21                                 │
  # ├─────────┼─────────────────────────┼──────────────────────────────────────────────┤
  # │ 20      │ Neutrophil              │ HVCN1, SLC22A4, SLC22A5, PFKM                │
  # └─────────┴─────────────────────────┴──────────────────────────────────────────────┘

myeloid.obj$myeloid_subtype <- unname(myeloid_labels[as.character(myeloid.obj$seurat_clusters)])

DimPlot(myeloid.obj, group.by = 'myeloid_subtype', label = TRUE, repel = TRUE,reduction='umap_orig') +
  theme(legend.position = 'none')

myeloid.obj.clean <- subset(myeloid.obj,subset = myeloid_subtype %in% c('Macrophage_Resident','Macrophage_C1q','Monocyte','Macrophage_Inflammatory','Myeloid_Proliferating','Neutrophil','Basophil'))
tnk.obj.clean <- subset(myeloid.obj,subset = myeloid_subtype %in% c('TNK_Cell'))

myeloid.obj.clean[["pca"]]     <- NULL  
myeloid.obj.clean[["harmony"]]     <- NULL
myeloid.obj.clean[["umap"]]     <- NULL

umap_cols <- grep("^umap_\\d+$|^UMAP_\\d+$|^umaporig_\\d+$",
                  colnames(myeloid.obj.clean@meta.data), value = TRUE)
if (length(umap_cols)) myeloid.obj.clean@meta.data[umap_cols] <- NULL

myeloid.obj.clean <- JoinLayers(myeloid.obj.clean)
myeloid.obj.clean <- FindVariableFeatures(myeloid.obj.clean)
myeloid.obj.clean <- ScaleData(myeloid.obj.clean,split.by='patient')
myeloid.obj.clean <- RunPCA(myeloid.obj.clean, npcs = 30)
myeloid.obj.clean <- RunHarmony(myeloid.obj.clean, "patient")

myeloid.obj.clean <- RunUMAP(myeloid.obj.clean, reduction = "harmony", dims = 1:30,reduction.name='umap')
myeloid.obj.clean <- RunUMAP(myeloid.obj.clean, reduction = "pca", dims = 1:30,reduction.name='umap_orig')

DimPlot(myeloid.obj.clean,group.by='myeloid_subtype',label=T,reduction='umap_orig')

pdf(paste0('./output/Xenium/', 'Xenium_myeloid_snUMAP_new.pdf'), width=5, height=5)
PlotEmbedding(myeloid.obj.clean,group.by='myeloid_subtype',reduction='umap_orig',point_size=.05,plot_under=TRUE,plot_theme=umap_theme()+NoLegend(),raster_dpi=400,raster_scale=0.5)
dev.off()

pdf(paste0('./output/Xenium/', 'Xenium_myeloid_snUMAP_new_pt.pdf'), width=5, height=5)
PlotEmbedding(myeloid.obj.clean,group.by='patient',reduction='umap_orig',point_size=.05,plot_under=TRUE,plot_theme=umap_theme()+NoLegend(),raster_dpi=400,raster_scale=0.5)
dev.off()






# CCR2+ marks
FeaturePlot(myeloid.obj.clean,
              features = c("FCN1","VCAN","PLAC8","SPP1","S100A8","CCR2"),
              ncol = 3)                                                                               
DotPlot(myeloid.obj.clean,                                                                          
          features = c("FCN1","VCAN","PLAC8","SPP1","S100A8","CCR2"),                                 
          group.by = "myeloid_subtype") + RotatedAxis() 
# Some of the monocytes are actually CCR2+ macrophages



myeloid.obj.clean <- FindNeighbors(myeloid.obj.clean, reduction = "harmony", dims = 1:30)
myeloid.obj.clean <- FindClusters(myeloid.obj.clean, resolution = 0.5)


myeloid.obj.clean.mem <- myeloid.obj.clean
myeloid.obj.clean.mem[["Xenium"]] <- CreateAssayObject(                                                  
    counts = as(GetAssayData(myeloid.obj.clean, layer = "counts"), "dgCMatrix")                            
)
myeloid.marks <- FindAllMarkers(myeloid.obj.clean.mem)                                                        
rm(myeloid.obj.clean.mem)  
write.csv(myeloid.marks, './output/Xenium/myeloid.clean.marks.csv')

# ── Annotate myeloid clean clusters ──────────────────────────────────────────
#
#   ┌─────────┬───────────────────────────┬───────────────────────────────────────┐
#   │ Cluster │         Cell Type         │             Key Markers               │
#   ├─────────┼───────────────────────────┼───────────────────────────────────────┤
#   │ 0       │ Macrophage_Resident       │ MARCO, MS4A6A, MS4A4A, CD14           │
#   │ 1       │ Macrophage_C1q            │ C1QA, CEBPD                           │
#   │ 2       │ Monocyte                  │ FCN1, VCAN, S100A12, MNDA, IL1R2      │
#   │ 3       │ Dendritic_Cell            │ CD209, CCL18, IKZF2                   │
#   │ 4       │ Neutrophil                │ HVCN1, ALS2CL, CPNE5                  │
#   │ 5       │ Macrophage_Resident_LYVE1 │ MRC1, LYVE1                           │
#   │ 6       │ Cardiomyocyte             │ ACTA1, MYL2, ACTC1, TNNI3, MYH7      │
#   │ 7       │ Schwann_Cell              │ SCN1B, GALC, SCRG1, EPHB2             │
#   │ 8       │ Dendritic_Cell            │ CD1C, CLEC10A, FCER1A, CD1E (cDC2)   │
#   │ 9       │ Myeloid_Proliferating     │ MKI67, TOP2A, CENPF, PCNA             │
#   │ 10      │ Macrophage_Inflammatory   │ CXCL8, CXCL2, CCL2, OLR1, NR4A1     │
#   │ 11      │ Mast_Cell                 │ CPA3, KIT, CTSG, MS4A2               │
#   │ 12      │ Macrophage_TREM2          │ HAVCR2, CTSK, LPL, LILRB4, PLA2G7   │
#   │ 13      │ CD8_T                     │ CD8A, TRAC, GZMK, PRF1, CD28         │
#   └─────────┴───────────────────────────┴───────────────────────────────────────┘
#   True myeloid: 0,1,2,3,5,8,9,10,11,12
#   Contamination/Doublets: 4 (Doublet: Neutrophil+CM), 5 (Doublet: LYVE1-Mac+CM), 13 (CD8_T)
#   Note: clusters 6+7 (MBP+/MRC1+/CD68+ with sarcomeric/myelin DE genes) = phagocytic
#         Macrophage_Resident that have engulfed CM/myelin debris, merged with cluster 0

myeloid_clean_labels <- c(
  '0'  = 'Macrophage_Resident',
  '1'  = 'Macrophage_C1q',
  '2'  = 'Monocyte',
  '3'  = 'Dendritic_Cell',
  '4'  = 'Doublet',
  '5'  = 'Doublet',
  '6'  = 'Macrophage_Resident',
  '7'  = 'Macrophage_Resident',
  '8'  = 'Dendritic_Cell',
  '9'  = 'Myeloid_Proliferating',
  '10' = 'Macrophage_Inflammatory',
  '11' = 'Mast_Cell',
  '12' = 'Macrophage_TREM2',
  '13' = 'CD8_T'
)

myeloid.obj.clean$myeloid_subtype <- unname(myeloid_clean_labels[as.character(myeloid.obj.clean$seurat_clusters)])



saveRDS(myeloid.obj.clean,'./output/Xenium/myeloid_clean.rds')


myeloid.obj.clean.clean <- subset(myeloid.obj.clean,
  subset = myeloid_subtype %in% c('Macrophage_Resident', 'Macrophage_C1q', 'Monocyte',
                                   'Dendritic_Cell', 'Myeloid_Proliferating',
                                   'Macrophage_Inflammatory', 'Mast_Cell', 'Macrophage_TREM2'))

myeloid.obj.clean.clean[["pca"]]      <- NULL
myeloid.obj.clean.clean[["harmony"]]  <- NULL
myeloid.obj.clean.clean[["umap"]]     <- NULL
myeloid.obj.clean.clean[["umap_orig"]] <- NULL

umap_cols <- grep("^umap_\\d+$|^UMAP_\\d+$|^umaporig_\\d+$",
                  colnames(myeloid.obj.clean.clean@meta.data), value = TRUE)
if (length(umap_cols)) myeloid.obj.clean.clean@meta.data[umap_cols] <- NULL

myeloid.obj.clean.clean <- JoinLayers(myeloid.obj.clean.clean)
myeloid.obj.clean.clean <- FindVariableFeatures(myeloid.obj.clean.clean)
myeloid.obj.clean.clean <- ScaleData(myeloid.obj.clean.clean)
myeloid.obj.clean.clean <- RunPCA(myeloid.obj.clean.clean, npcs = 30)
myeloid.obj.clean.clean <- RunHarmony(myeloid.obj.clean.clean, "patient")
myeloid.obj.clean.clean <- RunUMAP(myeloid.obj.clean.clean, reduction = "harmony", dims = 1:30, reduction.name = 'umap')
myeloid.obj.clean.clean <- RunUMAP(myeloid.obj.clean.clean, reduction = "pca",     dims = 1:30, reduction.name = 'umap_orig')

myeloid.obj.clean.clean <- FindNeighbors(myeloid.obj.clean.clean, reduction = "harmony", dims = 1:30)
myeloid.obj.clean.clean <- FindClusters(myeloid.obj.clean.clean, resolution = 2)


DimPlot(myeloid.obj.clean.clean, group.by = 'myeloid_subtype', label = TRUE, repel = TRUE, reduction = 'umap_orig') +
  theme(legend.position = 'none')
DimPlot(myeloid.obj.clean.clean, group.by = 'patient', reduction = 'umap_orig')


myeloid.obj.clean.clean.mem <- myeloid.obj.clean.clean
myeloid.obj.clean.clean.mem[["Xenium"]] <- CreateAssayObject(                                                  
    counts = as(GetAssayData(myeloid.obj.clean.clean, layer = "counts"), "dgCMatrix")                            
)
myeloid.marks <- FindAllMarkers(myeloid.obj.clean.clean.mem)                                                        
rm(myeloid.obj.clean.clean.mem)  
write.csv(myeloid.marks, './output/Xenium/myeloid.clean.clean.marks.csv')


# 1. Pan-myeloid vs contamination check for ambiguous clusters   

DotPlot(myeloid.obj.clean.clean,                                                                    
features = c("CD68","MRC1","PTPRC","MARCO","CSF1R",  # macrophage                                 
             "ACTC1","MYH7","TNNI3",                  # CM                                        
             "SOX10","MPZ","S100B",                    # Schwann                                  
             "ADIPOQ","FABP4",                         # adipocyte                                
             "CNN1","ACTA2",                           # smooth muscle                            
             "DCN","COL1A1"),                          # fibroblast                               
idents = c(0,1,5,6,7,14,17,19,20,21,22,26,27,30),                                                 
group.by = "seurat_clusters") + RotatedAxis() + ggtitle("Ambiguous clusters: myeloid vs           
contamination")     
          


# 2. T/NK check                                                                                     
DotPlot(myeloid.obj.clean.clean,                                                                    
features = c("CD3E","CD3D","CD8A","CD4","GZMK","NKG7","GNLY",                                     
             "CCL5","GIMAP5","THEMIS"),                                                           
idents = c(1,19,29),                                    
group.by = "seurat_clusters") + RotatedAxis() + ggtitle("T/NK check") 

# 3. Resident mac subtypes                                                                          
DotPlot(myeloid.obj.clean.clean,                                                                    
features = c("LYVE1","MRC1","VSIG4","MARCO","CD163","FOLR2","ALOX15B",
             "C1QA","C1QB","SIGLEC1"),                                                            
idents = c(0,5,7,9,10,22,27),                           
group.by = "seurat_clusters") + RotatedAxis() + ggtitle("Resident mac subtypes")                  



# 4. Cluster 26 — TREM2+inflammatory vs fibroblast doublet                                          
DotPlot(myeloid.obj.clean.clean,                                                                    
features = c("TREM2","CTSK","RGCC","PLA2G7","CXCL8",                                              
             "DCN","COL1A1","FBLN1","LUM",                                                        
             "CD68","PTPRC"),                                                                     
idents = c(25,26),                                                                                
group.by = "seurat_clusters") + RotatedAxis() + ggtitle("Cluster 26 vs 25 (TREM2)")  

                                                                                                  
# 5. Cluster 6 — neuronal or IFN-stimulated mac?   
                                             
DotPlot(myeloid.obj.clean.clean,                                                                    
features = c("SYT1","SNAP25","RBFOX3",   # neuronal                                               
             "ISG15","MX1","IFI44L","IFIT1", # ISG                                                
             "CD68","PTPRC","CSF1R"),                                                             
idents = c(6,20),                                                                                 
group.by = "seurat_clusters") + RotatedAxis() + ggtitle("Cluster 6 & 20: neuronal vs ISG-mac")    

# 6. Cluster 30 — endothelial-DC doublet?                                                           
DotPlot(myeloid.obj.clean.clean,                                                                    
features = c("PECAM1","VWF","CDH5","FLT1",  # endothelial
             "CD209","IKZF2","HLA-DRA",       # DC                                                
             "PTPRC","CD68"),                           
idents = c(12,30),                                                                                
group.by = "seurat_clusters") + RotatedAxis() + ggtitle("Cluster 30 vs 12 (DC)") 

pdf('./output/Xenium/test.pdf', width = 10, height = 8)                                                                                               

DotPlot(myeloid.obj.clean.clean,
features = c("CD1C","FCER1A","CLEC10A","HLA-DRA",   # DC                                          
             "CD68","MRC1","CSF1R","MARCO","LYVE1",  # macrophage
             "PECAM1","EMCN","VWF",                  # endothelial                                
             "PTPRC"),                                                                            
idents = c(12,16,18),                                                                             
group.by = "seurat_clusters") + RotatedAxis() + ggtitle("Cluster 12: DC or macrophage?")  
        dev.off()                                                                                
 


#   ┌─────────┬───────────────────────────┬──────────────────────────────────────────┐
#   │ Cluster │         Cell Type         │             Key Markers / Evidence       │
#   ├─────────┼───────────────────────────┼──────────────────────────────────────────┤
#   │ 0       │ Macrophage_Resident       │ LYVE1/MRC1/VSIG4/CD163 (DotPlot)        │
#   │ 1       │ Macrophage_Resident       │ Mac markers+, no T/NK (DotPlot)         │
#   │ 2       │ Macrophage_C1q            │ C1QA top marker                          │
#   │ 3       │ Macrophage_Resident       │ MARCO/LYVE1/VSIG4 tissue-resident        │
#   │ 4       │ Macrophage_Resident       │ GALC/TCN2 (phagocytic)                   │
#   │ 5       │ Macrophage_Resident       │ Full resident panel (DotPlot)            │
#   │ 6       │ Macrophage_Resident       │ CSF1R-dominant (DotPlot)                 │
#   │ 7       │ Macrophage_Resident       │ Resident panel confirmed (DotPlot)       │
#   │ 8       │ Macrophage_Monocyte_Derived│ IL1R2/EREG/AREG; MoMF reparative         │
#   │ 9       │ Macrophage_Resident       │ ALOX15B+ M2 (DotPlot)                   │
#   │ 10      │ Macrophage_Resident       │ SCRG1/ZC3H12D, weak resident (DotPlot)  │
#   │ 11      │ Monocyte                  │ VCAN/S100A12/CCR2/SELL                   │
#   │ 12      │ Macrophage_Resident       │ CD209/CSF1R+; no CD1C/FCER1A (not DC)    │
#   │ 13      │ Adipocyte                 │ ADIPOQ/PPARG/FABP4 — REMOVE             │
#   │ 14      │ Macrophage_Resident       │ CM debris + MRC1/CSF1R (DotPlot)        │
#   │ 15      │ Myeloid_Proliferating     │ MKI67/TOP2A/PCNA                         │
#   │ 16      │ Dendritic_Cell            │ CD1C/FCER1A/CLEC10A (cDC2)               │
#   │ 17      │ Macrophage_Resident       │ CNN1 top but MRC1/CSF1R (DotPlot)       │
#   │ 18      │ Dendritic_Cell            │ IRF8/CLEC9A/CADM1 (cDC1)                │
#   │ 19      │ Doublet                   │ Mac+T (CD4 bright) + DCN — REMOVE       │
#   │ 20      │ Doublet                   │ Neuronal (SYT1/DLGAP2/GPM6A) — REMOVE   │
#   │ 21      │ Macrophage_Resident       │ Weak mac markers, generic DE genes       │
#   │ 22      │ Macrophage_Resident       │ LYVE1/SIGLEC1 bright (DotPlot)          │
#   │ 23      │ Mast_Cell                 │ CPA3/KIT/MS4A2/GATA2                     │
#   │ 24      │ Macrophage_Inflammatory   │ NR4A1/CCL3/CCL8/CD83                    │
#   │ 25      │ Macrophage_TREM2          │ CTSK/TREM2 (DotPlot confirmed)           │
#   │ 26      │ Doublet                   │ Mac+Fibroblast (DCN/COL1A1) — REMOVE    │
#   │ 27      │ Macrophage_Resident       │ Strongest resident panel (DotPlot)       │
#   │ 28      │ Low_Quality               │ NegControlCodeword top marker — REMOVE   │
#   │ 29      │ Doublet                   │ Neuronal (SEMA3D/NTRK3) — REMOVE        │
#   │ 30      │ Doublet                   │ EC+DC (PECAM1/VWF+CD209) — REMOVE       │
#   │ 31      │ Macrophage_Resident       │ CCL18/TEK                                │
#   └─────────┴───────────────────────────┴──────────────────────────────────────────┘

myeloid_clean_clean_labels <- c(
  '0'  = 'Macrophage_Resident',
  '1'  = 'Macrophage_Resident',
  '2'  = 'Macrophage_C1q',
  '3'  = 'Macrophage_Resident',
  '4'  = 'Macrophage_Resident',
  '5'  = 'Macrophage_Resident',
  '6'  = 'Macrophage_Resident',
  '7'  = 'Macrophage_Resident',
  '8'  = 'Macrophage_Monocyte_Derived',
  '9'  = 'Macrophage_Resident',
  '10' = 'Macrophage_Resident',
  '11' = 'Monocyte',
  '12' = 'Macrophage_Resident',
  '13' = 'Adipocyte',
  '14' = 'Macrophage_Resident',
  '15' = 'Myeloid_Proliferating',
  '16' = 'Dendritic_Cell',
  '17' = 'Macrophage_Resident',
  '18' = 'Dendritic_Cell',
  '19' = 'Doublet',
  '20' = 'Doublet',
  '21' = 'Macrophage_Resident',
  '22' = 'Macrophage_Resident',
  '23' = 'Mast_Cell',
  '24' = 'Macrophage_Inflammatory',
  '25' = 'Macrophage_TREM2',
  '26' = 'Doublet',
  '27' = 'Macrophage_Resident',
  '28' = 'Low_Quality',
  '29' = 'Doublet',
  '30' = 'Doublet',
  '31' = 'Macrophage_Resident'
)

myeloid.obj.clean.clean$myeloid_clean_subtype <- unname(
  myeloid_clean_clean_labels[as.character(myeloid.obj.clean.clean$seurat_clusters)]
)

# Remove contamination / doublets / low-quality clusters
myeloid.obj.clean.clean <- subset(myeloid.obj.clean.clean,
  subset = myeloid_clean_subtype %in% c('Macrophage_Resident', 'Macrophage_C1q', 'Monocyte',
                                         'Dendritic_Cell', 'Myeloid_Proliferating',
                                         'Macrophage_Inflammatory', 'Macrophage_Monocyte_Derived',
                                         'Mast_Cell', 'Macrophage_TREM2'))




pdf('./output/Xenium/test.pdf', width = 10, height = 8)
DimPlot(myeloid.obj.clean.clean, group.by = 'myeloid_subtype', label = TRUE, repel = TRUE, reduction = 'umap_orig') +
  theme(legend.position = 'none')
dev.off()





# Label transfer from snRV reference
myeloid.ref <- JoinLayers(myeloid.ref)
myeloid.anchors <- FindTransferAnchors(
  reference   = myeloid.ref,
  query       = myeloid.obj.clean.clean,
  dims        = 1:30,
  reference.reduction = 'pca'
)
myeloid.predictions <- TransferData(
  anchorset  = myeloid.anchors,
  refdata    = myeloid.ref$Subnames,
  dims       = 1:30
)
myeloid.obj.clean.clean <- AddMetaData(myeloid.obj.clean.clean, metadata = myeloid.predictions)

saveRDS(myeloid.obj.clean.clean, './output/Xenium/myeloid_clean_clean.rds')


pdf('./output/Xenium/test.pdf', width = 10, height = 8)
DimPlot(myeloid.obj.clean.clean, group.by = 'seurat_clusters', label = TRUE, repel = TRUE, reduction = 'umap_orig') +
  theme(legend.position = 'none') +
  ggtitle('Reference predicted labels')
dev.off()

# Reviewer DotPlot: canonical markers by myeloid subtype
myeloid_marker_genes <- c(
  # Dendritic cell
  'CD1C', 'FCER1A', 'CLEC10A', 'IRF8', 'CLEC9A',
  # C1q macrophage
  'C1QA', 'C1QB',
  # Inflammatory macrophage
  'CCL3', 'CCL8', 'NR4A1', 'CD83',
  # Monocyte-derived macrophage
  'IL1R2', 'EREG', 'AREG', 'SERPINB2',
  # Tissue-resident macrophage
  'LYVE1', 'MRC1', 'FOLR2', 'VSIG4', 'MARCO', 'CD163',
  # TREM2 disease-associated macrophage
  'TREM2', 'CTSK', 'HAVCR2', 'RGCC',
  # Mast cell
  'CPA3', 'KIT', 'MS4A2',
  # Monocyte
  'VCAN', 'S100A12', 'CCR2', 'FCN1',
  # Proliferating
  'MKI67', 'TOP2A',
  # Pan-myeloid
  'CD68', 'PTPRC', 'CSF1R'
)

myeloid_subtype_order <- c(
  'Dendritic_Cell', 'Macrophage_C1q', 'Macrophage_Inflammatory',
  'Macrophage_Monocyte_Derived', 'Macrophage_Resident', 'Macrophage_TREM2',
  'Mast_Cell', 'Monocyte', 'Myeloid_Proliferating'
)
myeloid.obj.clean.clean$myeloid_clean_subtype <- factor(
  myeloid.obj.clean.clean$myeloid_clean_subtype, levels = myeloid_subtype_order
)

pdf('./output/Xenium/myeloid_clean_clean_markers_dotplot.pdf', width = 16, height = 6)
DotPlot(myeloid.obj.clean.clean,
  features = myeloid_marker_genes,
  group.by = 'myeloid_clean_subtype',col.min=0,col.max=2) +
  RotatedAxis() +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 9)) +
  ggtitle('Myeloid subtype canonical markers')
dev.off()









# Myeloid subtype proportions by disease state
myeloid_meta <- myeloid.obj.clean.clean@meta.data %>%
  filter(!is.na(myeloid_clean_subtype)) %>%
  select(patient, group, myeloid_clean_subtype)

myeloid_totals <- myeloid_meta %>%
  group_by(patient, group) %>%
  summarise(total = n(), .groups = 'drop')

myeloid_prop <- myeloid_meta %>%
  group_by(patient, group, myeloid_clean_subtype) %>%
  summarise(n = n(), .groups = 'drop') %>%
  left_join(myeloid_totals, by = c('patient', 'group')) %>%
  mutate(proportion = n / total) %>%
  complete(nesting(patient, group), myeloid_clean_subtype, fill = list(n = 0, proportion = 0)) %>%
  mutate(group = factor(group, levels = c('NF', 'pRV', 'RVF')))

pairwise_comparisons <- list(c('NF', 'pRV'), c('NF', 'RVF'), c('pRV', 'RVF'))

p_myeloid_prop <- ggplot(myeloid_prop, aes(x = group, y = proportion, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.5) +
  stat_compare_means(comparisons = pairwise_comparisons, method = 'wilcox.test', label = 'p.signif') +
  stat_compare_means(method = 'kruskal.test', label.y.npc = 0.95, size = 3) +
  facet_wrap(~ myeloid_clean_subtype, scales = 'free_y') +
  scale_x_discrete(limits = c('NF', 'pRV', 'RVF')) +
  theme_bw() +
  theme(legend.position = 'none',
        strip.text = element_text(face = 'bold'),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = 'Proportion of myeloid cells', title = 'Myeloid subtype distribution by group')

pdf('./output/Xenium/myeloid_proportions.pdf', width = 10, height = 8)
print(p_myeloid_prop)
dev.off()


##############################################
#### Pseudobulk DESeq2 by myeloid subtype
##############################################

library(DESeq2)

# Build pseudobulk metadata directly from cell-level metadata (avoids column name parsing)
pseudo_meta <- unique(myeloid.obj.clean.clean@meta.data[, c('patient', 'group', 'myeloid_clean_subtype')])
pseudo_meta$patient <- as.character(pseudo_meta$patient)
pseudo_meta$myeloid_clean_subtype <- as.character(pseudo_meta$myeloid_clean_subtype)
pseudo_meta$pseudobulk_id <- paste0(pseudo_meta$patient, '..', pseudo_meta$myeloid_clean_subtype)
pseudo_meta$group <- factor(pseudo_meta$group, levels = c('NF', 'pRV', 'RVF'))
rownames(pseudo_meta) <- pseudo_meta$pseudobulk_id

# Pseudobulk: aggregate raw counts by patient × cell type (use '..' separator to avoid conflicts)
myeloid.obj.clean.clean$pseudobulk_id <- paste0(
  as.character(myeloid.obj.clean.clean$patient), '..',
  as.character(myeloid.obj.clean.clean$myeloid_clean_subtype)
)

pseudo_counts <- AggregateExpression(
  myeloid.obj.clean.clean,
  assays = 'Xenium',
  group.by = 'pseudobulk_id',
  slot = 'counts',
  return.seurat = FALSE
)$Xenium

# AggregateExpression prepends 'g' and converts '_' to '-' — normalize both
colnames(pseudo_counts) <- gsub('-', '_', sub('^g', '', colnames(pseudo_counts)))

# Run DESeq2 per cell type: RVF vs NF
cell_types <- unique(pseudo_meta$myeloid_clean_subtype)
deseq_results <- list()

for (ct in cell_types) {
  ct_samples <- pseudo_meta$pseudobulk_id[pseudo_meta$myeloid_clean_subtype == ct]
  # Keep only samples that exist in pseudo_counts (some patient×celltype combos may have no cells)
  ct_samples <- ct_samples[ct_samples %in% colnames(pseudo_counts)]
  if (length(ct_samples) < 2) next

  ct_counts <- pseudo_counts[, ct_samples, drop = FALSE]
  ct_meta   <- pseudo_meta[ct_samples, , drop = FALSE]

  # Need at least 2 groups represented
  if (length(unique(ct_meta$group)) < 2) next

  # Filter low-count genes
  keep <- rowSums(ct_counts >= 5) >= 2
  ct_counts <- ct_counts[keep, , drop = FALSE]
  if (nrow(ct_counts) < 10) next

  dds <- DESeqDataSetFromMatrix(
    countData = ct_counts,
    colData   = ct_meta,
    design    = ~ group
  )
  dds <- tryCatch(DESeq(dds), error = function(e) NULL)
  if (is.null(dds)) next

  # RVF vs NF
  res_rvf <- tryCatch(
    results(dds, contrast = c('group', 'RVF', 'NF')) %>%
      as.data.frame() %>%
      tibble::rownames_to_column('gene') %>%
      mutate(cell_type = ct, comparison = 'RVF_vs_NF') %>%
      arrange(padj),
    error = function(e) NULL
  )

  # pRV vs NF
  res_prv <- tryCatch(
    results(dds, contrast = c('group', 'pRV', 'NF')) %>%
      as.data.frame() %>%
      tibble::rownames_to_column('gene') %>%
      mutate(cell_type = ct, comparison = 'pRV_vs_NF') %>%
      arrange(padj),
    error = function(e) NULL
  )

  # RVF vs pRV
  res_rvf_prv <- tryCatch(
    results(dds, contrast = c('group', 'RVF', 'pRV')) %>%
      as.data.frame() %>%
      tibble::rownames_to_column('gene') %>%
      mutate(cell_type = ct, comparison = 'RVF_vs_pRV') %>%
      arrange(padj),
    error = function(e) NULL
  )

  if (!is.null(res_rvf))     deseq_results[[paste0(ct, '_RVF_vs_NF')]]  <- res_rvf
  if (!is.null(res_prv))     deseq_results[[paste0(ct, '_pRV_vs_NF')]]  <- res_prv
  if (!is.null(res_rvf_prv)) deseq_results[[paste0(ct, '_RVF_vs_pRV')]] <- res_rvf_prv
}

# Combine and save
all_deseq <- bind_rows(deseq_results)
write.csv(all_deseq, './output/Xenium/myeloid_pseudobulk_deseq2.csv', row.names = FALSE)

# Summary: significant DEGs per cell type and comparison
sig_summary <- all_deseq %>%
  filter(padj < 0.05, abs(log2FoldChange) > 0.5) %>%
  group_by(cell_type, comparison) %>%
  summarise(n_up = sum(log2FoldChange > 0),
            n_down = sum(log2FoldChange < 0),
            .groups = 'drop')
print(sig_summary)


all_deseq %>% filter(padj < 0.05, abs(log2FoldChange) > 0.5) %>%
    filter(cell_type == 'Dendritic_Cell')

all_deseq %>% filter(padj < 0.05, abs(log2FoldChange) > 0.5) %>%
    filter(cell_type == 'Monocyte')

all_deseq %>% filter(padj < 0.05, abs(log2FoldChange) > 0.5) %>%
    filter(cell_type == 'Macrophage_Monocyte_Derived')

all_deseq %>% filter(padj < 0.05, abs(log2FoldChange) > 0.5) %>%
    filter(cell_type == 'Macrophage_Resident')


##############################################
#### Myeloid populations: Volcano + GO/Reactome (all cell types)
##############################################

library(EnhancedVolcano)
library(enrichR)

# Load original panel genes from corrected (non-imputed) object
xenium.base <- readRDS('~/Documents/XeniumWorkflow/functions/Xenium_resegmented_corrected.rds')
panel_genes <- rownames(xenium.base)
rm(xenium.base); gc()

# Helper: cap -log10(padj) at 5 by flooring padj at 1e-5
cap_pval <- function(df) mutate(df, padj = ifelse(is.na(padj), NA_real_, pmax(padj, 1e-5)))

# Helper: build custom color vector — red for panel genes, blue for imputed, grey for NS
make_volcano_colors <- function(df, panel_genes) {
  is_panel <- df$gene %in% panel_genes
  is_sig   <- !is.na(df$padj) & df$padj < 0.05 & abs(df$log2FoldChange) > 0.5
  cols <- case_when(
    is_sig & is_panel  ~ 'red3',
    is_sig & !is_panel ~ 'royalblue3',
    !is_sig & is_panel ~ 'salmon',
    TRUE               ~ 'grey80'
  )
  names(cols) <- df$gene
  cols
}

# Helper: top N significant gene labels
sig_labs <- function(df, n = 30) {
  df %>%
    filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 0.5) %>%
    arrange(padj) %>%
    head(n) %>%
    pull(gene)
}

volcano_args <- list(
  x = 'log2FoldChange', y = 'padj',
  pCutoff = 0.05, FCcutoff = 0.5,
  pointSize = 2, labSize = 3.5,
  drawConnectors = TRUE, widthConnectors = 0.5,
  xlim = c(-5, 5), ylim = c(0, 5),
  subtitle = 'Red = original panel gene | Blue = imputed gene',
  legendPosition = 'none'
)

enrichR_dbs <- c('GO_Biological_Process_2023', 'Reactome_2022')
setEnrichrSite('Enrichr')

combine_enrich <- function(enrich_res, direction, comparison) {
  if (is.null(enrich_res)) return(NULL)
  lapply(names(enrich_res), function(db) {
    df <- enrich_res[[db]]
    if (nrow(df) == 0) return(NULL)
    df %>%
      mutate(Term = as.character(Term),
             database = db, direction = direction, comparison = comparison)
  }) %>% bind_rows()
}

plot_enrich <- function(enrich_df, title_str) {
  if (is.null(enrich_df) || nrow(enrich_df) == 0) return(NULL)
  enrich_df %>%
    arrange(Adjusted.P.value) %>%
    head(10) %>%
    mutate(Term = factor(Term, levels = rev(Term))) %>%
    ggplot(aes(x = -log10(Adjusted.P.value), y = Term)) +
    geom_col(fill = 'steelblue') +
    theme_bw() +
    labs(x = '-log10(adj. P)', y = NULL, title = title_str)
}

# Loop over all myeloid cell types
myeloid_cell_types <- unique(all_deseq$cell_type)
all_enrich_combined <- list()

for (ct in myeloid_cell_types) {
  ct_slug <- gsub('_', '-', ct)   # safe filename prefix
  cat('Processing:', ct, '\n')

  dat_rvf     <- all_deseq %>% filter(cell_type == ct, comparison == 'RVF_vs_NF')  %>% cap_pval()
  dat_prv     <- all_deseq %>% filter(cell_type == ct, comparison == 'pRV_vs_NF')  %>% cap_pval()
  dat_rvf_prv <- all_deseq %>% filter(cell_type == ct, comparison == 'RVF_vs_pRV') %>% cap_pval()

  # --- Volcano plots ---
  make_volcano_pdf <- function(dat, filepath, title_str) {
    if (nrow(dat) == 0) { cat('  Skipping', title_str, '— no data\n'); return(invisible(NULL)) }
    pdf(filepath, width = 10, height = 8)
    tryCatch(
      print(do.call(EnhancedVolcano, c(list(toptable = dat, lab = dat$gene,
        title = title_str,
        selectLab = sig_labs(dat),
        colCustom = make_volcano_colors(dat, panel_genes)), volcano_args))),
      error = function(e) { cat('  Volcano error for', title_str, ':', conditionMessage(e), '\n') }
    )
    dev.off()
  }

  make_volcano_pdf(dat_rvf,     paste0('./output/Xenium/', ct_slug, '_volcano_RVF_vs_NF.pdf'),  paste(ct, ': RVF vs NF'))
  make_volcano_pdf(dat_prv,     paste0('./output/Xenium/', ct_slug, '_volcano_pRV_vs_NF.pdf'),  paste(ct, ': pRV vs NF'))
  make_volcano_pdf(dat_rvf_prv, paste0('./output/Xenium/', ct_slug, '_volcano_RVF_vs_pRV.pdf'), paste(ct, ': RVF vs pRV'))

  # --- EnrichR ---
  run_enrich <- function(dat, comp) {
    sig <- dat %>% filter(padj < 0.05, abs(log2FoldChange) > 0.5)
    up   <- sig %>% filter(log2FoldChange > 0) %>% pull(gene)
    down <- sig %>% filter(log2FoldChange < 0) %>% pull(gene)
    list(
      up   = if (length(up)   > 0) enrichr(up,   enrichR_dbs) else NULL,
      down = if (length(down) > 0) enrichr(down, enrichR_dbs) else NULL
    )
  }

  er_rvf     <- run_enrich(dat_rvf,     'RVF_vs_NF')
  er_prv     <- run_enrich(dat_prv,     'pRV_vs_NF')
  er_rvf_prv <- run_enrich(dat_rvf_prv, 'RVF_vs_pRV')

  ct_enrich <- bind_rows(
    combine_enrich(er_rvf$up,      'up',   'RVF_vs_NF'),
    combine_enrich(er_rvf$down,    'down', 'RVF_vs_NF'),
    combine_enrich(er_prv$up,      'up',   'pRV_vs_NF'),
    combine_enrich(er_prv$down,    'down', 'pRV_vs_NF'),
    combine_enrich(er_rvf_prv$up,  'up',   'RVF_vs_pRV'),
    combine_enrich(er_rvf_prv$down,'down', 'RVF_vs_pRV')
  )
  if (!is.null(ct_enrich) && nrow(ct_enrich) > 0) {
    ct_enrich$cell_type <- ct
    all_enrich_combined[[ct]] <- ct_enrich
  }

  pdf(paste0('./output/Xenium/', ct_slug, '_enrichr_plots.pdf'), width = 12, height = 7)
  for (comp_label in c('RVF_vs_NF', 'pRV_vs_NF', 'RVF_vs_pRV')) {
    er <- switch(comp_label, RVF_vs_NF = er_rvf, pRV_vs_NF = er_prv, RVF_vs_pRV = er_rvf_prv)
    for (dir_label in c('up', 'down')) {
      er_dir <- er[[dir_label]]
      if (!is.null(er_dir)) {
        p1 <- plot_enrich(er_dir$GO_Biological_Process_2023,
          paste(ct, comp_label, toupper(dir_label), 'GO BP'))
        p2 <- plot_enrich(er_dir$Reactome_2022,
          paste(ct, comp_label, toupper(dir_label), 'Reactome'))
        if (!is.null(p1)) print(p1)
        if (!is.null(p2)) print(p2)
      }
    }
  }
  dev.off()
}

# Save combined enrichment table
all_enrich_all_ct <- bind_rows(all_enrich_combined)
write.csv(all_enrich_all_ct, './output/Xenium/myeloid_all_enrichr.csv', row.names = FALSE)


##############################################
#### CM proximity to Macrophage_Resident: pseudobulk DESeq2
##############################################

library(RANN)




# install.packages('presto')  # run once if not installed
library(presto)

# Downsample each cell_type_rctd_doublet to at most 5000 cells, then densify and run wilcoxauc
Idents(xenium.obj) <- 'cell_type_rctd_doublet'
xenium.ds <- subset(xenium.obj, downsample = 5000)
mat_ds <- as(GetAssayData(xenium.ds, assay = 'Xenium', layer = 'data'), 'dgCMatrix')
all.marks <- wilcoxauc(mat_ds, y = xenium.ds$cell_type_rctd_doublet)
rm(xenium.ds, mat_ds); gc()
write.csv(all.marks, './output/Xenium/all.marks.csv')


near_cutoff <- 15   # µm — CMs within this distance of a Macrophage_Resident
far_cutoff  <- 100   # µm — CMs beyond this distance from ANY macrophage

# --- Extract CM coordinates from xenium.obj ---
cm_cells <- WhichCells(xenium.obj, expression = cell_type_rctd_doublet == "CM")
cm_cells_seurat <- WhichCells(xenium.obj, expression = cell_type_seurat %in% c('Cm1','Cm2','Cm3','Cm4','Cm5','Cm6','Cm7','Cm8','Cm9'))
cm_cells <- intersect(cm_cells,cm_cells_seurat)

cm_coords <- as.matrix(xenium.obj@meta.data[cm_cells, c("x_centroid", "y_centroid")])

# --- Macrophage_Resident coordinates (refined annotation) for "near" group ---
#mac_res_cells <- WhichCells(myeloid.obj.clean.clean,
#  expression = myeloid_clean_subtype == "Macrophage_Resident")
#mac_res_coords <- as.matrix(myeloid.obj.clean.clean@meta.data[mac_res_cells, c("x_centroid", "y_centroid")])

# --- All macrophage coordinates (conservative union) for "far" exclusion ---
# Source 1: rMac_1 and rMac_2 from whole xenium.obj (cell_type_seurat)
broad_mac_cells <- WhichCells(xenium.obj,
  expression = cell_type_seurat %in% c("CCR2- rMac2", "CCR2- rMac1"))
# Source 2: Macrophage_Resident from myeloid.obj.clean.clean
# Union of barcodes, deduplicated
all_mac_barcodes <- broad_mac_cells#unique(c(broad_mac_cells, mac_res_cells))
all_mac_barcodes <- all_mac_barcodes[all_mac_barcodes %in% rownames(xenium.obj@meta.data)]
all_mac_coords <- as.matrix(xenium.obj@meta.data[all_mac_barcodes, c("x_centroid", "y_centroid")])
mac_res_coords <- all_mac_coords

cat("CMs:", nrow(cm_coords), "| Macrophage_Resident (refined):", nrow(mac_res_coords),
    "| All macrophages (union):", nrow(all_mac_coords), "\n")

# --- Compute nearest-neighbor distances ---
# For "near" group: distance from each CM to nearest Macrophage_Resident
nn_near <- nn2(mac_res_coords, cm_coords, k = 1)
dist_to_mac_res <- as.numeric(nn_near$nn.dists)

# For "far" group: distance from each CM to nearest macrophage (any source)
nn_far <- nn2(all_mac_coords, cm_coords, k = 1)
dist_to_any_mac <- as.numeric(nn_far$nn.dists)

# --- Assign groups ---
near_cm_cells <- cm_cells[dist_to_mac_res <= near_cutoff]
far_cm_cells  <- cm_cells[dist_to_any_mac > far_cutoff]

cat("Near CMs (<=", near_cutoff, "µm from Mac_Resident):", length(near_cm_cells), "\n")
cat("Far CMs (>", far_cutoff, "µm from any macrophage):", length(far_cm_cells), "\n")

# Per-patient breakdown
near_patients <- table(xenium.obj@meta.data[near_cm_cells, "patient"])
far_patients  <- table(xenium.obj@meta.data[far_cm_cells, "patient"])
cat("Near CMs per patient:\n"); print(near_patients)
cat("Far CMs per patient:\n");  print(far_patients)

# --- Annotate and subset ---
xenium.obj@meta.data$cm_mac_proximity <- NA_character_
xenium.obj@meta.data[near_cm_cells, "cm_mac_proximity"] <- "near_macrophage"
xenium.obj@meta.data[far_cm_cells, "cm_mac_proximity"]  <- "far_from_macrophage"

cm_prox <- subset(xenium.obj, cells = c(near_cm_cells, far_cm_cells))

# --- Pseudobulk DESeq2: near vs far, controlling for patient ---
cm_prox$pseudobulk_id <- paste0(
  as.character(cm_prox$patient), '..',
  as.character(cm_prox$cm_mac_proximity)
)

prox_meta <- unique(cm_prox@meta.data[, c('patient', 'group', 'cm_mac_proximity')])
prox_meta$patient <- as.character(prox_meta$patient)
prox_meta$pseudobulk_id <- paste0(prox_meta$patient, '..', prox_meta$cm_mac_proximity)
prox_meta$group <- factor(prox_meta$group, levels = c('NF', 'pRV', 'RVF'))
prox_meta$cm_mac_proximity <- factor(prox_meta$cm_mac_proximity,
  levels = c('far_from_macrophage', 'near_macrophage'))
rownames(prox_meta) <- prox_meta$pseudobulk_id

prox_counts <- AggregateExpression(
  cm_prox,
  assays = 'Xenium',
  group.by = 'pseudobulk_id',
  slot = 'counts',
  return.seurat = FALSE
)$Xenium

# Normalize column names (AggregateExpression quirks)
colnames(prox_counts) <- gsub('-', '_', sub('^g', '', colnames(prox_counts)))

# Filter to samples present in both counts and metadata
common_samples <- intersect(colnames(prox_counts), rownames(prox_meta))
prox_counts <- prox_counts[, common_samples, drop = FALSE]
prox_meta   <- prox_meta[common_samples, , drop = FALSE]

# --- Suppress myeloid-specific genes ---
# Load marker results from downsampled wilcoxauc run
all.marks <- read.csv('./output/Xenium/all.marks.csv')

# Myeloid cell types (cell_type_seurat labels)
myeloid_types <- c('Myeloid')

# CM types: any group starting with 'CM'
cm_types <- unique(all.marks$group[grepl('^CM', all.marks$group)])

# Myeloid marker genes: significant in any myeloid cell type (padj < 0.05, logFC > 0.5, auc > 0.6)
myeloid_marker_genes <- all.marks %>%
  filter(group %in% myeloid_types, padj < 0.05, logFC > 0.5, auc > 0.6) %>%
  pull(feature) %>% unique()

# CM marker genes: significant in any CM subtype
cm_marker_genes <- all.marks %>%
  filter(group %in% cm_types, padj < 0.05, logFC > 0.5, auc > 0.6) %>%
  pull(feature) %>% unique()

# Exclude genes that are myeloid markers but NOT CM markers
genes_to_exclude <- setdiff(myeloid_marker_genes, cm_marker_genes)
cat("Myeloid marker genes:", length(myeloid_marker_genes),
    "| CM marker genes:", length(cm_marker_genes),
    "| Genes to exclude:", length(genes_to_exclude), "\n")

prox_counts <- prox_counts[!rownames(prox_counts) %in% genes_to_exclude, , drop = FALSE]

# Filter low-count genes
keep <- rowSums(prox_counts >= 5) >= 2
prox_counts <- prox_counts[keep, , drop = FALSE]

cat("Pseudobulk samples:", ncol(prox_counts), "| Genes after filtering:", nrow(prox_counts), "\n")

# (a) Pooled analysis: ~ patient + cm_mac_proximity
dds_prox <- DESeqDataSetFromMatrix(
  countData = prox_counts,
  colData   = prox_meta,
  design    = ~ patient + cm_mac_proximity
)
dds_prox <- DESeq(dds_prox)

res_prox <- results(dds_prox, contrast = c('cm_mac_proximity', 'near_macrophage', 'far_from_macrophage')) %>%
  as.data.frame() %>%
  tibble::rownames_to_column('gene') %>%
  mutate(comparison = 'near_vs_far_pooled') %>%
  arrange(padj)

# (b) Per-group analysis
prox_deseq_pergroup <- list()
for (grp in c('NF', 'pRV', 'RVF')) {
  grp_samples <- prox_meta$pseudobulk_id[prox_meta$group == grp]
  grp_samples <- grp_samples[grp_samples %in% colnames(prox_counts)]
  if (length(grp_samples) < 2) next

  grp_counts <- prox_counts[, grp_samples, drop = FALSE]
  grp_meta   <- prox_meta[grp_samples, , drop = FALSE]

  if (length(unique(grp_meta$cm_mac_proximity)) < 2) next

  keep_grp <- rowSums(grp_counts >= 5) >= 2
  grp_counts <- grp_counts[keep_grp, , drop = FALSE]
  if (nrow(grp_counts) < 10) next

  dds_grp <- tryCatch({
    dds_tmp <- DESeqDataSetFromMatrix(
      countData = grp_counts,
      colData   = grp_meta,
      design    = ~ cm_mac_proximity
    )
    DESeq(dds_tmp)
  }, error = function(e) NULL)
  if (is.null(dds_grp)) next

  res_grp <- tryCatch(
    results(dds_grp, contrast = c('cm_mac_proximity', 'near_macrophage', 'far_from_macrophage')) %>%
      as.data.frame() %>%
      tibble::rownames_to_column('gene') %>%
      mutate(comparison = paste0('near_vs_far_', grp)) %>%
      arrange(padj),
    error = function(e) NULL
  )
  if (!is.null(res_grp)) prox_deseq_pergroup[[grp]] <- res_grp
}

# Combine all results and save
all_prox_deseq <- bind_rows(c(list(pooled = res_prox), prox_deseq_pergroup))
write.csv(all_prox_deseq, './output/Xenium/cm_proximity_deseq2.csv', row.names = FALSE)

sig_prox <- all_prox_deseq %>% filter(padj < 0.05, abs(log2FoldChange) > 0.5)
cat("Significant DEGs (pooled):", nrow(sig_prox %>% filter(comparison == 'near_vs_far_pooled')), "\n")

sig_prox <- all_prox_deseq %>% filter(padj < 0.05, abs(log2FoldChange) > 0.5)
cat("Significant DEGs (pRV):", nrow(sig_prox %>% filter(comparison == 'near_vs_far_pRV')), "\n")

sig_prox <- all_prox_deseq %>% filter(padj < 0.05, abs(log2FoldChange) > 0.5)
cat("Significant DEGs (RVF):", nrow(sig_prox %>% filter(comparison == 'near_vs_far_RVF')), "\n")



# --- Volcano plot: pooled near vs far ---
res_prox_capped <- res_prox %>% cap_pval()


volcano_args <- list(
  x = 'log2FoldChange', y = 'padj',
  pCutoff = 0.05, FCcutoff = 0.5,
  pointSize = 2, labSize = 3.5,
  drawConnectors = TRUE, widthConnectors = 0.5,
  xlim = c(-1, 1), ylim = c(0, 5),
  subtitle = 'Red = original panel gene | Blue = imputed gene',
  legendPosition = 'none'
)

pdf('./output/Xenium/cm_proximity_volcano.pdf', width = 10, height = 8)
do.call(EnhancedVolcano, c(list(
  toptable = res_prox_capped,
  lab = res_prox_capped$gene,
  title = 'CM near Macrophage_Resident vs far from macrophage',
  selectLab = sig_labs(res_prox_capped,100),
  colCustom = make_volcano_colors(res_prox_capped, panel_genes)),
  volcano_args))
dev.off()

# --- EnrichR: pooled near vs far ---
sig_pooled <- res_prox %>% filter(padj < 0.05, abs(log2FoldChange) > 0.5)
up_genes_prox   <- sig_pooled %>% filter(log2FoldChange > 0) %>% pull(gene)
down_genes_prox <- sig_pooled %>% filter(log2FoldChange < 0) %>% pull(gene)

enrich_prox_up   <- if (length(up_genes_prox) > 0)   enrichr(up_genes_prox, enrichR_dbs)   else NULL
enrich_prox_down <- if (length(down_genes_prox) > 0) enrichr(down_genes_prox, enrichR_dbs) else NULL

all_prox_enrich <- bind_rows(
  combine_enrich(enrich_prox_up,   'up',   'near_vs_far'),
  combine_enrich(enrich_prox_down, 'down', 'near_vs_far')
)
write.csv(all_prox_enrich, './output/Xenium/cm_proximity_enrichr.csv', row.names = FALSE)

pdf('./output/Xenium/cm_proximity_enrichr_plots.pdf', width = 12, height = 7)
if (!is.null(enrich_prox_up))   print(plot_enrich(enrich_prox_up$GO_Biological_Process_2023,   'CM near Mac - UP GO BP'))
if (!is.null(enrich_prox_up))   print(plot_enrich(enrich_prox_up$Reactome_2022,                'CM near Mac - UP Reactome'))
if (!is.null(enrich_prox_down)) print(plot_enrich(enrich_prox_down$GO_Biological_Process_2023, 'CM near Mac - DOWN GO BP'))
if (!is.null(enrich_prox_down)) print(plot_enrich(enrich_prox_down$Reactome_2022,              'CM near Mac - DOWN Reactome'))
dev.off()

xenium.obj$cell_type_seurat[far_cm_cells]
xenium.obj$cell_type_seurat[near_cm_cells]

##############################################
#### CM proximity to iMac: pseudobulk DESeq2
##############################################

library(RANN)




# install.packages('presto')  # run once if not installed
library(presto)

# Downsample each cell_type_rctd_doublet to at most 5000 cells, then densify and run wilcoxauc
Idents(xenium.obj) <- 'cell_type_rctd_doublet'
xenium.ds <- subset(xenium.obj, downsample = 5000)
mat_ds <- as(GetAssayData(xenium.ds, assay = 'Xenium', layer = 'data'), 'dgCMatrix')
all.marks <- wilcoxauc(mat_ds, y = xenium.ds$cell_type_rctd_doublet)
rm(xenium.ds, mat_ds); gc()
write.csv(all.marks, './output/Xenium/all.marks.csv')


near_cutoff <- 15   # µm — CMs within this distance of a Macrophage_Resident
far_cutoff  <- 50   # µm — CMs beyond this distance from ANY macrophage

# --- Extract CM coordinates from xenium.obj ---
cm_cells <- WhichCells(xenium.obj, expression = cell_type_rctd_doublet == "CM")
cm_cells_seurat <- WhichCells(xenium.obj, expression = cell_type_seurat %in% c('Cm1','Cm2','Cm3','Cm4','Cm5','Cm6','Cm7','Cm8','Cm9'))
cm_cells <- intersect(cm_cells,cm_cells_seurat)

cm_coords <- as.matrix(xenium.obj@meta.data[cm_cells, c("x_centroid", "y_centroid")])

broad_mac_cells <- WhichCells(xenium.obj,
  expression = cell_type_seurat %in% c("iMac"))

all_mac_barcodes <- broad_mac_cells#unique(c(broad_mac_cells, mac_res_cells))
all_mac_barcodes <- all_mac_barcodes[all_mac_barcodes %in% rownames(xenium.obj@meta.data)]
all_mac_coords <- as.matrix(xenium.obj@meta.data[all_mac_barcodes, c("x_centroid", "y_centroid")])
mac_res_coords <- all_mac_coords

cat("CMs:", nrow(cm_coords), "| iMac (refined):", nrow(mac_res_coords),
    "| iMac (union):", nrow(all_mac_coords), "\n")

# --- Compute nearest-neighbor distances ---
# For "near" group: distance from each CM to nearest Macrophage_Resident
nn_near <- nn2(mac_res_coords, cm_coords, k = 1)
dist_to_mac_res <- as.numeric(nn_near$nn.dists)

# For "far" group: distance from each CM to nearest macrophage (any source)
nn_far <- nn2(all_mac_coords, cm_coords, k = 1)
dist_to_any_mac <- as.numeric(nn_far$nn.dists)

# --- Assign groups ---
near_cm_cells <- cm_cells[dist_to_mac_res <= near_cutoff]
far_cm_cells  <- cm_cells[dist_to_any_mac > far_cutoff]

cat("Near CMs (<=", near_cutoff, "µm from iMac):", length(near_cm_cells), "\n")
cat("Far CMs (>", far_cutoff, "µm from any iMac):", length(far_cm_cells), "\n")

# Per-patient breakdown
near_patients <- table(xenium.obj@meta.data[near_cm_cells, "patient"])
far_patients  <- table(xenium.obj@meta.data[far_cm_cells, "patient"])
cat("Near CMs per patient:\n"); print(near_patients)
cat("Far CMs per patient:\n");  print(far_patients)

# --- Annotate and subset ---
xenium.obj@meta.data$cm_mac_proximity <- NA_character_
xenium.obj@meta.data[near_cm_cells, "cm_mac_proximity"] <- "near_macrophage"
xenium.obj@meta.data[far_cm_cells, "cm_mac_proximity"]  <- "far_from_macrophage"

cm_prox <- subset(xenium.obj, cells = c(near_cm_cells, far_cm_cells))

# --- Pseudobulk DESeq2: near vs far, controlling for patient ---
cm_prox$pseudobulk_id <- paste0(
  as.character(cm_prox$patient), '..',
  as.character(cm_prox$cm_mac_proximity)
)

prox_meta <- unique(cm_prox@meta.data[, c('patient', 'group', 'cm_mac_proximity')])
prox_meta$patient <- as.character(prox_meta$patient)
prox_meta$pseudobulk_id <- paste0(prox_meta$patient, '..', prox_meta$cm_mac_proximity)
prox_meta$group <- factor(prox_meta$group, levels = c('NF', 'pRV', 'RVF'))
prox_meta$cm_mac_proximity <- factor(prox_meta$cm_mac_proximity,
  levels = c('far_from_macrophage', 'near_macrophage'))
rownames(prox_meta) <- prox_meta$pseudobulk_id

prox_counts <- AggregateExpression(
  cm_prox,
  assays = 'Xenium',
  group.by = 'pseudobulk_id',
  slot = 'counts',
  return.seurat = FALSE
)$Xenium

# Normalize column names (AggregateExpression quirks)
colnames(prox_counts) <- gsub('-', '_', sub('^g', '', colnames(prox_counts)))

# Filter to samples present in both counts and metadata
common_samples <- intersect(colnames(prox_counts), rownames(prox_meta))
prox_counts <- prox_counts[, common_samples, drop = FALSE]
prox_meta   <- prox_meta[common_samples, , drop = FALSE]

# --- Suppress myeloid-specific genes ---
# Load marker results from downsampled wilcoxauc run
all.marks <- read.csv('./output/Xenium/all.marks.csv')

# Myeloid cell types (cell_type_seurat labels)
myeloid_types <- c('Myeloid')

# CM types: any group starting with 'CM'
cm_types <- unique(all.marks$group[grepl('^CM', all.marks$group)])

# Myeloid marker genes: significant in any myeloid cell type (padj < 0.05, logFC > 0.5, auc > 0.6)
myeloid_marker_genes <- all.marks %>%
  filter(group %in% myeloid_types, padj < 0.05, logFC > 0.5, auc > 0.6) %>%
  pull(feature) %>% unique()

# CM marker genes: significant in any CM subtype
cm_marker_genes <- all.marks %>%
  filter(group %in% cm_types, padj < 0.05, logFC > 0.5, auc > 0.6) %>%
  pull(feature) %>% unique()

# Exclude genes that are myeloid markers but NOT CM markers
genes_to_exclude <- setdiff(myeloid_marker_genes, cm_marker_genes)
cat("Myeloid marker genes:", length(myeloid_marker_genes),
    "| CM marker genes:", length(cm_marker_genes),
    "| Genes to exclude:", length(genes_to_exclude), "\n")

prox_counts <- prox_counts[!rownames(prox_counts) %in% genes_to_exclude, , drop = FALSE]

# Filter low-count genes
keep <- rowSums(prox_counts >= 5) >= 2
prox_counts <- prox_counts[keep, , drop = FALSE]

cat("Pseudobulk samples:", ncol(prox_counts), "| Genes after filtering:", nrow(prox_counts), "\n")

# (a) Pooled analysis: ~ patient + cm_mac_proximity
dds_prox <- DESeqDataSetFromMatrix(
  countData = prox_counts,
  colData   = prox_meta,
  design    = ~ patient + cm_mac_proximity
)
dds_prox <- DESeq(dds_prox)

res_prox <- results(dds_prox, contrast = c('cm_mac_proximity', 'near_macrophage', 'far_from_macrophage')) %>%
  as.data.frame() %>%
  tibble::rownames_to_column('gene') %>%
  mutate(comparison = 'near_vs_far_pooled') %>%
  arrange(padj)

# (b) Per-group analysis
prox_deseq_pergroup <- list()
for (grp in c('NF', 'pRV', 'RVF')) {
  grp_samples <- prox_meta$pseudobulk_id[prox_meta$group == grp]
  grp_samples <- grp_samples[grp_samples %in% colnames(prox_counts)]
  if (length(grp_samples) < 2) next

  grp_counts <- prox_counts[, grp_samples, drop = FALSE]
  grp_meta   <- prox_meta[grp_samples, , drop = FALSE]

  if (length(unique(grp_meta$cm_mac_proximity)) < 2) next

  keep_grp <- rowSums(grp_counts >= 5) >= 2
  grp_counts <- grp_counts[keep_grp, , drop = FALSE]
  if (nrow(grp_counts) < 10) next

  dds_grp <- tryCatch({
    dds_tmp <- DESeqDataSetFromMatrix(
      countData = grp_counts,
      colData   = grp_meta,
      design    = ~ cm_mac_proximity
    )
    DESeq(dds_tmp)
  }, error = function(e) NULL)
  if (is.null(dds_grp)) next

  res_grp <- tryCatch(
    results(dds_grp, contrast = c('cm_mac_proximity', 'near_macrophage', 'far_from_macrophage')) %>%
      as.data.frame() %>%
      tibble::rownames_to_column('gene') %>%
      mutate(comparison = paste0('near_vs_far_', grp)) %>%
      arrange(padj),
    error = function(e) NULL
  )
  if (!is.null(res_grp)) prox_deseq_pergroup[[grp]] <- res_grp
}

# Combine all results and save
all_prox_deseq <- bind_rows(c(list(pooled = res_prox), prox_deseq_pergroup))
write.csv(all_prox_deseq, './output/Xenium/cm_proximity_imac_deseq2.csv', row.names = FALSE)

sig_prox <- all_prox_deseq %>% filter(padj < 0.05, abs(log2FoldChange) > 0.5)
cat("Significant DEGs (pooled):", nrow(sig_prox %>% filter(comparison == 'near_vs_far_pooled')), "\n")

# --- Volcano plot: pooled near vs far ---
res_prox_capped <- res_prox %>% cap_pval()


volcano_args <- list(
  x = 'log2FoldChange', y = 'padj',
  pCutoff = 0.05, FCcutoff = 0.5,
  pointSize = 2, labSize = 3.5,
  drawConnectors = TRUE, widthConnectors = 0.5,
  xlim = c(-1, 1), ylim = c(0, 5),
  subtitle = 'Red = original panel gene | Blue = imputed gene',
  legendPosition = 'none'
)

pdf('./output/Xenium/cm_proximity_imac_volcano.pdf', width = 10, height = 8)
do.call(EnhancedVolcano, c(list(
  toptable = res_prox_capped,
  lab = res_prox_capped$gene,
  title = 'CM near iMac vs far from macrophage',
  selectLab = sig_labs(res_prox_capped,100),
  colCustom = make_volcano_colors(res_prox_capped, panel_genes)),
  volcano_args))
dev.off()

# --- EnrichR: pooled near vs far ---
sig_pooled <- res_prox %>% filter(padj < 0.05, abs(log2FoldChange) > 0.5)
up_genes_prox   <- sig_pooled %>% filter(log2FoldChange > 0) %>% pull(gene)
down_genes_prox <- sig_pooled %>% filter(log2FoldChange < 0) %>% pull(gene)

enrich_prox_up   <- if (length(up_genes_prox) > 0)   enrichr(up_genes_prox, enrichR_dbs)   else NULL
enrich_prox_down <- if (length(down_genes_prox) > 0) enrichr(down_genes_prox, enrichR_dbs) else NULL

all_prox_enrich <- bind_rows(
  combine_enrich(enrich_prox_up,   'up',   'near_vs_far'),
  combine_enrich(enrich_prox_down, 'down', 'near_vs_far')
)
write.csv(all_prox_enrich, './output/Xenium/cm_proximity_imac_enrichr.csv', row.names = FALSE)

pdf('./output/Xenium/cm_proximity_enrichr_plots.pdf', width = 12, height = 7)
if (!is.null(enrich_prox_up))   print(plot_enrich(enrich_prox_up$GO_Biological_Process_2023,   'CM near Mac - UP GO BP'))
if (!is.null(enrich_prox_up))   print(plot_enrich(enrich_prox_up$Reactome_2022,                'CM near Mac - UP Reactome'))
if (!is.null(enrich_prox_down)) print(plot_enrich(enrich_prox_down$GO_Biological_Process_2023, 'CM near Mac - DOWN GO BP'))
if (!is.null(enrich_prox_down)) print(plot_enrich(enrich_prox_down$Reactome_2022,              'CM near Mac - DOWN Reactome'))
dev.off()

xenium.obj$cell_type_seurat[far_cm_cells]
xenium.obj$cell_type_seurat[near_cm_cells]


##############################################
#### CM proximity to capEC: pseudobulk DESeq2
##############################################

library(RANN)




# install.packages('presto')  # run once if not installed
library(presto)

# Downsample each cell_type_rctd_doublet to at most 5000 cells, then densify and run wilcoxauc
Idents(xenium.obj) <- 'cell_type_rctd_doublet'
xenium.ds <- subset(xenium.obj, downsample = 5000)
mat_ds <- as(GetAssayData(xenium.ds, assay = 'Xenium', layer = 'data'), 'dgCMatrix')
all.marks <- wilcoxauc(mat_ds, y = xenium.ds$cell_type_rctd_doublet)
rm(xenium.ds, mat_ds); gc()
write.csv(all.marks, './output/Xenium/all.marks.csv')


near_cutoff <- 15   # µm — CMs within this distance of a Macrophage_Resident
far_cutoff  <- 50   # µm — CMs beyond this distance from ANY macrophage

# --- Extract CM coordinates from xenium.obj ---
cm_cells <- WhichCells(xenium.obj, expression = cell_type_rctd_doublet == "CM")
cm_cells_seurat <- WhichCells(xenium.obj, expression = cell_type_seurat %in% c('Cm1','Cm2','Cm3','Cm4','Cm5','Cm6','Cm7','Cm8','Cm9'))
cm_cells <- intersect(cm_cells,cm_cells_seurat)

cm_coords <- as.matrix(xenium.obj@meta.data[cm_cells, c("x_centroid", "y_centroid")])

broad_mac_cells <- WhichCells(xenium.obj,
  expression = cell_type_seurat %in% c("EC_Capillary"))

all_mac_barcodes <- broad_mac_cells#unique(c(broad_mac_cells, mac_res_cells))
all_mac_barcodes <- all_mac_barcodes[all_mac_barcodes %in% rownames(xenium.obj@meta.data)]
all_mac_coords <- as.matrix(xenium.obj@meta.data[all_mac_barcodes, c("x_centroid", "y_centroid")])
mac_res_coords <- all_mac_coords

cat("CMs:", nrow(cm_coords), "| cap EC (refined):", nrow(mac_res_coords),
    "| cap EC (union):", nrow(all_mac_coords), "\n")

# --- Compute nearest-neighbor distances ---
# For "near" group: distance from each CM to nearest Macrophage_Resident
nn_near <- nn2(mac_res_coords, cm_coords, k = 1)
dist_to_mac_res <- as.numeric(nn_near$nn.dists)

# For "far" group: distance from each CM to nearest macrophage (any source)
nn_far <- nn2(all_mac_coords, cm_coords, k = 1)
dist_to_any_mac <- as.numeric(nn_far$nn.dists)

# --- Assign groups ---
near_cm_cells <- cm_cells[dist_to_mac_res <= near_cutoff]
far_cm_cells  <- cm_cells[dist_to_any_mac > far_cutoff]

cat("Near CMs (<=", near_cutoff, "µm from cap EC):", length(near_cm_cells), "\n")
cat("Far CMs (>", far_cutoff, "µm from any cap EC):", length(far_cm_cells), "\n")

# Per-patient breakdown
near_patients <- table(xenium.obj@meta.data[near_cm_cells, "patient"])
far_patients  <- table(xenium.obj@meta.data[far_cm_cells, "patient"])
cat("Near CMs per patient:\n"); print(near_patients)
cat("Far CMs per patient:\n");  print(far_patients)

# --- Annotate and subset ---
xenium.obj@meta.data$cm_mac_proximity <- NA_character_
xenium.obj@meta.data[near_cm_cells, "cm_ec_proximity"] <- "near_cap_ec"
xenium.obj@meta.data[far_cm_cells, "cm_ec_proximity"]  <- "far_from_cap_ec"

cm_prox <- subset(xenium.obj, cells = c(near_cm_cells, far_cm_cells))

# --- Pseudobulk DESeq2: near vs far, controlling for patient ---
cm_prox$pseudobulk_id <- paste0(
  as.character(cm_prox$patient), '..',
  as.character(cm_prox$cm_ec_proximity)
)

prox_meta <- unique(cm_prox@meta.data[, c('patient', 'group', 'cm_ec_proximity')])
prox_meta$patient <- as.character(prox_meta$patient)
prox_meta$pseudobulk_id <- paste0(prox_meta$patient, '..', prox_meta$cm_ec_proximity)
prox_meta$group <- factor(prox_meta$group, levels = c('NF', 'pRV', 'RVF'))
prox_meta$cm_mac_proximity <- factor(prox_meta$cm_ec_proximity,
  levels = c('far_from_cap_ec', 'near_cap_ec'))
rownames(prox_meta) <- prox_meta$pseudobulk_id

prox_counts <- AggregateExpression(
  cm_prox,
  assays = 'Xenium',
  group.by = 'pseudobulk_id',
  slot = 'counts',
  return.seurat = FALSE
)$Xenium

# Normalize column names (AggregateExpression quirks)
colnames(prox_counts) <- gsub('-', '_', sub('^g', '', colnames(prox_counts)))

# Filter to samples present in both counts and metadata
common_samples <- intersect(colnames(prox_counts), rownames(prox_meta))
prox_counts <- prox_counts[, common_samples, drop = FALSE]
prox_meta   <- prox_meta[common_samples, , drop = FALSE]

# --- Suppress myeloid-specific genes ---
# Load marker results from downsampled wilcoxauc run
all.marks <- read.csv('./output/Xenium/all.marks.csv')

# EC cell types (cell_type_seurat labels)
ec_types <- c('EC')

# CM types: any group starting with 'CM'
cm_types <- unique(all.marks$group[grepl('^CM', all.marks$group)])

# EC marker genes: significant in any myeloid cell type (padj < 0.05, logFC > 0.5, auc > 0.6)
ec_marker_genes <- all.marks %>%
  filter(group %in% ec_types, padj < 0.05, logFC > 0.5, auc > 0.6) %>%
  pull(feature) %>% unique()

# CM marker genes: significant in any CM subtype
cm_marker_genes <- all.marks %>%
  filter(group %in% cm_types, padj < 0.05, logFC > 0.5, auc > 0.6) %>%
  pull(feature) %>% unique()

# Exclude genes that are myeloid markers but NOT CM markers
genes_to_exclude <- setdiff(ec_marker_genes, cm_marker_genes)
cat("EC marker genes:", length(ec_marker_genes),
    "| CM marker genes:", length(cm_marker_genes),
    "| Genes to exclude:", length(genes_to_exclude), "\n")

prox_counts <- prox_counts[!rownames(prox_counts) %in% genes_to_exclude, , drop = FALSE]

# Filter low-count genes
keep <- rowSums(prox_counts >= 5) >= 2
prox_counts <- prox_counts[keep, , drop = FALSE]

cat("Pseudobulk samples:", ncol(prox_counts), "| Genes after filtering:", nrow(prox_counts), "\n")

# (a) Pooled analysis: ~ patient + cm_mac_proximity
dds_prox <- DESeqDataSetFromMatrix(
  countData = prox_counts,
  colData   = prox_meta,
  design    = ~ patient + cm_ec_proximity
)
dds_prox <- DESeq(dds_prox)

res_prox <- results(dds_prox, contrast = c('cm_ec_proximity', 'near_cap_ec', 'far_from_cap_ec')) %>%
  as.data.frame() %>%
  tibble::rownames_to_column('gene') %>%
  mutate(comparison = 'near_vs_far_pooled') %>%
  arrange(padj)

# (b) Per-group analysis
prox_deseq_pergroup <- list()
for (grp in c('NF', 'pRV', 'RVF')) {
  grp_samples <- prox_meta$pseudobulk_id[prox_meta$group == grp]
  grp_samples <- grp_samples[grp_samples %in% colnames(prox_counts)]
  if (length(grp_samples) < 2) next

  grp_counts <- prox_counts[, grp_samples, drop = FALSE]
  grp_meta   <- prox_meta[grp_samples, , drop = FALSE]

  if (length(unique(grp_meta$cm_ec_proximity)) < 2) next

  keep_grp <- rowSums(grp_counts >= 5) >= 2
  grp_counts <- grp_counts[keep_grp, , drop = FALSE]
  if (nrow(grp_counts) < 10) next

  dds_grp <- tryCatch({
    dds_tmp <- DESeqDataSetFromMatrix(
      countData = grp_counts,
      colData   = grp_meta,
      design    = ~ cm_ec_proximity
    )
    DESeq(dds_tmp)
  }, error = function(e) NULL)
  if (is.null(dds_grp)) next

  res_grp <- tryCatch(
    results(dds_grp, contrast = c('cm_ec_proximity', 'near_cap_ec', 'far_from_cap_ec')) %>%
      as.data.frame() %>%
      tibble::rownames_to_column('gene') %>%
      mutate(comparison = paste0('near_vs_far_', grp)) %>%
      arrange(padj),
    error = function(e) NULL
  )
  if (!is.null(res_grp)) prox_deseq_pergroup[[grp]] <- res_grp
}

# Combine all results and save
all_prox_deseq <- bind_rows(c(list(pooled = res_prox), prox_deseq_pergroup))
write.csv(all_prox_deseq, './output/Xenium/cm_proximity_ec_deseq2.csv', row.names = FALSE)

sig_prox <- all_prox_deseq %>% filter(padj < 0.05, abs(log2FoldChange) > 0.5)
cat("Significant DEGs (pooled):", nrow(sig_prox %>% filter(comparison == 'near_vs_far_pooled')), "\n")

sig_prox <- all_prox_deseq %>% filter(padj < 0.05, abs(log2FoldChange) > 0.5)
cat("Significant DEGs (NF):", nrow(sig_prox %>% filter(comparison == 'near_vs_far_NF')), "\n")

sig_prox <- all_prox_deseq %>% filter(padj < 0.05, abs(log2FoldChange) > 0.5)
cat("Significant DEGs (pRV):", nrow(sig_prox %>% filter(comparison == 'near_vs_far_pRV')), "\n")

sig_prox <- all_prox_deseq %>% filter(padj < 0.05, abs(log2FoldChange) > 0.5)
cat("Significant DEGs (RVF):", nrow(sig_prox %>% filter(comparison == 'near_vs_far_RVF')), "\n")


# --- Volcano plot: pooled near vs far ---
res_prox_capped <- res_prox %>% cap_pval()


volcano_args <- list(
  x = 'log2FoldChange', y = 'padj',
  pCutoff = 0.05, FCcutoff = 0.5,
  pointSize = 2, labSize = 3.5,
  drawConnectors = TRUE, widthConnectors = 0.5,
  xlim = c(-1, 1), ylim = c(0, 5),
  subtitle = 'Red = original panel gene | Blue = imputed gene',
  legendPosition = 'none'
)

pdf('./output/Xenium/cm_proximity_ec_volcano.pdf', width = 10, height = 8)
do.call(EnhancedVolcano, c(list(
  toptable = res_prox_capped,
  lab = res_prox_capped$gene,
  title = 'CM near iMac vs far from macrophage',
  selectLab = sig_labs(res_prox_capped,100),
  colCustom = make_volcano_colors(res_prox_capped, panel_genes)),
  volcano_args))
dev.off()

# --- EnrichR: pooled near vs far ---
sig_pooled <- res_prox %>% filter(padj < 0.05, abs(log2FoldChange) > 0.5)
up_genes_prox   <- sig_pooled %>% filter(log2FoldChange > 0) %>% pull(gene)
down_genes_prox <- sig_pooled %>% filter(log2FoldChange < 0) %>% pull(gene)

enrich_prox_up   <- if (length(up_genes_prox) > 0)   enrichr(up_genes_prox, enrichR_dbs)   else NULL
enrich_prox_down <- if (length(down_genes_prox) > 0) enrichr(down_genes_prox, enrichR_dbs) else NULL

all_prox_enrich <- bind_rows(
  combine_enrich(enrich_prox_up,   'up',   'near_vs_far'),
  combine_enrich(enrich_prox_down, 'down', 'near_vs_far')
)
write.csv(all_prox_enrich, './output/Xenium/cm_proximity_ec_enrichr.csv', row.names = FALSE)

pdf('./output/Xenium/cm_proximity_ec_enrichr_plots.pdf', width = 12, height = 7)
if (!is.null(enrich_prox_up))   print(plot_enrich(enrich_prox_up$GO_Biological_Process_2023,   'CM near Mac - UP GO BP'))
if (!is.null(enrich_prox_up))   print(plot_enrich(enrich_prox_up$Reactome_2022,                'CM near Mac - UP Reactome'))
if (!is.null(enrich_prox_down)) print(plot_enrich(enrich_prox_down$GO_Biological_Process_2023, 'CM near Mac - DOWN GO BP'))
if (!is.null(enrich_prox_down)) print(plot_enrich(enrich_prox_down$Reactome_2022,              'CM near Mac - DOWN Reactome'))
dev.off()

xenium.obj$cell_type_seurat[far_cm_cells]
xenium.obj$cell_type_seurat[near_cm_cells]











nkt.obj <- subset(xenium.obj, subset = cell_type_rctd_doublet == "NKT")
nkt.obj <- merge(nkt.obj, tnk.obj.clean)
nkt.obj <- JoinLayers(nkt.obj)

umap_cols <- grep("^umap_\\d+$|^UMAP_\\d+$|^umaporig_\\d+$",
                  colnames(nkt.obj@meta.data), value = TRUE)
if (length(umap_cols)) nkt.obj@meta.data[umap_cols] <- NULL

nkt.obj[["pca"]]     <- NULL
nkt.obj[["harmony"]]     <- NULL
nkt.obj[["umap"]]     <- NULL
nkt.obj[["umap_orig"]]     <- NULL

nkt.obj <- FindVariableFeatures(nkt.obj)
nkt.obj <- ScaleData(nkt.obj)
nkt.obj <- RunPCA(nkt.obj, npcs = 30)
nkt.obj <- RunHarmony(nkt.obj, "patient")

nkt.obj <- RunUMAP(nkt.obj, reduction = "harmony", dims = 1:30,reduction.name='umap')
nkt.obj <- RunUMAP(nkt.obj, reduction = "pca", dims = 1:30,reduction.name='umap_orig')

nkt.obj <- FindNeighbors(nkt.obj, reduction = "harmony", dims = 1:30)
nkt.obj <- FindClusters(nkt.obj, resolution = 0.5)

nkt.marks <- FindAllMarkers(nkt.obj)
write.csv(nkt.marks,'./output/Xenium/nkt.marks.csv')


DimPlot(nkt.obj, group.by = 'patient',reduction='umap_orig')

nkt_labels <- c(
  '0'  = 'CD8_T',
  '1'  = 'Fibroblast',
  '2'  = 'Cardiomyocyte',
  '3'  = 'CD4_T',
  '4'  = 'NK',
  '5'  = 'Endothelial',
  '6'  = 'Myeloid',
  '7'  = 'Neutrophil',
  '8'  = 'Unknown',
  '9'  = 'Pericyte',
  '10' = 'Pericyte',
  '11' = 'Basophil'
)

# True lymphoid clusters: 0 (CD8_T), 3 (CD4_T), 4 (NK)
# Clusters 1,2,5,6,7,9,10,11 are contamination
nkt.obj$nkt_subtype <- unname(nkt_labels[as.character(nkt.obj$seurat_clusters)])

DimPlot(nkt.obj, group.by = 'nkt_subtype', label = TRUE, repel = TRUE) +
  theme(legend.position = 'none')

nkt.obj.clean <- subset(nkt.obj, subset = nkt_subtype %in% c('CD8_T', 'CD4_T', 'NK'))

nkt.obj.clean[["pca"]]     <- NULL
nkt.obj.clean[["harmony"]] <- NULL
nkt.obj.clean[["umap"]]    <- NULL
nkt.obj.clean[["umap_orig"]] <- NULL

umap_cols <- grep("^umap_\\d+$|^UMAP_\\d+$|^umaporig_\\d+$",
                  colnames(nkt.obj.clean@meta.data), value = TRUE)
if (length(umap_cols)) nkt.obj.clean@meta.data[umap_cols] <- NULL

nkt.obj.clean <- JoinLayers(nkt.obj.clean)
nkt.obj.clean <- FindVariableFeatures(nkt.obj.clean)
nkt.obj.clean <- ScaleData(nkt.obj.clean)
nkt.obj.clean <- RunPCA(nkt.obj.clean, npcs = 30)
nkt.obj.clean <- RunHarmony(nkt.obj.clean, "patient")
nkt.obj.clean <- RunUMAP(nkt.obj.clean, reduction = "harmony", dims = 1:30, reduction.name = 'umap')
nkt.obj.clean <- RunUMAP(nkt.obj.clean, reduction = "pca",     dims = 1:30, reduction.name = 'umap_orig')

DimPlot(nkt.obj.clean, group.by = 'nkt_subtype', label = TRUE, repel = TRUE, reduction = 'umap_orig') +
  theme(legend.position = 'none')
DimPlot(nkt.obj.clean, group.by = 'patient', reduction = 'umap_orig')

saveRDS(nkt.obj.clean, './output/Xenium/nkt_clean.rds')

pdf(paste0('./output/Xenium/', 'Xenium_nkt_snUMAP_new.pdf'), width=5, height=5)
PlotEmbedding(nkt.obj.clean,group.by='nkt_subtype',reduction='umap_orig',point_size=.05,plot_under=TRUE,plot_theme=umap_theme()+NoLegend(),raster_dpi=400,raster_scale=0.5)
dev.off()

pdf(paste0('./output/Xenium/', 'Xenium_nkt_snUMAP_new_pt.pdf'), width=5, height=5)
PlotEmbedding(nkt.obj.clean,group.by='patient',reduction='umap_orig',point_size=.05,plot_under=TRUE,plot_theme=umap_theme()+NoLegend(),raster_dpi=400,raster_scale=0.5)
dev.off()

# NKT subtype proportions by disease state
nkt_meta <- nkt.obj.clean@meta.data %>%
  filter(!is.na(nkt_subtype)) %>%
  select(patient, group, nkt_subtype)

nkt_totals <- nkt_meta %>%
  group_by(patient, group) %>%
  summarise(total = n(), .groups = 'drop')

nkt_prop <- nkt_meta %>%
  group_by(patient, group, nkt_subtype) %>%
  summarise(n = n(), .groups = 'drop') %>%
  left_join(nkt_totals, by = c('patient', 'group')) %>%
  mutate(proportion = n / total) %>%
  complete(nesting(patient, group), nkt_subtype, fill = list(n = 0, proportion = 0)) %>%
  mutate(group = factor(group, levels = c('NF', 'pRV', 'RVF')))

pairwise_comparisons <- list(c('NF', 'pRV'), c('NF', 'RVF'), c('pRV', 'RVF'))

p_nkt_prop <- ggplot(nkt_prop, aes(x = group, y = proportion, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.5) +
  stat_compare_means(comparisons = pairwise_comparisons, method = 'wilcox.test', label = 'p.signif') +
  stat_compare_means(method = 'kruskal.test', label.y.npc = 0.95, size = 3) +
  facet_wrap(~ nkt_subtype, scales = 'free_y') +
  scale_x_discrete(limits = c('NF', 'pRV', 'RVF')) +
  theme_bw() +
  theme(legend.position = 'none',
        strip.text = element_text(face = 'bold'),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = 'Proportion of NKT cells', title = 'NKT subtype distribution by group')

pdf('./output/Xenium/nkt_proportions.pdf', width = 6, height = 4)
print(p_nkt_prop)
dev.off()


##############################################
#### Pseudobulk DESeq2 by NKT subtype
##############################################

# Build pseudobulk metadata
nkt_pseudo_meta <- unique(nkt.obj.clean@meta.data[, c('patient', 'group', 'nkt_subtype')])
nkt_pseudo_meta$patient     <- as.character(nkt_pseudo_meta$patient)
nkt_pseudo_meta$nkt_subtype <- as.character(nkt_pseudo_meta$nkt_subtype)
nkt_pseudo_meta$pseudobulk_id <- paste0(nkt_pseudo_meta$patient, '..', nkt_pseudo_meta$nkt_subtype)
nkt_pseudo_meta$group <- factor(nkt_pseudo_meta$group, levels = c('NF', 'pRV', 'RVF'))
rownames(nkt_pseudo_meta) <- nkt_pseudo_meta$pseudobulk_id

nkt.obj.clean$pseudobulk_id <- paste0(
  as.character(nkt.obj.clean$patient), '..',
  as.character(nkt.obj.clean$nkt_subtype)
)

nkt_pseudo_counts <- AggregateExpression(
  nkt.obj.clean,
  assays = 'Xenium',
  group.by = 'pseudobulk_id',
  slot = 'counts',
  return.seurat = FALSE
)$Xenium

colnames(nkt_pseudo_counts) <- gsub('-', '_', sub('^g', '', colnames(nkt_pseudo_counts)))

nkt_cell_types  <- unique(nkt_pseudo_meta$nkt_subtype)
nkt_deseq_results <- list()

for (ct in nkt_cell_types) {
  ct_samples <- nkt_pseudo_meta$pseudobulk_id[nkt_pseudo_meta$nkt_subtype == ct]
  ct_samples <- ct_samples[ct_samples %in% colnames(nkt_pseudo_counts)]
  if (length(ct_samples) < 2) next

  ct_counts <- nkt_pseudo_counts[, ct_samples, drop = FALSE]
  ct_meta   <- nkt_pseudo_meta[ct_samples, , drop = FALSE]

  if (length(unique(ct_meta$group)) < 2) next

  keep <- rowSums(ct_counts >= 5) >= 2
  ct_counts <- ct_counts[keep, , drop = FALSE]
  if (nrow(ct_counts) < 10) next

  dds <- DESeqDataSetFromMatrix(countData = ct_counts, colData = ct_meta, design = ~ group)
  dds <- tryCatch(DESeq(dds), error = function(e) NULL)
  if (is.null(dds)) next

  for (comp in list(c('RVF', 'NF'), c('pRV', 'NF'), c('RVF', 'pRV'))) {
    comp_label <- paste0(comp[1], '_vs_', comp[2])
    res <- tryCatch(
      results(dds, contrast = c('group', comp[1], comp[2])) %>%
        as.data.frame() %>%
        tibble::rownames_to_column('gene') %>%
        mutate(cell_type = ct, comparison = comp_label) %>%
        arrange(padj),
      error = function(e) NULL
    )
    if (!is.null(res)) nkt_deseq_results[[paste0(ct, '_', comp_label)]] <- res
  }
}

nkt_all_deseq <- bind_rows(nkt_deseq_results)
write.csv(nkt_all_deseq, './output/Xenium/nkt_pseudobulk_deseq2.csv', row.names = FALSE)

nkt_sig_summary <- nkt_all_deseq %>%
  filter(padj < 0.05, abs(log2FoldChange) > 0.5) %>%
  group_by(cell_type, comparison) %>%
  summarise(n_up = sum(log2FoldChange > 0), n_down = sum(log2FoldChange < 0), .groups = 'drop')
print(nkt_sig_summary)


##############################################
#### NKT subtypes: Volcano + GO/Reactome
##############################################

nkt_enrich_combined <- list()

for (ct in unique(nkt_all_deseq$cell_type)) {
  ct_slug <- gsub('_', '-', ct)
  cat('NKT processing:', ct, '\n')

  dat_rvf     <- nkt_all_deseq %>% filter(cell_type == ct, comparison == 'RVF_vs_NF')  %>% cap_pval()
  dat_prv     <- nkt_all_deseq %>% filter(cell_type == ct, comparison == 'pRV_vs_NF')  %>% cap_pval()
  dat_rvf_prv <- nkt_all_deseq %>% filter(cell_type == ct, comparison == 'RVF_vs_pRV') %>% cap_pval()

  make_volcano_pdf(dat_rvf,     paste0('./output/Xenium/nkt_', ct_slug, '_volcano_RVF_vs_NF.pdf'),  paste(ct, ': RVF vs NF'))
  make_volcano_pdf(dat_prv,     paste0('./output/Xenium/nkt_', ct_slug, '_volcano_pRV_vs_NF.pdf'),  paste(ct, ': pRV vs NF'))
  make_volcano_pdf(dat_rvf_prv, paste0('./output/Xenium/nkt_', ct_slug, '_volcano_RVF_vs_pRV.pdf'), paste(ct, ': RVF vs pRV'))

  run_enrich <- function(dat, comp) {
    sig  <- dat %>% filter(padj < 0.05, abs(log2FoldChange) > 0.5)
    up   <- sig %>% filter(log2FoldChange > 0) %>% pull(gene)
    down <- sig %>% filter(log2FoldChange < 0) %>% pull(gene)
    list(
      up   = if (length(up)   > 0) enrichr(up,   enrichR_dbs) else NULL,
      down = if (length(down) > 0) enrichr(down, enrichR_dbs) else NULL
    )
  }

  er_rvf     <- run_enrich(dat_rvf,     'RVF_vs_NF')
  er_prv     <- run_enrich(dat_prv,     'pRV_vs_NF')
  er_rvf_prv <- run_enrich(dat_rvf_prv, 'RVF_vs_pRV')

  ct_enrich <- bind_rows(
    combine_enrich(er_rvf$up,      'up',   'RVF_vs_NF'),
    combine_enrich(er_rvf$down,    'down', 'RVF_vs_NF'),
    combine_enrich(er_prv$up,      'up',   'pRV_vs_NF'),
    combine_enrich(er_prv$down,    'down', 'pRV_vs_NF'),
    combine_enrich(er_rvf_prv$up,  'up',   'RVF_vs_pRV'),
    combine_enrich(er_rvf_prv$down,'down', 'RVF_vs_pRV')
  )
  if (!is.null(ct_enrich) && nrow(ct_enrich) > 0) {
    ct_enrich$cell_type <- ct
    nkt_enrich_combined[[ct]] <- ct_enrich
  }

  pdf(paste0('./output/Xenium/nkt_', ct_slug, '_enrichr_plots.pdf'), width = 12, height = 7)
  for (comp_label in c('RVF_vs_NF', 'pRV_vs_NF', 'RVF_vs_pRV')) {
    er <- switch(comp_label, RVF_vs_NF = er_rvf, pRV_vs_NF = er_prv, RVF_vs_pRV = er_rvf_prv)
    for (dir_label in c('up', 'down')) {
      er_dir <- er[[dir_label]]
      if (!is.null(er_dir)) {
        p1 <- plot_enrich(er_dir$GO_Biological_Process_2023, paste(ct, comp_label, toupper(dir_label), 'GO BP'))
        p2 <- plot_enrich(er_dir$Reactome_2022,              paste(ct, comp_label, toupper(dir_label), 'Reactome'))
        if (!is.null(p1)) print(p1)
        if (!is.null(p2)) print(p2)
      }
    }
  }
  dev.off()
}

nkt_all_enrich <- bind_rows(nkt_enrich_combined)
write.csv(nkt_all_enrich, './output/Xenium/nkt_all_enrichr.csv', row.names = FALSE)


##############################################
#### EC subclustering
##############################################
ref <- readRDS('../snRV_ref.rds')                          
ec.ref <- subset(ref,subset = Names %in% c('EC','Endo','LEC'))
ec.ref <- subset(ref,subset = Subnames %in% c('EC_Arterial','EC_Capillary','EC_Endocardial','EC_Lymph','EC_Venous'))


ec.obj <- subset(xenium.obj, subset = cell_type_rctd_doublet %in% c('EC', 'Endo', 'LEC'))
ec.obj <- JoinLayers(ec.obj)

umap_cols <- grep("^umap_\\d+$|^UMAP_\\d+$|^umaporig_\\d+$",
                  colnames(ec.obj@meta.data), value = TRUE)
if (length(umap_cols)) ec.obj@meta.data[umap_cols] <- NULL

ec.obj[["pca"]]      <- NULL
ec.obj[["harmony"]]  <- NULL
ec.obj[["umap"]]     <- NULL
ec.obj[["umap_orig"]] <- NULL

ec.obj <- FindVariableFeatures(ec.obj)
ec.obj <- ScaleData(ec.obj)
ec.obj <- RunPCA(ec.obj, npcs = 30)
ec.obj <- RunHarmony(ec.obj, "patient")

ec.obj <- RunUMAP(ec.obj, reduction = "harmony", dims = 1:30, reduction.name = 'umap')
ec.obj <- RunUMAP(ec.obj, reduction = "pca",     dims = 1:30, reduction.name = 'umap_orig')

ec.obj <- FindNeighbors(ec.obj, reduction = "harmony", dims = 1:30)
ec.obj <- FindClusters(ec.obj, resolution = 0.5)

ec.marks <- FindAllMarkers(ec.obj)
write.csv(ec.marks, './output/Xenium/ec.marks.csv')

DimPlot(ec.obj, group.by = 'seurat_clusters', label = TRUE, repel = TRUE, reduction = 'umap_orig') +
  theme(legend.position = 'none')
DimPlot(ec.obj, group.by = 'patient', reduction = 'umap_orig')

# ── Annotate clusters after reviewing ec.marks.csv + canonical marker DotPlot ──
#
#   ┌─────────┬──────────────────┬────────────────────────────────────────┐
#   │ Cluster │    Cell Type     │             Key Markers                │
#   ├─────────┼──────────────────┼────────────────────────────────────────┤
#   │ 0       │ Capillary_EC     │ CD300LG, EMCN+                         │
#   │ 1       │ Cardiomyocyte    │ HSPB3, CAVIN4, KCNJ2                   │
#   │ 2       │ Arterial_EC      │ DLL1, SEMA4C                           │
#   │ 3       │ Fibroblast       │ PI16, CILP, FNDC1, SERPINA3            │
#   │ 4       │ Capillary_EC     │ TM4SF1, SEMA3G, EMCN+                  │
#   │ 5       │ Neutrophil       │ HVCN1, SLC22A4, SLC22A5               │
#   │ 6       │ Myeloid          │ C1QB, TYROBP, C5AR2                    │
#   │ 7       │ Venous_EC        │ ACKR1 (bright), EMCN+                  │
#   │ 8       │ Pericyte         │ KCNJ8, KCNE4, CCN4                     │
#   │ 9       │ Endocardial      │ NPR3+, CDH11 (highest), EMCN+          │
#   │ 10      │ Pericyte         │ RGS5, ENPEP, NOTCH3                    │
#   │ 11      │ Endocardial      │ CDH11+, EMCN+                          │
#   │ 12      │ Unknown          │ —                                      │
#   │ 13      │ LEC              │ PROX1, PDPN, LYVE1, CCL21, STAB2, FLT4│
#   │ 14      │ TNK_Cell         │ ZAP70, IL2RB                           │
#   │ 15      │ Pericyte         │ RGS5, NOTCH3, CARMN                    │
#   │ 16      │ Unknown          │ too few cells for markers              │
#   │ 17      │ Venous_EC        │ EMCN+, NR2F2+                          │
#   │ 18      │ Venous_EC        │ EMCN+ (highest), NR2F2+                │
#   └─────────┴──────────────────┴────────────────────────────────────────┘

ec_labels <- c(
  '0'  = 'Capillary_EC',
  '1'  = 'Cardiomyocyte',
  '2'  = 'Arterial_EC',
  '3'  = 'Fibroblast',
  '4'  = 'Capillary_EC',
  '5'  = 'Neutrophil',
  '6'  = 'Myeloid',
  '7'  = 'Venous_EC',
  '8'  = 'Pericyte',
  '9'  = 'Endocardial',
  '10' = 'Pericyte',
  '11' = 'Endocardial',
  '12' = 'Unknown',
  '13' = 'LEC',
  '14' = 'TNK_Cell',
  '15' = 'Pericyte',
  '16' = 'Unknown',
  '17' = 'Venous_EC',
  '18' = 'Venous_EC'
)

ec.obj$ec_subtype <- unname(ec_labels[as.character(ec.obj$seurat_clusters)])

DimPlot(ec.obj, group.by = 'ec_subtype', label = TRUE, repel = TRUE, reduction = 'umap_orig') +
  theme(legend.position = 'none')

ec.obj.clean.clean <- subset(ec.obj, subset = ec_subtype %in% c('Capillary_EC', 'Arterial_EC', 'Venous_EC', 'Endocardial', 'LEC'))

ec.obj.clean.clean[["pca"]]      <- NULL
ec.obj.clean.clean[["harmony"]]  <- NULL
ec.obj.clean.clean[["umap"]]     <- NULL
ec.obj.clean.clean[["umap_orig"]] <- NULL

umap_cols <- grep("^umap_\\d+$|^UMAP_\\d+$|^umaporig_\\d+$",
                  colnames(ec.obj.clean.clean@meta.data), value = TRUE)
if (length(umap_cols)) ec.obj.clean.clean@meta.data[umap_cols] <- NULL

ec.obj.clean.clean <- JoinLayers(ec.obj.clean.clean)
ec.obj.clean.clean <- FindVariableFeatures(ec.obj.clean.clean)
ec.obj.clean.clean <- ScaleData(ec.obj.clean.clean)
ec.obj.clean.clean <- RunPCA(ec.obj.clean.clean, npcs = 30)
ec.obj.clean.clean <- RunHarmony(ec.obj.clean.clean, "patient")
ec.obj.clean.clean <- RunUMAP(ec.obj.clean.clean, reduction = "harmony", dims = 1:30, reduction.name = 'umap')
ec.obj.clean.clean <- RunUMAP(ec.obj.clean.clean, reduction = "pca",     dims = 1:30, reduction.name = 'umap_orig')

ec.obj.clean.clean <- FindNeighbors(ec.obj.clean.clean, reduction = "harmony", dims = 1:30)
ec.obj.clean.clean <- FindClusters(ec.obj.clean.clean, resolution = 0.5)

ec.obj.clean.clean.mem <- ec.obj.clean.clean
ec.obj.clean.clean.mem[["Xenium"]] <- CreateAssayObject(                                                  
    counts = as(GetAssayData(ec.obj.clean.clean, layer = "counts"), "dgCMatrix")                            
)
ec.marks <- FindAllMarkers(ec.obj.clean.clean.mem)                                                        
rm(ec.obj.clean.clean.mem)  
write.csv(ec.marks, './output/Xenium/ec.clean.marks.csv')

# ── Annotate EC clean clusters ───────────────────────────────────────────────
#
#   ┌─────────┬────────────────┬──────────────────────────────────────────┐
#   │ Cluster │   Cell Type    │             Key Markers                  │
#   ├─────────┼────────────────┼──────────────────────────────────────────┤
#   │ 0       │ Capillary_EC   │ FLT1, CD300LG, PTGDS, TINAGL1           │
#   │ 1       │ Capillary_EC   │ PECAM1+, CDH5+, VWF+, FLT1+             │
#   │ 2       │ Arterial_EC    │ EFNB2, TMEM100, PECAM1, ACTA2           │
#   │ 3       │ Arterial_EC    │ DLL1, TNC                                │
#   │ 4       │ Capillary_EC   │ SMAGP, CCM2L, MICALL1, DUSP5            │
#   │ 5       │ Venous_EC      │ ACKR1, SELE, GATA2                       │
#   │ 6       │ Endocardial    │ NPR3+, CDH11+ (bright)                   │
#   │ 7       │ Endocardial    │ NPPB, MYL2 + PECAM1/CDH5/VWF+           │
#   │ 8       │ LEC            │ PROX1, LYVE1, PDPN, MMRN1, EDNRB        │
#   │ 9       │ Endocardial    │ CDH11 dominant                           │
#   │ 10      │ Cardiomyocyte  │ DES, TNNT2, MYH7                        │
#   │ 11      │ Schwann_Cell   │ GRIK3, ERBB3, ASPA                      │
#   │ 12      │ Proliferating_EC│ TTK, CHEK2, CENPU, PARPBP              │
#   │ 13      │ Activated_EC   │ IFNAR1, TRIM5, TNFSF10 + PECAM1/CDH5+   │
#   └─────────┴────────────────┴──────────────────────────────────────────┘

ec_clean_labels <- c(
  '0'  = 'Capillary_EC',
  '1'  = 'Capillary_EC',
  '2'  = 'Arterial_EC',
  '3'  = 'Arterial_EC',
  '4'  = 'Capillary_EC',
  '5'  = 'Venous_EC',
  '6'  = 'Endocardial',
  '7'  = 'Endocardial',
  '8'  = 'LEC',
  '9'  = 'Endocardial',
  '10' = 'Cardiomyocyte',
  '11' = 'Schwann_Cell',
  '12' = 'Proliferating_EC',
  '13' = 'Activated_EC'
)

ec.obj.clean.clean$ec_subtype <- unname(ec_clean_labels[as.character(ec.obj.clean.clean$seurat_clusters)])

DimPlot(ec.obj.clean.clean, group.by = 'ec_subtype', label = TRUE, repel = TRUE, reduction = 'umap_orig') +
  theme(legend.position = 'none')

ec.obj.clean.clean.clean <- subset(ec.obj.clean.clean,
  subset = ec_subtype %in% c('Capillary_EC', 'Arterial_EC', 'Venous_EC',
                              'Endocardial', 'LEC', 'Proliferating_EC', 'Activated_EC'))

DimPlot(ec.obj.clean.clean.clean, group.by = 'ec_subtype', label = TRUE, repel = TRUE, reduction = 'umap_orig') +
  theme(legend.position = 'none')
DimPlot(ec.obj.clean.clean.clean, group.by = 'patient', reduction = 'umap_orig')

saveRDS(ec.obj.clean.clean.clean, './output/Xenium/ec_clean_clean.rds')

# Label transfer from snRV reference
ec.ref <- JoinLayers(ec.ref)
ec.anchors <- FindTransferAnchors(
  reference   = ec.ref,
  query       = ec.obj.clean.clean.clean,
  dims        = 1:30,
  reference.reduction = 'pca'
)
ec.predictions <- TransferData(
  anchorset  = ec.anchors,
  refdata    = ec.ref$Subnames,
  dims       = 1:30
)
ec.obj.clean.clean.clean <- AddMetaData(ec.obj.clean.clean.clean, metadata = ec.predictions)

DimPlot(ec.obj.clean.clean.clean, group.by = 'predicted.id', label = TRUE, repel = TRUE, reduction = 'umap_orig') +
  theme(legend.position = 'none') +
  ggtitle('Reference predicted labels')

# Cross-tab: predicted label vs seurat cluster
table(ec.obj.clean.clean$predicted.id, ec.obj.clean.clean$seurat_clusters)

saveRDS(ec.obj.clean.clean.clean, './output/Xenium/ec_clean_clean.rds')

DimPlot(ec.obj.clean.clean.clean, group.by = 'seurat_clusters', label = TRUE, repel = TRUE, reduction = 'umap_orig') +
  theme(legend.position = 'none') +
  ggtitle('Reference predicted labels')

macro.prolif <- subset(ec.obj.clean.clean.clean,subset = seurat_clusters == 12)
saveRDS(macro.prolif, './output/Xenium/macro.prolif.rds')

ec.obj.clean.clean.clean <- subset(ec.obj.clean.clean.clean,subset = seurat_clusters != 12)

saveRDS(ec.obj.clean.clean.clean, './output/Xenium/ec_clean_clean.rds')


pdf('./output/Xenium/Xenium_ec_snUMAP_new.pdf', width = 5, height = 5)
PlotEmbedding(ec.obj.clean.clean.clean, group.by = 'predicted.id', reduction = 'umap_orig', point_size = .05,
              plot_under = TRUE, plot_theme = umap_theme() + NoLegend(),
              raster_dpi = 400, raster_scale = 0.5)
dev.off()

pdf('./output/Xenium/Xenium_ec_snUMAP_new_pt.pdf', width = 5, height = 5)
PlotEmbedding(ec.obj.clean.clean.clean, group.by = 'patient', reduction = 'umap_orig', point_size = .05,
              plot_under = TRUE, plot_theme = umap_theme() + NoLegend(),
              raster_dpi = 400, raster_scale = 0.5)
dev.off()

pdf('./output/Xenium/ec_dot.pdf', width = 8, height = 5)
# Canonical EC subtype marker DotPlot
DotPlot(ec.obj.clean.clean.clean,
        features = c(
          # Arterial
          "EFNB2", "TMEM100",
          # Capillary
          "CD300LG", "RGCC", "CA4",
          # Endocardial
          "NPR3", "CDH11",
                    # LEC
          "PROX1", "LYVE1", "PDPN", "CCL21",
                    # Venous
          "ACKR1", "NR2F2"
        ),
        group.by = "predicted.id",col.min=0,col.max=2) +
  RotatedAxis() +
  theme(axis.text.x = element_text(size = 9)) +
  ggtitle("Canonical EC subtype markers")
dev.off()


##############################################
#### EC proportion by disease state
##############################################

ec_meta <- ec.obj.clean.clean.clean@meta.data %>%
  select(patient, group, predicted.id)

ec_totals <- xenium.obj@meta.data %>%
  group_by(patient, group) %>%
  summarise(total = n(), .groups = 'drop')

ec_prop <- ec_meta %>%
  group_by(patient, group, predicted.id) %>%
  summarise(n = n(), .groups = 'drop') %>%
  left_join(ec_totals, by = c('patient', 'group')) %>%
  mutate(proportion = n / total) %>%
  complete(nesting(patient, group), predicted.id, fill = list(n = 0, proportion = 0)) %>%
  mutate(group = factor(group, levels = c('NF', 'pRV', 'RVF')))

pairwise_comparisons <- list(c('NF', 'pRV'), c('NF', 'RVF'), c('pRV', 'RVF'))

p_ec_prop <- ggplot(ec_prop, aes(x = group, y = proportion, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.5) +
  stat_compare_means(comparisons = pairwise_comparisons, method = 'wilcox.test', label = 'p.signif') +
  stat_compare_means(method = 'kruskal.test', label.y.npc = 0.95, size = 3) +
  facet_wrap(~ predicted.id, scales = 'free_y') +
  scale_x_discrete(limits = c('NF', 'pRV', 'RVF')) +
  theme_bw() +
  theme(legend.position = 'none',
        strip.text = element_text(face = 'bold'),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = 'Proportion of all cells', title = 'EC subtype distribution by group')

pdf('./output/Xenium/ec_proportions.pdf', width = 8, height = 6)
print(p_ec_prop)
dev.off()

##############################################
#### EC fractional area coverage by disease state
##############################################

cell_areas <- read.csv('./output/Xenium/cell_areas_2d.csv', row.names = 1)

# EC cell metadata with areas
ec_area_meta <- ec.obj.clean.clean.clean@meta.data %>%
  tibble::rownames_to_column('cell') %>%
  select(cell, patient, group, ec_subtype) %>%
  left_join(cell_areas, by = 'cell')

# Total tissue area per patient (denominator = all cells)
total_area <- xenium.obj@meta.data %>%
  tibble::rownames_to_column('cell') %>%
  select(cell, patient, group) %>%
  left_join(cell_areas, by = 'cell') %>%
  group_by(patient, group) %>%
  summarise(total_area = sum(cell_area_2d, na.rm = TRUE), .groups = 'drop')

# Fractional area per EC subtype per patient
ec_area_prop <- ec_area_meta %>%
  group_by(patient, group, ec_subtype) %>%
  summarise(ec_area = sum(cell_area_2d, na.rm = TRUE), .groups = 'drop') %>%
  left_join(total_area, by = c('patient', 'group')) %>%
  mutate(frac_area = ec_area / total_area) %>%
  complete(nesting(patient, group), ec_subtype, fill = list(ec_area = 0, frac_area = 0)) %>%
  mutate(group = factor(group, levels = c('NF', 'pRV', 'RVF')))

pairwise_comparisons <- list(c('NF', 'pRV'), c('NF', 'RVF'), c('pRV', 'RVF'))

p_ec_area <- ggplot(ec_area_prop, aes(x = group, y = frac_area, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.5) +
  stat_compare_means(comparisons = pairwise_comparisons, method = 'wilcox.test', label = 'p.signif') +
  stat_compare_means(method = 'kruskal.test', label.y.npc = 0.95, size = 3) +
  facet_wrap(~ ec_subtype, scales = 'free_y') +
  scale_x_discrete(limits = c('NF', 'pRV', 'RVF')) +
  theme_bw() +
  theme(legend.position = 'none',
        strip.text = element_text(face = 'bold'),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = 'Fractional area coverage', title = 'EC subtype fractional area by group')

pdf('./output/Xenium/ec_area_coverage.pdf', width = 10, height = 8)
print(p_ec_area)
dev.off()




##############################################
#### Pseudobulk DESeq2 by EC subtype
##############################################

ec_pseudo_meta <- unique(ec.obj.clean.clean@meta.data[, c('patient', 'group', 'ec_subtype')])
ec_pseudo_meta$patient    <- as.character(ec_pseudo_meta$patient)
ec_pseudo_meta$ec_subtype <- as.character(ec_pseudo_meta$ec_subtype)
ec_pseudo_meta$pseudobulk_id <- paste0(ec_pseudo_meta$patient, '..', ec_pseudo_meta$ec_subtype)
ec_pseudo_meta$group <- factor(ec_pseudo_meta$group, levels = c('NF', 'pRV', 'RVF'))
rownames(ec_pseudo_meta) <- ec_pseudo_meta$pseudobulk_id

ec.obj.clean.clean$pseudobulk_id <- paste0(
  as.character(ec.obj.clean.clean$patient), '..',
  as.character(ec.obj.clean.clean$ec_subtype)
)

ec_pseudo_counts <- AggregateExpression(
  ec.obj.clean.clean,
  assays = 'Xenium',
  group.by = 'pseudobulk_id',
  slot = 'counts',
  return.seurat = FALSE
)$Xenium

colnames(ec_pseudo_counts) <- gsub('-', '_', sub('^g', '', colnames(ec_pseudo_counts)))

ec_cell_types <- unique(ec_pseudo_meta$ec_subtype)
ec_deseq_results <- list()

for (ct in ec_cell_types) {
  ct_samples <- ec_pseudo_meta$pseudobulk_id[ec_pseudo_meta$ec_subtype == ct]
  ct_samples <- ct_samples[ct_samples %in% colnames(ec_pseudo_counts)]
  if (length(ct_samples) < 2) next

  ct_counts <- ec_pseudo_counts[, ct_samples, drop = FALSE]
  ct_meta   <- ec_pseudo_meta[ct_samples, , drop = FALSE]

  if (length(unique(ct_meta$group)) < 2) next

  keep <- rowSums(ct_counts >= 5) >= 2
  ct_counts <- ct_counts[keep, , drop = FALSE]
  if (nrow(ct_counts) < 10) next

  dds <- DESeqDataSetFromMatrix(countData = ct_counts, colData = ct_meta, design = ~ group)
  dds <- tryCatch(DESeq(dds), error = function(e) NULL)
  if (is.null(dds)) next

  for (comp in list(c('RVF', 'NF'), c('pRV', 'NF'), c('RVF', 'pRV'))) {
    comp_label <- paste0(comp[1], '_vs_', comp[2])
    res <- tryCatch(
      results(dds, contrast = c('group', comp[1], comp[2])) %>%
        as.data.frame() %>%
        tibble::rownames_to_column('gene') %>%
        mutate(cell_type = ct, comparison = comp_label) %>%
        arrange(padj),
      error = function(e) NULL
    )
    if (!is.null(res)) ec_deseq_results[[paste0(ct, '_', comp_label)]] <- res
  }
}

ec_all_deseq <- bind_rows(ec_deseq_results)
write.csv(ec_all_deseq, './output/Xenium/ec_pseudobulk_deseq2.csv', row.names = FALSE)

ec_sig_summary <- ec_all_deseq %>%
  filter(padj < 0.05, abs(log2FoldChange) > 0.5) %>%
  group_by(cell_type, comparison) %>%
  summarise(n_up = sum(log2FoldChange > 0), n_down = sum(log2FoldChange < 0), .groups = 'drop')
print(ec_sig_summary)


##############################################
#### EC subtypes: Volcano + GO/Reactome
##############################################

ec_enrich_combined <- list()

for (ct in unique(ec_all_deseq$cell_type)) {
  ct_slug <- gsub('_', '-', ct)
  cat('EC processing:', ct, '\n')

  dat_rvf     <- ec_all_deseq %>% filter(cell_type == ct, comparison == 'RVF_vs_NF')  %>% cap_pval()
  dat_prv     <- ec_all_deseq %>% filter(cell_type == ct, comparison == 'pRV_vs_NF')  %>% cap_pval()
  dat_rvf_prv <- ec_all_deseq %>% filter(cell_type == ct, comparison == 'RVF_vs_pRV') %>% cap_pval()

  make_volcano_pdf(dat_rvf,     paste0('./output/Xenium/ec_', ct_slug, '_volcano_RVF_vs_NF.pdf'),  paste(ct, ': RVF vs NF'))
  make_volcano_pdf(dat_prv,     paste0('./output/Xenium/ec_', ct_slug, '_volcano_pRV_vs_NF.pdf'),  paste(ct, ': pRV vs NF'))
  make_volcano_pdf(dat_rvf_prv, paste0('./output/Xenium/ec_', ct_slug, '_volcano_RVF_vs_pRV.pdf'), paste(ct, ': RVF vs pRV'))

  run_enrich <- function(dat, comp) {
    sig  <- dat %>% filter(padj < 0.05, abs(log2FoldChange) > 0.5)
    up   <- sig %>% filter(log2FoldChange > 0) %>% pull(gene)
    down <- sig %>% filter(log2FoldChange < 0) %>% pull(gene)
    list(
      up   = if (length(up)   > 0) enrichr(up,   enrichR_dbs) else NULL,
      down = if (length(down) > 0) enrichr(down, enrichR_dbs) else NULL
    )
  }

  er_rvf     <- run_enrich(dat_rvf,     'RVF_vs_NF')
  er_prv     <- run_enrich(dat_prv,     'pRV_vs_NF')
  er_rvf_prv <- run_enrich(dat_rvf_prv, 'RVF_vs_pRV')

  ct_enrich <- bind_rows(
    combine_enrich(er_rvf$up,       'up',   'RVF_vs_NF'),
    combine_enrich(er_rvf$down,     'down', 'RVF_vs_NF'),
    combine_enrich(er_prv$up,       'up',   'pRV_vs_NF'),
    combine_enrich(er_prv$down,     'down', 'pRV_vs_NF'),
    combine_enrich(er_rvf_prv$up,   'up',   'RVF_vs_pRV'),
    combine_enrich(er_rvf_prv$down, 'down', 'RVF_vs_pRV')
  )
  if (!is.null(ct_enrich) && nrow(ct_enrich) > 0) {
    ct_enrich$cell_type <- ct
    ec_enrich_combined[[ct]] <- ct_enrich
  }

  pdf(paste0('./output/Xenium/ec_', ct_slug, '_enrichr_plots.pdf'), width = 12, height = 7)
  for (comp_label in c('RVF_vs_NF', 'pRV_vs_NF', 'RVF_vs_pRV')) {
    er <- switch(comp_label, RVF_vs_NF = er_rvf, pRV_vs_NF = er_prv, RVF_vs_pRV = er_rvf_prv)
    for (dir_label in c('up', 'down')) {
      er_dir <- er[[dir_label]]
      if (!is.null(er_dir)) {
        p1 <- plot_enrich(er_dir$GO_Biological_Process_2023, paste(ct, comp_label, toupper(dir_label), 'GO BP'))
        p2 <- plot_enrich(er_dir$Reactome_2022,              paste(ct, comp_label, toupper(dir_label), 'Reactome'))
        if (!is.null(p1)) print(p1)
        if (!is.null(p2)) print(p2)
      }
    }
  }
  dev.off()
}

ec_all_enrich <- bind_rows(ec_enrich_combined)
write.csv(ec_all_enrich, './output/Xenium/ec_all_enrichr.csv', row.names = FALSE)


##############################################
#### FB subclustering
##############################################
ref <- readRDS('../snRV_ref.rds')                          

fb.ref <- subset(ref, subset = Names %in% c('FB'))

fb.obj <- subset(xenium.obj, subset = cell_type_rctd_doublet %in% c('FB'))
fb.obj <- JoinLayers(fb.obj)

# Strip stale UMAP coordinate columns inherited from xenium.obj@meta.data.
# DimPlot's FetchData checks meta.data before @reductions, so leftover
# umap_1/umap_2 shadow any freshly injected reduction and DimPlot silently
# plots the old global UMAP instead of the subcluster integration.
umap_cols <- grep("^umap_\\d+$|^UMAP_\\d+$|^umaporig_\\d+$",
                  colnames(fb.obj@meta.data), value = TRUE)
if (length(umap_cols)) fb.obj@meta.data[umap_cols] <- NULL

fb.obj[["pca"]]       <- NULL
fb.obj[["harmony"]]   <- NULL
fb.obj[["umap"]]      <- NULL
fb.obj[["umap_orig"]] <- NULL

fb.obj <- FindVariableFeatures(fb.obj)
fb.obj <- ScaleData(fb.obj)
fb.obj <- RunPCA(fb.obj, npcs = 30)
fb.obj <- RunHarmony(fb.obj, group.by.vars = "patient",
                     theta            = 8,
                     lambda           = 0.5,
                     sigma            = 0.1,
                     max_iter         = 30,
                     epsilon.cluster  = -Inf,
                     epsilon.harmony  = -Inf,
                     plot_convergence = TRUE)


# RunUMAP silently ignores the requested reduction slot in some Seurat v5 /
# v3-assay combinations and builds from default-assay variable features instead.
# Build UMAPs directly with uwot and inject as DimReducObjects to guarantee
# the harmony embedding is actually the input.
library(uwot)
.build_umap_reduc <- function(obj, reduction, name, key_prefix, dims = NULL) {
  emb <- Embeddings(obj, reduction)
  if (!is.null(dims)) emb <- emb[, dims, drop = FALSE]
  um <- uwot::umap(emb, n_neighbors = 30, min_dist = 0.3, metric = "cosine")
  rownames(um) <- rownames(emb)
  colnames(um) <- paste0(key_prefix, 1:2)
  obj[[name]] <- CreateDimReducObject(embeddings = um,
                                      key = key_prefix,
                                      assay = DefaultAssay(obj))
  obj
}
fb.obj <- .build_umap_reduc(fb.obj, "harmony", "umap",      "umap_",     dims = 1:30)
fb.obj <- .build_umap_reduc(fb.obj, "pca",     "umap_orig", "umaporig_", dims = 1:30)

fb.obj <- FindNeighbors(fb.obj, reduction = "harmony", dims = 1:30)
fb.obj <- FindClusters(fb.obj, resolution = 0.5)

fb.obj.mem <- fb.obj
fb.obj.mem[["Xenium"]] <- CreateAssayObject(
  counts = as(GetAssayData(fb.obj, layer = "counts"), "dgCMatrix")
)
fb.marks <- FindAllMarkers(fb.obj.mem)
rm(fb.obj.mem)
write.csv(fb.marks, './output/Xenium/fb.marks.csv')

DimPlot(fb.obj, group.by = 'seurat_clusters', label = TRUE, repel = TRUE, reduction = 'umap') +
  theme(legend.position = 'none')
DimPlot(fb.obj, group.by = 'patient', reduction = 'umap')

# ── Annotate clusters after reviewing fb.marks.csv ───────────────────────────
#
#   ┌─────────┬────────────────────┬────────────────────────────────────────┐
#   │ Cluster │     Cell Type      │          Key Markers (top 20)          │
#   ├─────────┼────────────────────┼────────────────────────────────────────┤
#   │ 0       │ FB_PCOLCE2         │ PCOLCE2, F13A1, PDPN, COCH, CXCL6      │
#   │ 1       │ FB_Epicardial      │ WT1, GPX3, LUM, STC1, KAZN             │
#   │ 2       │ FB_Homeostatic     │ DCN, FBLN1, ELN, ASPN, DPT, SFRP4, OGN │
#   │ 3       │ drop_CM            │ ACTC1, MYH7, NPPA, NPPB, TNNI3, MYBPC3 │
#   │ 4       │ FB_Activated       │ HAPLN1, FNDC1, TGFB2, F2R, AEBP1       │
#   │ 5       │ FB_Resting  │ CILP, SERPINA3, HGF, WNT5B, PTGIR      │
#   │ 6       │ drop_LowConf       │ junk markers — no clear identity       │
#   │ 7       │ drop_Ambig         │ IL7R, COL28A1, PTGFR — no FB signature │
#   │ 8       │ drop_Macrophage    │ LYVE1, OLR1, PLA2G2A, PLA2G7, RNASE1   │
#   │ 9       │ drop_Mast          │ CPA3, KIT, CTSG, CXCL8, CCL2           │
#   │ 10      │ drop_Pericyte      │ KCNJ8, CCN4, EPOR, KCNE4               │
#   │ 11      │ drop_Mural         │ RGS5, NOTCH3, ANGPT2, ENPEP, BMP5      │
#   │ 12      │ drop_Ambig         │ FHL2, ITIH4, NTF3 — heterogeneous      │
#   │ 13      │ drop_EC_Art        │ GJA5, DKK2, JAG2, CLDN11               │
#   │ 14      │ drop_cDC_Mono      │ CD1C, FCER1A, HLA-DMA, HLA-DRB5, FCN1  │
#   │ 15      │ drop_Adipocyte     │ PLIN4, ADIPOQ, LPL, PPARG, CIDEC       │
#   │ 16      │ drop_Tcell         │ LCK, ZAP70, IL2RB, ITK, CD6, CXCR4     │
#   │ 17      │ drop_CM2           │ DES, TNNT1, TNNT2                      │
#   │ 18      │ drop_CM_stressed   │ ANKRD1, TMEM65, HMGCS2                 │
#   │ 19      │ drop_Schwann       │ SOX10, ERBB3, ALK, PTPRZ1              │
#   │ 20      │ drop_LEC           │ MMRN1, CCL21, STAB2, PRG4              │
#   │ 21      │ drop_Mural2        │ RGS5, NOTCH3, ANGPT2 — dup of c11      │
#   │ 22      │ drop_Myeloid_Prolif│ MKI67+/TOP2A+ but CSF1R+, DCN/PDGFRA−  │
#   └─────────┴────────────────────┴────────────────────────────────────────┘

fb_labels <- c(
  '0'  = 'FB_PCOLCE2',
  '1'  = 'FB_Epicardial',
  '2'  = 'FB_Homeostatic',
  '3'  = 'drop_CM',
  '4'  = 'FB_Activated',
  '5'  = 'FB_Resting',
  '6'  = 'drop_LowConf',
  '7'  = 'drop_Ambig',
  '8'  = 'drop_Macrophage',
  '9'  = 'drop_Mast',
  '10' = 'drop_Pericyte',
  '11' = 'drop_Mural',
  '12' = 'drop_Ambig',
  '13' = 'drop_EC_Art',
  '14' = 'drop_cDC_Mono',
  '15' = 'drop_Adipocyte',
  '16' = 'drop_Tcell',
  '17' = 'drop_CM2',
  '18' = 'drop_CM_stressed',
  '19' = 'drop_Schwann',
  '20' = 'drop_LEC',
  '21' = 'drop_Mural2',
  '22' = 'drop_Myeloid_Prolif'
)

fb.obj$fb_subtype <- unname(fb_labels[as.character(fb.obj$seurat_clusters)])

DimPlot(fb.obj, group.by = 'fb_subtype', label = TRUE, repel = TRUE, reduction = 'umap') +
  theme(legend.position = 'none')

fb.obj.clean <- subset(fb.obj, subset = fb_subtype %in% c(
  'FB_PCOLCE2', 'FB_Epicardial', 'FB_Homeostatic', 'FB_Activated',
  'FB_Resting'
))

fb.obj.clean[["pca"]]       <- NULL
fb.obj.clean[["harmony"]]   <- NULL
fb.obj.clean[["umap"]]      <- NULL
fb.obj.clean[["umap_orig"]] <- NULL

umap_cols <- grep("^umap_\\d+$|^UMAP_\\d+$|^umaporig_\\d+$",
                  colnames(fb.obj.clean@meta.data), value = TRUE)
if (length(umap_cols)) fb.obj.clean@meta.data[umap_cols] <- NULL

fb.obj.clean <- JoinLayers(fb.obj.clean)
fb.obj.clean <- FindVariableFeatures(fb.obj.clean)
fb.obj.clean <- ScaleData(fb.obj.clean)
fb.obj.clean <- RunPCA(fb.obj.clean, npcs = 30)
fb.obj.clean <- RunHarmony(fb.obj.clean, "patient")
fb.obj.clean <- RunUMAP(fb.obj.clean, reduction = "harmony", dims = 1:30, reduction.name = 'umap')
fb.obj.clean <- RunUMAP(fb.obj.clean, reduction = "pca",     dims = 1:30, reduction.name = 'umap_orig')

fb.obj.clean <- FindNeighbors(fb.obj.clean, reduction = "harmony", dims = 1:30)
fb.obj.clean <- FindClusters(fb.obj.clean, resolution = 0.5)

fb.obj.clean.mem <- fb.obj.clean
fb.obj.clean.mem[["Xenium"]] <- CreateAssayObject(
  counts = as(GetAssayData(fb.obj.clean, layer = "counts"), "dgCMatrix")
)
fb.clean.marks <- FindAllMarkers(fb.obj.clean.mem)
rm(fb.obj.clean.mem)
write.csv(fb.clean.marks, './output/Xenium/fb.clean.marks.csv')

# ── Annotate FB clean clusters ───────────────────────────────────────────────
#
#   ┌─────────┬────────────────────┬────────────────────────────────────────┐
#   │ Cluster │     Cell Type      │           Key Markers                  │
#   ├─────────┼────────────────────┼────────────────────────────────────────┤
#   │ 0       │ FB_PCOLCE2         │ PCOLCE2, F13A1, COCH, PDPN, CXCL6      │
#   │ 1       │ FB_Resting  │ CILP, SERPINA3, HGF, WNT5B, NFATC4     │
#   │ 2       │ FB_NOX4            │ NOX4, THBS4, PDGFRA, POSTN, FAP, F2R   │
#   │ 3       │ FB_Activated       │ HAPLN1, FNDC1, TGFB2, COL16A1, AEBP1   │
#   │ 4       │ FB_Homeostatic     │ DCN, FBLN1, ASPN, SFRP4, ELN, DPT, OGN │
#   │ 5       │ drop_EC            │ CDH5, PECAM1, VWF, ESAM, ARHGEF15      │
#   │ 6       │ FB_Resting  │ CILP, SERPINA3, HGF (merged with c1)   │
#   │ 7       │ drop_Ambig         │ KIT/CPA3 p1<0.03 — not mast, no FB id  │
#   │ 8       │ drop_Mural         │ RGS5, NOTCH3, MYL9, ITGA7, CARMN       │
#   │ 9       │ FB_Stress          │ CXCL8, CXCL6, FOSL1, MT1E/X/M, HK2, PKM│
#   │ 10      │ FB_NTN4            │ NTN4, HAS1, STK32B, SLC4A4, CREB5      │
#   └─────────┴────────────────────┴────────────────────────────────────────┘
#
# ── Literature support for these cardiac FB subtypes ────────────────────────
#
# General human cardiac FB taxonomy (healthy + heart failure atlases):
#   • Litvinukova et al. Nature 2020 — "Cells of the adult human heart"
#     (Heart Cell Atlas; FB1–FB7 populations across LV/RV/atria)
#   • Tucker et al. Circulation 2020 — "Transcriptional and Cellular Diversity
#     of the Human Heart" (snRNA-seq atlas; DCN+/PDGFRA+ homeostatic FBs)
#   • Koenig et al. Nat Cardiovasc Res 2022 — "Single-cell transcriptomics
#     reveals cell-type-specific diversification in human heart failure"
#     (DCM; defines FB-Homeostatic, FB-Activated, FB-POSTN, PDPN+/F13A1+ FB)
#   • Chaffin et al. Nature 2022 — "Single-nucleus profiling of human dilated
#     and hypertrophic cardiomyopathy" (confirms activated FB expansion in HF)
#   • Reichart et al. Science 2022 — DCM/ARVC snRNA-seq; FB state transitions
#   • Amrute et al. Nature 2024 — multi-disease human HF atlas
#
# FB_Homeostatic (DCN+/FBLN1+/DPT+/SFRP4+/ELN+/OGN+):
#   • "FB1" in Litvinukova 2020; "FB-Homeostatic" in Koenig 2022 and
#     Chaffin 2022 — the dominant resting/quiescent cardiac FB population
#
# FB_Activated (POSTN+/FAP+/FNDC1+/HAPLN1+/TGFB2+/COL16A1+):
#   • Kanisicak et al. Nat Commun 2016 — POSTN+ myofibroblast lineage tracing
#     in injured mouse heart
#   • Ruiz-Villalba et al. Circulation 2020 — CTHRC1+ activated FBs post-MI
#   • "FB-Activated" / "FB-POSTN" state in Koenig 2022 and Chaffin 2022
#
# FB_NOX4 (NOX4+/THBS4+/POSTN+ profibrotic activated FB):
#   • Cucoranu et al. Circ Res 2005 — NOX4 mediates TGFβ1-induced cardiac
#     FB → myofibroblast differentiation
#   • Kuroda et al. PNAS 2010 — NOX4 as a major oxidative-stress source in
#     the failing heart
#   • Zhao et al. Redox Biol 2015 — NOX4 in cardiac fibrosis pathogenesis
#
# FB_Resting (CILP+/SERPINA3+/HGF+/LARP6+):
#   • Homeostatic CILP+/HGF+ fibroblast co-expressing collagen biosynthetic
#     machinery (LARP6) and anti-fibrotic mediators. Corresponds to the CILP+
#     resting fibroblast state enriched in non-failing human hearts (Koenig
#     et al. 2022; Amrute et al. 2024). Previously labeled "matrifibrocyte"
#     but lacks canonical mouse matrifibrocyte markers (COMP, CHAD, THBS4;
#     Fu et al. 2018) and is enriched in NF rather than scar tissue.
#
# FB_PCOLCE2 (PCOLCE2+/F13A1+/PDPN+/COCH+ interstitial FB):
#   • PDPN+/F13A1+ interstitial/inflammatory FB state in Koenig 2022
#   • Comparable F13A1+ FB population in Litvinukova 2020 (FB4)
#
# FB_Stress (CXCL8+/CXCL6+/FOSL1+/MT1E/X/M+/HK2+/PKM+ IEG + hypoxia response):
#   • Farbehi et al. eLife 2019 — injury-response cardiac FBs with IEG and
#     metallothionein/glycolysis stress signatures
#   • "Stressed" / "FB-activated-stress" state reported in Koenig 2022 and
#     Chaffin 2022 DCM snRNA-seq
#
# FB_NTN4 (NTN4+/HAS1+ hyaluronan-producing FB):
#   • Less canonical in cardiac atlases — HAS1+ matrix-remodeling FBs have
#     been described in pulmonary fibrosis (Tsukui et al. Nat Commun 2020)
#     and skin fibrosis (Tabib et al. Nat Commun 2021); NTN4 marks a
#     basement-membrane-producing FB subset in multiple tissues. Flagged
#     here as a putative novel cardiac state for downstream validation.

fb_clean_labels <- c(
  '0'  = 'FB_PCOLCE2',
  '1'  = 'FB_Resting',
  '2'  = 'FB_NOX4',
  '3'  = 'FB_Activated',
  '4'  = 'FB_Homeostatic',
  '5'  = 'drop_EC',
  '6'  = 'FB_Resting',
  '7'  = 'drop_Ambig',
  '8'  = 'drop_Mural',
  '9'  = 'FB_Stress',
  '10' = 'FB_NTN4'
)

fb.obj.clean$fb_subtype <- unname(fb_clean_labels[as.character(fb.obj.clean$seurat_clusters)])

DimPlot(fb.obj.clean, group.by = 'fb_subtype', label = TRUE, repel = TRUE, reduction = 'umap') +
  theme(legend.position = 'none')

fb.obj.clean.clean <- subset(fb.obj.clean, subset = fb_subtype %in% c(
  'FB_PCOLCE2', 'FB_Resting', 'FB_NOX4', 'FB_Activated',
  'FB_Homeostatic', 'FB_Stress', 'FB_NTN4'
))


fb.obj.clean.clean[["pca"]]       <- NULL
fb.obj.clean.clean[["harmony"]]   <- NULL
fb.obj.clean.clean[["umap"]]      <- NULL
fb.obj.clean.clean[["umap_orig"]] <- NULL

umap_cols <- grep("^umap_\\d+$|^UMAP_\\d+$|^umaporig_\\d+$",
                  colnames(fb.obj.clean.clean@meta.data), value = TRUE)
if (length(umap_cols)) fb.obj.clean.clean@meta.data[umap_cols] <- NULL

fb.obj.clean.clean <- FindVariableFeatures(fb.obj.clean.clean)
fb.obj.clean.clean <- ScaleData(fb.obj.clean.clean)
fb.obj.clean.clean <- RunPCA(fb.obj.clean.clean, npcs = 30)
fb.obj.clean.clean <- RunHarmony(fb.obj.clean.clean, "patient")
fb.obj.clean.clean <- RunUMAP(fb.obj.clean.clean, reduction = "harmony", dims = 1:30, reduction.name = 'umap')
fb.obj.clean.clean <- RunUMAP(fb.obj.clean.clean, reduction = "pca",     dims = 1:30, reduction.name = 'umap_orig')



DimPlot(fb.obj.clean.clean, group.by = 'fb_subtype', label = TRUE, repel = TRUE, reduction = 'umap') +
  theme(legend.position = 'none')
DimPlot(fb.obj.clean.clean, group.by = 'patient', reduction = 'umap')

saveRDS(fb.obj.clean.clean, './output/Xenium/fb_clean_clean.rds')

# Label transfer from snRV reference
fb.ref <- JoinLayers(fb.ref)
fb.anchors <- FindTransferAnchors(
  reference           = fb.ref,
  query               = fb.obj.clean.clean,
  dims                = 1:30,
  reference.reduction = 'pca'
)
fb.predictions <- TransferData(
  anchorset = fb.anchors,
  refdata   = fb.ref$Subnames,
  dims      = 1:30
)
fb.obj.clean.clean <- AddMetaData(fb.obj.clean.clean, metadata = fb.predictions)

DimPlot(fb.obj.clean.clean, group.by = 'predicted.id', label = TRUE, repel = TRUE, reduction = 'umap_orig') +
  theme(legend.position = 'none') +
  ggtitle('Reference predicted labels')

table(fb.obj.clean.clean$predicted.id, fb.obj.clean.clean$seurat_clusters)

# ── Cross-reference de novo fb_subtype vs snRNA reference predicted.id ───────
fb_cross <- as.data.frame(table(
  fb_subtype   = fb.obj.clean.clean$fb_subtype,
  predicted.id = fb.obj.clean.clean$predicted.id
))
# Row-normalize within fb_subtype (fraction of each de novo cluster that
# maps to each reference label)
fb_cross <- fb_cross %>%
  group_by(fb_subtype) %>%
  mutate(frac = Freq / sum(Freq)) %>%
  ungroup()

pdf('./output/Xenium/fb_subtype_vs_predictedid.pdf', width = 7, height = 5)
ggplot(fb_cross, aes(x = predicted.id, y = fb_subtype, fill = frac)) +
  geom_tile(color = 'grey90') +
  geom_text(aes(label = ifelse(frac > 0.02, sprintf('%.0f%%', 100 * frac), '')),
            size = 3) +
  scale_fill_gradient(low = 'white', high = '#b2182b',
                      limits = c(0, 1), name = 'frac of\nfb_subtype') +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = 'snRNA reference (predicted.id)',
       y = 'De novo Xenium fb_subtype',
       title = 'De novo FB clusters vs snRNA reference label transfer')
dev.off()

saveRDS(fb.obj.clean.clean, './output/Xenium/fb_clean_clean.rds')

pdf('./output/Xenium/Xenium_fb_snUMAP_new.pdf', width = 5, height = 5)
PlotEmbedding(fb.obj.clean.clean, group.by = 'fb_subtype', reduction = 'umap', point_size = .05,
              plot_under = TRUE, plot_theme = umap_theme() + NoLegend(),
              raster_dpi = 400, raster_scale = 0.5)
dev.off()

pdf('./output/Xenium/Xenium_fb_snUMAP_new_pt.pdf', width = 5, height = 5)
PlotEmbedding(fb.obj.clean.clean, group.by = 'patient', reduction = 'umap', point_size = .05,
              plot_under = TRUE, plot_theme = umap_theme() + NoLegend(),
              raster_dpi = 400, raster_scale = 0.5)
dev.off()

fb.obj.clean.clean$fb_subtype <- factor(
  fb.obj.clean.clean$fb_subtype,
  levels = c('FB_Homeostatic', 'FB_PCOLCE2', 'FB_Activated', 'FB_NOX4',
             'FB_Resting', 'FB_Stress', 'FB_NTN4')
)

pdf('./output/Xenium/fb_dot.pdf', width = 13, height = 5)
DotPlot(fb.obj.clean.clean,
        features = c(
          # FB_Activated (myofibroblast / POSTN+ lineage)
          "HAPLN1", "FNDC1", "TGFB2", "COL16A1", "AEBP1", "EDNRA",
          # FB_Homeostatic (resting, FB1-like)
          "ASPN", "SFRP4", "ELN", "DPT", "OGN", "LTBP2", "FBN1",
          # FB_Resting (Fu 2018; CILP+)
          "CILP", "SERPINA3", "HGF", "WNT5B", "NFATC4", "PTGIR",
          # FB_NOX4 (profibrotic NOX4+, peri-myocyte)
          "NOX4", "THBS4", "POSTN", "FAP", "F2R",
          # FB_NTN4 (HAS1+ matrix-producing)
          "NTN4", "HAS1", "SLC4A4", "CREB5",
          # FB_PCOLCE2 (PDPN+/F13A1+ interstitial)
          "PCOLCE2", "F13A1", "COCH", "PDPN",
          # FB_Stress (IEG + hypoxia / glycolysis)
          "CXCL8", "CXCL6", "FOSL1", "MT1E", "MT1X", "HK2", "PKM", "NAMPT",
          # Pan-fibroblast
          "DCN", "COL1A1", "FBLN1", "PDGFRA"
        ),
        group.by = "fb_subtype", col.min = 0, col.max = 2) +
  RotatedAxis() +
  theme(axis.text.x = element_text(size = 8)) +
  ggtitle("FB subtype markers (de novo fb_subtype)")
dev.off()


##############################################
#### FB proportion by disease state
##############################################

pairwise_comparisons <- list(c('NF', 'pRV'), c('NF', 'RVF'), c('pRV', 'RVF'))

fb_meta <- fb.obj.clean.clean@meta.data %>%
  select(patient, group, fb_subtype)

fb_totals <- xenium.obj@meta.data %>%
  group_by(patient, group) %>%
  summarise(total = n(), .groups = 'drop')

fb_prop <- fb_meta %>%
  group_by(patient, group, fb_subtype) %>%
  summarise(n = n(), .groups = 'drop') %>%
  left_join(fb_totals, by = c('patient', 'group')) %>%
  mutate(proportion = n / total) %>%
  complete(nesting(patient, group), fb_subtype, fill = list(n = 0, proportion = 0)) %>%
  mutate(group = factor(group, levels = c('NF', 'pRV', 'RVF')))

p_fb_prop <- ggplot(fb_prop, aes(x = group, y = proportion, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.5) +
  stat_compare_means(comparisons = pairwise_comparisons, method = 'wilcox.test', label = 'p.signif') +
  stat_compare_means(method = 'kruskal.test', label.y.npc = 0.95, size = 3) +
  facet_wrap(~ fb_subtype, scales = 'free_y') +
  scale_x_discrete(limits = c('NF', 'pRV', 'RVF')) +
  theme_bw() +
  theme(legend.position = 'none',
        strip.text = element_text(face = 'bold'),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = 'Proportion of all cells', title = 'FB subtype distribution by group')

pdf('./output/Xenium/fb_proportions.pdf', width = 8, height = 6)
print(p_fb_prop)
dev.off()

##############################################
#### FB fractional area coverage by disease state
##############################################

cell_areas <- read.csv('./output/Xenium/cell_areas_2d.csv', row.names = 1)

total_area <- xenium.obj@meta.data %>%
  tibble::rownames_to_column('cell') %>%
  left_join(cell_areas, by = 'cell') %>%
  group_by(patient, group) %>%
  summarise(total_area = sum(cell_area_2d, na.rm = TRUE), .groups = 'drop')

fb_area_meta <- fb.obj.clean.clean@meta.data %>%
  tibble::rownames_to_column('cell') %>%
  select(cell, patient, group, fb_subtype) %>%
  left_join(cell_areas, by = 'cell')

fb_area_prop <- fb_area_meta %>%
  group_by(patient, group, fb_subtype) %>%
  summarise(fb_area = sum(cell_area_2d, na.rm = TRUE), .groups = 'drop') %>%
  left_join(total_area, by = c('patient', 'group')) %>%
  mutate(frac_area = fb_area / total_area) %>%
  complete(nesting(patient, group), fb_subtype, fill = list(fb_area = 0, frac_area = 0)) %>%
  mutate(group = factor(group, levels = c('NF', 'pRV', 'RVF')))

p_fb_area <- ggplot(fb_area_prop, aes(x = group, y = frac_area, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.5) +
  stat_compare_means(comparisons = pairwise_comparisons, method = 'wilcox.test', label = 'p.signif') +
  stat_compare_means(method = 'kruskal.test', label.y.npc = 0.95, size = 3) +
  facet_wrap(~ fb_subtype, scales = 'free_y') +
  scale_x_discrete(limits = c('NF', 'pRV', 'RVF')) +
  theme_bw() +
  theme(legend.position = 'none',
        strip.text = element_text(face = 'bold'),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = 'Fractional area coverage', title = 'FB subtype fractional area by group')

pdf('./output/Xenium/fb_area_coverage.pdf', width = 10, height = 8)
print(p_fb_area)
dev.off()


##############################################
#### Pseudobulk DESeq2 by FB subtype
##############################################

fb_pseudo_meta <- unique(fb.obj.clean.clean@meta.data[, c('patient', 'group', 'fb_subtype')])
fb_pseudo_meta$patient    <- as.character(fb_pseudo_meta$patient)
fb_pseudo_meta$fb_subtype <- as.character(fb_pseudo_meta$fb_subtype)
fb_pseudo_meta$pseudobulk_id <- paste0(fb_pseudo_meta$patient, '..', fb_pseudo_meta$fb_subtype)
fb_pseudo_meta$group <- factor(fb_pseudo_meta$group, levels = c('NF', 'pRV', 'RVF'))
rownames(fb_pseudo_meta) <- fb_pseudo_meta$pseudobulk_id

fb.obj.clean.clean$pseudobulk_id <- paste0(
  as.character(fb.obj.clean.clean$patient), '..',
  as.character(fb.obj.clean.clean$fb_subtype)
)

fb_pseudo_counts <- AggregateExpression(
  fb.obj.clean.clean,
  assays = 'Xenium',
  group.by = 'pseudobulk_id',
  slot = 'counts',
  return.seurat = FALSE
)$Xenium

colnames(fb_pseudo_counts) <- gsub('-', '_', sub('^g', '', colnames(fb_pseudo_counts)))

fb_cell_types <- unique(fb_pseudo_meta$fb_subtype)
fb_deseq_results <- list()

for (ct in fb_cell_types) {
  ct_samples <- fb_pseudo_meta$pseudobulk_id[fb_pseudo_meta$fb_subtype == ct]
  ct_samples <- ct_samples[ct_samples %in% colnames(fb_pseudo_counts)]
  if (length(ct_samples) < 2) next

  ct_counts <- fb_pseudo_counts[, ct_samples, drop = FALSE]
  ct_meta   <- fb_pseudo_meta[ct_samples, , drop = FALSE]

  if (length(unique(ct_meta$group)) < 2) next

  keep <- rowSums(ct_counts >= 5) >= 2
  ct_counts <- ct_counts[keep, , drop = FALSE]
  if (nrow(ct_counts) < 10) next

  dds <- DESeqDataSetFromMatrix(countData = ct_counts, colData = ct_meta, design = ~ group)
  dds <- tryCatch(DESeq(dds), error = function(e) NULL)
  if (is.null(dds)) next

  for (comp in list(c('RVF', 'NF'), c('pRV', 'NF'), c('RVF', 'pRV'))) {
    comp_label <- paste0(comp[1], '_vs_', comp[2])
    res <- tryCatch(
      results(dds, contrast = c('group', comp[1], comp[2])) %>%
        as.data.frame() %>%
        tibble::rownames_to_column('gene') %>%
        mutate(cell_type = ct, comparison = comp_label) %>%
        arrange(padj),
      error = function(e) NULL
    )
    if (!is.null(res)) fb_deseq_results[[paste0(ct, '_', comp_label)]] <- res
  }
}

fb_all_deseq <- bind_rows(fb_deseq_results)
write.csv(fb_all_deseq, './output/Xenium/fb_pseudobulk_deseq2.csv', row.names = FALSE)

fb_sig_summary <- fb_all_deseq %>%
  filter(padj < 0.05, abs(log2FoldChange) > 0.5) %>%
  group_by(cell_type, comparison) %>%
  summarise(n_up = sum(log2FoldChange > 0), n_down = sum(log2FoldChange < 0), .groups = 'drop')
print(fb_sig_summary)


##############################################
#### FB subtypes: Volcano + GO/Reactome
##############################################

# Load original panel genes from corrected (non-imputed) object
xenium.base <- readRDS('./Xenium_resegmented_corrected.rds')
panel_genes <- rownames(xenium.base)
rm(xenium.base); gc()

# Helper: cap -log10(padj) at 5 by flooring padj at 1e-5
cap_pval <- function(df) mutate(df, padj = ifelse(is.na(padj), NA_real_, pmax(padj, 1e-5)))

# Helper: build custom color vector — red for panel genes, blue for imputed, grey for NS
make_volcano_colors <- function(df, panel_genes) {
  is_panel <- df$gene %in% panel_genes
  is_sig   <- !is.na(df$padj) & df$padj < 0.05 & abs(df$log2FoldChange) > 0.5
  cols <- case_when(
    is_sig & is_panel  ~ 'red3',
    is_sig & !is_panel ~ 'royalblue3',
    !is_sig & is_panel ~ 'salmon',
    TRUE               ~ 'grey80'
  )
  names(cols) <- df$gene
  cols
}

# Helper: top N significant gene labels
sig_labs <- function(df, n = 30) {
  df %>%
    filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 0.5) %>%
    arrange(padj) %>%
    head(n) %>%
    pull(gene)
}

volcano_args <- list(
  x = 'log2FoldChange', y = 'padj',
  pCutoff = 0.05, FCcutoff = 0.5,
  pointSize = 2, labSize = 3.5,
  drawConnectors = TRUE, widthConnectors = 0.5,
  xlim = c(-5, 5), ylim = c(0, 5),
  subtitle = 'Red = original panel gene | Blue = imputed gene',
  legendPosition = 'none'
)

enrichR_dbs <- c('GO_Biological_Process_2023', 'Reactome_2022')
setEnrichrSite('Enrichr')

combine_enrich <- function(enrich_res, direction, comparison) {
  if (is.null(enrich_res)) return(NULL)
  lapply(names(enrich_res), function(db) {
    df <- enrich_res[[db]]
    if (nrow(df) == 0) return(NULL)
    df %>%
      mutate(Term = as.character(Term),
             database = db, direction = direction, comparison = comparison)
  }) %>% bind_rows()
}

plot_enrich <- function(enrich_df, title_str) {
  if (is.null(enrich_df) || nrow(enrich_df) == 0) return(NULL)
  enrich_df %>%
    arrange(Adjusted.P.value) %>%
    head(10) %>%
    mutate(Term = factor(Term, levels = rev(Term))) %>%
    ggplot(aes(x = -log10(Adjusted.P.value), y = Term)) +
    geom_col(fill = 'steelblue') +
    theme_bw() +
    labs(x = '-log10(adj. P)', y = NULL, title = title_str)
}

  # --- Volcano plots ---
  make_volcano_pdf <- function(dat, filepath, title_str) {
    if (nrow(dat) == 0) { cat('  Skipping', title_str, '— no data\n'); return(invisible(NULL)) }
    pdf(filepath, width = 10, height = 8)
    tryCatch(
      print(do.call(EnhancedVolcano, c(list(toptable = dat, lab = dat$gene,
        title = title_str,
        selectLab = sig_labs(dat),
        colCustom = make_volcano_colors(dat, panel_genes)), volcano_args))),
      error = function(e) { cat('  Volcano error for', title_str, ':', conditionMessage(e), '\n') }
    )
    dev.off()
  }


fb_enrich_combined <- list()

for (ct in unique(fb_all_deseq$cell_type)) {
  ct_slug <- gsub('_', '-', ct)
  cat('FB processing:', ct, '\n')

  dat_rvf     <- fb_all_deseq %>% filter(cell_type == ct, comparison == 'RVF_vs_NF')  %>% cap_pval()
  dat_prv     <- fb_all_deseq %>% filter(cell_type == ct, comparison == 'pRV_vs_NF')  %>% cap_pval()
  dat_rvf_prv <- fb_all_deseq %>% filter(cell_type == ct, comparison == 'RVF_vs_pRV') %>% cap_pval()

  make_volcano_pdf(dat_rvf,     paste0('./output/Xenium/fb_', ct_slug, '_volcano_RVF_vs_NF.pdf'),  paste(ct, ': RVF vs NF'))
  make_volcano_pdf(dat_prv,     paste0('./output/Xenium/fb_', ct_slug, '_volcano_pRV_vs_NF.pdf'),  paste(ct, ': pRV vs NF'))
  make_volcano_pdf(dat_rvf_prv, paste0('./output/Xenium/fb_', ct_slug, '_volcano_RVF_vs_pRV.pdf'), paste(ct, ': RVF vs pRV'))

  run_enrich <- function(dat, comp) {
    sig  <- dat %>% filter(padj < 0.05, abs(log2FoldChange) > 0.5)
    up   <- sig %>% filter(log2FoldChange > 0) %>% pull(gene)
    down <- sig %>% filter(log2FoldChange < 0) %>% pull(gene)
    list(
      up   = if (length(up)   > 0) enrichr(up,   enrichR_dbs) else NULL,
      down = if (length(down) > 0) enrichr(down, enrichR_dbs) else NULL
    )
  }

  er_rvf     <- run_enrich(dat_rvf,     'RVF_vs_NF')
  er_prv     <- run_enrich(dat_prv,     'pRV_vs_NF')
  er_rvf_prv <- run_enrich(dat_rvf_prv, 'RVF_vs_pRV')

  ct_enrich <- bind_rows(
    combine_enrich(er_rvf$up,       'up',   'RVF_vs_NF'),
    combine_enrich(er_rvf$down,     'down', 'RVF_vs_NF'),
    combine_enrich(er_prv$up,       'up',   'pRV_vs_NF'),
    combine_enrich(er_prv$down,     'down', 'pRV_vs_NF'),
    combine_enrich(er_rvf_prv$up,   'up',   'RVF_vs_pRV'),
    combine_enrich(er_rvf_prv$down, 'down', 'RVF_vs_pRV')
  )
  if (!is.null(ct_enrich) && nrow(ct_enrich) > 0) {
    ct_enrich$cell_type <- ct
    fb_enrich_combined[[ct]] <- ct_enrich
  }

  pdf(paste0('./output/Xenium/fb_', ct_slug, '_enrichr_plots.pdf'), width = 12, height = 7)
  for (comp_label in c('RVF_vs_NF', 'pRV_vs_NF', 'RVF_vs_pRV')) {
    er <- switch(comp_label, RVF_vs_NF = er_rvf, pRV_vs_NF = er_prv, RVF_vs_pRV = er_rvf_prv)
    for (dir_label in c('up', 'down')) {
      er_dir <- er[[dir_label]]
      if (!is.null(er_dir)) {
        p1 <- plot_enrich(er_dir$GO_Biological_Process_2023, paste(ct, comp_label, toupper(dir_label), 'GO BP'))
        p2 <- plot_enrich(er_dir$Reactome_2022,              paste(ct, comp_label, toupper(dir_label), 'Reactome'))
        if (!is.null(p1)) print(p1)
        if (!is.null(p2)) print(p2)
      }
    }
  }
  dev.off()
}

fb_all_enrich <- bind_rows(fb_enrich_combined)
write.csv(fb_all_enrich, './output/Xenium/fb_all_enrichr.csv', row.names = FALSE)


##############################################
#### Compensation vs decompensation signature score
##############################################
# Biological rationale:
#   The pseudobulk DESeq2 contrasts (fb_pseudobulk_deseq2.csv) reveal a striking two-phase
#   pattern across all 7 FB subtypes: genes that have been classically linked to
#   cardioprotection (ACE2, the counter-regulatory arm of the RAAS; IL1RL1/ST2, the
#   soluble IL-33 decoy receptor; GPX3, extracellular glutathione peroxidase; TIMP4,
#   inhibitor of matrix metalloproteinases) are progressively LOST from compensated
#   pRV into overt RVF, while a distinct fibrogenic programme (SFRP2/SFRP4 Wnt
#   modulators, BMP6, the POSTN-lineage matricellular genes TNC/THBS4/OGN, and
#   NTN4) is specifically GAINED at the pRV->RVF transition rather than already
#   present in pRV. Scoring each cell for both axes simultaneously lets us visualise
#   whether individual fibroblasts are walking a diagonal trajectory (losing
#   protection while gaining fibrogenesis) or whether the two transitions are
#   decoupled in time (protection collapses in pRV, fibrogenesis appears only at RVF).
#   References:
#     ACE2/RAAS counter-regulation:  Oudit et al. Circ Res 2003; Tikellis & Thomas 2012
#     IL-33/ST2 in HF:               Kakkar & Lee, Nat Rev Cardiol 2008
#     SFRP2 as fibrosis biomarker:   Mastri et al. AJP Heart Circ 2014; Lin et al. PNAS 2016
#     Matrifibrocyte programme:      Fu et al. JCI 2018 (TNC/THBS4/OGN)
#     POSTN lineage:                 Kanisicak et al. Nat Commun 2016

library(ggpubr)

compensation_genes   <- c('ACE2','IL1RL1','GPX3','TIMP4','CYP4B1','HMGCS2','PLA2G2A','EDNRB','CRISPLD2')
decompensation_genes <- c('SFRP2','SFRP4','BMP6','OGN','THBS4','TNC','NTN4','TCF4','ANKRD1','EDNRA')

# Keep only genes present in the (panel + imputed) Xenium assay
compensation_genes   <- intersect(compensation_genes,   rownames(fb.obj.clean.clean))
decompensation_genes <- intersect(decompensation_genes, rownames(fb.obj.clean.clean))
cat('Compensation genes kept:  ', paste(compensation_genes,   collapse=', '), '\n')
cat('Decompensation genes kept:', paste(decompensation_genes, collapse=', '), '\n')

DefaultAssay(fb.obj.clean.clean) <- 'Xenium'
fb.obj.clean.clean <- AddModuleScore(fb.obj.clean.clean,
                                     features = list(compensation_genes),
                                     name     = 'comp_score',
                                     assay    = 'Xenium')
fb.obj.clean.clean <- AddModuleScore(fb.obj.clean.clean,
                                     features = list(decompensation_genes),
                                     name     = 'decomp_score',
                                     assay    = 'Xenium')
fb.obj.clean.clean$comp_score1   <- scale(fb.obj.clean.clean$comp_score1)[,1]
fb.obj.clean.clean$decomp_score1 <- scale(fb.obj.clean.clean$decomp_score1)[,1]
fb.obj.clean.clean$delta_score   <- fb.obj.clean.clean$decomp_score1 - fb.obj.clean.clean$comp_score1

score_df <- fb.obj.clean.clean@meta.data[, c('fb_subtype','group','comp_score1','decomp_score1','delta_score')]
score_df$group <- factor(score_df$group, levels = c('NF','pRV','RVF'))

# --- Scatter: compensation x decompensation, faceted by subtype, colored by group ---
p_scatter <- ggplot(score_df,
                    aes(x = comp_score1, y = decomp_score1, color = group)) +
  geom_point(size = 0.15, alpha = 0.25) +
  stat_ellipse(level = 0.68, linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey50') +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey50') +
  facet_wrap(~ fb_subtype, nrow = 2) +
  scale_color_manual(values = c(NF = '#1f77b4', pRV = '#ff7f0e', RVF = '#d62728')) +
  labs(x = 'Compensation score (ACE2/IL1RL1/GPX3/TIMP4/...)',
       y = 'Decompensation score (SFRP2/SFRP4/BMP6/TNC/...)',
       title = 'Two-phase FB trajectory: protection loss precedes fibrogenic commitment') +
  theme_bw(base_size = 10) +
  theme(legend.position = 'bottom')

ggsave('./output/Xenium/fb_comp_decomp_scatter.pdf', p_scatter, width = 12, height = 6)

# --- Boxplot: delta score (decomp - comp) by group within each subtype ---
# Interpretation: increasingly positive delta = cell has shifted from protective toward fibrogenic state.
pairwise_comparisons <- list(c('NF','pRV'), c('NF','RVF'), c('pRV','RVF'))

p_delta <- ggplot(score_df, aes(x = group, y = delta_score, fill = group)) +
  geom_boxplot(outlier.size = 0.2, alpha = 0.8) +
  facet_wrap(~ fb_subtype, nrow = 2, scales = 'free_y') +
  stat_compare_means(comparisons = pairwise_comparisons, method = 'wilcox.test',
                     size = 2.5, label = 'p.signif') +
  scale_fill_manual(values = c(NF = '#1f77b4', pRV = '#ff7f0e', RVF = '#d62728')) +
  labs(x = NULL, y = 'Decompensation - Compensation (z)',
       title = 'Protective-to-fibrogenic shift per FB subtype across disease stage') +
  theme_bw(base_size = 10) +
  theme(legend.position = 'none')

ggsave('./output/Xenium/fb_comp_decomp_delta_boxplot.pdf', p_delta, width = 12, height = 6)

# --- Per-group, per-subtype mean scores for the manuscript table ---
score_summary <- score_df %>%
  group_by(fb_subtype, group) %>%
  summarise(mean_comp   = mean(comp_score1,   na.rm = TRUE),
            mean_decomp = mean(decomp_score1, na.rm = TRUE),
            mean_delta  = mean(delta_score,   na.rm = TRUE),
            n_cells     = dplyr::n(),
            .groups = 'drop')
write.csv(score_summary, './output/Xenium/fb_comp_decomp_summary.csv', row.names = FALSE)


##############################################
#### Wnt / ACE2 cross-subtype heatmap
##############################################
# Biological rationale:
#   The ACE2-RAAS counter-regulatory axis and the secreted Wnt modulator family
#   (SFRP2/SFRP4/DKK2) sit at two opposing nodes of the cardiac fibrosis network:
#     * ACE2 converts Ang II to Ang-(1-7), braking the pro-fibrotic AT1R arm of RAAS.
#       ACE2 knockout accelerates cardiac fibrosis (Oudit et al. Cardiovasc Res 2007),
#       and recombinant ACE2 reduces RV fibrosis in pulmonary hypertension models
#       (Shenoy et al. Am J Respir Crit Care Med 2010).
#     * Secreted Frizzled-Related Proteins (SFRP2, SFRP4) were initially described
#       as Wnt ANTAGONISTS but are now recognised as bona fide fibrosis EFFECTORS:
#       SFRP2 has BMP-like activity, drives myofibroblast persistence, and is a
#       validated serum biomarker of cardiac fibrosis (Mastri et al. 2014;
#       Kobayashi et al. Nat Cell Biol 2009; Lin et al. PNAS 2016). SFRP4 is in
#       the top panel of circulating fibrosis markers in HF cohorts.
#     * The endothelin axis (EDNRA pro-fibrotic, EDNRB protective) and TIMP4
#       (uniquely cardiac-enriched MMP inhibitor) complete the protective/
#       fibrogenic switchboard.
#   A subtype x gene log2FC heatmap, split by NF->pRV, pRV->RVF, NF->RVF, reveals
#   (a) whether protective gene loss is uniform across all FB states or concentrated
#   in specific subtypes, and (b) which subtype carries the SFRP2/SFRP4 signal —
#   i.e. which FB state is the cellular source of the circulating biomarker.

fb_deg <- read.csv('./output/Xenium/fb_pseudobulk_deseq2.csv', stringsAsFactors = FALSE)

ace2_axis_genes <- c('ACE2','IL1RL1','AGTR1','EDNRA','EDNRB','TIMP4','GPX3')
wnt_axis_genes  <- c('SFRP2','SFRP4','DKK2','DACT2','WNT5B','WT1','BMP6')
hm_genes <- c(ace2_axis_genes, wnt_axis_genes)

hm_df <- fb_deg %>%
  filter(gene %in% hm_genes) %>%
  mutate(axis    = ifelse(gene %in% ace2_axis_genes, 'ACE2 / protective', 'Wnt / fibrogenic'),
         lfc_cap = pmax(pmin(log2FoldChange, 3), -3),                         # cap for display
         sig     = ifelse(is.na(padj), '',
                    ifelse(padj < 0.001, '***',
                    ifelse(padj < 0.01,  '**',
                    ifelse(padj < 0.05,  '*', '')))),
         comparison = factor(comparison, levels = c('pRV_vs_NF','RVF_vs_pRV','RVF_vs_NF')),
         gene       = factor(gene, levels = hm_genes),
         axis       = factor(axis, levels = c('ACE2 / protective','Wnt / fibrogenic')))

p_heat <- ggplot(hm_df, aes(x = cell_type, y = gene, fill = lfc_cap)) +
  geom_tile(color = 'white', linewidth = 0.2) +
  geom_text(aes(label = sig), size = 3, vjust = 0.78) +
  facet_grid(axis ~ comparison, scales = 'free_y', space = 'free_y') +
  scale_fill_gradient2(low = '#2166ac', mid = 'white', high = '#b2182b',
                       midpoint = 0, limits = c(-3, 3),
                       name = expression(log[2]*' FC (capped)')) +
  labs(x = NULL, y = NULL,
       title = 'Protective (ACE2/RAAS) vs fibrogenic (Wnt) axes across FB subtypes and disease stages',
       subtitle = '*p<0.05  **p<0.01  ***p<0.001 (pseudobulk DESeq2)') +
  theme_bw(base_size = 10) +
  theme(axis.text.x      = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.y     = element_text(angle = 0))

ggsave('./output/Xenium/fb_wnt_ace2_heatmap.pdf', p_heat, width = 11, height = 6)


##############################################
#### FB subtype spatial organisation
####   A. FB-vs-FB co-localization
####   B. Transmural (endo <-> epi) preference
##############################################
# Rationale:
#   Pseudobulk DESeq2 has shown that the 7 FB subtypes carry distinct
#   programmes (protective loss vs fibrogenic commitment). The open question
#   is whether these subtypes are distributed uniformly through the RV wall
#   or whether they occupy specific spatial niches that explain their
#   transcriptional divergence.
#   Part A asks whether subtypes co-localize in tight clusters (same-subtype
#   enrichment on diagonal) or whether specific subtype pairs neighbour each
#   other more than chance (off-diagonal).
#   Part B asks whether any subtype has an anatomical bias along the
#   transmural (endocardium -> epicardium) axis, and whether that bias is
#   reshaped by disease stage.

library(RANN)
library(ggridges)

fb_groups_levels <- c('NF','pRV','RVF')
group_palette    <- c(NF = '#1f77b4', pRV = '#ff7f0e', RVF = '#d62728')
fb.obj.clean.clean$group <- factor(fb.obj.clean.clean$group, levels = fb_groups_levels)

# =========================================================
# PART A - FB subtype co-localization (patient-wise enrichment)
# =========================================================
# Reuses the existing compute_spatial_proximity() helper (line ~763) on an
# FB-only slice so that neighbour edges run exclusively between FB cells.
# Enrichment > 1 = subtype pair adjacent more often than a random shuffle
# of labels.

fb_levels <- levels(fb.obj.clean.clean$fb_subtype)
patients  <- unique(fb.obj.clean.clean$patient)

per_patient_enrich <- list()
for (pat in patients) {
  pat_cells <- WhichCells(fb.obj.clean.clean, expression = patient == pat)
  if (length(pat_cells) < 200) {
    cat(pat, ': too few FB cells (', length(pat_cells), ') - skipping\n', sep = '')
    next
  }
  sub  <- subset(fb.obj.clean.clean, cells = pat_cells)
  prox <- compute_spatial_proximity(sub, label_col = 'fb_subtype', n_neighbors = 10)
  enr  <- prox$enrichment

  # Pad to the common 7x7 frame (a patient may lack some rare subtypes)
  full <- matrix(NA_real_, nrow = length(fb_levels), ncol = length(fb_levels),
                 dimnames = list(fb_levels, fb_levels))
  common_r <- intersect(rownames(enr), fb_levels)
  common_c <- intersect(colnames(enr), fb_levels)
  full[common_r, common_c] <- enr[common_r, common_c]
  per_patient_enrich[[pat]] <- full
}

mean_enrich <- function(mats) {
  if (length(mats) == 0) return(NULL)
  arr <- simplify2array(mats)
  if (length(dim(arr)) == 2) return(arr)   # only 1 patient
  apply(arr, c(1, 2), mean, na.rm = TRUE)
}
pooled_enrich <- mean_enrich(per_patient_enrich)

patient_meta <- unique(fb.obj.clean.clean@meta.data[, c('patient','group')])
group_enrich <- lapply(fb_groups_levels, function(g) {
  pats_g <- as.character(patient_meta$patient[patient_meta$group == g])
  pats_g <- intersect(pats_g, names(per_patient_enrich))
  if (length(pats_g) == 0) return(NULL)
  mean_enrich(per_patient_enrich[pats_g])
})
names(group_enrich) <- fb_groups_levels

# ---- Heatmap: pooled + per-group ----
library(reshape2)
melt_mat <- function(m, tag) {
  if (is.null(m)) return(NULL)
  df <- melt(m, varnames = c('subtype_A','subtype_B'), value.name = 'enrichment')
  df$panel <- tag
  df
}
hm_df <- rbind(
  melt_mat(pooled_enrich, 'Pooled (all patients)'),
  do.call(rbind,
          Map(melt_mat, group_enrich, paste0('Group: ', names(group_enrich))))
)
hm_df$panel <- factor(hm_df$panel,
                      levels = c('Pooled (all patients)',
                                 paste0('Group: ', fb_groups_levels)))
hm_df$subtype_A <- factor(hm_df$subtype_A, levels = fb_levels)
hm_df$subtype_B <- factor(hm_df$subtype_B, levels = fb_levels)

p_coloc <- ggplot(hm_df,
                  aes(x = subtype_A, y = subtype_B,
                      fill = log2(pmax(enrichment, 1e-3)))) +
  geom_tile(color = 'white', linewidth = 0.2) +
  facet_wrap(~ panel, nrow = 1) +
  scale_fill_gradient2(low = '#2166ac', mid = 'white', high = '#b2182b',
                       midpoint = 0, limits = c(-2, 2), oob = scales::squish,
                       name = expression(log[2]*' obs/exp')) +
  labs(x = NULL, y = NULL,
       title = 'FB subtype co-localization (k=10 nearest-FB neighbour enrichment)') +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid  = element_blank())

ggsave('./output/Xenium/fb_subtype_colocalization.pdf',
       p_coloc, width = 16, height = 5)

write.csv(pooled_enrich,
          './output/Xenium/fb_subtype_colocalization_pooled.csv')
for (g in fb_groups_levels) {
  if (!is.null(group_enrich[[g]]))
    write.csv(group_enrich[[g]],
              paste0('./output/Xenium/fb_subtype_colocalization_', g, '.csv'))
}

# =========================================================
# PART B - Transmural preference (endocardium <-> epicardium)
# =========================================================
# Self-contained: compute depth locally for FB cells using Endo/Epi anchors
# from the parent xenium.obj. The script-wide transmural_depth block runs
# later, so we cannot rely on it here.
#
#   depth = d_endo / (d_endo + d_epi)
#     0 = on the endocardium, 1 = on the epicardium, 0.5 = mid-wall.

fb.obj.clean.clean$fb_transmural_depth <- NA_real_

for (pat in patients) {
  fb_cells_p <- WhichCells(fb.obj.clean.clean, expression = patient == pat)
  if (length(fb_cells_p) == 0) next

  parent_meta <- xenium.obj@meta.data
  endo_idx <- parent_meta$patient == pat &
              parent_meta$cell_type_rctd_doublet == 'Endo'
  epi_idx  <- parent_meta$patient == pat &
              parent_meta$cell_type_rctd_doublet == 'Epi'
  cat(pat, ': Endo =', sum(endo_idx), ' Epi =', sum(epi_idx), '\n')
  if (sum(endo_idx) < 5 || sum(epi_idx) < 5) {
    cat('  Skipping - insufficient Endo/Epi anchors\n'); next
  }

  endo_xy <- as.matrix(parent_meta[endo_idx, c('x_centroid','y_centroid')])
  epi_xy  <- as.matrix(parent_meta[epi_idx,  c('x_centroid','y_centroid')])
  fb_xy   <- as.matrix(fb.obj.clean.clean@meta.data[fb_cells_p,
                       c('x_centroid','y_centroid')])

  d_endo <- nn2(endo_xy, fb_xy, k = 1)$nn.dists[, 1]
  d_epi  <- nn2(epi_xy,  fb_xy, k = 1)$nn.dists[, 1]

  fb.obj.clean.clean$fb_transmural_depth[fb_cells_p] <-
    d_endo / (d_endo + d_epi)
}

depth_df <- fb.obj.clean.clean@meta.data[,
              c('fb_subtype','group','patient','fb_transmural_depth',
                'x_centroid','y_centroid')]
depth_df <- depth_df[!is.na(depth_df$fb_transmural_depth), ]
depth_df$group <- factor(depth_df$group, levels = fb_groups_levels)

# ---- B1. Pooled ridge plot ----
p_ridge_pool <- ggplot(depth_df,
                       aes(x = fb_transmural_depth, y = fb_subtype, fill = fb_subtype)) +
  geom_density_ridges(alpha = 0.8, scale = 1.1, linewidth = 0.3) +
  scale_x_continuous(limits = c(0, 1),
                     breaks = c(0, 0.33, 0.67, 1),
                     labels = c('Endo','','','Epi')) +
  labs(x = 'Transmural depth (0 = endocardium, 1 = epicardium)',
       y = NULL,
       title = 'FB subtype distribution across the RV wall') +
  theme_bw(base_size = 10) +
  theme(legend.position = 'none')
ggsave('./output/Xenium/fb_subtype_transmural_ridges.pdf',
       p_ridge_pool, width = 7, height = 5)

# ---- B2. Ridge plot by disease group ----
p_ridge_grp <- ggplot(depth_df,
                      aes(x = fb_transmural_depth, y = fb_subtype, fill = group)) +
  geom_density_ridges(alpha = 0.6, scale = 1.1, linewidth = 0.3) +
  scale_x_continuous(limits = c(0, 1),
                     breaks = c(0, 0.33, 0.67, 1),
                     labels = c('Endo','','','Epi')) +
  scale_fill_manual(values = group_palette) +
  facet_wrap(~ group, nrow = 1) +
  labs(x = 'Transmural depth',
       y = NULL,
       title = 'FB subtype transmural distribution across disease stages') +
  theme_bw(base_size = 10) +
  theme(legend.position = 'none')
ggsave('./output/Xenium/fb_subtype_transmural_ridges_by_group.pdf',
       p_ridge_grp, width = 12, height = 5)

# ---- B3. Summary table ----
depth_summary <- depth_df %>%
  mutate(third = cut(fb_transmural_depth,
                     breaks = c(-1e-6, 1/3, 2/3, 1 + 1e-6),
                     labels = c('endo','mid','epi'))) %>%
  group_by(fb_subtype, group) %>%
  summarise(n_cells      = dplyr::n(),
            mean_depth   = mean(fb_transmural_depth),
            median_depth = median(fb_transmural_depth),
            frac_endo    = mean(third == 'endo'),
            frac_mid     = mean(third == 'mid'),
            frac_epi     = mean(third == 'epi'),
            .groups = 'drop')
write.csv(depth_summary,
          './output/Xenium/fb_subtype_transmural_summary.csv',
          row.names = FALSE)

# ---- B4. Statistical tests ----
# (a) Within each group, are the 7 subtypes distributed differently across depth?
kw_per_group <- depth_df %>%
  group_by(group) %>%
  summarise(p_kruskal = kruskal.test(fb_transmural_depth ~ fb_subtype)$p.value,
            .groups = 'drop')
write.csv(kw_per_group,
          './output/Xenium/fb_subtype_transmural_kruskal.csv',
          row.names = FALSE)

# (b) Per subtype, does depth shift between disease groups?
wilcox_shift <- do.call(rbind, lapply(levels(depth_df$fb_subtype), function(st) {
  sub <- depth_df[depth_df$fb_subtype == st, ]
  pairs <- list(c('NF','pRV'), c('NF','RVF'), c('pRV','RVF'))
  do.call(rbind, lapply(pairs, function(pp) {
    x <- sub$fb_transmural_depth[sub$group == pp[1]]
    y <- sub$fb_transmural_depth[sub$group == pp[2]]
    if (length(x) < 20 || length(y) < 20) return(NULL)
    data.frame(fb_subtype = st,
               comparison = paste(pp, collapse = '_vs_'),
               p_wilcox   = suppressWarnings(wilcox.test(x, y)$p.value),
               delta_mean = mean(y) - mean(x))
  }))
}))
wilcox_shift$padj <- p.adjust(wilcox_shift$p_wilcox, method = 'BH')
write.csv(wilcox_shift,
          './output/Xenium/fb_subtype_transmural_group_shift.csv',
          row.names = FALSE)

# ---- B5. Spatial sanity scatter: FB coloured by subtype + Endo/Epi anchors ----
anchor_df <- xenium.obj@meta.data[
  xenium.obj$cell_type_rctd_doublet %in% c('Endo','Epi'),
  c('x_centroid','y_centroid','patient','cell_type_rctd_doublet')]
colnames(anchor_df)[4] <- 'anchor'

p_space <- ggplot() +
  geom_point(data = depth_df,
             aes(x = x_centroid, y = -y_centroid, color = fb_subtype),
             size = 0.15, alpha = 0.6) +
  geom_point(data = anchor_df,
             aes(x = x_centroid, y = -y_centroid, shape = anchor),
             size = 0.4, color = 'black') +
  facet_wrap(~ patient, nrow = 3, scales = 'free') +
  scale_shape_manual(values = c(Endo = 3, Epi = 4)) +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1)),
         shape = guide_legend(override.aes = list(size = 2))) +
  labs(title = 'FB subtypes in situ with Endo (+) / Epi (x) anchors',
       x = NULL, y = NULL) +
  theme_bw(base_size = 9) +
  theme(axis.text  = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank())
ggsave('./output/Xenium/fb_subtype_spatial_with_anchors.pdf',
       p_space, width = 14, height = 12)


##############################################
#### CM subclustering
##############################################
ref <- readRDS('../snRV_ref.rds')                          

cm.ref <- subset(ref, subset = Subnames_manual %in% c('Cm1', 'Cm2', 'Cm3', 'Cm4','Cm5','Cm6','Cm7','Cm8','Cm9','Cm10'))
# update Subnames after checking ref$Subnames

cm.obj <- subset(xenium.obj, subset = cell_type_rctd_doublet == 'CM')
cm.obj <- JoinLayers(cm.obj)

umap_cols <- grep("^umap_\\d+$|^UMAP_\\d+$|^umaporig_\\d+$",
                  colnames(cm.obj@meta.data), value = TRUE)
if (length(umap_cols)) cm.obj@meta.data[umap_cols] <- NULL

cm.obj[["pca"]]       <- NULL
cm.obj[["harmony"]]   <- NULL
cm.obj[["umap"]]      <- NULL
cm.obj[["umap_orig"]] <- NULL

cm.obj <- FindVariableFeatures(cm.obj)
cm.obj <- ScaleData(cm.obj)
cm.obj <- RunPCA(cm.obj, npcs = 30)
cm.obj <- RunHarmony(cm.obj, "patient")

cm.obj <- RunUMAP(cm.obj, reduction = "harmony", dims = 1:30, reduction.name = 'umap')
cm.obj <- RunUMAP(cm.obj, reduction = "pca",     dims = 1:30, reduction.name = 'umap_orig')

cm.obj <- FindNeighbors(cm.obj, reduction = "harmony", dims = 1:30)
cm.obj <- FindClusters(cm.obj, resolution = 0.5)

cm.obj.mem <- cm.obj
cm.obj.mem[["Xenium"]] <- CreateAssayObject(
  counts = as(GetAssayData(cm.obj, layer = "counts"), "dgCMatrix")
)
cm.marks <- FindAllMarkers(cm.obj.mem)
rm(cm.obj.mem)
write.csv(cm.marks, './output/Xenium/cm.marks.csv')

DimPlot(cm.obj, group.by = 'seurat_clusters', label = TRUE, repel = TRUE, reduction = 'umap') +
  theme(legend.position = 'none')
DimPlot(cm.obj, group.by = 'patient', reduction = 'umap')

# ── Annotate clusters after reviewing cm.marks.csv + canonical marker DotPlot ──
#
#   ┌─────────┬──────────────────────────┬──────────────────────────────────────┐
#   │ Cluster │       Cell Type          │           Key Markers                │
#   ├─────────┼──────────────────────────┼──────────────────────────────────────┤
#   │ 0       │ CM_Metabolic             │ COQ10A, ECSIT, DLAT                  │
#   │ 1       │ CM_Stressed              │ TNNI3, THBS4, IGF1, FBN1            │
#   │ 2       │ Fibroblast               │ DCN, TNC, ELN, ASPN, VCAN           │
#   │ 3       │ CM_Ventricular           │ HVCN1, CPNE5, EYA4, SLC22A5         │
#   │ 4       │ CM_Contractile           │ TNNT2, MYH6, NPPB, MYBPC3           │
#   │ 5       │ Endothelial              │ CD34, CLEC14A, BTNL9, DLL1           │
#   │ 6       │ CM_Ventricular_2         │ GALNTL6, EXPH5                       │
#   │ 7       │ Macrophage               │ LYVE1, IRF8, MS4A4A, MRC1            │
#   │ 8       │ Fibroblast_2             │ ANKRD1, COL6A2, COL3A1               │
#   │ 9       │ Fibroblast_3             │ COL6A2, COL3A1, COL1A2, COL6A1      │
#   │ 10      │ Endothelial_2            │ STC1, CDH5, TIE1, TEK, EMCN          │
#   │ 11      │ Pericyte                 │ COL27A1, RGS5, GUCY1A2              │
#   │ 12      │ Fibroblast_4             │ COL1A1, COL6A2, COL6A1, HEY2        │
#   │ 13      │ Mast_Cell                │ CTSG, GATA2, MS4A2, HDC, CPA3, KIT  │
#   │ 14      │ CM_Ventricular_4         │ ALS2CL, HVCN1, EPCAM                │
#   │ 15      │ T_NK_Cell                │ IL7R, PRF1, ZAP70, CD48              │
#   │ 16      │ Dendritic                │ HLA-DQB1, ITGAX, TNFRSF13B          │
#   │ 17      │ Adipocyte                │ PLIN4, ADIPOQ, TIMP4, LGALS12       │
#   │ 18      │ Macrophage_2             │ CSF1R, MERTK, SIGLEC1                │
#   │ 19      │ Fibroblast_5             │ COL1A1, COL6A1, COL1A2, FAP         │
#   │ 20      │ Lymphatic_EC             │ MMRN1, CCL21                         │
#   │ 21      │ Schwann                  │ SOX10, ERBB3, SLC5A7                 │
#   │ 22      │ Proliferating            │ TPX2, KIF11, BUB1B, NCAPG           │
#   │ 23      │ CM_Ventricular_5         │ ABRA, NPR3, ABCG2                   │
#   └─────────┴──────────────────────────┴──────────────────────────────────────┘

cm_labels <- c(
  '0'  = 'CM_Metabolic',
  '1'  = 'CM_Stressed',
  '2'  = 'Fibroblast',
  '3'  = 'CM_Ventricular',
  '4'  = 'CM_Contractile',
  '5'  = 'Endothelial',
  '6'  = 'CM_Ventricular_2',
  '7'  = 'Macrophage',
  '8'  = 'Fibroblast_2',
  '9'  = 'Fibroblast_3',
  '10' = 'Endothelial_2',
  '11' = 'Pericyte',
  '12' = 'Fibroblast_4',
  '13' = 'Mast_Cell',
  '14' = 'CM_Ventricular_4',
  '15' = 'T_NK_Cell',
  '16' = 'Dendritic',
  '17' = 'Adipocyte',
  '18' = 'Macrophage_2',
  '19' = 'Fibroblast_5',
  '20' = 'Lymphatic_EC',
  '21' = 'Schwann',
  '22' = 'Proliferating',
  '23' = 'CM_Ventricular_5'
)

cm.obj$cm_subtype <- unname(cm_labels[as.character(cm.obj$seurat_clusters)])

DimPlot(cm.obj.clean, group.by = 'cm_subtype', label = TRUE, repel = TRUE, reduction = 'umap_orig') +
  theme(legend.position = 'none')

cm.obj.clean <- subset(cm.obj, subset = cm_subtype %in% c(
  'CM_Metabolic', 'CM_Stressed', 'CM_Ventricular', 'CM_Contractile',
  'CM_Ventricular_2', 'CM_Ventricular_4', 'CM_Ventricular_5'))

cm.obj.clean[["pca"]]       <- NULL
cm.obj.clean[["harmony"]]   <- NULL
cm.obj.clean[["umap"]]      <- NULL
cm.obj.clean[["umap_orig"]] <- NULL

umap_cols <- grep("^umap_\\d+$|^UMAP_\\d+$|^umaporig_\\d+$",
                  colnames(cm.obj.clean@meta.data), value = TRUE)
if (length(umap_cols)) cm.obj.clean@meta.data[umap_cols] <- NULL

cm.obj.clean <- JoinLayers(cm.obj.clean)
cm.obj.clean <- FindVariableFeatures(cm.obj.clean)
cm.obj.clean <- ScaleData(cm.obj.clean)
cm.obj.clean <- RunPCA(cm.obj.clean, npcs = 30)
cm.obj.clean <- RunHarmony(cm.obj.clean, "patient")
cm.obj.clean <- RunUMAP(cm.obj.clean, reduction = "harmony", dims = 1:30, reduction.name = 'umap')
cm.obj.clean <- RunUMAP(cm.obj.clean, reduction = "pca",     dims = 1:30, reduction.name = 'umap_orig')

cm.obj.clean <- FindNeighbors(cm.obj.clean, reduction = "harmony", dims = 1:30)
cm.obj.clean <- FindClusters(cm.obj.clean, resolution = 0.5)

cm.obj.clean.mem <- cm.obj.clean
cm.obj.clean.mem[["Xenium"]] <- CreateAssayObject(
  counts = as(GetAssayData(cm.obj.clean, layer = "counts"), "dgCMatrix")
)
cm.clean.marks <- FindAllMarkers(cm.obj.clean.mem)
rm(cm.obj.clean.mem)
write.csv(cm.clean.marks, './output/Xenium/cm.clean.marks.csv')

# ── Annotate CM clean clusters ───────────────────────────────────────────────
#
#   ┌─────────┬────────────────────┬──────────────────────────────────────────┐
#   │ Cluster │     Cell Type      │             Key Markers                  │
#   ├─────────┼────────────────────┼──────────────────────────────────────────┤
#   │ 0       │ CM_Ventricular     │ PKP2, ATP2A2, ACTN2, PLN                │
#   │ 1       │ CM_Ventricular_2   │ RYR2, ATP2A2, TTN; CXCL8, FGL2         │
#   │ 2       │ CM_RV              │ HAND2, MYL7, CASQ2, GATA4, TBX5        │
#   │ 3       │ CM_RV_2            │ HAND2, MYL7, CASQ2, TBX5, DSP          │
#   │ 4       │ CM_Ventricular_3   │ DSP, TBX5, PLN, GRIK1                  │
#   │ 5       │ CM_Ventricular_4   │ RYR2, ATP2A2, TTN, SCN5A               │
#   │ 6       │ CM_RV_3            │ MYL7, CASQ2, SCN5A, TBX5, ACTN2       │
#   │ 7       │ CM_Ventricular_5   │ HEY2, ATP2A2, RYR2, PLN, DSP           │
#   │ 8       │ CM_Contractile     │ NPPA, TNNI3, DES, TNNT2, ACTC1, MYH6   │
#   │ 9       │ CM_Ventricular_6   │ TBX5, ACTN2, CASQ2, RYR2, DSP          │
#   │ 10      │ CM_Ventricular_7   │ CASQ2, SCN5A, TBX5, ACTN2, MYL3        │
#   │ 11      │ CM_Stressed        │ MYLK, IL1RL1, FBN1, BMP4; PKP2, DSP    │
#   └─────────┴────────────────────┴──────────────────────────────────────────┘

cm_clean_labels <- c(
  '0'  = 'CM_Ventricular',
  '1'  = 'CM_Ventricular_2',
  '2'  = 'CM_RV',
  '3'  = 'CM_RV_2',
  '4'  = 'CM_Ventricular_3',
  '5'  = 'CM_Ventricular_4',
  '6'  = 'CM_RV_3',
  '7'  = 'CM_Ventricular_5',
  '8'  = 'CM_Contractile',
  '9'  = 'CM_Ventricular_6',
  '10' = 'CM_Ventricular_7',
  '11' = 'CM_Stressed'
)

cm.obj.clean$cm_subtype <- unname(cm_clean_labels[as.character(cm.obj.clean$seurat_clusters)])

DimPlot(cm.obj.clean, group.by = 'cm_subtype', label = TRUE, repel = TRUE, reduction = 'uamp') +
  theme(legend.position = 'none')

cm.obj.clean.clean <- cm.obj.clean  # all clusters are true CMs — no further filtering

DimPlot(cm.obj.clean.clean, group.by = 'cm_subtype', label = TRUE, repel = TRUE, reduction = 'umap') +
  theme(legend.position = 'none')
DimPlot(cm.obj.clean.clean, group.by = 'patient', reduction = 'umap')

saveRDS(cm.obj.clean.clean, './output/Xenium/cm_clean_clean.rds')
#156866 nuclei

# Label transfer from snRV reference
cm.ref <- JoinLayers(cm.ref)
cm.anchors <- FindTransferAnchors(
  reference           = cm.ref,
  query               = cm.obj.clean.clean,
  dims                = 1:30,
  reference.reduction = 'pca'
)
cm.predictions <- TransferData(
  anchorset = cm.anchors,
  refdata   = cm.ref$Subnames_manual,
  dims      = 1:30
)
cm.obj.clean.clean <- AddMetaData(cm.obj.clean.clean, metadata = cm.predictions)

DimPlot(cm.obj.clean.clean, group.by = 'predicted.id', label = TRUE, repel = TRUE, reduction = 'umap_orig') +
  theme(legend.position = 'none') +
  ggtitle('Reference predicted labels')

table(cm.obj.clean.clean$predicted.id, cm.obj.clean.clean$cm_subtype)

saveRDS(cm.obj.clean.clean, './output/Xenium/cm_clean_clean.rds')

pdf('./output/Xenium/Xenium_cm_snUMAP_new.pdf', width = 5, height = 5)
PlotEmbedding(cm.obj.clean.clean, group.by = 'cm_subtype', reduction = 'umap', point_size = .05,
              plot_under = TRUE, plot_theme = umap_theme() + NoLegend(),
              raster_dpi = 400, raster_scale = 0.5)
dev.off()

pdf('./output/Xenium/Xenium_cm_snUMAP_new_pt.pdf', width = 5, height = 5)
PlotEmbedding(cm.obj.clean.clean, group.by = 'patient', reduction = 'umap', point_size = .05,
              plot_under = TRUE, plot_theme = umap_theme() + NoLegend(),
              raster_dpi = 400, raster_scale = 0.5)
dev.off()

pdf('./output/Xenium/cm_dot.pdf', width = 14, height = 5)
DotPlot(cm.obj.clean.clean,
        features = c(
          # CM_Contractile — subendocardial natriuretic peptide-high
          "NPPA", "NPPB", "MYH6", "ANKRD1",
          # CM_RV — HAND2+ right ventricular identity
          "HAND2", "MYL7", "GATA4",
          # CM_RV_2 — CACNA1E+/FGF7+/NREP+
          "CACNA1E", "FGF7", "NREP",
          # CM_RV_3 — DOT1L+/STEAP4+/EEF1A2+
          "DOT1L", "STEAP4", "EEF1A2",
          # CM_Stressed — MYLK, IL1RL1, FBN1, C7
          "MYLK", "IL1RL1", "FBN1", "C7",
          # CM_Ventricular — TIMP4+/PDCD4+/SDHD+
          "TIMP4", "PDCD4", "SDHD",
          # CM_Ventricular_2 — CXCL8+/MEDAG+/SLC45A4+
          "CXCL8", "MEDAG", "SLC45A4",
          # CM_Ventricular_3 — PRDX6+/ESR2+/SLC19A2+
          "PRDX6", "ESR2", "SLC19A2",
          # CM_Ventricular_4 — GUCY1A2+/EGFLAM+/NFASC+
          "GUCY1A2", "EGFLAM", "NFASC",
          # CM_Ventricular_5 — HEY2+/C1QTNF1+/CAMK1D+
          "HEY2", "C1QTNF1", "CAMK1D",
          # CM_Ventricular_6 — CELSR1+/LPCAT4+/SPTBN4+
          "CELSR1", "LPCAT4", "SPTBN4",
          # CM_Ventricular_7 — HPR+/FGF13+/ABRA+/NPR3+
          "HPR", "FGF13", "ABRA", "NPR3",
          # Pan-CM sarcomeric (last)
          "MYH7", "TNNT2", "TNNI3", "ACTC1", "MYL2", "MYL3", "DES", "TTN", "RYR2"
        ),
        group.by = "cm_subtype", col.min = 0, col.max = 2) +
  RotatedAxis() +
  theme(axis.text.x = element_text(size = 9)) +
  ggtitle("Canonical CM subtype markers")
dev.off()


##############################################
#### CM proportion by disease state
##############################################

pairwise_comparisons <- list(c('NF', 'pRV'), c('NF', 'RVF'), c('pRV', 'RVF'))

cm_meta <- cm.obj.clean.clean@meta.data %>%
  select(patient, group, cm_subtype)

cm_totals <- xenium.obj@meta.data %>%
  group_by(patient, group) %>%
  summarise(total = n(), .groups = 'drop')

cm_prop <- cm_meta %>%
  group_by(patient, group, cm_subtype) %>%
  summarise(n = n(), .groups = 'drop') %>%
  left_join(cm_totals, by = c('patient', 'group')) %>%
  mutate(proportion = n / total) %>%
  complete(nesting(patient, group), cm_subtype, fill = list(n = 0, proportion = 0)) %>%
  mutate(group = factor(group, levels = c('NF', 'pRV', 'RVF')))

p_cm_prop <- ggplot(cm_prop, aes(x = group, y = proportion, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.5) +
  stat_compare_means(comparisons = pairwise_comparisons, method = 'wilcox.test', label = 'p.signif') +
  stat_compare_means(method = 'kruskal.test', label.y.npc = 0.95, size = 3) +
  facet_wrap(~ cm_subtype, scales = 'free_y') +
  scale_x_discrete(limits = c('NF', 'pRV', 'RVF')) +
  theme_bw() +
  theme(legend.position = 'none',
        strip.text = element_text(face = 'bold'),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = 'Proportion of all cells', title = 'CM subtype distribution by group')

pdf('./output/Xenium/cm_proportions.pdf', width = 8, height = 6)
print(p_cm_prop)
dev.off()

##############################################
#### CM fractional area coverage by disease state
##############################################

cell_areas <- read.csv('./output/Xenium/cell_areas_2d.csv', row.names = 1)
total_area <- xenium.obj@meta.data %>%
  tibble::rownames_to_column('cell') %>%
  select(cell, patient, group) %>%
  left_join(cell_areas, by = 'cell') %>%
  group_by(patient, group) %>%
  summarise(total_area = sum(cell_area_2d, na.rm = TRUE), .groups = 'drop')

cm_area_meta <- cm.obj.clean.clean@meta.data %>%
  tibble::rownames_to_column('cell') %>%
  select(cell, patient, group, cm_subtype) %>%
  left_join(cell_areas, by = 'cell')

cm_area_prop <- cm_area_meta %>%
  group_by(patient, group, cm_subtype) %>%
  summarise(cm_area = sum(cell_area_2d, na.rm = TRUE), .groups = 'drop') %>%
  left_join(total_area, by = c('patient', 'group')) %>%
  mutate(frac_area = cm_area / total_area) %>%
  complete(nesting(patient, group), cm_subtype, fill = list(cm_area = 0, frac_area = 0)) %>%
  mutate(group = factor(group, levels = c('NF', 'pRV', 'RVF')))

p_cm_area <- ggplot(cm_area_prop, aes(x = group, y = frac_area, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.5) +
  stat_compare_means(comparisons = pairwise_comparisons, method = 'wilcox.test', label = 'p.signif') +
  stat_compare_means(method = 'kruskal.test', label.y.npc = 0.95, size = 3) +
  facet_wrap(~ cm_subtype, scales = 'free_y') +
  scale_x_discrete(limits = c('NF', 'pRV', 'RVF')) +
  theme_bw() +
  theme(legend.position = 'none',
        strip.text = element_text(face = 'bold'),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = 'Fractional area coverage', title = 'CM subtype fractional area by group')

pdf('./output/Xenium/cm_area_coverage.pdf', width = 10, height = 8)
print(p_cm_area)
dev.off()


##############################################
#### Pseudobulk DESeq2 by CM subtype
##############################################

cm_pseudo_meta <- unique(cm.obj.clean.clean@meta.data[, c('patient', 'group', 'cm_subtype')])
cm_pseudo_meta$patient    <- as.character(cm_pseudo_meta$patient)
cm_pseudo_meta$cm_subtype <- as.character(cm_pseudo_meta$cm_subtype)
cm_pseudo_meta$pseudobulk_id <- paste0(cm_pseudo_meta$patient, '..', cm_pseudo_meta$cm_subtype)
cm_pseudo_meta$group <- factor(cm_pseudo_meta$group, levels = c('NF', 'pRV', 'RVF'))
rownames(cm_pseudo_meta) <- cm_pseudo_meta$pseudobulk_id

cm.obj.clean.clean$pseudobulk_id <- paste0(
  as.character(cm.obj.clean.clean$patient), '..',
  as.character(cm.obj.clean.clean$cm_subtype)
)

cm_pseudo_counts <- AggregateExpression(
  cm.obj.clean.clean,
  assays = 'Xenium',
  group.by = 'pseudobulk_id',
  slot = 'counts',
  return.seurat = FALSE
)$Xenium

colnames(cm_pseudo_counts) <- gsub('-', '_', sub('^g', '', colnames(cm_pseudo_counts)))

cm_cell_types <- unique(cm_pseudo_meta$cm_subtype)
cm_deseq_results <- list()

for (ct in cm_cell_types) {
  ct_samples <- cm_pseudo_meta$pseudobulk_id[cm_pseudo_meta$cm_subtype == ct]
  ct_samples <- ct_samples[ct_samples %in% colnames(cm_pseudo_counts)]
  if (length(ct_samples) < 2) next

  ct_counts <- cm_pseudo_counts[, ct_samples, drop = FALSE]
  ct_meta   <- cm_pseudo_meta[ct_samples, , drop = FALSE]

  if (length(unique(ct_meta$group)) < 2) next

  keep <- rowSums(ct_counts >= 5) >= 2
  ct_counts <- ct_counts[keep, , drop = FALSE]
  if (nrow(ct_counts) < 10) next

  dds <- DESeqDataSetFromMatrix(countData = ct_counts, colData = ct_meta, design = ~ group)
  dds <- tryCatch(DESeq(dds), error = function(e) NULL)
  if (is.null(dds)) next

  for (comp in list(c('RVF', 'NF'), c('pRV', 'NF'), c('RVF', 'pRV'))) {
    comp_label <- paste0(comp[1], '_vs_', comp[2])
    res <- tryCatch(
      results(dds, contrast = c('group', comp[1], comp[2])) %>%
        as.data.frame() %>%
        tibble::rownames_to_column('gene') %>%
        mutate(cell_type = ct, comparison = comp_label) %>%
        arrange(padj),
      error = function(e) NULL
    )
    if (!is.null(res)) cm_deseq_results[[paste0(ct, '_', comp_label)]] <- res
  }
}

cm_all_deseq <- bind_rows(cm_deseq_results)
write.csv(cm_all_deseq, './output/Xenium/cm_pseudobulk_deseq2.csv', row.names = FALSE)

cm_sig_summary <- cm_all_deseq %>%
  filter(padj < 0.05, abs(log2FoldChange) > 0.5) %>%
  group_by(cell_type, comparison) %>%
  summarise(n_up = sum(log2FoldChange > 0), n_down = sum(log2FoldChange < 0), .groups = 'drop')
print(cm_sig_summary)


##############################################
#### CM subtypes: Volcano + GO/Reactome
##############################################

cap_pval <- function(df) mutate(df, padj = ifelse(is.na(padj), NA_real_, pmax(padj, 1e-5)))

make_volcano_colors <- function(df, panel_genes) {
  is_panel <- df$gene %in% panel_genes
  is_sig   <- !is.na(df$padj) & df$padj < 0.05 & abs(df$log2FoldChange) > 0.5
  cols <- case_when(
    is_sig & is_panel  ~ 'red3',
    is_sig & !is_panel ~ 'royalblue3',
    !is_sig & is_panel ~ 'salmon',
    TRUE               ~ 'grey80'
  )
  names(cols) <- df$gene
  cols
}

sig_labs <- function(df, n = 30) {
  df %>%
    filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 0.5) %>%
    arrange(padj) %>%
    head(n) %>%
    pull(gene)
}

volcano_args <- list(
  x = 'log2FoldChange', y = 'padj',
  pCutoff = 0.05, FCcutoff = 0.5,
  pointSize = 2, labSize = 3.5,
  drawConnectors = TRUE, widthConnectors = 0.5,
  xlim = c(-5, 5), ylim = c(0, 5),
  subtitle = 'Red = original panel gene | Blue = imputed gene',
  legendPosition = 'none'
)

make_volcano_pdf <- function(dat, filepath, title_str) {
  if (nrow(dat) == 0) { cat('  Skipping', title_str, '— no data\n'); return(invisible(NULL)) }
  pdf(filepath, width = 10, height = 8)
  tryCatch(
    print(do.call(EnhancedVolcano, c(list(toptable = dat, lab = dat$gene,
      title = title_str,
      selectLab = sig_labs(dat),
      colCustom = make_volcano_colors(dat, panel_genes)), volcano_args))),
    error = function(e) { cat('  Volcano error for', title_str, ':', conditionMessage(e), '\n') }
  )
  dev.off()
}

enrichR_dbs <- c('GO_Biological_Process_2023', 'Reactome_2022')
setEnrichrSite('Enrichr')

combine_enrich <- function(enrich_res, direction, comparison) {
  if (is.null(enrich_res)) return(NULL)
  lapply(names(enrich_res), function(db) {
    df <- enrich_res[[db]]
    if (nrow(df) == 0) return(NULL)
    df %>%
      mutate(Term = as.character(Term),
             database = db, direction = direction, comparison = comparison)
  }) %>% bind_rows()
}

plot_enrich <- function(enrich_df, title_str) {
  if (is.null(enrich_df) || nrow(enrich_df) == 0) return(NULL)
  enrich_df %>%
    arrange(Adjusted.P.value) %>%
    head(10) %>%
    mutate(Term = factor(Term, levels = rev(Term))) %>%
    ggplot(aes(x = -log10(Adjusted.P.value), y = Term)) +
    geom_col(fill = 'steelblue') +
    theme_bw() +
    labs(x = '-log10(adj. P)', y = NULL, title = title_str)
}

cm_enrich_combined <- list()

for (ct in unique(cm_all_deseq$cell_type)) {
  ct_slug <- gsub('_', '-', ct)
  cat('CM processing:', ct, '\n')

  dat_rvf     <- cm_all_deseq %>% filter(cell_type == ct, comparison == 'RVF_vs_NF')  %>% cap_pval()
  dat_prv     <- cm_all_deseq %>% filter(cell_type == ct, comparison == 'pRV_vs_NF')  %>% cap_pval()
  dat_rvf_prv <- cm_all_deseq %>% filter(cell_type == ct, comparison == 'RVF_vs_pRV') %>% cap_pval()

  make_volcano_pdf(dat_rvf,     paste0('./output/Xenium/cm_', ct_slug, '_volcano_RVF_vs_NF.pdf'),  paste(ct, ': RVF vs NF'))
  make_volcano_pdf(dat_prv,     paste0('./output/Xenium/cm_', ct_slug, '_volcano_pRV_vs_NF.pdf'),  paste(ct, ': pRV vs NF'))
  make_volcano_pdf(dat_rvf_prv, paste0('./output/Xenium/cm_', ct_slug, '_volcano_RVF_vs_pRV.pdf'), paste(ct, ': RVF vs pRV'))

  run_enrich <- function(dat, comp) {
    sig  <- dat %>% filter(padj < 0.05, abs(log2FoldChange) > 0.5)
    up   <- sig %>% filter(log2FoldChange > 0) %>% pull(gene)
    down <- sig %>% filter(log2FoldChange < 0) %>% pull(gene)
    list(
      up   = if (length(up)   > 0) enrichr(up,   enrichR_dbs) else NULL,
      down = if (length(down) > 0) enrichr(down, enrichR_dbs) else NULL
    )
  }

  er_rvf     <- run_enrich(dat_rvf,     'RVF_vs_NF')
  er_prv     <- run_enrich(dat_prv,     'pRV_vs_NF')
  er_rvf_prv <- run_enrich(dat_rvf_prv, 'RVF_vs_pRV')

  ct_enrich <- bind_rows(
    combine_enrich(er_rvf$up,       'up',   'RVF_vs_NF'),
    combine_enrich(er_rvf$down,     'down', 'RVF_vs_NF'),
    combine_enrich(er_prv$up,       'up',   'pRV_vs_NF'),
    combine_enrich(er_prv$down,     'down', 'pRV_vs_NF'),
    combine_enrich(er_rvf_prv$up,   'up',   'RVF_vs_pRV'),
    combine_enrich(er_rvf_prv$down, 'down', 'RVF_vs_pRV')
  )
  if (!is.null(ct_enrich) && nrow(ct_enrich) > 0) {
    ct_enrich$cell_type <- ct
    cm_enrich_combined[[ct]] <- ct_enrich
  }

  pdf(paste0('./output/Xenium/cm_', ct_slug, '_enrichr_plots.pdf'), width = 12, height = 7)
  for (comp_label in c('RVF_vs_NF', 'pRV_vs_NF', 'RVF_vs_pRV')) {
    er <- switch(comp_label, RVF_vs_NF = er_rvf, pRV_vs_NF = er_prv, RVF_vs_pRV = er_rvf_prv)
    for (dir_label in c('up', 'down')) {
      er_dir <- er[[dir_label]]
      if (!is.null(er_dir)) {
        p1 <- plot_enrich(er_dir$GO_Biological_Process_2023, paste(ct, comp_label, toupper(dir_label), 'GO BP'))
        p2 <- plot_enrich(er_dir$Reactome_2022,              paste(ct, comp_label, toupper(dir_label), 'Reactome'))
        if (!is.null(p1)) print(p1)
        if (!is.null(p2)) print(p2)
      }
    }
  }
  dev.off()
}

cm_all_enrich <- bind_rows(cm_enrich_combined)
write.csv(cm_all_enrich, './output/Xenium/cm_all_enrichr.csv', row.names = FALSE)





##############################################
#### 9. Nuclei distribution
##############################################

library(tidyr)
library(ggpubr)

nuclei_meta <- list(
  list(file = 'output-XETG00217__0038213__Region_1__20241206__182124_top_nuclei_assignments.csv',    patient = '1697', group = 'NF'),
  list(file = 'output-XETG00217__0038213__Region_1__20241206__182124_middle_nuclei_assignments.csv', patient = '1691', group = 'NF'),
  list(file = 'output-XETG00217__0038213__Region_1__20241206__182124_bottom_nuclei_assignments.csv', patient = '1618', group = 'pRV'),
  list(file = 'output-XETG00217__0038216__Region_1__20241206__182124_top_nuclei_assignments.csv',    patient = '1567', group = 'pRV'),
  list(file = 'output-XETG00217__0038216__Region_1__20241206__182124_bottom_nuclei_assignments.csv', patient = '1692', group = 'pRV'),
  list(file = 'output-XETG00217__0038290__Region_1__20241212__142808_nuclei_assignments.csv',        patient = '1561', group = 'NF'),
  list(file = 'output-XETG00217__0038290__Region_2__20241212__142808_nuclei_assignments.csv',        patient = '1343', group = 'RVF'),
  list(file = 'output-XETG00217__0038291__Region_1__20241212__142808_nuclei_assignments.csv',        patient = '1467', group = 'RVF'),
  list(file = 'output-XETG00217__0038291__Region_2__20241212__142808_nuclei_assignments.csv',        patient = '1632', group = 'RVF')
)

nuclei_raw <- bind_rows(lapply(nuclei_meta, function(m) {
  df <- read.csv(m$file)
  df$patient <- m$patient
  df$group   <- m$group
  df
}))

# Compute proportions (Unassigned excluded from denominator and plot)
nuclei_assigned <- nuclei_raw[!nuclei_raw$transferred_identity %in% c('Unassigned', 'None') &
                               !is.na(nuclei_raw$transferred_identity), ]

nuclei_totals <- nuclei_assigned %>%
  group_by(patient, group) %>%
  summarise(total = n(), .groups = 'drop')

nuclei_counts <- nuclei_assigned %>%
  group_by(patient, group, transferred_identity) %>%
  summarise(n = n(), .groups = 'drop') %>%
  left_join(nuclei_totals, by = c('patient', 'group')) %>%
  mutate(proportion = n / total)

# Fill missing cell types with 0 for each patient
nuclei_prop <- nuclei_counts %>%
  complete(nesting(patient, group), transferred_identity, fill = list(n = 0, proportion = 0)) %>%
  mutate(group = factor(group, levels = c('NF', 'pRV', 'RVF')))

pairwise_comparisons <- list(c('NF', 'pRV'), c('NF', 'RVF'), c('pRV', 'RVF'))

p_nuclei_dist <- ggplot(nuclei_prop, aes(x = group, y = proportion, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.5) +
  stat_compare_means(comparisons = pairwise_comparisons, method = 'wilcox.test', label = 'p.signif') +
  stat_compare_means(method = 'kruskal.test', label.y.npc = 0.95, size = 3) +
  facet_wrap(~ transferred_identity, scales = 'free_y') +
  theme_bw() +
  theme(legend.position = 'none',
        strip.text = element_text(face = 'bold'),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = 'Proportion of nuclei', title = 'Nuclei cell-type distribution by group')

pdf('./output/Xenium/nuclei_distribution.pdf', width = 6, height = 6)
print(p_nuclei_dist)
dev.off()

# Nuclei per proseg cell by disease state
nuclei_per_cell <- bind_rows(lapply(nuclei_meta, function(m) {
  rev_file <- sub('\\.csv$', '_rev.csv', m$file)
  rev_df <- read.csv(rev_file)
  rev_df$n_nuclei <- lengths(strsplit(as.character(rev_df$enclosed_nuclei), ';'))

  # Get cell type for each proseg cell from the forward CSV
  fwd_df <- read.csv(m$file)
  cell_types <- fwd_df %>%
    filter(!is.na(reference_cell_id)) %>%
    distinct(reference_cell_id, transferred_identity) %>%
    filter(!transferred_identity %in% c('Unassigned', 'None'))

  rev_df <- rev_df %>%
    left_join(cell_types, by = 'reference_cell_id') %>%
    filter(!is.na(transferred_identity))

  rev_df$patient <- m$patient
  rev_df$group   <- m$group
  rev_df
})) %>%
  mutate(group = factor(group, levels = c('NF', 'pRV', 'RVF')))

nuclei_per_cell_pb <- nuclei_per_cell %>%
  group_by(patient, group, transferred_identity) %>%
  summarise(mean_n_nuclei = mean(n_nuclei), .groups = 'drop')

p_nuclei_per_cell <- ggplot(nuclei_per_cell_pb, aes(x = group, y = mean_n_nuclei, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.5) +
  stat_compare_means(comparisons = pairwise_comparisons, method = 'wilcox.test', label = 'p.signif') +
  stat_compare_means(method = 'kruskal.test', label.y.npc = 0.95, size = 3) +
  facet_wrap(~ transferred_identity, scales = 'free_y') +
  theme_bw() +
  theme(legend.position = 'none',
        strip.text = element_text(face = 'bold'),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = 'Mean nuclei per proseg cell', title = 'Nuclei per proseg cell by disease state and cell type')

pdf('./output/Xenium/nuclei_per_cell.pdf', width = 6, height = 6)
print(p_nuclei_per_cell)
dev.off()


##############################################
#### 10. Neuron subclustering
##############################################

neuron.obj <- subset(xenium.obj, subset = cell_type_rctd_doublet == 'Neuron')
neuron.obj <- JoinLayers(neuron.obj)

umap_cols <- grep("^umap_\\d+$|^UMAP_\\d+$|^umaporig_\\d+$",
                  colnames(neuron.obj@meta.data), value = TRUE)
if (length(umap_cols)) neuron.obj@meta.data[umap_cols] <- NULL

neuron.obj[["pca"]]      <- NULL
neuron.obj[["harmony"]]  <- NULL
neuron.obj[["umap"]]     <- NULL
neuron.obj[["umap_orig"]] <- NULL

neuron.obj <- FindVariableFeatures(neuron.obj)
neuron.obj <- ScaleData(neuron.obj)
neuron.obj <- RunPCA(neuron.obj, npcs = 30)
neuron.obj <- RunHarmony(neuron.obj, "patient")

neuron.obj <- RunUMAP(neuron.obj, reduction = "harmony", dims = 1:30, reduction.name = 'umap')
neuron.obj <- RunUMAP(neuron.obj, reduction = "pca",     dims = 1:30, reduction.name = 'umap_orig')

neuron.obj <- FindNeighbors(neuron.obj, reduction = "harmony", dims = 1:30)
neuron.obj <- FindClusters(neuron.obj, resolution = 0.5)

neuron.marks <- FindAllMarkers(neuron.obj)
write.csv(neuron.marks, './output/Xenium/neuron.marks.csv')

DimPlot(neuron.obj, group.by = 'seurat_clusters', label = TRUE, repel = TRUE, reduction = 'umap') +
  theme(legend.position = 'none')
DimPlot(neuron.obj, group.by = 'patient', reduction = 'umap')

# ── Annotate clusters after reviewing neuron.marks.csv ───────────────────────
#
#   ┌─────────┬──────────────────────────┬──────────────────────────────────────────┐
#   │ Cluster │        Cell Type         │             Key Markers                  │
#   ├─────────┼──────────────────────────┼──────────────────────────────────────────┤
#   │ 0       │ Schwann_Myelinating      │ NTRK1, GFRA1, DIO2, SERPINA3, CYP4X1     │
#   │ 1       │ drop_CM                  │ CRYAB, HSPB3, NEB — sarcomeric/HSP        │
#   │ 2       │ Neuron_Autonomic         │ NRCAM, KCNIP4, ZFHX4, SCRN1              │
#   │ 3       │ Schwann_Nonmyelinating   │ MAL, S100B, CNTN1, NTNG1, SLC5A7         │
#   │ 4       │ drop_EC_Art              │ CLEC1A, AQP1, ICAM2, DLL4, BCL6B         │
#   │ 5       │ drop_Neutrophil          │ HVCN1, SLC22A4, SLC22A5, CPNE5           │
#   │ 6       │ drop_Macrophage          │ TLR1, CD200R1, C1QB, FCGR2B, P2RY12      │
#   │ 7       │ drop_Pericyte            │ RGS5, NOTCH3, ENPEP, BMP5, GUCY1A1       │
#   │ 8       │ drop_Pericyte2           │ KCNJ8, CCN4, KCNE4, EPOR                 │
#   │ 9       │ drop_FB                  │ NTF3, CSRP1, ITIH4, WFDC1                │
#   │ 10      │ drop_TNK                 │ SLAMF6, LEF1, KLRK1, CD6, TCF7           │
#   │ 11      │ drop_Basophil            │ HDC, HPGD, CLNK + NegControl probes      │
#   │ 12      │ drop_LowConf             │ NegControlCodeword top marker, no identity│
#   │ 13      │ drop_Unknown             │ COLEC11, PTGS1, SLC6A4, CLDN11            │
#   │ 14      │ drop_Adipocyte           │ GPD1, CIDEA, PFKFB1, PECR                 │
#   │ 15      │ drop_Proliferating       │ ASPM, KIF18B, ANLN, TPX2, KIF11           │
#   │ 16      │ drop_LEC                 │ FLT4, STAB2, CCL21, NR2F1                 │
#   └─────────┴──────────────────────────┴──────────────────────────────────────────┘

neuron_labels <- c(
  '0'  = 'Schwann_Myelinating',
  '1'  = 'drop_CM',
  '2'  = 'Neuron_Autonomic',
  '3'  = 'Schwann_Nonmyelinating',
  '4'  = 'drop_EC_Art',
  '5'  = 'drop_Neutrophil',
  '6'  = 'drop_Macrophage',
  '7'  = 'drop_Pericyte',
  '8'  = 'drop_Pericyte2',
  '9'  = 'drop_FB',
  '10' = 'drop_TNK',
  '11' = 'drop_Basophil',
  '12' = 'drop_LowConf',
  '13' = 'drop_Unknown',
  '14' = 'drop_Adipocyte',
  '15' = 'drop_Proliferating',
  '16' = 'drop_LEC'
)

neuron.obj$neuron_subtype <- unname(neuron_labels[as.character(neuron.obj$seurat_clusters)])

DimPlot(neuron.obj, group.by = 'neuron_subtype', label = TRUE, repel = TRUE, reduction = 'umap_orig') +
  theme(legend.position = 'none')

neuron.obj.clean <- subset(neuron.obj,
  subset = neuron_subtype %in% c('Schwann_Myelinating', 'Neuron_Autonomic', 'Schwann_Nonmyelinating'))

neuron.obj.clean[["pca"]]      <- NULL
neuron.obj.clean[["harmony"]]  <- NULL
neuron.obj.clean[["umap"]]     <- NULL
neuron.obj.clean[["umap_orig"]] <- NULL

umap_cols <- grep("^umap_\\d+$|^UMAP_\\d+$|^umaporig_\\d+$",
                  colnames(neuron.obj.clean@meta.data), value = TRUE)
if (length(umap_cols)) neuron.obj.clean@meta.data[umap_cols] <- NULL

neuron.obj.clean <- JoinLayers(neuron.obj.clean)
neuron.obj.clean <- FindVariableFeatures(neuron.obj.clean)
neuron.obj.clean <- ScaleData(neuron.obj.clean)
neuron.obj.clean <- RunPCA(neuron.obj.clean, npcs = 30)
neuron.obj.clean <- RunHarmony(neuron.obj.clean, "patient")
neuron.obj.clean <- RunUMAP(neuron.obj.clean, reduction = "harmony", dims = 1:30, reduction.name = 'umap')
neuron.obj.clean <- RunUMAP(neuron.obj.clean, reduction = "pca",     dims = 1:30, reduction.name = 'umap_orig')

neuron.obj.clean <- FindNeighbors(neuron.obj.clean, reduction = "harmony", dims = 1:30)
neuron.obj.clean <- FindClusters(neuron.obj.clean, resolution = 0.5)

neuron.obj.clean.mem <- neuron.obj.clean
neuron.obj.clean.mem[["Xenium"]] <- CreateAssayObject(
  counts = as(GetAssayData(neuron.obj.clean, layer = "counts"), "dgCMatrix")
)
neuron.clean.marks <- FindAllMarkers(neuron.obj.clean.mem)
rm(neuron.obj.clean.mem)
write.csv(neuron.clean.marks, './output/Xenium/neuron.clean.marks.csv')

# ── Annotate clean clusters after reviewing neuron.clean.marks.csv ──────────
#
#   ┌─────────┬──────────────────────────┬──────────────────────────────────────────────────┐
#   │ Cluster │        Cell Type         │             Key Markers / Evidence               │
#   ├─────────┼──────────────────────────┼──────────────────────────────────────────────────┤
#   │ 0       │ drop_Mesenchymal         │ ACTC1/FN1/DCN/PDGFRA mix; SOX10=35%; pct1≈pct2  │
#   │ 1       │ drop_FB                  │ ELN/VCAN/ASPN/OGN; MBP=85% but SOX10=0%         │
#   │ 2       │ drop_LowConf             │ 302 pos vs 1585 neg markers; SOX10=0%            │
#   │ 3       │ Schwann_Nonmyelinating   │ SOX10=88%, PMP22(6.3), MAL(4.5), S100B(4.5),    │
#   │         │                          │ CNTN1; MBP=0% — classic Remak Schwann            │
#   │ 4       │ drop_Ambig               │ Top markers in 5-10% cells only; SOX10=0%        │
#   │ 5       │ Neuron_Cholinergic       │ NTNG1(81%), SLC5A7(73%), RELN(100%),             │
#   │         │                          │ KCNH8(100%), NRG3, AK5, NKAIN2                   │
#   └─────────┴──────────────────────────┴──────────────────────────────────────────────────┘
#
# ── Literature support for cardiac neuronal/glial populations ────────────────
#
#   Two populations survive cleaning: one glial (Remak Schwann) and one neuronal
#   (parasympathetic ganglionic). This is expected for RV spatial transcriptomics.
#
#   Schwann_Nonmyelinating (Remak Schwann cells):
#     The heart is innervated predominantly by unmyelinated C-fibers —
#     postganglionic sympathetic and parasympathetic axons. Nonmyelinating
#     Schwann cells wrap these in Remak bundles, making them the dominant
#     glial cell type in cardiac tissue. The SOX10+/PMP22+/S100B+/MAL+/MBP-
#     profile is textbook nonmyelinating Schwann. The absence of a clean
#     myelinating Schwann population (MBP+/SOX10+) is consistent — myelinated
#     fibers are rare in the myocardium (mostly limited to large nerve trunks
#     near the epicardium).
#     Refs: Litvinukova et al. Nature 2020 (Heart Cell Atlas);
#           Kanemaru et al. Nature 2023 (spatial heart atlas);
#           Jessen & Mirsky, J Anat 2019 (Schwann cell repair roles)
#
#   Neuron_Cholinergic (intrinsic cardiac ganglionic neurons):
#     The heart contains intrinsic cardiac ganglia (ICG) — clusters of
#     parasympathetic postganglionic neurons in epicardial fat pads and
#     ventricular walls. SLC5A7 (high-affinity choline transporter/CHT1)
#     is a definitive cholinergic neuron marker required for ACh synthesis.
#     NTNG1, RELN, NRG3, NKAIN2 are neuronal adhesion/signaling molecules
#     consistent with ganglionic neurons. The RV receives vagal
#     parasympathetic innervation, so cholinergic neuron cell bodies are
#     expected.
#
#   No sympathetic neurons detected: TH/DBH absent from the marker list
#     (likely not in the 541-gene Xenium panel, or sympathetic neuron cell
#     bodies reside in stellate/paravertebral ganglia outside the heart —
#     only their axons project into the myocardium).
#
#   Clinical relevance: cardiac denervation and autonomic remodeling are
#     hallmarks of heart failure. Schwann cell loss and cholinergic
#     withdrawal in RVF could be significant in NF/pRV/RVF comparisons.

neuron_clean_labels <- c(
  '0' = 'drop_Mesenchymal',
  '1' = 'drop_FB',
  '2' = 'drop_LowConf',
  '3' = 'Schwann_Nonmyelinating',
  '4' = 'drop_Ambig',
  '5' = 'Neuron_Cholinergic'
)

neuron.obj.clean$neuron_subtype <- unname(neuron_clean_labels[as.character(neuron.obj.clean$seurat_clusters)])

DimPlot(neuron.obj.clean, group.by = 'neuron_subtype', label = TRUE, repel = TRUE, reduction = 'umap_orig') +
  theme(legend.position = 'none')

neuron.obj.clean.clean <- subset(neuron.obj.clean,
  subset = neuron_subtype %in% c('Schwann_Nonmyelinating', 'Neuron_Cholinergic'))

DimPlot(neuron.obj.clean.clean, group.by = 'neuron_subtype', label = TRUE, repel = TRUE, reduction = 'umap_orig') +
  theme(legend.position = 'none')
DimPlot(neuron.obj.clean.clean, group.by = 'patient', reduction = 'umap_orig')

saveRDS(neuron.obj.clean.clean, './output/Xenium/neuron_clean.rds')

pdf('./output/Xenium/Xenium_neuron_snUMAP.pdf', width = 5, height = 5)
PlotEmbedding(neuron.obj.clean.clean, group.by = 'neuron_subtype', reduction = 'umap_orig',
              point_size = .1, plot_under = TRUE, plot_theme = umap_theme() + NoLegend(),
              raster_dpi = 400, raster_scale = 0.5)
dev.off()

pdf('./output/Xenium/Xenium_neuron_snUMAP_pt.pdf', width = 5, height = 5)
PlotEmbedding(neuron.obj.clean.clean, group.by = 'patient', reduction = 'umap_orig',
              point_size = .1, plot_under = TRUE, plot_theme = umap_theme() + NoLegend(),
              raster_dpi = 400, raster_scale = 0.5)
dev.off()

##############################################
#### Neuron subtype proportions by disease state
##############################################

neuron_meta <- neuron.obj.clean.clean@meta.data %>%
  tibble::rownames_to_column('cell') %>%
  select(cell, patient, group, neuron_subtype)

neuron_totals <- xenium.obj@meta.data %>%
  tibble::rownames_to_column('cell') %>%
  select(cell, patient, group) %>%
  group_by(patient, group) %>%
  summarise(total = n(), .groups = 'drop')

neuron_prop <- neuron_meta %>%
  group_by(patient, group, neuron_subtype) %>%
  summarise(n = n(), .groups = 'drop') %>%
  left_join(neuron_totals, by = c('patient', 'group')) %>%
  mutate(proportion = n / total) %>%
  complete(nesting(patient, group), neuron_subtype, fill = list(n = 0, proportion = 0)) %>%
  mutate(group = factor(group, levels = c('NF', 'pRV', 'RVF')))

p_neuron_prop <- ggplot(neuron_prop, aes(x = group, y = proportion, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.5) +
  stat_compare_means(comparisons = pairwise_comparisons, method = 'wilcox.test', label = 'p.signif') +
  stat_compare_means(method = 'kruskal.test', label.y.npc = 0.95, size = 3) +
  facet_wrap(~ neuron_subtype, scales = 'free_y') +
  scale_x_discrete(limits = c('NF', 'pRV', 'RVF')) +
  theme_bw() +
  theme(legend.position = 'none',
        strip.text = element_text(face = 'bold'),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = 'Proportion of all cells', title = 'Neuron subtype distribution by group')

pdf('./output/Xenium/neuron_proportions.pdf', width = 8, height = 6)
print(p_neuron_prop)
dev.off()

##############################################
#### Neuron fractional area coverage by disease state
##############################################

cell_areas <- read.csv('./output/Xenium/cell_areas_2d.csv', row.names = 1)

total_area <- xenium.obj@meta.data %>%
  tibble::rownames_to_column('cell') %>%
  select(cell, patient, group) %>%
  left_join(cell_areas, by = 'cell') %>%
  group_by(patient, group) %>%
  summarise(total_area = sum(cell_area_2d, na.rm = TRUE), .groups = 'drop')

neuron_area_meta <- neuron.obj.clean.clean@meta.data %>%
  tibble::rownames_to_column('cell') %>%
  select(cell, patient, group, neuron_subtype) %>%
  left_join(cell_areas, by = 'cell')

neuron_area_prop <- neuron_area_meta %>%
  group_by(patient, group, neuron_subtype) %>%
  summarise(neuron_area = sum(cell_area_2d, na.rm = TRUE), .groups = 'drop') %>%
  left_join(total_area, by = c('patient', 'group')) %>%
  mutate(frac_area = neuron_area / total_area) %>%
  complete(nesting(patient, group), neuron_subtype, fill = list(neuron_area = 0, frac_area = 0)) %>%
  mutate(group = factor(group, levels = c('NF', 'pRV', 'RVF')))

p_neuron_area <- ggplot(neuron_area_prop, aes(x = group, y = frac_area, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.5) +
  stat_compare_means(comparisons = pairwise_comparisons, method = 'wilcox.test', label = 'p.signif') +
  stat_compare_means(method = 'kruskal.test', label.y.npc = 0.95, size = 3) +
  facet_wrap(~ neuron_subtype, scales = 'free_y') +
  scale_x_discrete(limits = c('NF', 'pRV', 'RVF')) +
  theme_bw() +
  theme(legend.position = 'none',
        strip.text = element_text(face = 'bold'),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = 'Fractional area coverage', title = 'Neuron subtype fractional area by group')

pdf('./output/Xenium/neuron_area_coverage.pdf', width = 8, height = 6)
print(p_neuron_area)
dev.off()

##############################################
#### Pseudobulk DESeq2 by Neuron subtype
##############################################

neuron_pseudo_meta <- unique(neuron.obj.clean.clean@meta.data[, c('patient', 'group', 'neuron_subtype')])
neuron_pseudo_meta$patient       <- as.character(neuron_pseudo_meta$patient)
neuron_pseudo_meta$neuron_subtype <- as.character(neuron_pseudo_meta$neuron_subtype)
neuron_pseudo_meta$pseudobulk_id <- paste0(neuron_pseudo_meta$patient, '..', neuron_pseudo_meta$neuron_subtype)
neuron_pseudo_meta$group <- factor(neuron_pseudo_meta$group, levels = c('NF', 'pRV', 'RVF'))
rownames(neuron_pseudo_meta) <- neuron_pseudo_meta$pseudobulk_id

neuron.obj.clean.clean$pseudobulk_id <- paste0(
  as.character(neuron.obj.clean.clean$patient), '..',
  as.character(neuron.obj.clean.clean$neuron_subtype)
)

neuron_pseudo_counts <- AggregateExpression(
  neuron.obj.clean.clean,
  assays = 'Xenium',
  group.by = 'pseudobulk_id',
  slot = 'counts',
  return.seurat = FALSE
)$Xenium

colnames(neuron_pseudo_counts) <- gsub('-', '_', sub('^g', '', colnames(neuron_pseudo_counts)))

neuron_cell_types <- unique(neuron_pseudo_meta$neuron_subtype)
neuron_deseq_results <- list()

for (ct in neuron_cell_types) {
  ct_samples <- neuron_pseudo_meta$pseudobulk_id[neuron_pseudo_meta$neuron_subtype == ct]
  ct_samples <- ct_samples[ct_samples %in% colnames(neuron_pseudo_counts)]
  if (length(ct_samples) < 2) next

  ct_counts <- neuron_pseudo_counts[, ct_samples, drop = FALSE]
  ct_meta   <- neuron_pseudo_meta[ct_samples, , drop = FALSE]

  if (length(unique(ct_meta$group)) < 2) next

  keep <- rowSums(ct_counts >= 5) >= 2
  ct_counts <- ct_counts[keep, , drop = FALSE]
  if (nrow(ct_counts) < 10) next

  dds <- DESeqDataSetFromMatrix(countData = ct_counts, colData = ct_meta, design = ~ group)
  dds <- tryCatch(DESeq(dds), error = function(e) NULL)
  if (is.null(dds)) next

  for (comp in list(c('RVF', 'NF'), c('pRV', 'NF'), c('RVF', 'pRV'))) {
    comp_label <- paste0(comp[1], '_vs_', comp[2])
    res <- tryCatch(
      results(dds, contrast = c('group', comp[1], comp[2])) %>%
        as.data.frame() %>%
        tibble::rownames_to_column('gene') %>%
        mutate(cell_type = ct, comparison = comp_label) %>%
        arrange(padj),
      error = function(e) NULL
    )
    if (!is.null(res)) neuron_deseq_results[[paste0(ct, '_', comp_label)]] <- res
  }
}

neuron_all_deseq <- bind_rows(neuron_deseq_results)
write.csv(neuron_all_deseq, './output/Xenium/neuron_pseudobulk_deseq2.csv', row.names = FALSE)

neuron_sig_summary <- neuron_all_deseq %>%
  filter(padj < 0.05, abs(log2FoldChange) > 0.5) %>%
  group_by(cell_type, comparison) %>%
  summarise(n_up = sum(log2FoldChange > 0), n_down = sum(log2FoldChange < 0), .groups = 'drop')
print(neuron_sig_summary)

##############################################
#### Neuron subtypes: Volcano + GO/Reactome
##############################################

neuron_enrich_combined <- list()

for (ct in unique(neuron_all_deseq$cell_type)) {
  ct_slug <- gsub('_', '-', ct)
  cat('Neuron processing:', ct, '\n')

  dat_rvf     <- neuron_all_deseq %>% filter(cell_type == ct, comparison == 'RVF_vs_NF')  %>% cap_pval()
  dat_prv     <- neuron_all_deseq %>% filter(cell_type == ct, comparison == 'pRV_vs_NF')  %>% cap_pval()
  dat_rvf_prv <- neuron_all_deseq %>% filter(cell_type == ct, comparison == 'RVF_vs_pRV') %>% cap_pval()

  make_volcano_pdf(dat_rvf,     paste0('./output/Xenium/neuron_', ct_slug, '_volcano_RVF_vs_NF.pdf'),  paste(ct, ': RVF vs NF'))
  make_volcano_pdf(dat_prv,     paste0('./output/Xenium/neuron_', ct_slug, '_volcano_pRV_vs_NF.pdf'),  paste(ct, ': pRV vs NF'))
  make_volcano_pdf(dat_rvf_prv, paste0('./output/Xenium/neuron_', ct_slug, '_volcano_RVF_vs_pRV.pdf'), paste(ct, ': RVF vs pRV'))

  run_enrich <- function(dat, comp) {
    sig  <- dat %>% filter(padj < 0.05, abs(log2FoldChange) > 0.5)
    up   <- sig %>% filter(log2FoldChange > 0) %>% pull(gene)
    down <- sig %>% filter(log2FoldChange < 0) %>% pull(gene)
    list(
      up   = if (length(up)   > 0) enrichr(up,   enrichR_dbs) else NULL,
      down = if (length(down) > 0) enrichr(down, enrichR_dbs) else NULL
    )
  }

  er_rvf     <- run_enrich(dat_rvf,     'RVF_vs_NF')
  er_prv     <- run_enrich(dat_prv,     'pRV_vs_NF')
  er_rvf_prv <- run_enrich(dat_rvf_prv, 'RVF_vs_pRV')

  ct_enrich <- bind_rows(
    combine_enrich(er_rvf$up,       'up',   'RVF_vs_NF'),
    combine_enrich(er_rvf$down,     'down', 'RVF_vs_NF'),
    combine_enrich(er_prv$up,       'up',   'pRV_vs_NF'),
    combine_enrich(er_prv$down,     'down', 'pRV_vs_NF'),
    combine_enrich(er_rvf_prv$up,   'up',   'RVF_vs_pRV'),
    combine_enrich(er_rvf_prv$down, 'down', 'RVF_vs_pRV')
  )
  if (!is.null(ct_enrich) && nrow(ct_enrich) > 0) {
    ct_enrich$cell_type <- ct
    neuron_enrich_combined[[ct]] <- ct_enrich
  }

  pdf(paste0('./output/Xenium/neuron_', ct_slug, '_enrichr_plots.pdf'), width = 12, height = 7)
  for (comp_label in c('RVF_vs_NF', 'pRV_vs_NF', 'RVF_vs_pRV')) {
    er <- switch(comp_label, RVF_vs_NF = er_rvf, pRV_vs_NF = er_prv, RVF_vs_pRV = er_rvf_prv)
    for (dir_label in c('up', 'down')) {
      er_dir <- er[[dir_label]]
      if (!is.null(er_dir)) {
        p1 <- plot_enrich(er_dir$GO_Biological_Process_2023, paste(ct, comp_label, toupper(dir_label), 'GO BP'))
        p2 <- plot_enrich(er_dir$Reactome_2022,              paste(ct, comp_label, toupper(dir_label), 'Reactome'))
        if (!is.null(p1)) print(p1)
        if (!is.null(p2)) print(p2)
      }
    }
  }
  dev.off()
}

neuron_all_enrich <- bind_rows(neuron_enrich_combined)
write.csv(neuron_all_enrich, './output/Xenium/neuron_all_enrichr.csv', row.names = FALSE)


##############################################
#### 11. Mural cell (SM + PC) subclustering
##############################################

mural.obj <- subset(xenium.obj, subset = cell_type_rctd_doublet %in% c('SM', 'PC'))
mural.obj <- JoinLayers(mural.obj)

umap_cols <- grep("^umap_\\d+$|^UMAP_\\d+$|^umaporig_\\d+$",
                  colnames(mural.obj@meta.data), value = TRUE)
if (length(umap_cols)) mural.obj@meta.data[umap_cols] <- NULL

mural.obj[["pca"]]      <- NULL
mural.obj[["harmony"]]  <- NULL
mural.obj[["umap"]]     <- NULL
mural.obj[["umap_orig"]] <- NULL

mural.obj <- FindVariableFeatures(mural.obj)
mural.obj <- ScaleData(mural.obj)
mural.obj <- RunPCA(mural.obj, npcs = 30)
mural.obj <- RunHarmony(mural.obj, "patient")

mural.obj <- RunUMAP(mural.obj, reduction = "harmony", dims = 1:30, reduction.name = 'umap')
mural.obj <- RunUMAP(mural.obj, reduction = "pca",     dims = 1:30, reduction.name = 'umap_orig')

mural.obj <- FindNeighbors(mural.obj, reduction = "harmony", dims = 1:30)
mural.obj <- FindClusters(mural.obj, resolution = 0.5)

mural.obj.mem <- mural.obj
mural.obj.mem[["Xenium"]] <- CreateAssayObject(
  counts = as(GetAssayData(mural.obj, layer = "counts"), "dgCMatrix")
)
mural.marks <- FindAllMarkers(mural.obj.mem)
rm(mural.obj.mem)
write.csv(mural.marks, './output/Xenium/mural.marks.csv')

DimPlot(mural.obj, group.by = 'seurat_clusters', label = TRUE, repel = TRUE, reduction = 'umap') +
  theme(legend.position = 'none')
DimPlot(mural.obj, group.by = 'patient', reduction = 'umap')
# Also visualise original annotation to see SM vs PC distribution
DimPlot(mural.obj, group.by = 'cell_type_rctd_doublet', reduction = 'umap')

# ── Annotate clusters after reviewing mural.marks.csv ─────────────────────────
#
#   ┌─────────┬──────────────────────┬──────────────────────────────────────────────────────┐
#   │ Cluster │      Cell Type       │            Key Markers / Evidence                    │
#   ├─────────┼──────────────────────┼──────────────────────────────────────────────────────┤
#   │ 0       │ Pericyte_KCNJ8       │ RGS5+2.8, KCNJ8+7.2, NOTCH3+1.6, ABCC9+2.9,       │
#   │         │                      │ ENPEP+1.9, FABP4+3.3; ACTA2/MYH11 depleted          │
#   │ 1       │ drop_FB              │ COL1A1+3.2, COL1A2+1.8; ACTG2+49, FBN1+37,          │
#   │         │                      │ FBLN5+34; all mural markers depleted                 │
#   │ 2       │ VSMC_Arterial        │ MYH11+174, ACTA2+75, CNN1+16, MYLK+19, TAGLN+3.5,   │
#   │         │                      │ HEY2+6.9, CARMN+2.3, TNC+47 — arterial SM            │
#   │ 3       │ drop_CM              │ NPPA+703, MYL2+305, NPPB+292, MYBPC3+87, TNNI3+47   │
#   │ 4       │ drop_EC              │ VWF+27, AQP1+12.7, DLL4+12.4, BCL6B+12              │
#   │ 5       │ drop_Neutrophil      │ HVCN1+10.3, SLC22A4+4.9, SLC22A5+5.4                │
#   │ 6       │ Pericyte              │ RGS5(0.99), NOTCH3(1.00), ENPEP(0.96), CARMN(1.00); │
#   │         │                      │ THBS4+11 top marker — generic pericyte                │
#   │ 7       │ VSMC                  │ MYH11(0.81), ACTA2(0.77), TAGLN(0.73), CARMN(1.00), │
#   │         │                      │ MYOCD(0.98), NOTCH3(1.00); CXCL9+14.6                │
#   │ 8       │ drop_Macrophage      │ LYVE1+23, ACE2+20, BLNK+14, CD200R1+14, C1QB+14     │
#   │ 9       │ Pericyte_Activated   │ STEAP4+34, NR4A1+24, NR2F2+12, SFRP2+12, THY1+9;    │
#   │         │                      │ RGS5+0.5; balanced 2204 up / 2087 down — stress/      │
#   │         │                      │ activation signature in venous mural cells             │
#   │ 10      │ Mural_Unspecified     │ MYOCD(0.83), CARMN(0.79); no distinctive subtype      │
#   │         │                      │ markers but real mural identity                        │
#   │ 11      │ Pericyte2             │ CARMN(1.00), RGS5(0.99); real pericyte, no unique     │
#   │         │                      │ subtype markers                                        │
#   │ 12      │ drop_TNK             │ IKZF3+14, IL7R+17, LCK+15, RASGRP1+16, CD3E+2.9     │
#   │ 13      │ drop_Unknown         │ COLEC11+14, SLC6A4+13, CLDN11+12, PTGS1+14;          │
#   │         │                      │ heterogeneous, no coherent mural identity              │
#   │ 14      │ drop_Adipocyte       │ ADIPOQ+82, PLIN4+55, CIDEA+16, GPD1+15, FABP4+3.6   │
#   │ 15      │ drop_Neuron          │ SOX10+14, SLC5A7+13, NTNG1+9, ERBB3+15               │
#   │ 16      │ drop_LEC             │ CCL21+15, CDH5+2.9 — lymphatic endothelial            │
#   │ 17      │ drop_Proliferating   │ MKI67+22, TPX2+15, TOP2A+9, CENPE/KIF11/NCAPG+16    │
#   └─────────┴──────────────────────┴──────────────────────────────────────────────────────┘

mural_labels <- c(
  '0'  = 'Pericyte_KCNJ8',
  '1'  = 'drop_FB',
  '2'  = 'VSMC_Arterial',
  '3'  = 'drop_CM',
  '4'  = 'drop_EC',
  '5'  = 'drop_Neutrophil',
  '6'  = 'Pericyte',
  '7'  = 'VSMC',
  '8'  = 'drop_Macrophage',
  '9'  = 'Pericyte_Activated',
  '10' = 'Mural_Unspecified',
  '11' = 'Pericyte2',
  '12' = 'drop_TNK',
  '13' = 'drop_Unknown',
  '14' = 'drop_Adipocyte',
  '15' = 'drop_Neuron',
  '16' = 'drop_LEC',
  '17' = 'drop_Proliferating'
)

mural.obj$mural_subtype <- unname(mural_labels[as.character(mural.obj$seurat_clusters)])

DimPlot(mural.obj, group.by = 'mural_subtype', label = TRUE, repel = TRUE, reduction = 'umap') +
  theme(legend.position = 'none')

# ── Clean subset: keep only true mural subtypes ──
keep_cells <- !is.na(mural.obj$mural_subtype) & !grepl('^drop_', mural.obj$mural_subtype)
mural.obj.clean <- subset(mural.obj, cells = colnames(mural.obj)[keep_cells])

mural.obj.clean[["pca"]]      <- NULL
mural.obj.clean[["harmony"]]  <- NULL
mural.obj.clean[["umap"]]     <- NULL
mural.obj.clean[["umap_orig"]] <- NULL

umap_cols <- grep("^umap_\\d+$|^UMAP_\\d+$|^umaporig_\\d+$",
                  colnames(mural.obj.clean@meta.data), value = TRUE)
if (length(umap_cols)) mural.obj.clean@meta.data[umap_cols] <- NULL

mural.obj.clean <- JoinLayers(mural.obj.clean)
mural.obj.clean <- FindVariableFeatures(mural.obj.clean)
mural.obj.clean <- ScaleData(mural.obj.clean)
mural.obj.clean <- RunPCA(mural.obj.clean, npcs = 30)
mural.obj.clean <- RunHarmony(mural.obj.clean, "patient")
mural.obj.clean <- RunUMAP(mural.obj.clean, reduction = "harmony", dims = 1:30, reduction.name = 'umap')
mural.obj.clean <- RunUMAP(mural.obj.clean, reduction = "pca",     dims = 1:30, reduction.name = 'umap_orig')

mural.obj.clean <- FindNeighbors(mural.obj.clean, reduction = "harmony", dims = 1:30)
mural.obj.clean <- FindClusters(mural.obj.clean, resolution = 0.5)

mural.obj.clean.mem <- mural.obj.clean
mural.obj.clean.mem[["Xenium"]] <- CreateAssayObject(
  counts = as(GetAssayData(mural.obj.clean, layer = "counts"), "dgCMatrix")
)
mural.clean.marks <- FindAllMarkers(mural.obj.clean.mem)
rm(mural.obj.clean.mem)
write.csv(mural.clean.marks, './output/Xenium/mural.clean.marks.csv')

# ── Annotate mural clean clusters after reviewing mural.clean.marks.csv ──────
#   All 12 clusters are real mural subtypes — no drops needed.
#
#   Literature support:
#   - Litvinukova et al. (2020) Nature 588:466 — human heart cell atlas; defined
#     KCNJ8/ABCC9 pericytes and arterial/basic SMC subtypes
#   - Koenig et al. (2022) Nat Commun 13:3027 — cardiac vascular niche in HF;
#     3 VSMC + 4 pericyte subclusters incl. stressed & IFN-responsive mural cells
#   - Muhl et al. (2022) Front Cardiovasc Med 9:876591 — organ-specific pericyte
#     markers; cardiac pericytes defined by KCNJ8, ABCC9, RGS5
#   - Fischer et al. (2004) Genes Dev 18:901 — Hey1/Hey2 (Notch targets) required
#     for arterial vascular identity and VSMC specification
#   - Doi et al. (2006) J Biol Chem 281:28555 — HEY2 antagonises MYOCD/SRF at
#     CArG boxes, modulating contractile gene expression in arterial VSMC
#   - Frangogiannis et al. (2012) FASEB J 26:2363 — THBS4 regulates cardiac
#     fibrosis/remodeling under pressure overload; expressed by ventricular pericytes
#   - Dobnikar et al. (2018) Nat Commun 9:4567 — disease-relevant transcriptional
#     signatures in individual SMC; contractile-to-synthetic switching via KLF4
#   - Liu & Gomez (2019) JAHA 8:e031121 — VSMC phenotypic switching: contractile
#     (MYH11/ACTA2/TAGLN/MYOCD) → synthetic (ECM genes, ELN, PRG4) spectrum
#   - van Kuijk et al. (2022) Cell Commun Signal 20:180 — CXCL9/10/11 induced by
#     IFNg in VSMC; inflammatory VSMC phenotype in vascular disease
#
#   ┌─────────┬──────────────────────┬──────────────────────────────────────────────────────┐
#   │ Cluster │      Cell Type       │            Key Markers / Evidence                    │
#   ├─────────┼──────────────────────┼──────────────────────────────────────────────────────┤
#   │ 0       │ Pericyte_KCNJ8       │ KCNJ8+6.0, RGS5+1.7, ABCC9+2.2, FABP4+3.5,        │
#   │         │                      │ ENPEP+0.71 — classical KCNJ8+ pericyte               │
#   │         │                      │ [Litvinukova 2020; Muhl 2022]                        │
#   │ 1       │ Pericyte              │ RGS5+0.42, KCNJ8+2.2, ABCC9+1.2, COL1A2+0.51;     │
#   │         │                      │ weaker KCNJ8 than c0, generic pericyte               │
#   │         │                      │ [Litvinukova 2020 — "basic pericyte"]                │
#   │ 2       │ VSMC_Arterial        │ MYH11+286, ACTA2+111, CNN1+28, MYLK+74, TNC+70,    │
#   │         │                      │ CARMN+1.2 — strong contractile arterial VSMC          │
#   │         │                      │ [Litvinukova 2020 — "arterial SMC"; Koenig 2022]    │
#   │ 3       │ Pericyte_THBS4       │ RGS5+1.2, NOTCH3+0.51, ABCC9+0.90; THBS4+37.5      │
#   │         │                      │ top marker — thrombospondin-4 enriched pericyte       │
#   │         │                      │ [Frangogiannis 2012 — THBS4 in ventricular PCs]     │
#   │ 4       │ VSMC_Synthetic       │ MYOCD+0.85, TAGLN+1.6, TNC+67, PRG4+17, MKX+10,    │
#   │         │                      │ ELN+8.8 — dedifferentiated/ECM-remodeling VSMC        │
#   │         │                      │ [Liu 2019; Dobnikar 2018 — contractile→synthetic]    │
#   │ 5       │ VSMC                  │ MYLK+28, MYOCD+0.80, TAGLN+1.3, CARMN+1.2;         │
#   │         │                      │ generic SM with MYLK as top distinguishing marker     │
#   │         │                      │ [Koenig 2022 — "VSMC1/VSMC2" differing by Acta2]    │
#   │ 6       │ Pericyte_Mixed       │ RGS5+0.88, NOTCH3+0.42, THBS4+37.5; NPPA+42,       │
#   │         │                      │ LYVE1+19 — pericyte with ambient CM signal            │
#   │         │                      │ [Koenig 2022 — stromal pericyte with ambient noise]  │
#   │ 7       │ VSMC_Inflamed        │ MYOCD+0.82, TAGLN+1.2, CARMN+1.15; CXCL9+14        │
#   │         │                      │ — IFNg-responsive / inflamed VSMC                     │
#   │         │                      │ [Koenig 2022 "IntPeri"; van Kuijk 2022 — CXCL9/IFNg]│
#   │ 8       │ VSMC_Arterial2       │ MYH11+131, ACTA2+90, FN1+91, HEY2 pct1=0.50,       │
#   │         │                      │ MYOCD+1.5, TAGLN+1.6, CARMN+0.77 — HEY2+ arterial   │
#   │         │                      │ [Fischer 2004; Doi 2006 — HEY2/Notch arterial ID]   │
#   │ 9       │ Pericyte_Quiescent   │ RGS5 pct1=0.97, ABCC9 pct1=0.98, CARMN pct1=1.00;  │
#   │         │                      │ top: DNMT1, ARHGAP22, SLCO2B1 — sparse quiescent PC   │
#   │         │                      │ [Muhl 2022 — quiescent cardiac PC signature]         │
#   │ 10      │ Mural_Unspecified    │ MYOCD+0.75, CARMN pct1=0.75, NOTCH3 pct1=0.28;      │
#   │         │                      │ NPPB+28 top marker — ambiguous with CM-adjacent noise  │
#   │         │                      │ [intermediate identity; possible ambient NPPB from CM]│
#   │ 11      │ VSMC_Transitional    │ MYOCD+0.97, TAGLN+1.49, MYH11 pct1=0.96;            │
#   │         │                      │ top: GRIK1, CADM2 — SM lineage, partially dediffer.   │
#   │         │                      │ [Liu 2019 — contractile→synthetic transition state]  │
#   └─────────┴──────────────────────┴──────────────────────────────────────────────────────┘

mural_clean_labels <- c(
  '0'  = 'Pericyte_KCNJ8',
  '1'  = 'Pericyte',
  '2'  = 'VSMC_Arterial',
  '3'  = 'Pericyte_THBS4',
  '4'  = 'VSMC_Synthetic',
  '5'  = 'VSMC',
  '6'  = 'Pericyte_Mixed',
  '7'  = 'VSMC_Inflamed',
  '8'  = 'VSMC_Arterial2',
  '9'  = 'Pericyte_Quiescent',
  '10' = 'Mural_Unspecified',
  '11' = 'VSMC_Transitional'
)

mural.obj.clean$mural_subtype <- unname(mural_clean_labels[as.character(mural.obj.clean$seurat_clusters)])

DimPlot(mural.obj.clean, group.by = 'mural_subtype', label = TRUE, repel = TRUE, reduction = 'umap') +
  theme(legend.position = 'none')
DimPlot(mural.obj.clean, group.by = 'patient', reduction = 'umap')

saveRDS(mural.obj.clean, './output/Xenium/mural_clean.rds')

pdf('./output/Xenium/Xenium_mural_snUMAP.pdf', width = 5, height = 5)
PlotEmbedding(mural.obj.clean, group.by = 'mural_subtype', reduction = 'umap',
              point_size = .05, plot_under = TRUE, plot_theme = umap_theme() + NoLegend(),
              raster_dpi = 400, raster_scale = 0.5)
dev.off()

pdf('./output/Xenium/Xenium_mural_snUMAP_pt.pdf', width = 5, height = 5)
PlotEmbedding(mural.obj.clean, group.by = 'patient', reduction = 'umap',
              point_size = .05, plot_under = TRUE, plot_theme = umap_theme() + NoLegend(),
              raster_dpi = 400, raster_scale = 0.5)
dev.off()

##############################################
#### Mural subtype proportions by disease state
##############################################

mural_meta <- mural.obj.clean@meta.data %>%
  tibble::rownames_to_column('cell') %>%
  select(cell, patient, group, mural_subtype)

mural_totals <- xenium.obj@meta.data %>%
  tibble::rownames_to_column('cell') %>%
  select(cell, patient, group) %>%
  group_by(patient, group) %>%
  summarise(total = n(), .groups = 'drop')

mural_prop <- mural_meta %>%
  group_by(patient, group, mural_subtype) %>%
  summarise(n = n(), .groups = 'drop') %>%
  left_join(mural_totals, by = c('patient', 'group')) %>%
  mutate(proportion = n / total) %>%
  complete(nesting(patient, group), mural_subtype, fill = list(n = 0, proportion = 0)) %>%
  mutate(group = factor(group, levels = c('NF', 'pRV', 'RVF')))

p_mural_prop <- ggplot(mural_prop, aes(x = group, y = proportion, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.5) +
  stat_compare_means(comparisons = pairwise_comparisons, method = 'wilcox.test', label = 'p.signif') +
  stat_compare_means(method = 'kruskal.test', label.y.npc = 0.95, size = 3) +
  facet_wrap(~ mural_subtype, scales = 'free_y') +
  scale_x_discrete(limits = c('NF', 'pRV', 'RVF')) +
  theme_bw() +
  theme(legend.position = 'none',
        strip.text = element_text(face = 'bold'),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = 'Proportion of all cells', title = 'Mural cell subtype distribution by group')

pdf('./output/Xenium/mural_proportions.pdf', width = 10, height = 6)
print(p_mural_prop)
dev.off()

##############################################
#### Mural fractional area coverage by disease state
##############################################

cell_areas <- read.csv('./output/Xenium/cell_areas_2d.csv', row.names = 1)

total_area <- xenium.obj@meta.data %>%
  tibble::rownames_to_column('cell') %>%
  select(cell, patient, group) %>%
  left_join(cell_areas, by = 'cell') %>%
  group_by(patient, group) %>%
  summarise(total_area = sum(cell_area_2d, na.rm = TRUE), .groups = 'drop')

mural_area_meta <- mural.obj.clean@meta.data %>%
  tibble::rownames_to_column('cell') %>%
  select(cell, patient, group, mural_subtype) %>%
  left_join(cell_areas, by = 'cell')

mural_area_prop <- mural_area_meta %>%
  group_by(patient, group, mural_subtype) %>%
  summarise(mural_area = sum(cell_area_2d, na.rm = TRUE), .groups = 'drop') %>%
  left_join(total_area, by = c('patient', 'group')) %>%
  mutate(frac_area = mural_area / total_area) %>%
  complete(nesting(patient, group), mural_subtype, fill = list(mural_area = 0, frac_area = 0)) %>%
  mutate(group = factor(group, levels = c('NF', 'pRV', 'RVF')))

p_mural_area <- ggplot(mural_area_prop, aes(x = group, y = frac_area, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.5) +
  stat_compare_means(comparisons = pairwise_comparisons, method = 'wilcox.test', label = 'p.signif') +
  stat_compare_means(method = 'kruskal.test', label.y.npc = 0.95, size = 3) +
  facet_wrap(~ mural_subtype, scales = 'free_y') +
  scale_x_discrete(limits = c('NF', 'pRV', 'RVF')) +
  theme_bw() +
  theme(legend.position = 'none',
        strip.text = element_text(face = 'bold'),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = 'Fractional area coverage', title = 'Mural cell subtype fractional area by group')

pdf('./output/Xenium/mural_area_coverage.pdf', width = 10, height = 6)
print(p_mural_area)
dev.off()

##############################################
#### Pseudobulk DESeq2 by mural subtype
##############################################

mural_pseudo_meta <- unique(mural.obj.clean@meta.data[, c('patient', 'group', 'mural_subtype')])
mural_pseudo_meta$patient       <- as.character(mural_pseudo_meta$patient)
mural_pseudo_meta$mural_subtype <- as.character(mural_pseudo_meta$mural_subtype)
mural_pseudo_meta$pseudobulk_id <- paste0(mural_pseudo_meta$patient, '..', mural_pseudo_meta$mural_subtype)
mural_pseudo_meta$group <- factor(mural_pseudo_meta$group, levels = c('NF', 'pRV', 'RVF'))
rownames(mural_pseudo_meta) <- mural_pseudo_meta$pseudobulk_id

mural.obj.clean$pseudobulk_id <- paste0(
  as.character(mural.obj.clean$patient), '..',
  as.character(mural.obj.clean$mural_subtype)
)

mural_pseudo_counts <- AggregateExpression(
  mural.obj.clean,
  assays = 'Xenium',
  group.by = 'pseudobulk_id',
  slot = 'counts',
  return.seurat = FALSE
)$Xenium

colnames(mural_pseudo_counts) <- gsub('-', '_', sub('^g', '', colnames(mural_pseudo_counts)))

mural_cell_types <- unique(mural_pseudo_meta$mural_subtype)
mural_deseq_results <- list()

for (ct in mural_cell_types) {
  ct_samples <- mural_pseudo_meta$pseudobulk_id[mural_pseudo_meta$mural_subtype == ct]
  ct_samples <- ct_samples[ct_samples %in% colnames(mural_pseudo_counts)]
  if (length(ct_samples) < 2) next

  ct_counts <- mural_pseudo_counts[, ct_samples, drop = FALSE]
  ct_meta   <- mural_pseudo_meta[ct_samples, , drop = FALSE]

  if (length(unique(ct_meta$group)) < 2) next

  keep <- rowSums(ct_counts >= 5) >= 2
  ct_counts <- ct_counts[keep, , drop = FALSE]
  if (nrow(ct_counts) < 10) next

  dds <- DESeqDataSetFromMatrix(countData = ct_counts, colData = ct_meta, design = ~ group)
  dds <- tryCatch(DESeq(dds), error = function(e) NULL)
  if (is.null(dds)) next

  for (comp in list(c('RVF', 'NF'), c('pRV', 'NF'), c('RVF', 'pRV'))) {
    comp_label <- paste0(comp[1], '_vs_', comp[2])
    res <- tryCatch(
      results(dds, contrast = c('group', comp[1], comp[2])) %>%
        as.data.frame() %>%
        tibble::rownames_to_column('gene') %>%
        mutate(cell_type = ct, comparison = comp_label) %>%
        arrange(padj),
      error = function(e) NULL
    )
    if (!is.null(res)) mural_deseq_results[[paste0(ct, '_', comp_label)]] <- res
  }
}

mural_all_deseq <- bind_rows(mural_deseq_results)
write.csv(mural_all_deseq, './output/Xenium/mural_pseudobulk_deseq2.csv', row.names = FALSE)

mural_sig_summary <- mural_all_deseq %>%
  filter(padj < 0.05, abs(log2FoldChange) > 0.5) %>%
  group_by(cell_type, comparison) %>%
  summarise(n_up = sum(log2FoldChange > 0), n_down = sum(log2FoldChange < 0), .groups = 'drop')
print(mural_sig_summary)

##############################################
#### Mural subtypes: Volcano + GO/Reactome
##############################################

mural_enrich_combined <- list()

for (ct in unique(mural_all_deseq$cell_type)) {
  ct_slug <- gsub('_', '-', ct)
  cat('Mural processing:', ct, '\n')

  dat_rvf     <- mural_all_deseq %>% filter(cell_type == ct, comparison == 'RVF_vs_NF')  %>% cap_pval()
  dat_prv     <- mural_all_deseq %>% filter(cell_type == ct, comparison == 'pRV_vs_NF')  %>% cap_pval()
  dat_rvf_prv <- mural_all_deseq %>% filter(cell_type == ct, comparison == 'RVF_vs_pRV') %>% cap_pval()

  make_volcano_pdf(dat_rvf,     paste0('./output/Xenium/mural_', ct_slug, '_volcano_RVF_vs_NF.pdf'),  paste(ct, ': RVF vs NF'))
  make_volcano_pdf(dat_prv,     paste0('./output/Xenium/mural_', ct_slug, '_volcano_pRV_vs_NF.pdf'),  paste(ct, ': pRV vs NF'))
  make_volcano_pdf(dat_rvf_prv, paste0('./output/Xenium/mural_', ct_slug, '_volcano_RVF_vs_pRV.pdf'), paste(ct, ': RVF vs pRV'))

  run_enrich <- function(dat, comp) {
    sig  <- dat %>% filter(padj < 0.05, abs(log2FoldChange) > 0.5)
    up   <- sig %>% filter(log2FoldChange > 0) %>% pull(gene)
    down <- sig %>% filter(log2FoldChange < 0) %>% pull(gene)
    list(
      up   = if (length(up)   > 0) enrichr(up,   enrichR_dbs) else NULL,
      down = if (length(down) > 0) enrichr(down, enrichR_dbs) else NULL
    )
  }

  er_rvf     <- run_enrich(dat_rvf,     'RVF_vs_NF')
  er_prv     <- run_enrich(dat_prv,     'pRV_vs_NF')
  er_rvf_prv <- run_enrich(dat_rvf_prv, 'RVF_vs_pRV')

  ct_enrich <- bind_rows(
    combine_enrich(er_rvf$up,       'up',   'RVF_vs_NF'),
    combine_enrich(er_rvf$down,     'down', 'RVF_vs_NF'),
    combine_enrich(er_prv$up,       'up',   'pRV_vs_NF'),
    combine_enrich(er_prv$down,     'down', 'pRV_vs_NF'),
    combine_enrich(er_rvf_prv$up,   'up',   'RVF_vs_pRV'),
    combine_enrich(er_rvf_prv$down, 'down', 'RVF_vs_pRV')
  )
  if (!is.null(ct_enrich) && nrow(ct_enrich) > 0) {
    ct_enrich$cell_type <- ct
    mural_enrich_combined[[ct]] <- ct_enrich
  }

  pdf(paste0('./output/Xenium/mural_', ct_slug, '_enrichr_plots.pdf'), width = 12, height = 7)
  for (comp_label in c('RVF_vs_NF', 'pRV_vs_NF', 'RVF_vs_pRV')) {
    er <- switch(comp_label, RVF_vs_NF = er_rvf, pRV_vs_NF = er_prv, RVF_vs_pRV = er_rvf_prv)
    for (dir_label in c('up', 'down')) {
      er_dir <- er[[dir_label]]
      if (!is.null(er_dir)) {
        p1 <- plot_enrich(er_dir$GO_Biological_Process_2023, paste(ct, comp_label, toupper(dir_label), 'GO BP'))
        p2 <- plot_enrich(er_dir$Reactome_2022,              paste(ct, comp_label, toupper(dir_label), 'Reactome'))
        if (!is.null(p1)) print(p1)
        if (!is.null(p2)) print(p2)
      }
    }
  }
  dev.off()
}

mural_all_enrich <- bind_rows(mural_enrich_combined)
write.csv(mural_all_enrich, './output/Xenium/mural_all_enrichr.csv', row.names = FALSE)


##############################################
#### 12. Epi (Epicardium) subclustering
##############################################

epi.obj <- subset(xenium.obj, subset = cell_type_rctd_doublet == 'Epi')
epi.obj <- JoinLayers(epi.obj)

umap_cols <- grep("^umap_\\d+$|^UMAP_\\d+$|^umaporig_\\d+$",
                  colnames(epi.obj@meta.data), value = TRUE)
if (length(umap_cols)) epi.obj@meta.data[umap_cols] <- NULL

epi.obj[["pca"]]      <- NULL
epi.obj[["harmony"]]  <- NULL
epi.obj[["umap"]]     <- NULL
epi.obj[["umap_orig"]] <- NULL

epi.obj <- FindVariableFeatures(epi.obj)
epi.obj <- ScaleData(epi.obj)
epi.obj <- RunPCA(epi.obj, npcs = 30)
epi.obj <- RunHarmony(epi.obj, "patient")

epi.obj <- RunUMAP(epi.obj, reduction = "harmony", dims = 1:30, reduction.name = 'umap')
epi.obj <- RunUMAP(epi.obj, reduction = "pca",     dims = 1:30, reduction.name = 'umap_orig')

epi.obj <- FindNeighbors(epi.obj, reduction = "harmony", dims = 1:30)
epi.obj <- FindClusters(epi.obj, resolution = 0.5)

epi.marks <- FindAllMarkers(epi.obj)
write.csv(epi.marks, './output/Xenium/epi.marks.csv')

DimPlot(epi.obj, group.by = 'seurat_clusters', label = TRUE, repel = TRUE, reduction = 'umap') +
  theme(legend.position = 'none')
DimPlot(epi.obj, group.by = 'patient', reduction = 'umap')

# ── Annotate clusters after reviewing epi.marks.csv ──────────────────────────
#   19 clusters (c0-c18). Heavy contamination: 14 non-epicardial, 5 retained.
#
#   ┌─────────┬──────────────────────┬──────────────────────────────────────────────────────┐
#   │ Cluster │      Cell Type       │            Key Markers / Evidence                    │
#   ├─────────┼──────────────────────┼──────────────────────────────────────────────────────┤
#   │ 0       │ Epi_Mesothelial      │ WT1+1.4, TBX18+1.6, ALDH1A2+1.3, COL1A2+1.1,      │
#   │         │                      │ DCN+1.5, PDGFRA+1.8 — classic epicardial/mesothelial │
#   │ 1       │ drop_EC              │ CDH5+6.5, TIE1+6.7, VWF+2.3 — endothelial           │
#   │ 2       │ drop_CM              │ TNNT2+, RYR2+0.8, TTN+1.0, MYOCD+1.9 — CM           │
#   │ 3       │ drop_Myeloid         │ C1QB+10.2, TYROBP+10.1, CSF1R+4.9, CD163+1.5        │
#   │ 4       │ Epi_EPDCf            │ TBX18+1.6, ALDH1A2+1.6, COL1A1+3.4, FAP+5.1,       │
#   │         │                      │ POSTN+1.6, PI16+9.2 — epicardial-derived fibroblast  │
#   │ 5       │ drop_CM2             │ RYR2+1.9, TTN+2.0, MYH7+1.5, TNNT2+1.2 — CM         │
#   │ 6       │ drop_Pericyte        │ RGS5+3.0, NOTCH3+2.7, CARMN+2.1 — pericyte           │
#   │ 7       │ Epi_ALDH1A2          │ TBX18+1.7, ALDH1A2+2.4, COL1A2+1.6 — epicardial     │
#   │         │                      │ with strong retinoic acid signaling                   │
#   │ 8       │ drop_Pericyte2       │ RGS5+3.1, KCNJ8+10.7, NOTCH3+3.3, ACTA2+2.7         │
#   │ 9       │ drop_VSMC            │ TAGLN+5.5, NOTCH3+3.1, CARMN+2.8, MYH11+3.8          │
#   │ 10      │ drop_TNK             │ IL7R+1.9, CD6+13.2, TOP2A+2.3 — T/NK cells           │
#   │ 11      │ drop_Adipocyte       │ GPD1+14.4, CIDEC+14.2, PLIN4+6.1, ADIPOQ+6.7        │
#   │ 12      │ drop_Mural           │ RGS5+2.6, MYOCD+2.1, CARMN+1.7; NegControl top      │
#   │ 13      │ Epi_Neural           │ SLC8A3+13.6, HDC+12.4, HPGD+9.2; RGS5+0.9 (low),   │
#   │         │                      │ CARMN+0.6 (low) — epicardial/neural crest features    │
#   │ 14      │ Epi_Mesothelial2     │ WT1+2.9, UPK3B+7.2, BMP4+2.7, CLDN1+14.4,           │
#   │         │                      │ LRP2+14.3 — strongest epicardial; mature mesothelium  │
#   │ 15      │ drop_Neuron          │ NRXN1+7.5, NRXN3+4.5, SOX10+11.9 — Schwann/neuron   │
#   │ 16      │ drop_Proliferating   │ KIF11+14.2, BUB1B+14.2, TPX2+11.1 — cycling cells    │
#   │ 17      │ drop_EC2             │ CDH5+2.4, TIE1+3.1, STAB2+14.3 — sinusoidal EC       │
#   │ 18      │ drop_CM3             │ TNNT2+6.1, MYH7+6.0, MYH6+6.4, NPPA+7.2, NPPB+8.0  │
#   └─────────┴──────────────────────┴──────────────────────────────────────────────────────┘

epi_labels <- c(
  '0'  = 'Epi_Mesothelial',
  '1'  = 'drop_EC',
  '2'  = 'drop_CM',
  '3'  = 'drop_Myeloid',
  '4'  = 'Epi_EPDCf',
  '5'  = 'drop_CM2',
  '6'  = 'drop_Pericyte',
  '7'  = 'Epi_ALDH1A2',
  '8'  = 'drop_Pericyte2',
  '9'  = 'drop_VSMC',
  '10' = 'drop_TNK',
  '11' = 'drop_Adipocyte',
  '12' = 'drop_Mural',
  '13' = 'Epi_Neural',
  '14' = 'Epi_Mesothelial2',
  '15' = 'drop_Neuron',
  '16' = 'drop_Proliferating',
  '17' = 'drop_EC2',
  '18' = 'drop_CM3'
)

epi.obj$epi_subtype <- unname(epi_labels[as.character(epi.obj$seurat_clusters)])

DimPlot(epi.obj, group.by = 'epi_subtype', label = TRUE, repel = TRUE, reduction = 'umap') +
  theme(legend.position = 'none')

keep_cells <- !is.na(epi.obj$epi_subtype) & !grepl('^drop_', epi.obj$epi_subtype)
epi.obj.clean <- subset(epi.obj, cells = colnames(epi.obj)[keep_cells])

epi.obj.clean[["pca"]]      <- NULL
epi.obj.clean[["harmony"]]  <- NULL
epi.obj.clean[["umap"]]     <- NULL
epi.obj.clean[["umap_orig"]] <- NULL

umap_cols <- grep("^umap_\\d+$|^UMAP_\\d+$|^umaporig_\\d+$",
                  colnames(epi.obj.clean@meta.data), value = TRUE)
if (length(umap_cols)) epi.obj.clean@meta.data[umap_cols] <- NULL

epi.obj.clean <- JoinLayers(epi.obj.clean)
epi.obj.clean <- FindVariableFeatures(epi.obj.clean)
epi.obj.clean <- ScaleData(epi.obj.clean)
epi.obj.clean <- RunPCA(epi.obj.clean, npcs = 30)
epi.obj.clean <- RunHarmony(epi.obj.clean, "patient")
epi.obj.clean <- RunUMAP(epi.obj.clean, reduction = "harmony", dims = 1:30, reduction.name = 'umap')
epi.obj.clean <- RunUMAP(epi.obj.clean, reduction = "pca",     dims = 1:30, reduction.name = 'umap_orig')

epi.obj.clean <- FindNeighbors(epi.obj.clean, reduction = "harmony", dims = 1:30)
epi.obj.clean <- FindClusters(epi.obj.clean, resolution = 0.5)

epi.obj.clean.mem <- epi.obj.clean
epi.obj.clean.mem[["Xenium"]] <- CreateAssayObject(
  counts = as(GetAssayData(epi.obj.clean, layer = "counts"), "dgCMatrix")
)
epi.clean.marks <- FindAllMarkers(epi.obj.clean.mem)
rm(epi.obj.clean.mem)
write.csv(epi.clean.marks, './output/Xenium/epi.clean.marks.csv')

# ── Annotate Epi clean clusters after reviewing epi.clean.marks.csv ─────────
#   8 clusters (c0-c7). 7 retained, 1 dropped (c6 = multi-lineage contamination).
#
#   Nomenclature note: only c7 (Epi_Mesothelial) retains canonical epicardial
#   identity (WT1+, UPK3B+, BMP4+). Clusters c0-c5 are WT1-negative and have
#   undergone epithelial-to-mesenchymal transition (EMT), making them
#   epicardial-DERIVED cells (EPDCs) rather than true epicardial mesothelium.
#   WT1 is required for epicardial EMT (von Gise et al. 2011) and is rapidly
#   downregulated in invading EPDCs (Hesse et al. 2021). The UMAP separation
#   between c7 and the EPDC clusters reflects this genuine transcriptional gulf.
#   We retain the "Epi_" prefix for all subtypes to denote lineage origin, not
#   current mesothelial identity.
#
#   Literature support:
#   - Hesse et al. (2021) eLife 10:e65921 — scRNA-seq of WT1+ epicardial
#     stromal cells; 11 EpiSC populations in 3 groups post-MI; EMT hierarchy
#   - Knight-Schrijver et al. (2022) Nat Cardiovasc Res 1:1200 — adult vs fetal
#     human epicardium; adult epicardium solely mesothelial (UPK3B+); fibroblast-
#     like EPDCs only in fetal hearts; PRG4/ITLN1 mesothelial-specific markers
#   - von Gise et al. (2011) Circ Res 109:e1 — WT1 required for epicardial EMT
#     via β-catenin and retinoic acid signaling
#   - Cao & Bhatt (2023) Nat Biotechnol 41:1534 — epicardioid scRNA-seq; ALDH1A2+
#     epicardial cells produce retinoic acid received by neighboring CMs;
#     EPDCs differentiate into fibroblasts (TNC, FN1, COL1A1) and mural cells
#   - Bochmann et al. (2010) PLoS One 5:e11834 — ciliated epicardial/mesothelial
#     cells on the cardiac surface; serous pericardium contains ciliated mesothelium
#
#   ┌─────────┬──────────────────────┬──────────────────────────────────────────────────────┐
#   │ Cluster │      Cell Type       │            Key Markers / Evidence                    │
#   ├─────────┼──────────────────────┼──────────────────────────────────────────────────────┤
#   │ 0       │ Epi_EPDCf            │ COL1A1+1.4, FAP+2.5, POSTN+2.5, PI16+6.4,           │
#   │         │                      │ THBS4+27.9; WT1-neg — epicardial-derived fibroblast   │
#   │         │                      │ [Hesse 2021; Knight-Schrijver 2022 — fetal EPDC-fb]  │
#   │ 1       │ Epi_Activated        │ TIMP3+93.9, GPX3+91.4, CXCL8+45.3, CXCL2+32.1,      │
#   │         │                      │ DCN+7.8, PDGFRA+7.7, THY1+14.8; TBX18+ but WT1-neg  │
#   │         │                      │ — chemokine-secreting activated EPDCs                  │
#   │         │                      │ [CXCL8 = key neutrophil attractant in cardiac injury] │
#   │ 2       │ Epi_Quiescent        │ No strong lineage markers +/-; minimal contamination  │
#   │         │                      │ — low-expression baseline quiescent epicardial state   │
#   │ 3       │ Epi_ALDH1A2          │ FSTL3+9.7; 130 up vs 2996 down; sparse transcription  │
#   │         │                      │ — retinoic acid-producing epicardial progenitor state  │
#   │         │                      │ [Cao & Bhatt 2023 — ALDH1A2+ epi produce RA for CMs] │
#   │ 4       │ Epi_Ciliated         │ FOXJ1+2.5 (master cilia TF), PON3+3.1; only 151 up   │
#   │         │                      │ genes; WT1-neg, BMP4-neg — rare ciliated mesothelial   │
#   │         │                      │ [Bochmann 2010 — ciliated cells on cardiac surface]   │
#   │ 5       │ Epi_Metabolic        │ HK2+6.3 (glycolysis), PPRC1+7.8 (PGC-1 coactivator), │
#   │         │                      │ SLC25A44+6.2 — metabolically active epicardial state   │
#   │ 6       │ drop_Mixed           │ RGS5+4.4, NOTCH3+4.4 (mural); CDH5+7.3, TIE1+6.7    │
#   │         │                      │ (EC); HDC+9.3 (mast); multi-lineage contamination     │
#   │ 7       │ Epi_Mesothelial      │ WT1+24.4, UPK3B+18.9, BMP4+23.4, CLDN1+12.1,        │
#   │         │                      │ LRP2+11.8, PRG4+153; strongest epicardial identity    │
#   │         │                      │ [Knight-Schrijver 2022 — UPK3B/PRG4 mesothelial-only] │
#   └─────────┴──────────────────────┴──────────────────────────────────────────────────────┘

epi_clean_labels <- c(
  '0' = 'Epi_EPDCf',
  '1' = 'Epi_Activated',
  '2' = 'Epi_Quiescent',
  '3' = 'Epi_ALDH1A2',
  '4' = 'Epi_Ciliated',
  '5' = 'Epi_Metabolic',
  '6' = 'drop_Mixed',
  '7' = 'Epi_Mesothelial'
)

epi.obj.clean$epi_subtype <- unname(epi_clean_labels[as.character(epi.obj.clean$seurat_clusters)])

# Remove residual contamination cluster
keep_clean <- !is.na(epi.obj.clean$epi_subtype) & !grepl('^drop_', epi.obj.clean$epi_subtype)
epi.obj.clean <- subset(epi.obj.clean, cells = colnames(epi.obj.clean)[keep_clean])

epi.obj.clean[["pca"]]      <- NULL
epi.obj.clean[["harmony"]]  <- NULL
epi.obj.clean[["umap"]]     <- NULL
epi.obj.clean[["umap_orig"]] <- NULL

umap_cols <- grep("^umap_\\d+$|^UMAP_\\d+$|^umaporig_\\d+$",
                  colnames(epi.obj.clean@meta.data), value = TRUE)
if (length(umap_cols)) epi.obj.clean@meta.data[umap_cols] <- NULL

epi.obj.clean <- JoinLayers(epi.obj.clean)
epi.obj.clean <- FindVariableFeatures(epi.obj.clean)
epi.obj.clean <- ScaleData(epi.obj.clean)
epi.obj.clean <- RunPCA(epi.obj.clean, npcs = 30)
epi.obj.clean <- RunHarmony(epi.obj.clean, "patient")
epi.obj.clean <- RunUMAP(epi.obj.clean, reduction = "harmony", dims = 1:30, reduction.name = 'umap')
epi.obj.clean <- RunUMAP(epi.obj.clean, reduction = "pca",     dims = 1:30, reduction.name = 'umap_orig')


DimPlot(epi.obj.clean, group.by = 'epi_subtype', label = TRUE, repel = TRUE, reduction = 'umap') +
  theme(legend.position = 'none')
DimPlot(epi.obj.clean, group.by = 'patient', reduction = 'umap')

saveRDS(epi.obj.clean, './output/Xenium/epi_clean_clean.rds')

pdf('./output/Xenium/Xenium_epi_snUMAP.pdf', width = 5, height = 5)
PlotEmbedding(epi.obj.clean, group.by = 'epi_subtype', reduction = 'umap',
              point_size = .1, plot_under = TRUE, plot_theme = umap_theme() + NoLegend(),
              raster_dpi = 400, raster_scale = 0.5)
dev.off()

pdf('./output/Xenium/Xenium_epi_snUMAP_pt.pdf', width = 5, height = 5)
PlotEmbedding(epi.obj.clean, group.by = 'patient', reduction = 'umap',
              point_size = .1, plot_under = TRUE, plot_theme = umap_theme() + NoLegend(),
              raster_dpi = 400, raster_scale = 0.5)
dev.off()

##############################################
#### Epi subtype proportions by disease state
##############################################

epi_meta <- epi.obj.clean@meta.data %>%
  tibble::rownames_to_column('cell') %>%
  select(cell, patient, group, epi_subtype)

epi_totals <- xenium.obj@meta.data %>%
  tibble::rownames_to_column('cell') %>%
  select(cell, patient, group) %>%
  group_by(patient, group) %>%
  summarise(total = n(), .groups = 'drop')

epi_prop <- epi_meta %>%
  group_by(patient, group, epi_subtype) %>%
  summarise(n = n(), .groups = 'drop') %>%
  left_join(epi_totals, by = c('patient', 'group')) %>%
  mutate(proportion = n / total) %>%
  complete(nesting(patient, group), epi_subtype, fill = list(n = 0, proportion = 0)) %>%
  mutate(group = factor(group, levels = c('NF', 'pRV', 'RVF')))

p_epi_prop <- ggplot(epi_prop, aes(x = group, y = proportion, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.5) +
  stat_compare_means(comparisons = pairwise_comparisons, method = 'wilcox.test', label = 'p.signif') +
  stat_compare_means(method = 'kruskal.test', label.y.npc = 0.95, size = 3) +
  facet_wrap(~ epi_subtype, scales = 'free_y') +
  scale_x_discrete(limits = c('NF', 'pRV', 'RVF')) +
  theme_bw() +
  theme(legend.position = 'none',
        strip.text = element_text(face = 'bold'),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = 'Proportion of all cells', title = 'Epicardium subtype distribution by group')

pdf('./output/Xenium/epi_proportions.pdf', width = 8, height = 6)
print(p_epi_prop)
dev.off()

##############################################
#### Epi fractional area coverage by disease state
##############################################

cell_areas <- read.csv('./output/Xenium/cell_areas_2d.csv', row.names = 1)

total_area <- xenium.obj@meta.data %>%
  tibble::rownames_to_column('cell') %>%
  select(cell, patient, group) %>%
  left_join(cell_areas, by = 'cell') %>%
  group_by(patient, group) %>%
  summarise(total_area = sum(cell_area_2d, na.rm = TRUE), .groups = 'drop')

epi_area_meta <- epi.obj.clean@meta.data %>%
  tibble::rownames_to_column('cell') %>%
  select(cell, patient, group, epi_subtype) %>%
  left_join(cell_areas, by = 'cell')

epi_area_prop <- epi_area_meta %>%
  group_by(patient, group, epi_subtype) %>%
  summarise(epi_area = sum(cell_area_2d, na.rm = TRUE), .groups = 'drop') %>%
  left_join(total_area, by = c('patient', 'group')) %>%
  mutate(frac_area = epi_area / total_area) %>%
  complete(nesting(patient, group), epi_subtype, fill = list(epi_area = 0, frac_area = 0)) %>%
  mutate(group = factor(group, levels = c('NF', 'pRV', 'RVF')))

p_epi_area <- ggplot(epi_area_prop, aes(x = group, y = frac_area, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.5) +
  stat_compare_means(comparisons = pairwise_comparisons, method = 'wilcox.test', label = 'p.signif') +
  stat_compare_means(method = 'kruskal.test', label.y.npc = 0.95, size = 3) +
  facet_wrap(~ epi_subtype, scales = 'free_y') +
  scale_x_discrete(limits = c('NF', 'pRV', 'RVF')) +
  theme_bw() +
  theme(legend.position = 'none',
        strip.text = element_text(face = 'bold'),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = 'Fractional area coverage', title = 'Epicardium subtype fractional area by group')

pdf('./output/Xenium/epi_area_coverage.pdf', width = 8, height = 6)
print(p_epi_area)
dev.off()

##############################################
#### Pseudobulk DESeq2 by Epi subtype
##############################################

epi_pseudo_meta <- unique(epi.obj.clean@meta.data[, c('patient', 'group', 'epi_subtype')])
epi_pseudo_meta$patient     <- as.character(epi_pseudo_meta$patient)
epi_pseudo_meta$epi_subtype <- as.character(epi_pseudo_meta$epi_subtype)
epi_pseudo_meta$pseudobulk_id <- paste0(epi_pseudo_meta$patient, '..', epi_pseudo_meta$epi_subtype)
epi_pseudo_meta$group <- factor(epi_pseudo_meta$group, levels = c('NF', 'pRV', 'RVF'))
rownames(epi_pseudo_meta) <- epi_pseudo_meta$pseudobulk_id

epi.obj.clean$pseudobulk_id <- paste0(
  as.character(epi.obj.clean$patient), '..',
  as.character(epi.obj.clean$epi_subtype)
)

epi_pseudo_counts <- AggregateExpression(
  epi.obj.clean,
  assays = 'Xenium',
  group.by = 'pseudobulk_id',
  slot = 'counts',
  return.seurat = FALSE
)$Xenium

colnames(epi_pseudo_counts) <- gsub('-', '_', sub('^g', '', colnames(epi_pseudo_counts)))

epi_cell_types <- unique(epi_pseudo_meta$epi_subtype)
epi_deseq_results <- list()

for (ct in epi_cell_types) {
  ct_samples <- epi_pseudo_meta$pseudobulk_id[epi_pseudo_meta$epi_subtype == ct]
  ct_samples <- ct_samples[ct_samples %in% colnames(epi_pseudo_counts)]
  if (length(ct_samples) < 2) next

  ct_counts <- epi_pseudo_counts[, ct_samples, drop = FALSE]
  ct_meta   <- epi_pseudo_meta[ct_samples, , drop = FALSE]

  if (length(unique(ct_meta$group)) < 2) next

  keep <- rowSums(ct_counts >= 5) >= 2
  ct_counts <- ct_counts[keep, , drop = FALSE]
  if (nrow(ct_counts) < 10) next

  dds <- DESeqDataSetFromMatrix(countData = ct_counts, colData = ct_meta, design = ~ group)
  dds <- tryCatch(DESeq(dds), error = function(e) NULL)
  if (is.null(dds)) next

  for (comp in list(c('RVF', 'NF'), c('pRV', 'NF'), c('RVF', 'pRV'))) {
    comp_label <- paste0(comp[1], '_vs_', comp[2])
    res <- tryCatch(
      results(dds, contrast = c('group', comp[1], comp[2])) %>%
        as.data.frame() %>%
        tibble::rownames_to_column('gene') %>%
        mutate(cell_type = ct, comparison = comp_label) %>%
        arrange(padj),
      error = function(e) NULL
    )
    if (!is.null(res)) epi_deseq_results[[paste0(ct, '_', comp_label)]] <- res
  }
}

epi_all_deseq <- bind_rows(epi_deseq_results)
write.csv(epi_all_deseq, './output/Xenium/epi_pseudobulk_deseq2.csv', row.names = FALSE)

epi_sig_summary <- epi_all_deseq %>%
  filter(padj < 0.05, abs(log2FoldChange) > 0.5) %>%
  group_by(cell_type, comparison) %>%
  summarise(n_up = sum(log2FoldChange > 0), n_down = sum(log2FoldChange < 0), .groups = 'drop')
print(epi_sig_summary)

##############################################
#### Epi subtypes: Volcano + GO/Reactome
##############################################

epi_enrich_combined <- list()

for (ct in unique(epi_all_deseq$cell_type)) {
  ct_slug <- gsub('_', '-', ct)
  cat('Epi processing:', ct, '\n')

  dat_rvf     <- epi_all_deseq %>% filter(cell_type == ct, comparison == 'RVF_vs_NF')  %>% cap_pval()
  dat_prv     <- epi_all_deseq %>% filter(cell_type == ct, comparison == 'pRV_vs_NF')  %>% cap_pval()
  dat_rvf_prv <- epi_all_deseq %>% filter(cell_type == ct, comparison == 'RVF_vs_pRV') %>% cap_pval()

  make_volcano_pdf(dat_rvf,     paste0('./output/Xenium/epi_', ct_slug, '_volcano_RVF_vs_NF.pdf'),  paste(ct, ': RVF vs NF'))
  make_volcano_pdf(dat_prv,     paste0('./output/Xenium/epi_', ct_slug, '_volcano_pRV_vs_NF.pdf'),  paste(ct, ': pRV vs NF'))
  make_volcano_pdf(dat_rvf_prv, paste0('./output/Xenium/epi_', ct_slug, '_volcano_RVF_vs_pRV.pdf'), paste(ct, ': RVF vs pRV'))

  run_enrich <- function(dat, comp) {
    sig  <- dat %>% filter(padj < 0.05, abs(log2FoldChange) > 0.5)
    up   <- sig %>% filter(log2FoldChange > 0) %>% pull(gene)
    down <- sig %>% filter(log2FoldChange < 0) %>% pull(gene)
    list(
      up   = if (length(up)   > 0) enrichr(up,   enrichR_dbs) else NULL,
      down = if (length(down) > 0) enrichr(down, enrichR_dbs) else NULL
    )
  }

  er_rvf     <- run_enrich(dat_rvf,     'RVF_vs_NF')
  er_prv     <- run_enrich(dat_prv,     'pRV_vs_NF')
  er_rvf_prv <- run_enrich(dat_rvf_prv, 'RVF_vs_pRV')

  ct_enrich <- bind_rows(
    combine_enrich(er_rvf$up,       'up',   'RVF_vs_NF'),
    combine_enrich(er_rvf$down,     'down', 'RVF_vs_NF'),
    combine_enrich(er_prv$up,       'up',   'pRV_vs_NF'),
    combine_enrich(er_prv$down,     'down', 'pRV_vs_NF'),
    combine_enrich(er_rvf_prv$up,   'up',   'RVF_vs_pRV'),
    combine_enrich(er_rvf_prv$down, 'down', 'RVF_vs_pRV')
  )
  if (!is.null(ct_enrich) && nrow(ct_enrich) > 0) {
    ct_enrich$cell_type <- ct
    epi_enrich_combined[[ct]] <- ct_enrich
  }

  pdf(paste0('./output/Xenium/epi_', ct_slug, '_enrichr_plots.pdf'), width = 12, height = 7)
  for (comp_label in c('RVF_vs_NF', 'pRV_vs_NF', 'RVF_vs_pRV')) {
    er <- switch(comp_label, RVF_vs_NF = er_rvf, pRV_vs_NF = er_prv, RVF_vs_pRV = er_rvf_prv)
    for (dir_label in c('up', 'down')) {
      er_dir <- er[[dir_label]]
      if (!is.null(er_dir)) {
        p1 <- plot_enrich(er_dir$GO_Biological_Process_2023, paste(ct, comp_label, toupper(dir_label), 'GO BP'))
        p2 <- plot_enrich(er_dir$Reactome_2022,              paste(ct, comp_label, toupper(dir_label), 'Reactome'))
        if (!is.null(p1)) print(p1)
        if (!is.null(p2)) print(p2)
      }
    }
  }
  dev.off()
}

epi_all_enrich <- bind_rows(epi_enrich_combined)
write.csv(epi_all_enrich, './output/Xenium/epi_all_enrichr.csv', row.names = FALSE)


##############################################
#### 13. Adipo (Adipocyte) subclustering
##############################################

adipo.obj <- subset(xenium.obj, subset = cell_type_rctd_doublet == 'Adipo')
adipo.obj <- JoinLayers(adipo.obj)

umap_cols <- grep("^umap_\\d+$|^UMAP_\\d+$|^umaporig_\\d+$",
                  colnames(adipo.obj@meta.data), value = TRUE)
if (length(umap_cols)) adipo.obj@meta.data[umap_cols] <- NULL

adipo.obj[["pca"]]      <- NULL
adipo.obj[["harmony"]]  <- NULL
adipo.obj[["umap"]]     <- NULL
adipo.obj[["umap_orig"]] <- NULL

adipo.obj <- FindVariableFeatures(adipo.obj)
adipo.obj <- ScaleData(adipo.obj)
adipo.obj <- RunPCA(adipo.obj, npcs = 30)
adipo.obj <- RunHarmony(adipo.obj, "patient")

adipo.obj <- RunUMAP(adipo.obj, reduction = "harmony", dims = 1:30, reduction.name = 'umap')
adipo.obj <- RunUMAP(adipo.obj, reduction = "pca",     dims = 1:30, reduction.name = 'umap_orig')

adipo.obj <- FindNeighbors(adipo.obj, reduction = "harmony", dims = 1:30)
adipo.obj <- FindClusters(adipo.obj, resolution = 0.5)

adipo.marks <- FindAllMarkers(adipo.obj)
write.csv(adipo.marks, './output/Xenium/adipo.marks.csv')

DimPlot(adipo.obj, group.by = 'seurat_clusters', label = TRUE, repel = TRUE, reduction = 'umap') +
  theme(legend.position = 'none')
DimPlot(adipo.obj, group.by = 'patient', reduction = 'umap')

# ── Annotate clusters after reviewing adipo.marks.csv ────────────────────────
#   20 clusters (c0-c19). Heavy contamination: 15 non-adipocyte, 5 retained.
#
#   ┌─────────┬──────────────────────┬──────────────────────────────────────────────────────┐
#   │ Cluster │      Cell Type       │            Key Markers / Evidence                    │
#   ├─────────┼──────────────────────┼──────────────────────────────────────────────────────┤
#   │ 0       │ drop_CM              │ TNNT2+3.0, MYH7+2.9, NPPA+4.3, NPPB+5.0 — CM       │
#   │ 1       │ Adipocyte            │ All 11 adipo markers+: ADIPOQ+2.8, PLIN4+3.2,       │
#   │         │                      │ LEP+5.3, GPD1+3.2, PCK1+3.8 — strongest adipocyte   │
#   │ 2       │ drop_FB              │ COL1A2+1.7, DCN+1.5, PDGFRA+1.8 — fibroblast        │
#   │ 3       │ drop_EC              │ CDH5+6.4, TIE1+6.9, VWF+2.6 — endothelial           │
#   │ 4       │ drop_Pericyte        │ RGS5+3.2, NOTCH3+2.8, CARMN+1.6 — pericyte          │
#   │ 5       │ drop_FB2             │ COL1A1+4.1, FAP+6.4, POSTN+2.5 — activated FB       │
#   │ 6       │ drop_Myeloid         │ C1QB+10.3, TYROBP+12.0, CSF1R+5.4 — macrophage      │
#   │ 7       │ Adipocyte2           │ PLIN1+2.4, CIDEC+2.5, GPD1+2.1, LIPE+2.4;           │
#   │         │                      │ lower ADIPOQ/PPARG than c1 — mature adipocyte variant│
#   │ 8       │ Adipo_Fibro          │ 10 adipo markers+ (moderate) + DCN+1.1, PDGFRA+1.4, │
#   │         │                      │ TAGLN+1.9 — fibro-adipogenic progenitor features     │
#   │ 9       │ drop_CM2             │ RYR2+2.3, TTN+2.4, MYOCD+3.4 — cardiomyocyte        │
#   │ 10      │ drop_Pericyte2       │ KCNJ8+13.6, RGS5+3.7, NOTCH3+3.8 — KCNJ8+ pericyte │
#   │ 11      │ drop_FB3             │ COL1A1+1.8, COL1A2+1.7, DCN+0.6 — fibroblast        │
#   │ 12      │ drop_VSMC            │ TAGLN+4.8, MYH11+4.7, NOTCH3+3.5, CARMN+2.3 — VSMC  │
#   │ 13      │ Adipo_Vascular       │ All 11 adipo markers+ + PECAM1+2.6, VWF+2.2 —       │
#   │         │                      │ adipocyte with vascular EC features                   │
#   │ 14      │ Adipo_Immune         │ 10 adipo markers+ + CD163+3.0, CSF1R+2.2 —           │
#   │         │                      │ adipocyte with lipid-associated macrophage features   │
#   │ 15      │ drop_TNK             │ RASGRP1+16.4, RUNX3+16.2, IL7R+5.1 — T/NK cells     │
#   │ 16      │ drop_Neuron          │ NRXN1+8.2, SOX10+15.2 — Schwann/neuron               │
#   │ 17      │ drop_EC2             │ CDH5+2.5, TIE1+2.7 — endothelial                     │
#   │ 18      │ drop_CM3             │ TNNT2+1.3, MYH7+1.2; NegControls top — low-quality   │
#   │ 19      │ drop_Proliferating   │ TPX2+16.3, TOP2A+5.3, MKI67+4.7 — cycling cells      │
#   └─────────┴──────────────────────┴──────────────────────────────────────────────────────┘

adipo_labels <- c(
  '0'  = 'drop_CM',
  '1'  = 'Adipocyte',
  '2'  = 'drop_FB',
  '3'  = 'drop_EC',
  '4'  = 'drop_Pericyte',
  '5'  = 'drop_FB2',
  '6'  = 'drop_Myeloid',
  '7'  = 'Adipocyte2',
  '8'  = 'Adipo_Fibro',
  '9'  = 'drop_CM2',
  '10' = 'drop_Pericyte2',
  '11' = 'drop_FB3',
  '12' = 'drop_VSMC',
  '13' = 'Adipo_Vascular',
  '14' = 'Adipo_Immune',
  '15' = 'drop_TNK',
  '16' = 'drop_Neuron',
  '17' = 'drop_EC2',
  '18' = 'drop_CM3',
  '19' = 'drop_Proliferating'
)

adipo.obj$adipo_subtype <- unname(adipo_labels[as.character(adipo.obj$seurat_clusters)])

DimPlot(adipo.obj, group.by = 'adipo_subtype', label = TRUE, repel = TRUE, reduction = 'umap') +
  theme(legend.position = 'none')

keep_cells <- !is.na(adipo.obj$adipo_subtype) & !grepl('^drop_', adipo.obj$adipo_subtype)
adipo.obj.clean <- subset(adipo.obj, cells = colnames(adipo.obj)[keep_cells])

adipo.obj.clean[["pca"]]      <- NULL
adipo.obj.clean[["harmony"]]  <- NULL
adipo.obj.clean[["umap"]]     <- NULL
adipo.obj.clean[["umap_orig"]] <- NULL

umap_cols <- grep("^umap_\\d+$|^UMAP_\\d+$|^umaporig_\\d+$",
                  colnames(adipo.obj.clean@meta.data), value = TRUE)
if (length(umap_cols)) adipo.obj.clean@meta.data[umap_cols] <- NULL

adipo.obj.clean <- JoinLayers(adipo.obj.clean)
adipo.obj.clean <- FindVariableFeatures(adipo.obj.clean)
adipo.obj.clean <- ScaleData(adipo.obj.clean)
adipo.obj.clean <- RunPCA(adipo.obj.clean, npcs = 30)
adipo.obj.clean <- RunHarmony(adipo.obj.clean, "patient")
adipo.obj.clean <- RunUMAP(adipo.obj.clean, reduction = "harmony", dims = 1:30, reduction.name = 'umap')
adipo.obj.clean <- RunUMAP(adipo.obj.clean, reduction = "pca",     dims = 1:30, reduction.name = 'umap_orig')

adipo.obj.clean <- FindNeighbors(adipo.obj.clean, reduction = "harmony", dims = 1:30)
adipo.obj.clean <- FindClusters(adipo.obj.clean, resolution = 0.5)

adipo.obj.clean.mem <- adipo.obj.clean
adipo.obj.clean.mem[["Xenium"]] <- CreateAssayObject(
  counts = as(GetAssayData(adipo.obj.clean, layer = "counts"), "dgCMatrix")
)
adipo.clean.marks <- FindAllMarkers(adipo.obj.clean.mem)
rm(adipo.obj.clean.mem)
write.csv(adipo.clean.marks, './output/Xenium/adipo.clean.marks.csv')

# ── Annotate Adipo clean clusters after reviewing adipo.clean.marks.csv ─────
#   11 clusters (c0-c10). 5 retained, 6 dropped (residual FB/EC/mural/macrophage).
#   Initial c13 (Adipo_Vascular) and c14 (Adipo_Immune) turned out to be
#   contamination — re-clustered into pure EC (c8) and macrophage (c6/c10).
#
#   ┌─────────┬──────────────────────┬──────────────────────────────────────────────────────┐
#   │ Cluster │      Cell Type       │            Key Markers / Evidence                    │
#   ├─────────┼──────────────────────┼──────────────────────────────────────────────────────┤
#   │ 0       │ Adipo_Mature         │ All 11 adipo markers: ADIPOQ+79.5, LPL+30.5,        │
#   │         │                      │ PLIN4+39.1, PCK1+0.9; CM ambient (MYH7+125)          │
#   │ 1       │ Adipo_Quiescent      │ Moderate PLIN1+0.4, FABP4+0.4, GPD1+0.6, LIPE+0.6; │
#   │         │                      │ ADIPOQ−113, PLIN4−200 — low-expression resting state  │
#   │ 2       │ Adipo_Sparse         │ 284 up genes; TBL1Y+2.9, APLF+2.8; all major adipo  │
#   │         │                      │ markers negative — possible pre-adipocyte             │
#   │ 3       │ Adipo_Sparse2        │ 419 up genes; TYRO3+4.0; ADIPOQ/PLIN4/PPARG all     │
#   │         │                      │ negative — possible pre-adipocyte                     │
#   │ 4       │ drop_FB              │ DCN+114.4, FN1+32.0, PDGFRA+19.2 — fibroblast        │
#   │ 5       │ drop_Mixed           │ SERPINE1+88.7; CDH5+5.9, TIE1+6.3 (EC), RGS5 (mural)│
#   │         │                      │ — all adipo markers negative; multi-lineage           │
#   │ 6       │ drop_Macrophage      │ CSF1R+5.8, MPEG1+9.6, LILRB2+7.2, MAF+6.5           │
#   │ 7       │ drop_Mural           │ RGS5+3.9; ADIPOQ−155, PPARG−45 — pericyte/CM         │
#   │ 8       │ drop_EC              │ PECAM1+13.6, VWF+23.4, FLT1+9.7 — endothelial        │
#   │ 9       │ Adipo_Lipogenic      │ ELOVL6+6.9 (fatty acid elongase), SAA1+7.1;          │
#   │         │                      │ moderate adipo reduction — active lipogenesis state    │
#   │ 10      │ drop_LAM             │ LYVE1+35.7, CD163+15.9, FOLR2+12.5, CSF1R+5.3 —     │
#   │         │                      │ lipid-associated macrophage                            │
#   └─────────┴──────────────────────┴──────────────────────────────────────────────────────┘

adipo_clean_labels <- c(
  '0'  = 'Adipo_Mature',
  '1'  = 'Adipo_Quiescent',
  '2'  = 'Adipo_Sparse',
  '3'  = 'Adipo_Sparse2',
  '4'  = 'drop_FB',
  '5'  = 'drop_Mixed',
  '6'  = 'drop_Macrophage',
  '7'  = 'drop_Mural',
  '8'  = 'drop_EC',
  '9'  = 'Adipo_Lipogenic',
  '10' = 'drop_LAM'
)

adipo.obj.clean$adipo_subtype <- unname(adipo_clean_labels[as.character(adipo.obj.clean$seurat_clusters)])

# Remove residual contamination — subset to adipo.obj.clean.clean
keep_clean2 <- !is.na(adipo.obj.clean$adipo_subtype) & !grepl('^drop_', adipo.obj.clean$adipo_subtype)
adipo.obj.clean.clean <- subset(adipo.obj.clean, cells = colnames(adipo.obj.clean)[keep_clean2])

adipo.obj.clean.clean[["pca"]]      <- NULL
adipo.obj.clean.clean[["harmony"]]  <- NULL
adipo.obj.clean.clean[["umap"]]     <- NULL
adipo.obj.clean.clean[["umap_orig"]] <- NULL

umap_cols <- grep("^umap_\\d+$|^UMAP_\\d+$|^umaporig_\\d+$",
                  colnames(adipo.obj.clean.clean@meta.data), value = TRUE)
if (length(umap_cols)) adipo.obj.clean.clean@meta.data[umap_cols] <- NULL

adipo.obj.clean.clean <- JoinLayers(adipo.obj.clean.clean)
adipo.obj.clean.clean <- FindVariableFeatures(adipo.obj.clean.clean)
adipo.obj.clean.clean <- ScaleData(adipo.obj.clean.clean)
adipo.obj.clean.clean <- RunPCA(adipo.obj.clean.clean, npcs = 30)
adipo.obj.clean.clean <- RunHarmony(adipo.obj.clean.clean, "patient")
adipo.obj.clean.clean <- RunUMAP(adipo.obj.clean.clean, reduction = "harmony", dims = 1:30, reduction.name = 'umap')
adipo.obj.clean.clean <- RunUMAP(adipo.obj.clean.clean, reduction = "pca",     dims = 1:30, reduction.name = 'umap_orig')


DimPlot(adipo.obj.clean.clean, group.by = 'patient', reduction = 'umap')
DimPlot(adipo.obj.clean.clean, group.by = 'adipo_subtype', label = TRUE, repel = TRUE, reduction = 'umap') +
  theme(legend.position = 'none')


# ── Adipocyte subtype annotation (after triple-clean) ─────────────────────────
#   5 retained clusters from adipo.obj.clean.clean. Labels carried forward from
#   round-2 (adipo_clean_labels) because the triple-clean only removed contamination
#   clusters — the 5 adipocyte identities remain unchanged.
#
#   Literature support:
#
#   ┌─────────┬──────────────────────┬──────────────────────────────────────────────────────┐
#   │ Cluster │      Cell Type       │            Key Markers / Evidence                    │
#   ├─────────┼──────────────────────┼──────────────────────────────────────────────────────┤
#   │ 0       │ Adipo_Mature         │ ADIPOQ+79.5, PLIN4+39.1, LPL+30.5, CD36+117,        │
#   │         │                      │ MEST+26.2, TIMP4+59.3; canonical mature adipocyte.    │
#   │         │                      │ Maps to AdipoPLIN / hAd3 (Emont et al. 2022) with    │
#   │         │                      │ highest ADIPOQ and CFD expression [1,2].              │
#   │ 1       │ Adipo_Quiescent      │ Near-universal PLIN1/FABP4/GPD1/LIPE (pct1≈1.0) but  │
#   │         │                      │ max FC <2.7; ADIPOQ−113, PLIN4−200 — low-expression   │
#   │         │                      │ resting state. Corresponds to hAd1/hAd2 in Emont      │
#   │         │                      │ et al. 2022, which have "relatively few specific       │
#   │         │                      │ markers" [1]. Quiescent adipocyte states described     │
#   │         │                      │ in Sárvári et al. 2021 [3].                            │
#   │ 2       │ Adipo_Sparse         │ TBL1Y+2.9, KLHL23+3.3, CACNA1C+2.1 (calcium channel); │
#   │         │                      │ all major adipo markers strongly negative (ADIPOQ−137,  │
#   │         │                      │ PLIN4−247, PPARG−37) but housekeeping lipid genes       │
#   │         │                      │ retained (PLIN1/FABP4/CIDEC pct1=1.0). Transcript-     │
#   │         │                      │ sparse committed adipocyte, not preadipocyte (which     │
#   │         │                      │ would lack PLIN1/FABP4). Ion channel signature          │
#   │         │                      │ (CACNA1C/CACNB2) distinguishes from c3 [4,5].          │
#   │ 3       │ Adipo_Sparse2        │ TYRO3+4.0 (TAM receptor), RND3+2.6, HIPK1+2.5;        │
#   │         │                      │ ADIPOQ−135, PLIN4−210, PPARG−35 but FABP4/LIPE         │
#   │         │                      │ retained (pct1=1.0). Only 6/30 top genes shared with    │
#   │         │                      │ c2 — distinct kinase/signaling program. TYRO3 (TAM      │
#   │         │                      │ family) involved in efferocytosis and anti-inflammatory  │
#   │         │                      │ signaling [6]. Both c2 and c3 may represent              │
#   │         │                      │ dedifferentiated or small adipocyte states [7].          │
#   │ 9       │ Adipo_Lipogenic      │ ELOVL6+6.9 (fatty acid elongase), SAA1+7.1, SAA2+6.8;  │
#   │         │                      │ maps to AdipoSAA subpopulation (Emont et al. 2022),     │
#   │         │                      │ which groups with lipid-scavenging adipocytes (hAd3/4).  │
#   │         │                      │ SAA1/2 are established inflammatory adipokines (Yang     │
#   │         │                      │ et al. 2006) [1,8,9].                                    │
#   └─────────┴──────────────────────┴──────────────────────────────────────────────────────┘
#
#   Adipo_Sparse vs Adipo_Sparse2: these are NOT the same population. They share
#   the "low adipocyte marker" phenotype (reduced ADIPOQ/PLIN4/LPL/PPARG with
#   retained PLIN1/FABP4/CIDEC/LIPE), confirming committed adipocyte identity, but
#   their top marker programs are distinct: c2 = TBL1Y/KLHL23/CACNA1C (Y-linked +
#   ion channels), c3 = TYRO3/RND3/HIPK1 (receptor kinase + cytoskeletal signaling).
#
#   References:
#   [1] Emont et al. (2022) Nature 603:926-933 — single-cell atlas of human/mouse
#       white adipose tissue; defines AdipoPLIN, AdipoLEP, AdipoSAA subtypes and
#       hAd1-hAd7 subpopulations.
#   [2] Sun et al. (2024) ATVB — pericoronary epicardial adipose tissue snRNA-seq;
#       5 adipocyte subclusters (Ad0-Ad4): Ad0/Ad1 lipogenic, Ad3 lipolytic, Ad4
#       senescent.
#   [3] Sárvári et al. (2021) Cell Metab 33:437-453 — spatial mapping of human
#       adipocyte subpopulations with distinct insulin sensitivities; quiescent and
#       active states.
#   [4] Corvera (2021) Nat Rev Endocrinol 17:726-735 — cellular heterogeneity in
#       adipose tissues; preadipocyte-to-adipocyte commitment hierarchy.
#   [5] Min et al. (2019) PNAS 116:17970-17979 — diverse repertoire of human
#       adipocyte subtypes from transcriptionally distinct mesenchymal progenitors.
#   [6] Lemke (2013) Cold Spring Harb Perspect Biol 5:a009076 — biology of TAM
#       receptors (TYRO3/AXL/MERTK) in efferocytosis and inflammation.
#   [7] Hildreth et al. (2021) Nat Immunol 22:639-653 — single-cell sequencing of
#       human WAT identifies preadipocytes, committed progenitors, and immune states.
#   [8] Yang et al. (2006) PLoS Med 3:e287 — acute-phase SAA as inflammatory
#       adipokine linking obesity and metabolic complications.
#   [9] Landgraeber et al. (2024) JAHA — conserved human adipose subpopulations
#       using targeted snRNA-seq; AdipoPLIN, AdipoLEP, AdipoSAA conservation.

saveRDS(adipo.obj.clean.clean, './output/Xenium/adipo_clean_clean.rds')

pdf('./output/Xenium/Xenium_adipo_snUMAP.pdf', width = 5, height = 5)
PlotEmbedding(adipo.obj.clean.clean, group.by = 'adipo_subtype', reduction = 'umap',
              point_size = .1, plot_under = TRUE, plot_theme = umap_theme() + NoLegend(),
              raster_dpi = 400, raster_scale = 0.5)
dev.off()

pdf('./output/Xenium/Xenium_adipo_snUMAP_pt.pdf', width = 5, height = 5)
PlotEmbedding(adipo.obj.clean.clean, group.by = 'patient', reduction = 'umap',
              point_size = .1, plot_under = TRUE, plot_theme = umap_theme() + NoLegend(),
              raster_dpi = 400, raster_scale = 0.5)
dev.off()

##############################################
#### Adipo subtype proportions by disease state
##############################################

adipo_meta <- adipo.obj.clean.clean@meta.data %>%
  tibble::rownames_to_column('cell') %>%
  select(cell, patient, group, adipo_subtype)

adipo_totals <- xenium.obj@meta.data %>%
  tibble::rownames_to_column('cell') %>%
  select(cell, patient, group) %>%
  group_by(patient, group) %>%
  summarise(total = n(), .groups = 'drop')

adipo_prop <- adipo_meta %>%
  group_by(patient, group, adipo_subtype) %>%
  summarise(n = n(), .groups = 'drop') %>%
  left_join(adipo_totals, by = c('patient', 'group')) %>%
  mutate(proportion = n / total) %>%
  complete(nesting(patient, group), adipo_subtype, fill = list(n = 0, proportion = 0)) %>%
  mutate(group = factor(group, levels = c('NF', 'pRV', 'RVF')))

p_adipo_prop <- ggplot(adipo_prop, aes(x = group, y = proportion, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.5) +
  stat_compare_means(comparisons = pairwise_comparisons, method = 'wilcox.test', label = 'p.signif') +
  stat_compare_means(method = 'kruskal.test', label.y.npc = 0.95, size = 3) +
  facet_wrap(~ adipo_subtype, scales = 'free_y') +
  scale_x_discrete(limits = c('NF', 'pRV', 'RVF')) +
  theme_bw() +
  theme(legend.position = 'none',
        strip.text = element_text(face = 'bold'),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = 'Proportion of all cells', title = 'Adipocyte subtype distribution by group')

pdf('./output/Xenium/adipo_proportions.pdf', width = 8, height = 6)
print(p_adipo_prop)
dev.off()

##############################################
#### Adipo fractional area coverage by disease state
##############################################

cell_areas <- read.csv('./output/Xenium/cell_areas_2d.csv', row.names = 1)

total_area <- xenium.obj@meta.data %>%
  tibble::rownames_to_column('cell') %>%
  select(cell, patient, group) %>%
  left_join(cell_areas, by = 'cell') %>%
  group_by(patient, group) %>%
  summarise(total_area = sum(cell_area_2d, na.rm = TRUE), .groups = 'drop')

adipo_area_meta <- adipo.obj.clean.clean@meta.data %>%
  tibble::rownames_to_column('cell') %>%
  select(cell, patient, group, adipo_subtype) %>%
  left_join(cell_areas, by = 'cell')

adipo_area_prop <- adipo_area_meta %>%
  group_by(patient, group, adipo_subtype) %>%
  summarise(adipo_area = sum(cell_area_2d, na.rm = TRUE), .groups = 'drop') %>%
  left_join(total_area, by = c('patient', 'group')) %>%
  mutate(frac_area = adipo_area / total_area) %>%
  complete(nesting(patient, group), adipo_subtype, fill = list(adipo_area = 0, frac_area = 0)) %>%
  mutate(group = factor(group, levels = c('NF', 'pRV', 'RVF')))

p_adipo_area <- ggplot(adipo_area_prop, aes(x = group, y = frac_area, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.5) +
  stat_compare_means(comparisons = pairwise_comparisons, method = 'wilcox.test', label = 'p.signif') +
  stat_compare_means(method = 'kruskal.test', label.y.npc = 0.95, size = 3) +
  facet_wrap(~ adipo_subtype, scales = 'free_y') +
  scale_x_discrete(limits = c('NF', 'pRV', 'RVF')) +
  theme_bw() +
  theme(legend.position = 'none',
        strip.text = element_text(face = 'bold'),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = 'Fractional area coverage', title = 'Adipocyte subtype fractional area by group')

pdf('./output/Xenium/adipo_area_coverage.pdf', width = 8, height = 6)
print(p_adipo_area)
dev.off()

##############################################
#### Pseudobulk DESeq2 by Adipo subtype
##############################################

adipo_pseudo_meta <- unique(adipo.obj.clean.clean@meta.data[, c('patient', 'group', 'adipo_subtype')])
adipo_pseudo_meta$patient      <- as.character(adipo_pseudo_meta$patient)
adipo_pseudo_meta$adipo_subtype <- as.character(adipo_pseudo_meta$adipo_subtype)
adipo_pseudo_meta$pseudobulk_id <- paste0(adipo_pseudo_meta$patient, '..', adipo_pseudo_meta$adipo_subtype)
adipo_pseudo_meta$group <- factor(adipo_pseudo_meta$group, levels = c('NF', 'pRV', 'RVF'))
rownames(adipo_pseudo_meta) <- adipo_pseudo_meta$pseudobulk_id

adipo.obj.clean.clean$pseudobulk_id <- paste0(
  as.character(adipo.obj.clean.clean$patient), '..',
  as.character(adipo.obj.clean.clean$adipo_subtype)
)

adipo_pseudo_counts <- AggregateExpression(
  adipo.obj.clean.clean,
  assays = 'Xenium',
  group.by = 'pseudobulk_id',
  slot = 'counts',
  return.seurat = FALSE
)$Xenium

colnames(adipo_pseudo_counts) <- gsub('-', '_', sub('^g', '', colnames(adipo_pseudo_counts)))

adipo_cell_types <- unique(adipo_pseudo_meta$adipo_subtype)
adipo_deseq_results <- list()

for (ct in adipo_cell_types) {
  ct_samples <- adipo_pseudo_meta$pseudobulk_id[adipo_pseudo_meta$adipo_subtype == ct]
  ct_samples <- ct_samples[ct_samples %in% colnames(adipo_pseudo_counts)]
  if (length(ct_samples) < 2) next

  ct_counts <- adipo_pseudo_counts[, ct_samples, drop = FALSE]
  ct_meta   <- adipo_pseudo_meta[ct_samples, , drop = FALSE]

  if (length(unique(ct_meta$group)) < 2) next

  keep <- rowSums(ct_counts >= 5) >= 2
  ct_counts <- ct_counts[keep, , drop = FALSE]
  if (nrow(ct_counts) < 10) next

  dds <- DESeqDataSetFromMatrix(countData = ct_counts, colData = ct_meta, design = ~ group)
  dds <- tryCatch(DESeq(dds), error = function(e) NULL)
  if (is.null(dds)) next

  for (comp in list(c('RVF', 'NF'), c('pRV', 'NF'), c('RVF', 'pRV'))) {
    comp_label <- paste0(comp[1], '_vs_', comp[2])
    res <- tryCatch(
      results(dds, contrast = c('group', comp[1], comp[2])) %>%
        as.data.frame() %>%
        tibble::rownames_to_column('gene') %>%
        mutate(cell_type = ct, comparison = comp_label) %>%
        arrange(padj),
      error = function(e) NULL
    )
    if (!is.null(res)) adipo_deseq_results[[paste0(ct, '_', comp_label)]] <- res
  }
}

adipo_all_deseq <- bind_rows(adipo_deseq_results)
write.csv(adipo_all_deseq, './output/Xenium/adipo_pseudobulk_deseq2.csv', row.names = FALSE)

adipo_sig_summary <- adipo_all_deseq %>%
  filter(padj < 0.05, abs(log2FoldChange) > 0.5) %>%
  group_by(cell_type, comparison) %>%
  summarise(n_up = sum(log2FoldChange > 0), n_down = sum(log2FoldChange < 0), .groups = 'drop')
print(adipo_sig_summary)

##############################################
#### Adipo subtypes: Volcano + GO/Reactome
##############################################

adipo_enrich_combined <- list()

for (ct in unique(adipo_all_deseq$cell_type)) {
  ct_slug <- gsub('_', '-', ct)
  cat('Adipo processing:', ct, '\n')

  dat_rvf     <- adipo_all_deseq %>% filter(cell_type == ct, comparison == 'RVF_vs_NF')  %>% cap_pval()
  dat_prv     <- adipo_all_deseq %>% filter(cell_type == ct, comparison == 'pRV_vs_NF')  %>% cap_pval()
  dat_rvf_prv <- adipo_all_deseq %>% filter(cell_type == ct, comparison == 'RVF_vs_pRV') %>% cap_pval()

  make_volcano_pdf(dat_rvf,     paste0('./output/Xenium/adipo_', ct_slug, '_volcano_RVF_vs_NF.pdf'),  paste(ct, ': RVF vs NF'))
  make_volcano_pdf(dat_prv,     paste0('./output/Xenium/adipo_', ct_slug, '_volcano_pRV_vs_NF.pdf'),  paste(ct, ': pRV vs NF'))
  make_volcano_pdf(dat_rvf_prv, paste0('./output/Xenium/adipo_', ct_slug, '_volcano_RVF_vs_pRV.pdf'), paste(ct, ': RVF vs pRV'))

  run_enrich <- function(dat, comp) {
    sig  <- dat %>% filter(padj < 0.05, abs(log2FoldChange) > 0.5)
    up   <- sig %>% filter(log2FoldChange > 0) %>% pull(gene)
    down <- sig %>% filter(log2FoldChange < 0) %>% pull(gene)
    list(
      up   = if (length(up)   > 0) enrichr(up,   enrichR_dbs) else NULL,
      down = if (length(down) > 0) enrichr(down, enrichR_dbs) else NULL
    )
  }

  er_rvf     <- run_enrich(dat_rvf,     'RVF_vs_NF')
  er_prv     <- run_enrich(dat_prv,     'pRV_vs_NF')
  er_rvf_prv <- run_enrich(dat_rvf_prv, 'RVF_vs_pRV')

  ct_enrich <- bind_rows(
    combine_enrich(er_rvf$up,       'up',   'RVF_vs_NF'),
    combine_enrich(er_rvf$down,     'down', 'RVF_vs_NF'),
    combine_enrich(er_prv$up,       'up',   'pRV_vs_NF'),
    combine_enrich(er_prv$down,     'down', 'pRV_vs_NF'),
    combine_enrich(er_rvf_prv$up,   'up',   'RVF_vs_pRV'),
    combine_enrich(er_rvf_prv$down, 'down', 'RVF_vs_pRV')
  )
  if (!is.null(ct_enrich) && nrow(ct_enrich) > 0) {
    ct_enrich$cell_type <- ct
    adipo_enrich_combined[[ct]] <- ct_enrich
  }

  pdf(paste0('./output/Xenium/adipo_', ct_slug, '_enrichr_plots.pdf'), width = 12, height = 7)
  for (comp_label in c('RVF_vs_NF', 'pRV_vs_NF', 'RVF_vs_pRV')) {
    er <- switch(comp_label, RVF_vs_NF = er_rvf, pRV_vs_NF = er_prv, RVF_vs_pRV = er_rvf_prv)
    for (dir_label in c('up', 'down')) {
      er_dir <- er[[dir_label]]
      if (!is.null(er_dir)) {
        p1 <- plot_enrich(er_dir$GO_Biological_Process_2023, paste(ct, comp_label, toupper(dir_label), 'GO BP'))
        p2 <- plot_enrich(er_dir$Reactome_2022,              paste(ct, comp_label, toupper(dir_label), 'Reactome'))
        if (!is.null(p1)) print(p1)
        if (!is.null(p2)) print(p2)
      }
    }
  }
  dev.off()
}

adipo_all_enrich <- bind_rows(adipo_enrich_combined)
write.csv(adipo_all_enrich, './output/Xenium/adipo_all_enrichr.csv', row.names = FALSE)


##############################################
#### Merge subclustering labels into xenium.obj
##############################################

# Load final clean objects and transfer subtype labels into xenium.obj
# Each RDS carries a *_subtype column with the final annotation


subcluster_rds <- list(
  list(file = '~/Documents/XeniumWorkflow/functions/output/Xenium/myeloid_clean_clean.rds', col = 'myeloid_subtype'),
  list(file = '~/Documents/XeniumWorkflow/functions/output/Xenium/nkt_clean_clean.rds',           col = 'nkt_subtype'),
  list(file = '~/Documents/XeniumWorkflow/functions/output/Xenium/ec_clean_clean.rds',      col = 'ec_subtype'),
  list(file = '~/Documents/XeniumWorkflow/functions/output/Xenium/fb_clean_clean.rds',      col = 'fb_subtype'),
  list(file = '~/Documents/XeniumWorkflow/functions/output/Xenium/cm_clean_clean.rds',      col = 'cm_subtype'),
  list(file = '~/Documents/XeniumWorkflow/functions/output/Xenium/neuron_clean_clean.rds',        col = 'neuron_subtype'),
  list(file = '~/Documents/XeniumWorkflow/functions/output/Xenium/mural_clean_clean.rds',         col = 'mural_subtype'),
  list(file = '~/Documents/XeniumWorkflow/functions/output/Xenium/epi_clean_clean.rds',     col = 'epi_subtype'),
  list(file = '~/Documents/XeniumWorkflow/functions/output/Xenium/adipo_clean_clean.rds',   col = 'adipo_subtype')
)

xenium.obj$cell_types_subclustering <- 'Unassigned'

for (entry in subcluster_rds) {
  if (!file.exists(entry$file)) {
    cat('Skipping (not found):', entry$file, '\n')
    next
  }
  obj <- readRDS(entry$file)
  sub_labels <- obj@meta.data[[entry$col]]
  names(sub_labels) <- colnames(obj)

  # Only overwrite cells that are in both xenium.obj and the subcluster object
  shared_cells <- intersect(names(sub_labels), colnames(xenium.obj))
  cat(entry$col, ':', length(shared_cells), 'cells transferred from', entry$file, '\n')
  xenium.obj$cell_types_subclustering[match(shared_cells, colnames(xenium.obj))] <- sub_labels[shared_cells]
}

cat('Unique subtypes:', length(unique(xenium.obj$cell_types_subclustering)), '\n')
table(xenium.obj$cell_types_subclustering) %>% sort(decreasing = TRUE) %>% print()

xenium.ref <- subset(xenium.obj,cell_types_subclustering != 'Unassigned')

saveRDS(xenium.ref,'~/Documents/XeniumWorkflow/functions/output/Xenium/xenium.ref.rds')

# Fill NA subclustering labels as Unassigned
xenium.obj$cell_types_subclustering[is.na(xenium.obj$cell_types_subclustering)] <- 'Unassigned'
xenium.ref$cell_types_subclustering[is.na(xenium.ref$cell_types_subclustering)] <- 'Unassigned'

DimPlot(xenium.obj,group.by='cell_types_subclustering',reduction='umap_orig',label=T)  + theme(legend.position = 'none')
DimPlot(xenium.ref,group.by='cell_types_subclustering',reduction='umap_orig',label=T)  + theme(legend.position = 'none')



pdf('~/Documents/XeniumWorkflow/functions//output/Xenium/Xenium_ref.pdf', width = 5, height = 5)
PlotEmbedding(xenium.obj, group.by = 'cell_types_subclustering', reduction = 'umap_orig',
              point_size = .1, plot_under = TRUE, plot_theme = umap_theme() + NoLegend(),
              raster_dpi = 400, raster_scale = 0.5)
dev.off()

pdf('~/Documents/XeniumWorkflow/functions//output/Xenium/Xenium_obj_w_ref.pdf', width = 5, height = 5)
PlotEmbedding(xenium.ref, group.by = 'cell_types_subclustering', reduction = 'umap_orig',
              point_size = .1, plot_under = TRUE, plot_theme = umap_theme() + NoLegend(),
              raster_dpi = 400, raster_scale = 0.5)
dev.off()

saveRDS(xenium.obj,'~/Documents/XeniumWorkflow/functions/output/Xenium/xenium_obj_subclustered.rds')


##############################################
#### Validation: recluster on panel genes only
#### (confirm subclustering was not an artifact of imputed genes)
##############################################

panel_genes_val <- panel_genes

# ── xenium.obj — full dataset, panel genes only ───────────────────────────────
xenium.obj.panel <- xenium.obj[panel_genes_val, ]
umap_cols <- grep("^umap_\\d+$|^UMAP_\\d+$|^umaporig_\\d+$",
                  colnames(xenium.obj.panel@meta.data), value = TRUE)
if (length(umap_cols)) xenium.obj.panel@meta.data[umap_cols] <- NULL
xenium.obj.panel[["pca"]]      <- NULL
xenium.obj.panel[["harmony"]]  <- NULL
xenium.obj.panel[["umap"]]     <- NULL
xenium.obj.panel[["umap_orig"]] <- NULL
xenium.obj.panel <- JoinLayers(xenium.obj.panel)
xenium.obj.panel <- NormalizeData(xenium.obj.panel)
xenium.obj.panel <- FindVariableFeatures(xenium.obj.panel, nfeatures = length(panel_genes_val))
xenium.obj.panel <- ScaleData(xenium.obj.panel, features = panel_genes_val)
xenium.obj.panel <- RunPCA(xenium.obj.panel, npcs = 30, features = panel_genes_val)
xenium.obj.panel <- RunHarmony(xenium.obj.panel, group.by.vars = 'patient',
                               reduction.save = 'harmony',
                               kmeans_init_nstart = 20, kmeans_init_iter_max = 100)
xenium.obj.panel <- FindNeighbors(xenium.obj.panel, dims = 1:30, reduction = 'harmony')
xenium.obj.panel <- FindClusters(xenium.obj.panel, resolution = 1)
xenium.obj.panel <- RunUMAP(xenium.obj.panel, reduction = 'harmony', dims = 1:30,
                             reduction.name = 'umap')

# Check 1 - do we get spatial separation consistent with original clustering
DimPlot(xenium.obj.panel,group.by='cell_types_subclustering',reduction='umap',label=T)  + theme(legend.position = 'none')

pdf('~/Documents/XeniumWorkflow/functions/output/Xenium/Xenium_obj_panel_only.pdf', width = 5, height = 5)
PlotEmbedding(xenium.obj.panel, group.by = 'cell_types_subclustering', reduction = 'umap',
              point_size = .1, plot_under = TRUE, plot_theme = umap_theme() + NoLegend(),
              raster_dpi = 400, raster_scale = 0.5)
dev.off()

# Check 2 - do manually annotated cells have superior clustering?
DimPlot(subset(xenium.obj.panel,cell_types_subclustering != 'Unassigned'),group.by='cell_types_subclustering',reduction='umap',label=T)  + theme(legend.position = 'none')

pdf('~/Documents/XeniumWorkflow/functions/output/Xenium/Xenium_obj_singlets_panel_only.pdf', width = 5, height = 5)
PlotEmbedding(subset(xenium.obj.panel,cell_types_subclustering != 'Unassigned'), group.by = 'cell_types_subclustering', reduction = 'umap',
              point_size = .1, plot_under = TRUE, plot_theme = umap_theme() + NoLegend(),
              raster_dpi = 400, raster_scale = 0.5)
dev.off()


# ── xenium.ref — subclustered cells only, panel genes ────────────────────────
xenium.ref.panel <- xenium.ref[panel_genes_val, ]
umap_cols <- grep("^umap_\\d+$|^UMAP_\\d+$|^umaporig_\\d+$",
                  colnames(xenium.ref.panel@meta.data), value = TRUE)
if (length(umap_cols)) xenium.ref.panel@meta.data[umap_cols] <- NULL
xenium.ref.panel[["pca"]]      <- NULL
xenium.ref.panel[["harmony"]]  <- NULL
xenium.ref.panel[["umap"]]     <- NULL
xenium.ref.panel[["umap_orig"]] <- NULL
xenium.ref.panel <- JoinLayers(xenium.ref.panel)
xenium.ref.panel <- NormalizeData(xenium.ref.panel)
xenium.ref.panel <- FindVariableFeatures(xenium.ref.panel, nfeatures = length(panel_genes_val))
xenium.ref.panel <- ScaleData(xenium.ref.panel, features = panel_genes_val)
xenium.ref.panel <- RunPCA(xenium.ref.panel, npcs = 30, features = panel_genes_val)
xenium.ref.panel <- RunHarmony(xenium.ref.panel, group.by.vars = 'patient',
                               reduction.save = 'harmony',
                               kmeans_init_nstart = 20, kmeans_init_iter_max = 100)
xenium.ref.panel <- FindNeighbors(xenium.ref.panel, dims = 1:30, reduction = 'harmony')
xenium.ref.panel <- FindClusters(xenium.ref.panel, resolution = 1)
xenium.ref.panel <- RunUMAP(xenium.ref.panel, reduction = 'harmony', dims = 1:30,
                             reduction.name = 'umap')





##############################################
#### RCTD: assign cell types to unassigned cells using xenium.ref.panel as reference
##############################################

library(spacexr)

# ── Build reference from xenium.ref.panel ─────────────────────────────────────
ref_counts <- GetAssayData(xenium.ref.panel, layer = 'counts')
ref_counts <- as(ref_counts, 'dgCMatrix')

ref_cell_types <- factor(xenium.ref.panel$cell_types_subclustering)
names(ref_cell_types) <- colnames(xenium.ref.panel)

# Reference() drops cells with < min_UMI (default 100) BEFORE checking the 25-cell
# minimum — so we must pre-filter to the same cells Reference() will keep, then
# drop any types that fall below 25 after that filter.
ref_nUMI    <- colSums(ref_counts)
keep_umi    <- ref_nUMI >= 10
ref_counts_f  <- ref_counts[, keep_umi]
ref_cell_types_f <- ref_cell_types[keep_umi]

type_counts_post_umi <- table(ref_cell_types_f)
cat('Types below 25 cells after UMI filter:\n')
print(type_counts_post_umi[type_counts_post_umi < 25])
valid_types <- names(type_counts_post_umi[type_counts_post_umi >= 25])
keep_ref    <- names(ref_cell_types_f)[ref_cell_types_f %in% valid_types]

ref_counts_f     <- ref_counts_f[, keep_ref]
ref_cell_types_f <- setNames(factor(as.character(ref_cell_types_f[keep_ref])), keep_ref)

cat('After filter:', ncol(ref_counts_f), 'cells,', nlevels(ref_cell_types_f), 'types\n')
cat('Min cells per type:', min(table(ref_cell_types_f)), '\n')
stopifnot(min(table(ref_cell_types_f)) >= 25)

# Pass min_UMI=0 — we already filtered above; don't let Reference() drop more cells
rctd_reference <- Reference(ref_counts_f, ref_cell_types_f, n_max_cells = 10000, min_UMI = 10)

# ── Build query from unassigned cells in xenium.obj.panel ─────────────────────
unassigned_cells <- colnames(xenium.obj.panel)[xenium.obj$cell_types_subclustering[colnames(xenium.obj.panel)] == 'Unassigned']
cat('Query: ', length(unassigned_cells), 'unassigned cells\n')

query_counts <- GetAssayData(xenium.obj.panel[, unassigned_cells], layer = 'counts')
query_counts <- as(query_counts, 'dgCMatrix')

query_coords <- xenium.obj@meta.data[unassigned_cells, c('x_centroid', 'y_centroid')]
colnames(query_coords) <- c('x', 'y')
rownames(query_coords) <- unassigned_cells

query_spatial <- SpatialRNA(query_coords, query_counts, colSums(query_counts))

# ── Run RCTD (doublet mode) ────────────────────────────────────────────────────
rctd_query <- create.RCTD(query_spatial, rctd_reference, max_cores = 16, test_mode = FALSE,UMI_min = 10)
rctd_query <- run.RCTD(rctd_query, doublet_mode = 'doublet')

saveRDS(rctd_query,'~/Documents/XeniumWorkflow/functions/output/Xenium/rctd_res_subclust.rds')

# ── Extract first-type labels and confidence ───────────────────────────────────
rctd_res_df <- rctd_query@results$results_df
rctd_cells <- rownames(rctd_res_df)
rctd_labels <- as.character(rctd_res_df$first_type)
names(rctd_labels) <- rctd_cells

cat('RCTD returned results for', length(rctd_cells), 'of', length(unassigned_cells), 'query cells\n')

# Spot types: singlet (high confidence), doublet_certain, doublet_uncertain, reject
rctd_spot_class <- as.character(rctd_res_df$spot_class)
names(rctd_spot_class) <- rctd_cells

cat('\nRCTD spot class distribution:\n')
print(table(rctd_spot_class))

# Assign labels to cells that were not rejected
assigned_mask <- rctd_spot_class != 'reject'
cat('Assigning', sum(assigned_mask), 'of', length(rctd_cells), 'RCTD-processed cells\n')

xenium.obj$cell_types_manual <- xenium.obj$cell_types_subclustering
xenium.obj$cell_types_manual[rctd_cells[assigned_mask]] <- rctd_labels[assigned_mask]
# Cells that RCTD rejects stay as Unassigned

# Propagate to xenium.obj.panel meta.data
xenium.obj.panel$cell_types_manual <- xenium.obj$cell_types_manual[colnames(xenium.obj.panel)]

cat('\nFinal cell_types_subclustering distribution:\n')
table(xenium.obj$cell_types_manual) %>% sort(decreasing = TRUE) %>% print()

# ── UMAP visualisation after RCTD fill-in ─────────────────────────────────────
DimPlot(subset(xenium.obj,cell_types_manual != 'Unassigned'), group.by = 'cell_types_manual', reduction = 'umap_orig',
        label = TRUE) + NoLegend() + ggtitle('Panel-only UMAP: subclustering + RCTD fill-in')

saveRDS(xenium.obj, '~/Documents/XeniumWorkflow/functions/output/Xenium/xenium_obj_subclustered.rds')


##############################################
#### Diagnose SOPA patch-boundary artifact
#### Patches are 1000µm non-overlapping tiles; cells at edges may lose transcripts
##############################################

# ── Step 1: Estimate grid origin ──────────────────────────────────────────────
# Brute-force search for the x/y offset that maximizes cells within 25µm of a
# boundary (i.e., finds where the real grid lines are).

all_x <- xenium.obj$x_centroid
all_y <- xenium.obj$y_centroid
patch_w <- 2000  # µm — must match the patch_width used in create_patches()

# Search offsets in 10µm steps
offsets <- seq(0, patch_w - 10, by = 10)
boundary_thr <- 25  # µm threshold for "near boundary"

best_ox <- 0; best_oy <- 0; best_n <- 0
for (ox in offsets) {
  dx <- pmin((all_x - ox) %% patch_w, patch_w - (all_x - ox) %% patch_w)
  n_x <- sum(dx < boundary_thr)
  if (n_x > best_n) { best_n <- n_x; best_ox <- ox }
}
for (oy in offsets) {
  dy <- pmin((all_y - oy) %% patch_w, patch_w - (all_y - oy) %% patch_w)
  n_y <- sum(dy < boundary_thr)
  if (n_y > best_n) { best_n <- n_y; best_oy <- oy }
}
cat('Estimated grid origin: x_offset =', best_ox, 'µm, y_offset =', best_oy, 'µm\n')

# ── Step 2: Distance to nearest patch boundary ───────────────────────────────
dist_x <- pmin((all_x - best_ox) %% patch_w, patch_w - (all_x - best_ox) %% patch_w)
dist_y <- pmin((all_y - best_oy) %% patch_w, patch_w - (all_y - best_oy) %% patch_w)
xenium.obj$dist_to_patch_boundary <- pmin(dist_x, dist_y)
xenium.obj$is_boundary <- xenium.obj$dist_to_patch_boundary < 50  # within 50µm

cat('Boundary cells (< 50µm):', sum(xenium.obj$is_boundary),
    '(', round(100 * mean(xenium.obj$is_boundary), 1), '%)\n')

# ── Step 3: Chi-square — CM subtype × boundary status ────────────────────────
cm_mask <- grepl('^CM_', xenium.obj$cell_types_subclustering)
cm_tab  <- table(
  subtype   = xenium.obj$cell_types_subclustering[cm_mask],
  boundary  = xenium.obj$is_boundary[cm_mask]
)
cat('\nCM subtype × boundary contingency table:\n')
print(cm_tab)
cm_chi <- chisq.test(cm_tab)
cat('Chi-square p =', format.pval(cm_chi$p.value, digits = 3), '\n')

# Per-subtype odds ratios (boundary enrichment)
cm_or <- apply(cm_tab, 1, function(row) {
  a <- row['TRUE']; b <- row['FALSE']
  total_b <- sum(cm_tab[, 'TRUE']); total_i <- sum(cm_tab[, 'FALSE'])
  or <- (a / total_b) / (b / total_i)
  return(or)
})
cat('\nCM subtype boundary enrichment (OR > 1 = enriched at boundary):\n')
print(sort(round(cm_or, 2), decreasing = TRUE))

# ── Repeat for ALL cell types (broad) ────────────────────────────────────────
broad_tab <- table(
  subtype  = xenium.obj$cell_types_subclustering,
  boundary = xenium.obj$is_boundary
)
broad_chi <- chisq.test(broad_tab)
cat('\nAll subtypes × boundary chi-square p =', format.pval(broad_chi$p.value, digits = 3), '\n')

# ── Step 4: Visualisation ────────────────────────────────────────────────────

# 4a: CM subtype proportions by distance bin
cm_df <- data.frame(
  subtype = xenium.obj$cell_types_subclustering[cm_mask],
  dist    = xenium.obj$dist_to_patch_boundary[cm_mask]
)
cm_df$dist_bin <- cut(cm_df$dist, breaks = c(0, 25, 50, 100, 200, 500),
                      include.lowest = TRUE, labels = c('0-25', '25-50', '50-100', '100-200', '200-500'))

cm_prop_by_dist <- cm_df %>%
  group_by(dist_bin, subtype) %>%
  summarise(n = n(), .groups = 'drop') %>%
  group_by(dist_bin) %>%
  mutate(prop = n / sum(n))

p_cm_boundary <- ggplot(cm_prop_by_dist, aes(x = dist_bin, y = prop, fill = subtype)) +
  geom_bar(stat = 'identity', position = 'stack') +
  theme_bw() +
  labs(x = 'Distance to nearest patch boundary (µm)', y = 'Proportion',
       title = 'CM subtype composition by distance to patch boundary') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.text = element_text(size = 7))

# 4b: Spatial plot — boundary cells highlighted
boundary_df <- data.frame(
  x = xenium.obj$x_centroid,
  y = -xenium.obj$y_centroid,
  boundary = xenium.obj$is_boundary
)

p_boundary_spatial <- ggplot(boundary_df, aes(x = x, y = y, color = boundary)) +
  geom_point(size = 0.01, alpha = 0.3) +
  scale_color_manual(values = c('FALSE' = 'grey80', 'TRUE' = 'red')) +
  coord_fixed() + theme_void() +
  labs(title = 'Patch boundary cells (red, < 50µm from edge)') +
  theme(legend.position = 'bottom')

pdf('./output/Xenium/patch_boundary_diagnostic.pdf', width = 12, height = 6)
print(p_cm_boundary)
print(p_boundary_spatial)
dev.off()

# ── Step 5: Save cell-level boundary data ────────────────────────────────────
boundary_out <- data.frame(
  cell = colnames(xenium.obj),
  dist_to_patch_boundary = xenium.obj$dist_to_patch_boundary,
  is_boundary = xenium.obj$is_boundary,
  cell_type = xenium.obj$cell_types_subclustering
)
write.csv(boundary_out, './output/Xenium/cell_patch_boundary_distances.csv', row.names = FALSE)
cat('Saved: patch_boundary_diagnostic.pdf, cell_patch_boundary_distances.csv\n')


'''
Bash

python seurat2cellnest.py Xenium_resegmented_imputed_final.rds --list_columns  


python seurat2cellnest.py Xenium_resegmented_imputed_final.rds --output_dir cellnest_input/ --fov_col patient --cell_type_col cell_type_rctd_doublet


'''



### Niches

M1_list <- SplitObject(subset(xenium.obj,cell_types_manual != 'Unassigned'), split.by = "patient")
fov_list <- c('X1467','X1632','X1561','X1343','X1697','X1691','X1618','X1567','X1692')

# Recalculate just for broad celltypes

# Calculate neighbors

for ( i in seq_along(M1_list) ){ message(i)
	object = M1_list[[i]]
	fov = fov_list[[i]]
	object <- BuildNicheAssay(object = object, fov = fov, group.by = "cell_types_manual",niches.k = 10, neighbors.k = 100)
	M1_list[[i]] = object
}


# Concatenate neighbors across all samples

col.unique = lapply( seq_along(M1_list), function(i){
    colnames(t( M1_list[[i]][['niche']]@counts ))
}  )
shared_feats <- Reduce(intersect,col.unique)


DAT.counts = lapply( seq_along(M1_list), function(i){
    t( M1_list[[i]][['niche']]@counts )[,shared_feats]
}  )
DAT.counts <- do.call("rbind", DAT.counts)

DAT.data = lapply( seq_along(M1_list), function(i){
    t( M1_list[[i]][['niche']]@data )[,shared_feats]
}  )
DAT.data <- do.call("rbind", DAT.data)

DAT.scale.data = lapply( seq_along(M1_list), function(i){
    t( M1_list[[i]][['niche']]@scale.data )[,shared_feats]
}  )
DAT.scale.data <- do.call("rbind", DAT.scale.data)

# Pad niche matrices with zeros for Unassigned cells so dimensions match xenium.obj
all_cells <- colnames(xenium.obj)
niche_cells <- rownames(DAT.counts)
missing_cells <- setdiff(all_cells, niche_cells)

if (length(missing_cells) > 0) {
  zero_pad <- matrix(0, nrow = length(missing_cells), ncol = ncol(DAT.counts),
                     dimnames = list(missing_cells, colnames(DAT.counts)))
  DAT.counts.full <- rbind(DAT.counts, zero_pad)[all_cells, ]
  DAT.data.full   <- rbind(DAT.data, zero_pad)[all_cells, ]
  DAT.scale.full  <- rbind(DAT.scale.data, zero_pad)[all_cells, ]
} else {
  DAT.counts.full <- DAT.counts[all_cells, ]
  DAT.data.full   <- DAT.data[all_cells, ]
  DAT.scale.full  <- DAT.scale.data[all_cells, ]
}

niche.assay <- CreateAssayObject(counts = t(DAT.counts.full))
niche.assay@data <- t(DAT.data.full)
niche.assay@scale.data <- t(DAT.scale.full)

xenium.obj[['niche_broad']] <- niche.assay
DefaultAssay(xenium.obj) <- 'niche_broad'

niches.k.range = c(10,15,20)

# CLR transform niche proportions to amplify differences in rare cell type fractions
# (raw proportions are dominated by the majority cell type, collapsing kmeans clusters)
DAT.clr <- DAT.data + 1e-6  # small pseudocount to avoid log(0)
DAT.clr <- log(DAT.clr) - rowMeans(log(DAT.clr))
DAT.clr <- scale(DAT.clr)  # z-score after CLR for balanced kmeans distances

res.clusters = data.frame(rownames = rownames(DAT.clr))

for ( k in niches.k.range ){ message("k=", k)
    # new column name
    newCol = paste0("kmeans_", k)
    # get centroids
    km_mb = ClusterR::MiniBatchKmeans(
        "data" = DAT.clr
        , "clusters" = k # the number of clusters
        , "batch_size" = 20 # the size of the mini batches
        , "num_init" = 20 # number of times the algorithm will be run with different centroid seeds
        , "max_iters" = 100 # the maximum number of clustering iterations. 
        , "init_fraction" = 0.2 # percentage of data to use for the initialization centroids (applies if initializer is kmeans++ or optimal_init). Should be a float number between 0.0 and 1.0.
        , "initializer" = "kmeans++" # the method of initialization. One of, optimal_init, quantile_init, kmeans++ and random. See details for more information
        , "early_stop_iter" = 10 # continue that many iterations after calculation of the best within-cluster-sum-of-squared-error
        , "verbose" = F
        , "CENTROIDS" = NULL
        , "tol" = 1e-04
        , "tol_optimal_init" = 0.3
        , "seed" = 1
    )
    
    # use centroids to get clusters

    res.clusters[,newCol] = ClusterR::predict_MBatchKMeans( # This function takes the data and the output centroids and returns the clusters.
        "data" = DAT.clr
        , "CENTROIDS" = km_mb$centroids
    )
    res.clusters[,newCol] = as.factor( res.clusters[,newCol] ) # change clusters to factors
    
}

rownames(res.clusters) <- res.clusters$rownames
res.clusters <- res.clusters[,-1,drop=F]

write.table(res.clusters,'./output/Xenium/Niche_bulk_clusters.csv',sep=',')

# Assign niche clusters back to xenium.obj; Unassigned cells get NA
niche_cells <- rownames(res.clusters)
for(i in colnames(res.clusters))
  {
    xenium.obj[[i]] <- NA
    xenium.obj[[i]][niche_cells,] <- res.clusters[niche_cells, i]
  }

niche.patient <- table(xenium.obj$kmeans_15,xenium.obj$patient)
niche.patient <- t(t(niche.patient) / colSums(niche.patient))



disease <- rep(c('NF','pRV','RVF'),9)[t(table(xenium.obj$patient,xenium.obj$group) > 0)]
disease <- c(t(replicate(15,disease)))

niche.patient <- data.frame(niche = disease,niche.patient)

niche.patient$niche <- factor(niche.patient$niche, levels=c('NF','pRV','RVF'))

ggplot(niche.patient,aes(Var1,Freq,color = niche))+geom_boxplot() + theme_classic()



niche.names <- table(xenium.obj$kmeans_15,xenium.obj$cell_type_rctd_doublet)
niche.names <- niche.names / rowSums(niche.names)

niche.names_alt <- table(xenium.obj$kmeans_15,xenium.obj$cell_type_seurat)
niche.names_alt <- niche.names_alt / rowSums(niche.names_alt)


##############################################
#### Niche annotation
##############################################

# --- 1. Heatmap of niche composition (row-normalized) ---
niche_df <- melt(niche.names)
colnames(niche_df) <- c("niche", "cell_type", "fraction")

pdf('./output/Xenium/niche_composition_heatmap.pdf', width = 12, height = 6)
ggplot(niche_df, aes(x = cell_type, y = factor(niche), fill = fraction)) +
  geom_tile() +
  scale_fill_gradient2(low = "white", high = "darkred", midpoint = 0.15,
                       name = "Fraction") +
  geom_text(aes(label = ifelse(fraction > 0.05, sprintf("%.0f%%", fraction * 100), "")),
            size = 2.5) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        axis.text.y = element_text(size = 10)) +
  labs(x = NULL, y = "Niche cluster", title = "Cell type composition per niche (k=15)")
dev.off()

# --- 2. Print top 3 cell types per niche ---
for (k in 1:nrow(niche.names)) {
  top3 <- sort(niche.names[k, ], decreasing = TRUE)[1:3]
  cat(sprintf("Niche %2d: %s (%.0f%%), %s (%.0f%%), %s (%.0f%%)\n",
              k,
              names(top3)[1], top3[1]*100,
              names(top3)[2], top3[2]*100,
              names(top3)[3], top3[3]*100))
}

# --- 3. Assign niche labels ---
niche_labels_15 <- c(
  '1'  = 'Adipose-myocardial',               # CM_Contractile 11%, Adipo_Quiescent 8%, CM_RV 7%, CM_Ventricular_4 6%
  '2'  = 'Mixed myocardium',                   # CM_RV 13%, CM_Ventricular 11%, CM_Contractile 9%, Capillary_EC 8%
  '3'  = 'Epicardial-fibrotic stroma',        # FB_Homeostatic 34%, Macrophage_Inflammatory 8%, Epi_Activated 6%
  '4'  = 'Interstitial myocardium',           # CM_Contractile 14%, CM_RV 8%, CM_RV_3 7%, FB_NTN4 6%, Arterial_EC 6%
  '5'  = 'Sub-epicardial myocardium',         # CM_RV 13%, CM_Ventricular_4 9%, CM_Contractile 9%, Epi_Quiescent 7%
  '6'  = 'Peri-arterial myocardium',          # CM_RV 12%, CM_Contractile 11%, CM_Ventricular_4 8%, Arterial_EC 8%, Capillary_EC 6%
  '7'  = 'Immune-infiltrated stroma',         # FB_Homeostatic 12%, Monocyte 8%, Macrophage_Inflammatory 7%, Venous_EC 7%, Macrophage_C1q 6%
  '8'  = 'Remodeling stroma',                 # FB_Resting 10%, FB_PCOLCE2 9%, Capillary_EC 7%, Macrophage_C1q 7%
  '9'  = 'Perivascular',                      # Arterial_EC 9%, Capillary_EC 8%, FB_Resting 7%, VSMC_Arterial 6%, CM_Ventricular 6%
  '10' = 'Capillary-rich myocardium',         # Capillary_EC 10%, CM_Contractile 8%, CM_RV 8%, CM_Ventricular_4 6%
  '11' = 'Compact myocardium',                # CM_RV 14%, CM_Contractile 10%, CM_Ventricular_4 8%, CM_Ventricular 7%, CM_Ventricular_2 6%
  '12' = 'Arterial-adjacent myocardium',      # Arterial_EC 11%, CM_Contractile 7%, CM_RV 7%, CM_Ventricular_4 7%, Capillary_EC 6%
  '13' = 'Fibrotic myocardium',               # FB_Resting 9%, CM_Ventricular 9%, Capillary_EC 8%, FB_PCOLCE2 8%
  '14' = 'Arterial wall',                     # VSMC_Arterial 10%, VSMC 7%, FB_PCOLCE2 7%, Macrophage_C1q 7%, Venous_EC 6%
  '15' = 'Venous-perivascular'                # Venous_EC 9%, Macrophage_C1q 6%, Capillary_EC 6%, FB_PCOLCE2 6%, Adipo_Mature 6%
)


xenium.obj$niche_label <- unname(niche_labels_15[as.character(xenium.obj$kmeans_15)])

Idents(xenium.obj) <- "kmeans_15"

# --- 4. Subtype-level composition per niche ---
niche.names_alt <- table(xenium.obj$kmeans_15, xenium.obj$cell_type_seurat)
niche.names_alt <- niche.names_alt / rowSums(niche.names_alt)

for (k in 1:nrow(niche.names_alt)) {
  top5 <- sort(niche.names_alt[k, ], decreasing = TRUE)[1:5]
  cat(sprintf("Niche %2d: %s (%.0f%%), %s (%.0f%%), %s (%.0f%%), %s (%.0f%%), %s (%.0f%%)\n",
              k,
              names(top5)[1], top5[1]*100,
              names(top5)[2], top5[2]*100,
              names(top5)[3], top5[3]*100,
              names(top5)[4], top5[4]*100,
              names(top5)[5], top5[5]*100))
}

# Subtype-level heatmap
niche_alt_df <- melt(niche.names_alt)
colnames(niche_alt_df) <- c("niche", "cell_subtype", "fraction")

pdf('./output/Xenium/niche_composition_subtype_heatmap.pdf', width = 18, height = 6)
ggplot(niche_alt_df, aes(x = cell_subtype, y = factor(niche), fill = fraction)) +
  geom_tile() +
  scale_fill_gradient2(low = "white", high = "darkred", midpoint = 0.08,
                       name = "Fraction") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.text.y = element_text(size = 10)) +
  labs(x = NULL, y = "Niche cluster",
       title = "Cell subtype composition per niche (k=15)")
dev.off()

##############################################
#### Niche proportions by disease state
##############################################

library(ggpubr)
library(tidyr)

niche_meta <- xenium.obj@meta.data %>%
  filter(!is.na(niche_label)) %>%
  select(patient, group, niche_label)

niche_totals <- niche_meta %>%
  group_by(patient, group) %>%
  summarise(total = n(), .groups = 'drop')

niche_prop <- niche_meta %>%
  group_by(patient, group, niche_label) %>%
  summarise(n = n(), .groups = 'drop') %>%
  left_join(niche_totals, by = c('patient', 'group')) %>%
  mutate(proportion = n / total) %>%
  complete(nesting(patient, group), niche_label, fill = list(n = 0, proportion = 0)) %>%
  mutate(group = factor(group, levels = c('NF', 'pRV', 'RVF')))

pairwise_comparisons <- list(c('NF', 'pRV'), c('NF', 'RVF'), c('pRV', 'RVF'))

p_niche_prop <- ggplot(niche_prop, aes(x = group, y = proportion, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.5) +
  stat_compare_means(comparisons = pairwise_comparisons, method = 'wilcox.test', label = 'p.signif') +
  stat_compare_means(method = 'kruskal.test', label.y.npc = 0.95, size = 3) +
  facet_wrap(~ niche_label, scales = 'free_y', ncol = 5) +
  scale_x_discrete(limits = c('NF', 'pRV', 'RVF')) +
  theme_bw() +
  theme(legend.position = 'none',
        strip.text = element_text(face = 'bold', size = 7),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = 'Proportion of cells',
       title = 'Spatial niche distribution by disease state (k=15)')

pdf('./output/Xenium/niche_proportions_by_group.pdf', width = 14, height = 10)
print(p_niche_prop)
dev.off()

# All niches on a single plot, colored by disease group, sorted by overall frequency
niche_order <- niche_prop %>%
  group_by(niche_label) %>%
  summarise(mean_prop = mean(proportion), .groups = 'drop') %>%
  arrange(mean_prop) %>%
  pull(niche_label)
niche_prop$niche_label <- factor(niche_prop$niche_label, levels = niche_order)

# Kruskal-Wallis per niche (BH-adjusted), then pairwise Wilcoxon for significant niches
niche_stats <- niche_prop %>%
  group_by(niche_label) %>%
  summarise(kw_p = kruskal.test(proportion ~ group)$p.value, .groups = 'drop') %>%
  mutate(kw_padj = p.adjust(kw_p, method = 'BH'),
         sig = ifelse(kw_padj < 1, sprintf('p=%.3f', kw_padj), 'ns'))

cat('\nNiche Kruskal-Wallis results (BH-adjusted):\n')
print(as.data.frame(niche_stats %>% arrange(kw_padj)), row.names = FALSE)

niche_prop$group <- factor(niche_prop$group, levels = c('RVF', 'pRV', 'NF'))
niche_stats$niche_label <- factor(niche_stats$niche_label, levels = niche_order)

p_niche_combined <- ggplot(niche_prop, aes(x = niche_label, y = proportion, color = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15), size = 1.5) +
  scale_color_manual(values = c(NF = 'steelblue', pRV = 'darkorange', RVF = 'firebrick')) +
  coord_flip() +
  theme_classic() +
  theme(axis.text.y = element_text(size = 8)) +
  geom_text(data = niche_stats, aes(x = niche_label, y = Inf, label = sig),
            inherit.aes = FALSE, hjust = 1.1, size = 2.5, color = 'black') +
  labs(x = NULL, y = 'Proportion of cells', color = 'Group',
       title = 'Spatial niche distribution by disease state (k=15)')

pdf('./output/Xenium/niche_proportions_by_group_combined.pdf', width = 10, height = 7)
print(p_niche_combined)
dev.off()

# ── Spatial niche map ─────────────────────────────────────────────────────────

# Tiling helper: arranges multiple patients per group into compact layouts
# (2 patients = side-by-side; 3 = 2 bottom + 1 centered on top; 4+ = 2-col grid)
offset_coords <- function(meta, gap_um = 500, row_height_um = 12000) {
  for (g in unique(meta$group)) {
    idx_g <- which(meta$group == g)
    patients <- unique(meta$patient[idx_g])
    for (p in patients) {
      idx_p <- which(meta$group == g & meta$patient == p)
      meta$x[idx_p] <- meta$x[idx_p] - min(meta$x[idx_p])
      meta$y[idx_p] <- meta$y[idx_p] - min(meta$y[idx_p])
    }
    if (length(patients) <= 1) next
    bboxes <- lapply(patients, function(p) {
      idx <- which(meta$group == g & meta$patient == p)
      list(w = diff(range(meta$x[idx])), h = diff(range(meta$y[idx])))
    })
    names(bboxes) <- patients
    areas <- sapply(bboxes, function(b) b$w * b$h)
    patients <- patients[order(-areas)]
    bboxes <- bboxes[as.character(patients)]
    if (length(patients) == 2) {
      idx1 <- which(meta$group == g & meta$patient == patients[1])
      idx2 <- which(meta$group == g & meta$patient == patients[2])
      meta$x[idx2] <- meta$x[idx2] + bboxes[[1]]$w + gap_um
    } else if (length(patients) == 3) {
      bot1 <- patients[1]; bot2 <- patients[2]; top1 <- patients[3]
      idx_bot1 <- which(meta$group == g & meta$patient == bot1)
      idx_bot2 <- which(meta$group == g & meta$patient == bot2)
      idx_top  <- which(meta$group == g & meta$patient == top1)
      bot_total_w <- bboxes[[1]]$w + gap_um + bboxes[[2]]$w
      meta$x[idx_bot2] <- meta$x[idx_bot2] + bboxes[[1]]$w + gap_um
      top_x_offset <- (bot_total_w - bboxes[[3]]$w) / 2
      meta$x[idx_top] <- meta$x[idx_top] + top_x_offset
      meta$y[idx_top] <- meta$y[idx_top] + row_height_um
    } else {
      ncols <- 2
      nrows <- ceiling(length(patients) / ncols)
      col_widths <- rep(0, ncols)
      row_heights_v <- rep(0, nrows)
      for (i in seq_along(patients)) {
        r <- ceiling(i / ncols)
        cc <- ((i - 1) %% ncols) + 1
        col_widths[cc] <- max(col_widths[cc], bboxes[[i]]$w)
        row_heights_v[r] <- max(row_heights_v[r], bboxes[[i]]$h)
      }
      x_starts <- c(0, cumsum(col_widths[-ncols] + gap_um))
      y_starts <- c(0, cumsum(row_heights_v[-nrows] + gap_um))
      for (i in seq_along(patients)) {
        p <- patients[i]
        r <- ceiling(i / ncols)
        cc <- ((i - 1) %% ncols) + 1
        idx_p <- which(meta$group == g & meta$patient == p)
        meta$x[idx_p] <- meta$x[idx_p] + x_starts[cc]
        meta$y[idx_p] <- meta$y[idx_p] + y_starts[r]
      }
    }
  }
  meta
}

# Build spatial metadata — only cells with niche assignments
niche_cells_mask <- !is.na(xenium.obj$niche_label)
meta_niche <- data.frame(
  x = xenium.obj$x_centroid[niche_cells_mask],
  y = -xenium.obj$y_centroid[niche_cells_mask],
  group = xenium.obj$group[niche_cells_mask],
  patient = xenium.obj$patient[niche_cells_mask],
  niche = xenium.obj$niche_label[niche_cells_mask]
)
meta_niche$group <- factor(meta_niche$group, levels = c('NF', 'pRV', 'RVF'))
meta_niche <- offset_coords(meta_niche)

p_niche_spatial <- ggplot(meta_niche, aes(x = x, y = y, color = niche)) +
  ggrastr::rasterise(geom_point(size = 0.1, stroke = 0), dpi = 300) +
  facet_wrap(~ group, nrow = 1) +
  labs(title = 'Spatial niche map (k=15)', color = 'Niche') +
  theme_void(base_size = 14) +
  theme(strip.text = element_text(face = 'bold', size = 16),
        legend.position = 'right',
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 12, face = 'bold')) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  coord_fixed()

pdf('./output/Xenium/niche_spatial_map.pdf', width = 20, height = 8)
print(p_niche_spatial)
dev.off()
cat('Saved niche_spatial_map.pdf\n')

# ── Niche composition by disease group for key niches ─────────────────────────

focus_niches <- c('Capillary-rich myocardium', 'Mixed myocardium')

# Per-patient cell subtype proportions within each focus niche
niche_comp <- data.frame(
  patient = xenium.obj$patient[niche_cells_mask],
  group = xenium.obj$group[niche_cells_mask],
  niche = xenium.obj$niche_label[niche_cells_mask],
  cell_type = xenium.obj$cell_types_manual[niche_cells_mask]
)
niche_comp <- niche_comp[niche_comp$niche %in% focus_niches, ]
niche_comp$group <- factor(niche_comp$group, levels = c('NF', 'pRV', 'RVF'))

# Compute per-patient proportions
niche_comp_prop <- niche_comp %>%
  group_by(niche, patient, group, cell_type) %>%
  summarise(n = n(), .groups = 'drop') %>%
  group_by(niche, patient) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# Keep top cell types (>3% in any patient×niche combo) to avoid clutter
top_types <- niche_comp_prop %>%
  group_by(niche, cell_type) %>%
  summarise(max_prop = max(prop), .groups = 'drop') %>%
  filter(max_prop > 0.03) %>%
  pull(cell_type) %>%
  unique()

niche_comp_prop_filt <- niche_comp_prop %>%
  filter(cell_type %in% top_types)

# (A) Grouped bar chart: mean proportion per group, faceted by niche
niche_comp_mean <- niche_comp_prop_filt %>%
  group_by(niche, group, cell_type) %>%
  summarise(mean_prop = mean(prop), se = sd(prop) / sqrt(n()), .groups = 'drop')

p_niche_comp_bar <- ggplot(niche_comp_mean, aes(x = reorder(cell_type, -mean_prop), y = mean_prop, fill = group)) +
  geom_bar(stat = 'identity', position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = mean_prop - se, ymax = mean_prop + se),
                position = position_dodge(width = 0.8), width = 0.25, linewidth = 0.3) +
  scale_fill_manual(values = c(NF = 'steelblue', pRV = 'darkorange', RVF = 'firebrick')) +
  facet_wrap(~ niche, scales = 'free_x') +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        strip.text = element_text(face = 'bold', size = 12)) +
  labs(x = NULL, y = 'Proportion of cells within niche', fill = 'Group',
       title = 'Cell subtype composition of key niches by disease group')

pdf('./output/Xenium/niche_composition_focus.pdf', width = 16, height = 6)
print(p_niche_comp_bar)
dev.off()
cat('Saved niche_composition_focus.pdf\n')

# (B) Dot plot: per-patient proportions with group coloring
p_niche_comp_dot <- ggplot(niche_comp_prop_filt, aes(x = reorder(cell_type, -prop), y = prop, color = group)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6),
              size = 2, alpha = 0.8) +
  stat_summary(aes(group = group), fun = mean, geom = 'crossbar',
               position = position_dodge(width = 0.6), width = 0.4, linewidth = 0.4) +
  scale_color_manual(values = c(NF = 'steelblue', pRV = 'darkorange', RVF = 'firebrick')) +
  facet_wrap(~ niche, scales = 'free_x') +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        strip.text = element_text(face = 'bold', size = 12)) +
  labs(x = NULL, y = 'Proportion of cells within niche', color = 'Group',
       title = 'Per-patient cell subtype composition of key niches')

pdf('./output/Xenium/niche_composition_focus_dotplot.pdf', width = 16, height = 6)
print(p_niche_comp_dot)
dev.off()
cat('Saved niche_composition_focus_dotplot.pdf\n')

# ── Niche analysis with collapsed CM subtypes ─────────────────────────────────
# CM subtypes show patient-specific variation that doesn't map to disease state,
# causing niches to split on CM subtype ratios rather than microenvironment.
# Collapse all CM_* subtypes into a single "CM" label and re-run niche analysis.

xenium.obj$cell_types_cm_collapsed <- as.character(xenium.obj$cell_types_manual)
cm_mask <- grepl('^CM_', xenium.obj$cell_types_cm_collapsed)
xenium.obj$cell_types_cm_collapsed[cm_mask] <- 'CM'
cat('Collapsed', sum(cm_mask), 'CM subtype cells into single CM label\n')
cat('Unique cell types after collapse:', length(unique(xenium.obj$cell_types_cm_collapsed)), '\n')

M1_list_cm <- SplitObject(subset(xenium.obj, cell_types_manual != 'Unassigned'), split.by = 'patient')
fov_list <- c('X1467','X1632','X1561','X1343','X1697','X1691','X1618','X1567','X1692')

for (i in seq_along(M1_list_cm)) {
  message('BuildNicheAssay (collapsed CM) - patient ', i)
  object <- M1_list_cm[[i]]
  fov <- fov_list[[i]]
  object <- BuildNicheAssay(object = object, fov = fov, group.by = 'cell_types_cm_collapsed',
                            niches.k = 10, neighbors.k = 100)
  M1_list_cm[[i]] <- object
}

# Concatenate niche matrices
col.unique.cm <- lapply(seq_along(M1_list_cm), function(i) {
  colnames(t(M1_list_cm[[i]][['niche']]@counts))
})
shared_feats_cm <- Reduce(intersect, col.unique.cm)

DAT.data.cm <- do.call('rbind', lapply(seq_along(M1_list_cm), function(i) {
  t(M1_list_cm[[i]][['niche']]@data)[, shared_feats_cm]
}))

# CLR transform + scale
DAT.clr.cm <- DAT.data.cm + 1e-6
DAT.clr.cm <- log(DAT.clr.cm) - rowMeans(log(DAT.clr.cm))
DAT.clr.cm <- scale(DAT.clr.cm)

# MiniBatchKmeans with k=15
km_cm <- ClusterR::MiniBatchKmeans(
  data = DAT.clr.cm, clusters = 15, batch_size = 20, num_init = 20,
  max_iters = 100, init_fraction = 0.2, initializer = 'kmeans++',
  early_stop_iter = 10, verbose = FALSE, CENTROIDS = NULL,
  tol = 1e-04, tol_optimal_init = 0.3, seed = 1
)
km_cm_clusters <- ClusterR::predict_MBatchKMeans(data = DAT.clr.cm, CENTROIDS = km_cm$centroids)

# Assign back to xenium.obj
xenium.obj$kmeans_15_cm <- NA
xenium.obj$kmeans_15_cm[rownames(DAT.clr.cm)] <- as.factor(km_cm_clusters)

# Print niche compositions
niche_comp_cm <- table(xenium.obj$kmeans_15_cm, xenium.obj$cell_types_cm_collapsed)
niche_comp_cm <- niche_comp_cm / rowSums(niche_comp_cm)

cat('\n=== Niche composition (collapsed CM, k=15) ===\n')
for (k in 1:nrow(niche_comp_cm)) {
  top5 <- sort(niche_comp_cm[k, ], decreasing = TRUE)[1:5]
  cat(sprintf('Niche %2d: %s (%.0f%%), %s (%.0f%%), %s (%.0f%%), %s (%.0f%%), %s (%.0f%%)\n',
              k, names(top5)[1], top5[1]*100, names(top5)[2], top5[2]*100,
              names(top5)[3], top5[3]*100, names(top5)[4], top5[4]*100,
              names(top5)[5], top5[5]*100))
}

# Assign niche labels (collapsed CM, k=15)
niche_labels_cm15 <- c(
  '1'  = 'Epicardial-fibrotic stroma',
  '2'  = 'Arterial wall',
  '3'  = 'Peri-arterial myocardium',        # grouped with 10
  '4'  = 'Venous-immune stroma',
  '5'  = 'Adipose-myocardial',
  '6'  = 'Sub-epicardial myocardium',
  '7'  = 'Remodeling stroma',
  '8'  = 'Homeostatic-myocardial',
  '9'  = 'Fibrotic myocardium',
  '10' = 'Peri-arterial myocardium',         # grouped with 3
  '11' = 'Adipose stroma',                   # grouped with 14, 15
  '12' = 'Capillary-rich myocardium',
  '13' = 'Immune-surveilled myocardium',
  '14' = 'Adipose stroma',                   # grouped with 11, 15
  '15' = 'Adipose stroma'                    # grouped with 11, 14
)

xenium.obj$niche_label_cm <- unname(niche_labels_cm15[as.character(xenium.obj$kmeans_15_cm)])
cat('Collapsed-CM niche labels (12 unique):\n')
print(sort(table(xenium.obj$niche_label_cm), decreasing = TRUE))

# Heatmap (using original 15 clusters before grouping)
niche_cm_df <- reshape2::melt(niche_comp_cm)
colnames(niche_cm_df) <- c('niche', 'cell_type', 'fraction')

pdf('./output/Xenium/niche_composition_collapsed_CM_heatmap.pdf', width = 14, height = 6)
ggplot(niche_cm_df, aes(x = cell_type, y = factor(niche), fill = fraction)) +
  geom_tile() +
  scale_fill_gradient2(low = 'white', high = 'darkred', midpoint = 0.1, name = 'Fraction') +
  geom_text(aes(label = ifelse(fraction > 0.05, sprintf('%.0f%%', fraction * 100), '')),
            size = 2.5) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        axis.text.y = element_text(size = 10)) +
  labs(x = NULL, y = 'Niche cluster',
       title = 'Niche composition (collapsed CM subtypes, k=15)')
dev.off()
cat('Saved niche_composition_collapsed_CM_heatmap.pdf\n')

# Niche proportions by disease state (using grouped labels = 12 niches)
niche_cells_cm_mask <- !is.na(xenium.obj$niche_label_cm)
niche_meta_cm <- data.frame(
  patient = xenium.obj$patient[niche_cells_cm_mask],
  group = xenium.obj$group[niche_cells_cm_mask],
  niche = xenium.obj$niche_label_cm[niche_cells_cm_mask]
)
niche_meta_cm$group <- factor(niche_meta_cm$group, levels = c('NF', 'pRV', 'RVF'))

niche_totals_cm <- niche_meta_cm %>%
  group_by(patient, group) %>%
  summarise(total = n(), .groups = 'drop')

niche_prop_cm <- niche_meta_cm %>%
  group_by(patient, group, niche) %>%
  summarise(n = n(), .groups = 'drop') %>%
  left_join(niche_totals_cm, by = c('patient', 'group')) %>%
  mutate(proportion = n / total) %>%
  complete(nesting(patient, group), niche, fill = list(n = 0, proportion = 0)) %>%
  mutate(group = factor(group, levels = c('NF', 'pRV', 'RVF')))

# KW stats
niche_stats_cm <- niche_prop_cm %>%
  group_by(niche) %>%
  summarise(kw_p = kruskal.test(proportion ~ group)$p.value, .groups = 'drop') %>%
  mutate(kw_padj = p.adjust(kw_p, method = 'BH'),
         sig = ifelse(kw_padj < 1, sprintf('p=%.3f', kw_padj), 'ns'))

cat('\nCollapsed-CM niche Kruskal-Wallis results (BH-adjusted):\n')
print(as.data.frame(niche_stats_cm %>% arrange(kw_padj)), row.names = FALSE)

# Faceted boxplots
pairwise_comparisons <- list(c('NF', 'pRV'), c('NF', 'RVF'), c('pRV', 'RVF'))

p_niche_cm <- ggplot(niche_prop_cm, aes(x = group, y = proportion, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, size = 1.5) +
  stat_compare_means(comparisons = pairwise_comparisons, method = 'wilcox.test',
                     label = 'p.signif', size = 3) +
  stat_compare_means(method = 'kruskal.test', label.y.npc = 0.95, size = 2.5) +
  scale_fill_manual(values = c(NF = 'steelblue', pRV = 'darkorange', RVF = 'firebrick')) +
  facet_wrap(~ niche, scales = 'free_y', ncol = 5) +
  theme_bw(base_size = 10) +
  theme(legend.position = 'none',
        strip.text = element_text(face = 'bold', size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = 'Proportion of cells',
       title = 'Niche proportions by disease state (collapsed CM, k=15)')

pdf('./output/Xenium/niche_proportions_collapsed_CM.pdf', width = 14, height = 10)
print(p_niche_cm)
dev.off()
cat('Saved niche_proportions_collapsed_CM.pdf\n')

# Spatial map
meta_niche_cm <- data.frame(
  x = xenium.obj$x_centroid[niche_cells_cm_mask],
  y = -xenium.obj$y_centroid[niche_cells_cm_mask],
  group = xenium.obj$group[niche_cells_cm_mask],
  patient = xenium.obj$patient[niche_cells_cm_mask],
  niche = xenium.obj$niche_label_cm[niche_cells_cm_mask]
)
meta_niche_cm$group <- factor(meta_niche_cm$group, levels = c('NF', 'pRV', 'RVF'))
meta_niche_cm <- offset_coords(meta_niche_cm)

p_niche_cm_spatial <- ggplot(meta_niche_cm, aes(x = x, y = y, color = niche)) +
  ggrastr::rasterise(geom_point(size = 0.1, stroke = 0), dpi = 300) +
  facet_wrap(~ group, nrow = 1) +
  labs(title = 'Spatial niche map (collapsed CM, 12 niches)', color = 'Niche') +
  theme_void(base_size = 14) +
  theme(strip.text = element_text(face = 'bold', size = 16),
        legend.position = 'right',
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 12, face = 'bold')) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  coord_fixed()

pdf('./output/Xenium/niche_spatial_map_collapsed_CM.pdf', width = 20, height = 8)
print(p_niche_cm_spatial)
dev.off()
cat('Saved niche_spatial_map_collapsed_CM.pdf\n')


# ── Cell type frequency by disease state / patient ────────────────────────────

# Compute per-patient cell type proportions (excluding Unassigned)
ct_freq <- data.frame(
  patient = xenium.obj$patient,
  group = xenium.obj$group,
  cell_type = xenium.obj$cell_types_manual
)
ct_freq <- ct_freq[ct_freq$cell_type != 'Unassigned', ]
ct_freq$group <- factor(ct_freq$group, levels = c('NF', 'pRV', 'RVF'))

ct_prop <- ct_freq %>%
  group_by(patient, group, cell_type) %>%
  summarise(n = n(), .groups = 'drop') %>%
  group_by(patient) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# Order cell types by overall mean proportion (descending)
ct_order <- ct_prop %>%
  group_by(cell_type) %>%
  summarise(mean_prop = mean(prop), .groups = 'drop') %>%
  arrange(desc(mean_prop)) %>%
  pull(cell_type)
ct_prop$cell_type <- factor(ct_prop$cell_type, levels = ct_order)

# (A) Combined boxplot: all cell types, colored by group
p_ct_freq <- ggplot(ct_prop, aes(x = cell_type, y = prop * 100, color = group)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.4) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75),
              size = 1.2, alpha = 0.7) +
  scale_color_manual(values = c(NF = 'steelblue', pRV = 'darkorange', RVF = 'firebrick')) +
  coord_flip() +
  theme_classic(base_size = 11) +
  theme(axis.text.y = element_text(size = 7)) +
  labs(x = NULL, y = 'Proportion of cells (%)', color = 'Group',
       title = 'Cell type frequency by disease state')

pdf('./output/Xenium/celltype_frequency_by_group.pdf', width = 10, height = 14)
print(p_ct_freq)
dev.off()
cat('Saved celltype_frequency_by_group.pdf\n')

# (B) Per-patient stacked bar chart
ct_prop_stacked <- ct_prop %>%
  mutate(cell_type = factor(cell_type, levels = rev(ct_order)))

p_ct_stacked <- ggplot(ct_prop_stacked, aes(x = patient, y = prop * 100, fill = cell_type)) +
  geom_bar(stat = 'identity', width = 0.8) +
  facet_wrap(~ group, scales = 'free_x', nrow = 1) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.3, 'cm'),
        strip.text = element_text(face = 'bold')) +
  labs(x = NULL, y = 'Proportion of cells (%)', fill = 'Cell type',
       title = 'Cell type composition per patient')

pdf('./output/Xenium/celltype_composition_per_patient.pdf', width = 12, height = 7)
print(p_ct_stacked)
dev.off()
cat('Saved celltype_composition_per_patient.pdf\n')

# (C) Faceted boxplots for top disease-variable cell types
# Select cell types with highest between-group variance relative to within-group
ct_signal <- ct_prop %>%
  group_by(cell_type, group) %>%
  summarise(group_mean = mean(prop), group_sd = sd(prop), .groups = 'drop') %>%
  group_by(cell_type) %>%
  summarise(between_cv = sd(group_mean) / mean(group_mean),
            within_cv_mean = mean(group_sd / (group_mean + 1e-6)),
            .groups = 'drop') %>%
  mutate(signal_ratio = between_cv / (within_cv_mean + 0.01)) %>%
  arrange(desc(signal_ratio))

top_signal_cts <- head(ct_signal$cell_type, 15)

ct_prop_top <- ct_prop %>% filter(cell_type %in% top_signal_cts)
ct_prop_top$cell_type <- factor(ct_prop_top$cell_type, levels = top_signal_cts)

pairwise_comparisons <- list(c('NF','pRV'), c('pRV','RVF'), c('NF','RVF'))

p_ct_top <- ggplot(ct_prop_top, aes(x = group, y = prop * 100, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, size = 1.5) +
  stat_compare_means(comparisons = pairwise_comparisons, method = 'wilcox.test',
                     label = 'p.signif', size = 3) +
  stat_compare_means(method = 'kruskal.test', label.y.npc = 0.95, size = 2.5) +
  scale_fill_manual(values = c(NF = 'steelblue', pRV = 'darkorange', RVF = 'firebrick')) +
  facet_wrap(~ cell_type, scales = 'free_y', ncol = 5) +
  theme_bw(base_size = 10) +
  theme(legend.position = 'none',
        strip.text = element_text(face = 'bold', size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = 'Proportion of cells (%)',
       title = 'Top disease-variable cell types (ranked by signal-to-noise ratio)')

pdf('./output/Xenium/celltype_top_disease_variable.pdf', width = 14, height = 10)
print(p_ct_top)
dev.off()
cat('Saved celltype_top_disease_variable.pdf\n')

# ── QC: Nuclei-level cell type proportions vs cell-level proportions ──────────
# Map each Xenium nucleus to its cell_types_manual annotation and verify that
# nucleus-level frequencies recapitulate cell-level disease trends.

# Re-use nuclei_meta from line 5070 (patient/group/file mapping)
nuclei_meta_qc <- list(
  list(file = 'output-XETG00217__0038213__Region_1__20241206__182124_top_nuclei_assignments.csv',    patient = '1697', group = 'NF'),
  list(file = 'output-XETG00217__0038213__Region_1__20241206__182124_middle_nuclei_assignments.csv', patient = '1691', group = 'NF'),
  list(file = 'output-XETG00217__0038213__Region_1__20241206__182124_bottom_nuclei_assignments.csv', patient = '1618', group = 'pRV'),
  list(file = 'output-XETG00217__0038216__Region_1__20241206__182124_top_nuclei_assignments.csv',    patient = '1567', group = 'pRV'),
  list(file = 'output-XETG00217__0038216__Region_1__20241206__182124_bottom_nuclei_assignments.csv', patient = '1692', group = 'pRV'),
  list(file = 'output-XETG00217__0038290__Region_1__20241212__142808_nuclei_assignments.csv',        patient = '1561', group = 'NF'),
  list(file = 'output-XETG00217__0038290__Region_2__20241212__142808_nuclei_assignments.csv',        patient = '1343', group = 'RVF'),
  list(file = 'output-XETG00217__0038291__Region_1__20241212__142808_nuclei_assignments.csv',        patient = '1467', group = 'RVF'),
  list(file = 'output-XETG00217__0038291__Region_2__20241212__142808_nuclei_assignments.csv',        patient = '1632', group = 'RVF')
)

# Load all nuclei assignments and join with cell_types_manual from xenium.obj
nuclei_qc <- bind_rows(lapply(nuclei_meta_qc, function(m) {
  df <- read.csv(m$file)
  df$patient <- m$patient
  df$group   <- m$group
  df
}))

# reference_cell_id = proseg cell ID; colnames(xenium.obj) = proseg ID + Seurat merge suffix
# Strip the merge suffix (_1, _2, ..._9) to recover original proseg cell ID
cell_annot <- data.frame(
  seurat_barcode = colnames(xenium.obj),
  cell_types_manual = as.character(xenium.obj$cell_types_manual),
  patient = as.character(xenium.obj$patient),
  stringsAsFactors = FALSE
)
cell_annot$proseg_id <- sub('_\\d+$', '', cell_annot$seurat_barcode)

# Match via reference_cell_id (proseg cell) + patient to avoid cross-patient collisions
nuclei_qc$reference_cell_id_str <- as.character(nuclei_qc$reference_cell_id)
nuclei_qc <- merge(nuclei_qc, cell_annot[, c('proseg_id', 'cell_types_manual', 'patient')],
                   by.x = c('reference_cell_id_str', 'patient'),
                   by.y = c('proseg_id', 'patient'),
                   all.x = TRUE)

cat('Nuclei matched to cell_types_manual:', sum(!is.na(nuclei_qc$cell_types_manual)),
    'of', nrow(nuclei_qc), '\n')
cat('Nuclei without match (not in xenium.obj):', sum(is.na(nuclei_qc$cell_types_manual)), '\n')

# Filter to assigned nuclei (exclude Unassigned from both coarse and fine annotations)
nuclei_qc <- nuclei_qc %>%
  filter(!is.na(cell_types_manual),
         cell_types_manual != 'Unassigned',
         !transferred_identity %in% c('Unassigned', 'None'))

nuclei_qc$group <- factor(nuclei_qc$group, levels = c('NF', 'pRV', 'RVF'))

# ── (A) Concordance: coarse nuclei label vs fine cell_types_manual ────────────
# Build mapping from fine cell_types_manual -> expected coarse compartment
compartment_map <- c(
  CM = 'CM', CM_RV = 'CM', CM_RV_2 = 'CM', CM_RV_3 = 'CM',
  CM_Contractile = 'CM', CM_Ventricular = 'CM', CM_Ventricular_2 = 'CM',
  CM_Ventricular_3 = 'CM', CM_Ventricular_4 = 'CM', CM_Ventricular_5 = 'CM',
  CM_Ventricular_6 = 'CM', CM_Ventricular_7 = 'CM', CM_Stressed = 'CM',
  CM_Atrial = 'CM',
  FB_Resting = 'FB', FB_Activated = 'FB', FB_NOX4 = 'FB',
  FB_PCOLCE2 = 'FB', FB_Homeostatic = 'FB', FB_NTN4 = 'FB', FB_Stress = 'FB',
  Capillary_EC = 'EC', Activated_EC = 'EC', Venous_EC = 'EC',
  Arterial_EC = 'EC', Endocardial = 'EC',
  LEC = 'LEC',
  VSMC_Synthetic = 'SM', VSMC_Arterial = 'SM', VSMC_Arterial2 = 'SM',
  VSMC_Inflamed = 'SM', VSMC = 'SM',
  Pericyte = 'PC', Pericyte_KCNJ8 = 'PC', Pericyte_THBS4 = 'PC',
  Pericyte_Quiescent = 'PC', Mural_Unspecified = 'SM',
  Macrophage_Resident = 'Myeloid', Macrophage_C1q = 'Myeloid',
  Macrophage_Inflammatory = 'Myeloid', Macrophage_TREM2 = 'Myeloid',
  Macrophage_Monocyte_Derived = 'Myeloid', Monocyte = 'Myeloid',
  Myeloid_Proliferating = 'Myeloid', Dendritic_Cell = 'Myeloid',
  Mast_Cell = 'Myeloid',
  CD8_T = 'NKT', CD4_T = 'NKT', NK = 'NKT',
  Adipo_Mature = 'Adipo', Adipo_Quiescent = 'Adipo',
  Adipo_Sparse = 'Adipo', Adipo_Sparse2 = 'Adipo', Adipo_Lipogenic = 'Adipo',
  Stromal = 'Adipo',
  Neuron = 'Neuron', Schwann_Cell = 'Neuron',
  Epicardial = 'Epi', Mesothelial = 'Epi'
)

nuclei_qc$expected_coarse <- compartment_map[nuclei_qc$cell_types_manual]
nuclei_qc$concordant <- nuclei_qc$transferred_identity == nuclei_qc$expected_coarse

concordance_rate <- mean(nuclei_qc$concordant, na.rm = TRUE)
cat(sprintf('Overall nuclei-cell compartment concordance: %.1f%%\n', concordance_rate * 100))

# Per-compartment concordance
concordance_by_compartment <- nuclei_qc %>%
  filter(!is.na(expected_coarse)) %>%
  group_by(expected_coarse) %>%
  summarise(n = n(),
            concordant = sum(concordant),
            rate = concordant / n,
            .groups = 'drop') %>%
  arrange(desc(n))
cat('\nConcordance by compartment:\n')
print(as.data.frame(concordance_by_compartment))

# Confusion matrix: transferred_identity (rows) vs expected_coarse (cols)
confusion <- table(nuclei_qc$transferred_identity, nuclei_qc$expected_coarse)
cat('\nConfusion matrix (rows = nuclei label, cols = expected from cell_types_manual):\n')
print(confusion)

# ── (B) Nucleus-level cell type proportions by disease group ──────────────────
# Compute proportions using cell_types_manual (fine-grained)
nuc_prop <- nuclei_qc %>%
  group_by(patient, group, cell_types_manual) %>%
  summarise(n = n(), .groups = 'drop') %>%
  group_by(patient) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# ── (C) Compare nucleus vs cell proportions side-by-side ──────────────────────
# cell-level proportions (ct_prop already computed above)
cell_prop_summary <- ct_prop %>%
  group_by(cell_type, group) %>%
  summarise(cell_mean_prop = mean(prop), .groups = 'drop')

nuc_prop_summary <- nuc_prop %>%
  group_by(cell_types_manual, group) %>%
  summarise(nuc_mean_prop = mean(prop), .groups = 'drop') %>%
  dplyr::rename(cell_type = cell_types_manual)

prop_compare <- merge(cell_prop_summary, nuc_prop_summary,
                      by = c('cell_type', 'group'), all = TRUE) %>%
  mutate(cell_mean_prop = ifelse(is.na(cell_mean_prop), 0, cell_mean_prop),
         nuc_mean_prop  = ifelse(is.na(nuc_mean_prop), 0, nuc_mean_prop))

# Correlation between cell and nucleus proportions
cor_val <- cor(prop_compare$cell_mean_prop, prop_compare$nuc_mean_prop, use = 'complete.obs')
cat(sprintf('\nCorrelation between cell-level and nucleus-level proportions: r = %.3f\n', cor_val))

# Scatter: cell vs nucleus proportions colored by group
p_prop_scatter <- ggplot(prop_compare, aes(x = cell_mean_prop * 100, y = nuc_mean_prop * 100,
                                           color = group)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'grey50') +
  geom_text(data = prop_compare %>% filter(cell_mean_prop > 0.03 | nuc_mean_prop > 0.03),
            aes(label = cell_type), size = 2.2, nudge_y = 0.3, show.legend = FALSE) +
  scale_color_manual(values = c(NF = 'steelblue', pRV = 'darkorange', RVF = 'firebrick')) +
  theme_bw(base_size = 11) +
  labs(x = 'Cell-level mean proportion (%)', y = 'Nucleus-level mean proportion (%)',
       color = 'Group',
       title = sprintf('Cell vs nucleus proportions (r = %.3f)', cor_val))

pdf('./output/Xenium/qc_cell_vs_nucleus_proportions.pdf', width = 8, height = 7)
print(p_prop_scatter)
dev.off()
cat('Saved qc_cell_vs_nucleus_proportions.pdf\n')

# ── (D) Nucleus-level faceted boxplots for top disease-variable cell types ────
# Re-use top_signal_cts from cell-level analysis
nuc_prop_top <- nuc_prop %>%
  filter(cell_types_manual %in% top_signal_cts) %>%
  mutate(cell_types_manual = factor(cell_types_manual, levels = top_signal_cts))

p_nuc_top <- ggplot(nuc_prop_top, aes(x = group, y = prop * 100, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, size = 1.5) +
  stat_compare_means(comparisons = pairwise_comparisons, method = 'wilcox.test',
                     label = 'p.signif', size = 3) +
  stat_compare_means(method = 'kruskal.test', label.y.npc = 0.95, size = 2.5) +
  scale_fill_manual(values = c(NF = 'steelblue', pRV = 'darkorange', RVF = 'firebrick')) +
  facet_wrap(~ cell_types_manual, scales = 'free_y', ncol = 5) +
  theme_bw(base_size = 10) +
  theme(legend.position = 'none',
        strip.text = element_text(face = 'bold', size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = 'Proportion of nuclei (%)',
       title = 'Top disease-variable cell types - NUCLEUS-level proportions')

pdf('./output/Xenium/qc_nucleus_top_disease_variable.pdf', width = 14, height = 10)
print(p_nuc_top)
dev.off()
cat('Saved qc_nucleus_top_disease_variable.pdf\n')

# ── (E) Side-by-side comparison: cell vs nucleus for key cell types ───────────
# Combine cell and nucleus proportions into one long dataframe for paired plotting
ct_prop_long <- ct_prop %>%
  filter(cell_type %in% top_signal_cts) %>%
  select(patient, group, cell_type, prop) %>%
  mutate(source = 'Cell-level')

nuc_prop_long <- nuc_prop %>%
  filter(cell_types_manual %in% top_signal_cts) %>%
  select(patient, group, cell_type = cell_types_manual, prop) %>%
  mutate(source = 'Nucleus-level')

combined_prop <- bind_rows(ct_prop_long, nuc_prop_long) %>%
  mutate(cell_type = factor(cell_type, levels = top_signal_cts),
         source = factor(source, levels = c('Cell-level', 'Nucleus-level')))

p_cell_vs_nuc <- ggplot(combined_prop, aes(x = group, y = prop * 100, fill = source)) +
  geom_boxplot(outlier.shape = NA, width = 0.7, alpha = 0.7,
               position = position_dodge(width = 0.75)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75),
              size = 1, alpha = 0.6) +
  scale_fill_manual(values = c('Cell-level' = '#4DAF4A', 'Nucleus-level' = '#984EA3')) +
  facet_wrap(~ cell_type, scales = 'free_y', ncol = 5) +
  theme_bw(base_size = 10) +
  theme(strip.text = element_text(face = 'bold', size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'bottom') +
  labs(x = NULL, y = 'Proportion (%)', fill = NULL,
       title = 'Cell vs nucleus proportions: top disease-variable types')

pdf('./output/Xenium/qc_cell_vs_nucleus_sidebyside.pdf', width = 14, height = 10)
print(p_cell_vs_nuc)
dev.off()
cat('Saved qc_cell_vs_nucleus_sidebyside.pdf\n')















#xenium.obj <- IntegrateLayers(
#  object = xenium.obj, method = HarmonyIntegration,
#  orig.reduction = "pca", new.reduction = "harmony",
#  verbose = TRUE
#)


#xenium.obj <- FindNeighbors(xenium.obj, dims = 1:30, reduction = "harmony")
#xenium.obj <- FindClusters(xenium.obj, resolution = 1)
#xenium.obj <- RunUMAP(xenium.obj, reduction = "pca", dims = 1:30,reduction.name='umap_orig')
#xenium.obj <- RunUMAP(xenium.obj, reduction = "harmony", dims = 1:30,reduction.name='umap')



#p1 <- DimPlot(xenium.obj,group.by='cell_type_rctd_doublet',reduction = 'umap', label=T)
#p2 <- DimPlot(xenium.obj,group.by='cell_type_rctd_doublet',reduction = 'umap_orig', label=T)
#wrap_plots(c(p1, p2), ncol = 2, byrow = F)


############################################## ####
#### Transmural endo -> epi gradient analysis  ####
############################################## ####

xenium.base <- readRDS('./Xenium_resegmented_corrected.rds')
panel_genes <- rownames(xenium.base)
rm(xenium.base); gc()

sn <- readRDS('../snRV_ref.rds')
Idents(sn) <- sn$Subnames_manual
all.marks <- FindAllMarkers(sn, only.pos = TRUE, min.pct = 0.25,
                            logfc.threshold = 0.5)

write.csv(all.marks,'./output/Xenium/ref.marks.csv')

## Restrict to genes on the real Xenium panel (panel_genes) — exclude imputed
## genes, which have no real spatial signal.
all.marks.panel <- all.marks[all.marks$gene %in% panel_genes, ]

## Top-ranked Endo and Epi markers restricted to the panel
endo_top <- all.marks.panel %>%
  filter(cluster == 'EC_Endocardial') %>%
  arrange(desc(avg_log2FC)) %>%
  head(20)
epi_top <- all.marks.panel %>%
  filter(cluster == 'Epi') %>%
  arrange(desc(avg_log2FC)) %>%
  head(20)

cat('--- Top Endo markers (in Xenium panel) ---\n'); print(endo_top)
cat('--- Top Epi  markers (in Xenium panel) ---\n'); print(epi_top)


library(RANN)

## Helper: retain only the `top_frac` of anchor cells nearest to the convex-hull
## boundary of all cells in the patient section, discarding internal Endo/Epi
## cells that would corrupt the transmural depth calculation.
## `min_keep` prevents over-pruning when anchors are already sparse.
.boundary_anchors <- function(anchor_cells, all_xy, anchor_xy,
                               top_frac = 0.30, min_keep = 10) {
  if (length(anchor_cells) == 0 || nrow(anchor_xy) == 0 || nrow(all_xy) < 10)
    return(anchor_cells)
  hull_verts <- all_xy[chull(all_xy), , drop = FALSE]
  d          <- RANN::nn2(hull_verts, anchor_xy, k = 1)$nn.dists[, 1]
  n_keep     <- max(min_keep, round(top_frac * length(anchor_cells)))
  n_keep     <- min(n_keep, length(anchor_cells))
  anchor_cells[order(d)[seq_len(n_keep)]]
}

## ---------------------------------------------------------------------------
## Manual polygon tracing helpers for Endo / Epi region refinement.
## The user draws polygons on a per-patient spatial plot to define which tissue
## regions are endocardial vs epicardial.  Traced polygons are cached to an RDS
## file so the interactive step only happens once.
## ---------------------------------------------------------------------------

manual_trace_cache <- './output/Xenium/transmural_manual_traces.rds'

## Ray-casting point-in-polygon test (vectorised over points).
## Returns logical vector of length(px).
.point_in_poly <- function(px, py, vx, vy) {
  n   <- length(vx)
  inside <- rep(FALSE, length(px))
  j <- n
  for (i in seq_len(n)) {
    cond <- ((vy[i] > py) != (vy[j] > py)) &
            (px < (vx[j] - vx[i]) * (py - vy[i]) / (vy[j] - vy[i]) + vx[i])
    inside <- xor(inside, cond)
    j <- i
  }
  inside
}

## Test membership against a list of polygons (union).
## `polys` is a list of data.frames each with columns x, y.
## Returns logical vector of length(x).
.cells_in_polys <- function(x, y, polys) {
  hit <- rep(FALSE, length(x))
  for (poly in polys) {
    hit <- hit | .point_in_poly(x, y, poly$x, poly$y)
  }
  hit
}

## Interactive tracing for one patient.
## Shows the spatial plot with current automated anchors, then asks the user to
## draw polygon(s) for Endo and Epi regions using locator().
## Returns list(endo = list_of_polys, epi = list_of_polys).
.trace_patient <- function(pat, pat_meta) {
  cat(sprintf('\n===== Tracing patient: %s =====\n', pat))

  ## Plot all cells grey, overlay automated anchors
  plot(pat_meta$x_centroid, -pat_meta$y_centroid,
       pch = '.', col = 'grey80', xlab = '', ylab = '',
       main = sprintf('%s — click polygons for Endo (red) then Epi (blue)', pat),
       asp = 1)
  endo_idx <- which(pat_meta$transmural_anchor == 'Endo')
  epi_idx  <- which(pat_meta$transmural_anchor == 'Epi')
  if (length(endo_idx) > 0)
    points(pat_meta$x_centroid[endo_idx], -pat_meta$y_centroid[endo_idx],
           pch = 20, cex = 0.4, col = '#e31a1c')
  if (length(epi_idx) > 0)
    points(pat_meta$x_centroid[epi_idx], -pat_meta$y_centroid[epi_idx],
           pch = 20, cex = 0.4, col = '#1f78b4')
  legend('topright', legend = c('Endo (auto)', 'Epi (auto)', 'Other'),
         col = c('#e31a1c', '#1f78b4', 'grey80'), pch = 20, cex = 0.8)

  .draw_polys <- function(label, col) {
    polys <- list()
    repeat {
      cat(sprintf('  Draw %s polygon #%d (right-click / Esc to finish this polygon):\n',
                  label, length(polys) + 1))
      loc <- locator(type = 'l', col = col, lwd = 2)
      if (is.null(loc) || length(loc$x) < 3) {
        if (length(polys) == 0) {
          cat(sprintf('  No %s polygons drawn.\n', label))
        }
        break
      }
      ## Close the polygon visually
      lines(c(loc$x, loc$x[1]), c(loc$y, loc$y[1]), col = col, lwd = 2)
      polys[[length(polys) + 1]] <- data.frame(x = loc$x, y = loc$y)
      cat(sprintf('  %s polygon with %d vertices saved. Draw another? (right-click/Esc to stop)\n',
                  label, length(loc$x)))
    }
    polys
  }

  cat('  --- Step A: draw ENDO region(s) in red ---\n')
  cat('  Click vertices to outline the endocardial surface region.\n')
  cat('  Right-click (or press Esc) to close a polygon. Repeat for additional polygons.\n')
  cat('  Right-click/Esc with <3 points to finish Endo and move to Epi.\n')
  endo_polys <- .draw_polys('Endo', '#e31a1c')

  cat('  --- Step B: draw EPI region(s) in blue ---\n')
  cat('  Click vertices to outline the epicardial surface region.\n')
  cat('  Right-click/Esc with <3 points to finish Epi.\n')
  epi_polys <- .draw_polys('Epi', '#1f78b4')

  list(endo = endo_polys, epi = epi_polys)
}

## Wrapper: trace selected patients, merge with existing cache, save.
## `patients`  — character vector of patient IDs to trace (default: all)
## `retrace`   — if TRUE, re-trace even if cached polygons already exist
## Run interactively: trace_patients()          # trace all missing
##                    trace_patients('P1')       # trace one patient
##                    trace_patients(retrace=T)  # redo all from scratch
trace_patients <- function(patients = unique(xenium.obj$patient),
                           retrace = FALSE) {
  ## Load existing cache
  if (file.exists(manual_trace_cache) && !retrace) {
    manual_traces <- readRDS(manual_trace_cache)
    cat(sprintf('Loaded %d cached traces from %s\n',
                length(manual_traces), manual_trace_cache))
  } else {
    manual_traces <- list()
  }

  for (pat in patients) {
    if (!retrace && !is.null(manual_traces[[pat]])) {
      cat(sprintf('Skipping %s (already cached). Use retrace=TRUE to redo.\n', pat))
      next
    }
    pat_cells <- WhichCells(xenium.obj, expression = patient == pat)
    pat_meta  <- xenium.obj@meta.data[pat_cells, ]
    traces <- .trace_patient(pat, pat_meta)
    manual_traces[[pat]] <- traces
  }

  saveRDS(manual_traces, manual_trace_cache)
  cat(sprintf('Saved %d patient traces to %s\n',
              length(manual_traces), manual_trace_cache))
  invisible(manual_traces)
}

## Load existing manual traces (empty list if none cached)
if (file.exists(manual_trace_cache)) {
  manual_traces <- readRDS(manual_trace_cache)
  cat(sprintf('Loaded %d manual trace(s) from cache.\n', length(manual_traces)))
} else {
  manual_traces <- list()
  cat('No manual traces cached. Run trace_patients() to create them.\n')
}

## --- Step 1: per-cell transmural depth = d_endo / (d_endo + d_epi) ---

## Compute per-cell marker scores once across all cells.
##
## Endo: NPPA and NPPB are natriuretic peptides expressed in subendocardial /
## trabecular cardiomyocytes (NPPA-high) rather than in the NPPA-low compact
## myocardium near the epi surface. They are NOT endocardial-EC markers, but
## because trabecular CMs form a broad, dense layer along the inner ventricular
## wall, they robustly trace the endocardial side. They also do not cross-react
## with coronary vasculature (unlike PECAM1/VWF/FLT1), which was the failure
## mode of previous pan-endothelial sets.
##
## Epi: top-ranked FindAllMarkers hits for the snRV Epi cluster, restricted to
## panel_genes and filtered for specificity.
endo_markers <- c('NPPA', 'NPPB')
epi_markers  <- c('PRG4', 'CFTR', 'WT1', 'CFB', 'ROR2', 'PLA2G2A')
.score_mat <- GetAssayData(xenium.obj, assay = 'Xenium', layer = 'data')
xenium.obj$endo_score <- Matrix::colMeans(.score_mat[endo_markers, , drop = FALSE])
xenium.obj$epi_score  <- Matrix::colMeans(.score_mat[epi_markers,  , drop = FALSE])
rm(.score_mat); gc()

xenium.obj$transmural_depth  <- NA_real_
xenium.obj$transmural_anchor <- 'Other'   # tracks final anchors for QC plot

anchor_summary <- list()
for (pat in unique(xenium.obj$patient)) {
  pat_cells <- WhichCells(xenium.obj, expression = patient == pat)
  pat_meta  <- xenium.obj@meta.data[pat_cells, ]

  if (!is.null(manual_traces[[pat]])) {
    ## ---- Manual polygon mode: use cached hand-drawn regions ----
    cat(sprintf('  Using manual traces for %s\n', pat))
    traces <- manual_traces[[pat]]
    px <- pat_meta$x_centroid
    ## NOTE: polygons are drawn on the plot where y = -y_centroid,
    ## so test against -y_centroid to match the drawn coordinates.
    py <- -pat_meta$y_centroid

    endo_mask <- if (length(traces$endo) > 0) .cells_in_polys(px, py, traces$endo) else rep(FALSE, length(px))
    epi_mask  <- if (length(traces$epi)  > 0) .cells_in_polys(px, py, traces$epi)  else rep(FALSE, length(px))

    ## Resolve any overlap: if a cell falls in both, assign to whichever marker
    ## score is higher (endo_score vs epi_score)
    both <- endo_mask & epi_mask
    if (any(both)) {
      endo_mask[both] <- pat_meta$endo_score[both] >= pat_meta$epi_score[both]
      epi_mask[both]  <- !endo_mask[both]
    }

    endo_cells <- pat_cells[endo_mask]
    epi_cells  <- pat_cells[epi_mask]

  } else {
    ## ---- Automatic mode: marker scores + smoothing + mutual exclusion ----

    ## Spatially smooth the marker scores within this patient so that isolated
    ## vessel cells drop out and only coherent regions of co-expression survive.
    ## Each cell's score becomes the mean over its k nearest spatial neighbours.
    all_xy   <- as.matrix(pat_meta[, c('x_centroid', 'y_centroid')])
    k_smooth <- min(30, nrow(all_xy) - 1)
    nn_idx   <- RANN::nn2(all_xy, all_xy, k = k_smooth)$nn.idx
    endo_smooth <- rowMeans(matrix(pat_meta$endo_score[nn_idx], ncol = k_smooth))
    epi_smooth  <- rowMeans(matrix(pat_meta$epi_score [nn_idx], ncol = k_smooth))

    ## Select anchors by smoothed score differential to enforce mutual exclusivity:
    ## a cell must score higher for one marker set than the other
    endo_diff  <- endo_smooth - epi_smooth
    epi_diff   <- epi_smooth  - endo_smooth
    endo_cells <- pat_cells[endo_diff >= quantile(endo_diff, 0.95)]
    epi_cells  <- pat_cells[epi_diff  >= quantile(epi_diff,  0.95)]

    ## Spatial mutual exclusion: drop any anchor within `excl_radius_um` of an
    ## opposite-class anchor. These are contested/transitional regions that would
    ## otherwise corrupt the d_endo / d_epi ratio.
    excl_radius_um <- 150
    if (length(endo_cells) > 0 && length(epi_cells) > 0) {
      endo_xy_f <- as.matrix(xenium.obj@meta.data[endo_cells, c('x_centroid', 'y_centroid')])
      epi_xy_f  <- as.matrix(xenium.obj@meta.data[epi_cells,  c('x_centroid', 'y_centroid')])
      d_endo <- RANN::nn2(epi_xy_f,  endo_xy_f, k = 1)$nn.dists[, 1]
      d_epi  <- RANN::nn2(endo_xy_f, epi_xy_f,  k = 1)$nn.dists[, 1]
      endo_cells <- endo_cells[d_endo > excl_radius_um]
      epi_cells  <- epi_cells [d_epi  > excl_radius_um]
    }
  }

  xenium.obj$transmural_anchor[endo_cells] <- 'Endo'
  xenium.obj$transmural_anchor[epi_cells]  <- 'Epi'

  anchor_summary[[pat]] <- data.frame(
    patient = pat,
    n_cells = length(pat_cells),
    n_endo  = length(endo_cells),
    n_epi   = length(epi_cells),
    stringsAsFactors = FALSE
  )
  cat(sprintf('%s: n=%d  Endo=%d  Epi=%d\n',
              pat, length(pat_cells), length(endo_cells), length(epi_cells)))

  if (length(endo_cells) < 5 || length(epi_cells) < 5) {
    cat('  Skipping - insufficient Endo/Epi anchors\n'); next
  }

  all_coords  <- as.matrix(pat_meta[, c('x_centroid', 'y_centroid')])
  endo_coords <- as.matrix(xenium.obj@meta.data[endo_cells, c('x_centroid', 'y_centroid')])
  epi_coords  <- as.matrix(xenium.obj@meta.data[epi_cells,  c('x_centroid', 'y_centroid')])

  d_endo <- RANN::nn2(endo_coords, all_coords, k = 1)$nn.dists[, 1]
  d_epi  <- RANN::nn2(epi_coords,  all_coords, k = 1)$nn.dists[, 1]

  xenium.obj$transmural_depth[pat_cells] <- d_endo / (d_endo + d_epi)
}
anchor_summary <- do.call(rbind, anchor_summary)
write.csv(anchor_summary, './output/Xenium/transmural_anchor_summary.csv', row.names = FALSE)

## --- Step 1b: QC - visualise Endo/Epi anchor locations per patient ---
# All cells shown as light grey background; Endo and Epi anchors overlaid in
# distinct colours so the user can verify they trace the correct walls before
# trusting the transmural depth computation.
all_meta <- xenium.obj@meta.data[, c('patient','x_centroid','y_centroid',
                                      'transmural_anchor')]
all_meta$layer <- factor(all_meta$transmural_anchor, levels = c('Other','Endo','Epi'))

p_anchors <- ggplot(all_meta, aes(x = x_centroid, y = -y_centroid,
                                   color = layer, size = layer, alpha = layer)) +
  geom_point() +
  scale_color_manual(values = c(Other = 'grey85', Endo = '#e31a1c', Epi = '#1f78b4')) +
  scale_size_manual(values  = c(Other = 0.05,     Endo = 0.6,       Epi = 0.6)) +
  scale_alpha_manual(values = c(Other = 0.3,      Endo = 1,         Epi = 1)) +
  facet_wrap(~ patient, scales = 'free', ncol = 3) +
  labs(title = 'Endo (red) and Epi (blue) anchor cells per patient - QC',
       subtitle = 'Verify red traces the endocardial surface and blue the epicardial surface',
       x = NULL, y = NULL, color = NULL) +
  guides(size = 'none', alpha = 'none',
         color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme_bw(base_size = 10) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank())

ggsave('./output/Xenium/transmural_anchor_qc.pdf',
       ggrastr::rasterise(p_anchors, layers = 'Point', dpi = 300),
       width = 14, height = 12, device = cairo_pdf)

## --- Step 2: sanity-check spatial scatter colored by depth ---
depth_df <- xenium.obj@meta.data[, c('patient', 'x_centroid', 'y_centroid',
                                     'transmural_depth', 'cell_type_rctd_doublet')]
depth_df <- depth_df[!is.na(depth_df$transmural_depth), ]

pdf('./output/Xenium/transmural_depth_per_patient.pdf', width = 14, height = 12)
print(
  ggplot(depth_df, aes(x = x_centroid, y = y_centroid, color = transmural_depth)) +
    geom_point(size = 0.15, alpha = 0.6) +
    scale_color_gradientn(colors = c('navy', 'deepskyblue', 'lightyellow', 'orange', 'firebrick')) +
    facet_wrap(~ patient, scales = 'free', ncol = 3) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(title = 'Transmural depth (0 = endocardium, 1 = epicardium)',
         color = 'depth')
)
# depth histogram per patient
print(
  ggplot(depth_df, aes(x = transmural_depth)) +
    geom_histogram(bins = 40, fill = 'steelblue', color = 'white') +
    facet_wrap(~ patient, scales = 'free_y', ncol = 3) +
    theme_bw() +
    labs(title = 'Transmural depth distribution per patient')
)
dev.off()
rm(depth_df); gc()

## --- Step 3: per-gene Spearman rho of bin-averaged expression per patient ---
n_bins   <- 20
expr_mat <- GetAssayData(xenium.obj, assay = 'Xenium', layer = 'data')   # genes x cells
all_genes <- rownames(expr_mat)
patients  <- sort(unique(xenium.obj$patient))

patient_rho <- matrix(NA_real_,
                      nrow = length(all_genes),
                      ncol = length(patients),
                      dimnames = list(all_genes, patients))

# Also store bin-averaged expression per patient for heatmap/line plots
bin_expr_by_patient <- vector('list', length(patients))
names(bin_expr_by_patient) <- patients

for (pat in patients) {
  cells <- WhichCells(xenium.obj, expression = patient == pat)
  depth <- xenium.obj$transmural_depth[cells]
  valid <- !is.na(depth)
  cells <- cells[valid]; depth <- depth[valid]
  if (length(cells) < 1000) { cat(pat, '- too few cells, skipping\n'); next }

  bin_idx <- cut(depth, breaks = seq(0, 1, length.out = n_bins + 1),
                 include.lowest = TRUE, labels = FALSE)

  bin_expr <- sapply(seq_len(n_bins), function(b) {
    cb <- cells[bin_idx == b]
    if (length(cb) < 10) return(rep(NA_real_, length(all_genes)))
    Matrix::rowMeans(expr_mat[, cb, drop = FALSE])
  })
  colnames(bin_expr) <- paste0('bin', sprintf('%02d', seq_len(n_bins)))
  bin_expr_by_patient[[pat]] <- bin_expr

  rho <- apply(bin_expr, 1, function(g) {
    if (all(is.na(g))) return(NA_real_)
    suppressWarnings(cor(g, seq_len(n_bins), method = 'spearman', use = 'complete.obs'))
  })
  patient_rho[, pat] <- rho
  cat(pat, '- done (', length(cells), 'cells )\n')
}

saveRDS(list(patient_rho = patient_rho,
             bin_expr_by_patient = bin_expr_by_patient),
        './output/Xenium/transmural_gradient_raw.rds')

## --- Step 4: pooled aggregation across all patients ---
patient_meta <- unique(xenium.obj@meta.data[, c('patient', 'group')])
rownames(patient_meta) <- patient_meta$patient

gradient_summary <- data.frame(
  gene          = all_genes,
  median_rho    = apply(patient_rho, 1, median, na.rm = TRUE),
  mean_rho      = apply(patient_rho, 1, mean,   na.rm = TRUE),
  n_pos         = apply(patient_rho, 1, function(r) sum(r >  0.3, na.rm = TRUE)),
  n_neg         = apply(patient_rho, 1, function(r) sum(r < -0.3, na.rm = TRUE)),
  n_patients    = apply(patient_rho, 1, function(r) sum(!is.na(r))),
  is_panel_gene = all_genes %in% panel_genes,
  stringsAsFactors = FALSE
)
gradient_summary$consistency <- with(gradient_summary,
                                     pmax(n_pos, n_neg) / pmax(n_patients, 1))
gradient_summary <- gradient_summary[order(-abs(gradient_summary$median_rho)), ]
write.csv(gradient_summary, './output/Xenium/transmural_gradient_pooled.csv', row.names = FALSE)

## --- Step 5: per-disease-group aggregation ---
for (grp in c('NF', 'pRV', 'RVF')) {
  grp_pats <- patient_meta$patient[patient_meta$group == grp]
  grp_pats <- intersect(grp_pats, colnames(patient_rho))
  if (length(grp_pats) == 0) next
  grp_rho  <- patient_rho[, grp_pats, drop = FALSE]
  grp_summary <- data.frame(
    gene          = all_genes,
    median_rho    = apply(grp_rho, 1, median, na.rm = TRUE),
    mean_rho      = apply(grp_rho, 1, mean,   na.rm = TRUE),
    n_pos         = apply(grp_rho, 1, function(r) sum(r >  0.3, na.rm = TRUE)),
    n_neg         = apply(grp_rho, 1, function(r) sum(r < -0.3, na.rm = TRUE)),
    group         = grp,
    is_panel_gene = all_genes %in% panel_genes,
    stringsAsFactors = FALSE
  )
  grp_summary <- grp_summary[order(-abs(grp_summary$median_rho)), ]
  write.csv(grp_summary,
            paste0('./output/Xenium/transmural_gradient_', grp, '.csv'),
            row.names = FALSE)
}

## --- Step 6: disease-altered gradient genes (delta rho between groups) ---
rho_by_group <- sapply(c('NF', 'pRV', 'RVF'), function(grp) {
  gp <- intersect(patient_meta$patient[patient_meta$group == grp], colnames(patient_rho))
  if (length(gp) == 0) return(rep(NA_real_, length(all_genes)))
  apply(patient_rho[, gp, drop = FALSE], 1, median, na.rm = TRUE)
})
colnames(rho_by_group) <- paste0('rho_', colnames(rho_by_group))
rho_by_group <- as.data.frame(rho_by_group)
rho_by_group$gene          <- all_genes
rho_by_group$delta_pRV_NF  <- rho_by_group$rho_pRV - rho_by_group$rho_NF
rho_by_group$delta_RVF_NF  <- rho_by_group$rho_RVF - rho_by_group$rho_NF
rho_by_group$delta_RVF_pRV <- rho_by_group$rho_RVF - rho_by_group$rho_pRV
rho_by_group$is_panel_gene <- all_genes %in% panel_genes
rho_by_group <- rho_by_group[order(-abs(rho_by_group$delta_RVF_NF)), ]
write.csv(rho_by_group, './output/Xenium/transmural_gradient_by_group.csv', row.names = FALSE)

## --- Step 7: visualizations ---

# (a) Heatmap of top gradient genes: rows = top genes, cols = bins x patients
top_pos <- gradient_summary %>%
  filter(consistency >= 0.67, median_rho >  0.3) %>%
  arrange(desc(median_rho)) %>% head(20) %>% pull(gene)
top_neg <- gradient_summary %>%
  filter(consistency >= 0.67, median_rho < -0.3) %>%
  arrange(median_rho) %>% head(20) %>% pull(gene)
top_genes <- c(top_neg, top_pos)

# stack bin_expr_by_patient into a wide matrix: genes x (patient.bin)
stacked <- do.call(cbind, lapply(patients, function(pat) {
  be <- bin_expr_by_patient[[pat]]
  if (is.null(be)) return(NULL)
  colnames(be) <- paste0(pat, '.', colnames(be))
  be
}))
# z-score per gene across all bins
z_rows <- t(scale(t(stacked[top_genes, , drop = FALSE])))
z_rows[is.na(z_rows)] <- 0

col_split <- sub('\\..*$', '', colnames(z_rows))
col_split <- factor(col_split, levels = patients)

pdf('./output/Xenium/transmural_gradient_heatmap.pdf', width = 16, height = 10)
tryCatch({
  library(ComplexHeatmap)
  ht <- Heatmap(z_rows,
                name = 'z-score',
                cluster_rows = FALSE, cluster_columns = FALSE,
                column_split = col_split,
                show_column_names = FALSE,
                row_names_side = 'left',
                row_split = factor(c(rep('endo-biased', length(top_neg)),
                                     rep('epi-biased',  length(top_pos))),
                                   levels = c('endo-biased', 'epi-biased')),
                column_title_gp = gpar(fontsize = 9),
                row_title_gp = gpar(fontsize = 10, fontface = 'bold'))
  draw(ht)
}, error = function(e) {
  cat('ComplexHeatmap failed:', conditionMessage(e), '- using pheatmap\n')
  pheatmap::pheatmap(z_rows, cluster_rows = FALSE, cluster_cols = FALSE,
                     show_colnames = FALSE, fontsize_row = 8)
})
dev.off()

# (b) Per-group line plots: bin-averaged expression vs depth, top genes, faceted
line_df <- do.call(rbind, lapply(patients, function(pat) {
  be <- bin_expr_by_patient[[pat]]
  if (is.null(be)) return(NULL)
  df <- reshape2::melt(be[top_genes, , drop = FALSE],
                       varnames = c('gene', 'bin'), value.name = 'expr')
  df$bin_idx <- as.integer(sub('bin', '', df$bin))
  df$patient <- pat
  df$group   <- patient_meta[pat, 'group']
  df
}))
line_summary <- line_df %>%
  group_by(gene, bin_idx, group) %>%
  summarise(mean_expr = mean(expr, na.rm = TRUE),
            se        = sd(expr, na.rm = TRUE) / sqrt(sum(!is.na(expr))),
            .groups   = 'drop') %>%
  mutate(gene = factor(gene, levels = top_genes),
         group = factor(group, levels = c('NF', 'pRV', 'RVF')))

pdf('./output/Xenium/transmural_gradient_top_lineplots.pdf', width = 16, height = 12)
print(
  ggplot(line_summary, aes(x = bin_idx, y = mean_expr, color = group, fill = group)) +
    geom_ribbon(aes(ymin = mean_expr - se, ymax = mean_expr + se),
                alpha = 0.2, color = NA) +
    geom_line(linewidth = 0.7) +
    facet_wrap(~ gene, scales = 'free_y', ncol = 5) +
    scale_color_manual(values = c(NF = 'steelblue', pRV = 'darkorange', RVF = 'firebrick')) +
    scale_fill_manual(values  = c(NF = 'steelblue', pRV = 'darkorange', RVF = 'firebrick')) +
    theme_bw() +
    labs(x = 'Transmural depth bin (endo -> epi)', y = 'Mean log expression',
         title = 'Top transmural gradient genes by disease group')
)
dev.off()

# (c) Spatial maps of top 6 gradient genes across patients
top_vis <- c(head(top_neg, 3), head(top_pos, 3))
vis_meta <- xenium.obj@meta.data[, c('patient', 'x_centroid', 'y_centroid')]
vis_meta$cell <- rownames(vis_meta)

pdf('./output/Xenium/transmural_top_genes_spatial.pdf', width = 16, height = 10)
for (g in top_vis) {
  if (!g %in% rownames(expr_mat)) next
  vals <- as.numeric(as.matrix(expr_mat[g, , drop = FALSE]))
  names(vals) <- colnames(expr_mat)
  vis_meta$expr <- vals[vis_meta$cell]
  p <- ggplot(vis_meta, aes(x = x_centroid, y = y_centroid, color = expr)) +
    geom_point(size = 0.15, alpha = 0.6) +
    scale_color_gradientn(colors = c('grey90', 'orange', 'firebrick')) +
    facet_wrap(~ patient, scales = 'free', ncol = 3) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(title = paste0(g, ' - transmural spatial distribution'))
  print(p)
  print(g)
}
dev.off()

## --- Step 8: enrichR on top gradient gene sets ---
endo_biased <- gradient_summary %>%
  filter(median_rho < -0.5, consistency >= 0.67) %>%
  arrange(median_rho) %>% head(200) %>% pull(gene)
epi_biased  <- gradient_summary %>%
  filter(median_rho >  0.5, consistency >= 0.67) %>%
  arrange(desc(median_rho)) %>% head(200) %>% pull(gene)
rvf_altered <- rho_by_group %>%
  filter(abs(delta_RVF_NF) > 0.5) %>%
  arrange(desc(abs(delta_RVF_NF))) %>% head(200) %>% pull(gene)

grad_sets <- list(
  endo_biased = endo_biased,
  epi_biased  = epi_biased,
  rvf_altered = rvf_altered
)

grad_enrich_all <- list()
pdf('./output/Xenium/transmural_enrichr_plots.pdf', width = 10, height = 6)
for (sn in names(grad_sets)) {
  genes <- grad_sets[[sn]]
  cat(sn, ':', length(genes), 'genes\n')
  if (length(genes) < 5) next
  er <- tryCatch(enrichr(genes, enrichR_dbs), error = function(e) NULL)
  if (is.null(er)) next
  grad_enrich_all[[sn]] <- bind_rows(
    combine_enrich(er, direction = sn, comparison = 'transmural_gradient')
  )
  p1 <- plot_enrich(er$GO_Biological_Process_2023, paste0(sn, ' - GO BP'))
  p2 <- plot_enrich(er$Reactome_2022,              paste0(sn, ' - Reactome'))
  if (!is.null(p1)) print(p1)
  if (!is.null(p2)) print(p2)
}
dev.off()

grad_enrich_df <- bind_rows(grad_enrich_all)
write.csv(grad_enrich_df, './output/Xenium/transmural_enrichr.csv', row.names = FALSE)


