library(Seurat)
library(harmony)
library(spacexr)
library(FNN)
library(reticulate)

use_python("/Users/ikuz/.pyenv/versions/3.11.6/bin/python3", required = TRUE)

# Scrublet doublet detection (direct reticulate call)
scrub <- reticulate::import("scrublet")
scipy <- reticulate::import("scipy.sparse")

##########################
##### Myeloid labels #####
##########################

M1 <- readRDS(file = "./dependencies/shared/myeloid_subclust_new.rds")


df <- data.frame(colnames(M1),M1$cluster_id)

write.csv(df,'./dependencies/shared/myeloid_subclustering_new.csv')


#####################
##### EC labels #####
#####################

M1 <- readRDS(file = "./dependencies/shared/EC_hdWGCNA_by_celltype.rds")

df <- data.frame(colnames(M1),M1$Names)

write.csv(df,'./dependencies/shared/ec_subclustering.csv')

#####################
##### CM labels #####
#####################

M1 <- readRDS('./dependencies/shared/cm_subclust_new_new.rds')
df <- data.frame(colnames(M1),M1@active.ident)
write.csv(df,'./dependencies/shared/cm_subclustering.csv')


#####################
##### FB labels #####
#####################

M1 <- readRDS(file = "./dependencies/shared/fb_subclust.rds")

new.cluster.ids <- c("Fb_Resident","Fb_Adventitial","Fb_Elastogenic",
                     "Fb_Interstitial","Fb_Stressed","Fb_Pro-fibrotic",
                     "Fb_Anti-fibrotic")
names(new.cluster.ids) <- levels(M1)
M1 <- RenameIdents(M1, new.cluster.ids)

df <- data.frame(colnames(M1),M1@active.ident)

write.csv(df,'./dependencies/shared/fb_subclustering.csv')


########################
##### Mural labels #####
########################


M1 <- readRDS(file = "./dependencies/shared/pc_sm_subclust.rds")


new.cluster.ids <- c("Sm","Pc","Pc","Sm")

names(new.cluster.ids) <- levels(M1)
M1 <- RenameIdents(M1, new.cluster.ids)


df <- data.frame(colnames(M1),M1@active.ident)

write.csv(df,'./dependencies/shared/mural_subclustering.csv')


##################
##### Assign #####
##################

M1 <- readRDS('./dependencies/shared/snRV_ref.rds')

fb <- read.csv('./dependencies/shared/fb_subclustering.csv',header=T,row.names = 1)
cm <- read.csv('./dependencies/shared/cm_subclustering.csv',header=T,row.names = 1)
ec <- read.csv('./dependencies/shared/ec_subclustering.csv',header=T,row.names = 1)
myeloid <- read.csv('./dependencies/shared/myeloid_subclustering_new.csv',header=T,row.names = 1)
mural <- read.csv('./dependencies/shared/mural_subclustering.csv',header=T,row.names = 1)


M1$Subnames <- paste0(M1$Names,'_unc')

M1$Subnames[rownames(fb)] <- fb$'M1.active.ident'
M1$Subnames[rownames(cm)] <- cm$'M1.active.ident'
M1$Subnames[rownames(myeloid)] <- myeloid$'M1.cluster_id'
M1$Subnames[rownames(ec)] <- paste0('EC_',ec$'M1.Names')

M1$Subnames[M1$Subnames == 'LEC_unc'] = 'EC_Lymph'
M1$Subnames[M1$Subnames == 'Endo_unc'] = 'EC_Endocardial'
M1$Subnames[M1$Subnames == 'NKT_unc'] = 'NK_T'
M1$Subnames[M1$Subnames == 'Neuron_unc'] = 'Neuron'
M1$Subnames[M1$Subnames == 'PC_unc'] = 'PC'
M1$Subnames[M1$Subnames == 'SM_unc'] = 'SM'
M1$Subnames[M1$Subnames == 'Adipo_unc'] = 'Adipo'
M1$Subnames[M1$Subnames == 'Epi_unc'] = 'Epi'

saveRDS(M1,'./output/snRV_ref.rds')




#############################
##### Manual correction #####
#############################

M1 <- NormalizeData(M1)
M1 <- FindVariableFeatures(M1)
M1 <- ScaleData(M1)
M1 <- RunPCA(M1)
M1 <- RunHarmony(M1,'patient')
M1 <- RunUMAP(M1,reduction='harmony',dim=1:30)
M1$Subnames_manual <- M1$Subnames
M1$scrublet_call <- 'Singlet'
M1$scrublet_score <- 0


DimPlot(M1,label=T,group.by='Names')
DimPlot(M1,label=T,group.by='Subnames')

# EC
ec <- subset(M1,Names %in% c('EC','Endo','LEC'))
ec <- NormalizeData(ec)
ec <- FindVariableFeatures(ec)
ec <- ScaleData(ec)
ec <- RunPCA(ec)
ec <- RunHarmony(ec,'patient')
ec <- RunUMAP(ec,reduction='harmony',dim=1:30)
ec <- FindNeighbors(ec,dim=1:30)
ec <- FindClusters(ec,resolution=.2)


p1 <- DimPlot(ec,label=T,group.by='Subnames')
ec$Subnames_manual <- ec$Subnames

ec_known <- subset(ec, Subnames != 'EC_unc')

# KNN label transfer from ec_known to EC_unc cells
ec_unc <- subset(ec, Subnames == 'EC_unc')

ref_embeddings <- Embeddings(ec_known, reduction = 'harmony')[, 1:30]
query_embeddings <- Embeddings(ec_unc, reduction = 'harmony')[, 1:30]

k <- 15
nn <- FNN::get.knnx(ref_embeddings, query_embeddings, k = k)

ref_labels <- ec_known$Subnames
predicted_labels <- apply(nn$nn.index, 1, function(idx) {
  neighbor_labels <- ref_labels[idx]
  names(sort(table(neighbor_labels), decreasing = TRUE))[1]
})

ec$Subnames_manual[colnames(ec_unc)] <- predicted_labels

p2 <- DimPlot(ec, label = T, group.by = 'Subnames_manual')
p1 + p2



counts <- as(GetAssayData(ec, layer = "counts"), "dgCMatrix")
counts_t <- Matrix::t(counts)  # cells x genes for scrublet
counts_py <- scipy$csc_matrix(
  tuple(list(as.numeric(counts_t@x), as.integer(counts_t@i), as.integer(counts_t@p))),
  shape = tuple(as.integer(nrow(counts_t)), as.integer(ncol(counts_t)))
)

scrublet_obj <- scrub$Scrublet(counts_py)
result <- scrublet_obj$scrub_doublets()

ec[["scrublet_score"]] <- result[[1]]
ec[["scrublet_call"]] <- ifelse(result[[2]], "Doublet", "Singlet")

table(ec$scrublet_call)
FeaturePlot(ec, features = 'scrublet_score')
DimPlot(ec, group.by = 'scrublet_call')


M1$Subnames_manual[colnames(ec)] <- unlist(ec[["Subnames_manual"]])
M1$scrublet_score[colnames(ec)] <- unlist(ec[["scrublet_score"]])
M1$scrublet_call[colnames(ec)] <- unlist(ec[["scrublet_call"]])

# CM

ec <- subset(M1,Names %in% c('CM'))
ec <- NormalizeData(ec)
ec <- FindVariableFeatures(ec)
ec <- ScaleData(ec)
ec <- RunPCA(ec)
ec <- RunHarmony(ec,'patient')
ec <- RunUMAP(ec,reduction='harmony',dim=1:30)
ec <- FindNeighbors(ec,dim=1:30)
ec <- FindClusters(ec,resolution=.2)


p1 <- DimPlot(ec,label=T,group.by='Subnames')

ec_known <- subset(ec, Subnames != 'CM_unc')

# KNN label transfer from ec_known to CM_unc cells
ec_unc <- subset(ec, Subnames == 'CM_unc')

ref_embeddings <- Embeddings(ec_known, reduction = 'harmony')[, 1:30]
query_embeddings <- Embeddings(ec_unc, reduction = 'harmony')[, 1:30]

k <- 15
nn <- FNN::get.knnx(ref_embeddings, query_embeddings, k = k)

ref_labels <- ec_known$Subnames
predicted_labels <- apply(nn$nn.index, 1, function(idx) {
  neighbor_labels <- ref_labels[idx]
  names(sort(table(neighbor_labels), decreasing = TRUE))[1]
})

ec$Subnames_manual[colnames(ec_unc)] <- predicted_labels

p2 <- DimPlot(ec, label = T, group.by = 'Subnames_manual')
p1 + p2

# Scrublet doublet detection (direct reticulate call)
counts <- as(GetAssayData(ec, layer = "counts"), "dgCMatrix")
counts_t <- Matrix::t(counts)  # cells x genes for scrublet
counts_py <- scipy$csc_matrix(
  tuple(list(as.numeric(counts_t@x), as.integer(counts_t@i), as.integer(counts_t@p))),
  shape = tuple(as.integer(nrow(counts_t)), as.integer(ncol(counts_t)))
)

scrublet_obj <- scrub$Scrublet(counts_py)
result <- scrublet_obj$scrub_doublets()

ec[["scrublet_score"]] <- result[[1]]
ec[["scrublet_call"]] <- ifelse(result[[2]], "Doublet", "Singlet")

table(ec$scrublet_call)
FeaturePlot(ec, features = 'scrublet_score')
DimPlot(ec, group.by = 'scrublet_call')

M1$Subnames_manual[colnames(ec)] <- unlist(ec[["Subnames_manual"]])
M1$scrublet_score[colnames(ec)] <- unlist(ec[["scrublet_score"]])
M1$scrublet_call[colnames(ec)] <- unlist(ec[["scrublet_call"]])


# Myeloid

ec <- subset(M1,Names %in% c('Myeloid'))
ec <- NormalizeData(ec)
ec <- FindVariableFeatures(ec)
ec <- ScaleData(ec)
ec <- RunPCA(ec)
ec <- RunHarmony(ec,'patient')
ec <- RunUMAP(ec,reduction='harmony',dim=1:30)
ec <- FindNeighbors(ec,dim=1:30)
ec <- FindClusters(ec,resolution=.2)


p1 <- DimPlot(ec,label=T,group.by='Subnames')

ec_known <- subset(ec, Subnames != 'Myeloid_unc')

# KNN label transfer from ec_known to Myeloid_unc cells
ec_unc <- subset(ec, Subnames == 'Myeloid_unc')

ref_embeddings <- Embeddings(ec_known, reduction = 'harmony')[, 1:30]
query_embeddings <- Embeddings(ec_unc, reduction = 'harmony')[, 1:30]

k <- 15
nn <- FNN::get.knnx(ref_embeddings, query_embeddings, k = k)

ref_labels <- ec_known$Subnames
predicted_labels <- apply(nn$nn.index, 1, function(idx) {
  neighbor_labels <- ref_labels[idx]
  names(sort(table(neighbor_labels), decreasing = TRUE))[1]
})

ec$Subnames_manual[colnames(ec_unc)] <- predicted_labels

p2 <- DimPlot(ec, label = T, group.by = 'Subnames_manual')
p1 + p2

# Scrublet doublet detection (direct reticulate call)
counts <- as(GetAssayData(ec, layer = "counts"), "dgCMatrix")
counts_t <- Matrix::t(counts)  # cells x genes for scrublet
counts_py <- scipy$csc_matrix(
  tuple(list(as.numeric(counts_t@x), as.integer(counts_t@i), as.integer(counts_t@p))),
  shape = tuple(as.integer(nrow(counts_t)), as.integer(ncol(counts_t)))
)

scrublet_obj <- scrub$Scrublet(counts_py)
result <- scrublet_obj$scrub_doublets()

ec[["scrublet_score"]] <- result[[1]]
ec[["scrublet_call"]] <- ifelse(result[[2]], "Doublet", "Singlet")

table(ec$scrublet_call)
FeaturePlot(ec, features = 'scrublet_score')
DimPlot(ec, group.by = 'scrublet_call')

M1$Subnames_manual[colnames(ec)] <- unlist(ec[["Subnames_manual"]])
M1$scrublet_score[colnames(ec)] <- unlist(ec[["scrublet_score"]])
M1$scrublet_call[colnames(ec)] <- unlist(ec[["scrublet_call"]])



# Doublets


counts <- as(GetAssayData(M1, layer = "counts"), "dgCMatrix")
counts_t <- Matrix::t(counts)  # cells x genes for scrublet
counts_py <- scipy$csc_matrix(
  tuple(list(as.numeric(counts_t@x), as.integer(counts_t@i), as.integer(counts_t@p))),
  shape = tuple(as.integer(nrow(counts_t)), as.integer(ncol(counts_t)))
)

scrublet_obj <- scrub$Scrublet(counts_py)
result <- scrublet_obj$scrub_doublets()

M1[["scrublet_score"]] <- result[[1]]
M1[["scrublet_call"]] <- ifelse(result[[2]], "Doublet", "Singlet")

table(M1$scrublet_call)
FeaturePlot(M1, features = 'scrublet_score')
DimPlot(M1, group.by = 'scrublet_call')



saveRDS(M1,'./output/snRV_ref.rds')




