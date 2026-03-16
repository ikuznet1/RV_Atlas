##########################
##### Myeloid labels #####
##########################

M1 <- readRDS(file = "./dependencies/shared/myeloid_subclust.rds")

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



df <- data.frame(colnames(M1),M1$Subsubnames)

write.csv(df,'./dependencies/shared/myeloid_subclustering.csv')


#####################
##### EC labels #####
#####################

M1 <- readRDS(file = "./dependencies/shared/EC_hdWGCNA_by_celltype.rds")

df <- data.frame(colnames(M1),M1$Names)

write.csv(df,'./dependencies/shared/ec_subclustering.csv')

#####################
##### CM labels #####
#####################

M1 <- readRDS('./dependencies/shared/cm_new_subclust.rds')
df <- data.frame(colnames(M1),M1@active.ident)
write.csv(df,'./dependencies/shared/cm_subclustering.csv')


#####################
##### FB labels #####
#####################

M1 <- readRDS(file = "./dependencies/shared/fb_subclust.rds")

new.cluster.ids <- c("Fb1","Fb2","Fb3","Fb4","Fb5","Fb6","Fb7")
names(new.cluster.ids) <- levels(M1)
M1 <- RenameIdents(M1, new.cluster.ids)

df <- data.frame(colnames(M1),M1@active.ident)

write.csv(df,'./dependencies/shared/fb_subclustering.csv')


########################
##### Mural labels #####
########################


M1 <- readRDS(file = "./dependencies/Supplementary_Figure_4/pc_sm_subclust.rds")


new.cluster.ids <- c("Sm","Pc","Pc","Sm")

names(new.cluster.ids) <- levels(M1)
M1 <- RenameIdents(M1, new.cluster.ids)


df <- data.frame(colnames(M1),M1@active.ident)

write.csv(df,'./dependencies/shared/mural_subclustering.csv')


##################
##### Assign #####
##################

M1 <- readRDS('~/Documents/XeniumWorkflow/snRV_ref.rds')

fb <- read.csv('./dependencies/shared/fb_subclustering.csv',header=T,row.names = 1)
cm <- read.csv('./dependencies/shared/cm_subclustering.csv',header=T,row.names = 1)
ec <- read.csv('./dependencies/shared/ec_subclustering.csv',header=T,row.names = 1)
myeloid <- read.csv('./dependencies/shared/myeloid_subclustering.csv',header=T,row.names = 1)
mural <- read.csv('./dependencies/shared/mural_subclustering.csv',header=T,row.names = 1)


M1$Subnames <- paste0(M1$Names,'_unc')

M1$Subnames[rownames(fb)] <- fb$'M1.active.ident'
M1$Subnames[rownames(cm)] <- cm$'M1.active.ident'
M1$Subnames[rownames(myeloid)] <- myeloid$'M1.Subsubnames'
M1$Subnames[rownames(ec)] <- paste0('EC_',ec$'M1.Names')

M1$Subnames[M1$Subnames == 'LEC_unc'] = 'EC_Lymph'
M1$Subnames[M1$Subnames == 'Endo_unc'] = 'EC_Endocardial'
M1$Subnames[M1$Subnames == 'NKT_unc'] = 'NK_T'
M1$Subnames[M1$Subnames == 'Neuron_unc'] = 'Neuron'
M1$Subnames[M1$Subnames == 'PC_unc'] = 'PC'
M1$Subnames[M1$Subnames == 'SM_unc'] = 'SM'
M1$Subnames[M1$Subnames == 'Adipo_unc'] = 'Adipo'
M1$Subnames[M1$Subnames == 'Epi_unc'] = 'Epi'

saveRDS(M1,'~/Documents/XeniumWorkflow/snRV_ref.rds')


