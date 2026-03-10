#Take doublet removed data and decontaminate

library(singleCellTK)
library(Seurat)
library(celda)
library(SeuratWrappers)
library(monocle3)



M1<- readRDS(file = "./output/CellTypes/integrated_no_doublets.rds")

new.cluster.ids <- c("CM","FB","Myeloid","EC","Endo","PC","Adipo","SM","CM","LEC","Myeloid","EC","NKT","FB","Neuron","CM","Myeloid","Epi")
names(new.cluster.ids) <- levels(M1)
M1 <- RenameIdents(M1, new.cluster.ids)
M1$Names <- M1@active.ident

DefaultAssay(M1) <- "RNA"

cds <- as.cell_data_set(M1)
cds$cluster <- M1$Names
cds <- decontX(cds,z=M1$Names,batch=M1$patient)
#plot_cells(cds[,colData(cds) %>% subset(patient == "1392") %>% row.names],color_cells_by="decontX_clusters",show_trajectory_graph=F)
#plot_cells(cds[,colData(cds) %>% subset(patient == "1392") %>% row.names],color_cells_by="decontX_contamination",show_trajectory_graph=F)

markers <- list(adipo = "PLIN1",cm = "RYR2",endo = "LEPR",ec = "VWF",epi = "HAS1",fb = "DCN",lec = "CCL21",myeloid = "CSF1R",neuron = "NRXN1",nkt = "CD2",pc = "PDGFRB",sm = "MYH11")

#plotDecontXMarkerPercentage(cds[,colData(cds) %>% subset(patient == "1392") %>% row.names],markers = markers,assayName = c("counts","decontXcounts"))

#library(scater)
#a <- cds[,colData(cds) %>% subset(patient == "1392") %>% row.names]
#a <- logNormCounts(a,exprs_values = "decontXcounts",name = "decontXlogcounts")

#plotDecontXMarkerExpression(a,markers = markers[["cm"]],ncol = 3)
#plotDecontXMarkerExpression(a,markers = markers[["fb"]],ncol = 3)
#plotDecontXMarkerExpression(a,markers = markers[["ec"]],ncol = 3)
#plotDecontXMarkerExpression(a,markers = markers[["nkt"]],ncol = 3)

#plotDecontXMarkerExpression(a,markers = markers[["cm"]],ncol = 3, assayName = c("logcounts", "decontXlogcounts"))
#plotDecontXMarkerExpression(cds[,colData(cds) %>% subset(patient == "1392") %>% row.names],markers = markers[["fb"]],ncol = 3, assayName = c("logcounts", "decontXlogcounts"))
#plotDecontXMarkerExpression(cds[,colData(cds) %>% subset(patient == "1392") %>% row.names],markers = markers[["ec"]],ncol = 3, assayName = c("logcounts", "decontXlogcounts"))
#plotDecontXMarkerExpression(cds[,colData(cds) %>% subset(patient == "1392") %>% row.names],markers = markers[["nkt"]],ncol = 3, assayName = c("logcounts", "decontXlogcounts"))

M1[["decontXcounts"]] <- CreateAssayObject(counts = decontXcounts(cds))
DefaultAssay(M1) <- "decontXcounts"

rm(cds)
gc()
M2 <- SplitObject(M1,split.by="patient")
rm(M1)
gc()

M2 <- lapply(X = M2, FUN = function(x) {
	x[['SCT']] <- NULL
	x
	#x<-SCTransform(x,assay='decontXcounts',vst.flavor="v2",verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')
})


saveRDS(M2,file = "./output/BroadData_Post_Decontam.rds")

M2 <- lapply(X = M2, FUN = function(x) {
	DefaultAssay(x) <- "SCT"
	x
})


#M2 <- lapply(X = M2, FUN = function(x) {
#	x[['decontXcounts']] <- NULL
#	x
#})


data1<-M2[[1]]
data2<-M2[[2]]
data3<-M2[[3]]
data4<-M2[[4]]
data5<-M2[[5]]
data6<-M2[[6]]
data7<-M2[[7]]
data8<-M2[[8]]
data9<-M2[[9]]
data10<-M2[[10]]
data11<-M2[[11]]
var.feats <- SelectIntegrationFeatures(object.list = M2, nfeatures = 2000)
rm(M2)
gc()

M1 <- merge(data1,data2)
M1 <- merge(M1,data3)
rm(data1)
rm(data2)
rm(data3)
gc()
M1 <- merge(M1,data4)
M1 <- merge(M1,data5)
M1 <- merge(M1,data6)
rm(data4)
rm(data5)
rm(data6)
gc()
M1 <- merge(M1,data7)
rm(data7)
gc()
M1 <- merge(M1,data8)
rm(data8)
gc()
M1 <- merge(M1,data9)
rm(data9)
gc()
M1 <- merge(M1,data10)
rm(data10)
gc()
M1 <- merge(M1,data11)
rm(data11)
gc()


saveRDS(M1,file = "./output/CellTypes/temp.rds")
VariableFeatures(M1[["SCT"]]) <- var.feats

library(harmony)
library(tidyverse)


M1<-RunPCA(M1)
M1 <- RunHarmony(M1,'patient')
M1 <- M1%>% 
  RunUMAP(reduction = "harmony", dims = 1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:50) %>% 
  FindClusters() %>% 
  identity()
saveRDS(M1,file = "./output/temp.rds")



##########
cm <- subset(M1,subset=Names=="CM")
cm <- cm%>%  
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()
DotPlot(cm,features=c("PLIN1","RYR2","LEPR","VWF","HAS1","DCN","CCL21","CSF1R","NRXN1","CD2","PDGFRB","MYH11"),col.min=0,col.max=2)+xlab('Marker Score')+xlab('Cluster')
cm <- subset(cm,ident=14,invert=T)

saveRDS(cm,file = "./output/CellTypes/cm.rds")
write.table(colnames(cm), file="./output/CellTypes/cm.tsv", sep="\t", quote=F, col.names=NA)
write.table(colnames(cm), file="./output/CellTypes/all.tsv", sep="\t", quote=F, col.names=NA)  

fb <- subset(M1,subset=Names=="FB")
fb <- fb%>%  
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()
DotPlot(cm,features=c("PLIN1","RYR2","LEPR","VWF","HAS1","DCN","CCL21","CSF1R","NRXN1","CD2","PDGFRB","MYH11"),col.min=0,col.max=2)+xlab('Marker Score')+xlab('Cluster')
fb <- subset(fb,ident=12,invert=T)
p1<-DimPlot(fb,label=T)
cells.located.1 <- CellSelector(plot = p1)
cells.located.2 <- CellSelector(plot = p1)
to_remove <- c(cells.located.1,cells.located.2)
fb <- subset(fb,cell=to_remove,invert=T)


saveRDS(fb,file = "./output/CellTypes/fb.rds")
write.table(colnames(fb), file="./output/CellTypes/fb.tsv", sep="\t", quote=F, col.names=NA)
all<-read.table("./output/CellTypes/all.tsv", sep="\t")
write.table(c(all[,2],colnames(fb)), file="./output/CellTypes/all.tsv", sep="\t", quote=F, col.names=NA)  

pc <- subset(M1,subset=Names=="PC")
pc <- pc%>%  
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()
p1<-DimPlot(pc,label=T)
cells.located.1 <- CellSelector(plot = p1)
cells.located.2 <- CellSelector(plot = p1)
cells.located.3 <- CellSelector(plot = p1)

to_remove <- c(cells.located.1,cells.located.2,cells.located.3)
pc <- subset(pc,cell=to_remove,invert=T)
pc <- subset(pc,ident=4,invert=T)

saveRDS(pc,file = "./output/CellTypes/pc.rds")
write.table(colnames(pc), file="./output/CellTypes/pc.tsv", sep="\t", quote=F, col.names=NA)
all<-read.table("./output/CellTypes/all.tsv", sep="\t")
write.table(c(all[,2],colnames(pc)), file="./output/CellTypes/all.tsv", sep="\t", quote=F, col.names=NA)  


myeloid <- subset(M1,subset=Names=="myeloid")
myeloid <- myeloid %>%  
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()
p1<-DimPlot(myeloid,label=T)
cells.located.1 <- CellSelector(plot = p1)
cells.located.2 <- CellSelector(plot = p1)
to_remove <- c(cells.located.1,cells.located.2)
myeloid <- subset(myeloid,cell=to_remove,invert=T)
myeloid <- subset(myeloid,ident=6,invert=T)

saveRDS(myeloid,file = "./output/CellTypes/myeloid.rds")
write.table(colnames(myeloid), file="./output/CellTypes/myeloid.tsv", sep="\t", quote=F, col.names=NA)
all<-read.table("./output/CellTypes/all.tsv", sep="\t")
write.table(c(all[,2],colnames(myeloid)), file="./output/CellTypes/all.tsv", sep="\t", quote=F, col.names=NA)  



nkt <- subset(M1,subset=Names=="NKT")
nkt <- nkt %>%  
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()
p1<-DimPlot(nkt,label=T)
cells.located.1 <- CellSelector(plot = p1)
cells.located.2 <- CellSelector(plot = p1)
cells.located.3 <- CellSelector(plot = p1)


to_remove <- c(cells.located.1,cells.located.2)
nkt <- subset(nkt,cell=to_remove,invert=T)

saveRDS(nkt,file = "./output/CellTypes/nkt.rds")
write.table(colnames(nkt), file="./output/CellTypes/nkt.tsv", sep="\t", quote=F, col.names=NA)
all<-read.table("./output/CellTypes/all.tsv", sep="\t")
write.table(c(all[,2],colnames(nkt)), file="./output/CellTypes/all.tsv", sep="\t", quote=F, col.names=NA)  


epi <- subset(M1,subset=Names=="Epi")
epi <- epi %>%  
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()
epi <- subset(epi,ident=3,invert=T)

saveRDS(epi,file = "./output/CellTypes/epi.rds")
write.table(colnames(epi), file="./output/CellTypes/epi.tsv", sep="\t", quote=F, col.names=NA)
all<-read.table("./output/CellTypes/all.tsv", sep="\t")
write.table(c(all[,2],colnames(epi)), file="./output/CellTypes/all.tsv", sep="\t", quote=F, col.names=NA)  

adipo <- subset(M1,subset=Names=="Adipo")
adipo <- adipo %>%  
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()
adipo <- subset(adipo,DCN>1.75,invert=T)

  
p1<-DimPlot(adipo,label=T)
cells.located.1 <- CellSelector(plot = p1)


to_remove <- c(cells.located.1)
adipo <- subset(adipo,cell=to_remove,invert=T)


saveRDS(adipo,file = "./output/CellTypes/adipo.rds")
write.table(colnames(adipo), file="./output/CellTypes/adipo.tsv", sep="\t", quote=F, col.names=NA)
all<-read.table("./output/CellTypes/all.tsv", sep="\t")
write.table(c(all[,2],colnames(adipo)), file="./output/CellTypes/all.tsv", sep="\t", quote=F, col.names=NA)  

ec <- subset(M1,subset=Names=="EC")
ec <- ec %>%  
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()
p1<-DimPlot(ec,label=T)
cells.located.1 <- CellSelector(plot = p1)
cells.located.2 <- CellSelector(plot = p1)
cells.located.3 <- CellSelector(plot = p1)


to_remove <- c(cells.located.1,cells.located.2,cells.located.3)
ec <- subset(ec,cell=to_remove,invert=T)

saveRDS(ec,file = "./output/CellTypes/ec")
write.table(colnames(ec), file="./output/CellTypes/ec", sep="\t", quote=F, col.names=NA)
all<-read.table("./output/CellTypes/all.tsv", sep="\t")
write.table(c(all[,2],colnames(ec)), file="./output/CellTypes/all.tsv", sep="\t", quote=F, col.names=NA)  

lec <- subset(M1,subset=Names=="LEC")
lec <- lec %>%  
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()
p1<-DimPlot(lec,label=T)
cells.located.1 <- CellSelector(plot = p1)
cells.located.2 <- CellSelector(plot = p1)
cells.located.3 <- CellSelector(plot = p1)


to_remove <- c(cells.located.1,cells.located.2,cells.located.3)
lec <- subset(lec,cell=to_remove,invert=T)

saveRDS(lec,file = "./output/CellTypes/lec")
write.table(colnames(lec), file="./output/CellTypes/lec", sep="\t", quote=F, col.names=NA)
all<-read.table("./output/CellTypes/all.tsv", sep="\t")
write.table(c(all[,2],colnames(lec)), file="./output/CellTypes/all.tsv", sep="\t", quote=F, col.names=NA)  


endo <- subset(M1,subset=Names=="Endo")
endo <- endo %>%  
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()
p1<-DimPlot(endo,label=T)
cells.located.1 <- CellSelector(plot = p1)
cells.located.2 <- CellSelector(plot = p1)
cells.located.3 <- CellSelector(plot = p1)


to_remove <- c(cells.located.1,cells.located.2,cells.located.3)
endo <- subset(endo,cell=to_remove,invert=T)

saveRDS(endo,file = "./output/CellTypes/endo")
write.table(colnames(endo), file="./output/CellTypes/endo", sep="\t", quote=F, col.names=NA)
all<-read.table("./output/CellTypes/all.tsv", sep="\t")
write.table(c(all[,2],colnames(endo)), file="./output/CellTypes/all.tsv", sep="\t", quote=F, col.names=NA)  

sm <- subset(M1,subset=Names=="SM")
sm <- sm %>%  
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()
p1<-DimPlot(sm,label=T)
cells.located.1 <- CellSelector(plot = p1)
cells.located.2 <- CellSelector(plot = p1)
cells.located.3 <- CellSelector(plot = p1)


to_remove <- c(cells.located.1,cells.located.2,cells.located.3)
sm <- subset(sm,cell=to_remove,invert=T)

saveRDS(sm,file = "./output/CellTypes/sm")
write.table(colnames(sm), file="./output/CellTypes/sm", sep="\t", quote=F, col.names=NA)
all<-read.table("./output/CellTypes/all.tsv", sep="\t")
write.table(c(all[,2],colnames(sm)), file="./output/CellTypes/all.tsv", sep="\t", quote=F, col.names=NA)  




neuro <- subset(M1,subset=Names=="Neuron")
neuro <- neuro %>%  
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()
p1<-DimPlot(neuro,label=T)
cells.located.1 <- CellSelector(plot = p1)
cells.located.2 <- CellSelector(plot = p1)
cells.located.3 <- CellSelector(plot = p1)


to_remove <- c(cells.located.1,cells.located.2,cells.located.3)
neuro <- subset(neuro,cell=to_remove,invert=T)

saveRDS(neuro,file = "./output/CellTypes/neuro.rds")
write.table(colnames(neuro), file="./output/CellTypes/neuro.tsv", sep="\t", quote=F, col.names=NA)
all<-read.table("./output/CellTypes/all.tsv", sep="\t")
write.table(c(all[,2],colnames(neuro)), file="./output/CellTypes/all.tsv", sep="\t", quote=F, col.names=NA)  


####SUBSET
all<-read.table("./output/CellTypes/all.tsv", sep="\t")
all <- unique(all[,2])
M1 <- subset(M1,cells=all)



M1 <- M1%>% 
  RunUMAP(reduction = "harmony", dims = 1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:50) %>% 
  FindClusters() %>% 
  identity()

saveRDS(M1,file = "./output/CellTypes/Post_R3.rds")

#Clean up

p1<-DimPlot(M1,label=T)
cells.located.1 <- CellSelector(plot = p1)
cells.located.2 <- CellSelector(plot = p1)
cells.located.3 <- CellSelector(plot = p1)
to_remove <- c(cells.located.1,cells.located.2,cells.located.3)
M1 <- subset(M1,cell=to_remove,invert=T)
M1 <- M1%>% 
  RunUMAP(reduction = "harmony", dims = 1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:50) %>% 
  FindClusters() %>% 
  identity()
  
p1<-DimPlot(M1,label=T)
cells.located.1 <- CellSelector(plot = p1)
to_remove <- c(cells.located.1)

M1 <- subset(M1,cell=to_remove,invert=T)

saveRDS(M1,file = "./output/Post_R3_FINAL.rds")


