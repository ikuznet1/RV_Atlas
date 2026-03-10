library(Seurat)


M1<-readRDS(file = "./output/BroadData_Harmony_SCTransform_with_Doublets.rds")

M1 <- M1%>% 
  RunUMAP(reduction = "harmony", dims = 1:40, verbose = F) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:40) %>% 
  FindClusters(resolution=1) %>% 
  identity()
  
#Extract CM clusters
cm <- subset(M1,ident=c("0","4","9","12","21","25","26"))


cm <- cm%>% 
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(dims=1:80)%>%
  RunUMAP(reduction = "harmony", dims=1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()
DotPlot(cm,features=c("PLIN1","RYR2","LEPR","VWF","HAS1","DCN","CCL21","CSF1R","NRXN1","CD2","PDGFRB","MYH11"),col.min=0,col.max=2)+xlab('Marker Score')+xlab('Cluster')
#11,13 has VWF
#7,10,13 has DCN
#14 has PDGFRB

#Throw out 10,11,13,14
cm <- subset(cm,ident=c('10','11','13','14'),invert=T)
cm <- cm%>% 
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(dims=1:80)%>%
  RunUMAP(reduction = "harmony", dims=1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()
  
cm7 <- subset(cm,ident=c('7'))
cm7 <- cm7%>% 
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(dims=1:80)%>%
  RunUMAP(reduction = "harmony", dims=1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()

to_remove <- colnames(subset(cm7,ident=c("0","2","3","5")))

cm <- subset(cm,cells=to_remove,invert=T)
cm <- cm%>% 
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(dims=1:80)%>%
  RunUMAP(reduction = "harmony", dims=1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()

saveRDS(cm,file = "./output/CellTypes/cm.rds")
write.table(colnames(cm), file="./output/CellTypes/cm.tsv", sep="\t", quote=F, col.names=NA)

#Extract FB clusters
fb <- subset(M1,ident=c("1","2","11","22"))
fb <- fb%>% 
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(dims=1:80)%>%
  RunUMAP(reduction = "harmony", dims=1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()
DotPlot(fb,features=c("PLIN1","RYR2","LEPR","VWF","HAS1","DCN","CCL21","CSF1R","NRXN1","CD2","PDGFRB","MYH11"),col.min=0,col.max=2)+xlab('Marker Score')+xlab('Cluster')

#9 is RYR2
fb <- subset(fb,ident=c('9'),invert=T)
fb <- fb%>% 
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(dims=1:80)%>%
  RunUMAP(reduction = "harmony", dims=1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()
  
fb4 <- subset(fb,ident=c('4'))
fb4 <- fb4 %>% 
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(dims=1:80)%>%
  RunUMAP(reduction = "harmony", dims=1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()
to_remove <- colnames(subset(fb4,ident=c("1","2","4")))

fb <- subset(fb,cells=to_remove,invert=T)
fb <- fb%>% 
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(dims=1:80)%>%
  RunUMAP(reduction = "harmony", dims=1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()
  
saveRDS(fb,file = "./output/CellTypes/fb.rds")
write.table(colnames(fb), file="./output/CellTypes/fb.tsv", sep="\t", quote=F, col.names=NA)
  
write.table(c(colnames(cm),colnames(fb)), file="./output/CellTypes/all.tsv", sep="\t", quote=F, col.names=NA)

#Extract EC clusters
ec <- subset(M1,ident=c("8","13","15","19","24"))
ec <- ec%>% 
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(dims=1:80)%>%
  RunUMAP(reduction = "harmony", dims=1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()
DotPlot(ec,features=c("PLIN1","RYR2","LEPR","VWF","HAS1","DCN","CCL21","CSF1R","NRXN1","CD2","PDGFRB","MYH11"),col.min=0,col.max=2)+xlab('Marker Score')+xlab('Cluster')

#0,8 with a lot of DCN
#7,9 with a lot of RYR2

ec <- subset(ec,ident=c("0","8","7","9"),invert=T)

ec <- ec %>% 
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(dims=1:80)%>%
  RunUMAP(reduction = "harmony", dims=1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()
  
ec <- subset(ec,ident=c("9"),invert=T)
ec6 <- subset(ec,ident=c("6"))
subset(ec6,ident=c('0','1'))


ec6<-ec6 %>% 
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(dims=1:80)%>%
  RunUMAP(reduction = "harmony", dims=1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()
to_remove <- colnames(subset(ec6,ident=c('0','1')))

ec <- subset(ec,cells=to_remove,invert=T)
ec <- ec%>% 
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(dims=1:80)%>%
  RunUMAP(reduction = "harmony", dims=1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()

saveRDS(fb,file = "./output/CellTypes/ec.rds")
write.table(colnames(ec), file="./output/CellTypes/ec.tsv", sep="\t", quote=F, col.names=NA)
  
write.table(c(colnames(cm),colnames(fb),colnames(ec)), file="./output/CellTypes/all.tsv", sep="\t", quote=F, col.names=NA)
  
#########Had to remove a lot of EC cells

#Extract Endocardial clusters
endo <- subset(M1,ident=c("5"))
endo <- endo%>% 
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(dims=1:80)%>%
  RunUMAP(reduction = "harmony", dims=1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()
  
#7,9 are distant probably doublet

endo <- subset(endo,ident=c("7","9"),invert=T)
endo0 <- subset(endo,ident=c("0"))


endo0<-endo0 %>% 
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(dims=1:80)%>%
  RunUMAP(reduction = "harmony", dims=1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()
to_remove <- colnames(subset(endo0,ident=c('3')))
endo <- subset(endo,cells=to_remove,invert=T)

endo <- endo%>% 
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(dims=1:80)%>%
  RunUMAP(reduction = "harmony", dims=1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()
#2 with a lot of CM
#0 with a lot of DCN
endo2 <- subset(endo,ident="2")
endo2<-endo2 %>% 
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(dims=1:80)%>%
  RunUMAP(reduction = "harmony", dims=1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()
to_remove <- colnames(subset(endo2,ident=c('3')))
endo <- subset(endo,cells=to_remove,invert=T)
endo <- endo%>% 
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(dims=1:80)%>%
  RunUMAP(reduction = "harmony", dims=1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()
saveRDS(endo,file = "./output/CellTypes/endo.rds")
write.table(colnames(endo), file="./output/CellTypes/endo.tsv", sep="\t", quote=F, col.names=NA)
  
write.table(c(colnames(cm),colnames(fb),colnames(ec),colnames(endo)), file="./output/CellTypes/all.tsv", sep="\t", quote=F, col.names=NA)


#Extract adipo clusters
adipo <- subset(M1,ident=c("7"))
adipo <- adipo%>% 
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(dims=1:80)%>%
  RunUMAP(reduction = "harmony", dims=1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()
#0 with a lot of DCN
#8 with a lot of RYR2
#11,12 look like junk
#10 looks like junk
  
adipo <- subset(adipo,ident=c("11","12","10"),invert=T)
adipo <- adipo%>% 
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(dims=1:80)%>%
  RunUMAP(reduction = "harmony", dims=1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()
adipo6789 <- subset(adipo,ident=c("6","7","8","9"))
adipo6789 <- adipo6789%>% 
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(dims=1:80)%>%
  RunUMAP(reduction = "harmony", dims=1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()
#Remove 1,5
to_remove <- colnames(subset(adipo6789,ident=c("5","1")))
adipo <- subset(adipo,cells=to_remove,invert=T)
adipo <- adipo%>% 
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(dims=1:80)%>%
  RunUMAP(reduction = "harmony", dims=1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()
  
saveRDS(adipo,file = "./output/CellTypes/adipo.rds")
write.table(colnames(adipo), file="./output/CellTypes/adipo.tsv", sep="\t", quote=F, col.names=NA)
  
write.table(c(colnames(cm),colnames(fb),colnames(ec),colnames(endo),colnames(adipo)), file="./output/CellTypes/all.tsv", sep="\t", quote=F, col.names=NA)

#Extract myeloid clusters
myeloid <- subset(M1,ident=c("3","16","17","20","29"))
myeloid <- myeloid%>% 
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(dims=1:80)%>%
  RunUMAP(reduction = "harmony", dims=1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()
#11,10,9 is bad, maybe 7,8
#4 is all EC
#0 is DCN rich
myeloid_subclust <- subset(myeloid,ident=c("0","4","7","8","9","10","11"))
myeloid_subclust <- myeloid_subclust %>% 
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(dims=1:80)%>%
  RunUMAP(reduction = "harmony", dims=1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()
#Toss 1,11
to_remove <- colnames(subset(myeloid_subclust,ident=c("11","1")))
myeloid <- subset(myeloid,cells=to_remove,invert=T)
myeloid <- myeloid%>% 
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(dims=1:80)%>%
  RunUMAP(reduction = "harmony", dims=1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()
#10,11,12 looks bad, maybe 9
myeloid<-subset(myeloid,ident=c("10","11","12"),invert=T)
myeloid <- myeloid%>% 
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(dims=1:80)%>%
  RunUMAP(reduction = "harmony", dims=1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()
  
myeloid_subclust <- subset(myeloid,ident=c("2"))
myeloid_subclust <- myeloid_subclust %>% 
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(dims=1:80)%>%
  RunUMAP(reduction = "harmony", dims=1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()
#Toss2
to_remove <- colnames(subset(myeloid_subclust,ident=c("2")))
myeloid <- subset(myeloid,cells=to_remove,invert=T)
myeloid <- myeloid%>% 
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(dims=1:80)%>%
  RunUMAP(reduction = "harmony", dims=1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()
  
saveRDS(myeloid,file = "./output/CellTypes/myeloid.rds")
write.table(colnames(myeloid), file="./output/CellTypes/myeloid.tsv", sep="\t", quote=F, col.names=NA)
  
all<-read.table("./output/CellTypes/all.tsv", sep="\t")
write.table(c(colnames(cm),colnames(fb),colnames(ec),colnames(endo),colnames(adipo),colnames(myeloid)), file="./output/CellTypes/all.tsv", sep="\t", quote=F, col.names=NA)

#1,2,3,5,6,7,9,10 are all res
#0,4 are other
#0 is ITGAX (CD11c)+ and FLT3+
#4 might be a fibro/residnet mix
#9 might be a doublet pop - RGS5+/PDGFRB+ (pericyte) <-- remove?

#0 is a mix of DC and monos
#8 is Mac2 - nonresident

#Extract SM clusters
sm <- subset(M1,ident=c("10"))
sm <- sm%>% 
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(dims=1:80)%>%
  RunUMAP(reduction = "harmony", dims=1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()

#7 is VWF, 6 is CM, 5 is DCN
sm567<-subset(sm,ident=c('5','6','7'))
sm567 <- sm567%>% 
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(dims=1:80)%>%
  RunUMAP(reduction = "harmony", dims=1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()
to_remove<-colnames(subset(sm567,ident=c("0","3")))
sm <- subset(sm,cells=to_remove,invert=T)
sm <- sm%>% 
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(dims=1:80)%>%
  RunUMAP(reduction = "harmony", dims=1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()
sm<-subset(sm,ident=c('8'),invert=T)
sm7<-subset(sm,ident=c('7'))

sm7 <- sm7%>% 
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(dims=1:80)%>%
  RunUMAP(reduction = "harmony", dims=1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()
to_remove<-colnames(subset(sm7,ident=c("0")))
sm <- subset(sm,cells=to_remove,invert=T)
sm <- sm%>% 
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(dims=1:80)%>%
  RunUMAP(reduction = "harmony", dims=1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()
sm<-subset(sm,ident=c("8"),invert=T)
sm <- sm%>% 
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(dims=1:80)%>%
  RunUMAP(reduction = "harmony", dims=1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()
  
saveRDS(sm,file = "./output/CellTypes/sm.rds")
write.table(colnames(sm), file="./output/CellTypes/sm.tsv", sep="\t", quote=F, col.names=NA)
  
all<-read.table("./output/CellTypes/all.tsv", sep="\t")
write.table(c(all[,2],colnames(sm)), file="./output/CellTypes/all.tsv", sep="\t", quote=F, col.names=NA)

#Extract neuron clusters
neuro <- subset(M1,ident=c("23"))
neuro <- neuro%>% 
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(dims=1:80)%>%
  RunUMAP(reduction = "harmony", dims=1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()
  
#4,5,6,7
neuro <- subset(neuro,ident=c('4','5','6','7'),invert=T)
neuro <- neuro%>% 
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(dims=1:80)%>%
  RunUMAP(reduction = "harmony", dims=1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity() 
neuro <- subset(neuro,ident=c('3'),invert=T)
neuro <- neuro%>% 
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(dims=1:80)%>%
  RunUMAP(reduction = "harmony", dims=1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity() 
saveRDS(neuro,file = "./output/CellTypes/neuro.rds")
write.table(colnames(neuro), file="./output/CellTypes/neuro.tsv", sep="\t", quote=F, col.names=NA)
  
all<-read.table("./output/CellTypes/all.tsv", sep="\t")
write.table(c(all[,2],colnames(neuro)), file="./output/CellTypes/all.tsv", sep="\t", quote=F, col.names=NA)

#Extract epicardial clusters
epi <- subset(M1,ident=c("30"))
epi <- epi%>% 
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(dims=1:80)%>%
  RunUMAP(reduction = "harmony", dims=1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity() 
epi0 <- subset(epi,ident=0)
epi0 <- epi0%>% 
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(dims=1:80)%>%
  RunUMAP(reduction = "harmony", dims=1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity() 
to_remove <- colnames(subset(epi0,ident=1))
epi <- subset(epi,cells=to_remove,invert=T)
epi <- epi%>% 
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(dims=1:80)%>%
  RunUMAP(reduction = "harmony", dims=1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity() 
saveRDS(epi,file = "./output/CellTypes/epi.rds")
write.table(colnames(epi), file="./output/CellTypes/epi.tsv", sep="\t", quote=F, col.names=NA)
  
all<-read.table("./output/CellTypes/all.tsv", sep="\t")
write.table(c(all[,2],colnames(epi)), file="./output/CellTypes/all.tsv", sep="\t", quote=F, col.names=NA)  
  
#Extract ccl21 clusters
ccl <- subset(M1,ident=c("14"))
ccl <- ccl%>% 
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(dims=1:80)%>%
  RunUMAP(reduction = "harmony", dims=1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity() 
ccl <- subset(ccl,ident=c('8','9'),invert=T)
ccl <- ccl%>% 
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(dims=1:80)%>%
  RunUMAP(reduction = "harmony", dims=1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity() 
ccl <- subset(ccl,ident=c('6','7'),invert=T)
ccl <- ccl%>% 
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(dims=1:80)%>%
  RunUMAP(reduction = "harmony", dims=1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity() 
#ccl0<-subset(ccl,ident=0)
#ccl0 <- ccl0%>% 
#  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
#  RunPCA(dims=1:80)%>%
#  RunUMAP(reduction = "harmony", dims=1:50, verbose = F) %>% 
#  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
#  FindClusters(resolution=1) %>% 
#  identity() 
#to_remove <- colnames(subset(ccl0,ident=1))
#epi <- subset(epi,cells=to_remove,invert=T)
#epi <- epi%>% 
#  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
#  RunPCA(dims=1:80)%>%
#  RunUMAP(reduction = "harmony", dims=1:50, verbose = F) %>% 
#  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
#  FindClusters(resolution=1) %>% 
#  identity() 
saveRDS(ccl,file = "./output/CellTypes/lec.rds")
write.table(colnames(ccl), file="./output/CellTypes/lec.tsv", sep="\t", quote=F, col.names=NA)
  
all<-read.table("./output/CellTypes/all.tsv", sep="\t")
write.table(c(all[,2],colnames(ccl)), file="./output/CellTypes/all.tsv", sep="\t", quote=F, col.names=NA)  

#Extract pericyte clusters
pc <- subset(M1,ident=c("6","27"))
pc <- pc%>% 
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(dims=1:80)%>%
  RunUMAP(reduction = "harmony", dims=1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity() 
#1,2,3,5,6 are suspicious  
pc12356 <- subset(pc,ident=c('1','2','3','5','6'))
pc12356 <- pc12356%>% 
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(dims=1:80)%>%
  RunUMAP(reduction = "harmony", dims=1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity() 
to_remove <- colnames(subset(pc12356,ident=c("0","7","3")))
pc <- subset(pc,cells=to_remove,invert=T)  
pc <- pc%>% 
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(dims=1:80)%>%
  RunUMAP(reduction = "harmony", dims=1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity() 

pc8910 <- subset(pc,ident=c('8',"9","10"))

pc8910 <- pc8910%>% 
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(dims=1:80)%>%
  RunUMAP(reduction = "harmony", dims=1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity() 
to_remove <- colnames(subset(pc12356,ident=c("3")))
pc <- subset(pc,cells=to_remove,invert=T)  
pc <- pc%>% 
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(dims=1:80)%>%
  RunUMAP(reduction = "harmony", dims=1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity() 
saveRDS(pc,file = "./output/CellTypes/pc.rds")
write.table(colnames(pc), file="./output/CellTypes/pc.tsv", sep="\t", quote=F, col.names=NA)
  
all<-read.table("./output/CellTypes/all.tsv", sep="\t")
write.table(c(all[,2],colnames(pc)), file="./output/CellTypes/all.tsv", sep="\t", quote=F, col.names=NA)  

#Extract NK/T clusters
nkt <- subset(M1,ident=c("18"))
nkt <- nkt%>% 
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(dims=1:80)%>%
  RunUMAP(reduction = "harmony", dims=1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity() 
nkt <- subset(nkt,ident=c("5","6","8"),invert=T)
nkt <- nkt%>% 
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(dims=1:80)%>%
  RunUMAP(reduction = "harmony", dims=1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity() 
nkt <- subset(nkt,ident=c("6"),invert=T)
nkt <- nkt%>% 
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(dims=1:80)%>%
  RunUMAP(reduction = "harmony", dims=1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity() 
saveRDS(nkt,file = "./output/CellTypes/nkt.rds")
write.table(colnames(nkt), file="./output/CellTypes/nkt.tsv", sep="\t", quote=F, col.names=NA)
  
all<-read.table("./output/CellTypes/all.tsv", sep="\t")
write.table(c(all[,2],colnames(nkt)), file="./output/CellTypes/all.tsv", sep="\t", quote=F, col.names=NA)  

#Prolif cells
prolif <- subset(M1,ident=c("28"))
prolif <- prolif%>% 
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(dims=1:80)%>%
  RunUMAP(reduction = "harmony", dims=1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity() 
prolif <- subset(prolif,ident=c("3"),invert=T)
prolif <- prolif%>% 
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(dims=1:80)%>%
  RunUMAP(reduction = "harmony", dims=1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity() 
prolif12 <- subset(prolif,ident=c('1','2'))
prolif12 <- prolif12%>% 
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(dims=1:80)%>%
  RunUMAP(reduction = "harmony", dims=1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity() 
to_remove <- subset(prolif12,ident=c('2','3')) 
  
prolif <- subset(prolif,cells=to_remove,invert=T)  
prolif <- prolif%>% 
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(dims=1:80)%>%
  RunUMAP(reduction = "harmony", dims=1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity() 
prolif <- subset(prolif,ident=1,invert=T)
prolif <- prolif%>% 
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(dims=1:80)%>%
  RunUMAP(reduction = "harmony", dims=1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity() 
  
saveRDS(prolif,file = "./output/CellTypes/prolif.rds")
write.table(colnames(prolif), file="./output/CellTypes/prolif.tsv", sep="\t", quote=F, col.names=NA)
  
all<-read.table("./output/CellTypes/all.tsv", sep="\t")
write.table(c(all[,2],colnames(prolif)), file="./output/CellTypes/all.tsv", sep="\t", quote=F, col.names=NA)  


###Re-clust
all<-read.table("./output/CellTypes/all.tsv", sep="\t")
all <- all[,2]

M1 <- subset(M1,cells=all)  

M1 <- M1%>% 
  RunUMAP(reduction = "harmony", dims = 1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()

saveRDS(M1,file = "./output/CellTypes/dataset_post_r1.rds")
#Need to clean up VWF from T Cells
#Extra pericyte cluster

##############################Second round
library('harmony')  

#Extract CM clusters
cm <- subset(M1,ident=c("0","2","11","12","19","23"))
cm <- cm%>%  
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()
#DotPlot(cm,features=c("PLIN1","RYR2","LEPR","VWF","HAS1","DCN","CCL21","CSF1R","NRXN1","CD2","PDGFRB","MYH11"),col.min=0,col.max=2)+xlab('Marker Score')+xlab('Cluster')

cm <- cm%>%
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(assay="SCT",npcs=50)%>%
  RunUMAP(reduction = "harmony", dims = 1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()
  

 
saveRDS(cm,file = "./output/CellTypes/cm.rds")
write.table(colnames(cm), file="./output/CellTypes/cm.tsv", sep="\t", quote=F, col.names=NA)
write.table(colnames(cm), file="./output/CellTypes/all.tsv", sep="\t", quote=F, col.names=NA)  




#Extract FB clusters
fb <- subset(M1,ident=c("1","4","15","18"))
fb <- fb%>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()

p1<-DimPlot(fb,label=T)
cells.located.1 <- CellSelector(plot = p1)
cells.located.2 <- CellSelector(plot = p1)
cells.located.3 <- CellSelector(plot = p1)
to_remove <- c(cells.located.1,cells.located.2,cells.located.3)
fb <- subset(fb,cell=to_remove,invert=T)

fb <- fb%>%
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(assay="SCT",npcs=50)%>%
  RunUMAP(reduction = "harmony", dims = 1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()
  

 
saveRDS(fb,file = "./output/CellTypes/fb.rds")
write.table(colnames(fb), file="./output/CellTypes/fb.tsv", sep="\t", quote=F, col.names=NA)
  
all<-read.table("./output/CellTypes/all.tsv", sep="\t")
write.table(c(all[,2],colnames(fb)), file="./output/CellTypes/all.tsv", sep="\t", quote=F, col.names=NA)  

#Extract EC clusters
ec <- subset(M1,ident=c("7","16","17"))
ec <- ec%>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()

p1<-DimPlot(ec,label=T)
cells.located.1 <- CellSelector(plot = p1)
cells.located.2 <- CellSelector(plot = p1)
cells.located.3 <- CellSelector(plot = p1)
to_remove <- c(cells.located.1,cells.located.2,cells.located.3)
ec <- subset(ec,cell=to_remove,invert=T)

ec <- ec%>%
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(assay="SCT",npcs=50)%>%
  RunUMAP(reduction = "harmony", dims = 1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()
  

 
saveRDS(ec,file = "./output/CellTypes/ec.rds")
write.table(colnames(ec), file="./output/CellTypes/ec.tsv", sep="\t", quote=F, col.names=NA)
  
all<-read.table("./output/CellTypes/all.tsv", sep="\t")
write.table(c(all[,2],colnames(ec)), file="./output/CellTypes/all.tsv", sep="\t", quote=F, col.names=NA)  


#Extract Endocardial clusters
endo <- subset(M1,ident=c("5"))

endo <- endo%>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()

p1<-DimPlot(endo,label=T)
cells.located.1 <- CellSelector(plot = p1)
cells.located.2 <- CellSelector(plot = p1)
cells.located.3 <- CellSelector(plot = p1)
to_remove <- c(cells.located.1,cells.located.2,cells.located.3)
endo <- subset(endo,cell=to_remove,invert=T)

endo <- endo%>%
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(assay="SCT",npcs=50)%>%
  RunUMAP(reduction = "harmony", dims = 1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()
  

 
saveRDS(endo,file = "./output/CellTypes/endo.rds")
write.table(colnames(endo), file="./output/CellTypes/endo.tsv", sep="\t", quote=F, col.names=NA)
  
all<-read.table("./output/CellTypes/all.tsv", sep="\t")
write.table(c(all[,2],colnames(endo)), file="./output/CellTypes/all.tsv", sep="\t", quote=F, col.names=NA)  



#Extract adipo clusters
adipo <- subset(M1,ident=c("6"))

adipo <- adipo %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()

p1<-DimPlot(adipo,label=T)
cells.located.1 <- CellSelector(plot = p1)
to_remove <- c(cells.located.1)
adipo <- subset(adipo,cell=to_remove,invert=T)

adipo <- adipo%>%
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(assay="SCT",npcs=50)%>%
  RunUMAP(reduction = "harmony", dims = 1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()
 
saveRDS(adipo,file = "./output/CellTypes/adipo.rds")
write.table(colnames(adipo), file="./output/CellTypes/adipo.tsv", sep="\t", quote=F, col.names=NA)
  
all<-read.table("./output/CellTypes/all.tsv", sep="\t")
write.table(c(all[,2],colnames(adipo)), file="./output/CellTypes/all.tsv", sep="\t", quote=F, col.names=NA)  


#Extract myeloid clusters
myeloid <- subset(M1,ident=c("3","10","14"))
myeloid <- myeloid %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()

p1<-DimPlot(myeloid,label=T)
cells.located.1 <- CellSelector(plot = p1)
cells.located.2 <- CellSelector(plot = p1)
cells.located.3 <- CellSelector(plot = p1)

to_remove <- c(cells.located.1,cells.located.2,cells.located.3)
myeloid <- subset(myeloid,cell=to_remove,invert=T)

myeloid <- myeloid%>%
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(assay="SCT",npcs=50)%>%
  RunUMAP(reduction = "harmony", dims = 1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()
  
myeloid5 <- subset(myeloid,ident=5)

myeloid5 <- myeloid5%>%
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(assay="SCT",npcs=50)%>%
  RunUMAP(reduction = "harmony", dims = 1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution = 1) %>% 
  identity()
  
to_remove <- colnames(subset(myeloid5,ident=c("3","4")))
myeloid <- subset(myeloid,cell=to_remove,invert=T)
myeloid <- myeloid%>%
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(assay="SCT",npcs=50)%>%
  RunUMAP(reduction = "harmony", dims = 1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()
 
saveRDS(myeloid,file = "./output/CellTypes/myeloid.rds")
write.table(colnames(myeloid), file="./output/CellTypes/myeloid.tsv", sep="\t", quote=F, col.names=NA)
  
all<-read.table("./output/CellTypes/all.tsv", sep="\t")
write.table(c(all[,2],colnames(myeloid)), file="./output/CellTypes/all.tsv", sep="\t", quote=F, col.names=NA)  


#Extract SM clusters
sm <- subset(M1,ident=c("9"))
  
sm <- sm %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()

p1<-DimPlot(sm,label=T)
cells.located.1 <- CellSelector(plot = p1)


to_remove <- c(cells.located.1)
sm <- subset(sm,cell=to_remove,invert=T)

sm <-sm%>%
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(assay="SCT",npcs=50)%>%
  RunUMAP(reduction = "harmony", dims = 1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()
   
saveRDS(sm,file = "./output/CellTypes/sm.rds")
write.table(colnames(sm), file="./output/CellTypes/sm.tsv", sep="\t", quote=F, col.names=NA)
  
all<-read.table("./output/CellTypes/all.tsv", sep="\t")
write.table(c(all[,2],colnames(sm)), file="./output/CellTypes/all.tsv", sep="\t", quote=F, col.names=NA)  




#Neuro clusts
neuro <-subset(M1,ident=c("22"))
  
neuro <- neuro %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()

p1<-DimPlot(neuro,label=T)
cells.located.1 <- CellSelector(plot = p1)
cells.located.2 <- CellSelector(plot = p1)

to_remove <- c(cells.located.1,cells.located.2)
neuro <- subset(neuro,cell=to_remove,invert=T)

neuro <-neuro%>%
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(assay="SCT",npcs=50)%>%
  RunUMAP(reduction = "harmony", dims = 1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()
   
saveRDS(neuro,file = "./output/CellTypes/neuro.rds")
write.table(colnames(neuro), file="./output/CellTypes/neuro.tsv", sep="\t", quote=F, col.names=NA)
  
all<-read.table("./output/CellTypes/all.tsv", sep="\t")
write.table(c(all[,2],colnames(neuro)), file="./output/CellTypes/all.tsv", sep="\t", quote=F, col.names=NA)  


#Extract epicardial clusters
epi <- subset(M1,ident=c("25"))

epi <- epi %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()

p1<-DimPlot(neuro,label=T)
#cells.located.1 <- CellSelector(plot = p1)
#cells.located.2 <- CellSelector(plot = p1)

#to_remove <- c(cells.located.1,cells.located.2)
#neuro <- subset(neuro,cell=to_remove,invert=T)

epi <-epi%>%
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(assay="SCT",npcs=50)%>%
  RunUMAP(reduction = "harmony", dims = 1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()
   
saveRDS(epi,file = "./output/CellTypes/epi.rds")
write.table(colnames(epi), file="./output/CellTypes/epi.tsv", sep="\t", quote=F, col.names=NA)
  
all<-read.table("./output/CellTypes/all.tsv", sep="\t")
write.table(c(all[,2],colnames(epi)), file="./output/CellTypes/all.tsv", sep="\t", quote=F, col.names=NA)  
  
#Extract ccl21 clusters
ccl <- subset(M1,ident=c("13"))

ccl <- ccl %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()

p1<-DimPlot(ccl,label=T)
cells.located.1 <- CellSelector(plot = p1)
cells.located.2 <- CellSelector(plot = p1)

to_remove <- c(cells.located.1,cells.located.2)
ccl <- subset(ccl,cell=to_remove,invert=T)

ccl <-ccl%>%
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(assay="SCT",npcs=50)%>%
  RunUMAP(reduction = "harmony", dims = 1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()
   
saveRDS(ccl,file = "./output/CellTypes/ccl.rds")
write.table(colnames(ccl), file="./output/CellTypes/ccl.tsv", sep="\t", quote=F, col.names=NA)
  
all<-read.table("./output/CellTypes/all.tsv", sep="\t")
write.table(c(all[,2],colnames(ccl)), file="./output/CellTypes/all.tsv", sep="\t", quote=F, col.names=NA)  

#Extract pericyte clusters
pc <- subset(M1,ident=c("8","20"))
pc <- pc %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()

p1<-DimPlot(pc,label=T)
cells.located.1 <- CellSelector(plot = p1)
cells.located.2 <- CellSelector(plot = p1)

to_remove <- c(cells.located.1,cells.located.2)
pc <- subset(pc,cell=to_remove,invert=T)

pc <-pc%>%
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(assay="SCT",npcs=50)%>%
  RunUMAP(reduction = "harmony", dims = 1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()
   
saveRDS(pc,file = "./output/CellTypes/pc.rds")
write.table(colnames(pc), file="./output/CellTypes/pc.tsv", sep="\t", quote=F, col.names=NA)
  
all<-read.table("./output/CellTypes/all.tsv", sep="\t")
write.table(c(all[,2],colnames(pc)), file="./output/CellTypes/all.tsv", sep="\t", quote=F, col.names=NA) 


#Extract NK/T clusters
nkt <- subset(M1,ident=c("21"))

nkt <- nkt %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()

p1<-DimPlot(nkt,label=T)
cells.located.1 <- CellSelector(plot = p1)
cells.located.2 <- CellSelector(plot = p1)
cells.located.3 <- CellSelector(plot = p1)


to_remove <- c(cells.located.1,cells.located.2,cells.located.3)
nkt <- subset(nkt,cell=to_remove,invert=T)

nkt <-nkt%>%
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(assay="SCT",npcs=50)%>%
  RunUMAP(reduction = "harmony", dims = 1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution = 1) %>% 
  identity()
   
saveRDS(nkt,file = "./output/CellTypes/nkt.rds")
write.table(colnames(nkt), file="./output/CellTypes/nkt.tsv", sep="\t", quote=F, col.names=NA)
  
all<-read.table("./output/CellTypes/all.tsv", sep="\t")
write.table(c(all[,2],colnames(nkt)), file="./output/CellTypes/all.tsv", sep="\t", quote=F, col.names=NA) 


#Prolif cells
prolif <- subset(M1,ident=c("24"))

prolif <- prolif %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()

p1<-DimPlot(prolif,label=T)
cells.located.1 <- CellSelector(plot = p1)


to_remove <- c(cells.located.1)
prolif <- subset(prolif,cell=to_remove,invert=T)

prolif <-prolif%>%
  SCTransform(verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')%>%
  RunPCA(assay="SCT",npcs=50)%>%
  RunUMAP(reduction = "harmony", dims = 1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:50) %>% 
  FindClusters(resolution = 1) %>% 
  identity()
   
saveRDS(prolif,file = "./output/CellTypes/prolif.rds")
write.table(colnames(prolif), file="./output/CellTypes/prolif.tsv", sep="\t", quote=F, col.names=NA)
  
all<-read.table("./output/CellTypes/all.tsv", sep="\t")
write.table(c(all[,2],colnames(prolif)), file="./output/CellTypes/all.tsv", sep="\t", quote=F, col.names=NA) 


####Recluster

###Re-clust
all<-read.table("./output/CellTypes/all.tsv", sep="\t")
all <- all[,2]

M1 <- subset(M1,cells=all)  

M1 <- M1%>% 
  RunUMAP(reduction = "harmony", dims = 1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:50) %>% 
  FindClusters(resolution=0.5) %>% 
  identity()
  
p1<-DimPlot(M1,label=T)
cells.located.1 <- CellSelector(plot = p1)
cells.located.2 <- CellSelector(plot = p1)
cells.located.3 <- CellSelector(plot = p1)
to_remove <- c(cells.located.1,cells.located.2,cells.located.3)
M1 <- subset(M1,cell=to_remove,invert=T)
M1 <- M1%>% 
  RunUMAP(reduction = "harmony", dims = 1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:50) %>% 
  FindClusters(resolution=0.5) %>% 
  identity()
  
saveRDS(M1,file = "./output/CellTypes/integrated_no_doublets.rds")
M1 <- PrepSCTFindMarkers(M1)
cm.markers <- FindMarkers(M1,ident.1=c('0','8','15'))
write.table(cm.markers, file="./output/CellTypes/cm_genes.tsv", sep="\t", quote=F, col.names=NA)
