library(Seurat)
library(harmony)
library(sctransform)

ifnb.list <- readRDS('./output/dataset_post_clipping_qc.rds')

ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
	x<-SCTransform(x,verbose=FALSE,conserve.memory=TRUE,vars.to.regress='percent.mt')
})

ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
	x<-RunPCA(x,npcs=80)
	x<-RunUMAP(x,dims=1:80,verbose=FALSE)
})

names = list("1343", "1392", "1467", "1561", "1567", "1618", "1632", "1681", "1691", "1692", "1697")

saveRDS(ifnb.list, file = "./output/dataset_pre_harmony_after_sctransform.rds")
var.feats <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 2000)

ifnb.list[[1]]$orig.ident<- names[[1]]
ifnb.list[[2]]$orig.ident<- names[[2]]
ifnb.list[[3]]$orig.ident<- names[[3]]
ifnb.list[[4]]$orig.ident<- names[[4]]
ifnb.list[[5]]$orig.ident<- names[[5]]
ifnb.list[[6]]$orig.ident<- names[[6]]
ifnb.list[[7]]$orig.ident<- names[[7]]
ifnb.list[[8]]$orig.ident<- names[[8]]
ifnb.list[[9]]$orig.ident<- names[[9]]
ifnb.list[[10]]$orig.ident<- names[[10]]
ifnb.list[[11]]$orig.ident<- names[[11]]

data1<-RenameCells(ifnb.list[[1]],add.cell.id = names[[1]])
data2<-RenameCells(ifnb.list[[2]],add.cell.id = names[[2]])
data3<-RenameCells(ifnb.list[[3]],add.cell.id = names[[3]])
data4<-RenameCells(ifnb.list[[4]],add.cell.id = names[[4]])
data5<-RenameCells(ifnb.list[[5]],add.cell.id = names[[5]])
data6<-RenameCells(ifnb.list[[6]],add.cell.id = names[[6]])
data7<-RenameCells(ifnb.list[[7]],add.cell.id = names[[7]])
data8<-RenameCells(ifnb.list[[8]],add.cell.id = names[[8]])
data9<-RenameCells(ifnb.list[[9]],add.cell.id = names[[9]])
data10<-RenameCells(ifnb.list[[10]],add.cell.id = names[[10]])
data11<-RenameCells(ifnb.list[[11]],add.cell.id = names[[11]])
rm('ifnb.list')

M1<- merge(data1,data2)
rm('data1')
rm('data2')
gc()
M2 <- merge(data3,data4)
rm('data3')
rm('data4')
gc()
M1<- merge(M1,M2)
M2 <- merge(data5,data6)
rm('data5')
rm('data6')
gc()
M1<- merge(M1,M2)
M2 <- merge(data7,data8)
rm('data7')
rm('data8')
gc()
M1<- merge(M1,M2)
rm('M2')
gc()
M1<- merge(M1,data9)
rm('data9')
gc()
M1<- merge(M1,data10)
rm('data10')
gc()
M1<- merge(M1,data11)
rm('data11')
gc()
saveRDS(M1, file = "./output/BroadData_Merge_SCTransform_with_Doublets.rds")


M1[['patient']] <- sapply(strsplit(rownames(M1@meta.data),'_'), `[`, 1)
VariableFeatures(M1[["SCT"]]) <- var.feats
M1<-RunPCA(M1)
M1 <- RunHarmony(M1,'patient')
saveRDS(M1, file = "./output/BroadData_Harmony_SCTransform_with_Doublets.rds")

M1 <- M1%>% 
  RunUMAP(reduction = "harmony", dims = 1:40, verbose = F) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:40) %>% 
  FindClusters() %>% 
  identity()
