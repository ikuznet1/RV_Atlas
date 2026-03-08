library(Seurat)
library(singleCellTK)
library(reticulate)
library(patchwork)
library(dplyr)
library(hdf5r)
library(tidyverse)
library(harmony)
conda_create("scrublet")
conda_install(envname="scrublet", packages ="pip","git")
conda_install(envname='scrublet','git+https://github.com/AllonKleinLab/scrublet.git',pip=T)
conda_install(envname="scrublet", packages ="numpy=1.20.0")



#sctkPythonInstallConda(envname = "sctk-reticulate",conda = "auto",packages = c("scipy", "numpy", "astroid",# "six"),pipPackages = c("scrublet", "scanpy", "louvain", "leidenalg", "bbknn", "scanorama","anndata"))


use_miniconda("scrublet", required=T)
scrub <- import("scrublet")
np <- import("numpy")
ndd <- import("ndd")



###################LOADING FILES###################

files <- list.files("/Volumes/Extreme\ SSD/CellBender_Final")

for(i in files){
        #loc <- paste("/Volumes/Extreme\ Pro/Broad_snRNAseq/cellranger/",i,"/filtered_feature_bc_matrix",sep="")
	loc <- paste("/Volumes/Extreme\ SSD/CellBender_Final/",i,"/RV_",i,"_seurat.h5",sep="")
	print(loc)
	assign(paste0("RV_",i,".data"),Read10X_h5(loc))
}


#RV_1298 = CreateSeuratObject(counts=RV_1298.data,project="RVF",min.cells=3,min.features = 200)
RV_1343 = CreateSeuratObject(counts=RV_1343.data,project="RVF",min.cells=3,min.features = 200)
RV_1392 = CreateSeuratObject(counts=RV_1392.data,project="pRV",min.cells=3,min.features = 200)
RV_1467 = CreateSeuratObject(counts=RV_1467.data,project="RVF",min.cells=3,min.features = 200)
RV_1561 = CreateSeuratObject(counts=RV_1561.data,project="NF",min.cells=3,min.features = 200)
RV_1567 = CreateSeuratObject(counts=RV_1567.data,project="pRV",min.cells=3,min.features = 200)
RV_1618 = CreateSeuratObject(counts=RV_1618.data,project="pRV",min.cells=3,min.features = 200)
RV_1632 = CreateSeuratObject(counts=RV_1632.data,project="RVF",min.cells=3,min.features = 200)
RV_1681 = CreateSeuratObject(counts=RV_1681.data,project="NF",min.cells=3,min.features = 200)
RV_1691 = CreateSeuratObject(counts=RV_1691.data,project="NF",min.cells=3,min.features = 200)
RV_1692 = CreateSeuratObject(counts=RV_1692.data,project="pRV",min.cells=3,min.features = 200)
RV_1697 = CreateSeuratObject(counts=RV_1697.data,project="NF",min.cells=3,min.features = 200)

#RV_1298@meta.data$group = "RVF"
RV_1343@meta.data$group = "RVF"
RV_1392@meta.data$group = "pRV"
RV_1467@meta.data$group = "RVF"
RV_1561@meta.data$group = "NF"
RV_1567@meta.data$group = "pRV"
RV_1618@meta.data$group = "pRV"
RV_1632@meta.data$group = "RVF"
RV_1681@meta.data$group = "NF"
RV_1691@meta.data$group = "NF"
RV_1692@meta.data$group = "pRV"
RV_1697@meta.data$group = "NF"

#Exclude 1298



ifnb.list = list(RV_1343, RV_1392, RV_1467, RV_1561, RV_1567, RV_1618, RV_1632, RV_1681, RV_1691, RV_1692, RV_1697)

names = list("1343", "1392", "1467", "1561", "1567", "1618", "1632", "1681", "1691", "1692", "1697")

###################Calculate exonic RNA and per-nuclei entropy###################

#use_python("/Users/ikuz/python-microscopy/bin/python")
ifnb.list  <- lapply(seq_along(ifnb.list), FUN = function(x) {
	fname = names[[x]]
	print(fname)
	loc <- paste("/Volumes/Extreme\ SSD/cellranger_final/",fname,"/",fname,"_scrinvex.tsv",sep="")
	out <- read.csv(loc,sep="\t", header=TRUE)
	data = out %>% group_by(barcode) %>% summarise(introns=sum(introns),exons=sum(exons),junctions=sum(junctions))
	index_exons_in_seurat <- match(rownames(ifnb.list[[x]]@meta.data),data$barcode)
	ifnb.list[[x]][["percent.exon"]] <-
mapply(function(y) data$exons[[y]] / (data$exons[[y]] + data$introns[[y]] + data$junctions[[y]]), index_exons_in_seurat)
	entropy_estimate <- apply(ifnb.list[[x]]@assays$RNA@counts,2,function(x) ndd$entropy(x,k=35614))
	ifnb.list[[x]][["entropy"]] <- entropy_estimate
	ifnb.list[[x]]
})  


raw <- ifnb.list

###################Basic QC###################

#Clean data
#Plotting
temp <- lapply(seq_along(ifnb.list), FUN = function(x) {
	data <- ifnb.list[[x]]

	data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
	qc <- VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.exon","entropy"), ncol = 5)
	pdf(paste("/Volumes/Extreme\ SSD/Final_Analysis/QC_Metrics/",names[[x]],".pdf",sep=""))
	plot(qc)
	dev.off()
    x <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 100000 & percent.mt < 5)
})  

#Process
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
	x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-")
	x
#    x <- subset(x, subset = nFeature_RNA > 200 & nFeature_RNA < 100000 & percent.mt < 5 & entropy > 5 & percent.exon < 0.25)
})  


###################Identify doublets###################
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
	sce <- convertSeuratToSCE(x)
	sce <- runScrublet(sce)
	scores <- colData(sce)['scrublet_score']
	call <- colData(sce)['scrublet_call']
	singlet <- which(call[,1] == "Singlet")
	print(length(singlet))
	t <- format(Sys.time(), "%Y-%m-%d-%H-%M-%S-%Z")
	p <- plotScrubletResults(inSCE=sce,reducedDimName="scrublet_UMAP")
	pdf(paste("/Volumes/Extreme\ SSD/Final_Analysis/QC_Metrics/",t,".pdf",sep=""))
	plot(p)
	dev.off()

	x[["scrublet_score"]] <- list(scores$scrublet_score)
	x[["scrublet_call"]] <- list(call$scrublet_call)
	x
})


saveRDS(ifnb.list, file = "/Volumes/Extreme\ SSD/Final_Analysis/dataset_with_qc_corr.rds")



###################Per-sample clustering-based nuclei removal###################

# normalize and identify variable features for each dataset independently   
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# PCA
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
    x <- ScaleData(x, verbose = TRUE)
    x <- RunPCA(x, verbose = TRUE)
})


ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
    x <- FindNeighbors(x, dims = 1:50)
    x <- FindClusters(x, resolution = 2)
    x <- RunUMAP(x, dims = 1:50)
})

temp <- lapply(seq_along(ifnb.list), FUN = function(x) {
    data <- ifnb.list[[x]] 
    p <- VlnPlot(data, features = c("percent.mt","percent.exon","entropy"), ncol = 3)
    pdf(paste("/Volumes/Extreme\ SSD/Final_Analysis/QC_Metrics/cluster_qc_",names[[x]],".pdf",sep=""))
    plot(p)
    dev.off()
})


ifnb.list <- lapply(seq_along(ifnb.list), FUN = function(x) {
    data_all <- ifnb.list[[x]] 
    num_clusters <- max(as.numeric(as.character(
    data_all@meta.data$seurat_clusters)))

    for (i in 0:num_clusters) {
	    data <- subset(data_all, idents = i)
	    mito_quartile_third=quantile(unlist(data[["percent.mt"]]), prob=0.75) 
	    mito_quartile_first=quantile(unlist(data[["percent.mt"]]), prob=0.25) 
	    mito_iqr = mito_quartile_third - mito_quartile_first 
	    mito_cutoff = 1.5 * mito_iqr +  mito_quartile_third

	    exon_quartile_third=quantile(unlist(data[["percent.exon"]]), prob=0.75) 
	    exon_quartile_first=quantile(unlist(data[["percent.exon"]]), prob=0.25) 
	    exon_iqr = exon_quartile_third - exon_quartile_first 
	    exon_cutoff = 1.5 * exon_iqr +  exon_quartile_third

	    entropy_quartile_third=quantile(unlist(data[["entropy"]]), prob=0.75) 
	    entropy_quartile_first=quantile(unlist(data[["entropy"]]), prob=0.25) 
	    entropy_iqr = exon_quartile_third - exon_quartile_first 
	    entropy_cutoff = exon_quartile_first - 1.5 * entropy_iqr
	    data_all <- tryCatch({
	    	cells <- WhichCells(data,expression = percent.mt > mito_cutoff | entropy < 	entropy_cutoff | percent.exon > exon_cutoff)
	    	print(paste("Deleting",length(cells),"from cluster",i,"in",names[[x]]))
	    	#print(cells)
        	data_all <- subset(data_all,cells=cells,invert=TRUE)
	    	},
	    	error=function(cond) {
	    	print(paste("Deleting 0 from cluster",i,"in",names[[x]]))
	    	data_all <- data_all
	    	return(data_all)
	    	})
        
    }
    data_all
})

temp <- lapply(seq_along(ifnb.list), FUN = function(x) {
    data <- ifnb.list[[x]] 
    p <- VlnPlot(data, features = c("percent.mt","percent.exon","entropy"), ncol = 3)
    pdf(paste("/Volumes/Extreme\ SSD/Final_Analysis/QC_Metrics/cluster_qc_",names[[x]],"_post.pdf",sep=""))
    plot(p)
    dev.off()
})

#Step 5 ends here
saveRDS(ifnb.list, file = "/Volumes/Extreme\ SSD/Final_Analysis/dataset_post_cluster_qc.rds")

ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
	x<-BuildClusterTree(x)
})


skl <- import("sklearn.covariance")


ifnb.list <- lapply(seq_along(ifnb.list), FUN = function(x) {
    data_all <- ifnb.list[[x]] 
    num_clusters <- max(as.numeric(as.character(
    data_all@meta.data$seurat_clusters)))

    for (i in 0:num_clusters) {
    	print(length(colnames(data_all)))
	    data <- subset(data_all, idents = i)
	    if (length(colnames(data)) < 3){break}
	    print(length(colnames(data)))
	    mat <- as.matrix(data@meta.data[c("percent.exon","entropy","percent.mt")])
	  cov <- skl$EllipticEnvelope(contamination = 0.05)$fit(mat)
	  pred <- cov$predict(mat)	  
	    	data_all <- tryCatch({
	    		print(paste("Deleting",sum(pred == -1),"from cluster",i,"in",names[[x]]))
        	data_all <- subset(data_all,cells=colnames(data)[pred == -1],invert=TRUE)
	    	},
	    	error=function(cond) {
	    		print(paste("Deleting 0 from cluster",i,"in",names[[x]]))
	    		data_all <- data_all
	    		return(data_all)
	    	})
        
    }
    data_all
})


saveRDS(ifnb.list, file = "/Volumes/Extreme\ SSD/Final_Analysis/dataset_post_ellipsoid_qc.rds")

ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
	x
    x <- subset(x, subset = nFeature_RNA > 200 & nFeature_RNA < 100000 & percent.mt < 1 & entropy > 5 & percent.exon < 0.25)
})  

temp <- lapply(seq_along(ifnb.list), FUN = function(x) {
    data <- ifnb.list[[x]] 
    p <- VlnPlot(data, features = c("percent.mt","percent.exon","entropy"), ncol = 3)
    pdf(paste("/Volumes/Extreme\ SSD/Final_Analysis/QC_Metrics/final_pre_doublets_",names[[x]],".pdf",sep=""))
    plot(p)
    dev.off()
})

saveRDS(ifnb.list, file = "/Volumes/Extreme\ SSD/Final_Analysis/dataset_post_clipping_qc.rds")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("harmony", version = "3.16")
library(scCustomize)
immune.combined <- Merge_Seurat_List(ifnb.list,names)

immune.combined<- NormalizeData(immune.combined, verbose = F)
immune.combined <- FindVariableFeatures(immune.combined, selection.method = "vst", nfeatures = 2000, verbose = F)
immune.combined <- ScaleData(immune.combined, verbose = F)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = F)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30, verbose = F)
DimPlot(immune.combined, reduction = "umap")

immune.combined[['patient']] <- sapply(str_split(rownames(immune.combined@meta.data),'_'), `[`, 1)
immune.combined <- RunHarmony(immune.combined,'patient')
DimPlot(object = immune.combined, reduction = "harmony", pt.size = .1, group.by = "patient") + NoLegend()


immune.combined <- immune.combined%>% 
  RunUMAP(reduction = "harmony", dims = 1:30, verbose = F) %>% 
  FindNeighbors(reduction = "harmony", k.param = 10, dims = 1:30) %>% 
  FindClusters() %>% 
  identity()
  

  
immune.combined <- SetIdent(immune.combined,value = "patient")

DimPlot(immune.combined,reduction = "umap")
saveRDS(immune.combined, file = "/Volumes/Extreme\ SSD/Final_Analysis/harmony.rds")