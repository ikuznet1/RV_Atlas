# RV_Atlas
Scripts for generating figures for RV atlas paper


### DEPENDENCIES ###

First create a conda enviornment:

conda create -n RV_atlas -c conda-forge -c bioconda r-base=4.4 mamba
conda activate RV_atlas

mamba install -c conda-forge -c bioconda r-seurat r-hdf5r r-igraph r-tidyverse r-ggraph r-harmony r-enrichr r-devtools r-cairo r-sf
mamba install -c conda-forge -c bioconda bioconductor-ucell bioconductor-genomicranges bioconductor-geneoverlap 


Then open R and run:

install.packages("BiocManager")
BiocManager::install()
BiocManager::install("WGCNA")
BiocManager::install('EnsDb.Hsapiens.v86')
BiocManager::install("edgeR")
BiocManager::install('sva')
BiocManager::install('DESeq2')
BiocManager::install('scCustomize')
BiocManager::install('tximport')
BiocManager::install('tximportData')
BiocManager::install("biomaRt")
BiocManager::install("DESeq2")
BiocManager::install("DOSE")
BiocManager::install('EnhancedVolcano')
BiocManager::install('Nebulosa')
BiocManager::install('glmGamPoi')

install.packages('arrow')
install.packages('reticulate')
install.packages('colormap')
install.packages('cowplot')
install.packages('forcats')
install.packages('ggeasy')
install.packages('ggfortify')
install.packages('gplots')
install.packages('viridis')
install.packages('stringr')
install.packages('ggpubr')
install.packages('reshape2')
install.packages('readxl')
install.packages('RColorBrewer')
install.packages('pracma')
install.packages('patchwork')
install.packages('matrixStats')
install.packages('ashr')
install.packages("ClusterR")
install.packages('VennDiagram')

devtools::install_github('smorabit/hdWGCNA', ref='dev')
remotes::install_github('satijalab/seurat-wrappers')
devtools::install_github('immunogenomics/presto')
devtools::install_github('cole-trapnell-lab/monocle3')



