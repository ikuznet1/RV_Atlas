# RV_Atlas
Scripts for generating figures for RV atlas paper


### DEPENDENCIES ###

First create a conda enviornment:

conda create -n RV_atlas -c conda-forge -c bioconda r-base=4.4 mamba
conda activate RV_atlas

mamba install -c conda-forge -c bioconda r-seurat r-hdf5r r-igraph r-tidyverse r-ggraph r-harmony r-enrichr r-devtools 
mamba install -c conda-forge -c bioconda bioconductor-ucell bioconductor-genomicranges bioconductor-geneoverlap 
mamba install -c conda-forge r-cairo
mamba install -c conda-forge r-sf


Then open R and run:

install.packages("BiocManager")
BiocManager::install()

BiocManager::install("WGCNA")
devtools::install_github('smorabit/hdWGCNA', ref='dev')

install.packages('arrow')
install.packages('reticulate')
BiocManager::install('EnsDb.Hsapiens.v86')
BiocManager::install("edgeR")
BiocManager::install('sva')
BiocManager::install('DESeq2')
BiocManager::install('scCustomize')
remotes::install_github('satijalab/seurat-wrappers')

BiocManager::install('tximport')
BiocManager::install('tximportData')



BiocManager::install("biomaRt")
install.packages('colormap')
install.packages('cowplot')
BiocManager::install("DESeq2")
BiocManager::install("DOSE")
BiocManager::install('EnhancedVolcano')
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
BiocManager::install('Nebulosa')
install.packages('matrixStats')
install.packages('ashr')
devtools::install_github('immunogenomics/presto')
install.packages("ClusterR")
install.packages('VennDiagram')
BiocManager::install('glmGamPoi')



devtools::install_github('cole-trapnell-lab/monocle3')



