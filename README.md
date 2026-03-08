# RV_Atlas

**Human Adult and Pediatric Single-Nucleus Transcriptomic Atlas of Progression from Pressure Loaded to Failing Right Ventricle**

*Kuznetsov IA, Li K, Guedira Y, Simonson B, Chaffin M, Bedi KC Jr., Thome T, Zhao W, Zhu W, Zhou W, Yang Y, Kadyrov F, Amrute JM, Lai L, Griffin J, Li L, Li J, Miyamoto SD, Ellinor P, Margulies KB, Lavine KJ, Arany Z\*, Edwards JJ\**

---

## About

This repository contains the R analysis code used to generate all figures in the manuscript. The paper characterizes the transcriptional landscape of healthy, pressure-overloaded, and failing human right ventricles (RVs) using a multi-modal approach: bulk RNA-seq (n=142), single-nucleus RNA-seq (snRNA-seq; n=11 adult, n=14 pediatric), and spatial transcriptomics (10X Xenium). Key findings include:

- **Cardiomyocytes** downregulate nuclear-encoded mitochondrial transcripts and show reduced mitochondrial respiration in RV failure (RVF)
- **Myeloid cells** upregulate MHCII-associated genes, indicating a shift toward antigen presentation and a pro-inflammatory state in RVF
- **Endothelial cells** expand in RVF, driven by capillary and arterial subtypes in adults and venous subtypes in pediatric hearts — an expansion not seen in left ventricular failure
- A murine pulmonary artery-banding (PAB) model of RVF recapitulates EC expansion but diverges from human RVF in myeloid and cardiomyocyte transcriptional programs, cautioning against its uncritical use
- Pediatric RVF (hypoplastic left heart syndrome) largely mirrors adult RVF transcriptionally, with notable differences in mitochondrial and endothelial programs

## Repository purpose

This repository enables reproduction of all publication figures (`Figure_1.R` – `Figure_8.R`) and supplementary figures (`Supplementary_Figure_1.R` – `Supplementary_Figure_8.r`). Each script is self-contained and reads processed data objects from the `./dependencies/` directory.

---

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



