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

This repository enables reproduction of all publication figures (`Figure_1.R` – `Figure_8.R`) and supplementary figures (`Supplementary_Figure_1.R` – `Supplementary_Figure_8.r`). Each script reads processed data objects from `./dependencies/`; see *Data bundle (Zenodo)* below for layout and the one cross-script dependency.

---

### DATA BUNDLE (Zenodo) ###

Figure scripts read processed data objects from `./dependencies/`. The full bundle (≈42 GB) is archived on Zenodo at **[DOI placeholder]**. Download and extract into the repository root so that `dependencies/` sits alongside the `Figure_*.R` scripts.

**Layout**

- `dependencies/shared/` — DEG CSVs, metadata, snRNA-seq + bulk Seurat objects, the cached module eigengenes matrix (`scWGCNA_bulk2sn_MEs.rds`), bulk-WGCNA consensus module table, and the Kallisto bulk RNA-seq outputs (`BulkRNA/`).
- `dependencies/shared/Xenium/` — Xenium spatial objects (`Xenium_resegmented_corrected.rds`, `Xenium_resegmented_imputed_final.rds`, `xenium_obj_subclustered.rds`), nine per-lineage `*_clean_clean.rds` subcluster objects, `cellnest_output/` (CellNEST cell–cell communication tables), Xenium pseudobulk DEG CSVs, and Xenium metadata.
- `dependencies/Figure_8/` — Pediatric + adult co-embedded Seurat objects for cross-cohort comparison. The `*_annotated.rds` files are scale.data-stripped and xz-compressed (~85% smaller than the working copies); standard `readRDS()` loads them transparently.

**Run order**

Most scripts are independent and can run in any order. Two cases warrant attention:

1. `Figure_2.R` and `Figure_5.R` rebuild a small module-eigengenes cache on first run, then reuse it. The cache is shipped pre-built as `dependencies/shared/scWGCNA_bulk2sn_MEs.rds`, so first-run cost is zero on a fresh Zenodo download.
2. **`Supplementary_Figure_4.R` panels E/F+ (the legacy myeloid block, gated by `RUN_LEGACY_MYELOID=TRUE`) read `./output/Figure_6/myeloid_subclust_new.rds`.** If you want to render that gated block, run `Figure_6.R` first so it produces that file. The default run path (without the env var) skips the gated block.

**Known TODO panels** — placeholders in the current bundle:

- `Figure_8.R` panels G/H/I require `dependencies/shared/Koenig_LV_with_RV_modules.rds` (not bundled).
- `Figure_8.R` panels L/M require `dependencies/shared/rv_dcm_merged.rds` (not bundled).

The script renders informative placeholders in these positions when the files are absent.

**Optional re-derivation (not needed for figure reproduction)**

- `additional_scripts/` contains pipelines for re-deriving bundled processed objects from raw data: `XeniumReanalysis.r` (Xenium resegmentation + subclustering), `snRNAReanalysis.r` (snRNA re-clustering), and `AnalysisPAH.R` (PAH bulk DEG generation, which also writes the `pah.*` / `rv.*` CSVs referenced by `Supplementary_Figure_1.R`). These scripts reference paths outside `dependencies/` and are not required to reproduce the published figures.
- `scWGCNA_bulk2sn_projection.rds` (the full 20 GB Seurat object with bulk-to-snRNA WGCNA projection) is **not bundled**. `Figure_2.R` and `Figure_5.R` read the precomputed MEs cache (`scWGCNA_bulk2sn_MEs.rds`, 13 MB) instead. To re-derive the projection from scratch, build the bulk WGCNA via `Figure_1.R` and then run `hdWGCNA::ProjectModules()` against the snRNA Seurat object (`RV_data.rds`).
- Raw sequencing reads (bulk + snRNA-seq + Xenium) are deposited in GEO under accession **[GSE placeholder]**.

---

### PREPROCESSING PIPELINE ###

Raw CellBender-corrected H5 files are processed through a four-step pipeline in `preprocess_pipeline/` before figure generation. Steps must be run in order; each saves intermediate `.rds` files to `./output/`.

**STEP1_Preprocess.R** — Per-sample QC, doublet scoring, and initial clustering
Reads CellBender H5 files from `./dependencies/CellBender_Final/`. For each of the 11 adult samples (1343, 1392, 1467, 1561, 1567, 1618, 1632, 1681, 1691, 1692, 1697) it:
- Assigns disease group labels: NF (non-failing), pRV (pressure-overloaded RV), or RVF (RV failure)
- Computes per-nuclei QC metrics: mitochondrial fraction (`percent.mt`), exonic fraction (`percent.exon`), and transcriptional entropy using `scrinvex` output and the `ndd` Python library
- Scores doublets using Scrublet (via `reticulate`; requires the `scrublet` conda environment)
- Clusters each sample individually and removes low-quality nuclei using IQR-based outlier thresholds on `percent.mt`, `percent.exon`, and entropy — applied per cluster
- Applies a second outlier-removal pass using sklearn's `EllipticEnvelope` (contamination = 0.05) on the three QC metrics jointly
- Applies hard cutoffs: `nFeature_RNA > 200`, `percent.mt < 1`, `entropy > 5`, `percent.exon < 0.25`
- Outputs: `./output/dataset_post_clipping_qc.rds`

**STEP0_Integration.R** — Cross-sample integration with Harmony
Reads `./output/dataset_post_clipping_qc.rds`. For each sample applies SCTransform (regressing out `percent.mt`), PCA (80 PCs), and UMAP. Merges all 11 samples, selects 2,000 integration features, and runs Harmony batch correction on patient ID. Re-embeds with UMAP and Leiden clustering.
Outputs: `./output/BroadData_Harmony_SCTransform_with_Doublets.rds`

**STEP2_Doublet_removal.R** — Cell-type isolation and doublet purification
Reads the Harmony-integrated object. Iteratively extracts each of 12 cell types (CM, FB, EC, Endocardial, Adipo, Myeloid, SM, Neuron, Epicardial, LEC, PC, NKT, Proliferating) by their cluster identities, re-clusters each subset with SCTransform + Harmony, and removes contaminating nuclei identified by canonical marker expression (e.g., VWF for EC contamination in CMs, DCN for fibroblast contamination). A second round of refinement is performed after merging all clean cells. The process runs for three rounds total, with manual cell selection (`CellSelector`) used for fine-grained cleanup in round 2.
Outputs per cell type to `./output/CellTypes/` and final cleaned object to `./output/CellTypes/integrated_no_doublets.rds`

**STEP3_Decontaminate_ambient.R** — Ambient RNA removal and final integration
Reads `./output/CellTypes/integrated_no_doublets.rds`. Converts to SingleCellExperiment and runs `decontX` (from the `celda` package) with cell-type and patient-batch labels to model and remove ambient RNA contamination. The decontaminated count matrix is stored as a new `decontXcounts` assay. Samples are split, re-merged, and re-integrated with SCTransform + Harmony. A final three-round cell-type re-extraction and cleanup produces the publication-ready object.
Final output: `./output/Post_R3_FINAL.rds` (used as `./dependencies/shared/RV_data.rds` by figure scripts)

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


### VISUALIZATION ###

To quickly visualize single-nuclei gene-expression first make sure you have the dependencies folder downloaded.  

Then open R and run:  

shiny::runApp("shinyViz")

![PLIN1 UMAP](shinyViz/UMAP_PLIN1.png)


To quickly visualize Xenium gene-expression first make sure you have the dependencies folder downloaded.  

Then open R and run:  

shiny::runApp("XeniumExp")

![PLIN1 UMAP](XeniumExp/Xenium_P1697_niche_manual.png)

---

### FIGURES (v57 manuscript) ###

#### Main Figures

**Figure 1 — Study overview and bulk RNA-sequencing of human RVF** (`Figure_1.R`)
Cohort schematic across three platforms (bulk n=142, snRNA n=11, Xenium n=9), PCA of bulk transcriptomes, WGCNA UMAP of 29 modules with top hub genes and GO enrichments, module-score violins across NF/pRV/RVF for significant modules, and cross-cohort concordance dotplots between RVF and pulmonary arterial hypertension (PAH).

**Figure 2 — Single-nucleus transcriptomic atlas of the right ventricle** (`Figure_2.R`)
UMAP of 61,398 nuclei colored by 12 lineage clusters (CM, FB, EC, Endo, SM, Myeloid, NKT, PC, Adipo, Epi, LEC, Neuron) with patient-overlay inset for integration QC, stacked-violin canonical lineage markers, per-lineage frequency boxplots across disease states (capturing progressive EC expansion and CM contraction), and a two-row dotplot of bulk WGCNA module scores by lineage (i) and disease group (ii).

**Figure 3 — Xenium spatial transcriptomic atlas of the right ventricle** (`Figure_3.R`)
UMAP of 625,305 Xenium cells across 12 lineages with patient-overlay inset, canonical-marker dotplot, cross-platform pseudobulk log₂FC concordance scatter (snRNA-seq vs Xenium per lineage × contrast), spatial maps of all NF/pRV/RVF sections by lineage and by 14-niche assignment (BuildNicheAssay, k=100, CLR-transformed), and niche frequency boxplots.

**Figure 4 — Two-phase transcriptional trajectory of RV failure** (`Figure_4.R`)
Two-phase model of NF→pRV (Phase 1) and pRV→RVF (Phase 2) transitions integrated across bulk, snRNA-seq pseudobulk, and Xenium pseudobulk. Lineage-level key-gene heatmaps per phase, bulk WGCNA module × phase × platform × lineage heatmap, Phase 1 vs Phase 2 DEG asymmetry quantification, CellNEST ligand-receptor pairs across disease states, and cell-type communication chord diagrams.

**Figure 5 — RV cardiomyocyte transcriptional and metabolic remodeling** (`Figure_5.R`)
CM-enriched WGCNA module dotplot by subcluster × disease state, GO enrichments on intersected module + DEG gene lists, M2 module volcano (pRV vs NF), mitochondrial-module volcanoes (M10, M25, M26, M28), CollecTRI/decoupleR TF activity for mitochondrial regulators, snRNA-seq pseudobulk violins of *ESRRA*/*ESRRG*/*PPARA*/*PPARGC1A*, Oroboros 2k high-resolution respirometry on isolated RV mitochondria (adult NF/pRV/RVF and pediatric NF/HLHS-palliated; PAB mouse 2-week + 2-month), and a cross-cohort OXPHOS summary panel showing reduced maximal coupled respiration in adult RVF and murine PAB but preserved in pediatric HLHS-palliated.

**Figure 6 — Myeloid subclustering and disease-associated programs** (`Figure_6.R`)
Six myeloid subtypes: *CCR2*⁻ resident macrophage (rMac), inflammatory macrophages (iMac), Mono / Mono-derived, *TREM2*⁺ macrophages, dendritic cells, and proliferating myeloid. Marker dotplot, disease-state frequencies, four myeloid-relevant bulk WGCNA module scores (M1, M3, M4, M8), ChEA TF enrichments for up/down DEGs, and pseudobulk violins for the four transcriptional programs (GR-homeostatic and HIF/vascular drift in *CCR2*⁻ rMac; NF-κB MHCII/inflammasome in iMac; IFNγ antigen presentation pooled). Spatial localization of myeloid subtypes in the full Xenium dataset.

**Figure 7 — Endothelial subclustering and disease-associated programs** (`Figure_7.R`)
Five EC subtypes (Arterial, Capillary, Venous, Lymph, Endocardial); per-subtype frequency boxplots; hdWGCNA on EC snRNA-seq with seven modules; GO enrichments for capillary-core (ecM1), arterial (ecM4) and venous/inflamed (ecM7) modules; EC module-score dotplot and FeaturePlots; Xenium pan-EC pseudobulk volcano (RVF vs pRV) with Phase-2 candidate annotation; cross-platform concordance scatter for Phase-2 candidates; per-subtype expression of *NRG1* (endocardial) and *MECOM* (arterial).

**Figure 8 — Cross-species, cross-age, and cross-etiology comparisons** (`Figure_8.R`)
Co-embedding of HLHS (n=209,604 nuclei) and adult RV (n=61,398 nuclei) datasets; HLHS cell-type abundance by disease state; bulk WGCNA module expression by cell type and disease state across both cohorts; cross-cohort log₂FC scatter (adult-RVF vs ped-RVF); HLHS CM-specific module expression by disease state (M10/M26 rise with severity; M25/M28 flat); GO enrichments for pooled mitochondrial-module DEGs; LV (Koenig 2022 DCM) cross-comparison panels for M2 and pooled mitochondrial modules at both per-nuclei and per-subject pseudobulk resolution; cross-ventricle log₂FC scatter (LV DCM/NF vs RV RVF/NF).

---

#### Supplementary Figures

**Figure S1 — Bulk RNA-seq supplemental analyses** (`Supplementary_Figure_1.R`)
RVF vs NF volcano; NF-vs-RVF vs pRV-vs-NF fold-change scatter; GO biological-process enrichments for pRV-up and RVF-up genes; heatmaps of monotonically increasing and decreasing gene sets across NF/pRV/RVF; GO enrichments for stepwise-decreasing genes (M1/M25/M26, mitochondrial translation machinery); progressively-downregulated DEG volcanoes (NF vs pRV; NF vs RVF) highlighting mitochondrial ribosome components; PAH-comparison divergent patterns in *CD163* and mitochondrial-translation genes.

**Figure S2 — Cardiomyocyte subclustering and spatial validation** (`Supplementary_Figure_2.R`)
UMAP of 10 CM subpopulations with patient-overlay inset; marker dotplot; per-state cluster frequency; per-cluster GO enrichments; Xenium spatial distribution of CM subclusters; per-cluster WGCNA module-score dotplot; MitoCarta3.0 mitochondrial-gene dotplot by CM subtype × disease; *HAND2* and *EDNRA* pseudobulk violins; WGA-stained minimum-Feret diameter quantification of CM hypertrophy (NF/pRV/RVF; cell-level linear mixed model); per-patient Xenium *NPPA*/*NPPB* rasterized spatial tiles.

**Figure S3 — Fibroblast and mural cell subclustering and spatial validation** (`Supplementary_Figure_3.R`)
Seven FB subpopulations (UMAP, marker dotplot, per-patient frequencies); FB Phase 1 sn-vs-Xenium concordance heatmap and Phase 2 concordance scatter at FB lineage; "erosion of FB identity" Phase-1 (Koenig 2022 donor-FB signature) and "matrifibrocyte commitment" Phase-2 (Fu 2018 mature-scar signature) per-patient module scores; Xenium FB subtype dot plot with SpaGE-imputed values; Xenium↔snRNA FB cluster similarity heatmap via label transfer; FB–myeloid spatial colocalization (K=15 nearest-neighbor analysis) with per-patient co-abundance scatter and spatial vignettes; adult RV trichrome/Sirius-red fibrosis quantification; mural-cell subclustering (PC and SM) with Phase 1 and Phase 2 heatmaps, vascular-protective Phase-1 and vascular IEG/stress Phase-2 module scores, and sn-vs-Xenium concordance scatter.

**Figure S4 — Left ventricular failure (DCM) comparison** (`Supplementary_Figure_4.R`)
Integrated human RV fibroblast UMAP with reference-mapped LV; FB subtype concordance between ventricles; Koenig (LV DCM) per-patient pseudobulk module scores for FB Phase 1 (donor-FB identity erosion) and Phase 2 (matrifibrocyte commitment); LV-FB vs RV-FB log₂FC scatter; integrated RV myeloid UMAP; reference mapping of LV snRNA (i) and scRNA (ii) myeloid onto RV with original LV annotations (iii, iv); cross-chamber concordance scatter for the four myeloid bulk-WGCNA modules (M1, M3, M4, M8); LV (Koenig 2022) myeloid dot plot of *AddModuleScore* values for the four modules by subtype × Donor/DCM; per-patient pseudobulk module-score box plots for the four RV-derived myeloid programs (GR-homeostatic, HIF/vascular drift, NF-κB MHCII/inflammasome, IFNγ antigen presentation) in LV myeloid (Donor vs DCM).

**Figure S5 — Murine PAB model: immune and vascular characterization** (`Supplementary_Figure_5.R`)
Mouse PAB snRNA-seq atlas (n=55,420 nuclei, 11 cell types); marker-based annotation; per-subject stacked-bar cell-type proportions across disease states; mouse EC subclustering with adult RV ecM1/ecM4/ecM7 module-score overlay (capturing arterial vs venous bias); mouse-EC stacked-bar proportions; human myeloid WGCNA module scores in mouse PAB myeloid clusters; per-mouse pseudobulk scores for the four human RVF myeloid programs (Kruskal–Wallis + pairwise Wilcoxon); flow-cytometry validation of immune populations in PAB heart; orthogonal bulk RNA-seq validation in flow-sorted cardiac macrophages.

**Figure S6 — Mouse PAB cardiomyocyte comparison with human RV failure** (`Supplementary_Figure_6.R`)
Reference mapping of mouse PAB CMs onto human RV CM subtypes; per-subcluster expression-score validation and mapping-score quantification; cross-species fold-change scatter (human RVF vs NF against mouse PAB severe vs NF) showing poor correlation; CM WGCNA module scores in mouse PAB by disease state (human trends poorly conserved); per-subject MitoCarta3.0 module score (mouse vs human); per-patient pseudobulk DESeq2 volcano of human MitoCarta3.0 genes with apeglm shrinkage; within-mouse Severe-vs-Moderate fold-change trajectories indicating mild PAB is transcriptionally distinct from severe; Sirius-red fibrosis quantification across PAB stages; cross-species concordance summary for three CM programs (fatty-acid oxidation, failure/fetal, glucocorticoid receptor).

**Figure S7 — Pediatric snRNA-seq myeloid and endothelial comparison with adult RV** (`Supplementary_Figure_7.R`)
Adult-RV myeloid marker scores in HLHS myeloid clusters; HLHS bulk WGCNA module-score dotplot by cell type × disease; GO biological-process enrichment for upregulated M8 DEGs (failing single-ventricle vs NF bi-ventricle); ChEA TF enrichments for NF and failing single-ventricle vs NF bi-ventricle, subset to M1 and M8 (CIITA, NR3C1, IRF8); per-patient HLHS myeloid program scores (one-way ANOVA per program); adult-RV EC hdWGCNA module scores in HLHS EC clusters; HLHS EC stacked-bar proportions (venous expansion); Phase-1 EC vasoprotective and IFN-engagement scores (pediatric vs adult); subtype-restricted Phase-2 readouts including *MECOM*, arterialization-TF score, Notch-target score, *SMAD1* and *NR2F2* in pediatric vs adult ECs.

**Figure S8 — Adult vs pediatric non-failing RV comparison** (`Supplementary_Figure_8.r`)
Co-embedding of healthy adult and pediatric NF RVs; cell-type abundance by dataset origin (one-way ANOVA); bulk WGCNA module expression by cell type and origin; pseudobulk CM volcano, PCA embedding, and GO enrichments for adult-NF vs pediatric-NF cardiomyocytes (baseline differences independent of disease).
