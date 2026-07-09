###############################################################################
## Figure 8 — Cross-species, cross-age, and cross-etiology comparisons
##
## Panels (AUTHORITATIVE — v57 manuscript legend,
##         .figure_run_logs/v57_figure_legends.md). v57 F8 is A–J:
##   (A) UMAP co-embedding of HLHS (n=209,604) + adult RV (n=61,398);
##       shared cell types.
##   (B) Cell-type abundance by disease state, HLHS dataset (one-way ANOVA).
##   (C) Bulk WGCNA module expression by cell type (top) and disease
##       state (bottom) for HLHS and adult RV.
##   (D) Scatter of pseudobulk log2FC: adult-RVF-vs-adult-NF vs
##       HLHS ped-RVF-vs-ped-NF. Genes absent in either / low-abundance
##       (<5% nuclei both) removed.
##   (E) Cardiomyocyte-specific module expression by disease state, HLHS
##       (M10/M26 rise ped-NF→HLHS→ped-RVF; M25/M28 no clear trend).
##   (F) GO-BP enrichment of upregulated DEGs from pooled M10/M25/M26/M28
##       for ped-RVF-vs-ped-HLHS (top) and ped-HLHS-vs-ped-NF (bottom).
##   (G) Seurat WGCNA module scores by cell type (top) and disease state
##       (bottom) — RV vs LV conservation.
##   (H) Per-nuclei violins of Seurat WGCNA module scores: M2 (left) and
##       pooled M10/M25/M26/M28 (right), NF vs DCM LV cardiomyocytes.
##   (I) Dot plot of Seurat WGCNA module scores over LV cardiomyocytes
##       split by disease state (M2 up, M10/M25/M26/M28 down, as in RV).
##   (J) Scatter of log2FC LV(DCM vs NF) vs RV(RVF vs NF) by gene; well
##       correlated, notable outliers labelled.
##
## v57 CROSS-CHECK / STATUS: this script is the legacy v52 SKELETON and
## does NOT yet implement the v57 A–J layout (it was authored to a v52
## A–M legend). Aligning F8 to the v57 panels above is the DEFERRED F8
## integration rewrite (tasks #42/#44) — to be done LAST per the user.
## The dead-code purge belongs to that rewrite (many commented
## `#<- readRDS(...)` lines are intentional cache-toggle docs, not dead).
##
## Source: copied from v51 Figure_8.R on 2026-04-10
##
## Output: ./output/Figure_8/Figure_8_panel_<L>.pdf (per-panel, v57)
###############################################################################

source('./helper_scripts/_shared_helpers.R')

## Per-figure output directory (introduced for consistent output paths)
V52_FIG_DIR <- './output/Figure_8'
dir.create(V52_FIG_DIR, showWarnings = FALSE, recursive = TRUE)


## Suppress R's default Rplots.pdf in cwd when Rscript hits a plot call
## that's outside an explicit pdf() ... dev.off() envelope.
pdf(NULL)
## Composite figure dimensions (inches) — used by theme_v52() scaling
COMP_W <- 12
COMP_H <- 5

## Publication scales (geom widths, point sizes, text sizes) scaled to COMP_W.
PS <- pub_scales(COMP_W)

## ── No dedicated cartoon in Figure_8 (Panel A is a full UMAP of all cell types).
## If a species/chamber cross-comparison schematic is later commissioned, save as
## ./new_scripts/assets/Figure_8_panel_A_schematic.png and compose via:
##   insert_asset('Figure_8_panel_A_schematic.png')
## Example re-crop from the published Figure_8.pdf (top-left UMAP label), if ever needed:
##   library(magick)
##   img  <- image_read_pdf('~/Downloads/hdWGCNA_TOM/Manuscripts/Figure_8.pdf',
##                          pages = 1, density = 300)
##   cart <- image_crop(img, '400x400+40+20')   # 17 cm composite ≈ 2008 px wide
##   image_write(cart, './new_scripts/assets/Figure_8_panel_A_schematic.png',
##               format = 'png', density = 300)

library(Seurat)
library(hdWGCNA)
library(ggeasy)
library(dplyr)


source('./helper_scripts/spatial_functions.R')



#######################################
#############  FIGURE 8A  #############
#######################################
M1 <- readRDS('./dependencies/shared/all_peds_data.rds')
M1$Names <- M1$cell.type
M1$NewNames <- M1$Names
M1 <- SetIdent(M1, value = "NewNames")

cluster.ids.new <- c("EC", "FB", "CM", "Myeloid","PC","Endo","NKT","SM","LEC","Adipo","Neuron","B","Mast")
names(cluster.ids.new) <- levels(M1)
M1 <- RenameIdents(M1, cluster.ids.new)
reorder_levels <- c("Adipo","CM","EC","Endo","FB","LEC","Myeloid","Neuron","NKT","PC","SM","Mast","B")
levels(M1) <- reorder_levels
M1$NewNames <- M1@active.ident


#Do some extra doublet removal


feats <- c("PLIN1","RYR2","VWF","LEPR","DCN","CCL21","CSF1R","NRXN1","CD2","PDGFRB","MYH11","KIT","MS4A1")
p<-VlnPlot(M1, features = feats,ncol=1,pt.size=F,group.by="NewNames")



for(i in 1:13) {  
   p[[i]] <- p[[i]] + NoLegend() + easy_remove_axes(which="y",what = c("ticks", "text","line")) + ggtitle("") + ylab(feats[i])
   if(i<13){p[[i]]<-p[[i]]+easy_remove_axes(which="x")}
   
}


pdf(paste0('./output/Figure_8/', 'snPeds_Vln.pdf'), width=4, height=20)
print(p)
dev.off()



M2 <- readRDS('./dependencies/shared/RV_data.rds')

## merge.SCTAssay's rbind(new_residual, old_residuals) fails when the two
## SCT objects' residual/scale.data feature sets differ (our
## all_peds_data.rds vs RV_data.rds differ from the legacy inputs that this
## merge originally worked on). Restrict BOTH objects to their SHARED SCT
## features first, so the original merge -> VariableFeatures -> RunPCA
## workflow runs unchanged (and the common-feature scale.data is smaller).
.f1 <- tryCatch(rownames(SeuratObject::GetAssayData(M1, assay = 'SCT', layer = 'scale.data')),
                error = function(e) character(0))
.f2 <- tryCatch(rownames(SeuratObject::GetAssayData(M2, assay = 'SCT', layer = 'scale.data')),
                error = function(e) character(0))
if (length(.f1) == 0) .f1 <- rownames(M1[['SCT']])
if (length(.f2) == 0) .f2 <- rownames(M2[['SCT']])
.common_sct <- intersect(.f1, .f2)
message('Pre-merge: restricting M1/M2 SCT to ', length(.common_sct),
        ' shared features (M1 ', length(.f1), ', M2 ', length(.f2), ').')
M1 <- subset(M1, features = .common_sct)
M2 <- subset(M2, features = .common_sct)

M2 <- merge(M1,M2)
VariableFeatures(M2[["SCT"]]) <- rownames(M2[["SCT"]]@scale.data)
M2$origin <- M2$orig.ident == "SeuratProject"
M2$patient[is.na(M2$patient)] = M2$sample[is.na(M2$patient)]

## Free M1 (peds, ~6 GB) — already merged into M2; downstream uses M2 only
## until line ~135 where M1 is restored from the subset. Saves enough RAM
## to keep RunHarmony + RunUMAP + DotPlot calls under the 48 GB system cap.
rm(M1); invisible(gc(verbose = FALSE))

## OOM-bypass: the committed pre-integrated co-embed (already carries
## pca/harmony/umap + CombinedNames/origin) lets us skip the OOM-prone
## RunPCA/RunHarmony/RunUMAP recompute on the full ~270k-cell merge
## (the documented intent of the old commented `#M2 <- readRDS(...merge)`).
.f8_merge_shrunk    <- './dependencies/Figure_8/RV_Peds_merge_shrunk.rds'
.cache_M2_integrated <- './output/Figure_8/fig8_M2_integrated_cache.rds'
if (file.exists(.f8_merge_shrunk)) {
  message('Loading committed pre-integrated co-embed ',
          '(RV_Peds_merge_shrunk.rds) — skips OOM-prone integration recompute.')
  rm(M2); invisible(gc(verbose = FALSE))
  M2 <- readRDS(.f8_merge_shrunk)
} else if (file.exists(.cache_M2_integrated)) {
  message('Loading cached M2 integrated (PCA/Harmony/UMAP/clusters)...')
  M2 <- readRDS(.cache_M2_integrated)
} else {
  M2 <- RunPCA(M2)

  M2 <- RunHarmony(M2,'patient')

  M2 <- M2%>%
    RunUMAP(reduction = "harmony", dims = 1:50, verbose = F) %>%
    FindNeighbors(reduction = "harmony", dims = 1:50) %>%
    FindClusters(resolution=0.5) %>%
    identity()
  ## The merged SCTModel (from older-class input objects) lacks the
  ## `median_umi` slot the current SeuratObject saveRDS hook
  ## (containsOutOfMemoryData) expects, so caching the integrated object
  ## errors. The cache is only a speed optimisation — don't let it halt F8.
  tryCatch(saveRDS(M2, .cache_M2_integrated),
           error = function(e)
             message('Skipping integrated-M2 cache (', conditionMessage(e), ')'))
}

M2$CombinedNames <- M2$NewNames
M2$CombinedNames[is.na(M2$NewNames)] <- M2$Names[is.na(M2$NewNames)]

## (SCT assay preserved as-is — SCT features were intersected pre-merge so
## the legacy merge / VariableFeatures / RunPCA workflow above works
## unchanged; no scale.data demote needed.)

#saveRDS(M2,'./output/Figure_8/RV_Peds_merge.rds')
#M2 <- readRDS('./output/Figure_8/RV_Peds_merge.rds')
#saveRDS(M1,'./output/Figure_8/Peds_clean.rds')
#M2 <- readRDS('./dependencies/Figure_8/RV_Peds_merge.rds')

## Restore M1 (peds-only subset). L109 `rm(M1)` freed RAM during integration;
## the L107-108 comment promised M1 would be "restored from the subset" here,
## but the restore line was left commented — so L149's PlotEmbedding(M1,...)
## crashed with `object 'M1' not found` (independent of the memory issue).
M1 <- subset(M2, origin == TRUE)


pdf(paste0('./output/Figure_8/', 'sn_RV_Peds_UMAP.pdf'), width=5, height=5)
p_8A <- PlotEmbedding(M2,group.by='CombinedNames',point_size=PS$umap_pt,plot_under=TRUE,plot_theme=umap_theme()+NoLegend(),raster_dpi=400,raster_scale=0.5)
print(p_8A)
dev.off()
save_figure(p_8A, 'Figure_8_panel_A.pdf', width = 5, height = 5)

pdf(paste0('./output/Figure_8/', 'sn_RV_Peds_UMAP_origin.pdf'), width=5, height=5)
print(PlotEmbedding(M2,group.by='origin',point_size=PS$umap_pt,plot_under=TRUE,plot_theme=umap_theme()+NoLegend(),raster_dpi=400,raster_scale=0.5))
dev.off()


pdf(paste0('./output/Figure_8/', 'sn_Peds_UMAP_reprojected.pdf'), width=5, height=5)
print(PlotEmbedding(M1,group.by='NewNames',point_size=PS$umap_pt,plot_under=TRUE,plot_theme=umap_theme()+NoLegend(),raster_dpi=400,raster_scale=0.5))
dev.off()
rm(M1); invisible(gc(verbose = FALSE))   # peds subset no longer needed until L511

M2$condition[M2$condition == "NF"] = "pRV" 
M2$condition[M2$condition == "Donor"] = "NF" 
M2$condition[M2$condition == "SystolicHF"] = "RVF" 



M2$group[is.na(M2$group)] <- M2$condition[is.na(M2$group)]


#######################################
#############  FIGURE 8B  #############
#######################################

t(table(M2$CombinedNames,M2$group))/rowSums(t(table(M2$CombinedNames,M2$group)))


NF_percent_cell <- cbind(as.data.frame(table(subset(M2,group=="NF")@active.ident)/length(subset(M2,group=="NF")@active.ident)*100),type = "NF")
NF_percent_cell$sum <- (rev(cumsum(rev(NF_percent_cell$Freq))) - NF_percent_cell$Freq/2)/100
NF_percent_cell$Freq <- NF_percent_cell$Freq/100


pRV_percent_cell <- cbind(as.data.frame(table(subset(M2,group=="pRV")@active.ident)/length(subset(M2,group=="pRV")@active.ident)*100),type = "pRV")
pRV_percent_cell$sum <- (rev(cumsum(rev(pRV_percent_cell$Freq))) - pRV_percent_cell$Freq/2)/100
pRV_percent_cell$Freq <- pRV_percent_cell$Freq/100


RVF_percent_cell <- cbind(as.data.frame(table(subset(M2,group=="RVF")@active.ident)/length(subset(M2,group=="RVF")@active.ident)*100),type = "RVF")
RVF_percent_cell$sum <- (rev(cumsum(rev(RVF_percent_cell$Freq))) - RVF_percent_cell$Freq/2)/100
RVF_percent_cell$Freq <- RVF_percent_cell$Freq/100

percent_cell_df <- rbind(NF_percent_cell,pRV_percent_cell,RVF_percent_cell)
pdf('./output/Figure_8/RV_Peds_prev_stacked.pdf',width=6,height=2.5)
print(
ggplot(percent_cell_df, aes(fill=Var1, y=Freq, x=type,label=round(sum,1))) +
geom_bar(position="stack", stat="identity",width=0.6, linewidth = PS$geom_lw) + theme_v52(COMP_W) + coord_flip()+
scale_fill_lineage(name = "Cell type") +
xlab("Disease State") + ylab("Frequency") + labs(fill="Cell type",color='black') + scale_y_continuous(expand=c(0,0)) +
geom_label_repel(aes(type,sum,label=scales::percent(round(Freq,2))),fill=NA,nudge_x=0.5,direction="y",
                 size = PS$text_mm, family = FONT_FAMILY)
)
dev.off()

cells <- table(M2$CombinedNames,M2$patient)
cells <- cells[c('CM','FB','EC','Myeloid'),]
cells <- sweep(cells,2,colSums(cells),'/')
cells <- data.frame(cells)

cells$group<-M2$group[match(cells$Var2,M2$patient)]

cells$group<-as.factor(cells$group)

library(ggpubr)
pdf('./output/Figure_8/Peds_RV_clust_freq.pdf',width=8,height=5)
p_8B <- ggboxplot(cells[length(cells$group):1,],x="group",y="Freq",fill="group",group="group", size = PS$geom_lw)+
  theme_v52(COMP_W) +
  scale_fill_disease() +
  labs(color='Group',x="Disease",y='Frequency') +
  facet_wrap(~Var1,ncol=7) +
  scale_y_touch() +
  stat_compare_means(aes(group=group),method="anova", size = PS$text_mm, symnum.args = pub_signif_args)
print(p_8B)
dev.off()
save_figure(p_8B, 'Figure_8_panel_B.pdf', width = 8, height = 5)

cells <- table(M2$NewNames,M2$patient)
cells <- cells[,12:25]
cells <- cells[c('CM','FB','EC','Myeloid'),]
cells <- sweep(cells,2,colSums(cells),'/')
cells <- data.frame(cells)

cells$group<-M2$group[match(cells$Var2,M2$patient)]

cells$group<-as.factor(cells$group)

library(ggpubr)
pdf('./output/Figure_8/Peds_clust_freq.pdf',width=8,height=5)
p <- ggboxplot(cells[length(cells$group):1,],x="group",y="Freq",fill="group",group="group", size = PS$geom_lw)+
  theme_v52(COMP_W) +
  scale_fill_disease() +
  labs(color='Group',x="Disease",y='Frequency') +
  facet_wrap(~Var1,ncol=7) +
  scale_y_touch() +
  stat_compare_means(aes(group=group),method="anova", size = PS$text_mm, symnum.args = pub_signif_args)
print(p)
dev.off()

#######################################
#############  FIGURE 8C  #############
#######################################

########Embed bulk in single nuc peds


######Load module
consensus_modules <- read.csv("./dependencies/shared/bulk_heart_modules.csv")
consensus_modules <- consensus_modules[,1:3]
consensus_modules <- subset(consensus_modules, gene_name %in% rownames(M2))


# remove duplicate gene names
consensus_modules <- consensus_modules[match(unique(consensus_modules$gene_name), consensus_modules$gene_name),]



## Skip the heavy SetupForWGCNA/ScaleData/ProjectModules machinery when the
## committed bulk2sn MEs dependency is present (ProjectModules fails on
## harmony 2.x / Seurat-v5, and ScaleData on the ~210k-cell co-embed is an
## OOM risk). MEs are loaded directly from the dependency below.
.f8_MEs_dep <- './dependencies/shared/scWGCNA_bulk2sn_MEs.rds'
if (!file.exists(.f8_MEs_dep)) {
#A bit hacky since ref here should be NULL and wgcna_name should be None but function errors then. This will be overwritten by modules=consensus_modules anyway
DefaultAssay(M2) <- 'RNA'
M2<-FindVariableFeatures(M2)
M2 <- ScaleData(M2,block.size=1000)

M2 <- SetupForWGCNA(
  M2,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "Cardiomyocyte" # the name of the hdWGCNA experiment
)


.cache_M2_projected <- './output/Figure_8/fig8_M2_projected_cache.rds'
if (file.exists(.cache_M2_projected)) {
  message('Loading cached M2 ProjectModules result...')
  M2 <- readRDS(.cache_M2_projected)
} else {
  ## hdWGCNA's ProjectModules with group.by.vars internally calls
  ## RunHarmony(assay.use=...), which harmony 2.x rejects via
  ## check_legacy_args ("Argument assay.use is unhandled"). Drop
  ## group.by.vars — the merged M2 is already harmonised (same fix as
  ## F7 / S8). Guard the cache saveRDS for the legacy-SCTModel
  ## median_umi serialization issue.
  ## ProjectModules can fail inside hdWGCNA/.CreateStdAssay on this
  ## Seurat-v5 stack. Don't halt — fall back to the committed scWGCNA
  ## bulk->sn MEs dependency at the GetMEs step below (same approach
  ## that fixed F6 Panel D).
  M2 <- tryCatch(
    ProjectModules(M2, modules = consensus_modules, seurat_ref = M2,
                   wgcna_name = "Cardiomyocyte", wgcna_name_proj = 'bulk2sn'),
    error = function(e) {
      message('ProjectModules failed (', conditionMessage(e),
              ') — will use cached scWGCNA_bulk2sn MEs dependency.'); M2 })
  tryCatch(saveRDS(M2, .cache_M2_projected),
           error = function(e)
             message('Skipping M2-projected cache (', conditionMessage(e), ')'))
}
}   # end if(!file.exists(.f8_MEs_dep)) — heavy ProjectModules path

mapping <- labels2colors(1:100)
if (file.exists(.f8_MEs_dep)) {
  message('Loading committed scWGCNA bulk2sn MEs dependency (', .f8_MEs_dep, ')')
  MEs <- readRDS(.f8_MEs_dep)
} else {
  M2 <- SetActiveWGCNA(M2, 'bulk2sn')
  .cache_MEs_bulk2sn <- './output/Figure_8/fig8_MEs_bulk2sn_cache.rds'
  if (file.exists(.cache_MEs_bulk2sn)) {
    message('Loading cached GetMEs (bulk2sn harmonized)...')
    MEs <- readRDS(.cache_MEs_bulk2sn)
  } else {
    MEs <- GetMEs(M2, harmonized=TRUE)
    saveRDS(MEs, .cache_MEs_bulk2sn)
  }
}
mods <- colnames(MEs); mods <- mods[mods != 'grey']
mods_num <- paste0('M',match(mods,mapping))
prv_vs_rv_signif <- c('M3','M4','M5','M10','M11','M12','M14')
all_signif <- c('M1','M2','M3','M4','M5','M8','M10','M11','M12','M14','M20','M25','M26','M28')


colnames(MEs)<-paste0('M',match(colnames(MEs),mapping))
## Align MEs to M2 cells. The committed bulk2sn MEs dependency covers the
## adult-RV subset (61,398 cells); M2 is the full co-embed (209,604).
## cbind is POSITIONAL, so build an MEs matrix indexed to M2's barcodes
## (NA where a cell has no projected ME) — per-celltype/disease
## aggregation then simply ignores the NA (peds) cells.
if (nrow(MEs) != nrow(M2@meta.data)) {
  .me_full <- matrix(NA_real_, nrow = nrow(M2@meta.data), ncol = ncol(MEs),
                     dimnames = list(rownames(M2@meta.data), colnames(MEs)))
  .shared <- intersect(rownames(MEs), rownames(M2@meta.data))
  message(sprintf('Aligning MEs to M2: %d / %d M2 cells have projected MEs',
                  length(.shared), nrow(M2@meta.data)))
  if (length(.shared) > 0)
    .me_full[.shared, ] <- as.matrix(MEs)[.shared, , drop = FALSE]
  MEs <- as.data.frame(.me_full, check.names = FALSE)
}
M2@meta.data <- cbind(M2@meta.data, MEs)
M2 <- SetIdent(M2, value = "CombinedNames")


#consensus_modules <- read.csv("./dependencies/shared/bulk_heart_modules.csv")
#consensus_modules <- consensus_modules[,1:3]
#consensus_modules <- subset(consensus_modules, gene_name %in% rownames(M2))
# remove duplicate gene names
#consensus_modules <- consensus_modules[match(unique(consensus_modules$gene_name), consensus_modules$gene_name),]


library(dplyr)
score_calc <- consensus_modules %>% group_by(module) %>% group_split()
module_colors <- unique(unlist(lapply(score_calc,'[[','module')))
module_colors <- paste0('M',match(module_colors,mapping))
#saveRDS(M2, './output/Figure_8/scWGCNA_RV_Peds_bulk2sn_projection.rds')
#M2<- readRDS('./output/Figure_8/scWGCNA_RV_Peds_bulk2sn_projection.rds')

DefaultAssay(M2) <- 'SCT'



#rm(seurat_ref)
#gc()
#seurat_ref<-readRDS('./dependencies/shared/RV_data.rds')
#seurat_ref <- SetIdent(seurat_ref, value = "Names")
#seurat_ref@meta.data <- cbind(seurat_ref@meta.data, MEs)



.cache_M2_modulescore <- './output/Figure_8/fig8_M2_modulescore_cache.rds'
if (file.exists(.cache_M2_modulescore)) {
  message('Loading cached AddModuleScore columns (M2)...')
  .ms_cols <- readRDS(.cache_M2_modulescore)
  M2@meta.data[, colnames(.ms_cols)] <- .ms_cols
} else {
  M2 <- AddModuleScore(M2,lapply(score_calc,'[[','gene_name'),name="module_score")
  cols_current <- colnames(M2@meta.data)
  cols_current[startsWith(colnames(M2@meta.data),'module_score')] <- paste0('module_',module_colors)
  colnames(M2@meta.data) <- cols_current
  .ms_col_names <- paste0('module_', module_colors)
  .ms_cols <- M2@meta.data[, .ms_col_names, drop = FALSE]
  saveRDS(.ms_cols, .cache_M2_modulescore)
}

M2$origin[M2$origin] = 'Peds'
M2$origin[M2$origin == FALSE] = 'RV'

M2$CombinedNamesSplit <- paste0(M2$CombinedNames,'_',M2$origin)

pdf(paste0('./output/Figure_8/', 'RV_Peds_Dot.pdf'), width=7, height=5)
p <- DotPlot(M2,paste0('module_',all_signif),group.by='CombinedNamesSplit',dot.min=0,col.min=0,col.max=2,idents=c("CM","EC","FB","Myeloid","PC","SM")) +
  RotatedAxis() + ylab('')+ xlab('')+
  scale_color_gradient2(high='red', mid='grey95', low='blue') +
  theme(
    panel.border = element_rect(size=1,fill=NA, color='black'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()
) 
print(p)
dev.off()


pdf(paste0('./output/Figure_8/', 'RV_Peds_Dot_ordered.pdf'), width=7, height=5)

p_8C_top <- DotPlot(M2,paste0('module_',c('M20','M5','M1','M3','M4','M8','M2','M12','M25','M26','M10','M28','M14','M11')),group.by='CombinedNamesSplit',dot.min=0,col.min=0,col.max=2,idents=c("CM","EC","FB","Myeloid","PC","SM")) +
  RotatedAxis() + ylab('')+ xlab('')+
  scale_color_gradient2(high='red', mid='grey95', low='blue') +
  theme(
    panel.border = element_rect(size=1,fill=NA, color='black'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()
) 
print(p_8C_top)
dev.off()

M2$groupSplit <- paste0(M2$group,'_',M2$origin)


pdf(paste0('./output/Figure_8/', 'RV_Peds_Dot_disease_ordered.pdf'), width=7, height=4)

p_8C_bot <- DotPlot(M2,paste0('module_',c('M20','M5','M1','M3','M4','M8','M2','M12','M25','M26','M10','M28','M14','M11')),group.by='groupSplit',dot.min=0,col.min=0,col.max=2,idents=c("CM","EC","FB","Myeloid","PC","SM")) +
  RotatedAxis() + ylab('')+ xlab('')+
  scale_color_gradient2(high='red', mid='grey95', low='blue') +
  theme(
    panel.border = element_rect(size=1,fill=NA, color='black'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()
) 
print(p_8C_bot)
dev.off()
p_8C <- p_8C_top / p_8C_bot
save_figure(p_8C, 'Figure_8_panel_C.pdf', width = 7, height = 9)



M2$group_split <- paste0(M2$group,'_',M2$origin)
M2$group_split <- factor(M2$group_split,levels=c('NF_RV','pRV_RV','RVF_RV','NF_Peds','pRV_Peds','RVF_Peds'))

pdf(paste0('./output/Figure_8/', 'RV_Peds_Dot_CM_Peds.pdf'), width=7, height=3)
p <- DotPlot(subset(M2,CombinedNames=='CM' & origin=='Peds'),paste0('module_',c('M2','M12','M28','M10','M25','M26')),group.by='group_split',dot.min=0,col.min=-2,col.max=2,idents=c("CM","EC","FB","Myeloid","PC","SM")) +
  RotatedAxis() + ylab('')+ xlab('')+
  scale_color_gradient2(high='red', mid='grey95', low='blue') +
  theme(
    panel.border = element_rect(size=1,fill=NA, color='black'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()
) 
print(p)
dev.off()

pdf(paste0('./output/Figure_8/', 'RV_Peds_Dot_CM_RV.pdf'), width=7, height=3)
p <- DotPlot(subset(M2,CombinedNames=='CM' & origin=='RV'),paste0('module_',c('M2','M12','M28','M10','M25','M26')),group.by='group_split',dot.min=0,col.min=-2,col.max=2,idents=c("CM","EC","FB","Myeloid","PC","SM")) +
  RotatedAxis() + ylab('')+ xlab('')+
  scale_color_gradient2(high='red', mid='grey95', low='blue') +
  theme(
    panel.border = element_rect(size=1,fill=NA, color='black'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()
) 
print(p)
dev.off()

pdf(paste0('./output/Figure_8/', 'RV_Peds_Dot_Myeloid_Peds.pdf'), width=7, height=3)
p <- DotPlot(subset(M2,CombinedNames=='Myeloid' & origin=='Peds'),paste0('module_',c('M1','M3','M4','M8')),group.by='group_split',dot.min=0,col.min=-2,col.max=2,idents=c("CM","EC","FB","Myeloid","PC","SM")) +
  RotatedAxis() + ylab('')+ xlab('')+
  scale_color_gradient2(high='red', mid='grey95', low='blue') +
  theme(
    panel.border = element_rect(size=1,fill=NA, color='black'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()
) 
print(p)
dev.off()

pdf(paste0('./output/Figure_8/', 'RV_Peds_Dot_Myeloid_RV.pdf'), width=7, height=3)
p <- DotPlot(subset(M2,CombinedNames=='Myeloid' & origin=='RV'),paste0('module_',c('M1','M3','M4','M8')),group.by='group_split',dot.min=0,col.min=-2,col.max=2,idents=c("CM","EC","FB","Myeloid","PC","SM")) +
  RotatedAxis() + ylab('')+ xlab('')+
  scale_color_gradient2(high='red', mid='grey95', low='blue') +
  theme(
    panel.border = element_rect(size=1,fill=NA, color='black'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()
) 
print(p)
dev.off()


#######################################
#############  FIGURE 8D  #############
#######################################

M2 <- SetIdent(M2, value = "groupSplit")
## SCT was demoted to a plain assay (no SCTModel.list); this median_umi prep
## is only needed for SCT residual recorrection, which we skip
## (FindMarkers below uses recorrect_umi=FALSE on the data layer).
if ('SCTModel.list' %in% methods::slotNames(M2[['SCT']]) &&
    length(M2[['SCT']]@SCTModel.list) > 0)
  tryCatch(
    for (.i in seq_along(M2[['SCT']]@SCTModel.list))
      slot(M2$SCT@SCTModel.list[[.i]], 'median_umi') <-
        median(M2$SCT@SCTModel.list[[.i]]@cell.attributes$umi),
    error = function(e)
      message('median_umi prep skipped (', conditionMessage(e),
              ') — FindMarkers uses recorrect_umi=FALSE so this is fine'))

.cache_gene_set_RV <- './output/Figure_8/fig8_gene_set_RV_cache.rds'
if (file.exists(.cache_gene_set_RV)) {
  message('Loading cached FindMarkers gene_set_RV...')
  gene_set_RV <- readRDS(.cache_gene_set_RV)
} else {
  gene_set_RV <- FindMarkers(M2, ident.1 = "RVF_RV", ident.2 = "NF_RV",recorrect_umi=F)
  saveRDS(gene_set_RV, .cache_gene_set_RV)
}

.cache_gene_set_Peds <- './output/Figure_8/fig8_gene_set_Peds_cache.rds'
if (file.exists(.cache_gene_set_Peds)) {
  message('Loading cached FindMarkers gene_set_Peds...')
  gene_set_Peds <- readRDS(.cache_gene_set_Peds)
} else {
  gene_set_Peds <- FindMarkers(M2, ident.1 = "RVF_Peds", ident.2 = "NF_Peds",recorrect_umi=F)
  saveRDS(gene_set_Peds, .cache_gene_set_Peds)
}




shared <- intersect(rownames(gene_set_RV),rownames(gene_set_Peds))
dataset <- data.frame(Peds=gene_set_Peds[shared,]$avg_log2FC,RV=gene_set_RV[shared,]$avg_log2FC)
rownames(dataset) <- shared
#Cor 0.5288722

pdf(paste0('./output/Figure_8/', 'Peds_vs_RV.pdf'), width=8, height=8)
p_8D <- ggplot(dataset, aes(x = RV, y=Peds)) + geom_point(size = PS$scatter_pt) +
  geom_text_repel(label=rownames(dataset),max.overlaps = 50, size = PS$text_mm, family = FONT_FAMILY, fontface = "italic") + theme_v52(COMP_W)
print(p_8D)
dev.off()
save_figure(p_8D, 'Figure_8_panel_D.pdf', width = 8, height = 8)

shared <- intersect(rownames(subset(gene_set_RV,pct.1>0.05 & pct.2>0.05)),
  rownames(subset(gene_set_Peds,pct.1>0.05 & pct.2>0.05)))

dataset <- data.frame(Peds=gene_set_Peds[shared,]$avg_log2FC,RV=gene_set_RV[shared,]$avg_log2FC)
rownames(dataset) <- shared


pdf(paste0('./output/Figure_8/', 'Peds_vs_RV_5percent.pdf'), width=8, height=8)
print(
ggplot(dataset, aes(x = RV, y=Peds)) + geom_point(size = PS$scatter_pt) +
  geom_text_repel(label=rownames(dataset),max.overlaps = 20, size = PS$text_mm, family = FONT_FAMILY, fontface = "italic") + theme_v52(COMP_W)
)
dev.off()
#Cor 0.4706919

#######################################
#############  FIGURE 8E  #############
#######################################

M1 <- readRDS('./dependencies/Figure_8/cardiomyocyte_annotated.rds')
M1$Names <- M1$cell.type
M1$NewNames <- M1$Names
M1$Subnames <- M1$sub.type
M1$NewSubnames <- M1$Subnames
M1 <- SetIdent(M1, value = "NewSubnames")
DefaultAssay(M1) <- 'SCT'


consensus_modules <- read.csv("./dependencies/shared/bulk_heart_modules.csv")
consensus_modules <- consensus_modules[,1:3]
consensus_modules <- subset(consensus_modules, gene_name %in% rownames(M1))
# remove duplicate gene names
consensus_modules <- consensus_modules[match(unique(consensus_modules$gene_name), consensus_modules$gene_name),]
library(dplyr)
score_calc <- consensus_modules %>% group_by(module) %>% group_split()
module_colors <- unique(unlist(lapply(score_calc,'[[','module')))
mapping <- labels2colors(1:100)
module_colors <- paste0('M',match(module_colors,mapping))
.cache_M1_modulescore <- './output/Figure_8/fig8_M1_modulescore_cache.rds'
if (file.exists(.cache_M1_modulescore)) {
  message('Loading cached AddModuleScore columns (M1)...')
  .ms1_cols <- readRDS(.cache_M1_modulescore)
  M1@meta.data[, colnames(.ms1_cols)] <- .ms1_cols
} else {
  M1 <- AddModuleScore(M1,lapply(score_calc,'[[','gene_name'),name="module_score")
  cols_current <- colnames(M1@meta.data)
  cols_current[startsWith(colnames(M1@meta.data),'module_score')] <- paste0('module_',module_colors)
  colnames(M1@meta.data) <- cols_current
  .ms1_col_names <- paste0('module_', module_colors)
  .ms1_cols <- M1@meta.data[, .ms1_col_names, drop = FALSE]
  saveRDS(.ms1_cols, .cache_M1_modulescore)
}

Idents(M1) <- factor(x = Idents(M1), levels = sort(levels(M1)))


#Dot Plot of enrichment cell type of CM enriched modules

pdf(paste0('./output/Figure_8/', 'peds_mod_trend_subcluster_CM.pdf'), width=4.5, height=3)

p_8E_top <- DotPlot(M1,paste0('module_',
  c('M2','M12','M25','M26','M10','M28')),dot.min=0,col.min=0,col.max=2) +
  RotatedAxis() + ylab('')+ xlab('')+
  scale_color_gradient2(high='red', mid='grey95', low='blue') +
  theme(
    panel.border = element_rect(size=1,fill=NA, color='black'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()
) 
print(p_8E_top)
dev.off()


M1 <- SetIdent(M1, value = "condition")
Idents(M1) <- factor(x = Idents(M1), levels = c('Donor','NF','SystolicHF'))


pdf(paste0('./output/Figure_8/', 'peds_mod_trend_condition_CM.pdf'), width=5, height=2.5)

p_8E_bot <- DotPlot(M1,paste0('module_',
  c('M2','M12','M25','M26','M10','M28')),dot.min=0,col.min=0,col.max=2) +
  RotatedAxis() + ylab('')+ xlab('')+
  scale_color_gradient2(high='red', mid='grey95', low='blue') +
  theme(
    panel.border = element_rect(size=1,fill=NA, color='black'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()
) 
print(p_8E_bot)
dev.off()
p_8E <- p_8E_top / p_8E_bot
save_figure(p_8E, 'Figure_8_panel_E.pdf', width = 5, height = 5.5)



#######################################
#############  FIGURE 8F  #############
#######################################

## ===== OLD F8F (per-cell FindMarkers + enrichR GO on mito modules) DISABLED (v67). =====
## Superseded by the pseudobulk MitoCarta volcanoes built after this block: the per-cell
## Wilcoxon pseudoreplicates ~30k CM cells and the enrichR step is a circular GO of a
## mito-module gene input, so it cannot establish mitochondrial induction.
if (FALSE) {

M1 <- SetIdent(M1, value = "condition")
Idents(M1) <- factor(x = Idents(M1), levels = c('Donor','NF','SystolicHF'))

dbs <-c('GO_Biological_Process_2023','GO_Cellular_Component_2023','GO_Molecular_Function_2023','Reactome_2022', 'ChEA_2022',"LINCS_L1000_Chem_Pert_Consensus_Sigs")

library(enrichR)
#Run enrichment by cell type
.cache_fig8f_combined_set    <- './output/Figure_8/fig8_fig8F_combined_set_cache.rds'
.cache_fig8f_combined_output <- './output/Figure_8/fig8_fig8F_combined_output_cache.rds'
if (file.exists(.cache_fig8f_combined_set) && file.exists(.cache_fig8f_combined_output)) {
  message('Loading cached Figure 8F FindMarkers + enrichR loop results...')
  combined_set    <- readRDS(.cache_fig8f_combined_set)
  combined_output <- readRDS(.cache_fig8f_combined_output)
} else {
  combined_set <- data.frame()
  combined_output <- data.frame()
  bulk_modules <- consensus_modules
  bulk_modules$module <- match(consensus_modules$module,mapping)
  mods_idx <- list(2,12,10,25,26,28)
  cell_types <- unique(M1$NewNames)
  comparison <- list(c("SystolicHF","Donor"),c("SystolicHF","NF"))
  for (i in mods_idx){
    for (j in cell_types){
      for (k in comparison){
        key_genes <- subset(bulk_modules,module %in% c(i))$gene_name
        key_genes <- key_genes[key_genes %in% rownames(M1)]

        gene_set <- FindMarkers(M1, ident.1 = paste0(k[1]), ident.2 = paste0(k[2]),features=key_genes)

        gene_set<-subset(gene_set,p_val_adj<0.05)
        if (length(rownames(gene_set))==0){next}
        gene_set$module <- paste0('M',i)
        gene_set$color <- mapping[i]
        gene_set$comparison <- paste0(k[1],'_',k[2])
        gene_set$celltype <- j

        if (length(combined_set) == 0){
          combined_set <- gene_set
        }
        else {
          combined_set <- rbind(combined_set,gene_set)
        }

        gene_enrich <- subset(gene_set,avg_log2FC<0)
        enriched <- enrichR::enrichr(rownames(gene_enrich), dbs)
        Sys.sleep(5)
        for(db in names(enriched)){
            cur_df <- enriched[[db]]
            if (nrow(cur_df) > 1){
              cur_df$db <- db
              cur_df$module <- paste0('M',i)
              cur_df$celltype <- j
              cur_df$comparison <- paste0(k[1],'_',k[2])
              cur_df$color <- mapping[i]
              cur_df$direction <- 'down'
              combined_output <- rbind(combined_output, cur_df)
            }
        }

        gene_enrich <- subset(gene_set,avg_log2FC>0)
        enriched <- enrichR::enrichr(rownames(gene_enrich), dbs)
        Sys.sleep(5)
        for(db in names(enriched)){
            cur_df <- enriched[[db]]
            if (nrow(cur_df) > 1){
              cur_df$db <- db
              cur_df$module <- paste0('M',i)
              cur_df$celltype <- j
              cur_df$comparison <- paste0(k[1],'_',k[2])
              cur_df$color <- mapping[i]
              cur_df$direction <- 'up'
              combined_output <- rbind(combined_output, cur_df)
            }
        }
      }
    }
  }
  saveRDS(combined_set,    .cache_fig8f_combined_set)
  saveRDS(combined_output, .cache_fig8f_combined_output)
}

wrapText <- function(x, len) {
    sapply(x, function(y) paste(strwrap(y, len), collapse = "\n"), USE.NAMES = FALSE)
}


###Systolic HF vs Donor

#Up
selected_terms <- subset(combined_output,db=="GO_Biological_Process_2023")
selected_terms <- subset(selected_terms,direction=="up")
selected_terms <- subset(selected_terms,comparison=="SystolicHF_Donor")
selected_terms <- subset(selected_terms,color %in% mapping[c(2,12)])


# subset selected terms
selected_terms <- subset(selected_terms, P.value < 0.05)
selected_terms$module_celltype <- paste0(selected_terms$module,'_',selected_terms$celltype)
idx_top_1 <- match(unique(selected_terms$module_celltype),selected_terms$module_celltype)
idx_top_5 <- sort(c(idx_top_1,idx_top_1+1,idx_top_1+2,idx_top_1+3,idx_top_1+4))

selected_terms<-selected_terms[idx_top_5,]


selected_terms$group <- factor(
  as.character(selected_terms$module_celltype),
  levels = unique(selected_terms$module_celltype)
)


# set max pval
quantile(-log(selected_terms$P.value), 0.95)
max_p <- 10

selected_terms$logp <- -log(selected_terms$P.value)
selected_terms$logp <- ifelse(selected_terms$logp > max_p, max_p, selected_terms$logp)

# remove Reactome Term ID
library(stringr)
selected_terms$Term <- str_replace(selected_terms$Term, "\\(GO.*", "")

selected_terms <- selected_terms %>%
  arrange(group)


selected_terms$wrap <- wrapText(selected_terms$Term, 45)

selected_terms$Term <- factor(
  as.character(selected_terms$Term),
  levels = rev(unique(as.character(selected_terms$Term)))
)

selected_terms$wrap <- factor(
  as.character(selected_terms$wrap),
  levels = rev(unique(as.character(selected_terms$wrap)))
)

library(viridis)

# Reactome Term dot plot
p <- selected_terms %>%
  ggplot(aes(x = group, y = wrap, color =logp, size=log(Combined.Score))) +
  geom_point() +
  scale_color_stepsn(colors=rev(magma(256))) +
  RotatedAxis() + xlab('') + ylab('') +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_rect(size=1, color='black', fill=NA),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(0,0,0,0),
    panel.grid = element_line(size=0.25, color='lightgrey')
  )


mapping <- labels2colors(1:100)


color_df <- data.frame(group = selected_terms$group,colour = selected_terms$color)



color_df$group <- factor(
  as.character(color_df$group),
  levels = unique(color_df$group)
)


c_vect <- color_df$colour
names(c_vect) <- color_df$group



# make the colorbar as its own heatmap
color_df$var <- 1
colorbar <- color_df %>%
  ggplot(aes(x=group, y=var, fill=group)) +
  geom_tile() +
  scale_fill_manual(values=c_vect) +
  coord_equal() +
  NoLegend() + RotatedAxis() +
  theme(
    plot.title=element_blank(),
    axis.line=element_blank(),
    axis.ticks.y =element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    plot.margin=margin(0,0,0,0),
  )


pdf(paste0('./output/Figure_8/', 'peds_CM_by_cluster_terms_cell_type_up_SystolicHF_vs_Donor.pdf'), width=6, height=4)
print(p / colorbar)
dev.off()



###RVF vs NF

#Down
selected_terms <- subset(combined_output,db=="GO_Biological_Process_2023")
selected_terms <- subset(selected_terms,direction=="down")
selected_terms <- subset(selected_terms,comparison=="SystolicHF_Donor")
selected_terms <- subset(selected_terms,color %in% mapping[c(10,25,26,28)])


# subset selected terms
selected_terms <- subset(selected_terms, P.value < 0.05)
selected_terms$module_celltype <- paste0(selected_terms$module,'_',selected_terms$celltype)
idx_top_1 <- match(unique(selected_terms$module_celltype),selected_terms$module_celltype)
idx_top_5 <- sort(c(idx_top_1,idx_top_1+1,idx_top_1+2,idx_top_1+3,idx_top_1+4))

selected_terms<-selected_terms[idx_top_5,]


selected_terms$group <- factor(
  as.character(selected_terms$module_celltype),
  levels = unique(selected_terms$module_celltype)
)


# set max pval
quantile(-log(selected_terms$P.value), 0.95)
max_p <- 10

selected_terms$logp <- -log(selected_terms$P.value)
selected_terms$logp <- ifelse(selected_terms$logp > max_p, max_p, selected_terms$logp)

# remove Reactome Term ID
library(stringr)
selected_terms$Term <- str_replace(selected_terms$Term, "\\(GO.*", "")

selected_terms <- selected_terms %>%
  arrange(group)


selected_terms$wrap <- wrapText(selected_terms$Term, 45)

selected_terms$Term <- factor(
  as.character(selected_terms$Term),
  levels = rev(unique(as.character(selected_terms$Term)))
)

selected_terms$wrap <- factor(
  as.character(selected_terms$wrap),
  levels = rev(unique(as.character(selected_terms$wrap)))
)

library(viridis)

# Reactome Term dot plot
p <- selected_terms %>%
  ggplot(aes(x = group, y = wrap, color =logp, size=log(Combined.Score))) +
  geom_point() +
  scale_color_stepsn(colors=rev(magma(256))) +
  RotatedAxis() + xlab('') + ylab('') +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_rect(size=1, color='black', fill=NA),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(0,0,0,0),
    panel.grid = element_line(size=0.25, color='lightgrey')
  )


mapping <- labels2colors(1:100)


color_df <- data.frame(group = selected_terms$group,colour = selected_terms$color)



color_df$group <- factor(
  as.character(color_df$group),
  levels = unique(color_df$group)
)


c_vect <- color_df$colour
names(c_vect) <- color_df$group



# make the colorbar as its own heatmap
color_df$var <- 1
colorbar <- color_df %>%
  ggplot(aes(x=group, y=var, fill=group)) +
  geom_tile() +
  scale_fill_manual(values=c_vect) +
  coord_equal() +
  NoLegend() + RotatedAxis() +
  theme(
    plot.title=element_blank(),
    axis.line=element_blank(),
    axis.ticks.y =element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    plot.margin=margin(0,0,0,0),
  )


pdf(paste0('./output/Figure_8/', 'peds_CM_by_cluster_terms_cell_type_down_SystolicHF_vs_Donor.pdf'), width=6, height=4)
print(p / colorbar)
dev.off()



#Both
selected_terms <- subset(combined_output,db=="GO_Biological_Process_2023")
#selected_terms <- subset(selected_terms,direction=="up")
selected_terms <- subset(selected_terms,comparison=="SystolicHF_Donor")
selected_terms <- subset(selected_terms,color %in% mapping[c(2,12,10,25,26,28)])


# subset selected terms
selected_terms <- subset(selected_terms, P.value < 0.05)
selected_terms$module_celltype <- paste0(selected_terms$module,'_',selected_terms$celltype)
idx_top_1 <- match(unique(selected_terms$module_celltype),selected_terms$module_celltype)
idx_top_5 <- sort(c(idx_top_1,idx_top_1+1,idx_top_1+2,idx_top_1+3,idx_top_1+4))

selected_terms<-selected_terms[idx_top_5,]


selected_terms$group <- factor(
  as.character(selected_terms$module_celltype),
  levels = unique(selected_terms$module_celltype)
)


# set max pval
quantile(-log(selected_terms$P.value), 0.95)
max_p <- 10

selected_terms$logp <- -log(selected_terms$P.value)
selected_terms$logp <- ifelse(selected_terms$logp > max_p, max_p, selected_terms$logp)

# remove Reactome Term ID
library(stringr)
selected_terms$Term <- str_replace(selected_terms$Term, "\\(GO.*", "")

selected_terms <- selected_terms %>%
  arrange(group)


selected_terms$wrap <- wrapText(selected_terms$Term, 45)

selected_terms$Term <- factor(
  as.character(selected_terms$Term),
  levels = rev(unique(as.character(selected_terms$Term)))
)

selected_terms$wrap <- factor(
  as.character(selected_terms$wrap),
  levels = rev(unique(as.character(selected_terms$wrap)))
)

library(viridis)

# Reactome Term dot plot
p <- selected_terms %>%
  ggplot(aes(x = group, y = wrap, color =logp, size=log(Combined.Score))) +
  geom_point() +
  scale_color_stepsn(colors=rev(magma(256))) +
  RotatedAxis() + xlab('') + ylab('') +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_rect(size=1, color='black', fill=NA),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(0,0,0,0),
    panel.grid = element_line(size=0.25, color='lightgrey')
  )


mapping <- labels2colors(1:100)


color_df <- data.frame(group = selected_terms$group,colour = selected_terms$color)



color_df$group <- factor(
  as.character(color_df$group),
  levels = unique(color_df$group)
)


c_vect <- color_df$colour
names(c_vect) <- color_df$group



# make the colorbar as its own heatmap
color_df$var <- 1
colorbar <- color_df %>%
  ggplot(aes(x=group, y=var, fill=group)) +
  geom_tile() +
  scale_fill_manual(values=c_vect) +
  coord_equal() +
  NoLegend() + RotatedAxis() +
  theme(
    plot.title=element_blank(),
    axis.line=element_blank(),
    axis.ticks.y =element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    plot.margin=margin(0,0,0,0),
  )


pdf(paste0('./output/Figure_8/', 'peds_CM_by_cluster_terms_cell_type_both_SystolicHF_vs_Donor.pdf'), width=6, height=9)
print(p / colorbar)
dev.off()



#Deep dive M2

bulk_modules <- read.csv("./dependencies/shared/bulk_heart_modules.csv")
bulk_modules$module <- match(bulk_modules$module,mapping)


.cache_deep_SvD_markers <- './output/Figure_8/fig8_deep_SystolicHF_vs_Donor_markers_cache.rds'
if (file.exists(.cache_deep_SvD_markers)) {
  message('Loading cached FindMarkers SystolicHF vs Donor (deep dive)...')
  combined_set <- readRDS(.cache_deep_SvD_markers)
} else {
  combined_set <- data.frame()
  mods_idx <- c(2,12,28,10,25,26)
  for (i in mods_idx){
    key_genes <- subset(bulk_modules,module %in% c(i))$gene_name
    key_genes <- key_genes[key_genes %in% rownames(M1)]
    gene_set <- FindMarkers(M1, ident.1 = "SystolicHF", ident.2 = "Donor",features=key_genes)
    gene_set<-subset(gene_set,p_val_adj<0.05)
    gene_set$module <- paste0('M',i)
    gene_set$color <- mapping[i]
    if (length(combined_set) == 0){
      combined_set <- gene_set
    }
    else {
      combined_set <- rbind(combined_set,gene_set)
    }
  }
  saveRDS(combined_set, .cache_deep_SvD_markers)
}

#cat(rownames(subset(combined_set,module=='M2' & avg_log2FC>0)),sep='\n')

M2_genes_up <- rownames(subset(combined_set,module=="M2" & avg_log2FC>0))
M10_genes_up <- rownames(subset(combined_set,module=="M10" & avg_log2FC>0))
M26_genes_up <- rownames(subset(combined_set,module=="M26" & avg_log2FC>0))
M25_genes_up <- rownames(subset(combined_set,module=="M25" & avg_log2FC>0))
M28_genes_up <- rownames(subset(combined_set,module=="M28" & avg_log2FC>0))


M2_genes_down <- rownames(subset(combined_set,module=="M2" & avg_log2FC<0))
M10_genes_down <- rownames(subset(combined_set,module=="M10" & avg_log2FC<0))
M26_genes_down <- rownames(subset(combined_set,module=="M26" & avg_log2FC<0))


dbs <- c("ChEA_2022","WikiPathway_2023_Human","Reactome_2016","GO_Biological_Process_2023")

library(enrichR)
library(forcats)

parse_ratio <- function(ratio) {
    ratio <- sub("^\\s*", "", as.character(ratio))
    ratio <- sub("\\s*$", "", ratio)
    numerator <- as.numeric(sub("/\\d+$", "", ratio))
    denominator <- as.numeric(sub("^\\d+/", "", ratio))
    return(numerator/denominator)
}

wrapText <- function(x, len) {
    sapply(x, function(y) paste(strwrap(y, len), collapse = "\n"), USE.NAMES = FALSE)
}

.cache_enrichr_SvD_M2_up <- './output/Figure_8/fig8_enrichr_SystolicHF_vs_Donor_M2_up_cache.rds'
if (file.exists(.cache_enrichr_SvD_M2_up)) {
  message('Loading cached enrichr SystolicHF vs Donor M2 up...')
  enriched <- readRDS(.cache_enrichr_SvD_M2_up)
} else {
  enriched <- enrichr(M2_genes_up, dbs)
  saveRDS(enriched, .cache_enrichr_SvD_M2_up)
}
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('./output/Figure_8/CM_Peds_M2_enrichr_up.pdf',width=5,height=2.5)
p1<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),], 
  (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), 
  size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + 
  ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + 
  ggtitle('GO Biological Process Up') + 
  scale_y_discrete(labels= fct_inorder(
    wrapText(sapply(
      strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),]$Term," \\(GO"),
         `[`, 1),35))) + 
  theme(axis.text=element_text(colour="black"))+
  scale_color_stepsn(colors=rev(magma(256)))
print(p1)
dev.off()

.cache_enrichr_SvD_M2_down <- './output/Figure_8/fig8_enrichr_SystolicHF_vs_Donor_M2_down_cache.rds'
if (file.exists(.cache_enrichr_SvD_M2_down)) {
  message('Loading cached enrichr SystolicHF vs Donor M2 down...')
  enriched <- readRDS(.cache_enrichr_SvD_M2_down)
} else {
  enriched <- enrichr(M2_genes_down, dbs)
  saveRDS(enriched, .cache_enrichr_SvD_M2_down)
}
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('./output/Figure_8/CM_Peds_M2_enrichr_down.pdf',width=5,height=2.5)
p2<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),], 
  (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), 
  size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + 
  ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + 
  ggtitle('GO Biological Process Up') + 
  scale_y_discrete(labels= fct_inorder(
    wrapText(sapply(
      strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),]$Term," \\(GO"),
         `[`, 1),35))) + 
  theme(axis.text=element_text(colour="black"))+
  scale_color_stepsn(colors=rev(magma(256)))
print(p2)
dev.off()

pdf('./output/Figure_8/CM_Peds_M2_enrichr_up_down.pdf',width=6,height=4)
print(p1/p2)
dev.off()



.cache_enrichr_SvD_mito_up <- './output/Figure_8/fig8_enrichr_SystolicHF_vs_Donor_mito_up_cache.rds'
if (file.exists(.cache_enrichr_SvD_mito_up)) {
  message('Loading cached enrichr SystolicHF vs Donor mito up...')
  enriched <- readRDS(.cache_enrichr_SvD_mito_up)
} else {
  enriched <- enrichr(c(M10_genes_up,M25_genes_up,M26_genes_up,M28_genes_up), dbs)
  saveRDS(enriched, .cache_enrichr_SvD_mito_up)
}
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('./output/Figure_8/CM_Peds_mito_enrichr_up.pdf',width=6,height=3)
p1<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:5),], 
  (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), 
  size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + 
  ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + 
  ggtitle('GO Biological Process Up') + 
  scale_y_discrete(labels= fct_inorder(
    wrapText(sapply(
      strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:5),]$Term," \\(GO"),
         `[`, 1),35))) + 
  theme(axis.text=element_text(colour="black"))+
  scale_color_stepsn(colors=rev(magma(256)))
print(p1)
dev.off()




###PRV


bulk_modules <- read.csv("./dependencies/shared/bulk_heart_modules.csv")
bulk_modules$module <- match(bulk_modules$module,mapping)


.cache_deep_SvN_markers <- './output/Figure_8/fig8_deep_SystolicHF_vs_NF_markers_cache.rds'
if (file.exists(.cache_deep_SvN_markers)) {
  message('Loading cached FindMarkers SystolicHF vs NF (deep dive)...')
  combined_set <- readRDS(.cache_deep_SvN_markers)
} else {
  combined_set <- data.frame()
  mods_idx <- c(2,12,28,10,25,26)
  for (i in mods_idx){
    key_genes <- subset(bulk_modules,module %in% c(i))$gene_name
    key_genes <- key_genes[key_genes %in% rownames(M1)]
    gene_set <- FindMarkers(M1, ident.1 = "SystolicHF", ident.2 = "NF",features=key_genes)
    gene_set<-subset(gene_set,p_val_adj<0.05)
    gene_set$module <- paste0('M',i)
    gene_set$color <- mapping[i]
    if (length(combined_set) == 0){
      combined_set <- gene_set
    }
    else {
      combined_set <- rbind(combined_set,gene_set)
    }
  }
  saveRDS(combined_set, .cache_deep_SvN_markers)
}

#cat(rownames(subset(combined_set,module=='M2' & avg_log2FC>0)),sep='\n')

M2_genes_up <- rownames(subset(combined_set,module=="M2" & avg_log2FC>0))
M10_genes_up <- rownames(subset(combined_set,module=="M10" & avg_log2FC>0))
M26_genes_up <- rownames(subset(combined_set,module=="M26" & avg_log2FC>0))
M25_genes_up <- rownames(subset(combined_set,module=="M25" & avg_log2FC>0))
M28_genes_up <- rownames(subset(combined_set,module=="M28" & avg_log2FC>0))


M2_genes_down <- rownames(subset(combined_set,module=="M2" & avg_log2FC<0))
M10_genes_down <- rownames(subset(combined_set,module=="M10" & avg_log2FC<0))
M26_genes_down <- rownames(subset(combined_set,module=="M26" & avg_log2FC<0))


dbs <- c("ChEA_2022","WikiPathway_2023_Human","Reactome_2016","GO_Biological_Process_2023")

library(enrichR)
library(forcats)

parse_ratio <- function(ratio) {
    ratio <- sub("^\\s*", "", as.character(ratio))
    ratio <- sub("\\s*$", "", ratio)
    numerator <- as.numeric(sub("/\\d+$", "", ratio))
    denominator <- as.numeric(sub("^\\d+/", "", ratio))
    return(numerator/denominator)
}

wrapText <- function(x, len) {
    sapply(x, function(y) paste(strwrap(y, len), collapse = "\n"), USE.NAMES = FALSE)
}

.cache_enrichr_SvN_M2_up <- './output/Figure_8/fig8_enrichr_SystolicHF_vs_NF_M2_up_cache.rds'
if (file.exists(.cache_enrichr_SvN_M2_up)) {
  message('Loading cached enrichr SystolicHF vs NF M2 up...')
  enriched <- readRDS(.cache_enrichr_SvN_M2_up)
} else {
  enriched <- enrichr(M2_genes_up, dbs)
  saveRDS(enriched, .cache_enrichr_SvN_M2_up)
}
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('./output/Figure_8/CM_Peds_M2_enrichr_up_RVF_vs_pRV.pdf',width=5,height=2.5)
p1<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),], 
  (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), 
  size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + 
  ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + 
  ggtitle('GO Biological Process Up') + 
  scale_y_discrete(labels= fct_inorder(
    wrapText(sapply(
      strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),]$Term," \\(GO"),
         `[`, 1),35))) + 
  theme(axis.text=element_text(colour="black"))+
  scale_color_stepsn(colors=rev(magma(256)))
print(p1)
dev.off()

.cache_enrichr_SvN_M2_down <- './output/Figure_8/fig8_enrichr_SystolicHF_vs_NF_M2_down_cache.rds'
if (file.exists(.cache_enrichr_SvN_M2_down)) {
  message('Loading cached enrichr SystolicHF vs NF M2 down...')
  enriched <- readRDS(.cache_enrichr_SvN_M2_down)
} else {
  enriched <- enrichr(M2_genes_down, dbs)
  saveRDS(enriched, .cache_enrichr_SvN_M2_down)
}
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('./output/Figure_8/CM_Peds_M2_enrichr_down_RVF_vs_pRV.pdf',width=5,height=2.5)
p2<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),], 
  (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), 
  size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + 
  ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + 
  ggtitle('GO Biological Process Up') + 
  scale_y_discrete(labels= fct_inorder(
    wrapText(sapply(
      strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),]$Term," \\(GO"),
         `[`, 1),35))) + 
  theme(axis.text=element_text(colour="black"))+
  scale_color_stepsn(colors=rev(magma(256)))
print(p2)
dev.off()

pdf('./output/Figure_8/CM_Peds_M2_enrichr_up_down_RVF_vs_pRV.pdf',width=6,height=4)
p_8F_top <- p1/p2
print(p_8F_top)
dev.off()



.cache_enrichr_SvN_mito_up <- './output/Figure_8/fig8_enrichr_SystolicHF_vs_NF_mito_up_cache.rds'
if (file.exists(.cache_enrichr_SvN_mito_up)) {
  message('Loading cached enrichr SystolicHF vs NF mito up...')
  enriched <- readRDS(.cache_enrichr_SvN_mito_up)
} else {
  enriched <- enrichr(c(M10_genes_up,M25_genes_up,M26_genes_up,M28_genes_up), dbs)
  saveRDS(enriched, .cache_enrichr_SvN_mito_up)
}
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('./output/Figure_8/CM_Peds_mito_enrichr_up_RVF_vs_pRV.pdf',width=6,height=3)
p1<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:5),], 
  (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), 
  size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + 
  ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + 
  ggtitle('GO Biological Process Up') + 
  scale_y_discrete(labels= fct_inorder(
    wrapText(sapply(
      strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:5),]$Term," \\(GO"),
         `[`, 1),35))) + 
  theme(axis.text=element_text(colour="black"))+
  scale_color_stepsn(colors=rev(magma(256)))
print(p1)
dev.off()



###SV vs NF


bulk_modules <- read.csv("./dependencies/shared/bulk_heart_modules.csv")
bulk_modules$module <- match(bulk_modules$module,mapping)


.cache_deep_NvD_markers <- './output/Figure_8/fig8_deep_NF_vs_Donor_markers_cache.rds'
if (file.exists(.cache_deep_NvD_markers)) {
  message('Loading cached FindMarkers NF vs Donor (deep dive)...')
  combined_set <- readRDS(.cache_deep_NvD_markers)
} else {
  combined_set <- data.frame()
  mods_idx <- c(2,12,28,10,25,26)
  for (i in mods_idx){
    key_genes <- subset(bulk_modules,module %in% c(i))$gene_name
    key_genes <- key_genes[key_genes %in% rownames(M1)]
    gene_set <- FindMarkers(M1, ident.1 = "NF", ident.2 = "Donor",features=key_genes)
    gene_set<-subset(gene_set,p_val_adj<0.05)
    gene_set$module <- paste0('M',i)
    gene_set$color <- mapping[i]
    if (length(combined_set) == 0){
      combined_set <- gene_set
    }
    else {
      combined_set <- rbind(combined_set,gene_set)
    }
  }
  saveRDS(combined_set, .cache_deep_NvD_markers)
}

#cat(rownames(subset(combined_set,module=='M2' & avg_log2FC>0)),sep='\n')

M2_genes_up <- rownames(subset(combined_set,module=="M2" & avg_log2FC>0))
M10_genes_up <- rownames(subset(combined_set,module=="M10" & avg_log2FC>0))
M26_genes_up <- rownames(subset(combined_set,module=="M26" & avg_log2FC>0))
M25_genes_up <- rownames(subset(combined_set,module=="M25" & avg_log2FC>0))
M28_genes_up <- rownames(subset(combined_set,module=="M28" & avg_log2FC>0))


M2_genes_down <- rownames(subset(combined_set,module=="M2" & avg_log2FC<0))
M10_genes_down <- rownames(subset(combined_set,module=="M10" & avg_log2FC<0))
M26_genes_down <- rownames(subset(combined_set,module=="M26" & avg_log2FC<0))


dbs <- c("ChEA_2022","WikiPathway_2023_Human","Reactome_2016","GO_Biological_Process_2023")

library(enrichR)
library(forcats)

parse_ratio <- function(ratio) {
    ratio <- sub("^\\s*", "", as.character(ratio))
    ratio <- sub("\\s*$", "", ratio)
    numerator <- as.numeric(sub("/\\d+$", "", ratio))
    denominator <- as.numeric(sub("^\\d+/", "", ratio))
    return(numerator/denominator)
}

wrapText <- function(x, len) {
    sapply(x, function(y) paste(strwrap(y, len), collapse = "\n"), USE.NAMES = FALSE)
}

.cache_enrichr_NvD_M2_up <- './output/Figure_8/fig8_enrichr_NF_vs_Donor_M2_up_cache.rds'
if (file.exists(.cache_enrichr_NvD_M2_up)) {
  message('Loading cached enrichr NF vs Donor M2 up...')
  enriched <- readRDS(.cache_enrichr_NvD_M2_up)
} else {
  enriched <- enrichr(M2_genes_up, dbs)
  saveRDS(enriched, .cache_enrichr_NvD_M2_up)
}
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('./output/Figure_8/CM_Peds_M2_enrichr_up_pRV_vs_NF.pdf',width=5,height=2.5)
p1<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),], 
  (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), 
  size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + 
  ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + 
  ggtitle('GO Biological Process Up') + 
  scale_y_discrete(labels= fct_inorder(
    wrapText(sapply(
      strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),]$Term," \\(GO"),
         `[`, 1),35))) + 
  theme(axis.text=element_text(colour="black"))+
  scale_color_stepsn(colors=rev(magma(256)))
print(p1)
dev.off()

.cache_enrichr_NvD_M2_down <- './output/Figure_8/fig8_enrichr_NF_vs_Donor_M2_down_cache.rds'
if (file.exists(.cache_enrichr_NvD_M2_down)) {
  message('Loading cached enrichr NF vs Donor M2 down...')
  enriched <- readRDS(.cache_enrichr_NvD_M2_down)
} else {
  enriched <- enrichr(M2_genes_down, dbs)
  saveRDS(enriched, .cache_enrichr_NvD_M2_down)
}
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('./output/Figure_8/CM_Peds_M2_enrichr_down_pRV_vs_NF.pdf',width=5,height=2.5)
p2<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),], 
  (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), 
  size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + 
  ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + 
  ggtitle('GO Biological Process Up') + 
  scale_y_discrete(labels= fct_inorder(
    wrapText(sapply(
      strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:3),]$Term," \\(GO"),
         `[`, 1),35))) + 
  theme(axis.text=element_text(colour="black"))+
  scale_color_stepsn(colors=rev(magma(256)))
print(p2)
dev.off()

pdf('./output/Figure_8/CM_Peds_M2_enrichr_up_down_pRV_vs_NF.pdf',width=6,height=4)
print(p1/p2)
dev.off()
## (v57 Panel F is assembled after the NvD pooled-mito-up block below.)



.cache_enrichr_NvD_mito_up <- './output/Figure_8/fig8_enrichr_NF_vs_Donor_mito_up_cache.rds'
if (file.exists(.cache_enrichr_NvD_mito_up)) {
  message('Loading cached enrichr NF vs Donor mito up...')
  enriched <- readRDS(.cache_enrichr_NvD_mito_up)
} else {
  enriched <- enrichr(c(M10_genes_up,M25_genes_up,M26_genes_up,M28_genes_up), dbs)
  saveRDS(enriched, .cache_enrichr_NvD_mito_up)
}
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('./output/Figure_8/CM_Peds_mito_enrichr_up_pRV_vs_NF.pdf',width=6,height=3)
p1<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:5),],
  (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value),
  size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') +
  ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  +
  ggtitle('GO Biological Process Up') +
  scale_y_discrete(labels= fct_inorder(
    wrapText(sapply(
      strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:5),]$Term," \\(GO"),
         `[`, 1),35))) +
  theme(axis.text=element_text(colour="black"))+
  scale_color_stepsn(colors=rev(magma(256)))
print(p1)
dev.off()

## ── v57 Panel F: GO-BP enrichment of UPregulated pooled M10/M25/M26/M28
## DEGs — (top) ped-RVF vs ped-HLHS (SystolicHF vs NF); (bottom)
## ped-HLHS vs ped-NF (NF vs Donor). Built from the two cached enrichr
## results so it is independent of variable-overwrite order above.
.mk_8F <- function(rds, ttl) {
  e <- readRDS(rds)[[4]]
  e <- subset(e, Adjusted.P.value < 0.05)
  e <- e[order(e$Combined.Score, decreasing = TRUE), , drop = FALSE]
  e <- e[rev(seq_len(min(5, nrow(e)))), , drop = FALSE]
  ggplot(e, aes(x = Combined.Score, y = fct_inorder(Term),
                color = as.numeric(Adjusted.P.value),
                size = parse_ratio(Overlap))) +
    geom_point() + xlab('Combined Score') + ylab('Term') +
    labs(color = 'P value', size = 'Overlap') + ggtitle(ttl) +
    scale_y_discrete(labels = fct_inorder(
      wrapText(sapply(strsplit(e$Term, ' \\(GO'), `[`, 1), 35))) +
    theme_classic() + theme(axis.text = element_text(colour = 'black')) +
    scale_color_stepsn(colors = rev(magma(256)))
}
p_8F_top <- .mk_8F(.cache_enrichr_SvN_mito_up,
                   'ped-RVF vs ped-HLHS: pooled M10/M25/M26/M28 up (GO-BP)')
p_8F_bot <- .mk_8F(.cache_enrichr_NvD_mito_up,
                   'ped-HLHS vs ped-NF: pooled M10/M25/M26/M28 up (GO-BP)')
p_8F <- p_8F_top / p_8F_bot
}   # ===== end disabled OLD enrichR F8F =====


###############################################################################
##  FIGURE 8F (v67): pediatric CM MitoCarta volcanoes — rotated, stacked.
##  Per-patient pseudobulk from M1 (CM object) -> DESeq2 -> apeglm shrink (padj
##  from unshrunken Wald). MitoCarta3.0 genes coloured by top-level category.
##  Top = ped-HLHS (condition NF) vs Donor; bottom = ped-RVF (SystolicHF) vs Donor.
###############################################################################
suppressMessages({library(DESeq2); library(apeglm); library(readxl); library(stringr)
  library(ggrepel); library(patchwork); library(fgsea)})

M1$.pat   <- as.character(M1$sample)
.f8f_pb   <- tryCatch(AggregateExpression(M1, group.by='.pat', assays='RNA', slot='counts', return.seurat=FALSE),
                      error=function(e) AggregateExpression(M1, group.by='.pat', assays='RNA', return.seurat=FALSE))
.f8f_mat  <- round(as.matrix(.f8f_pb$RNA)); colnames(.f8f_mat) <- sub('^g','',colnames(.f8f_mat))
.f8f_meta <- unique(data.frame(patient=as.character(M1$sample), condition=as.character(M1$condition)))
.f8f_meta <- .f8f_meta[match(colnames(.f8f_mat), .f8f_meta$patient), ]; rownames(.f8f_meta) <- .f8f_meta$patient
.f8f_meta$condition <- factor(.f8f_meta$condition, levels=c('Donor','NF','SystolicHF'))
.f8f_dds  <- DESeq(DESeqDataSetFromMatrix(.f8f_mat, .f8f_meta, ~condition))

.f8f_Hmito <- read_excel('./dependencies/shared/Human.MitoCarta3.0.xls', sheet='A Human MitoCarta3.0')
.f8f_a <- trimws(unlist(lapply(lapply(lapply(.f8f_Hmito$MitoCarta3.0_MitoPathways, str_split, '>'), '[[', 1), '[[', 1)))
.f8f_a[.f8f_a=='Small molecule transport | Signaling'] <- 'Small molecule transport'
.f8f_remap <- function(x){x<-as.character(x);o<-character(length(x))
  o[is.na(x)|x=='']<-'No annotation';o[x=='Metabolism']<-'Metabolism';o[x=='OXPHOS']<-'Oxidative phos'
  o[x=='Protein homeostasis']<-'Protein import/sorting';o[x=='Mitochondrial central dogma']<-'Mito central dogma'
  o[x=='Mitochondrial dynamics and surveillance']<-'Mito dynamics';o[x=='Signaling']<-'Signaling'
  o[x=='Small molecule transport']<-'Transport';o[o=='']<-'No annotation';o}
.f8f_anno <- setNames(.f8f_remap(.f8f_a), .f8f_Hmito$Symbol)
.f8f_cols <- c('Oxidative phos'='#1B9E77','Metabolism'='#D95F02','No annotation'='#7570B3',
  'Mito central dogma'='#E7298A','Transport'='#66A61E','Signaling'='#E6AB02',
  'Mito dynamics'='#A6761D','Protein import/sorting'='#666666')
.f8f_hsym <- unique(.f8f_Hmito$Symbol); .f8f_nuc <- .f8f_hsym[!grepl('^MT-', .f8f_hsym, ignore.case=TRUE)]

.f8f_df <- function(coef, num){
  raw <- as.data.frame(results(.f8f_dds, contrast=c('condition', num, 'Donor')))
  shr <- as.data.frame(lfcShrink(.f8f_dds, coef=coef, type='apeglm')); shr$padj <- raw[rownames(shr),'padj']
  rk <- raw$stat; names(rk) <- rownames(raw); rk <- rk[!is.na(rk)]; rk <- rk[!duplicated(names(rk))]
  set.seed(1); nes <- fgsea::fgsea(list(M=intersect(.f8f_nuc, names(rk))), rk, nPermSimple=10000)$NES
  g <- intersect(.f8f_hsym, rownames(shr))
  d <- data.frame(gene=g, log2fc=shr[g,'log2FoldChange'], padj=shr[g,'padj'], category=.f8f_anno[g])
  d <- d[complete.cases(d[,c('log2fc','padj')]),]; d$neglogp <- pmin(-log10(d$padj), 6); list(d=d, nes=nes)
}
.f8f_H <- .f8f_df('condition_NF_vs_Donor', 'NF')
.f8f_R <- .f8f_df('condition_SystolicHF_vs_Donor', 'SystolicHF')
.f8f_yl <- max(4, ceiling(max(abs(c(.f8f_H$d$log2fc, .f8f_R$d$log2fc)))))
.f8f_FS <- 11.25
.f8f_volc <- function(x, title, nes, is_bottom){
  lab <- head(x[x$padj<0.05, ][order(x[x$padj<0.05,'padj']), ], 7)
  g <- ggplot(x, aes(neglogp, log2fc, color=category)) +
    geom_hline(yintercept=0, color='grey80') +
    geom_hline(yintercept=c(-0.5,0.5), linetype='dashed', color='grey60') +
    geom_vline(xintercept=-log10(0.05), linetype='dashed', color='grey60') +
    geom_point(size=1.6, alpha=0.85) +
    geom_text_repel(data=lab, aes(label=gene), size=.f8f_FS/.pt, fontface='italic', max.overlaps=20, show.legend=FALSE) +
    scale_color_manual(values=.f8f_cols, drop=FALSE, name=NULL) +
    scale_x_continuous(limits=c(0,6.3), expand=expansion(mult=c(0,0.02))) +
    scale_y_continuous(limits=c(-.f8f_yl, .f8f_yl)) +
    labs(title=sprintf('%s  (NES %.2f, n.s.)', title, nes),
         x=if(is_bottom) expression(-log[10]~adj.~P) else NULL, y=expression(log[2]~fold~change)) +
    theme_bw(base_size=.f8f_FS) +
    theme(text=element_text(size=.f8f_FS), plot.title=element_text(size=.f8f_FS),
          axis.title=element_text(size=.f8f_FS), axis.text=element_text(size=.f8f_FS),
          legend.text=element_text(size=.f8f_FS), panel.grid.minor=element_blank(),
          plot.margin=if(is_bottom) margin(0,3,3,3) else margin(3,3,0,3))
  if(!is_bottom) g <- g + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
  g
}
p_8F <- (.f8f_volc(.f8f_H$d, 'HLHS vs Donor', .f8f_H$nes, FALSE) /
         .f8f_volc(.f8f_R$d, 'RVF vs Donor',  .f8f_R$nes, TRUE)) +
        plot_layout(guides='collect') &
        theme(legend.position='right') & guides(color=guide_legend(ncol=1, override.aes=list(size=3)))
save_figure(p_8F, 'Figure_8_panel_F.pdf', width = 6.5, height = 4.3)

















# ##### OLD


# M3 <- subset(M2,CombinedNames=='CM')
# DefaultAssay(M3) <- 'RNA'
# M3[["RNA"]] <- split(M3[["RNA"]], f = M3$patient)


# M3 <- NormalizeData(M3)
# M3 <- FindVariableFeatures(M3)
# M3 <- ScaleData(M3)
# M3 <- RunPCA(M3)


# M3 <- IntegrateLayers(
#   object = M3, method = HarmonyIntegration,
#   orig.reduction = "pca", new.reduction = "integrated.cca",
#   verbose = FALSE
# )

# M3 <- FindNeighbors(M3, reduction = "integrated.cca", dims = 1:50)
# M3 <- FindClusters(M3, resolution = 2, cluster.name = "cca_clusters")
# M3 <- RunUMAP(M3, reduction = "integrated.cca", dims = 1:50, reduction.name = "umap.cca")

# p1 <- DimPlot(
#   M3,
#   reduction = "umap.cca",
#   group.by = c("origin"),
#   combine = FALSE, label.size = 2
# )


# pdf(paste0('./output/Figure_8/', 'sn_RV_Peds_CM_Harmonized_UMAP.pdf'), width=5, height=5)
# PlotEmbedding(M3,group.by='origin',reduction="umap.cca",point_size=0.2,plot_under=TRUE,plot_theme=umap_theme()+NoLegend(),raster_dpi=400,raster_scale=0.5)
# dev.off()


# M3 <- subset(M2,CombinedNames=='Myeloid')
# DefaultAssay(M3) <- 'RNA'
# M3[["RNA"]] <- split(M3[["RNA"]], f = M3$patient)


# M3 <- NormalizeData(M3)
# M3 <- FindVariableFeatures(M3)
# M3 <- ScaleData(M3)
# M3 <- RunPCA(M3)


# M3 <- IntegrateLayers(
#   object = M3, method = HarmonyIntegration,
#   orig.reduction = "pca", new.reduction = "integrated.cca",
#   verbose = FALSE
# )

# M3 <- FindNeighbors(M3, reduction = "integrated.cca", dims = 1:50)
# M3 <- FindClusters(M3, resolution = 2, cluster.name = "cca_clusters")
# M3 <- RunUMAP(M3, reduction = "integrated.cca", dims = 1:50, reduction.name = "umap.cca")

# p1 <- DimPlot(
#   M3,
#   reduction = "umap.cca",
#   group.by = c("origin"),
#   combine = FALSE, label.size = 2
# )


# pdf(paste0('./output/Figure_8/', 'sn_RV_Peds_Myeloid_Harmonized_UMAP.pdf'), width=5, height=5)
# PlotEmbedding(M3,group.by='origin',reduction="umap.cca",point_size=0.2,plot_under=TRUE,plot_theme=umap_theme()+NoLegend(),raster_dpi=400,raster_scale=0.5)
# dev.off()



















# DefaultAssay(M3) <- 'RNA'
# M3 <- NormalizeData(M3)
# M3 <- FindVariableFeatures(M3)
# M3 <- ScaleData(M3)
# M3 <- RunPCA(M3)
# M3 <- RunHarmony(M3,'origin')
# M3 <- RunUMAP(M3,reduction = "harmony", dims = 1:50, verbose = F) 

# RV_and_Peds <- intersect(rownames(subset(M3,origin=='RV')),rownames(subset(M3,origin=='Peds')))



# M3 <- SetIdent(M3, value = "origin")

# cm_compare <- FindMarkers(M3,ident.1='RV',ident.2='Peds')


# write.csv(cm_compare,'./output/Figure_8/cm_RV_vs_Peds_DEG.csv')



# M3 <- subset(M2,CombinedNames=='Myeloid')
# DefaultAssay(M3) <- 'RNA'
# M3 <- SetIdent(M3, value = "origin")
# myeloid_compare <- FindMarkers(M3,ident.1='RV',ident.2='Peds')









# RefMerge<-readRDS('./output/Figure_8/Kory_with_RV_modules_projected.rds')




# ####CM
# RefMerge<-SetActiveWGCNA(RefMerge, "CM_RV2LV")

# #RefMerge <- ModuleConnectivity(RefMerge,group.by = 'Names', group_name = 'Cardiomyocytes')


# modules <- GetModules(RefMerge)
# MEs <- GetMEs(RefMerge, T)
# genes_use <- as.character(modules$gene_name)
# params <- GetWGCNAParams(RefMerge)

# # get the assay
# assay <- DefaultAssay(RefMerge)

# cells.use <- RefMerge@meta.data %>% subset(get('Names') %in% 'Cardiomyocytes') %>% rownames
# MEs <- MEs[cells.use,]

# exp_mat <- SeuratObject::LayerData(RefMerge, assay=assay, layer='data')[genes_use,cells.use]

# kMEs <- corSparse(X = Matrix::t(exp_mat),Y = as.matrix(MEs))
# rownames(kMEs) <- genes_use
# kMEs <- as.data.frame(kMEs)

# modules <- modules[,1:3]
# mods <- levels(modules$module)
# colnames(kMEs) <- colnames(MEs)
# kMEs <- kMEs[,mods]
# colnames(kMEs) <- paste0("kME_", colnames(kMEs))
# kMEs <- cbind(modules, kMEs)
# RefMerge <- SetModules(RefMerge, kMEs, "CM_RV2LV")



# RefMerge <- ModuleExprScore(RefMerge,method='Seurat')


# mod_scores <-  GetModuleScores(RefMerge)
# mod_scores <- mod_scores[,colnames(mod_scores) != 'grey']

# # add hMEs to Seurat meta-data:
# RefMerge@meta.data <- cbind(
#   RefMerge@meta.data,
#   mod_scores
# )

# # plot with Seurat's DotPlot function
# p1 <- DotPlot(
#     RefMerge,
#     features = mixedsort(colnames(mod_scores)),
#     group.by = 'Names'
# )

# # flip the x/y axes, rotate the axis labels, and change color scheme:
# p1 <- p1 +
#   RotatedAxis() +
#   scale_color_gradient2(high='darkorchid4', mid='grey95', low='seagreen') +
#   xlab('') + ylab('')


# pdf(paste0('./output/Figure_8/', 'RV2LV_LV_CM_modules_dot.pdf'), width=7.5, height=5)
# p1
# dev.off()


# seurat_ref <- readRDS('./output/Figure_8/scWGCNA_all_celltypes.rds')
# seurat_ref<-SetActiveWGCNA(seurat_ref, "CM")
# seurat_ref <- ModuleConnectivity(seurat_ref,group.by = 'Names', group_name = 'CM')

# seurat_ref <- ModuleExprScore(seurat_ref,method='Seurat')


# mod_scores <-  GetModuleScores(seurat_ref)
# mod_scores <- mod_scores[,colnames(mod_scores) != 'grey']

# # add hMEs to Seurat meta-data:
# seurat_ref@meta.data <- cbind(
#   seurat_ref@meta.data,
#   mod_scores
# )

# # plot with Seurat's DotPlot function
# p2 <- DotPlot(
#     seurat_ref,
#     features = mixedsort(colnames(mod_scores)),
#     group.by = 'Names'
# )

# # flip the x/y axes, rotate the axis labels, and change color scheme:
# p2 <- p2 +
#   RotatedAxis() +
#   scale_color_gradient2(high='darkorchid4', mid='grey95', low='seagreen') +
#   xlab('') + ylab('')


# pdf(paste0('./output/Figure_8/', "LV_CM_modules_dot.pdf"), width=7.5, height=5)
# p2
# dev.off()

# pdf(paste0('./output/Figure_8/', "LV_CM_LV2RM_modules_dot.pdf"), width=7.5, height=8)
# wrap_plots(list(p1,p2),ncol=1)
# dev.off()

# ###EC

# RefMerge<-SetActiveWGCNA(RefMerge, "EC_RV2LV")

# #RefMerge <- ModuleConnectivity(RefMerge,group.by = 'Names', group_name = 'Cardiomyocytes')


# modules <- GetModules(RefMerge)
# MEs <- GetMEs(RefMerge, T)
# genes_use <- as.character(modules$gene_name)
# params <- GetWGCNAParams(RefMerge)

# # get the assay
# assay <- DefaultAssay(RefMerge)

# cells.use <- RefMerge@meta.data %>% subset(get('Names') %in% 'Endothelium') %>% rownames
# MEs <- MEs[cells.use,]

# exp_mat <- SeuratObject::LayerData(RefMerge, assay=assay, layer='data')[genes_use,cells.use]

# kMEs <- corSparse(X = Matrix::t(exp_mat),Y = as.matrix(MEs))
# rownames(kMEs) <- genes_use
# kMEs <- as.data.frame(kMEs)

# modules <- modules[,1:3]
# mods <- levels(modules$module)
# colnames(kMEs) <- colnames(MEs)
# kMEs <- kMEs[,mods]
# colnames(kMEs) <- paste0("kME_", colnames(kMEs))
# kMEs <- cbind(modules, kMEs)
# RefMerge <- SetModules(RefMerge, kMEs, "EC_RV2LV")



# RefMerge <- ModuleExprScore(RefMerge,method='Seurat')


# mod_scores <-  GetModuleScores(RefMerge)
# mod_scores <- mod_scores[,colnames(mod_scores) != 'grey']

# # add hMEs to Seurat meta-data:
# RefMerge@meta.data <- cbind(
#   RefMerge@meta.data,
#   mod_scores
# )

# # plot with Seurat's DotPlot function
# p1 <- DotPlot(
#     RefMerge,
#     features = mixedsort(colnames(mod_scores)),
#     group.by = 'Names'
# )

# # flip the x/y axes, rotate the axis labels, and change color scheme:
# p1 <- p1 +
#   RotatedAxis() +
#   scale_color_gradient2(high='darkorchid4', mid='grey95', low='seagreen') +
#   xlab('') + ylab('')


# pdf(paste0('./output/Figure_8/', 'RV2LV_LV_EC_modules_dot.pdf'), width=7.5, height=5)
# p1
# dev.off()


# seurat_ref<-SetActiveWGCNA(seurat_ref, "EC")
# seurat_ref <- ModuleConnectivity(seurat_ref,group.by = 'Names', group_name = 'EC')

# seurat_ref <- ModuleExprScore(seurat_ref,method='Seurat')


# mod_scores <-  GetModuleScores(seurat_ref)
# mod_scores <- mod_scores[,colnames(mod_scores) != 'grey']

# # add hMEs to Seurat meta-data:
# seurat_ref@meta.data <- cbind(
#   seurat_ref@meta.data,
#   mod_scores
# )

# # plot with Seurat's DotPlot function
# p2 <- DotPlot(
#     seurat_ref,
#     features = mixedsort(colnames(mod_scores)),
#     group.by = 'Names'
# )

# # flip the x/y axes, rotate the axis labels, and change color scheme:
# p2 <- p2 +
#   RotatedAxis() +
#   scale_color_gradient2(high='darkorchid4', mid='grey95', low='seagreen') +
#   xlab('') + ylab('')


# pdf(paste0('./output/Figure_8/', "LV_EC_modules_dot.pdf"), width=7.5, height=5)
# p2
# dev.off()

# pdf(paste0('./output/Figure_8/', "LV_EC_LV2RM_modules_dot.pdf"), width=7.5, height=8)
# wrap_plots(list(p1,p2),ncol=1)
# dev.off()



# ###FB

# RefMerge<-SetActiveWGCNA(RefMerge, "FB_RV2LV")

# #RefMerge <- ModuleConnectivity(RefMerge,group.by = 'Names', group_name = 'Cardiomyocytes')


# modules <- GetModules(RefMerge)
# MEs <- GetMEs(RefMerge, T)
# genes_use <- as.character(modules$gene_name)
# params <- GetWGCNAParams(RefMerge)

# # get the assay
# assay <- DefaultAssay(RefMerge)

# cells.use <- RefMerge@meta.data %>% subset(get('Names') %in% 'Fibroblasts') %>% rownames
# MEs <- MEs[cells.use,]

# exp_mat <- SeuratObject::LayerData(RefMerge, assay=assay, layer='data')[genes_use,cells.use]

# kMEs <- corSparse(X = Matrix::t(exp_mat),Y = as.matrix(MEs))
# rownames(kMEs) <- genes_use
# kMEs <- as.data.frame(kMEs)

# modules <- modules[,1:3]
# mods <- levels(modules$module)
# colnames(kMEs) <- colnames(MEs)
# kMEs <- kMEs[,mods]
# colnames(kMEs) <- paste0("kME_", colnames(kMEs))
# kMEs <- cbind(modules, kMEs)
# RefMerge <- SetModules(RefMerge, kMEs, "FB_RV2LV")



# RefMerge <- ModuleExprScore(RefMerge,method='Seurat')


# mod_scores <-  GetModuleScores(RefMerge)
# mod_scores <- mod_scores[,colnames(mod_scores) != 'grey']

# # add hMEs to Seurat meta-data:
# RefMerge@meta.data <- cbind(
#   RefMerge@meta.data,
#   mod_scores
# )

# # plot with Seurat's DotPlot function
# p1 <- DotPlot(
#     RefMerge,
#     features = mixedsort(colnames(mod_scores)),
#     group.by = 'Names'
# )

# # flip the x/y axes, rotate the axis labels, and change color scheme:
# p1 <- p1 +
#   RotatedAxis() +
#   scale_color_gradient2(high='darkorchid4', mid='grey95', low='seagreen') +
#   xlab('') + ylab('')


# pdf(paste0('./output/Figure_8/', 'RV2LV_LV_FB_modules_dot.pdf'), width=7.5, height=5)
# p1
# dev.off()


# seurat_ref<-SetActiveWGCNA(seurat_ref, "FB")
# seurat_ref <- ModuleConnectivity(seurat_ref,group.by = 'Names', group_name = 'FB')

# seurat_ref <- ModuleExprScore(seurat_ref,method='Seurat')


# mod_scores <-  GetModuleScores(seurat_ref)
# mod_scores <- mod_scores[,colnames(mod_scores) != 'grey']

# # add hMEs to Seurat meta-data:
# seurat_ref@meta.data <- cbind(
#   seurat_ref@meta.data,
#   mod_scores
# )

# # plot with Seurat's DotPlot function
# p2 <- DotPlot(
#     seurat_ref,
#     features = mixedsort(colnames(mod_scores)),
#     group.by = 'Names'
# )

# # flip the x/y axes, rotate the axis labels, and change color scheme:
# p2 <- p2 +
#   RotatedAxis() +
#   scale_color_gradient2(high='darkorchid4', mid='grey95', low='seagreen') +
#   xlab('') + ylab('')


# pdf(paste0('./output/Figure_8/', "LV_FB_modules_dot.pdf"), width=7.5, height=5)
# p2
# dev.off()

# pdf(paste0('./output/Figure_8/', "LV_FB_LV2RM_modules_dot.pdf"), width=7.5, height=8)
# wrap_plots(list(p1,p2),ncol=1)
# dev.off()



# ###Myeloid

# RefMerge<-SetActiveWGCNA(RefMerge, "Myeloid_RV2LV")

# #RefMerge <- ModuleConnectivity(RefMerge,group.by = 'Names', group_name = 'Cardiomyocytes')


# modules <- GetModules(RefMerge)
# MEs <- GetMEs(RefMerge, T)
# genes_use <- as.character(modules$gene_name)
# params <- GetWGCNAParams(RefMerge)

# # get the assay
# assay <- DefaultAssay(RefMerge)

# cells.use <- RefMerge@meta.data %>% subset(get('Names') %in% 'Myeloid') %>% rownames
# MEs <- MEs[cells.use,]

# exp_mat <- SeuratObject::LayerData(RefMerge, assay=assay, layer='data')[genes_use,cells.use]

# kMEs <- corSparse(X = Matrix::t(exp_mat),Y = as.matrix(MEs))
# rownames(kMEs) <- genes_use
# kMEs <- as.data.frame(kMEs)

# modules <- modules[,1:3]
# mods <- levels(modules$module)
# colnames(kMEs) <- colnames(MEs)
# kMEs <- kMEs[,mods]
# colnames(kMEs) <- paste0("kME_", colnames(kMEs))
# kMEs <- cbind(modules, kMEs)
# RefMerge <- SetModules(RefMerge, kMEs, "Myeloid_RV2LV")



# RefMerge <- ModuleExprScore(RefMerge,method='Seurat')


# mod_scores <-  GetModuleScores(RefMerge)
# mod_scores <- mod_scores[,colnames(mod_scores) != 'grey']

# # add hMEs to Seurat meta-data:
# RefMerge@meta.data <- cbind(
#   RefMerge@meta.data,
#   mod_scores
# )

# # plot with Seurat's DotPlot function
# p1 <- DotPlot(
#     RefMerge,
#     features = mixedsort(colnames(mod_scores)),
#     group.by = 'Names'
# )

# # flip the x/y axes, rotate the axis labels, and change color scheme:
# p1 <- p1 +
#   RotatedAxis() +
#   scale_color_gradient2(high='darkorchid4', mid='grey95', low='seagreen') +
#   xlab('') + ylab('')


# pdf(paste0('./output/Figure_8/', 'RV2LV_LV_Myeloid_modules_dot.pdf'), width=7.5, height=5)
# p1
# dev.off()


# seurat_ref<-SetActiveWGCNA(seurat_ref, "Myeloid")
# seurat_ref <- ModuleConnectivity(seurat_ref,group.by = 'Names', group_name = 'Myeloid')

# seurat_ref <- ModuleExprScore(seurat_ref,method='Seurat')


# mod_scores <-  GetModuleScores(seurat_ref)
# mod_scores <- mod_scores[,colnames(mod_scores) != 'grey']

# # add hMEs to Seurat meta-data:
# seurat_ref@meta.data <- cbind(
#   seurat_ref@meta.data,
#   mod_scores
# )

# # plot with Seurat's DotPlot function
# p2 <- DotPlot(
#     seurat_ref,
#     features = mixedsort(colnames(mod_scores)),
#     group.by = 'Names'
# )

# # flip the x/y axes, rotate the axis labels, and change color scheme:
# p2 <- p2 +
#   RotatedAxis() +
#   scale_color_gradient2(high='darkorchid4', mid='grey95', low='seagreen') +
#   xlab('') + ylab('')


# pdf(paste0('./output/Figure_8/', "LV_Myeloid_modules_dot.pdf"), width=7.5, height=5)
# p2
# dev.off()

# pdf(paste0('./output/Figure_8/', "LV_Myeloid_LV2RM_modules_dot.pdf"), width=7.5, height=8)
# wrap_plots(list(p1,p2),ncol=1)
# dev.off()



# #######################################
# #############  FIGURE 7B  #############
# #######################################


# RefMerge<-readRDS('./output/Figure_8/Kory_with_RV_modules_projected.rds')

# RefMerge<-SetActiveWGCNA(RefMerge, "CM_RV2LV")


# plot_list <- PlotModulePreservation(
#   RefMerge,
#   name="CM_pres",
#   statistics = "summary"
# )

# pdf(paste0('./output/Figure_8/', 'RV2LV_CM_pres.pdf'), width=5, height=5)
# wrap_plots(plot_list, ncol=2)
# dev.off()

###############################################################################
## FIGURE 8 PANELS G-J — LV (Koenig 2022) re-analysis
###############################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(ggrepel)
  library(ggrastr)
  library(gtools)
})

lv_obj_path <- './dependencies/Figure_8/Kory_with_RV_modules_projected_new.rds'

lv_modules_to_show <- paste0('module_', c('M2','M10','M12','M25','M26','M28'))
lv_cm_mito_modules <- c('M10','M25','M26','M28')

# Panel G
if (file.exists(lv_obj_path)) {
  message('Loading Koenig LV object with projected RV modules...')
  RefLV <- readRDS(lv_obj_path)

  if (!'group' %in% colnames(RefLV@meta.data) && 'condition' %in% colnames(RefLV@meta.data)) {
    RefLV$group <- RefLV$condition
  }

  ## Projected RV modules live in the hdWGCNA slot (GetMEs returns
  ## 'SM-M#' columns), NOT as module_M* meta columns. Extract and inject
  ## them as module_M<n> so the downstream G-J code (which references
  ## module_M2 / paste0('module_', ...) etc.) works unchanged.
  .lv_me <- tryCatch(hdWGCNA::GetMEs(RefLV), error = function(e) NULL)
  if (!is.null(.lv_me)) {
    colnames(.lv_me) <- paste0('module_', sub('^[^-]*-', '', colnames(.lv_me)))
    .lv_me <- .lv_me[rownames(RefLV@meta.data), , drop = FALSE]
    RefLV@meta.data[, colnames(.lv_me)] <- .lv_me
    message(sprintf('Injected %d projected RV module scores into RefLV',
                    ncol(.lv_me)))
  } else {
    message('GetMEs(RefLV) failed — G/H/I will fall back to placeholders.')
  }

  mod_score_cols <- grep('^module_M', colnames(RefLV@meta.data), value = TRUE)
  if (length(mod_score_cols) == 0) {
    mod_score_cols <- lv_modules_to_show
  }
  mod_score_cols <- intersect(mod_score_cols, colnames(RefLV@meta.data))
  mod_score_cols <- gtools::mixedsort(mod_score_cols)

  .g8_mods <- intersect(
    paste0('module_', c('M20','M5','M1','M3','M4','M8','M2','M12',
                        'M25','M26','M10','M28','M14','M11')),
    colnames(RefLV@meta.data))
  p_8G_top <- DotPlot(RefLV, features = .g8_mods, group.by = 'Names',
                       dot.min = 0, col.min = 0, col.max = 2) +
    RotatedAxis() + ylab('') + xlab('') +
    scale_color_gradient2(high = 'red', mid = 'grey95', low = 'blue') +
    theme_v52(COMP_W) +
    theme(panel.border = element_rect(size = 1, fill = NA, color = 'black'),
          axis.line.x = element_blank(),
          axis.line.y = element_blank())

  p_8G_bot <- DotPlot(RefLV, features = .g8_mods, group.by = 'group',
                       dot.min = 0, col.min = 0, col.max = 2) +
    RotatedAxis() + ylab('') + xlab('') +
    scale_color_gradient2(high = 'red', mid = 'grey95', low = 'blue') +
    theme_v52(COMP_W) +
    theme(panel.border = element_rect(size = 1, fill = NA, color = 'black'),
          axis.line.x = element_blank(),
          axis.line.y = element_blank())

  p_8G <- p_8G_top / p_8G_bot

  pdf('./output/Figure_8/LV_modules_celltype_and_disease.pdf', width = 7.5, height = 6)
  print(p_8G)
  dev.off()
  save_figure(p_8G, 'Figure_8_panel_G.pdf', width = 7.5, height = 6)

  # Panel H
  lv_cm <- subset(RefLV, Names == 'Cardiomyocytes' | Names == 'CM')
  pooled_cols <- intersect(paste0('module_', lv_cm_mito_modules), colnames(lv_cm@meta.data))
  if (length(pooled_cols) > 0) {
    lv_cm$pooled_mito <- rowMeans(as.matrix(lv_cm@meta.data[, pooled_cols, drop = FALSE]))
  } else {
    lv_cm$pooled_mito <- NA_real_
  }
  m2_col <- 'module_M2'
  if (!m2_col %in% colnames(lv_cm@meta.data)) {
    lv_cm$module_M2 <- NA_real_
  }

  vio_df <- data.frame(
    group       = as.character(lv_cm$group),
    M2          = as.numeric(lv_cm@meta.data[[m2_col]]),
    pooled_mito = as.numeric(lv_cm$pooled_mito)
  )
  vio_df <- vio_df[vio_df$group %in% c('Donor', 'NF', 'DCM'), , drop = FALSE]
  vio_df$group <- factor(ifelse(vio_df$group == 'Donor', 'NF', vio_df$group), levels = c('NF', 'DCM'))

  p_8H_M2 <- ggplot(vio_df, aes(x = group, y = M2, fill = group)) +
    geom_violin(scale = 'width', trim = TRUE, linewidth = PS$geom_lw) +
    scale_fill_manual(values = c(NF = 'grey70', DCM = '#984EA3'), guide = 'none') +
    ylab('M2 score') + xlab('') +
    theme_v52(COMP_W)

  p_8H_mito <- ggplot(vio_df, aes(x = group, y = pooled_mito, fill = group)) +
    geom_violin(scale = 'width', trim = TRUE, linewidth = PS$geom_lw) +
    scale_fill_manual(values = c(NF = 'grey70', DCM = '#984EA3'), guide = 'none') +
    ylab('Pooled M10/M25/M26/M28') + xlab('') +
    theme_v52(COMP_W)

  p_8H <- p_8H_M2 | p_8H_mito

  pdf('./output/Figure_8/LV_CM_M2_mito_violins.pdf', width = 5, height = 3)
  print(p_8H)
  dev.off()
  save_figure(p_8H, 'Figure_8_panel_H.pdf', width = 5, height = 3)

  # Panel I
  .i8_mods <- intersect(
    paste0('module_', c('M2','M12','M25','M26','M10','M28')),
    colnames(lv_cm@meta.data))
  p_8I <- DotPlot(lv_cm, features = .i8_mods, group.by = 'group',
                    dot.min = 0, col.min = 0, col.max = 2) +
    RotatedAxis() + ylab('') + xlab('') +
    scale_color_gradient2(high = 'red', mid = 'grey95', low = 'blue') +
    theme_v52(COMP_W) +
    theme(panel.border = element_rect(size = 1, fill = NA, color = 'black'),
          axis.line.x = element_blank(),
          axis.line.y = element_blank())

  pdf('./output/Figure_8/LV_CM_modules_dot_disease.pdf', width = 6, height = 3)
  print(p_8I)
  dev.off()
  save_figure(p_8I, 'Figure_8_panel_I.pdf', width = 6, height = 3)

  ## Panel-J input (LV side): CM DCM-vs-NF DEGs from the Koenig LV object,
  ## computed while lv_cm is still alive. Cached under dependencies/ (deps
  ## rule) so later runs skip the FindMarkers.
  .lv_deg_dep <- './dependencies/Figure_8/CM_LV_DCM_vs_NF_DEG.csv'
  if (!file.exists(.lv_deg_dep)) {
    ## The Kory LV RNA 'data' layer is raw counts (not log-normalised),
    ## so FindMarkers' expm1-based avg_log2FC explodes to +/-200. Log-
    ## normalise first so log2FC is bounded (fixes Panel J's blown y-scale).
    DefaultAssay(lv_cm) <- 'RNA'
    lv_cm <- NormalizeData(lv_cm, verbose = FALSE)
    .lv_grp <- unique(as.character(lv_cm$group))
    .lv_ref <- if ('Donor' %in% .lv_grp) 'Donor' else
               if ('NF' %in% .lv_grp) 'NF' else NA_character_
    .lv_d <- tryCatch(
      if ('DCM' %in% .lv_grp && !is.na(.lv_ref))
        FindMarkers(lv_cm, ident.1 = 'DCM', ident.2 = .lv_ref,
                    group.by = 'group') else NULL,
      error = function(e) { message('LV DEG FindMarkers failed: ',
                                    conditionMessage(e)); NULL })
    if (!is.null(.lv_d)) {
      dir.create('./dependencies/Figure_8', showWarnings = FALSE, recursive = TRUE)
      write.csv(.lv_d, .lv_deg_dep)
      message('Wrote LV CM DCM-vs-NF DEGs: ', .lv_deg_dep,
              ' (', nrow(.lv_d), ' genes)')
    }
  }

  rm(RefLV, lv_cm); gc()
} else {
  message('Koenig LV RDS not found: ', lv_obj_path, ' -- emitting placeholder panels G/H/I.')
  .ph8gj <- function(label)
    ggplot() +
      annotate('text', x = 0.5, y = 0.5, label = label,
               size = PS$text_mm, family = FONT_FAMILY, colour = 'grey50', hjust = 0.5, vjust = 0.5) +
      theme_void() +
      theme(panel.border = element_rect(colour = 'grey80', fill = NA, linewidth = PS$geom_lw))

  p_8G <- .ph8gj('[Panel G: LV module dot plots\n(Koenig_LV_with_RV_modules.rds pending)]')
  p_8H <- .ph8gj('[Panel H: LV CM M2 + pooled mito violins\n(Koenig_LV_with_RV_modules.rds pending)]')
  p_8I <- .ph8gj('[Panel I: LV CM module dot by disease\n(Koenig_LV_with_RV_modules.rds pending)]')

  save_figure(p_8G, 'Figure_8_panel_G.pdf', width = 7.5, height = 6)
  save_figure(p_8H, 'Figure_8_panel_H.pdf', width = 5, height = 3)
  save_figure(p_8I, 'Figure_8_panel_I.pdf', width = 6, height = 3)
}

# Panel J — LV vs RV CM log2FC scatter.
## RV side: CM RVF-vs-NF DEGs from the co-embed M2 (still in memory),
## computed + cached under dependencies/. LV side written above.
.rv_deg_dep <- './dependencies/Figure_8/CM_RV_RVF_vs_NF_DEG.csv'
if (!file.exists(.rv_deg_dep) && exists('M2')) {
  .rv_cm <- tryCatch(subset(M2, CombinedNames == 'CM' & origin == 'RV'),
                     error = function(e) NULL)
  if (!is.null(.rv_cm)) {
    .rv_d <- tryCatch(
      FindMarkers(.rv_cm, ident.1 = 'RVF', ident.2 = 'NF',
                  group.by = 'group', recorrect_umi = FALSE),
      error = function(e) { message('RV DEG FindMarkers failed: ',
                                    conditionMessage(e)); NULL })
    if (!is.null(.rv_d)) {
      dir.create('./dependencies/Figure_8', showWarnings = FALSE, recursive = TRUE)
      write.csv(.rv_d, .rv_deg_dep)
      message('Wrote RV CM RVF-vs-NF DEGs: ', .rv_deg_dep,
              ' (', nrow(.rv_d), ' genes)')
    }
    rm(.rv_cm); gc()
  }
}

## Prefer the dependencies copies (deps rule); fall back to ./output.
lv_deg_csv <- if (file.exists('./dependencies/Figure_8/CM_LV_DCM_vs_NF_DEG.csv'))
  './dependencies/Figure_8/CM_LV_DCM_vs_NF_DEG.csv' else
  './output/Figure_8/CM_LV_DCM_vs_NF_DEG.csv'
rv_deg_csv <- if (file.exists('./dependencies/Figure_8/CM_RV_RVF_vs_NF_DEG.csv'))
  './dependencies/Figure_8/CM_RV_RVF_vs_NF_DEG.csv' else
  './output/Figure_8/CM_RV_RVF_vs_NF_DEG.csv'

if (file.exists(lv_deg_csv) && file.exists(rv_deg_csv)) {
  lv_deg <- read.csv(lv_deg_csv, row.names = 1)
  rv_deg <- read.csv(rv_deg_csv, row.names = 1)

  lv_lfc_col <- intersect(c('avg_log2FC', 'log2FoldChange', 'logFC'), colnames(lv_deg))[1]
  rv_lfc_col <- intersect(c('avg_log2FC', 'log2FoldChange', 'logFC'), colnames(rv_deg))[1]

  shared <- intersect(rownames(lv_deg), rownames(rv_deg))
  scatter_df <- data.frame(
    gene = shared,
    LV   = lv_deg[shared, lv_lfc_col],
    RV   = rv_deg[shared, rv_lfc_col]
  )
  scatter_df <- scatter_df[is.finite(scatter_df$LV) & is.finite(scatter_df$RV), , drop = FALSE]

  fit <- lm(LV ~ RV, data = scatter_df)
  r2 <- summary(fit)$r.squared

  scatter_df$resid_abs <- abs(resid(fit))
  outliers <- scatter_df[order(-scatter_df$resid_abs), ][seq_len(min(25, nrow(scatter_df))), , drop = FALSE]

  p_8J <- ggplot(scatter_df, aes(x = RV, y = LV)) +
    ggrastr::rasterise(geom_point(size = PS$scatter_pt, alpha = 0.4, color = 'grey50'), dpi = 200) +
    geom_smooth(method = 'lm', se = FALSE, color = 'firebrick', linewidth = PS$geom_lw) +
    geom_text_repel(data = outliers, aes(label = gene), size = PS$text_mm,
                     family = FONT_FAMILY, fontface = 'italic', max.overlaps = 30) +
    annotate('text', x = -Inf, y = Inf,
             label = paste0('r^2 = ', round(r2, 3)),
             hjust = -0.1, vjust = 1.5,
             size = PS$text_mm, family = FONT_FAMILY) +
    labs(x = 'RV: RVF vs NF log2FC', y = 'LV: DCM vs NF log2FC') +
    theme_v52(COMP_W)

  pdf('./output/Figure_8/LV_vs_RV_log2FC_scatter.pdf', width = 7, height = 7)
  print(p_8J)
  dev.off()
  save_figure(p_8J, 'Figure_8_panel_J.pdf', width = 7, height = 7)
} else {
  message('LV/RV DEG CSV(s) missing -- emitting placeholder Panel J.')
  p_8J <- ggplot() +
    annotate('text', x = 0.5, y = 0.5,
             label = '[Panel J: LV vs RV CM log2FC scatter\n(CM_LV_DCM_vs_NF + CM_RV_RVF_vs_NF DEG CSVs pending)]',
             size = PS$text_mm, family = FONT_FAMILY, colour = 'grey50') +
    theme_void() +
    theme(panel.border = element_rect(colour = 'grey80', fill = NA, linewidth = PS$geom_lw))
  save_figure(p_8J, 'Figure_8_panel_J.pdf', width = 7, height = 7)
}

#Library EnrichR

#CM not preserved - M1, M5
#CM preserved - M7, M9, M10, M11, M14, M14

#EC not preserved - M1, M2, M4, M5, M10, M11, M12
#EC preserved - M8

# seurat_ref <- readRDS('./output/Figure_8/scWGCNA_all_celltypes.rds')
# seurat_ref<-SetActiveWGCNA(seurat_ref, "CM")


# dbs <- c('GO_Biological_Process_2023','GO_Cellular_Component_2023','GO_Molecular_Function_2023','ChEA_2022','Reactome_2022')

# # perform enrichment tests
# seurat_ref <- RunEnrichr(
#   seurat_ref,
#   dbs=dbs, # character vector of enrichr databases to test
#   max_genes = 2000 # number of genes per module to test. use max_genes = Inf to choose all genes!
# )

# # retrieve the output table
# enrich_df <- GetEnrichrTable(seurat_ref)

# EnrichrBarPlot(
#   seurat_ref,
#   outdir = "./output/Figure_8/sc_enrichr_plot_CM", 
#   n_terms = 5,
#   plot_size = c(5,4), # width, height of the output .pdfs
#   logscale=TRUE # do you want to show the enrichment as a log scale?
# )


# seurat_ref<-SetActiveWGCNA(seurat_ref, "EC")


# dbs <- c('GO_Biological_Process_2023','GO_Cellular_Component_2023','GO_Molecular_Function_2023','ChEA_2022','Reactome_2022')

# # perform enrichment tests
# seurat_ref <- RunEnrichr(
#   seurat_ref,
#   dbs=dbs, # character vector of enrichr databases to test
#   max_genes = 2000 # number of genes per module to test. use max_genes = Inf to choose all genes!
# )

# # retrieve the output table
# enrich_df <- GetEnrichrTable(seurat_ref)

# EnrichrBarPlot(
#   seurat_ref,
#   outdir = "./output/Figure_8/sc_enrichr_plot_EC", 
#   n_terms = 5,
#   plot_size = c(5,4), # width, height of the output .pdfs
#   logscale=TRUE # do you want to show the enrichment as a log scale?
# )


# #######################################
# #############  FIGURE 7D  #############
# #######################################
# library(GeneOverlap)

# seurat_ref <- readRDS('./output/Figure_8/scWGCNA_all_celltypes.rds')
# seurat_ref<-SetActiveWGCNA(seurat_ref, "CM")


# bulk_rv_vs_lv <- read.csv(paste0('./output/Figure_8/', 'RV_LV_align_NRVM_NRMV_ARVM.csv'))
# genes <- toupper(bulk_rv_vs_lv[,2])
# #nrvm <- bulk_rv_vs_lv[2:348,14]
# arvm <- bulk_rv_vs_lv[,20]
# arvm.p <-  bulk_rv_vs_lv[,22]
# arvm_rv <-genes[arvm.p < 0.05 & arvm < 0 & genes %in% rownames(seurat_ref)]
# arvm_lv <- genes[arvm.p < 0.05 & arvm > 0 & genes %in% rownames(seurat_ref)]
# arvm_share <-genes[arvm.p >= 1 & arvm < 0.5 & arvm > -0.5 & genes %in% rownames(seurat_ref)]


# # load modules
# modules <- GetModules(seurat_ref)
# mods <- levels(modules$module)
# genome.size <- nrow(modules)


# overlap_df_rv <- do.call(rbind, lapply(mods, function(cur_mod){

#   cur_genes <- modules %>% subset(module == cur_mod) %>% .$gene_name

#   cur_overlap <- testGeneOverlap(newGeneOverlap(
#       cur_genes,
#       arvm_rv,
#       genome.size=genome.size
#   ))

#   cur_overlap <- data.frame(
#     'odds.ratio' = cur_overlap@odds.ratio,
#     'pval' = cur_overlap@pval,
#     'Jaccard' = cur_overlap@Jaccard,
#     'size_intersection' = length(cur_overlap@intersection),
#     'module' = cur_mod
#   )

#   cur_overlap

# })) %>% as.data.frame()

# overlap_df_lv <- do.call(rbind, lapply(mods, function(cur_mod){

#   cur_genes <- modules %>% subset(module == cur_mod) %>% .$gene_name

#   cur_overlap <- testGeneOverlap(newGeneOverlap(
#       cur_genes,
#       arvm_lv,
#       genome.size=genome.size
#   ))

#   cur_overlap <- data.frame(
#     'odds.ratio' = cur_overlap@odds.ratio,
#     'pval' = cur_overlap@pval,
#     'Jaccard' = cur_overlap@Jaccard,
#     'size_intersection' = length(cur_overlap@intersection),
#     'module' = cur_mod
#   )

#   cur_overlap

# })) %>% as.data.frame()

# overlap_df_share <- do.call(rbind, lapply(mods, function(cur_mod){

#   cur_genes <- modules %>% subset(module == cur_mod) %>% .$gene_name

#   cur_overlap <- testGeneOverlap(newGeneOverlap(
#       cur_genes,
#       arvm_share,
#       genome.size=genome.size
#   ))

#   cur_overlap <- data.frame(
#     'odds.ratio' = cur_overlap@odds.ratio,
#     'pval' = cur_overlap@pval,
#     'Jaccard' = cur_overlap@Jaccard,
#     'size_intersection' = length(cur_overlap@intersection),
#     'module' = cur_mod
#   )

#   cur_overlap

# })) %>% as.data.frame()

# overlap_df_rv <- overlap_df_rv %>% mutate(fdr=p.adjust(pval, method='fdr'))
# overlap_df_rv <- overlap_df_rv %>% subset(module != 'grey')



# overlap_df_lv <- overlap_df_lv %>% mutate(fdr=p.adjust(pval, method='fdr'))
# overlap_df_lv <- overlap_df_lv %>% subset(module != 'grey')

# overlap_df_share <- overlap_df_share %>% mutate(fdr=p.adjust(pval, method='fdr'))
# overlap_df_share <- overlap_df_share %>% subset(module != 'grey')

# # Plot as a lollipop

# overlap_df_rv$shape <- ifelse(overlap_df_rv$fdr < 0.05, 21, 4)
# overlap_df_rv <- overlap_df_rv %>% arrange(odds.ratio, descending=TRUE)
# overlap_df_rv$module <- factor(as.character(overlap_df_rv$module), levels=as.character(overlap_df_rv$module))

# mod_colors <- dplyr::select(modules, c(module, color)) %>%
#   distinct
# cp <- mod_colors$color; names(cp) <- mod_colors$module

# p <- overlap_df_rv %>%
#   ggplot(aes(y=module, x=odds.ratio, size= size_intersection, color=module)) +
#   geom_segment(aes(y=module, yend=module, x=0, xend=odds.ratio), size=0.5, color='grey') +
#   geom_point() +
#   geom_point(shape=overlap_df_rv$shape, color='black', fill=NA) +
#   scale_color_manual(values=cp, guide='none') +
#   ylab('') + xlab("Odds ratio") +
#   scale_x_continuous(breaks = c(0, 1, 2,3)) +
#   labs(size='Size\nintersection') +
#   ggtitle('Overlap with SFARI genes') +

#   theme(
#     panel.border = element_rect(size=1, color='black', fill=NA),
#     axis.line.y = element_blank(),
#     axis.line.x = element_blank(),
#     plot.title = element_text(hjust=0.5, face='plain')
#   )

# pdf(paste0('./output/Figure_8/', 'RVbulk_R2LV_CM_overlap.pdf'), width=4.5, height=5)
# p
# dev.off()


# #######################################
# #############  FIGURE 7D  #############
# #######################################
# rm(RefMerge)
# rm(seurat_ref)
# seurat_ref <- readRDS('./dependencies/shared/RV_data.rds')
# load('./output/Figure_8/GSE183852_DCM_Integrated.Robj')


# sc_modules <- read.csv(paste0('./output/Figure_8/', 'sc_heart_modules.csv'))
# seurat_ref <- AddModuleScore(seurat_ref,list(
#   subset(sc_modules,module=="CM-M1")$gene_name,
#   subset(sc_modules,module=="CM-M2")$gene_name,
#   subset(sc_modules,module=="CM-M3")$gene_name,
#   subset(sc_modules,module=="CM-M4")$gene_name,
#   subset(sc_modules,module=="CM-M5")$gene_name,
#   subset(sc_modules,module=="CM-M6")$gene_name,
#   subset(sc_modules,module=="CM-M7")$gene_name
#   ),name='CM_score')
# RefMerge <- AddModuleScore(RefMerge,list(
#   subset(sc_modules,module=="CM-M1")$gene_name,
#   subset(sc_modules,module=="CM-M2")$gene_name,
#   subset(sc_modules,module=="CM-M3")$gene_name,
#   subset(sc_modules,module=="CM-M4")$gene_name,
#   subset(sc_modules,module=="CM-M5")$gene_name,
#   subset(sc_modules,module=="CM-M6")$gene_name,
#   subset(sc_modules,module=="CM-M7")$gene_name
#   )
#   ,name='CM_score')


# FeaturePlot(RefMerge,'CM_score1')
# DotPlot(RefMerge,c('CM_score1','CM_score2','CM_score3','CM_score4','CM_score5','CM_score6','CM_score7'),group.by='Names',split.by='condition')
# DotPlot(seurat_ref,c('CM_score1','CM_score2','CM_score3','CM_score4','CM_score5','CM_score6','CM_score7'),group.by='Names',split.by='group')


# cellsofint<- colnames(RefMerge)[RefMerge@meta.data$'CM-M6'>5 & RefMerge@meta.data$'EC-M7'>5 & RefMerge@meta.data$'CM-M11'>5]
# DimPlot(RefMerge, cells.highlight=cellsofint)
# FindMarkers(RefMerge,ident.1=cellsofint)



# #cellsofint<- colnames(seurat_ref)[seurat_ref@meta.data$'CM-M6'>0 & seurat_ref@meta.data$'EC-M7'>0 & seurat_ref@meta.data$'CM-M11'>0]
# cellsofint<- colnames(seurat_ref)[seurat_ref@meta.data$'CM-M6'>10]

# cellstocomp<- colnames(seurat_ref)[seurat_ref@meta.data$'CM-M6'<10 & seurat_ref@meta.data$Names == "CM"]

# DimPlot(seurat_ref, cells.highlight=cellsofint,raster=T)
# FindMarkers(seurat_ref,ident.1=cellsofint,ident.2=cellstocomp)

###############################################################################
## v52 NEW PANELS — L and M
## Cross-chamber pseudobulk DESeq2: DCM vs RVF
###############################################################################
suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(ggrepel)
  library(ggrastr)
  library(patchwork)
})

## DCM vs RVF requires a merged RV (RVF patients) + LV DCM Seurat object.
## Expected path: ./dependencies/shared/rv_dcm_merged.rds
## Metadata columns expected: ventricle ('RV'/'LV'), condition ('RVF'/'DCM'),
##   Names (lineage), patient
dcm_path <- './dependencies/shared/rv_dcm_merged.rds'

if (file.exists(dcm_path)) {
  dcm_obj <- readRDS(dcm_path)

  ## ── Panel L — PCA of pseudobulk samples (RVF vs DCM, CM lineage) ─────
  cm_dcm  <- subset(dcm_obj, Names == 'CM')
  .cache_dds_dcm <- './output/Figure_8/fig8_dds_dcm_cache.rds'
  if (file.exists(.cache_dds_dcm)) {
    message('Loading cached DESeq2 dds_dcm (DCM vs RVF pseudobulk)...')
    dds_dcm <- readRDS(.cache_dds_dcm)
    agg_dcm <- counts(dds_dcm)
    meta_dcm <- as.data.frame(colData(dds_dcm))
  } else {
    agg_dcm <- AggregateExpression(cm_dcm,
                                    assays   = 'RNA',
                                    group.by = c('patient', 'condition'),
                                    return.seurat = FALSE)[['RNA']]
    meta_dcm <- data.frame(sample = colnames(agg_dcm))
    meta_dcm$condition <- sub('^.*_', '', meta_dcm$sample)
    meta_dcm$patient   <- sub('_.*$', '', meta_dcm$sample)
    rownames(meta_dcm) <- meta_dcm$sample

    dds_dcm <- DESeqDataSetFromMatrix(
      countData = round(agg_dcm),
      colData   = meta_dcm,
      design    = ~ condition)
    dds_dcm <- DESeq(dds_dcm, quiet = TRUE)
    saveRDS(dds_dcm, .cache_dds_dcm)
  }

  ## PCA on VST counts
  vsd_dcm <- vst(dds_dcm, blind = TRUE)
  pca_dcm <- prcomp(t(assay(vsd_dcm)), scale. = FALSE)
  pca_df  <- as.data.frame(pca_dcm$x[, 1:2])
  pca_df$sample    <- rownames(pca_df)
  pca_df$condition <- meta_dcm$condition[match(pca_df$sample, meta_dcm$sample)]

  pct_var <- round(100 * pca_dcm$sdev^2 / sum(pca_dcm$sdev^2), 1)

  p_8L <- ggplot(pca_df, aes(x = PC1, y = PC2,
                              colour = condition, label = sample)) +
    geom_point(size = PS$scatter_pt * 3) +
    ggrepel::geom_label_repel(size = PS$text_mm, family = FONT_FAMILY, max.overlaps = 20, show.legend = FALSE) +
    scale_colour_manual(values = c(RVF = disease_pal['RVF'], DCM = '#984EA3'),
                        name = NULL) +
    labs(x = paste0('PC1 (', pct_var[1], '%)'),
         y = paste0('PC2 (', pct_var[2], '%)'),
         title = 'PCA: CM pseudobulk — RVF vs DCM') +
    theme_v52(COMP_W) +
    theme(legend.key.size = PS$legend_key)

  ## ── Panel M — Volcano: DCM vs RVF CM (shared vs divergent programs) ──
  .cache_res_dcm <- './output/Figure_8/fig8_res_dcm_cache.rds'
  if (file.exists(.cache_res_dcm)) {
    message('Loading cached DESeq2 results res_dcm (DCM vs RVF)...')
    res_dcm <- readRDS(.cache_res_dcm)
  } else {
    res_dcm  <- as.data.frame(results(dds_dcm, contrast = c('condition', 'DCM', 'RVF')))
    res_dcm$gene <- rownames(res_dcm)
    res_dcm      <- res_dcm[!is.na(res_dcm$padj), ]
    saveRDS(res_dcm, .cache_res_dcm)
  }

  ## Label contractile / mito programs
  contractile_genes <- c('MYH7','MYH6','TNNI3','TNNT2','MYBPC3','TPM1','ACTC1')
  mito_genes        <- c('ESRRA','ESRRG','PPARA','PPARGC1A','NRF1','NFE2L2',
                          'NDUFB5','COX5A','UQCRC1','ATP5F1A')
  res_dcm$category <- 'Other'
  res_dcm$category[res_dcm$gene %in% contractile_genes] <- 'Contractile'
  res_dcm$category[res_dcm$gene %in% mito_genes]        <- 'Mito TF/complex'

  cat_pal_M <- c(Contractile = '#377EB8', `Mito TF/complex` = '#FF7F00', Other = 'grey75')

  p_8M <- ggplot(res_dcm, aes(x = log2FoldChange,
                               y = -log10(padj + 1e-300),
                               colour = category,
                               size   = category != 'Other')) +
    ggrastr::rasterise(
      geom_point(data = subset(res_dcm, category == 'Other'), alpha = 0.35),
      dpi = 200) +
    geom_point(data = subset(res_dcm, category != 'Other'), alpha = 0.9) +
    ggrepel::geom_label_repel(
      data = subset(res_dcm, category != 'Other' & padj < 0.05),
      aes(label = gene), size = PS$text_mm, family = FONT_FAMILY, fontface = "italic", max.overlaps = 20, show.legend = FALSE) +
    scale_colour_manual(values = cat_pal_M, name = NULL) +
    scale_size_manual(values = c(`TRUE` = PS$scatter_pt * 2.2, `FALSE` = PS$scatter_pt * 0.7), guide = 'none') +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = 'dashed',
               linewidth = PS$geom_lw, colour = 'grey60') +
    labs(x = 'log\u2082FC (DCM vs RVF)',
         y = '-log\u2081\u2080(FDR)',
         title = 'CM: DCM vs RVF (contractile & mito programs)') +
    theme_v52(COMP_W) +
    theme(legend.position = 'top', legend.key.size = PS$legend_key)

  rm(dcm_obj, cm_dcm, agg_dcm, dds_dcm, vsd_dcm, res_dcm); gc()

} else {
  ## Placeholders
  .ph8 <- function(label)
    ggplot() +
      annotate('text', x = 0.5, y = 0.5, label = label,
               size = PS$text_mm, family = FONT_FAMILY, colour = 'grey50', hjust = 0.5, vjust = 0.5) +
      theme_void() +
      theme(panel.border = element_rect(colour = 'grey80', fill = NA, linewidth = PS$geom_lw))

  p_8L <- .ph8('[Panel L: DCM vs RVF pseudobulk PCA\n(rv_dcm_merged.rds pending)]')
  p_8M <- .ph8('[Panel M: DCM vs RVF volcano\n(rv_dcm_merged.rds pending)]')
}

save_figure(p_8L, 'Figure_8_panel_L.pdf', width = 5, height = 4.5)
save_figure(p_8M, 'Figure_8_panel_M.pdf', width = 6, height = 5)

## ── Assemble panels L-M ──────────────────────────────────────────────────
fig8_lm <- p_8L | p_8M +
  plot_annotation(tag_levels = list(c('L', 'M')),
    theme = theme(plot.tag = element_text(family = FONT_FAMILY, size = 16 * COMP_W / FINAL_WIDTH_IN, face = "bold")))

save_figure(fig8_lm, 'Figure_8_panels_L-M.pdf', width = COMP_W, height = COMP_H)
message('Figure 8 panels L-M complete.')




