###############################################################################
## Figure 4 (v52 draft) — RV cardiomyocyte transcriptional and metabolic remodeling
##
## Panels (from RV_snRNASeq_v52_draft.md figure legends):
##   (A) DotPlot of CM-enriched WGCNA module scores by CM subcluster/disease
##   (B) GO enrichment for CM-specific modules
##   (C) EnrichR for M2 module genes
##   (D) Volcano of M2 module DEGs RVF vs pRV (PRKAR1A, PDE1C, PDE4B, CFL2)
##   (E) Violin plots of key M12 genes (TNNI3, MYH7, BMPR2, TNNI3K)
##   (F) Module volcanos for mito modules M10/M25/M26/M28 (nuclear-encoded ETC)
##   (G) Xenium endocardial enrichment of NPPA/NPPB CMs
##   (H) Oroboros respirometry — adult RV (NF/pRV/RVF, top) AND pediatric RV
##       (NF vs HLHS-palliated failing, right) on a single composed panel.
##       (Note: assembled in Illustrator from two sub-plots rendered below.)
##   (I) Respirometry PAB 2-month (top) and 2-week (bottom) vs sham
##   (J) Cross-cohort summary: maximal coupled respiration normalized to
##       within-cohort controls across adult RVF, murine PAB, pediatric HLHS.
##   (K) Mito TF heatmap (ESRRA, ESRRG, PPARA, PPARGC1A, NRF1, NFE2L2) at bulk
##   (L) Single-nucleus violin plots of mito TFs in CM subtypes
##
## Source: copied from v51 Figure_4.R on 2026-04-10
## Status: SKELETON — v52 porting pending
##
## Output: ./output/Figure_5/v52_figures/Figure_4.pdf
###############################################################################

source('./helper_scripts/_shared_helpers.R')

## Per-figure output directory (introduced for consistent output paths)
V52_FIG_DIR <- './output/Figure_5'
dir.create(V52_FIG_DIR, showWarnings = FALSE, recursive = TRUE)


## Suppress R's default Rplots.pdf in cwd when Rscript hits a plot call
## that's outside an explicit pdf() ... dev.off() envelope.
pdf(NULL)
dir.create(file.path(V52_FIG_DIR, 'Xenium'), showWarnings = FALSE, recursive = TRUE)
COMP_W <- 14
COMP_H <- 10
PS <- pub_scales(COMP_W)

## -- Cartoon asset extraction for Figure 4 (run-once, commented) ------------
## Panel A of the published Figure_4.pdf has a small cardiomyocyte tissue
## cartoon icon (top-left) next to the "Cardiomyocyte" DotPlot header.
## Cropped from the original 17-cm-wide PDF at 300 DPI (~2008 px wide).
## Estimated crop (tune by eye):
##   widthxheight+x_offset+y_offset = 220x120+0+0
# library(magick)
# fig4 <- image_read_pdf("~/Downloads/hdWGCNA_TOM/Manuscripts/Figure_4.pdf", density=300)
# cm_cartoon <- image_crop(fig4, "220x120+0+0")
# image_write(cm_cartoon, "./new_scripts/assets/Figure_4_panel_A_cardiomyocyte.png")

library(Seurat)
library(hdWGCNA)
library(ggeasy)
library(harmony)
library(viridis)
library(stringr)



source('./helper_scripts/spatial_functions.R')

## Safe enrichR wrapper — returns empty named list if API is unreachable
.safe_enrichr <- function(genes, dbs) {
  tryCatch(
    enrichR::enrichr(genes, dbs),
    error = function(err) {
      msg <- tryCatch(err$message, error = function(e2) 'unknown error')
      message('enrichR unavailable (', msg, ') — skipping enrichment.')
      setNames(lapply(dbs, function(d)
        data.frame(Term='', Overlap='', P.value=NA, Adjusted.P.value=NA,
                   Combined.Score=NA, stringsAsFactors=FALSE)), dbs)
    }
  )
}


#######################################
#############  FIGURE 5A  #############
#######################################

mapping <- labels2colors(1:100)
.cache_MEs_bulk2sn <- './dependencies/shared/scWGCNA_bulk2sn_MEs.rds'

if (file.exists(.cache_MEs_bulk2sn)) {
  message('Loading cached MEs (bulk2sn) - skipping 20 GB scWGCNA load')
  MEs <- readRDS(.cache_MEs_bulk2sn)
} else {
  message('No MEs cache; loading scWGCNA_bulk2sn_projection.rds (one-time)')
  seurat_ref <- readRDS('./dependencies/shared/scWGCNA_bulk2sn_projection.rds')
  seurat_ref <- SetActiveWGCNA(seurat_ref, 'bulk2sn')
  dir.create(dirname(.cache_MEs_bulk2sn), showWarnings = FALSE, recursive = TRUE)
  MEs <- GetMEs(seurat_ref, harmonized = TRUE)
  saveRDS(MEs, .cache_MEs_bulk2sn)
  rm(seurat_ref); gc()
}
mods <- colnames(MEs); mods <- mods[mods != 'grey']
mods_num <- paste0('M',match(mods,mapping))
modules_up <- c('M2','M12','M28')
modules_down <- c('M10','M25','M26')



colnames(MEs) <- paste0('M', match(colnames(MEs), mapping))

consensus_modules <- read.csv("./dependencies/shared/bulk_heart_modules.csv")
consensus_modules <- consensus_modules[,1:3]


seurat_ref <- readRDS('./dependencies/shared/cm_subclust_new_new.rds')

# Gene-presence filter deferred until cm_subclust_new_new is loaded:
# AddModuleScore (called below) only scores genes in rownames(seurat_ref),
# so filtering here matches the gene universe actually used downstream.
consensus_modules <- subset(consensus_modules, gene_name %in% rownames(seurat_ref))
consensus_modules <- consensus_modules[match(unique(consensus_modules$gene_name), consensus_modules$gene_name),]

library(dplyr)
score_calc <- consensus_modules %>% group_by(module) %>% group_split()
module_colors <- unique(unlist(lapply(score_calc,'[[','module')))
module_colors <- paste0('M',match(module_colors,mapping))

# Idents and $Subnames already carry biological labels in this object:
# CM_Baseline / CM_BetaMHC / CM_XIRP / CM_RORA / CM_HAND2 / CM_HTR4 /
# CM_HMGCS2 / CM_KCNJ3 / CM_MYH6 / CM_NPP

idx <- rownames(MEs) %in% colnames(seurat_ref)
seurat_ref$SubNames_Groups <- paste(seurat_ref$Subnames,seurat_ref$group,sep='_')

seurat_ref <- SetIdent(seurat_ref, value = "Subnames")



seurat_ref@meta.data <- cbind(seurat_ref@meta.data, MEs[idx,])



.cache_module_score <- './output/Figure_5/fig5_module_score_cache.rds'
if (file.exists(.cache_module_score)) {
  message('Loading cached AddModuleScore columns...')
  .module_score_cols <- readRDS(.cache_module_score)
  seurat_ref@meta.data <- cbind(seurat_ref@meta.data, .module_score_cols)
} else {
  seurat_ref <- AddModuleScore(seurat_ref,lapply(score_calc,'[[','gene_name'),name="module_score")
  cols_current <- colnames(seurat_ref@meta.data)
  cols_current[startsWith(colnames(seurat_ref@meta.data),'module_score')] <- paste0('module_',module_colors)
  colnames(seurat_ref@meta.data) <- cols_current
  .module_score_cols <- seurat_ref@meta.data[, paste0('module_',module_colors), drop=FALSE]
  saveRDS(.module_score_cols, .cache_module_score)
}

#Dot Plot of enrichment cell type of CM enriched modules

## Shared dot-plot styling — matches reference Figure 5E (single-hue red,
## panel border, light grey gridlines for inter-dot spacing visibility).
.cm_dot_theme <- theme(
  panel.background = element_rect(fill = 'white', color = NA),
  panel.border     = element_rect(size = 1, fill = NA, color = 'black'),
  panel.grid.major = element_line(size = 0.25, color = 'lightgrey'),
  panel.grid.minor = element_blank(),
  axis.line.x      = element_blank(),
  axis.line.y      = element_blank(),
  legend.key.size  = PS$legend_key
)
.cm_dot_color <- scale_color_gradient(low = 'grey95', high = 'red',
                                      limits = c(0, 2), name = 'Avg Exp')

pdf(paste0('./output/Figure_5/', 'Figure_5_panel_A_dot_subclust_up.pdf'), width=5, height=5)

print({
p <- DotPlot(seurat_ref,paste0('module_',modules_up),group.by='Subnames',dot.min=0,col.min=0,col.max=2) +
  RotatedAxis() + ylab('')+ xlab('')+
  .cm_dot_color +
  scale_size(range = PS$dot_range, name = '% Exp') +
  .cm_dot_theme
p
})
dev.off()

pdf(paste0('./output/Figure_5/', 'Figure_5_panel_A_dot_subclust_down.pdf'), width=5, height=5)

print({
p <- DotPlot(seurat_ref,paste0('module_',modules_down),group.by='Subnames',dot.min=0,col.min=0,col.max=2) +
  RotatedAxis() + ylab('')+ xlab('')+
  .cm_dot_color +
  scale_size(range = PS$dot_range, name = '% Exp') +
  .cm_dot_theme
p
})
dev.off()


seurat_ref <- SetIdent(seurat_ref, value = "group")
my_levels <- c("NF","pRV","RVF")
Idents(seurat_ref) <- factor(Idents(seurat_ref), levels= my_levels)

pdf(paste0('./output/Figure_5/', 'Figure_5_supp_sc_trend_CM.pdf'), width=4.5, height=2)

print({
p <- DotPlot(seurat_ref,paste0('module_',
  c('M2','M12','M25','M26','M10','M28')),dot.min=0,col.min=0,col.max=2) +
  RotatedAxis() + ylab('')+ xlab('')+
  .cm_dot_color +
  scale_size(range = PS$dot_range, name = '% Exp') +
  .cm_dot_theme
p
})
dev.off()

seurat_ref <- SetIdent(seurat_ref, value = "Subnames")


pdf(paste0('./output/Figure_5/', 'Figure_5_supp_sc_trend_CM_subnames.pdf'), width=4.5, height=5)

print({
p <- DotPlot(seurat_ref,paste0('module_',
  c('M2','M12','M25','M26','M10','M28')),dot.min=0,col.min=0,col.max=2) +
  RotatedAxis() + ylab('')+ xlab('')+
  .cm_dot_color +
  scale_size(range = PS$dot_range, name = '% Exp') +
  .cm_dot_theme
p
})
dev.off()


#######################################
#############  FIGURE 5B  #############
#######################################
Idents(seurat_ref) <- "group"

dbs <-c('GO_Biological_Process_2023','GO_Cellular_Component_2023','GO_Molecular_Function_2023','Reactome_2022', 'ChEA_2022',"LINCS_L1000_Chem_Pert_Consensus_Sigs")


#Run enrichment by cell type
.cache_4B <- './output/Figure_5/fig5_4B_findmarkers_enrichr_cache.rds'
if (file.exists(.cache_4B)) {
  message('Loading cached Fig4B FindMarkers + enrichR results...')
  .cache_4B_data <- readRDS(.cache_4B)
  combined_set    <- .cache_4B_data$combined_set
  combined_output <- .cache_4B_data$combined_output
} else {
  combined_set <- data.frame()
  combined_output <- data.frame(
    Term=character(), Overlap=character(), P.value=numeric(),
    Adjusted.P.value=numeric(), Combined.Score=numeric(),
    db=character(), module=character(), celltype=character(),
    comparison=character(), color=character(), direction=character(),
    stringsAsFactors=FALSE)
  bulk_modules <- consensus_modules
  bulk_modules$module <- match(consensus_modules$module,mapping)
  mods_idx <- list(2,12,10,25,26,28)
  cell_types <- unique(seurat_ref$Names)
  comparison <- list(c("RVF","NF"),c("RVF","pRV"))
  for (i in mods_idx){
    for (j in cell_types){
      for (k in comparison){
        key_genes <- subset(bulk_modules,module %in% c(i))$gene_name
        key_genes <- key_genes[key_genes %in% rownames(seurat_ref)]

        gene_set <- FindMarkers(seurat_ref, ident.1 = paste0(k[1]), ident.2 = paste0(k[2]),features=key_genes)

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
        enriched <- .safe_enrichr(rownames(gene_enrich), dbs)
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
        enriched <- .safe_enrichr(rownames(gene_enrich), dbs)
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
  saveRDS(list(combined_set=combined_set, combined_output=combined_output), .cache_4B)
}
bulk_modules <- consensus_modules
bulk_modules$module <- match(consensus_modules$module,mapping)

wrapText <- function(x, len) {
    sapply(x, function(y) paste(strwrap(y, len), collapse = "\n"), USE.NAMES = FALSE)
}


if (nrow(combined_output) > 0) {  ## Panel 4B plots — skip when enrichR unavailable

###RVF vs NF

#Up
selected_terms <- subset(combined_output,db=="GO_Biological_Process_2023")
selected_terms <- subset(selected_terms,direction=="up")
selected_terms <- subset(selected_terms,comparison=="RVF_NF")
selected_terms <- subset(selected_terms,color %in% mapping[c(2,12)])


# subset selected terms
selected_terms <- subset(selected_terms, P.value < 0.05)
selected_terms$module_celltype <- paste0(selected_terms$module,'_',selected_terms$celltype)
idx_top_1 <- match(unique(selected_terms$module_celltype),selected_terms$module_celltype)
#idx_top_5 <- sort(c(idx_top_1,idx_top_1+1,idx_top_1+2,idx_top_1+3,idx_top_1+4))

selected_terms<-selected_terms[idx_top_1,]


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
selected_terms$Term <- str_replace(selected_terms$Term, "\\ (GO).*", "")

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
  scale_size(range = PS$dot_range) +
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
    panel.grid = element_line(size=0.25, color='lightgrey'),
    legend.key.size = PS$legend_key
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


pdf(paste0('./output/Figure_5/', 'Figure_5_panel_B_GO_modules_up.pdf'), width=6, height=4)
print({
p / colorbar
})
dev.off()



###RVF vs NF

#Down
selected_terms <- subset(combined_output,db=="GO_Biological_Process_2023")
selected_terms <- subset(selected_terms,direction=="down")
selected_terms <- subset(selected_terms,comparison=="RVF_NF")
selected_terms <- subset(selected_terms,color %in% mapping[c(10,25,26,28)])


# subset selected terms
selected_terms <- subset(selected_terms, P.value < 0.05)
selected_terms$module_celltype <- paste0(selected_terms$module,'_',selected_terms$celltype)
idx_top_1 <- match(unique(selected_terms$module_celltype),selected_terms$module_celltype)
#idx_top_5 <- sort(c(idx_top_1,idx_top_1+1,idx_top_1+2,idx_top_1+3,idx_top_1+4))

selected_terms<-selected_terms[idx_top_1,]


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
selected_terms$Term <- str_replace(selected_terms$Term, "\\ (GO).*", "")

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
  scale_size(range = PS$dot_range) +
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
    panel.grid = element_line(size=0.25, color='lightgrey'),
    legend.key.size = PS$legend_key
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


pdf(paste0('./output/Figure_5/', 'Figure_5_panel_B_GO_modules_down.pdf'), width=6, height=4)
print({
p / colorbar
})
dev.off()



#Both
selected_terms <- subset(combined_output,db=="GO_Biological_Process_2023")
#selected_terms <- subset(selected_terms,direction=="up")
selected_terms <- subset(selected_terms,comparison=="RVF_NF")
selected_terms <- subset(selected_terms,color %in% mapping[c(2,12,10,25,26,28)])


# subset selected terms
selected_terms <- subset(selected_terms, P.value < 0.05)
selected_terms$module_celltype <- paste0(selected_terms$module,'_',selected_terms$celltype)
idx_top_1 <- match(unique(selected_terms$module_celltype),selected_terms$module_celltype)
#idx_top_5 <- sort(c(idx_top_1,idx_top_1+1,idx_top_1+2,idx_top_1+3,idx_top_1+4))

selected_terms<-selected_terms[idx_top_1,]


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
selected_terms$Term <- str_replace(selected_terms$Term, "\\ (GO).*", "")

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
  scale_size(range = PS$dot_range) +
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
    panel.grid = element_line(size=0.25, color='lightgrey'),
    legend.key.size = PS$legend_key
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


pdf(paste0('./output/Figure_5/', 'Figure_5_panel_B_GO_modules_both.pdf'), width=6, height=4)
print({
p / colorbar
})
dev.off()

} else { message('Panel 4B: skipping enrichment plots (enrichR unavailable).') }


#######################################
#############  FIGURE 5C  #############
#######################################
#Deep dive M2

bulk_modules <- read.csv("./dependencies/shared/bulk_heart_modules.csv")
bulk_modules$module <- match(bulk_modules$module,mapping)

#RVF vs NF
Idents(seurat_ref) <- "group"

.cache_4C_NF <- './output/Figure_5/fig5_4C_findmarkers_RVF_vs_NF_cache.rds'
if (file.exists(.cache_4C_NF)) {
  message('Loading cached Fig4C FindMarkers (RVF vs NF)...')
  combined_set <- readRDS(.cache_4C_NF)
} else {
  combined_set <- data.frame()
  mods_idx <- c(2,12,28,10,25,26)
  for (i in mods_idx){
    key_genes <- subset(bulk_modules,module %in% c(i))$gene_name
    key_genes <- key_genes[key_genes %in% rownames(seurat_ref)]
    gene_set <- FindMarkers(seurat_ref, ident.1 = "RVF", ident.2 = "NF",features=key_genes)
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
  saveRDS(combined_set, .cache_4C_NF)
}

M2_genes_up <- rownames(subset(combined_set,module=="M2" & avg_log2FC>0))
M12_genes_up <- rownames(subset(combined_set,module=="M12" & avg_log2FC>0))
M2_genes_down <- rownames(subset(combined_set,module=="M2" & avg_log2FC<0))
M12_genes_down <- rownames(subset(combined_set,module=="M12" & avg_log2FC<0))

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

.cache_enrichr_M2_up_NF <- './output/Figure_5/fig5_enrichr_M2_genes_up_RVF_vs_NF_cache.rds'
if (file.exists(.cache_enrichr_M2_up_NF)) {
  message('Loading cached enrichR M2_genes_up (RVF vs NF)...')
  enriched <- readRDS(.cache_enrichr_M2_up_NF)
} else {
  enriched <- .safe_enrichr(M2_genes_up, dbs)
  saveRDS(enriched, .cache_enrichr_M2_up_NF)
}
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('./output/Figure_5/Figure_5_panel_C_M2_up_RVF_vs_NF.pdf',width=5,height=2.5)
print({
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
p1
})
dev.off()


.cache_enrichr_M12_up_NF <- './output/Figure_5/fig5_enrichr_M12_genes_up_RVF_vs_NF_cache.rds'
if (file.exists(.cache_enrichr_M12_up_NF)) {
  message('Loading cached enrichR M12_genes_up (RVF vs NF)...')
  enriched <- readRDS(.cache_enrichr_M12_up_NF)
} else {
  enriched <- .safe_enrichr(M12_genes_up, dbs)
  saveRDS(enriched, .cache_enrichr_M12_up_NF)
}
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('./output/Figure_5/Figure_5_supp_M12_enrichr_up_solo.pdf',width=5,height=2.5)
print({
p2<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:5),], 
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
p2
})
dev.off()

pdf('./output/Figure_5/Figure_5_supp_M12_enrichr_up.pdf',width=5,height=5)
print({
p1/p2
})
dev.off()

.cache_enrichr_M2_down_NF <- './output/Figure_5/fig5_enrichr_M2_genes_down_RVF_vs_NF_cache.rds'
if (file.exists(.cache_enrichr_M2_down_NF)) {
  message('Loading cached enrichR M2_genes_down (RVF vs NF)...')
  enriched <- readRDS(.cache_enrichr_M2_down_NF)
} else {
  enriched <- .safe_enrichr(M2_genes_down, dbs)
  saveRDS(enriched, .cache_enrichr_M2_down_NF)
}
pdf('./output/Figure_5/Figure_5_panel_C_M2_down_RVF_vs_NF.pdf',width=5,height=2.5)
print({
p3<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:5),], 
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
p3
})
dev.off()


.cache_enrichr_M12_down_NF <- './output/Figure_5/fig5_enrichr_M12_genes_down_RVF_vs_NF_cache.rds'
if (file.exists(.cache_enrichr_M12_down_NF)) {
  message('Loading cached enrichR M12_genes_down (RVF vs NF)...')
  enriched <- readRDS(.cache_enrichr_M12_down_NF)
} else {
  enriched <- .safe_enrichr(M12_genes_down, dbs)
  saveRDS(enriched, .cache_enrichr_M12_down_NF)
}
pdf('./output/Figure_5/Figure_5_supp_M12_enrichr_down.pdf',width=5,height=2.5)
print({
p4<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:5),], 
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
p4
})
dev.off()





bulk_modules <- read.csv("./dependencies/shared/bulk_heart_modules.csv")
bulk_modules$module <- match(bulk_modules$module,mapping)

#RVF vs pRV
Idents(seurat_ref) <- "group"

.cache_4C_pRV <- './output/Figure_5/fig5_4C_findmarkers_pRV_vs_NF_cache.rds'
if (file.exists(.cache_4C_pRV)) {
  message('Loading cached Fig4C FindMarkers (pRV vs NF)...')
  combined_set <- readRDS(.cache_4C_pRV)
} else {
  combined_set <- data.frame()
  mods_idx <- c(2,12,28,10,25,26)
  for (i in mods_idx){
    key_genes <- subset(bulk_modules,module %in% c(i))$gene_name
    key_genes <- key_genes[key_genes %in% rownames(seurat_ref)]
    gene_set <- FindMarkers(seurat_ref, ident.1 = "pRV", ident.2 = "NF",features=key_genes)
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
  saveRDS(combined_set, .cache_4C_pRV)
}

M2_genes_up <- rownames(subset(combined_set,module=="M2" & avg_log2FC>0))
M12_genes_up <- rownames(subset(combined_set,module=="M12" & avg_log2FC>0))
M2_genes_down <- rownames(subset(combined_set,module=="M2" & avg_log2FC<0))
M12_genes_down <- rownames(subset(combined_set,module=="M12" & avg_log2FC<0))

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

.cache_enrichr_M2_up_pRV <- './output/Figure_5/fig5_enrichr_M2_genes_up_pRV_vs_NF_cache.rds'
if (file.exists(.cache_enrichr_M2_up_pRV)) {
  message('Loading cached enrichR M2_genes_up (pRV vs NF)...')
  enriched <- readRDS(.cache_enrichr_M2_up_pRV)
} else {
  enriched <- .safe_enrichr(M2_genes_up, dbs)
  saveRDS(enriched, .cache_enrichr_M2_up_pRV)
}
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('./output/Figure_5/Figure_5_panel_C_M2_up_pRV_vs_NF.pdf',width=5,height=2.5)
print({
p5<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:5),], 
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
p5
})
dev.off()


.cache_enrichr_M2_down_pRV <- './output/Figure_5/fig5_enrichr_M2_genes_down_pRV_vs_NF_cache.rds'
if (file.exists(.cache_enrichr_M2_down_pRV)) {
  message('Loading cached enrichR M2_genes_down (pRV vs NF)...')
  enriched <- readRDS(.cache_enrichr_M2_down_pRV)
} else {
  enriched <- .safe_enrichr(M2_genes_down, dbs)
  saveRDS(enriched, .cache_enrichr_M2_down_pRV)
}
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('./output/Figure_5/Figure_5_panel_C_M2_down_pRV_vs_NF.pdf',width=5,height=2.5)
print({
p6<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:5),], 
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
p6
})
dev.off()

pdf('./output/Figure_5/Figure_5_panel_C_M2_up_composite.pdf',width=5,height=4.5)
print({
p1/p5
})
dev.off()



#######################################
###########################################################################
## LEGACY (was Panel 5D in v54) — single-cell-level M2 volcano.
## v57 Panel D is the per-patient pseudobulk DESeq2 + ashr volcano below;
## this single-cell version under-reports p-values by treating each nucleus
## as independent. Outputs go to *_supp_singlecell_volcano_*.
##
## Wrapped in tryCatch because in some Seurat/EnhancedVolcano combinations
## this block surfaces an opaque parse-error during plot construction; the
## failure must not block real Panel D / E / F / G / H / I / J downstream.
###########################################################################
tryCatch({
library(EnhancedVolcano)

temp_set <- combined_set[!grepl('MT-',rownames(combined_set)),]
temp_set <- subset(temp_set,module == "M2")
keyvals <- temp_set$color
names(keyvals)  <- temp_set$module
temp_set$p_val_adj[temp_set$p_val_adj < 1e-50] = 1e-50
temp_set$avg_log2FC[temp_set$avg_log2FC < -5] = -4.99
temp_set$avg_log2FC[temp_set$avg_log2FC > 5] = 5



pdf(paste0('./output/Figure_5/', 'Figure_5_supp_singlecell_volcano_M2_RVF_vs_pRV.pdf'), width=8, height=6)

print({
EnhancedVolcano(temp_set,lab=rownames(temp_set),
  x='avg_log2FC',y='p_val_adj',xlim=c(-5,5),
  FCcutoff = 0.1, labFace = "italic") + coord_flip()
})
dev.off()

pdf(paste0('./output/Figure_5/', 'Figure_5_supp_ACTA1_vln.pdf'), width=3, height=2.8)
print({
VlnPlot(seurat_ref,'ACTA1',group.by='group')
})
dev.off()
}, error = function(e) message('LEGACY Panel 5D single-cell volcano skipped: ', conditionMessage(e)))



#######################################
####  Panel E precomputation: M12 deep-dive (RVF vs NF FindMarkers)
####  Produces FindMarkers cache used downstream by mito-module volcano plot.
#######################################

#Deep dive M12

bulk_modules <- read.csv("./dependencies/shared/bulk_heart_modules.csv")
bulk_modules$module <- match(bulk_modules$module,mapping)

#RVF vs NF
Idents(seurat_ref) <- "group"

if (file.exists(.cache_4C_NF)) {
  message('Loading cached FindMarkers (RVF vs NF) from Fig4C cache...')
  combined_set <- readRDS(.cache_4C_NF)
} else {
  combined_set <- data.frame()
  mods_idx <- c(2,12,28,10,25,26)
  for (i in mods_idx){
    key_genes <- subset(bulk_modules,module %in% c(i))$gene_name
    key_genes <- key_genes[key_genes %in% rownames(seurat_ref)]
    gene_set <- FindMarkers(seurat_ref, ident.1 = "RVF", ident.2 = "NF",features=key_genes)
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
  saveRDS(combined_set, .cache_4C_NF)
}


temp_set <- combined_set[!grepl('MT-',rownames(combined_set)),]
temp_set <- subset(temp_set,module == "M12")
keyvals <- temp_set$color
names(keyvals)  <- temp_set$module
temp_set$p_val_adj[temp_set$p_val_adj < 1e-50] = 1e-50
temp_set$avg_log2FC[temp_set$avg_log2FC < -5] = -4.99
temp_set$avg_log2FC[temp_set$avg_log2FC > 5] = 5

#LRRC39,TNNI3,BMPR2,MYH7,TNNI3K

pdf(paste0('./output/Figure_5/', 'Figure_5_supp_M12_key_vln.pdf'), width=6, height=3)

print({
VlnPlot(seurat_ref,c('TNNI3','MYH7','BMPR2','TNNI3K'),group.by='group',ncol=4)
})
dev.off()


#######################################
####  Panel E precomputation: bulk × single-cell CM module gene overlap
####  Feeds the mito-module volcano below.
#######################################

library(GeneOverlap)

bulk_modules <- read.csv("./dependencies/shared/bulk_heart_modules.csv")
bulk_modules <- subset(bulk_modules,color %in% mods)

sc_modules <- read.csv("./dependencies/shared/sc_heart_modules.csv")
sc_modules <- subset(sc_modules,str_replace(module,"\\-.*", "")
 == "CM")

bulk_modules <- subset(bulk_modules,gene_name %in% sc_modules$gene_name)
sc_modules <- subset(sc_modules,gene_name %in% bulk_modules$gene_name)

idx<-match(bulk_modules$gene_name,sc_modules$gene_name)

bulk2sc_mods<-data.frame(gene_name=bulk_modules$gene_name,bulk_module=bulk_modules$module,sc_module=sc_modules[idx,2])
bulk2sc_mods$bulk_module <- factor(bulk2sc_mods$bulk_module)
bulk2sc_mods$sc_module <- factor(bulk2sc_mods$sc_module)


#table(bulk2sc_mods$sc_module,match(bulk2sc_mods$bulk_module,mapping))
cont_tbl <- table(bulk2sc_mods$sc_module,bulk2sc_mods$bulk_module)
cont_tbl <- t(cont_tbl) / colSums(cont_tbl)

#cat(subset(bulk2sc_mods,bulk_module == 'tan' & sc_module == 'CM-M1')$gene_name,sep='\n')
#cat(subset(bulk2sc_mods,bulk_module == 'tan' & sc_module == 'CM-M5')$gene_name,sep='\n')
#cat(subset(bulk2sc_mods,bulk_module == 'tan' & sc_module == 'CM-M6')$gene_name,sep='\n')

cat(subset(bulk2sc_mods,bulk_module == 'blue' & sc_module == 'CM-M10')$gene_name,sep='\n')
cat(subset(bulk2sc_mods,bulk_module == 'blue' & sc_module == 'CM-M4')$gene_name,sep='\n')
cat(subset(bulk2sc_mods,bulk_module == 'blue' & sc_module == 'CM-M5')$gene_name,sep='\n')
cat(subset(bulk2sc_mods,bulk_module == 'blue' & sc_module == 'CM-M6')$gene_name,sep='\n')

cat(subset(bulk2sc_mods,sc_module == 'CM-M6' & (bulk_module == 'blue' | bulk_module == 'tan' | bulk_module == 'skyblue'))$gene_name,sep='\n')
cat(subset(bulk2sc_mods,sc_module == 'CM-M5' & (bulk_module == 'blue' | bulk_module == 'tan' | bulk_module == 'skyblue'))$gene_name,sep='\n')

cat(subset(bulk2sc_mods,bulk_module == 'blue')$gene_name,sep='\n')

VlnPlot(seurat_ref,'module_M2',pt.size=0,group.by='group')

bulk_modules <- read.csv("./dependencies/shared/bulk_heart_modules.csv")
bulk_modules$module <- match(bulk_modules$module,mapping)

#key_genes <- subset(bulk_modules,module %in% c(10,25,26))$gene_name
#key_genes <- key_genes[key_genes %in% rownames(seurat_ref)]
#Idents(seurat_ref) <- "group"


#gene_set <- FindMarkers(seurat_ref, ident.1 = "NF", ident.2 = "RVF",features=key_genes,only.pos=T)
#gene_set<-subset(gene_set,p_val_adj<0.05)


#key_genes <- subset(bulk_modules,module %in% c(2))$gene_name
#key_genes <- key_genes[key_genes %in% rownames(seurat_ref)]
#Idents(seurat_ref) <- "group"


#gene_set <- FindMarkers(seurat_ref, ident.1 = "RVF", ident.2 = "NF",features=key_genes,only.pos=T)
#ene_set<-subset(gene_set,p_val_adj<0.05)

#key_genes <- subset(bulk_modules,module %in% c(12))$gene_name
#key_genes <- key_genes[key_genes %in% rownames(seurat_ref)]
#Idents(seurat_ref) <- "group"


#gene_set <- FindMarkers(seurat_ref, ident.1 = "RVF", ident.2 = "NF",features=key_genes,only.pos=T)
#gene_set<-subset(gene_set,p_val_adj<0.05)


#RVF vs NF
Idents(seurat_ref) <- "group"

if (file.exists(.cache_4C_NF)) {
  message('Loading cached FindMarkers (RVF vs NF) from Fig4C cache...')
  combined_set <- readRDS(.cache_4C_NF)
} else {
  combined_set <- data.frame()
  mods_idx <- c(2,12,28,10,25,26)
  for (i in mods_idx){
    key_genes <- subset(bulk_modules,module %in% c(i))$gene_name
    key_genes <- key_genes[key_genes %in% rownames(seurat_ref)]
    gene_set <- FindMarkers(seurat_ref, ident.1 = "RVF", ident.2 = "NF",features=key_genes)
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
  saveRDS(combined_set, .cache_4C_NF)
}

library(EnhancedVolcano)

combined_set <- combined_set[!grepl('MT-',rownames(combined_set)),]
keyvals <- combined_set$color
names(keyvals)  <- combined_set$module

pdf(paste0('./output/Figure_5/', 'Figure_5_supp_singlecell_volcano_M2_RVF_vs_NF.pdf'), width=8, height=6)

print({
EnhancedVolcano(combined_set,lab=rownames(combined_set),
  x='avg_log2FC',y='p_val_adj',
  FCcutoff = 0.1,
  colCustom=keyvals, labFace = "italic") + coord_flip()
})
dev.off()


sub_set <- subset(combined_set,combined_set$module != 'M2')
sub_set <- sub_set[!grepl('MT-',rownames(sub_set)),]
keyvals <- sub_set$color
names(keyvals)  <- sub_set$module

pdf(paste0('./output/Figure_5/', 'Figure_5_supp_singlecell_volcano_noM2_RVF_vs_NF.pdf'), width=8, height=6)

print({
EnhancedVolcano(sub_set,lab=rownames(sub_set),
  x='avg_log2FC',y='p_val_adj',
  FCcutoff = 0.1,
  colCustom=keyvals, labFace = "italic") + coord_flip()
})
dev.off()

subsub_set <- subset(sub_set,sub_set$module != 'M12')
keyvals <- subsub_set$color
names(keyvals)  <- subsub_set$module

pdf(paste0('./output/Figure_5/', 'Figure_5_supp_singlecell_volcano_mito_celllevel.pdf'), width=8, height=6)

print({
EnhancedVolcano(subsub_set,lab=rownames(subsub_set),
  x='avg_log2FC',y='p_val_adj',
  FCcutoff = 0.1,
  colCustom=keyvals, labFace = "italic") + coord_flip(xlim=c(-3, 1))
})
dev.off()

#RVF vs pRV
if (file.exists(.cache_4C_pRV)) {
  message('Loading cached FindMarkers (RVF vs pRV) from Fig4C cache...')
  combined_set <- readRDS(.cache_4C_pRV)
} else {
  combined_set <- data.frame()
  mods_idx <- c(2,12,28,10,25,26)
  for (i in mods_idx){
    key_genes <- subset(bulk_modules,module %in% c(i))$gene_name
    key_genes <- key_genes[key_genes %in% rownames(seurat_ref)]
    gene_set <- FindMarkers(seurat_ref, ident.1 = "pRV", ident.2 = "NF",features=key_genes)
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
  saveRDS(combined_set, .cache_4C_pRV)
}


library(EnhancedVolcano)

combined_set <- combined_set[!grepl('MT-',rownames(combined_set)),]
keyvals <- combined_set$color
names(keyvals)  <- combined_set$module

pdf(paste0('./output/Figure_5/', 'Figure_5_supp_singlecell_volcano_M2_RVF_vs_pRV.pdf'), width=8, height=6)

print({
EnhancedVolcano(combined_set,lab=rownames(combined_set),
  x='avg_log2FC',y='p_val_adj',
  FCcutoff = 0.1,
  colCustom=keyvals, labFace = "italic") + coord_flip()
})
dev.off()


sub_set <- sub_set[!grepl('MT-',rownames(sub_set)),]
sub_set <- subset(combined_set,combined_set$module != 'M2')
keyvals <- sub_set$color
names(keyvals)  <- sub_set$module

pdf(paste0('./output/Figure_5/', 'Figure_5_supp_singlecell_volcano_noM2_RVF_vs_pRV.pdf'), width=8, height=6)

print({
EnhancedVolcano(sub_set,lab=rownames(sub_set),
  x='avg_log2FC',y='p_val_adj',
  FCcutoff = 0.1,
  colCustom=keyvals, labFace = "italic") + coord_flip()
})
dev.off()

subsub_set <- subset(sub_set,sub_set$module != 'M12')
keyvals <- subsub_set$color
names(keyvals)  <- subsub_set$module

pdf(paste0('./output/Figure_5/', 'Figure_5_supp_singlecell_volcano_mito_RVF_vs_pRV.pdf'), width=8, height=6)

print({
EnhancedVolcano(subsub_set,lab=rownames(subsub_set),
  x='avg_log2FC',y='p_val_adj',
  FCcutoff = 0.1,
  colCustom=keyvals, labFace = "italic") + coord_flip(xlim=c(-3, 1))
})
dev.off()



##############################################
##############################################
#### Panel E — Mito-module volcanos (M10/M25/M26/M28)
##############################################
##############################################
Xenium.cm <- readRDS(path.expand(
  './dependencies/shared/Xenium/cm_clean_clean.rds'))
Xenium.cm <- UpdateSeuratObject(Xenium.cm)
Xenium.cm <- NormalizeData(Xenium.cm, verbose = FALSE)

pdf('./output/Figure_5/Xenium/CM_Xenium_Trends.pdf',width=20,height=5)

print({
# Drop the 'sct_' prefix — clean object has no SCT assay; default Xenium assay is fine.
p1 <- VlnPlot(Xenium.cm,'NPPA',group.by='group', pt.size = 0) + NoLegend()
p2 <- VlnPlot(Xenium.cm,'NPPB',group.by='group', pt.size = 0) + NoLegend()
p3 <- VlnPlot(Xenium.cm,'MYH6',group.by='group', pt.size = 0) + NoLegend()
p4 <- VlnPlot(Xenium.cm,'MYH7',group.by='group', pt.size = 0) + NoLegend()
p5 <- VlnPlot(Xenium.cm,'TNNI3',group.by='group', pt.size = 0) + NoLegend()
p6 <- VlnPlot(Xenium.cm,'ANKRD1',group.by='group', pt.size = 0) + NoLegend()
p7 <- VlnPlot(Xenium.cm,'BMPR2',group.by='group', pt.size = 0) + NoLegend()
p8 <- VlnPlot(Xenium.cm,'PANK1',group.by='group', pt.size = 0) + NoLegend()
p9 <- VlnPlot(Xenium.cm,'TMEM65',group.by='group', pt.size = 0) + NoLegend()


p1 | p2 | p4 | p5 | p6 | p7 | p8 | p9 | p3
})
dev.off()


###############################################################################
## NPPA / NPPB Xenium spatial tiles (Figure 5 Panel F regen).
## Uses the same per-patient rasterised tile pattern as Fig 3D, but with a
## continuous gray→red expression scale so low-expressing tissue greys out
## and NPPA/NPPB-high tissue (endocardium) lights up in red.
###############################################################################
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(ggrastr); library(cowplot)
})

## The minimalist Xenium.cm object does not carry spatial coordinates; load
## the full subclustered Xenium object to pull centroids per cell, then subset
## to CMs and join NPPA/NPPB expression.
##
## Wrapped in tryCatch because Seurat v5 sometimes rejects the cached FOV class
## (slots "coords_x_orientation" / "misc" absent on older objects). If the load
## fails we still want Panels D/E/F/G/H/I/J to run downstream.
tryCatch({
.xen_full_path <- './dependencies/shared/Xenium/xenium_obj_subclustered.rds'
if (!file.exists(.xen_full_path)) {
  message('Full Xenium object not found at ', .xen_full_path,
          '; skipping NPPA/NPPB spatial tiles.')
} else {
  message('Loading full Xenium object for NPPA/NPPB spatial tiles...')
  .xen_full <- readRDS(.xen_full_path)
  ## Identify CM cells using whatever lineage column the object carries.
  .lin_col <- intersect(c('Names','lineage','cell_type','Lineage','cell_type_rctd_doublet','cell_type_rctd'),
                         colnames(.xen_full@meta.data))[1]
  if (is.na(.lin_col)) stop('Cannot find lineage column on full Xenium object')
  .cm_cells <- rownames(.xen_full@meta.data)[
    grepl('^CM$|Cardiom', .xen_full@meta.data[[.lin_col]], ignore.case = TRUE)]
  message('  CM cells in full object: ', length(.cm_cells))
  .xen_cm <- subset(.xen_full, cells = .cm_cells)
  rm(.xen_full); gc(verbose = FALSE)

  ## Pull centroids from EVERY FOV (one per patient) and stack them.
  .imgs <- Images(.xen_cm)
  message('  FOVs: ', paste(.imgs, collapse = ', '))
  .coords <- do.call(rbind, lapply(.imgs, function(im) {
    cc <- tryCatch(GetTissueCoordinates(.xen_cm, which = 'centroids', image = im),
                   error = function(e) GetTissueCoordinates(.xen_cm[[im]]))
    if (!('cell' %in% colnames(cc))) {
      cc$cell <- rownames(cc)
    }
    cc
  }))
  .cell_col <- 'cell'
  .x_col <- intersect(c('x','x_centroid'), colnames(.coords))[1]
  .y_col <- intersect(c('y','y_centroid'), colnames(.coords))[1]
  message('  centroids retrieved: ', nrow(.coords),
          ' cells across ', length(.imgs), ' FOVs')

  .md <- .xen_cm@meta.data
  .pt_col <- intersect(c('patient','sample','orig.ident'), colnames(.md))[1]
  if (is.na(.pt_col)) stop('Cannot find patient/sample column on full Xenium object')
  .grp_col <- intersect(c('group','disease','disease_state'), colnames(.md))[1]

  ## FetchData handles BPCells / sparse / dense backends uniformly.
  .feat_df <- tryCatch(
    Seurat::FetchData(.xen_cm, vars = c('NPPA','NPPB'), layer = 'data'),
    error = function(e)
      Seurat::FetchData(.xen_cm, vars = c('NPPA','NPPB'))
  )
  .nppa <- if ('NPPA' %in% colnames(.feat_df)) .feat_df[, 'NPPA'] else NA_real_
  .nppb <- if ('NPPB' %in% colnames(.feat_df)) .feat_df[, 'NPPB'] else NA_real_

  .df <- data.frame(
    cell    = colnames(.xen_cm),
    patient = as.character(.md[[.pt_col]]),
    group   = if (!is.na(.grp_col)) as.character(.md[[.grp_col]]) else NA_character_,
    NPPA    = .nppa,
    NPPB    = .nppb,
    stringsAsFactors = FALSE
  )
  ## Combined natriuretic-peptide score = mean of NPPA + NPPB log-norm expression
  .df$NPP_score <- rowMeans(cbind(.df$NPPA, .df$NPPB), na.rm = TRUE)
  ## Merge in coordinates
  .coords_keep <- data.frame(cell = .coords[[.cell_col]],
                              x = .coords[[.x_col]],
                              y = .coords[[.y_col]])
  .df <- merge(.df, .coords_keep, by = 'cell', all.x = FALSE)
  rm(.xen_cm, .expr); gc(verbose = FALSE)

  .ums_per_inch <- 3000   # match Fig 3D physical scale
  .tile_dir <- file.path(V52_FIG_DIR, 'Figure_5_panel_F_NPPA_NPPB_tiles')
  dir.create(.tile_dir, showWarnings = FALSE, recursive = TRUE)

  ## Per-patient combined NPP-score tiles (independently rasterised, Illustrator-friendly).
  for (.pt in unique(.df$patient)) {
    .sub <- .df[.df$patient == .pt, ]
    .grp <- .sub$group[1]
    .xr  <- range(.sub$x); .yr <- range(.sub$y)
    .w_in <- diff(.xr) / .ums_per_inch
    .h_in <- diff(.yr) / .ums_per_inch
    if (all(is.na(.sub$NPP_score))) next
    ## Plot low-expressing cells first so high-expressing cells render on top.
    .sub2 <- .sub[order(.sub$NPP_score), ]
    .p <- ggplot(.sub2, aes(x = x, y = y, color = NPP_score)) +
      ggrastr::rasterise(geom_point(size = 0.45, stroke = 0), dpi = 400) +
      scale_color_gradient(low = 'grey90', high = 'red', name = 'NPPA/NPPB') +
      coord_fixed(xlim = .xr, ylim = .yr) +
      theme_void() +
      theme(legend.position = 'none', plot.margin = margin(0, 0, 0, 0))
    .fname <- sprintf('Figure_5_panel_F_NPP_%s_%s.pdf', .grp, .pt)
    ggsave(file.path(.tile_dir, .fname), plot = .p,
           width = .w_in, height = .h_in, device = cairo_pdf)
  }

  ## Composite (faceted) overview PDF — combined NPP score.
  if (!all(is.na(.df$NPP_score))) {
    .ord <- .df[order(.df$NPP_score), ]
    .p <- ggplot(.ord, aes(x = x, y = -y, color = NPP_score)) +
      ggrastr::rasterise(geom_point(size = 0.45, stroke = 0), dpi = 300) +
      facet_wrap(~ group + patient, scales = 'free') +
      scale_color_gradient(low = 'grey90', high = 'red',
                            name = 'NPPA/NPPB\n(mean log-norm)') +
      labs(title = 'Combined NPPA/NPPB score across Xenium CM cells',
           x = NULL, y = NULL) +
      theme_void() +
      theme(strip.text = element_text(size = 7),
            legend.position = 'right')
    save_figure(.p, 'Figure_5_supp_NPP_xenium_overview.pdf',
                width = 12, height = 8)
  }
  message('Combined NPPA/NPPB spatial tiles written to ', .tile_dir)
}
}, error = function(e) {
  message('Xenium NPP block skipped: ', conditionMessage(e))
})


###############################################################################
## PSEUDOBULK volcano regeneration for Figure 5 panels D and F.
## The volcanos saved earlier above use cell-level FindMarkers, which inflates
## p-values by treating each nucleus as independent. Here we re-render the same
## panels using sample-level pseudobulk DESeq2 (CM lineage from
## sn_pseudobulk_lineage_deseq2.csv) and overwrite the existing PDFs so the
## composite picks up the statistically defensible version.
###############################################################################
suppressPackageStartupMessages({
  library(EnhancedVolcano)
  library(dplyr)
})

## --- Recompute CM-lineage pseudobulk DESeq2 with lfcShrink (apeglm) -------
## The shared CSV stores unshrunken MLE log2FC, which produces extreme values
## for sparse genes (e.g., NTM lfc ~ -25). For volcano display we recompute
## with apeglm shrinkage (cached).
.cache_cm_pb_shrunk <- './output/Figure_5/fig5_cm_pseudobulk_shrunk.rds'
if (file.exists(.cache_cm_pb_shrunk)) {
  message('Loading cached CM pseudobulk DESeq2 (apeglm-shrunk)...')
  .sn_pb <- readRDS(.cache_cm_pb_shrunk)
} else {
  message('Computing CM pseudobulk DESeq2 with apeglm lfcShrink...')
  suppressPackageStartupMessages({
    library(Seurat); library(DESeq2)
  })
  if (!exists('sn') || !inherits(sn, 'Seurat')) {
    sn_full <- readRDS('./dependencies/shared/snRV_ref.rds')
  } else { sn_full <- sn }
  cm <- subset(sn_full, subset = Names == 'CM')
  agg <- AggregateExpression(cm, group.by = c('patient','group'),
                              assays = DefaultAssay(cm), slot = 'counts',
                              return.seurat = FALSE)[[DefaultAssay(cm)]]
  ## Build coldata from the column names AggregateExpression produced
  cn <- colnames(agg)
  cd <- data.frame(sample = cn,
                   group  = sub('.*_', '', cn),
                   patient = sub('_[^_]+$', '', cn),
                   stringsAsFactors = FALSE)
  cd$group <- factor(cd$group, levels = c('NF','pRV','RVF'))
  dds <- DESeqDataSetFromMatrix(countData = round(as.matrix(agg)),
                                  colData   = cd, design = ~ group)
  dds <- DESeq(dds)
  ## Use 'ashr' shrinkage (supports arbitrary contrasts without re-fitting).
  contrasts <- list(
    pRV_vs_NF  = c('group','pRV','NF'),
    RVF_vs_NF  = c('group','RVF','NF'),
    RVF_vs_pRV = c('group','RVF','pRV'))
  out_list <- lapply(names(contrasts), function(lab) {
    res <- tryCatch(
      lfcShrink(dds, contrast = contrasts[[lab]], type = 'ashr'),
      error = function(e) {
        message('  ashr failed for ', lab, '; falling back to normal: ', e$message)
        lfcShrink(dds, contrast = contrasts[[lab]], type = 'normal')
      })
    df  <- as.data.frame(res)
    df$gene <- rownames(df); df$contrast <- lab; df$subtype <- 'CM'
    df
  })
  .sn_pb <- do.call(rbind, out_list)
  saveRDS(.sn_pb, .cache_cm_pb_shrunk)
  rm(cm, agg, dds); gc()
}
.sn_pb <- subset(.sn_pb, subtype == 'CM')
.bulk_mods <- read.csv('./dependencies/shared/bulk_heart_modules.csv',
                        stringsAsFactors = FALSE)
.color2M <- function(color) paste0('M', match(color, WGCNA::labels2colors(1:100)))
.bulk_mods$M_id   <- .color2M(.bulk_mods$module)
.bulk_mods$mod_color <- .bulk_mods$module

## Helper: build a per-gene volcano frame for a CM contrast restricted to
## genes belonging to the requested set of M-numbers. Output mimics the
## EnhancedVolcano-friendly schema used elsewhere in the script (avg_log2FC,
## p_val_adj, module, color), so we can reuse the same plot call.
.pb_volcano_frame <- function(contrast_lab, m_ids) {
  ann <- subset(.bulk_mods, M_id %in% m_ids,
                select = c('gene_name','M_id','mod_color'))
  ann <- ann[match(unique(ann$gene_name), ann$gene_name), ]
  d   <- subset(.sn_pb, contrast == contrast_lab)
  d   <- d[!is.na(d$padj) & d$baseMean > 50, ]
  d   <- merge(d, ann, by.x = 'gene', by.y = 'gene_name')
  d   <- d[!grepl('^MT-', d$gene), ]
  data.frame(
    avg_log2FC = d$log2FoldChange,
    p_val_adj  = d$padj,
    module     = d$M_id,
    color      = d$mod_color,
    row.names  = d$gene,
    stringsAsFactors = FALSE)
}

## --- Panel D: M2 (blue) module, CM pseudobulk, pRV vs NF (Phase 1) -------
## M2 is induced primarily in Phase 1. Coloring: significant (padj<0.05 &
## |LFC|>0.5) red, NS gray. Specific biologically interesting genes labeled.
.M2_pb <- .pb_volcano_frame('pRV_vs_NF', 'M2')
if (nrow(.M2_pb) > 0) {
  .sig <- !is.na(.M2_pb$p_val_adj) &
          .M2_pb$p_val_adj < 0.05 &
          abs(.M2_pb$avg_log2FC) > 0.5
  .keyvals <- ifelse(.sig, '#E41A1C', '#BDBDBD')
  names(.keyvals) <- ifelse(.sig, 'Significant', 'NS')
  .label_genes <- intersect(
    c('EDNRA','KCNE1','ACTA1','MYL2','HTRA1','HSPB6','COL16A1','OXA1L'),
    rownames(.M2_pb))
  pdf('./output/Figure_5/Figure_5_panel_D_M2_pseudobulk_volcano_pRV_vs_NF.pdf', width = 8, height = 6)
  print(
    EnhancedVolcano(.M2_pb, lab = rownames(.M2_pb),
      x = 'avg_log2FC', y = 'p_val_adj',
      FCcutoff = 0.5, pCutoff = 0.05,
      xlim = c(-5, 5), ylim = c(0, 4),
      selectLab = .label_genes, max.overlaps = Inf, drawConnectors = TRUE,
      arrowheads = FALSE, widthConnectors = 0.3, caption = '',
      title = 'M2 module (CM pseudobulk, pRV vs NF)',
      subtitle = 'sample-level DESeq2 with ashr shrinkage',
      legendPosition = 'top',
      colCustom = .keyvals, labFace = 'italic') + coord_flip()
  )
  dev.off()
  message('Panel D: pseudobulk volcano written (', nrow(.M2_pb), ' genes; ',
          sum(.sig), ' significant; ', length(.label_genes), ' labeled).')
}

## --- Panel F: mitochondrial modules M10/M25/M26/M28, CM pseudobulk -------
## Phase 1 (NF→pRV) is the canonical suppression direction for these modules.
.mito_pb <- .pb_volcano_frame('pRV_vs_NF', c('M10','M25','M26','M28'))
if (nrow(.mito_pb) > 0) {
  ## Per-gene colors but legend-named by module (4 entries: M10/M25/M26/M28).
  .keyvals <- setNames(.mito_pb$color, .mito_pb$module)
  ## Label top 10 MitoCarta3.0 genes by p-value (any mito gene, not just ETC)
  .mc_path <- './dependencies/shared/Human.MitoCarta3.0.xls'
  .mc_genes <- tryCatch({
    suppressPackageStartupMessages(library(readxl))
    sheet <- readxl::read_excel(.mc_path, sheet = 'A Human MitoCarta3.0')
    sym_col <- intersect(c('Symbol', 'Gene Symbol', 'gene_symbol', 'Gene'), colnames(sheet))[1]
    unique(toupper(as.character(sheet[[sym_col]])))
  }, error = function(e) {
    message('  MitoCarta load failed (', e$message, '); falling back to ETC/mitoribosome regex.')
    grep('^(NDUF|COX[0-9]|SDH|UQCR|ATP5|MRPS|MRPL|COA|TIMM|TOMM)',
         rownames(.mito_pb), value = TRUE)
  })
  .mito_only <- intersect(rownames(.mito_pb), .mc_genes)
  .top10_p1 <- .mito_only[order(.mito_pb[.mito_only, 'p_val_adj'])][seq_len(min(10, length(.mito_only)))]
  message('  MitoCarta filter: ', length(.mito_only), ' / ', nrow(.mito_pb),
          ' module genes are MitoCarta3.0; labeling top 10 by padj.')
  pdf('./output/Figure_5/Figure_5_panel_E_mito_pseudobulk_volcano.pdf', width = 8, height = 6)
  print(
    EnhancedVolcano(.mito_pb, lab = rownames(.mito_pb),
      x = 'avg_log2FC', y = 'p_val_adj',
      FCcutoff = 0.5, pCutoff = 0.05,
      xlim = c(-1, 1), ylim = c(0, 2.5),
      selectLab = .top10_p1, max.overlaps = Inf, drawConnectors = TRUE,
      arrowheads = FALSE, widthConnectors = 0.3, caption = '',
      title = 'Mitochondrial modules (CM pseudobulk, pRV vs NF)',
      subtitle = 'M10/M25/M26/M28 — sample-level DESeq2',
      legendPosition = 'top',
      colCustom = .keyvals, labFace = 'italic') + coord_flip()
  )
  dev.off()
  message('Panel E: mito-module pseudobulk volcano written (', nrow(.mito_pb), ' Phase-1 genes).')
}


#M1 <- readRDS(file = "./dependencies/shared/cm_new_subclust.rds")

















####### OLD















#######################################
##  LEGACY (was Panel 5F in v52): MitoCarta3.0 module score per CM subcluster
##  In v57 this dotplot lives in Supplementary Figure 2 (panel G). Kept here
##  for reference; v57 Panel F is CollecTRI/decoupleR TF activity (TODO: add).
##  X = CM subclusters from Subnames_manual (CM_unc removed)
##  Y = disease group (NF / pRV / RVF)
##  Color = z(mean MitoCarta_All AddModuleScore), Size = % cells with score > 0
###############################################################################

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(ggplot2)
})

## Defensive: re-source helpers if running this block in isolation
if (!exists('theme_v52', mode = 'function')) {
  source('./helper_scripts/_shared_helpers.R')
}

.cache_panel_5F_mito <- './output/Figure_5/fig5_panel_5F_mitocarta_all_dotplot.rds'

if (file.exists(.cache_panel_5F_mito)) {
  message('Loading cached Panel 5F MitoCarta_All dotplot frame...')
  agg_mito <- readRDS(.cache_panel_5F_mito)
} else {
  message('Computing Panel 5F MitoCarta_All module score on snRV_ref CMs...')

  ## --- Load CM subset of full snRNA reference --------------------------------
  if (!exists('snRV_ref') || !inherits(snRV_ref, 'Seurat')) {
    snRV_ref <- readRDS('./dependencies/shared/snRV_ref.rds')
  }
  cm_obj <- subset(snRV_ref, subset = Names == 'CM')
  if (!'Subnames_manual' %in% colnames(cm_obj@meta.data)) {
    stop('snRV_ref CM subset is missing `Subnames_manual` — re-run annotate_subnames.r')
  }
  ## Drop unclassified CMs
  cm_obj <- subset(cm_obj, subset = Subnames_manual != 'CM_unc')

  ## --- Load MitoCarta3.0 (same loader as Panel F volcano) -------------------
  .mc_path <- './dependencies/shared/Human.MitoCarta3.0.xls'
  .mc_genes_all <- tryCatch({
    suppressPackageStartupMessages(library(readxl))
    sheet <- readxl::read_excel(.mc_path, sheet = 'A Human MitoCarta3.0')
    sym_col <- intersect(c('Symbol','Gene Symbol','gene_symbol','Gene'),
                         colnames(sheet))[1]
    unique(toupper(as.character(sheet[[sym_col]])))
  }, error = function(e) {
    message('  MitoCarta load failed (', e$message, '); falling back to ETC regex universe.')
    grep('^(NDUF|COX[0-9]|SDH|UQCR|ATP5|MRPS|MRPL|COA|TIMM|TOMM)',
         rownames(cm_obj), value = TRUE)
  })

  .universe   <- toupper(rownames(cm_obj))
  .mc_inscope <- intersect(.mc_genes_all, .universe)
  message('Panel 5F: ', length(.mc_inscope), ' / ', length(.mc_genes_all),
          ' MitoCarta3.0 genes present in CM object.')

  ## --- AddModuleScore (single set: all MitoCarta) ----------------------------
  cm_obj <- AddModuleScore(cm_obj, features = list(.mc_inscope),
                           name = 'mitoAll_')
  ## AddModuleScore appends as 'mitoAll_1'

  ## --- Aggregate per subtype × group ----------------------------------------
  md <- cm_obj@meta.data[, c('Subnames_manual','group','mitoAll_1')]
  md$group           <- factor(md$group,
                               levels = intersect(c('NF','pRV','RVF'),
                                                  unique(as.character(md$group))))
  md$Subnames_manual <- factor(md$Subnames_manual)

  agg_mito <- md %>%
    dplyr::group_by(group, Subnames_manual) %>%
    dplyr::summarise(mean_score = mean(mitoAll_1),
                     pct_pos    = 100 * mean(mitoAll_1 > 0),
                     n_cells    = dplyr::n(),
                     .groups    = 'drop') %>%
    dplyr::mutate(z = as.numeric(scale(mean_score)))

  saveRDS(agg_mito, .cache_panel_5F_mito)
  rm(cm_obj, md); gc(verbose = FALSE)
}

## --- Plot: subtype × group dotplot, MitoCarta_All only --------------------
p_5F_mito <- ggplot(agg_mito,
                    aes(x = Subnames_manual, y = group,
                        colour = z, size = pct_pos)) +
  geom_point() +
  scale_colour_gradient2(low = '#2166AC', mid = 'grey92', high = '#B2182B',
                         midpoint = 0, name = 'z(mean score)') +
  scale_size_area(max_size = 6, name = '% cells > 0') +
  labs(x = 'CM subcluster', y = NULL,
       title = 'MitoCarta3.0 module score (all genes) by CM subcluster × group') +
  theme_v52(COMP_W) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

save_figure(p_5F_mito, 'Figure_5_supp_mitocarta_all_dotplot.pdf',
            width = COMP_W, height = 3.0)
message('Figure 5 Panel 5F complete (MitoCarta_All dotplot, Subnames_manual × group).')


###############################################################################
## Panel F — Subject-level pan-CM mitochondrial-TF activity heatmap.
## Faithful reimplementation of additional_scripts/mito_TF.R "patient
## breakdown": CollecTRI ULM activity (run_ulm) over ALL CM cells, z-scored
## per TF (Seurat ScaleData equivalent), averaged PER PATIENT, shown for the
## 7 mito regulators, patients ordered NF -> pRV -> RVF. The expensive
## run_ulm is read from the committed cache dependencies/shared/TF_activity.rds
## (the exact mito_TF.R output); cell -> patient/group comes from seurat_ref
## (cm_subclust_new_new, already loaded; 100% barcode overlap).
###############################################################################
## Image row order (top -> bottom): NFE2L2, NRF1, PPARGC1A, PPARD, PPARA,
## ESRRG, ESRRA. ggplot puts the first y factor level at the bottom, so list
## bottom -> top.
.panel_F_mito_tfs <- c('NFE2L2','NRF1','PPARGC1A','PPARD','PPARA','ESRRG','ESRRA')
.panel_F_tf_levels <- rev(.panel_F_mito_tfs)        # ESRRA ... NFE2L2 (bottom->top)
.tf_activity_rds <- './dependencies/shared/TF_activity.rds'

panel_F_df <- NULL
if (!file.exists(.tf_activity_rds)) {
  message('Panel F: ', .tf_activity_rds, ' missing — skipping.')
} else if (!exists('seurat_ref')) {
  message('Panel F: seurat_ref not loaded — skipping.')
} else {
  message('Panel F: building subject-level pan-CM mito-TF heatmap from ',
          .tf_activity_rds, ' ...')
  panel_F_df <- tryCatch({
    acts <- readRDS(.tf_activity_rds)
    acts <- acts[acts$statistic == 'ulm' &
                   acts$source %in% .panel_F_mito_tfs, , drop = FALSE]
    md <- seurat_ref@meta.data
    md$.cell <- rownames(md)
    .pg <- unique(data.frame(.cell   = md$.cell,
                             patient = as.character(md$patient),
                             group   = as.character(md$group),
                             stringsAsFactors = FALSE))
    acts <- acts[acts$condition %in% .pg$.cell, , drop = FALSE]
    ## TF x cell score matrix (base R; small after the TF filter).
    .cells <- unique(acts$condition)
    M <- matrix(NA_real_, length(.panel_F_mito_tfs), length(.cells),
                dimnames = list(.panel_F_mito_tfs, .cells))
    M[cbind(match(acts$source, .panel_F_mito_tfs),
            match(acts$condition, .cells))] <- acts$score
    ## Per-TF z-score across cells (Seurat ScaleData equivalent).
    Mz <- t(scale(t(M)))
    Mz[!is.finite(Mz)] <- 0
    ## Per-patient mean of the scaled activity.
    .cellpat <- setNames(.pg$patient, .pg$.cell)[.cells]
    agg <- vapply(split(seq_along(.cells), .cellpat),
                  function(ix) rowMeans(Mz[, ix, drop = FALSE], na.rm = TRUE),
                  numeric(length(.panel_F_mito_tfs)))
    rownames(agg) <- .panel_F_mito_tfs
    ## Long frame + patient ordering NF -> pRV -> RVF (sorted within group).
    .p2g <- setNames(.pg$group, .pg$patient)
    long <- data.frame(
      tf      = rep(rownames(agg), times = ncol(agg)),
      patient = rep(colnames(agg), each  = nrow(agg)),
      value   = as.vector(agg), stringsAsFactors = FALSE)
    long$group <- factor(.p2g[long$patient], levels = c('NF','pRV','RVF'))
    .porder <- unlist(lapply(c('NF','pRV','RVF'),
                             function(g) sort(unique(long$patient[long$group == g]))))
    long$patient <- factor(long$patient, levels = .porder)
    long$tf      <- factor(long$tf,      levels = .panel_F_tf_levels)
    long[!is.na(long$patient) & !is.na(long$group), ]
  }, error = function(e) {
    message('Panel F build failed (', conditionMessage(e),
            '); skipping so downstream G/H/I/J still run.')
    NULL
  })
}

if (!is.null(panel_F_df) && nrow(panel_F_df) > 0) {
  p_5F <- ggplot(panel_F_df,
                 aes(x = patient, y = tf, fill = value)) +
          geom_tile(colour = 'white', linewidth = 0.3) +
          scale_fill_gradient2(low = '#2166ac', mid = 'white',
                                high = '#b2182b', midpoint = 0,
                                limits = c(-2, 2),
                                oob = scales::squish,
                                name = expression(-Log[10]*P)) +
          facet_grid(~ group, scales = 'free_x', space = 'free_x') +
          labs(x = NULL, y = NULL,
               title = 'Subject-level pan-CM mito-TF activity') +
          theme_v52(COMP_W) +
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                panel.spacing.x = unit(0.15, 'lines'))
  save_figure(p_5F, 'Figure_5_panel_F_collectri_TF_activity.pdf',
              width = 7, height = 3)
  message('Figure 5 Panel F complete (subject-level pan-CM mito-TF heatmap, ',
          nlevels(panel_F_df$patient), ' patients).')
} else {
  message('Panel F: no data — emitted no PDF.')
}


###############################################################################
## Panel G — Pseudo-bulk violin plots of ESRRA, ESRRG, PPARA, PPARGC1A,
##           NRF1, NFE2L2 expression in cardiomyocytes across NF / pRV / RVF.
## Per-patient DESeq2-VST aggregation (one point = one donor); pairwise Wilcoxon.
## Cache-aware: uses ./output/Figure_5/fig5_panel_G_pseudobulk.rds when present.
###############################################################################
.cache_panel_G <- './output/Figure_5/fig5_panel_G_pseudobulk.rds'
.panel_G_genes <- c('ESRRA','ESRRG','PPARA','PPARD','PPARGC1A','NRF1','NFE2L2')

if (file.exists(.cache_panel_G)) {
  message('Loading cached Panel G pseudobulk VST...')
  panel_G_df <- readRDS(.cache_panel_G)
} else if (exists('seurat_ref')) {
  message('Computing Panel G CM pseudobulk DESeq2-VST...')
  suppressPackageStartupMessages({ library(DESeq2); library(SummarizedExperiment) })
  cm_obj <- subset(seurat_ref, subset = Names == 'CM')
  DefaultAssay(cm_obj) <- 'RNA'
  agg <- Seurat::AggregateExpression(cm_obj, assays = 'RNA', slot = 'counts',
                                      group.by = c('patient','group'),
                                      return.seurat = FALSE)$RNA
  cn <- colnames(agg)
  cd <- data.frame(sample = cn,
                   patient = sub('_[^_]+$', '', cn),
                   group   = sub('.*_', '', cn))
  rownames(cd) <- cn
  cd$group <- factor(cd$group, levels = c('NF','pRV','RVF'))
  dds <- DESeqDataSetFromMatrix(countData = round(agg), colData = cd,
                                 design = ~ group)
  vsd <- vst(dds, blind = FALSE)
  vst_mat <- assay(vsd)
  use_genes <- intersect(.panel_G_genes, rownames(vst_mat))
  long <- do.call(rbind, lapply(use_genes, function(g) {
    data.frame(gene = g,
               sample = colnames(vst_mat),
               patient = cd[colnames(vst_mat),'patient'],
               group = cd[colnames(vst_mat),'group'],
               expr = vst_mat[g, ])
  }))
  panel_G_df <- long
  dir.create(dirname(.cache_panel_G), showWarnings = FALSE, recursive = TRUE)
  saveRDS(panel_G_df, .cache_panel_G)
  rm(cm_obj, agg, dds, vsd, vst_mat); gc(verbose = FALSE)
} else {
  message('seurat_ref not loaded; skipping Panel G.')
  panel_G_df <- NULL
}

if (!is.null(panel_G_df) && nrow(panel_G_df) > 0) {
  p_5G <- ggplot(panel_G_df,
                 aes(x = group, y = expr, fill = group)) +
          geom_violin(scale = 'width', trim = FALSE, colour = NA, alpha = 0.6) +
          geom_boxplot(width = 0.18, fill = 'white', outlier.shape = NA,
                        colour = 'grey20', linewidth = 0.3) +
          geom_jitter(width = 0.08, size = 0.7, alpha = 0.7, colour = 'grey25') +
          facet_wrap(~ gene, nrow = 1, scales = 'free_y') +
          scale_fill_manual(values = c(NF = '#4DAF4A', pRV = '#377EB8', RVF = '#E41A1C')) +
          labs(x = NULL, y = 'VST expression (CM pseudobulk)',
               title = 'CM mito TF expression (per-patient pseudobulk)') +
          theme_v52(COMP_W) +
          theme(legend.position = 'none')
  save_figure(p_5G, 'Figure_5_panel_G_mito_TF_pseudobulk_violins.pdf',
              width = COMP_W, height = 3.0)
  message('Figure 5 Panel G complete (CM mito-TF pseudobulk violins).')
} else {
  message('Panel G: no data — emitted no PDF.')
}


###############################################################################
## Figure 5 Panels H–J: respirometry across adult, PAB (2wk + 2mo), pediatric
##
##   G — Adult RV (NF/pRV/RVF):   CI, CI+II, CII, CIV; ROX-corrected;
##                                pmoles O2/sec/mg tissue; donor dots.
##   H — Murine PAB:              2-wk and 2-mo subpanels (ADP, OC, Succ + RCR);
##                                substrate-addition protocol (no ROX correction
##                                possible for this cohort); animal dots.
##   I — Pediatric (ped-NF / HLHS-failing): same 4-complex layout as G;
##                                          ROX-corrected; donor dots.
##   J — Cross-cohort summary: percent-of-control bar of the headline coupled
##       OXPHOS measure (CI+II for adult/peds, CI+II+FAO/Succ for PAB).
##
## Data source: dependencies/shared/mito_data.xlsx — sheets:
##   Adult_Data, PAB_2_Week_Data, PAB_8_Week_Data, Peds_Data
##
## v52-named output PDFs → v54 final-panel mapping (assembled in Illustrator):
##   Figure_5_panel_H_adult_respirometry.pdf       → top sub-plot of v54 5H
##   Figure_5_panel_H_pediatric_respirometry.pdf   → right sub-plot of v54 5H
##   Figure_5_panel_I_PAB_respirometry.pdf         → v54 5I
##   Figure_5_panel_J_cross_cohort_summary.pdf     → v54 5J
## (Internal variable names p_5G/p_5H_*/p_5I/p_5J retain the v52 layout.)
###############################################################################

suppressPackageStartupMessages({
  library(readxl); library(ggplot2); library(dplyr); library(tidyr)
  library(rstatix); library(ggpubr); library(patchwork)
})

.mito_xlsx <- './dependencies/shared/mito_data.xlsx'

GROUP_COLS <- c('NF' = '#4DAF4A',  'pRV' = '#FF7F00', 'RVF' = '#E41A1C',
                'Sham' = '#4DAF4A','PAB_Mod' = '#FF7F00','PAB_Severe' = '#E41A1C',
                'PAB_Sev'='#E41A1C',
                'ped_NF' = '#4DAF4A','ped_HLHS_fail' = '#E41A1C')

## ── helper: pull a long-form per-sample data.frame from a sheet ──────────────
.read_resp_sheet <- function(sheet) {
  raw <- suppressMessages(read_excel(.mito_xlsx, sheet = sheet, col_names = TRUE))
  ## Avoid `dplyr::rename(state = !!state_col)` — the !! unquote operator
  ## is being shadowed by another loaded package in this script's environment
  ## ("Error in !state_col : invalid argument type"). Rename column 1 directly.
  names(raw)[1] <- 'state'
  long <- raw %>%
    tidyr::pivot_longer(-state, names_to = 'sample', values_to = 'value') %>%
    dplyr::filter(!is.na(state) & !is.na(sample) & !is.na(value))
  long$state <- gsub(' ', ' ', long$state)   # NBSP → space ('Cyt C')
  long$state <- trimws(long$state)
  long
}

## ── helper: compute ROX-corrected complex values per sample ─────────────────
##   complex_map: state row → complex label (e.g. NADH → 'CI'). AmA is treated
##   as ROX baseline; all others are subtracted by AmA per matched sample.
.rox_correct <- function(long, complex_map) {
  rox <- long %>% filter(state == 'AmA') %>%
    select(sample, rox = value)
  long %>% inner_join(rox, by = 'sample') %>%
    filter(state %in% names(complex_map)) %>%
    mutate(complex = unname(complex_map[state]),
           value_corr = value - rox)
}

## ── helper: assign group label from sample-name prefix ──────────────────────
.classify_group <- function(samples, levels_map) {
  out <- rep(NA_character_, length(samples))
  for (lbl in names(levels_map)) {
    out[grepl(levels_map[[lbl]], samples)] <- lbl
  }
  out
}

## ── helper: per-complex KW + Mann-Whitney pairwise on a long data.frame ─────
.complex_stats <- function(df, value_col, complex_col, group_col, group_levels) {
  res <- list()
  for (cx in unique(df[[complex_col]])) {
    sub <- df[df[[complex_col]] == cx, ]
    grp_vals <- split(sub[[value_col]], factor(sub[[group_col]], levels = group_levels))
    grp_vals <- grp_vals[lengths(grp_vals) > 0]
    if (length(grp_vals) < 2) next
    kw <- tryCatch(kruskal.test(grp_vals)$p.value, error = function(e) NA_real_)
    pw <- list()
    grps <- names(grp_vals)
    for (i in seq_along(grps)) for (j in seq_along(grps)) if (i < j) {
      pn <- paste(grps[j], 'vs', grps[i])
      pw[[pn]] <- tryCatch(wilcox.test(grp_vals[[i]], grp_vals[[j]])$p.value,
                           error = function(e) NA_real_)
    }
    res[[cx]] <- list(kw = kw, pairwise = pw,
                      mean_per_grp = sapply(grp_vals, mean))
  }
  res
}

## ── helper: significance-stars from a p-value ───────────────────────────────
.pstar <- function(p) ifelse(is.na(p), 'ns',
                              ifelse(p < 0.001, '***',
                                     ifelse(p < 0.01,  '**',
                                            ifelse(p < 0.05, '*', 'ns'))))


#######################################
#############  FIGURE 5H (adult, top sub-plot)  #############
#######################################
## Adult RV respirometry (NF / pRV / RVF), ROX-corrected, four complexes.

message('Figure 5G: adult RV respirometry...')
adult_long <- .read_resp_sheet('Adult_Data')
adult_complex_map <- c('NADH' = 'CI', 'Succ' = 'CI+II',
                        'Rot'  = 'CII', 'Asc/TMPD' = 'CIV')
adult_corr <- .rox_correct(adult_long, adult_complex_map)
adult_corr$group <- .classify_group(adult_corr$sample,
                                    list(NF  = '^NF_',
                                         pRV = '^pRV_',
                                         RVF = '^RVF_'))
adult_corr <- adult_corr %>% filter(!is.na(group))
adult_corr$group   <- factor(adult_corr$group, levels = c('NF','pRV','RVF'))
adult_corr$complex <- factor(adult_corr$complex,
                              levels = c('CI','CII','CI+II','CIV'))

adult_stats <- .complex_stats(adult_corr, 'value_corr', 'complex', 'group',
                              c('NF','pRV','RVF'))
adult_kw_lab <- sapply(levels(adult_corr$complex),
                       function(cx) sprintf('%s\nKW p=%.3f', cx,
                                            adult_stats[[cx]]$kw))

.respiro_theme <- function() {
  theme_v52(COMP_W) +
  theme(legend.position  = 'none',
        axis.text.x      = element_text(size = 31),
        axis.text.y      = element_text(size = 31),
        axis.title.y     = element_text(size = 34),
        strip.text       = element_text(face = 'bold', size = 34),
        plot.title       = element_text(face = 'bold', size = 36),
        plot.tag         = element_text(face = 'bold', size = 43),
        panel.spacing.x  = unit(0.25, 'lines'))
}

p_5G <- ggplot(adult_corr,
               aes(x = group, y = value_corr, fill = group)) +
  geom_boxplot(width = 0.75, outlier.shape = NA, alpha = 0.55) +
  geom_jitter(width = 0.12, size = 2.625, alpha = 0.85, color = 'black') +
  facet_wrap(~ complex, nrow = 1, scales = 'free_y',
             labeller = as_labeller(adult_kw_lab)) +
  scale_fill_manual(values = GROUP_COLS) +
  ylab('Respiration (pmoles O₂ / sec / mg tissue)\nROX-corrected') +
  xlab(NULL) +
  .respiro_theme()
save_figure(p_5G, 'Figure_5_panel_H_adult_respirometry.pdf',
            width = 18, height = 7)
message('  wrote Figure_5_panel_H_adult_respirometry.pdf')


#######################################
#############  FIGURE 5I (PAB respirometry)  #############
#######################################
## Murine PAB respirometry (Sham / PAB_Mod / PAB_Severe) at 2-wk and 2-mo.
## Substrate-addition protocol: PM, ADP, Cyt C, OC, Succ.
## Headline display states: ADP (CI-coupled), OC (CI+FAO), Succ (CI+II+FAO).

message('Figure 5H: PAB respirometry (2-wk + 2-mo)...')
.pab_panel <- function(sheet, label) {
  long <- .read_resp_sheet(sheet)
  long$group <- .classify_group(long$sample,
                                list(Sham       = '^Sham_',
                                     PAB_Mod    = '^PAB_Mod_',
                                     PAB_Severe = '^PAB_Sev'))
  long <- long %>% filter(!is.na(group),
                          state %in% c('ADP','OC','Succ'))
  long$state <- factor(long$state,
                        levels = c('ADP','OC','Succ'),
                        labels = c('CI-coupled\n(ADP)','CI+FAO\n(OC)',
                                   'CI+II+FAO\n(Succ)'))
  long$group <- factor(long$group,
                        levels = c('Sham','PAB_Mod','PAB_Severe'))

  # KW per state
  kw <- long %>% group_by(state) %>%
    summarise(p = kruskal.test(value ~ group)$p.value, .groups='drop')
  kw_lab <- setNames(sprintf('%s\nKW p=%.3f', levels(long$state), kw$p),
                     levels(long$state))

  ggplot(long, aes(x = group, y = value, fill = group)) +
    geom_boxplot(width = 0.75, outlier.shape = NA, alpha = 0.55) +
    geom_jitter(width = 0.12, size = 2.44, alpha = 0.85, color = 'black') +
    facet_wrap(~ state, nrow = 1, scales = 'free_y',
               labeller = as_labeller(kw_lab)) +
    scale_fill_manual(values = GROUP_COLS) +
    labs(title = label, x = NULL,
         y = expression(atop('Respiration (pmoles O'[2]*' / sec / mg tissue)',
                              ''))) +
    .respiro_theme() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 31))
}
p_5H_2wk  <- .pab_panel('PAB_2_Week_Data', '2 weeks post-PAB')
p_5H_2mo  <- .pab_panel('PAB_8_Week_Data', '2 months post-PAB')
p_5H_full <- p_5H_2wk / p_5H_2mo

save_figure(p_5H_full, 'Figure_5_panel_I_PAB_respirometry.pdf',
            width = 15, height = 17)
message('  wrote Figure_5_panel_I_PAB_respirometry.pdf')


#######################################
#############  FIGURE 5H (pediatric, right sub-plot)  #############
#######################################
## Pediatric RV respirometry (ped-NF / HLHS-palliated failing), four complexes.

message('Figure 5I: pediatric RV respirometry...')
peds_long <- .read_resp_sheet('Peds_Data')
peds_complex_map <- c('NADH' = 'CI', 'Succ' = 'CI+II',
                       'Rot'  = 'CII', 'Asc/TMPD' = 'CIV')
peds_corr <- .rox_correct(peds_long, peds_complex_map)
peds_corr$group <- .classify_group(peds_corr$sample,
                                    list(ped_NF        = '^peds_NF_',
                                         ped_HLHS_fail = '^peds_HLHS_failing_'))
peds_corr <- peds_corr %>% filter(!is.na(group))
peds_corr$group <- factor(peds_corr$group,
                           levels = c('ped_NF','ped_HLHS_fail'),
                           labels = c('ped-NF','HLHS-failing'))
peds_corr$complex <- factor(peds_corr$complex,
                             levels = c('CI','CII','CI+II','CIV'))

peds_stats <- .complex_stats(peds_corr, 'value_corr', 'complex', 'group',
                             c('ped-NF','HLHS-failing'))
peds_kw_lab <- sapply(levels(peds_corr$complex), function(cx) {
  pwn <- names(peds_stats[[cx]]$pairwise)[1]
  sprintf('%s\nMW p=%.3f', cx, peds_stats[[cx]]$pairwise[[pwn]])
})

p_5I <- ggplot(peds_corr, aes(x = group, y = value_corr, fill = group)) +
  geom_boxplot(width = 0.75, outlier.shape = NA, alpha = 0.55) +
  geom_jitter(width = 0.12, size = 2.625, alpha = 0.85, color = 'black') +
  facet_wrap(~ complex, nrow = 1, scales = 'free_y',
             labeller = as_labeller(peds_kw_lab)) +
  scale_fill_manual(values = c('ped-NF' = '#4DAF4A',
                                'HLHS-failing' = '#E41A1C')) +
  ylab('Respiration (pmoles O₂ / sec / mg tissue)\nROX-corrected') +
  xlab(NULL) +
  .respiro_theme() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 31))
save_figure(p_5I, 'Figure_5_panel_H_pediatric_respirometry.pdf',
            width = 18, height = 7)
message('  wrote Figure_5_panel_H_pediatric_respirometry.pdf')


#######################################
#############  FIGURE 5J  #############
#######################################
## Cross-cohort summary: percent-of-control of the headline coupled OXPHOS
## measure per cohort/group, with within-cohort Mann-Whitney significance.
##
##   Cohort         | headline measure          | reference
##   adult-pRV/RVF  | CI+II (Succ; ROX-corr)    | adult NF mean
##   PAB-Sev/Mod 2w | CI+II+FAO (Succ raw)      | Sham 2w mean
##   PAB-Sev/Mod 2m | CI+II+FAO (Succ raw)      | Sham 2m mean
##   ped-HLHS-fail  | CI+II (Succ; ROX-corr)    | ped-NF mean

message('Figure 5J: cross-cohort percent-of-control summary...')

.cross_block <- function(df, value_col, group_col, ref_group, label_levels) {
  ref_mean <- mean(df[[value_col]][df[[group_col]] == ref_group], na.rm = TRUE)
  df$pct <- 100 * (df[[value_col]] / ref_mean - 1)
  df$label <- factor(df[[group_col]], levels = label_levels)
  df
}

# Adult: CI+II
adult_block <- adult_corr %>%
  filter(complex == 'CI+II') %>%
  mutate(cohort = 'Adult', tag = paste('Adult', as.character(group)),
         value = value_corr)
adult_block <- .cross_block(adult_block, 'value', 'group', 'NF',
                             c('NF','pRV','RVF'))
adult_block$display <- factor(paste('Adult', adult_block$group),
                               levels = c('Adult NF','Adult pRV','Adult RVF'))

# Pediatric: CI+II
peds_block <- peds_corr %>%
  filter(complex == 'CI+II') %>%
  mutate(cohort = 'Pediatric', value = value_corr)
peds_ref <- mean(peds_block$value[peds_block$group == 'ped-NF'], na.rm = TRUE)
peds_block$pct <- 100 * (peds_block$value / peds_ref - 1)
peds_block$display <- factor(ifelse(peds_block$group == 'ped-NF',
                                     'Ped NF', 'Ped HLHS-fail'),
                              levels = c('Ped NF','Ped HLHS-fail'))

# PAB 2-wk and 2-mo: Succ
.pab_block <- function(sheet, tp_label) {
  long <- .read_resp_sheet(sheet)
  long$group <- .classify_group(long$sample,
                                list(Sham       = '^Sham_',
                                     PAB_Mod    = '^PAB_Mod_',
                                     PAB_Severe = '^PAB_Sev'))
  long <- long %>% filter(!is.na(group), state == 'Succ')
  ref <- mean(long$value[long$group == 'Sham'], na.rm = TRUE)
  long$pct <- 100 * (long$value / ref - 1)
  long$display <- factor(paste0('PAB ', long$group, ' (', tp_label, ')'),
                          levels = c(paste0('PAB Sham (', tp_label, ')'),
                                     paste0('PAB PAB_Mod (', tp_label, ')'),
                                     paste0('PAB PAB_Severe (', tp_label, ')')))
  long
}
pab_2wk_block <- .pab_block('PAB_2_Week_Data', '2 wk')
pab_2mo_block <- .pab_block('PAB_8_Week_Data', '2 mo')

# Combine
.combine_cols <- c('display','pct')
cross_df <- bind_rows(
  adult_block[, .combine_cols],
  pab_2wk_block[, .combine_cols],
  pab_2mo_block[, .combine_cols],
  peds_block[, .combine_cols]
)
cross_df$display <- factor(
  cross_df$display,
  levels = c('Adult NF','Adult pRV','Adult RVF',
             'Ped NF','Ped HLHS-fail',
             'PAB Sham (2 wk)','PAB PAB_Mod (2 wk)','PAB PAB_Severe (2 wk)',
             'PAB Sham (2 mo)','PAB PAB_Mod (2 mo)','PAB PAB_Severe (2 mo)')
)

# Cohort coloring
cohort_pal <- c('Adult NF'='#4DAF4A','Adult pRV'='#FF7F00','Adult RVF'='#E41A1C',
                'PAB Sham (2 wk)'='#4DAF4A','PAB PAB_Mod (2 wk)'='#FF7F00','PAB PAB_Severe (2 wk)'='#E41A1C',
                'PAB Sham (2 mo)'='#4DAF4A','PAB PAB_Mod (2 mo)'='#FF7F00','PAB PAB_Severe (2 mo)'='#E41A1C',
                'Ped NF'='#4DAF4A','Ped HLHS-fail'='#E41A1C')

# Within-cohort Mann-Whitney p-values vs reference (NF / Sham / ped-NF)
.compare_to_ref <- function(values_disease, values_ref) {
  if (length(values_disease) < 2 || length(values_ref) < 2) return(NA_real_)
  tryCatch(wilcox.test(values_disease, values_ref)$p.value,
           error = function(e) NA_real_)
}
ref_lookup <- list(
  'Adult pRV'             = list(d = adult_block$value[adult_block$group=='pRV'],
                                  r = adult_block$value[adult_block$group=='NF']),
  'Adult RVF'             = list(d = adult_block$value[adult_block$group=='RVF'],
                                  r = adult_block$value[adult_block$group=='NF']),
  'PAB PAB_Mod (2 wk)'    = list(d = pab_2wk_block$value[pab_2wk_block$group=='PAB_Mod'],
                                  r = pab_2wk_block$value[pab_2wk_block$group=='Sham']),
  'PAB PAB_Severe (2 wk)' = list(d = pab_2wk_block$value[pab_2wk_block$group=='PAB_Severe'],
                                  r = pab_2wk_block$value[pab_2wk_block$group=='Sham']),
  'PAB PAB_Mod (2 mo)'    = list(d = pab_2mo_block$value[pab_2mo_block$group=='PAB_Mod'],
                                  r = pab_2mo_block$value[pab_2mo_block$group=='Sham']),
  'PAB PAB_Severe (2 mo)' = list(d = pab_2mo_block$value[pab_2mo_block$group=='PAB_Severe'],
                                  r = pab_2mo_block$value[pab_2mo_block$group=='Sham']),
  'Ped HLHS-fail'         = list(d = peds_block$value[peds_block$group=='HLHS-failing'],
                                  r = peds_block$value[peds_block$group=='ped-NF'])
)
ref_pvals <- sapply(ref_lookup, function(x) .compare_to_ref(x$d, x$r))
star_lookup <- sapply(ref_pvals, .pstar)

# Compute per-display median pct for star placement
star_df <- data.frame(
  display = names(ref_pvals),
  pct     = sapply(names(ref_pvals), function(d)
              median(cross_df$pct[cross_df$display == d], na.rm = TRUE)),
  star    = star_lookup
)
star_df$display <- factor(star_df$display, levels = levels(cross_df$display))

p_5J <- ggplot(cross_df, aes(x = display, y = pct, fill = display)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey50') +
  geom_boxplot(width = 0.78, outlier.shape = NA, alpha = 0.55) +
  geom_jitter(width = 0.12, size = 2.44, alpha = 0.85, color = 'black') +
  geom_text(data = star_df,
            aes(x = display, y = 30, label = star),
            inherit.aes = FALSE, size = 17, vjust = 0) +
  scale_fill_manual(values = cohort_pal) +
  ylab('Headline coupled OXPHOS\n(% change from cohort control)') +
  xlab(NULL) +
  ggtitle('Cross-cohort respiratory deficit — Adult CI+II / PAB CI+II+FAO / Pediatric CI+II') +
  .respiro_theme() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 31))
save_figure(p_5J, 'Figure_5_panel_J_cross_cohort_summary.pdf',
            width = 20, height = 9)
message('  wrote Figure_5_panel_J_cross_cohort_summary.pdf')


#######################################
####### COMBINED Figure 5 G-J #########
#######################################
## Layout (top→bottom, single column; each "boxplot row" gets equal height):
##   Row 1: G  (4-complex; adult)            — 1 unit
##   Row 2: I  (4-complex; pediatric)        — 1 unit
##   Row 3: H_2wk  (3-state PAB 2-week)      — 1 unit
##   Row 4: H_2mo  (3-state PAB 2-month)     — 1 unit
##   Row 5: J  (cross-cohort summary)        — 1 unit
##
## Each row is the same vertical height so boxplot scale is consistent across
## panels. Total height bumped 25% over the prior 15.75 in composite (→ 19.7 in).

p_5_combined <- p_5G / p_5I / p_5H_2wk / p_5H_2mo / p_5J +
  patchwork::plot_layout(heights = c(1, 1, 1, 1, 1)) +
  patchwork::plot_annotation(
    tag_levels = list(c('G', 'I', 'H', '', 'J'))
  ) &
  theme(plot.tag = element_text(face = 'bold', size = 13))
save_figure(p_5_combined, 'Figure_5_panels_GHIJ_combined.pdf',
            width = 20, height = 46)
message('  wrote Figure_5_panels_GHIJ_combined.pdf  (composite I, G, H_2wk, H_2mo, J — equal rows)')

message('Figure 5 respirometry panels (G–J + combined) complete.')

