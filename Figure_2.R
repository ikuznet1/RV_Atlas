###############################################################################
## Figure 2 (v53 final) — Single-nucleus transcriptomics atlas of 12 RV lineages
##
## Panels (final, derived from new_scripts/Figure_2.png):
##   (A) snRNA-seq UMAP of 61,398 nuclei coloured by 12 major lineage clusters
##       (CM, FB, EC, Endo, SM, Myeloid, NKT, PC, Adipo, Epi, LEC, Neuron)
##   (B) Stacked violin plots of canonical lineage markers
##       PLIN1, RYR2, LEPR, VWF, HAS1, DCN, CCL21, CSF1R, NRXN1, CD2, PDGFRB, MYH11
##   (C) Per-lineage frequency boxplots by NF/pRV/RVF with ANOVA p-values
##   (D) Two-row WGCNA module-score dotplot:
##       (i)  by lineage (CM, EC, FB, Myeloid, PC, SM)
##       (ii) by group   (NF, pRV, RVF)
##       Modules ordered M20, M5, M1, M3, M4, M8, M2, M12, M25, M26, M10,
##                       M28, M14, M11
##
## Outputs (./output/Figure_2/v52_figures/):
##   Figure_2_panel_A.pdf, Figure_2_panel_B_violins.pdf,
##   Figure_2_panel_C_freq_box.pdf,
##   Figure_2_panel_D_modules_by_lineage.pdf,
##   Figure_2_panel_D_modules_by_group.pdf
##   Supplementary QC: Figure_2_supp_panel_A_patient.pdf
###############################################################################

source('./helper_scripts/_shared_helpers.R')

## Per-figure output directory (introduced for consistent output paths)
V52_FIG_DIR <- './output/Figure_2'
dir.create(V52_FIG_DIR, showWarnings = FALSE, recursive = TRUE)


## Suppress R's default Rplots.pdf in cwd when Rscript hits a plot call
## that's outside an explicit pdf() ... dev.off() envelope.
pdf(NULL)
COMP_W <- 5
COMP_H <- 5
PS <- pub_scales(COMP_W)

## -- Cartoon asset extraction for Figure 2 (run-once, commented) ------------
## Inspection of ~/Downloads/hdWGCNA_TOM/Manuscripts/Figure_2.pdf shows no
## hand-drawn cartoons (all panels are generated data: UMAP, violin, boxplot,
## DotPlots). No insert_asset() calls are needed for Figure 2.
# library(magick)
# fig2 <- image_read_pdf("~/Downloads/hdWGCNA_TOM/Manuscripts/Figure_2.pdf", density=300)
# ## No hand-drawn cartoon was identified in the original.

library(Seurat)
library(hdWGCNA)
library(ggeasy)


source('./helper_scripts/spatial_functions.R')



M1<-readRDS('./dependencies/shared/RV_data.rds')

#######################################
#############  FIGURE 2B  #############
#######################################

feats <- c("PLIN1","RYR2","LEPR","VWF","HAS1","DCN","CCL21","CSF1R","NRXN1","CD2","PDGFRB","MYH11")
p_2B <- VlnPlot(M1, features = feats, ncol=1, pt.size=F, group.by="Names")

for(i in 1:12) {
   p_2B[[i]] <- p_2B[[i]] + NoLegend() +
                easy_remove_axes(which="y", what = c("ticks","text","line")) +
                ggtitle("") + ylab(feats[i])
   if(i < 12) p_2B[[i]] <- p_2B[[i]] + easy_remove_axes(which="x")
}

save_figure(p_2B, 'Figure_2_panel_B_violins.pdf', width = 4, height = 20)


#######################################
#############  FIGURE 2C  #############
#######################################
M1 <- SetIdent(M1, value = "Names")



all_percent_cell <- as.data.frame(table(M1@active.ident)/length(M1@active.ident)*100)

NF_percent_cell <- cbind(as.data.frame(table(subset(M1,group=="NF")@active.ident)/length(subset(M1,group=="NF")@active.ident)*100),type = "NF")
NF_percent_cell$sum <- (rev(cumsum(rev(NF_percent_cell$Freq))) - NF_percent_cell$Freq/2)/100
NF_percent_cell$Freq <- NF_percent_cell$Freq/100


pRV_percent_cell <- cbind(as.data.frame(table(subset(M1,group=="pRV")@active.ident)/length(subset(M1,group=="pRV")@active.ident)*100),type = "pRV")
pRV_percent_cell$sum <- (rev(cumsum(rev(pRV_percent_cell$Freq))) - pRV_percent_cell$Freq/2)/100
pRV_percent_cell$Freq <- pRV_percent_cell$Freq/100


RVF_percent_cell <- cbind(as.data.frame(table(subset(M1,group=="RVF")@active.ident)/length(subset(M1,group=="RVF")@active.ident)*100),type = "RVF")
RVF_percent_cell$sum <- (rev(cumsum(rev(RVF_percent_cell$Freq))) - RVF_percent_cell$Freq/2)/100
RVF_percent_cell$Freq <- RVF_percent_cell$Freq/100


percent_cell_df <- rbind(NF_percent_cell,pRV_percent_cell,RVF_percent_cell)

## Stacked-bar lineage composition (Supp QC; final Panel C uses the per-lineage
## boxplot view below).
p_2C_stack <- ggplot(percent_cell_df, aes(fill=Var1, y=Freq, x=type, label=round(sum,1))) +
  geom_bar(position="stack", stat="identity", width=0.6, linewidth = PS$geom_lw) +
  scale_fill_lineage(name = "Cell type") +
  theme_classic() +
  xlab("Disease State") + ylab("Frequency") +
  labs(fill="Cell type", color='black') +
  theme(text = element_text(size=20),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        legend.text = element_text(color="black"),
        legend.key.size = PS$legend_key) +
  scale_y_touch() +
  geom_label_repel(aes(type, sum, label = scales::percent(round(Freq,2))),
                   fill = NA, nudge_x = 0.5, direction = "y",
                   size = PS$text_mm, family = FONT_FAMILY)
save_figure(p_2C_stack, 'Figure_2_supp_panel_C_stacked.pdf', width = 5, height = 5)


#Prevalence comparisons
cells <- table(M1@active.ident,M1@meta.data$patient)
cells <- cells[c('CM','FB','Myeloid','EC','Adipo','PC','SM'),]
cells <- sweep(cells,2,colSums(cells),'/')
cells <- data.frame(cells)


#ggboxplot(cells[77:1,],x="Var2",y="Freq",fill="Var2",group="Var2")+theme_classic() + theme(axis.text.x=element_text(size=16),axis.text.y=element_text(size=16),axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),legend.title=element_text(size=16),legend.text=element_text(size=16),text=element_text(color='black'),axis.text=element_text(color='black')) + labs(color='Group',x="Disease",y='Frequency') + facet_wrap(~Var1,ncol=7)


cells$group = rep(c("RVF","pRV","RVF","NF","pRV","pRV","RVF","NF","NF","pRV","NF"),each=7)


#my_comparisons <- list( c("NF", "pRV"),c("pRV", "RVF"),c("NF", "RVF"))

library(ggpubr)
p_2C <- ggboxplot(cells[77:1,], x="group", y="Freq", fill="group", group="group",
                  size = PS$geom_lw) +
  scale_fill_disease(name = "Group") +
  theme_classic() +
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        legend.title = element_text(size=16),
        legend.text  = element_text(size=16),
        text = element_text(color='black'),
        axis.text = element_text(color='black'),
        legend.key.size = PS$legend_key) +
  labs(color='Group', x="Disease", y='Frequency') +
  facet_wrap(~ Var1, ncol = 7) +
  stat_compare_means(aes(group=group), method="anova",
                     size = PS$text_mm, family = FONT_FAMILY,
                     symnum.args = pub_signif_args)
save_figure(p_2C, 'Figure_2_panel_C_freq_box.pdf', width = 12.5, height = 5)






#######################################
#############  FIGURE 2D  #############
#######################################
mapping <- labels2colors(1:100)
.cache_MEs <- './dependencies/shared/scWGCNA_bulk2sn_MEs.rds'

if (file.exists(.cache_MEs)) {
  message('Loading cached MEs - skipping 20 GB scWGCNA load')
  MEs <- readRDS(.cache_MEs)
} else {
  message('No MEs cache; loading scWGCNA_bulk2sn_projection.rds (one-time)')
  seurat_ref <- readRDS('./dependencies/shared/scWGCNA_bulk2sn_projection.rds')
  seurat_ref <- SetActiveWGCNA(seurat_ref, 'bulk2sn')
  dir.create(dirname(.cache_MEs), showWarnings = FALSE, recursive = TRUE)
  MEs <- GetMEs(seurat_ref, harmonized = TRUE)
  saveRDS(MEs, .cache_MEs)
  rm(seurat_ref); gc()
}
mods <- colnames(MEs); mods <- mods[mods != 'grey']
mods_num <- paste0('M',match(mods,mapping))
prv_vs_rv_signif <- c('M3','M4','M5','M10','M11','M12','M14')
all_signif <- c('M1','M2','M3','M4','M5','M8','M10','M11','M12','M14','M20','M25','M26','M28')


colnames(MEs) <- paste0('M', match(colnames(MEs), mapping))

consensus_modules <- read.csv("./dependencies/shared/bulk_heart_modules.csv")
consensus_modules <- consensus_modules[,1:3]


seurat_ref <- readRDS('./dependencies/shared/RV_data.rds')

# Gene-presence filter deferred until RV_data is loaded:
# AddModuleScore (called below) only scores genes in rownames(seurat_ref),
# so filtering here matches the gene universe actually used downstream.
consensus_modules <- subset(consensus_modules, gene_name %in% rownames(seurat_ref))
consensus_modules <- consensus_modules[match(unique(consensus_modules$gene_name), consensus_modules$gene_name),]

library(dplyr)
score_calc <- consensus_modules %>% group_by(module) %>% group_split()
module_colors <- unique(unlist(lapply(score_calc,'[[','module')))
module_colors <- paste0('M',match(module_colors,mapping))
seurat_ref <- SetIdent(seurat_ref, value = "Names")
seurat_ref@meta.data <- cbind(seurat_ref@meta.data, MEs)



.cache_module_scores <- './output/Figure_2/fig2_module_scores_cache.rds'
if (file.exists(.cache_module_scores)) {
  message('Loading cached module scores...')
  module_score_df <- readRDS(.cache_module_scores)
  seurat_ref@meta.data <- cbind(seurat_ref@meta.data, module_score_df)
} else {
  seurat_ref <- AddModuleScore(seurat_ref,lapply(score_calc,'[[','gene_name'),name="module_score")
  cols_current <- colnames(seurat_ref@meta.data)
  cols_current[startsWith(colnames(seurat_ref@meta.data),'module_score')] <- paste0('module_',module_colors)
  colnames(seurat_ref@meta.data) <- cols_current
  module_score_df <- seurat_ref@meta.data[, paste0('module_', module_colors), drop = FALSE]
  saveRDS(module_score_df, .cache_module_scores)
}


## Module ordering used in both Panel D rows.
mods_D <- paste0('module_',
                 c('M20','M5','M1','M3','M4','M8','M2','M12',
                   'M25','M26','M10','M28','M14','M11'))

## Panel D (i): module score by lineage (CM, EC, FB, Myeloid, PC, SM).
p_2Di <- DotPlot(seurat_ref, mods_D, group.by = 'Names',
                 dot.min = 0, col.min = 0, col.max = 2,
                 idents = c('CM','EC','FB','Myeloid','PC','SM')) +
  RotatedAxis() + ylab('') + xlab('') +
  scale_color_gradient2(high = 'red', mid = 'grey95', low = 'blue') +
  theme(panel.border = element_rect(size = 1, fill = NA, color = 'black'),
        axis.line.x  = element_blank(),
        axis.line.y  = element_blank())
save_figure(p_2Di, 'Figure_2_panel_D_modules_by_lineage.pdf', width = 7, height = 4)

## Panel D (ii): module score by disease group (NF / pRV / RVF).
p_2Dii <- DotPlot(seurat_ref, mods_D, group.by = 'group',
                  dot.min = 0, col.min = 0, col.max = 2) +
  RotatedAxis() + ylab('') + xlab('') +
  scale_color_gradient2(high = 'red', mid = 'grey95', low = 'blue') +
  theme(panel.border = element_rect(size = 1, fill = NA, color = 'black'),
        axis.line.x  = element_blank(),
        axis.line.y  = element_blank())
save_figure(p_2Dii, 'Figure_2_panel_D_modules_by_group.pdf', width = 6.75, height = 2.5)

###############################################################################
#############  FIGURE 2A — v52: pixel-perfect UMAP (matches original)  ########
###############################################################################
## Uses PlotEmbedding (hdWGCNA) + umap_theme() to match the original style:
##   - simple text labels at cluster centroids (no white-box geom_label)
##   - UMAP1/UMAP2 axis arrows from umap_theme()
##   - "n = 61,398 nuclei" annotation (bottom-right, grey)
##   - no legend, no plot title
## Output: ./output/Figure_2/v52_figures/Figure_2_panel_A.pdf

Idents(seurat_ref) <- 'Names'
n_nuclei_A <- formatC(ncol(seurat_ref), format = 'd', big.mark = ',')

p_2A <- PlotEmbedding(seurat_ref,
                      group.by    = 'Names',
                      point_size  = PS$umap_pt,
                      plot_under  = TRUE,
                      raster_dpi  = 400,
                      raster_scale = 0.5,
                      plot_theme  = umap_theme() + NoLegend()) +
  annotate('text',
           x = Inf, y = -Inf,
           label    = paste0('n = ', n_nuclei_A, ' nuclei'),
           hjust    = 1.05, vjust = -0.4,
           size     = PS$text_mm,
           family   = FONT_FAMILY,
           colour   = 'grey40')

save_figure(p_2A, 'Figure_2_panel_A.pdf', width = 5, height = 5)

###############################################################################
##  Supp QC — UMAP coloured by patient (no patient-driven artefacts)
###############################################################################

n_nuclei_Apat <- formatC(ncol(seurat_ref), format = 'd', big.mark = ',')
n_pat <- length(unique(seurat_ref$patient))

p_2A_patient <- PlotEmbedding(seurat_ref,
                              group.by    = 'patient',
                              point_size  = PS$umap_pt,
                              plot_under  = TRUE,
                              raster_dpi  = 400,
                              raster_scale = 0.5,
                              label = FALSE,
                              plot_theme  = umap_theme()) +
  annotate('text',
           x = Inf, y = -Inf,
           label    = paste0('n = ', n_nuclei_Apat, ' nuclei, ', n_pat, ' patients'),
           hjust    = 1.05, vjust = -0.4,
           size     = PS$text_mm,
           family   = FONT_FAMILY,
           colour   = 'grey40') +
  guides(colour = guide_legend(override.aes = list(size = 3),
                               title = 'Patient',
                               ncol = 1)) +
  theme(legend.position = 'right',
        legend.key.size = PS$legend_key,
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9))

save_figure(p_2A_patient, 'Figure_2_supp_panel_A_patient.pdf', width = 6.5, height = 5)

message('Figure 2 (v53) per-panel PDFs written to ', V52_FIG_DIR)
