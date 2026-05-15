###############################################################################
## Figure 7 (v53 final) — Endothelial cells: Phase 1 antioxidant suppression
## with subtype-specific Phase 2 reactivation
##
## Panels (final, per v53 manuscript legend):
##   (A) EC UMAP, 5 subtypes (Arterial, Capillary, Venous, Lymph, Endocardial)
##   (B) Marker dotplot (curated canonical markers, diagonal-ordered)
##   (C) Subtype frequency boxplots NF/pRV/RVF with ANOVA
##   (D) hdWGCNA hub gene UMAPs + module identification
##   (E) GO enrichment for EC-specific modules (ecM1, ecM4, ecM7)
##   (F) EC module dotplots by subtype and by disease group
##   (G) FeaturePlot of ecM1, ecM4, ecM7 on UMAP
##   (H) Module volcano plots highlighting Phase 2 reactivation genes
##   (I) Xenium spatial: ANGPT2 and SNAI1 induction in RVF endothelium
##   (J) NRG1 expression in Endocardial EC across NF/pRV/RVF (violin)
##   (K) Angiogenesis biphasic score (combined positive + negative regulators)
##       across NF/pRV/RVF
##   (L) MECOM expression across EC subtypes × disease group
##
## Outputs (./output/Figure_7/v52_figures/):
##   Figure_7_panel_{A,B,C,D,E,F,G,H,I,J,K,L}.pdf
###############################################################################

source('./helper_scripts/_shared_helpers.R')

## Per-figure output directory (introduced for consistent output paths)
V52_FIG_DIR <- './output/Figure_7'
dir.create(V52_FIG_DIR, showWarnings = FALSE, recursive = TRUE)


## Suppress R's default Rplots.pdf in cwd when Rscript hits a plot call
## that's outside an explicit pdf() ... dev.off() envelope.
pdf(NULL)
dir.create('./output', recursive = TRUE, showWarnings = FALSE)

## Composite figure dimensions (inches) — used by theme_v52() scaling
COMP_W <- 14
COMP_H <- 15

## Publication scales (geom widths, point sizes, text sizes) scaled to COMP_W.
PS <- pub_scales(COMP_W)

## ── Panel A cartoon (Myeloid lineage schematic at top-left of Figure_6) ──
## To regenerate cropped cartoon from the published Figure_6.pdf, run in R:
##   library(magick)
##   img <- image_read_pdf('~/Downloads/hdWGCNA_TOM/Manuscripts/Figure_6.pdf',
##                         pages = 1, density = 300)
##   ## Approx. geometry for the Myeloid cartoon at top-left (17 cm ≈ 2008 px):
##   cart <- image_crop(img, '280x150+40+20')
##   image_write(cart, './new_scripts/assets/Figure_6_panel_A_myeloid_cartoon.png',
##               format = 'png', density = 300)
## Then compose with `insert_asset('Figure_7_panel_A_myeloid_cartoon.png')` in the layout.

library(Seurat)
library(hdWGCNA)
library(ggeasy)
library(harmony)


source('./helper_scripts/spatial_functions.R')



#######################################
#############  FIGURE 7A  #############
#######################################
M1 <- readRDS(file = "./dependencies/shared/EC_hdWGCNA_by_celltype.rds")


p_7A <- PlotEmbedding(M1, group.by = 'Names',
                      point_size = PS$umap_pt, plot_under = TRUE,
                      plot_theme = umap_theme() + NoLegend(),
                      raster_dpi = 400, raster_scale = 0.5)
save_figure(p_7A, 'Figure_7_panel_A.pdf', width = 5, height = 5)

## Supp QC — same UMAP coloured by disease group.
p_7A_group <- PlotEmbedding(M1, group.by = 'group',
                            point_size = PS$umap_pt, plot_under = TRUE,
                            plot_theme = umap_theme(),
                            raster_dpi = 400, raster_scale = 0.5)
save_figure(p_7A_group, 'Figure_7_supp_panel_A_group.pdf', width = 5, height = 5)

#######################################
#############  FIGURE 7B  #############
#######################################
## EC subtype marker dotplot — canonical lineage markers, ordered along the
## main diagonal (y-axis subtype order matches x-axis marker block order).
##
## Marker references:
##   Arterial:    GJA5 (Coppen 1999 Circ Res); DKK2 (Kalucka 2020 Cell);
##                EFNB2 (Wang 1998 Cell)
##   Capillary:   RGCC, CA4, CD36 (Litvinukova 2020 Nature; Kalucka 2020)
##   Venous:      NR2F2/COUP-TFII (You 2005 Nature); ACKR1/DARC (Thiriot 2017);
##                EPHB4 (Wang 1998 Cell)
##   Lymphatic:   PROX1 (Wigle 1999 Cell); PDPN (Schacht 2003)
##   Endocardial: NPR3 (Tian 2014 Cell Reports; Litvinukova 2020);
##                NRG3 (Tucker 2020 Circulation)
ec_subtype_levels <- c('Arterial','Capillary','Venous','Lymph','Endocardial')

ec_markers <- c(
  # Arterial
  'GJA5','DKK2','EFNB2',
  # Capillary
  'RGCC','CA4','CD36',
  # Venous
  'NR2F2','ACKR1','EPHB4',
  # Lymphatic
  'PROX1','PDPN',
  # Endocardial
  'NPR3','NRG3'
)
ec_markers <- intersect(ec_markers, rownames(M1))

M1$Names <- factor(M1$Names, levels = ec_subtype_levels)
Idents(M1) <- 'Names'

p_7B <- DotPlot(M1, features = ec_markers,
                col.min = 0, col.max = 2, dot.scale = 6) +
  RotatedAxis() +
  xlab(NULL) + ylab('EC subtype') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = 'italic'),
        panel.border = element_rect(colour = 'black', fill = NA, linewidth = PS$geom_lw),
        axis.line.x  = element_blank(),
        axis.line.y  = element_blank())

save_figure(p_7B, 'Figure_7_panel_B.pdf', width = 8, height = 2.7)


#######################################
#############  FIGURE 7C  #############
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

## Supp QC stacked-bar lineage composition (final Panel C uses per-subtype boxes).
p_7C_stack <- ggplot(percent_cell_df, aes(fill=Var1, y=Freq, x=type, label=round(sum,1))) +
  geom_bar(position="stack", stat="identity", width=0.6, linewidth = PS$geom_lw) +
  theme_v52(COMP_W) +
  scale_fill_lineage(name = "Cell type") +
  xlab("Disease State") + ylab("Frequency") + labs(fill="Cell type", color='black') +
  scale_y_touch() +
  geom_label_repel(aes(type, sum, label = scales::percent(round(Freq,2))),
                   fill = NA, nudge_x = 0.5, direction = "y",
                   size = PS$text_mm, family = FONT_FAMILY)
save_figure(p_7C_stack, 'Figure_7_supp_panel_C_stacked.pdf', width = 5, height = 5)


#Prevalence comparisons
cells <- table(M1@active.ident,M1@meta.data$patient)
cells <- sweep(cells,2,colSums(cells),'/')
cells <- data.frame(cells)


#ggboxplot(cells[77:1,],x="Var2",y="Freq",fill="Var2",group="Var2")+theme_classic() + theme(axis.text.x=element_text(size=16),axis.text.y=element_text(size=16),axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),legend.title=element_text(size=16),legend.text=element_text(size=16),text=element_text(color='black'),axis.text=element_text(color='black')) + labs(color='Group',x="Disease",y='Frequency') + facet_wrap(~Var1,ncol=7)


cells$group = rep(c("RVF","pRV","RVF","NF","pRV","pRV","RVF","NF","NF","pRV","NF"),each=5)


#my_comparisons <- list( c("NF", "pRV"),c("pRV", "RVF"),c("NF", "RVF"))

library(ggpubr)
p_7C <- ggboxplot(cells[55:1, ], x="group", y="Freq", fill="group", group="group",
                  size = PS$geom_lw) +
  theme_v52(COMP_W) +
  scale_fill_disease() +
  labs(color='Group', x="Disease", y='Frequency') +
  facet_wrap(~ Var1, ncol = 7) +
  stat_compare_means(aes(group = group), method = "anova",
                     size = PS$text_mm, symnum.args = pub_signif_args)
save_figure(p_7C, 'Figure_7_panel_C.pdf', width = 12.5, height = 4)


#######################################
#############  FIGURE 7D  #############
#######################################
DefaultAssay(M1) <- "RNA"

.cache_wgcna_fit <- './output/Figure_7/fig7_wgcna_fit_cache.rds'
if (file.exists(.cache_wgcna_fit)) {
  message('Loading cached hdWGCNA fit (SetupForWGCNA through ModuleConnectivity)...')
  M1 <- readRDS(.cache_wgcna_fit)
} else {
  M1 <- SetupForWGCNA(
    M1,
    gene_select = "fraction", # the gene selection approach
    fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
    wgcna_name = "ec" # the name of the hdWGCNA experiment
  )

  M1 <- MetacellsByGroups(seurat_obj = M1,
    group.by = c("Names"), # specify the columns in seurat_obj@meta.data to group by
    reduction = 'harmony', # select the dimensionality reduction to perform KNN on
    k = 25, # nearest-neighbors parameter
    max_shared = 10, # maximum number of shared cells between two metacells
    ident.group = 'Names' # set the Idents of the metacell seurat object
  )

  M1 <- NormalizeMetacells(M1)

  M1 <- SetDatExpr(
    M1,
    group_name = c("Capillary","Arterial","Venous"), # the name of the group of interest in the group.by column
    group.by='Names', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
    assay = 'RNA', # using RNA assay
    slot = 'data' # using normalized data
  )

  M1 <- TestSoftPowers(
    M1,
    networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
  )

  ## ConstructNetwork honours tom_dir, but the underlying WGCNA call still
  ## leaks `./TOM/<name>_TOM.rda` to getwd(). Run inside V52_FIG_DIR so the
  ## artifact stays under ./output/Figure_7/.
  .oldwd <- getwd()
  setwd(V52_FIG_DIR)
  on.exit(setwd(.oldwd), add = TRUE)
  M1 <- ConstructNetwork(overwrite_tom = TRUE,
    M1,
    tom_dir = 'hdWGCNA_TOM',
    tom_name = 'ec_net' # name of the topoligical overlap matrix written to disk
  )
  setwd(.oldwd)

  M1 <- ScaleData(M1)
  ## NOTE: hdWGCNA 0.4.11 + harmony 2.0.3 are incompatible — ModuleEigengenes
  ## internally passes `assay.use=...` to RunHarmony, which harmony 2.x rejects.
  ## M1 was loaded already-harmonised (EC_hdWGCNA_by_celltype.rds), so dropping
  ## group.by.vars here skips the redundant internal harmony call.
  M1 <- ModuleEigengenes(M1)

  M1 <- ModuleConnectivity(
    M1,
    group.by = 'Names', group_name =  c("Capillary","Arterial","Venous")
  )

  saveRDS(M1, .cache_wgcna_fit)
}

library(patchwork)
plot_list <- PlotSoftPowers(M1)
wrap_plots(plot_list, ncol=2)

MEs <- GetMEs(M1,harmonized=TRUE)

modules <- GetModules(M1)
mods <- levels(modules$module); mods <- mods[mods != 'grey']

#Blue is CM - M2
#Yellow module EC specific - M4
#Green is endocardial - M5
#Brown is SM - M3
#Turqouise is EC specific - M1
#Red is mural cells and fibroblast - M6
#Black is ECs - M7

# add hMEs to Seurat meta-data:
M1@meta.data <- cbind(M1@meta.data, MEs)

mapping <- labels2colors(1:100)
mods_num <- paste0('M',match(mods,mapping))

colNames <- colnames(M1@meta.data)
colNames[match(mods,colNames)] <- mods_num
colnames(M1@meta.data) <- colNames

# compute the module UMAPs
.cache_module_umap <- './output/Figure_7/fig7_module_umap_cache.rds'
if (file.exists(.cache_module_umap)) {
  message('Loading cached RunModuleUMAP (M1)...')
  M1 <- readRDS(.cache_module_umap)
} else {
  M1 <- RunModuleUMAP(
    M1,
    n_hubs = 5,
    n_neighbors=10,
    min_dist=0.3,
    spread=2,
    target_weight=0.1,
    supervised=TRUE
  )
  saveRDS(M1, .cache_module_umap)
}

# get the hub gene UMAP table from the seurat object
umap_df <- GetModuleUMAP(M1)

# plot with ggplot
plot_df <- umap_df

# compute coordinates for cluster labels
centroid_df <- data.frame()
for(cur_cluster in unique(plot_df[['module']])){
    cur_meta <- plot_df[plot_df[['module']] == cur_cluster,]
    df <- data.frame(
      cluster = cur_cluster,
      UMAP1 = mean(cur_meta$UMAP1),
      UMAP2 = mean(cur_meta$UMAP2)
    )
  centroid_df <- rbind(centroid_df, df)
}

# plot with ggplot
p <- ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
    ggrastr::rasterise(geom_point(
     color=umap_df$color,
     size=umap_df$kME*2
   ), dpi=500, scale=0.5) +
    umap_theme() +
    theme(
      plot.margin = margin(0,0,0,0),
      plot.title = element_text(hjust=0.5)
    ) + ggtitle('Bulk RNASeq') +
    ggrepel::geom_text_repel(data = centroid_df, label=centroid_df$cluster, color='black', max.overlaps=Inf, size = PS$text_mm, family = FONT_FAMILY)


hub_genes <- GetHubGenes(M1, 3)

  # add annotation
anno_genes <- hub_genes$gene_name
plot_df$anno <- ifelse(plot_df$gene %in% anno_genes, umap_df$gene, '')

plot_df_anno <- subset(plot_df, anno != '')
p <-  plot_df %>%
    ggplot(aes(x=UMAP1, y=UMAP2, color=module)) +
    ggrastr::rasterise(
      geom_point(
        inherit.aes=FALSE,
        data=plot_df,
        aes(x=UMAP1, y=UMAP2, color=module),
        color=plot_df$color,
        size=plot_df$kME*2,
      ), dpi=500, dpi_scale=0.5) +
    geom_point(
      inherit.aes = FALSE,
      data = plot_df_anno,
      shape=21, color='black',
      fill=plot_df_anno$color,
      size=plot_df_anno$kME*2,
      aes(x=UMAP1, y=UMAP2, fill=module)
    ) +
    # add labels
    ggrepel::geom_text_repel(data = centroid_df, label=paste0('M',match(centroid_df$cluster,mapping)), color='black', max.overlaps=Inf, size = PS$text_mm, family = FONT_FAMILY, fontface='bold') +
    geom_text_repel(label=plot_df$anno, max.overlaps=Inf, color='black', fontface='italic', size = PS$text_mm, family = FONT_FAMILY) +
    umap_theme() + NoLegend() +
    coord_equal() +
    theme(
      plot.margin = margin(0,0,0,0)
    )

pdf('./output/Figure_7/EC_hdWGCNA.pdf',width=6,height=6)
print(p)
dev.off()


library(enrichR)
dbs <- c('GO_Biological_Process_2023','GO_Cellular_Component_2023','GO_Molecular_Function_2023')

.cache_enrichr <- './output/Figure_7/fig7_enrichr_cache.rds'
if (file.exists(.cache_enrichr)) {
  message('Loading cached RunEnrichr (M1)...')
  M1 <- readRDS(.cache_enrichr)
} else {
  # perform enrichment tests
  M1 <- RunEnrichr(
    M1,
    dbs=dbs, # character vector of enrichr databases to test
    max_genes = 10000 # number of genes per module to test. use max_genes = Inf to choose all genes!
  )
  saveRDS(M1, .cache_enrichr)
}

# retrieve the output table
enrich_df <- GetEnrichrTable(M1)

EnrichrBarPlot(
  M1,
  outdir = "./output/Figure_7/scEC_subclust_enrichr_plot", 
  n_terms = 5,
  plot_size = c(5,4), # width, height of the output .pdfs
  logscale=TRUE # do you want to show the enrichment as a log scale?
)

#######################################
#############  FIGURE 7E  #############
#######################################


modules <- GetModules(M1)
color_df <- modules %>% subset(module!='grey') %>%
  select(c(module, color)) %>% distinct %>%
  dplyr::rename(group = module, colour = color)
mods <- levels(modules$module)
mods <- mods[mods!='grey']

# helper function to wrap text
wrapText <- function(x, len) {
    sapply(x, function(y) paste(strwrap(y, len), collapse = "\n"), USE.NAMES = FALSE)
}

combined_output <- GetEnrichrTable(M1)
selected_terms <- subset(combined_output,db=="GO_Biological_Process_2023")

# subset selected terms
selected_terms <- subset(selected_terms, P.value < 0.05)
idx_top_1 <- match(unique(selected_terms$module),selected_terms$module)
idx_top_3 <- sort(c(idx_top_1,idx_top_1+1,idx_top_1+2))

selected_terms<-selected_terms[idx_top_3,]

order_mod <- c(1,4,5,7,2,3,6)

selected_terms$group <- factor(
  as.character(selected_terms$module),
  levels = mapping[order_mod]
)


# set max pval
quantile(-log(selected_terms$P.value), 0.95)
max_p <- 10

selected_terms$logp <- -log(selected_terms$P.value)
selected_terms$logp <- ifelse(selected_terms$logp > max_p, max_p, selected_terms$logp)

# remove Reactome Term ID
library(stringr)
selected_terms$Term <- str_replace(selected_terms$Term, " \\s*\\([^\\)]+\\)", "")

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


color_df <- modules %>% subset(module!='grey') %>%
  select(c(module, color)) %>% distinct %>%
  dplyr::rename(group = module, colour = color)
color_df <- subset(color_df, colour %in% mods)


color_df$group<-paste0('M',match(color_df$group,mapping))

color_df$group <- factor(
  as.character(color_df$group),
  levels = paste0('M',order_mod)
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


pdf(paste0('./output/Figure_7/', 'EC_GO_terms.pdf'), width=8, height=7)
print(p / colorbar)
dev.off()




#######################################
#############  FIGURE 7F  #############
#######################################


pdf('./output/Figure_7/EC_celltype_mods.pdf',width=6,height=3)
p <- DotPlot(M1, features=sort(mods_num), group.by = 'Names',col.min=0)
p
dev.off()


pdf('./output/Figure_7/EC_group_mods.pdf',width=6,height=2)
p <- DotPlot(M1, features=sort(mods_num), group.by = 'group',col.min=0)
p
dev.off()

up_RVF = c('M1','M4','M5','M7')
down_RVF = c('M2','M3','M6')


#######################################
#############  FIGURE 7G  #############
#######################################

pdf('./output/Figure_7/EC_feature.pdf',width=4,height=6)
FeaturePlot(M1,c('M1','M4','M7'),label=T,min.cutoff=0,ncol=1)
dev.off()

#######################################$
############  FIGURE 7H - I ############
#######################################$
M1 <- SetIdent(M1, value = "group")

modules <- GetModules(M1)
modules$module <- match(modules$module,mapping)
.cache_find_markers <- './output/Figure_7/fig7_find_markers_cache.rds'
if (file.exists(.cache_find_markers)) {
  message('Loading cached FindMarkers (combined_set per module)...')
  combined_set <- readRDS(.cache_find_markers)
} else {
  combined_set <- data.frame()
  mods_idx <- c(1,4,5,7,2,3,6)
  for (i in mods_idx){
    key_genes <- subset(modules,module %in% c(i))$gene_name
    key_genes <- key_genes[key_genes %in% rownames(M1)]
    gene_set <- FindMarkers(M1, ident.1 = "RVF", ident.2 = "NF",features=key_genes)
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
  saveRDS(combined_set, .cache_find_markers)
}


library(EnhancedVolcano)

combined_set <- combined_set[!grepl('MT-',rownames(combined_set)),]
keyvals <- combined_set$color
names(keyvals)  <- combined_set$module



M1_subset <- subset(combined_set,color=='turquoise')
M1_subset_up <- subset(M1_subset,avg_log2FC>0)
M1_subset_down <- subset(M1_subset,avg_log2FC<0)
keyvals <- M1_subset$color
names(keyvals)  <- M1_subset$module


pdf(paste0('./output/Figure_7/', 'RV_EC_module_volcano_aM1.pdf'), width=8, height=6)

print(EnhancedVolcano(M1_subset,lab=rownames(M1_subset),
  x='avg_log2FC',y='p_val_adj',
  FCcutoff = 0.1,
  colCustom=keyvals, labFace = "italic") + coord_flip())
dev.off()




cat(rownames(M1_subset_up),sep='\n')
cat(rownames(M1_subset_down),sep='\n')


M4_subset <- subset(combined_set,module=='M4')
M4_subset_up <- subset(M4_subset,avg_log2FC>0)
M4_subset_down <- subset(M4_subset,avg_log2FC<0)


cat(rownames(M4_subset_up),sep='\n')
cat(rownames(M4_subset_down),sep='\n')

#saveRDS(M1,'./output/Figure_7/EC_hdWGCNA_by_celltype.rds')


DefaultAssay(M1) <- "SCT"

#VlnPlot(subset(M1,Names=='Arterial'),'MECOM',group.by='group')
#VlnPlot(M1,'VEGFC',group.by='group')
#VlnPlot(M1,'VASH1',group.by='group')
#VlnPlot(M1,'FLT1',group.by='group')
#VlnPlot(M1,'SP100',group.by='group')
#VlnPlot(M1,'SEMA6A',group.by='group')
#VlnPlot(M1,'PPP1R16B',group.by='group')
#VlnPlot(M1,'PLXND1',group.by='group')


#Get all terms for + regulation of angiogenesis from GO
library(biomaRt)
.cache_angio_scores <- './output/Figure_7/fig7_angio_module_scores_cache.rds'
if (file.exists(.cache_angio_scores)) {
  message('Loading cached biomaRt + AddModuleScore (angiogenesis scores)...')
  angio_scores_df <- readRDS(.cache_angio_scores)
  M1@meta.data[, colnames(angio_scores_df)] <- angio_scores_df
} else {
  ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl") #uses human ensembl annotations
  #gets gene symbol, transcript_id and go_id for all genes annotated with GO:0007507
  gene.data <- getBM(attributes=c('hgnc_symbol', 'ensembl_transcript_id', 'go_id'),
                     filters = 'go', values = 'GO:0045765', mart = ensembl)
  angio_genes <- unique(gene.data$hgnc_symbol)
  angio_genes <- intersect(angio_genes,rownames(M1))
  M1 <- AddModuleScore(M1,list(angio_genes),name="AngioGenes")

  gene.data <- getBM(attributes=c('hgnc_symbol', 'ensembl_transcript_id', 'go_id'),
                     filters = 'go', values = 'GO:0045766', mart = ensembl)
  angio_genes <- unique(gene.data$hgnc_symbol)
  angio_genes <- intersect(angio_genes,rownames(M1))
  M1 <- AddModuleScore(M1,list(angio_genes),name="PositiveAngioGenes")

  gene.data <- getBM(attributes=c('hgnc_symbol', 'ensembl_transcript_id', 'go_id'),
                     filters = 'go', values = 'GO:0016525', mart = ensembl)
  angio_genes <- unique(gene.data$hgnc_symbol)
  angio_genes <- intersect(angio_genes,rownames(M1))
  M1 <- AddModuleScore(M1,list(angio_genes),name="NegativeAngioGenes")

  gene.data <- getBM(attributes=c('hgnc_symbol', 'ensembl_transcript_id', 'go_id'),
                     filters = 'go', values = 'GO:0060978', mart = ensembl)
  angio_genes <- unique(gene.data$hgnc_symbol)
  angio_genes <- intersect(angio_genes,rownames(M1))
  M1 <- AddModuleScore(M1,list(angio_genes),name="CoronaryAngioGenes")

  gene.data <- getBM(attributes=c('hgnc_symbol', 'ensembl_transcript_id', 'go_id'),
                     filters = 'go', values = 'GO:0001525', mart = ensembl)
  angio_genes <- unique(gene.data$hgnc_symbol)
  angio_genes <- intersect(angio_genes,rownames(M1))
  M1 <- AddModuleScore(M1,list(angio_genes),name="AllAngioGenes")

  angio_score_cols <- c('AngioGenes1','PositiveAngioGenes1','NegativeAngioGenes1',
                        'CoronaryAngioGenes1','AllAngioGenes1')
  angio_scores_df <- M1@meta.data[, angio_score_cols, drop=FALSE]
  saveRDS(angio_scores_df, .cache_angio_scores)
}


###############################################################################
##  Panels J + K — Combined Phase-2 subtype-specific violins:
##              Panel J = NRG1 (Endocardial EC); Panel K = MECOM (Arterial EC).
##              Pseudobulk DESeq2-VST per patient (one point = one donor; Wilcoxon).
###############################################################################
suppressPackageStartupMessages({ library(DESeq2); library(SummarizedExperiment) })

.cache_panel_J_NRG1  <- './output/Figure_7/fig7_panel_J_endo_NRG1_pseudobulk.rds'
.cache_panel_J_MECOM <- './output/Figure_7/fig7_panel_L_arterial_MECOM_pseudobulk.rds'

.compute_subtype_pseudobulk_vst <- function(seurat_obj, subtype_label, gene) {
  .sub <- subset(seurat_obj, subset = Names == subtype_label)
  DefaultAssay(.sub) <- 'RNA'
  .agg <- Seurat::AggregateExpression(.sub, assays = 'RNA', slot = 'counts',
                                       group.by = c('patient','group'),
                                       return.seurat = FALSE)$RNA
  .cn <- colnames(.agg)
  .cd <- data.frame(sample = .cn,
                    patient = sub('_[^_]+$', '', .cn),
                    group   = sub('.*_', '', .cn),
                    stringsAsFactors = FALSE)
  .cd$group <- factor(.cd$group, levels = c('NF','pRV','RVF'))
  .dds <- DESeqDataSetFromMatrix(round(as.matrix(.agg)),
                                 colData = .cd, design = ~ group)
  .vsd  <- vst(.dds, blind = FALSE)
  .vmat <- as.matrix(SummarizedExperiment::assay(.vsd))
  data.frame(sample = .cd$sample, patient = .cd$patient,
             group  = .cd$group,  expr = .vmat[gene, ])
}

.nrg1_df <- if (file.exists(.cache_panel_J_NRG1)) {
  readRDS(.cache_panel_J_NRG1)
} else {
  d <- .compute_subtype_pseudobulk_vst(M1, 'Endocardial', 'NRG1')
  saveRDS(d, .cache_panel_J_NRG1); d
}
.mecom_df <- if (file.exists(.cache_panel_J_MECOM)) {
  readRDS(.cache_panel_J_MECOM)
} else {
  d <- .compute_subtype_pseudobulk_vst(M1, 'Arterial', 'MECOM')
  saveRDS(d, .cache_panel_J_MECOM); d
}

.violin_panel <- function(df, gene, subtype_lbl) {
  ggplot(df, aes(x = group, y = expr, fill = group)) +
    geom_violin(scale = 'width', alpha = 0.6, linewidth = PS$geom_lw) +
    geom_jitter(width = 0.08, size = 1.4, alpha = 0.9) +
    scale_fill_manual(values = disease_pal) +
    ggpubr::stat_compare_means(
      comparisons = list(c('NF','pRV'), c('pRV','RVF'), c('NF','RVF')),
      method = 'wilcox.test',
      size = PS$text_mm, family = FONT_FAMILY,
      symnum.args = pub_signif_args) +
    labs(x = NULL,
         y = bquote(italic(.(gene)) ~ 'expression (VST)'),
         title = paste0(subtype_lbl, ' (pseudobulk)')) +
    theme_v52(COMP_W) +
    theme(legend.position = 'none')
}

p_7J <- .violin_panel(.nrg1_df,  'NRG1',  'Endocardial')
p_7K <- .violin_panel(.mecom_df, 'MECOM', 'Arterial')
save_figure(p_7J, 'Figure_7_panel_J.pdf', width = 3.5, height = 4)
save_figure(p_7K, 'Figure_7_panel_K.pdf', width = 3.5, height = 4)


# 

#######################################
####  LEGACY: commented-out hdWGCNA setup (not produced for v57)
#######################################






















# DefaultAssay(M1) <- "RNA"


# M1 <- SetupForWGCNA(
#   M1,
#   gene_select = "fraction", # the gene selection approach
#   fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
#   wgcna_name = "ec_group" # the name of the hdWGCNA experiment
# )

# M1 <- MetacellsByGroups(seurat_obj = M1,
#   group.by = c("group"), # specify the columns in seurat_obj@meta.data to group by
#   reduction = 'harmony', # select the dimensionality reduction to perform KNN on
#   k = 25, # nearest-neighbors parameter
#   max_shared = 10, # maximum number of shared cells between two metacells
#   ident.group = 'group' # set the Idents of the metacell seurat object
# )

# M1 <- NormalizeMetacells(M1)


# M1 <- SetDatExpr(
#   M1,
#   group_name = c("NF","pRV","RVF"), # the name of the group of interest in the group.by column
#   group.by='group', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
#   assay = 'RNA', # using RNA assay
#   slot = 'data' # using normalized data
# )

# M1 <- TestSoftPowers(
#   M1,
#   networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
# )

# library(patchwork)
# plot_list <- PlotSoftPowers(M1)
# wrap_plots(plot_list, ncol=2)

# M1 <- ConstructNetwork(
#   M1,
#   tom_name = 'ec_net_group' # name of the topoligical overlap matrix written to disk
# )


# M1 <- ScaleData(M1)
# M1 <- ModuleEigengenes(
#  M1,
#  group.by.vars="patient"
# )


# MEs <- GetMEs(M1,harmonized=TRUE)

# M1 <- ModuleConnectivity(
#   M1,
#   group.by = 'group', group_name =  c("NF","pRV","RVF")
# )

# modules <- GetModules(M1)
# mods <- levels(modules$module); mods <- mods[mods != 'grey']
# saveRDS(M1,'./output/Figure_7/EC_hdWGCNA_by_group.rds')

###############################################################################
## v52 NEW PANELS — D through I
## hdWGCNA on EC + Phase 2 volcanos + Xenium spatial
###############################################################################
suppressPackageStartupMessages({
  library(hdWGCNA)
  library(WGCNA)
  library(DESeq2)
})

## ── Reuse the EC hdWGCNA object from earlier (already has 'ec' experiment) ──
ec_wgcna <- M1

.cache_mEs_ec <- './output/Figure_7/fig7_MEs_ec_cache.rds'
if (file.exists(.cache_mEs_ec)) {
  message('Loading cached GetMEs (ec_wgcna harmonized MEs)...')
  MEs_ec <- readRDS(.cache_mEs_ec)
} else {
  MEs_ec <- GetMEs(ec_wgcna, harmonized = TRUE)
  saveRDS(MEs_ec, .cache_mEs_ec)
}
mods_ec  <- colnames(MEs_ec)
mods_ec  <- mods_ec[mods_ec != 'grey']

## Map WGCNA color labels (turquoise, blue, brown, ...) to numeric M-labels.
## Convention: M_n where n = match(color, labels2colors(1:N)).
mapping_ec <- labels2colors(1:100)
mod_to_M   <- function(x) paste0('M', match(x, mapping_ec))

## Focus modules: M1 (capillary/Phase 1), M4 (arterial/angiogenic),
## M7 (activated/Phase 2). Map back to color names for module-table lookups
## but use the M-labels everywhere we render to a viewer.
focus_M    <- c(1, 4, 7)
focus_mods <- intersect(mapping_ec[focus_M], mods_ec)
if (length(focus_mods) == 0) focus_mods <- mods_ec[seq_len(min(3, length(mods_ec)))]
focus_mods_M <- mod_to_M(focus_mods)

## Add MEs to Seurat meta.data, plus M-named alias columns for plotting.
ec_wgcna@meta.data <- cbind(ec_wgcna@meta.data, MEs_ec)
for (.m in mods_ec) {
  ec_wgcna@meta.data[[mod_to_M(.m)]] <- ec_wgcna@meta.data[[.m]]
}

## ── Panel D — hdWGCNA Module UMAP (hub gene network) ─────────────────────
.cache_ec_wgcna_umap <- './output/Figure_7/fig7_ec_wgcna_module_umap_cache.rds'
tryCatch({
  if (file.exists(.cache_ec_wgcna_umap)) {
    message('Loading cached RunModuleUMAP (ec_wgcna)...')
    ec_wgcna <- readRDS(.cache_ec_wgcna_umap)
  } else {
    ec_wgcna <- RunModuleUMAP(ec_wgcna,
                               n_hubs      = 5,
                               n_neighbors = 10,
                               min_dist    = 0.3,
                               spread      = 2,
                               target_weight = 0.1,
                               supervised  = TRUE)
    saveRDS(ec_wgcna, .cache_ec_wgcna_umap)
  }
  p_7D <- ModuleUMAPPlot(ec_wgcna,
                          edge.alpha    = 0.25,
                          sample_edges  = TRUE,
                          edge_prop     = 0.1,
                          label_hubs    = 5,
                          keep_grey_edges = FALSE) +
    labs(title = 'EC hdWGCNA — module hub network') +
    theme_v52(COMP_W)
}, error = function(e) {
  message('ModuleUMAPPlot failed: ', conditionMessage(e), ' — using DimPlot fallback')
  p_7D <<- DimPlot(ec_wgcna, group.by = 'Names',
                   label = TRUE, repel = TRUE) +
    NoLegend() + labs(title = 'EC subtypes (UMAP)') + theme_v52(COMP_W)
})

save_figure(p_7D, 'Figure_7_panel_D.pdf', width = 6, height = 5)

## ── Panel E — GO enrichment for focus modules (clusterProfiler / enrichR) ─
modules_tbl <- GetModules(ec_wgcna)

.cache_panel_E_enrichr <- './output/Figure_7/fig7_panel_E_enrichr_cache.rds'
if (file.exists(.cache_panel_E_enrichr)) {
  message('Loading cached enrichr results (Panel E focus modules)...')
  p6E_enrichr_results <- readRDS(.cache_panel_E_enrichr)
} else {
  p6E_enrichr_results <- lapply(focus_mods, function(mod) {
    kme_col  <- paste0('kME_', mod)
    mod_rows <- modules_tbl$module == mod
    hub_genes <- modules_tbl$gene_name[mod_rows]
    if (kme_col %in% colnames(modules_tbl)) {
      hub_genes <- hub_genes[order(modules_tbl[[kme_col]][mod_rows], decreasing = TRUE)]
    }
    hub_genes <- hub_genes[1:min(50, length(hub_genes))]
    if (length(hub_genes) < 3) return(NULL)
    suppressPackageStartupMessages(library(enrichR))
    dbs_e <- c('GO_Biological_Process_2023')
    tryCatch(enrichr(hub_genes, databases = dbs_e)[[1]], error = function(e) NULL)
  })
  names(p6E_enrichr_results) <- focus_mods
  saveRDS(p6E_enrichr_results, .cache_panel_E_enrichr)
}

p_7E_list <- lapply(focus_mods, function(mod) {
  kme_col  <- paste0('kME_', mod)
  mod_rows <- modules_tbl$module == mod
  hub_genes <- modules_tbl$gene_name[mod_rows]
  if (kme_col %in% colnames(modules_tbl)) {
    hub_genes <- hub_genes[order(modules_tbl[[kme_col]][mod_rows], decreasing = TRUE)]
  }
  hub_genes <- hub_genes[1:min(50, length(hub_genes))]
  if (length(hub_genes) < 3) return(NULL)

  tryCatch({
    er  <- p6E_enrichr_results[[mod]]
    if (is.null(er)) stop('no result')
    er  <- er[er$Adjusted.P.value < 0.05, ]
    if (nrow(er) == 0) return(NULL)
    er$Term <- sub('\\s*\\(GO:\\d+\\)$', '', er$Term)
    er$Term <- factor(er$Term, levels = rev(er$Term[1:min(10, nrow(er))]))
    er <- er[1:min(10, nrow(er)), ]
    ggplot(er, aes(x = -log10(Adjusted.P.value), y = Term)) +
      geom_col(fill = lineage_pal['EC'], linewidth = PS$geom_lw) +
      labs(x = '-log\u2081\u2080(FDR)', y = NULL, title = mod_to_M(mod)) +
      theme_v52(COMP_W)
  }, error = function(e) {
    ggplot() +
      annotate('text', x = 0.5, y = 0.5,
               label = paste0('[', mod_to_M(mod), ' GO enrichment\n(enrichR unavailable)]'),
               size = PS$text_mm, family = FONT_FAMILY, colour = 'grey50', hjust = 0.5, vjust = 0.5) +
      theme_void()
  })
})
p_7E_list <- Filter(Negate(is.null), p_7E_list)
p_7E <- if (length(p_7E_list) > 0) {
  patchwork::wrap_plots(p_7E_list, ncol = length(p_7E_list))
} else {
  ggplot() + annotate('text', x = 0.5, y = 0.5,
                      label = '[Panel E: GO enrichment placeholder]',
                      size = PS$text_mm, family = FONT_FAMILY, colour = 'grey50', hjust = 0.5, vjust = 0.5) +
    theme_void()
}

save_figure(p_7E, 'Figure_7_panel_E.pdf', width = 10, height = 4)

## ── Panel F — Module DotPlot by EC subtype × disease group ────────────────
ec_wgcna$group <- factor(ec_wgcna$group, levels = c('NF', 'pRV', 'RVF'))
ec_wgcna$subtype_group <- paste(ec_wgcna$Names, ec_wgcna$group, sep = '_')

p_7F <- DotPlot(ec_wgcna,
                features   = focus_mods_M,
                group.by   = 'subtype_group',
                dot.scale  = 6) +
  scale_size(range = PS$dot_range) +
  coord_flip() +
  labs(x = 'Module', y = NULL, title = 'EC module scores by subtype × group') +
  theme_v52(COMP_W) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size = PS$legend_key)

save_figure(p_7F, 'Figure_7_panel_F.pdf', width = 9, height = 4)

## ── Panel G — FeaturePlot: module scores on EC UMAP ──────────────────────
p_7G_list <- lapply(focus_mods_M, function(mod_M) {
  tryCatch(
    FeaturePlot(ec_wgcna, features = mod_M, pt.size = PS$umap_pt) +
      scale_colour_gradient2(low = '#2166AC', mid = 'grey90', high = '#B2182B',
                             midpoint = 0, name = 'ME') +
      labs(title = mod_M) + NoAxes() + theme_v52(COMP_W),
    error = function(e) NULL)
})
p_7G_list <- Filter(Negate(is.null), p_7G_list)
p_7G <- patchwork::wrap_plots(p_7G_list, ncol = 3) +
  plot_annotation(title = 'EC module eigengene scores',
    theme = theme(plot.title = element_text(family = FONT_FAMILY, size = 7 * COMP_W / FINAL_WIDTH_IN, color = "black")))

save_figure(p_7G, 'Figure_7_panel_G.pdf', width = 10, height = 3.5)

## ── Panel H — EC EMT/angiogenic Phase-2 volcanos: snRNA-seq + Xenium ───────
## Stacked vertically (snRNA top, Xenium bottom) with a shared category legend.
## Xenium uses the **lineage** CSV restricted to EC (pan-pseudobulk), which is where
## the manuscript's quoted Phase-2 magnitudes (SNAI1 +2.94, MYC +2.58,
## ADAMTS1 +2.11, ANGPT2 +2.04) come from.
suppressPackageStartupMessages(library(EnhancedVolcano))

## Phase-2 EC categories (panel-filtered & coherent for adult EndMT-like
## activation rather than classical full EMT — see thread discussion).
##   Mesenchymal:    SNAI1 (master TF), FN1, ACTA2
##   Angiogenic:     ANGPT2, MYC, ADAMTS1
##   EC_activation:  MMRN1, PLA2G2A, ACKR1, PECAM1, VWF
##   Identity_TF:    MECOM
phase2_cats <- list(
  Mesenchymal   = c('SNAI1'),
  Angiogenic    = c('ANGPT2','MYC','ADAMTS1'),
  EC_activation = c('PLA2G2A','ACKR1','PECAM1','VWF'),
  EC_stress     = c('STC1','NR4A1','NR4A3'),
  Identity      = c('MECOM','CA4')
)
phase2_pal <- c(Mesenchymal   = '#B2182B', Angiogenic    = '#D6604D',
                EC_activation = '#2166AC', EC_stress     = '#1B9E77',
                Identity      = '#762A83', Other         = 'grey80')
phase2_lookup <- unlist(lapply(names(phase2_cats), function(k)
  setNames(rep(k, length(phase2_cats[[k]])), phase2_cats[[k]])))
phase2_genes <- names(phase2_lookup)

## Per the user's spec:
##   visual x-axis (post coord_flip) = -log10(padj) → 0 to 6 for both panels
##   visual y-axis (post coord_flip) = log2FC      → snRNA -3..3, Xenium -4..4
## Translate to EnhancedVolcano args (pre-flip):
##   ylim  = c(0, 6)        (will become x-axis after flip)
##   xlim  = c(-3, 3) sn    (will become y-axis after flip)
##   xlim  = c(-4, 4) xen   (will become y-axis after flip)
.make_ec_volcano <- function(df, label_text, fc_xlim, sig_ylim,
                              show_legend = TRUE) {
  df <- df[!is.na(df$padj) & !is.na(df$log2FoldChange), ]
  rownames(df) <- make.unique(as.character(df$gene))
  df$category  <- ifelse(df$gene %in% phase2_genes,
                         phase2_lookup[df$gene], 'Other')
  keyvals <- setNames(phase2_pal[df$category], df$category)
  sel     <- intersect(phase2_genes, df$gene)
  EnhancedVolcano(df, lab = rownames(df),
    x = 'log2FoldChange', y = 'padj',
    FCcutoff = 0.5, pCutoff = 0.05,
    xlim = fc_xlim, ylim = sig_ylim,
    selectLab = sel, max.overlaps = Inf, drawConnectors = TRUE,
    arrowheads = FALSE, widthConnectors = 0.3, caption = '',
    title    = label_text,
    subtitle = 'EC pseudobulk — RVF vs pRV (Phase 2)',
    legendPosition = if (show_legend) 'bottom' else 'none',
    colCustom = keyvals, labFace = 'italic') + coord_flip()
}

## --- snRNA-seq EC pan-pseudobulk (compute fresh, cache) ---------------------
.cache_panel_H_sn <- './output/Figure_7/fig7_panel_H_sn_ec_pan_pseudobulk_RVF_vs_pRV.rds'
if (file.exists(.cache_panel_H_sn)) {
  message('Loading cached snRNA pan-EC pseudobulk for Panel H...')
  sn_ec_de <- readRDS(.cache_panel_H_sn)
} else {
  message('Computing snRNA pan-EC pseudobulk DESeq2 (RVF vs pRV)...')
  suppressPackageStartupMessages({ library(DESeq2) })
  .agg <- Seurat::AggregateExpression(M1, assays = 'RNA', slot = 'counts',
                                       group.by = c('patient','group'),
                                       return.seurat = FALSE)$RNA
  .cn <- colnames(.agg)
  .cd <- data.frame(sample  = .cn,
                    patient = sub('_[^_]+$', '', .cn),
                    group   = sub('.*_', '', .cn),
                    stringsAsFactors = FALSE)
  .cd$group <- factor(.cd$group, levels = c('NF','pRV','RVF'))
  .dds <- DESeqDataSetFromMatrix(countData = round(as.matrix(.agg)),
                                  colData   = .cd, design = ~ group)
  .dds <- DESeq(.dds)
  .res <- lfcShrink(.dds, contrast = c('group','RVF','pRV'), type = 'ashr')
  sn_ec_de <- as.data.frame(.res); sn_ec_de$gene <- rownames(sn_ec_de)
  saveRDS(sn_ec_de, .cache_panel_H_sn)
  rm(.agg, .dds, .res); gc(verbose = FALSE)
}

## --- Xenium Activated_EC sublineage pseudobulk (precomputed CSV) ------------
xen_lin_de <- read.csv('./dependencies/shared/Xenium/xenium_pseudobulk_lineage_deseq2.csv')
xen_ec_de <- xen_lin_de[xen_lin_de$subtype == 'EC' &
                            xen_lin_de$contrast == 'RVF_vs_pRV', ]

p_7H_sn  <- .make_ec_volcano(sn_ec_de,    'snRNA-seq (EC pan-pseudobulk)',
                             fc_xlim = c(-3, 3), sig_ylim = c(0, 5),
                             show_legend = FALSE)
p_7H_xen <- .make_ec_volcano(xen_ec_de, 'Xenium (EC pan-pseudobulk, panel-only)',
                             fc_xlim = c(-2, 4), sig_ylim = c(0, 2),
                             show_legend = TRUE)
p_7H <- (p_7H_sn / p_7H_xen) + patchwork::plot_layout(guides = 'collect')

save_figure(p_7H, 'Figure_7_panel_H.pdf', width = 8, height = 10)

## ── Panel I — ANGPT2 and SNAI1 induction in RVF EC (pseudobulk violins) ───
## Replaces the prior Xenium spatial maps with patient-level pseudobulk
## DESeq2-VST violins (matches Panel J's design). EC = full M1 EC object
## (all 5 subtypes pooled).
.cache_panel_I <- './output/Figure_7/fig7_panel_I_ec_pseudobulk_angpt2_snai1.rds'

if (file.exists(.cache_panel_I)) {
  message('Loading cached EC pseudobulk frame for Panel I...')
  .pi_df <- readRDS(.cache_panel_I)
} else {
  message('Computing EC pseudobulk DESeq2-VST for ANGPT2/SNAI1...')
  ec_full <- M1
  DefaultAssay(ec_full) <- 'RNA'

  .agg <- Seurat::AggregateExpression(
    ec_full, assays = 'RNA', slot = 'counts',
    group.by = c('patient','group'), return.seurat = FALSE)$RNA
  .cn <- colnames(.agg)
  .cd <- data.frame(
    sample  = .cn,
    patient = sub('_[^_]+$', '', .cn),
    group   = sub('.*_', '', .cn),
    stringsAsFactors = FALSE)
  .cd$group <- factor(.cd$group, levels = c('NF','pRV','RVF'))

  .dds <- DESeqDataSetFromMatrix(
    countData = round(as.matrix(.agg)),
    colData   = .cd, design = ~ group)
  .vsd  <- vst(.dds, blind = FALSE)
  .vmat <- as.matrix(SummarizedExperiment::assay(.vsd))

  .pi_df <- data.frame(
    sample  = rep(.cd$sample, 2),
    patient = rep(.cd$patient, 2),
    group   = rep(.cd$group, 2),
    gene    = rep(c('ANGPT2','SNAI1'), each = nrow(.cd)),
    expr    = c(.vmat['ANGPT2', ], .vmat['SNAI1', ]),
    stringsAsFactors = FALSE)
  .pi_df$gene <- factor(.pi_df$gene, levels = c('ANGPT2','SNAI1'))
  saveRDS(.pi_df, .cache_panel_I)
  rm(.agg, .dds, .vsd, .vmat); gc(verbose = FALSE)
}

p_7I <- ggplot(.pi_df, aes(x = group, y = expr, fill = group)) +
  geom_violin(scale = 'width', alpha = 0.6, linewidth = PS$geom_lw) +
  geom_jitter(width = 0.08, size = 1.4, alpha = 0.9) +
  scale_fill_manual(values = disease_pal) +
  facet_wrap(~ gene, ncol = 2, scales = 'free_y') +
  ggpubr::stat_compare_means(
    comparisons = list(c('NF','pRV'), c('pRV','RVF'), c('NF','RVF')),
    method = 'wilcox.test',
    size = PS$text_mm, family = FONT_FAMILY,
    symnum.args = pub_signif_args) +
  labs(x = NULL, y = 'Expression (VST)',
       title = 'EC pseudobulk — ANGPT2 / SNAI1') +
  theme_v52(COMP_W) +
  theme(legend.position = 'none',
        strip.text = element_text(face = 'bold.italic'))

save_figure(p_7I, 'Figure_7_panel_I.pdf', width = 6, height = 4)
rm(ec_wgcna, MEs_ec, modules_tbl); gc()

## ── Assemble panels D-K (composite is a convenience; per-panel saves above) ─
tryCatch({
  fig7_dk <- (p_7D | p_7E | p_7F) /
             (p_7G | p_7H | p_7I) /
             (p_7J | p_7K) +
    plot_annotation(tag_levels = list(c('D','E','F','G','H','I','J','K')),
      theme = theme(plot.tag = element_text(family = FONT_FAMILY,
                                            size = 16 * COMP_W / FINAL_WIDTH_IN,
                                            face = "bold")))
  save_figure(fig7_dk, 'Figure_7_panels_D-K.pdf',
              width = COMP_W, height = COMP_H + 5)
}, error = function(e) {
  message('Composite assembly failed (likely patchwork viewport): ', e$message)
  message('Per-panel PDFs are still complete.')
})
message('Figure 7 (v53) per-panel + composite PDFs written to ', V52_FIG_DIR)



