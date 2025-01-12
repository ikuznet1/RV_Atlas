
library(Seurat)
library(hdWGCNA)
library(ggeasy)
library(harmony)
library('tximportData')
library('tximport')
library(viridis)



source('~/Downloads/hdWGCNA_TOM/spatial_functions.R')

#######################################
######### Consensus hdWGCNA #########
#######################################

M1 <- readRDS(file = "/Volumes/Extreme\ SSD/Final_Analysis/CellTypes/myeloid_subclust.rds")

M1 <- FindNeighbors(M1)
M1 <- FindClusters(M1,resolution=1)

new.cluster.ids <- c('CCR2- rMac2','CCR2- rMac1','CCR2- rMac2','CCR2- rMac2',
	'CCR2+ rMac','CCR2- rMac1','CCR2- rMac1',
	'CCR2- rMac1','DCs','iMac','TREM2 Mac','CCR2- rMac2')
names(new.cluster.ids) <- levels(M1)
M1 <- RenameIdents(M1, new.cluster.ids)

M2 <- subset(M1,idents=c('CCR2+ rMac','DCs'))
M2<-RunPCA(M2)
M2 <- M2%>% 
  RunUMAP(reduction = "pca", dims = 1:50, verbose = F) %>% 
  FindNeighbors(reduction = "pca", dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()

#0 is CD1c- CD163+ MERTK+
#1 is CCR2+ CD11c+
#2 is DC
#3 is mix of classical (FCN1+), non-classical (LILRB2, ITGAL, FCGR3A), and neutrophil-like (CSF3R)

M1$Subnames<-M1@active.ident
M1$Subsubnames <- M1$Subnames
cells_DC = colnames(M2)[M2$seurat_clusters==2]
cells_mono = colnames(M2)[M2$seurat_clusters==3]
cells_rMac = union(colnames(M2)[M2$seurat_clusters==0],colnames(M2)[M2$seurat_clusters==1])

levels(M1$Subsubnames) <- c(levels(M1$Subsubnames), 'Mono')
M1$Subsubnames[colnames(M1) %in% cells_DC] = 'DCs'
M1$Subsubnames[colnames(M1) %in% cells_mono] = 'Mono'
M1$Subsubnames[colnames(M1) %in% cells_rMac] = 'CCR2+ rMac'


M2 <- M1
M2<-RunPCA(M2)
M2 <- RunHarmony(M1,'patient')
M2 <- M2%>% 
  RunUMAP(reduction = "harmony", dims = 1:50, verbose = F) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:50) %>% 
  FindClusters(resolution=1) %>% 
  identity()

M1 <- M2
M1 <- SetIdent(M1,value='Subsubnames')

M1 <- AddModuleScore(M1, features=list(c("TREM2","GPNMB","MITF","SPP1")),assay="SCT",name="TREM2_Mac_Score")
M1 <- AddModuleScore(M1, features=list(c("CLEC9A","ZBTB46","CD1C","CD226")),assay="SCT",name="DC_Score")
M1 <- AddModuleScore(M1, features=list(c("FCN1","LILRB2","ITGAL","CSF3R")),assay="SCT",name="Mono_Score")
M1 <- AddModuleScore(M1, features=list(c("CCR2","CX3CR1","ITGAX")),assay="SCT",name="CCR2+_rMac_Score")
M1 <- AddModuleScore(M1, features=list(c("LYVE1","FOLR2","SIGLEC1","F13A1")),assay="SCT",name="CCR2-_rMac1_Score")
M1 <- AddModuleScore(M1, features=list(c("RBMS3","PLA2G5","EBF1")),assay="SCT",name="CCR2-_rMac2_Score")
M1 <- AddModuleScore(M1, features=list(c("IL1B","CCL3","CCL4","CXCL3","CXCL8")),assay="SCT",name="iMac_Score")





M2 <- readRDS('~/Downloads/hdWGCNA_TOM/PAB_data_clean.rds')

M2 <- SetIdent(M2, value = "Names")
M2$group <- M2$orig.ident


M3 <- subset(M2, Names %in% c("Myeloid"))

M3 <- RunPCA(M3)
M3 <- FindNeighbors(M3, dims = 1:50)
M3 <- FindClusters(M3, resolution = 0.5)
M3 <- RunUMAP(M3, dims = 1:10)

#0 is CCR2- rMac
#1 looks like cardiac contam
#2 is
#3 is DC and CCR, HLA high
#4 B cells
#5 is something non immune
#6 EC

M3$Subnames <- M3@active.ident
M3 <- subset(M3, Subnames %in% c('0','1','2','3') )
M3 <- RunUMAP(M3, dims = 1:50)

#0 rMac
#1
#2 
#3 HLA
M3 <- subset(M3, Subnames %in% c('0','3') )
M3 <- RunPCA(M3)
M3 <- RunUMAP(M3, dims = 1:10)

labels <- c('rMac','HLA')

names(labels) <- levels(M3)
M3 <- RenameIdents(M3, labels)
M3$Subnames <- M3@active.ident

human2mouse <- read.csv('~/Downloads/hdWGCNA_TOM/human2mouse.csv',header=F)
idx <- match(unique(human2mouse[,2]),human2mouse[,2])
human2mouse<-human2mouse[idx,]
colnames(human2mouse) <-c('human_name', 'mouse_name')

RNA <- M3@assays$RNA['counts']
newnames <- human2mouse$human_name[match(rownames(RNA),human2mouse$mouse_name)]
#newnames[is.na(newnames)] <- rownames(RNA)[is.na(newnames)]
rownames(RNA) <- newnames


M2 <- CreateAssayObject(RNA[!is.na(newnames),])
M3[['humanized']] <- M2
DefaultAssay(M3) <- "humanized"
M3[['SCT']] <- NULL

M3 <- SCTransform(M3, vst.flavor = "v2",assay='humanized',variable.features.n = 2000)
M3 <- RunPCA(M3, npcs = 50, verbose = FALSE)

M3$Species <- 'mouse'
M1$Species <- 'human'

DefaultAssay(M3) <- "humanized"
DefaultAssay(M1) <- "RNA"
M3[['SCT']] <- NULL
M1[['SCT']] <- NULL
M3$Subsubnames <- M3$Subnames


shared <- intersect(rownames(M1),rownames(M3))


M3[['RNA']] <- NULL
M3[['RNA']] <- M3[['humanized']]
DefaultAssay(M3) <- "RNA"
M3[['humanized']] <- NULL

M1$decontXcounts <- NULL

M2 <- merge(M1[shared,],M3[shared,])

M2 <- NormalizeData(M2)
M2<-FindVariableFeatures(M2,nfeatures=3000)
M2<-ScaleData(M2)
M2 <- RunPCA(M2)
M2 <- RunHarmony(M2, group.by.vars = c('Species','patient'), dims=1:50,lambda=c(1,0.01))

M2 <- RunUMAP(M2, reduction='harmony',dims = 1:30)

DimPlot(M2, group.by='patient', split.by='Species', raster=FALSE, label=TRUE) + umap_theme()
DimPlot(M2, group.by='patient', split.by='Subsubnames', raster=FALSE, label=TRUE) + umap_theme()

M2$Names <- M2$Subsubnames
M2$Names[M2$Names == 'CCR2- rMac2'] <- 'CCR2- rMac'
M2$Names[M2$Names == 'CCR2- rMac1'] <- 'CCR2- rMac'
M2$Names[M2$Names == 'rMac'] <- 'CCR2- rMac'
M2$Names[M2$Names == 'DCs'] <- 'HLA'
M2$Names[M2$Names == 'Mono'] <- 'HLA'
M2$Names[M2$Names == 'CCR2+ rMac'] <- 'HLA'
M2$Names[M2$Names == 'iMac'] <- 'HLA'




M2 <- SetupForWGCNA(
  M2,
  gene_select = "fraction",
  fraction = 0.05,
  wgcna_name = 'Myeloid_consensus'
)

# construct metacells:
M2 <- MetacellsByGroups(
  M2,
  group.by = c("Names", 'Species'),
  k = 25,
  max_shared = 12,
  min_cells = 50,
  target_metacells = 250,
  reduction = 'harmony',
  ident.group = 'Names'
)

M2 <- NormalizeMetacells(M2)

# setup expression matrices for each species in astrocytes
M2 <- SetMultiExpr(
  M2,
  group_name = "CCR2- rMac",
  group.by = "Names",
  multi.group.by ="Species",
  multi_groups = NULL
)


# identify soft power thresholds
M2 <- TestSoftPowersConsensus(M2)

# plot soft power results
plot_list <-  PlotSoftPowers(M2)
consensus_groups <- unique(M2$Species)
p_list <- lapply(1:length(consensus_groups), function(i){
  cur_group <- consensus_groups[[i]]
  plot_list[[i]][[1]] + ggtitle(paste0(cur_group)) + theme(plot.title=element_text(hjust=0.5))
})
library(patchwork)
wrap_plots(p_list, ncol=2)

# consensus network analysis
M2 <- ConstructNetwork(
  M2,soft_power=c(7,7),
  consensus=TRUE,
  tom_name = "Species_Consensus",
  overwrite_tom = TRUE
)

PlotDendrogram(M2, main='Resident macrophage cross species dendrogram')


M2 <- ModuleEigengenes(M2, group.by.vars=c("Species","patient"))
M2 <- ModuleConnectivity(M2, group_name ='CCR2- rMac', group.by='Names')

# re-name modules
M2 <- ResetModuleNames(M2, new_name = "rMac-CM")

# visualize network with UMAP
M2 <- RunModuleUMAP(
  M2,
  n_hubs = 5,
  n_neighbors=10,
  min_dist=0.3,
  spread=2,
  target_weight=0.1,
  supervised=TRUE
)

# get the hub gene UMAP table from the seurat object
umap_df <- GetModuleUMAP(M2)

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
    ) + ggtitle('HLA Consensus WGCNA') +
    ggrepel::geom_text_repel(data = centroid_df, label=centroid_df$cluster, color='black', max.overlaps=Inf, size=3)



hub_genes <- GetHubGenes(M2, 3)

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
    ggrepel::geom_text_repel(data = centroid_df, label=centroid_df$cluster, color='black', max.overlaps=Inf, size=3, fontface='bold') +
    geom_text_repel(label=plot_df$anno, max.overlaps=Inf, color='black', fontface='italic', size=3) +
    umap_theme() + NoLegend() +
    coord_equal() +
    theme(
      plot.margin = margin(0,0,0,0)
    )

pdf(paste0('~/Downloads/hdWGCNA_TOM/rMac_hubgene_umap_ggplot.pdf'), width=8, height=8)
print(p)
dev.off()

dbs <-c('GO_Biological_Process_2023','GO_Cellular_Component_2023','GO_Molecular_Function_2023','Reactome_2022', 'ChEA_2022',"LINCS_L1000_Chem_Pert_up")


M2 <- RunEnrichr(
  M2,
  dbs=dbs, # character vector of enrichr databases to test
  max_genes = Inf # number of genes per module to test. use max_genes = Inf to choose all genes!
)

enrich_df <- GetEnrichrTable(M2)

EnrichrBarPlot(
  M2,
  outdir = "~/Downloads/hdWGCNA_TOM/Myeloid_rMac_consensus_modules", # name of output directory
  n_terms = 10, # number of enriched terms to show (sometimes more show if there are ties!!!)
  plot_size = c(5,7), # width, height of the output .pdfs
  logscale=TRUE # do you want to show the enrichment as a log scale?
)

modules <- GetModules(M2) %>% subset(module != 'grey')

cat(rownames(subset(modules,module=='rMac-CM1')),sep='\n')
cat(rownames(subset(modules,module=='rMac-CM2')),sep='\n')
cat(rownames(subset(modules,module=='rMac-CM6')),sep='\n')
cat(rownames(subset(modules,module=='rMac-CM8')),sep='\n')
cat(rownames(subset(modules,module=='rMac-CM9')),sep='\n')
cat(rownames(subset(modules,module=='rMac-CM11')),sep='\n')
cat(rownames(subset(modules,module=='rMac-CM12')),sep='\n')

write.csv(modules,"~/Downloads/hdWGCNA_TOM/consensusWGCNA_rMac_modules.csv")

modules <- GetModules(M2)
color_df <- modules %>% subset(module!='grey') %>%
  select(c(module, color)) %>% distinct %>%
  rename(c(group=module, colour=color))
mods <- levels(modules$module)
mods <- mods[mods!='grey']

# helper function to wrap text
wrapText <- function(x, len) {
    sapply(x, function(y) paste(strwrap(y, len), collapse = "\n"), USE.NAMES = FALSE)
}

combined_output <- GetEnrichrTable(M2)
selected_terms <- subset(combined_output,db=="GO_Biological_Process_2023")

# subset selected terms
selected_terms <- subset(selected_terms, Term %in% selected_terms$Term & P.value < 0.05)
idx_top_1 <- match(unique(selected_terms$module),selected_terms$module)
idx_top_3 <- sort(c(idx_top_1,idx_top_1+1,idx_top_1+2))

selected_terms<-selected_terms[idx_top_1,]
#key_terms <- read.csv('~/Downloads/hdWGCNA_TOM/bulkRNA_GOterms_ofinterest.csv')
#selected_terms <- subset(selected_terms,Term %in% key_terms[[1]])
#selected_terms <- subset(combined_output, Term %in% selected_terms$Term & P.value < 0.05)


selected_terms$group <- factor(
  as.character(selected_terms$module),
  levels = mods
)


# set max pval
quantile(-log(selected_terms$P.value), 0.95)
max_p <- 10

selected_terms$logp <- -log(selected_terms$P.value)
selected_terms$logp <- ifelse(selected_terms$logp > max_p, max_p, selected_terms$logp)

# remove GO Term ID
library(stringr)
selected_terms$Term <- str_replace(selected_terms$Term, " \\s*\\([^\\)]+\\)", "")

selected_terms <- selected_terms %>%
  arrange(group)


selected_terms$wrap <- wrapText(selected_terms$Term, 35)

selected_terms$Term <- factor(
  as.character(selected_terms$Term),
  levels = rev(unique(as.character(selected_terms$Term)))
)

library(viridis)

# GO Term dot plot
p <- selected_terms %>%
  ggplot(aes(x = group, y = Term, color =logp, size=log(Combined.Score))) +
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
#color_df$group<-paste0('M',match(color_df$group,mapping))
lvls <- stringr::str_sort(unique(color_df$group), numeric = TRUE)
color_df$group <- factor(color_df$group, levels = lvls)
# make the colorbar as its own heatmap
color_df$var <- 1
colorbar <- color_df %>%
  ggplot(aes(x=group, y=var, fill=group)) +
  geom_tile() +
  scale_fill_manual(values=color_df$colour) +
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


pdf(paste0('~/Downloads/hdWGCNA_TOM/Consensus_WGNCA_rMac_selected_GO_terms.pdf'), width=13, height=8)
p / colorbar #+ plot_layout(heights=c(20,1))
dev.off()



### HLA

M2 <- SetMultiExpr(
  M2,
  group_name = "HLA",
  group.by = "Names",
  multi.group.by ="Species",
  multi_groups = NULL
)


# consensus network analysis
M2 <- ConstructNetwork(
  M2,soft_power=c(7,7),
  consensus=TRUE,
  tom_name = "Species_HLA_Consensus",
  overwrite_tom = TRUE
)

PlotDendrogram(M2, main='HLA macrophage cross species dendrogram')

M2 <- ModuleEigengenes(M2, group.by.vars=c("Species","patient"))
M2 <- ModuleConnectivity(M2, group_name ='HLA', group.by='Names')

# re-name modules
M2 <- ResetModuleNames(M2, new_name = "HLA-CM")

# visualize network with UMAP
M2 <- RunModuleUMAP(
  M2,
  n_hubs = 5,
  n_neighbors=10,
  min_dist=0.3,
  spread=2,
  target_weight=0.1,
  supervised=TRUE
)

# get the hub gene UMAP table from the seurat object
umap_df <- GetModuleUMAP(M2)

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
    ) + ggtitle('HLA Consensus WGCNA') +
    ggrepel::geom_text_repel(data = centroid_df, label=centroid_df$cluster, color='black', max.overlaps=Inf, size=3)



hub_genes <- GetHubGenes(M2, 3)

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
    ggrepel::geom_text_repel(data = centroid_df, label=centroid_df$cluster, color='black', max.overlaps=Inf, size=3, fontface='bold') +
    geom_text_repel(label=plot_df$anno, max.overlaps=Inf, color='black', fontface='italic', size=3) +
    umap_theme() + NoLegend() +
    coord_equal() +
    theme(
      plot.margin = margin(0,0,0,0)
    )

pdf(paste0('~/Downloads/hdWGCNA_TOM/HLA_hubgene_umap_ggplot.pdf'), width=8, height=8)
print(p)
dev.off()



dbs <-c('GO_Biological_Process_2023','GO_Cellular_Component_2023','GO_Molecular_Function_2023','Reactome_2022', 'ChEA_2022',"LINCS_L1000_Chem_Pert_up")


M2 <- RunEnrichr(
  M2,
  dbs=dbs, # character vector of enrichr databases to test
  max_genes = Inf # number of genes per module to test. use max_genes = Inf to choose all genes!
)

enrich_df <- GetEnrichrTable(M2)

EnrichrBarPlot(
  M2,
  outdir = "~/Downloads/hdWGCNA_TOM/Myeloid_HLA_consensus_modules", # name of output directory
  n_terms = 10, # number of enriched terms to show (sometimes more show if there are ties!!!)
  plot_size = c(5,7), # width, height of the output .pdfs
  logscale=TRUE # do you want to show the enrichment as a log scale?
)

modules <- GetModules(M2) %>% subset(module != 'grey')

write.csv(modules,"~/Downloads/hdWGCNA_TOM/consensusWGCNA_HLA_modules.csv")

cat(rownames(subset(modules,module=='HLA-CM10')),sep='\n')
cat(rownames(subset(modules,module=='HLA-CM20')),sep='\n')


modules <- GetModules(M2)
color_df <- modules %>% subset(module!='grey') %>%
  select(c(module, color)) %>% distinct %>%
  rename(c(group=module, colour=color))
mods <- levels(modules$module)
mods <- mods[mods!='grey']

# helper function to wrap text
wrapText <- function(x, len) {
    sapply(x, function(y) paste(strwrap(y, len), collapse = "\n"), USE.NAMES = FALSE)
}

combined_output <- GetEnrichrTable(M2)
selected_terms <- subset(combined_output,db=="GO_Biological_Process_2023")

# subset selected terms
selected_terms <- subset(selected_terms, Term %in% selected_terms$Term & P.value < 0.05)
idx_top_1 <- match(unique(selected_terms$module),selected_terms$module)
idx_top_3 <- sort(c(idx_top_1,idx_top_1+1,idx_top_1+2))

selected_terms<-selected_terms[idx_top_1,]
#key_terms <- read.csv('~/Downloads/hdWGCNA_TOM/bulkRNA_GOterms_ofinterest.csv')
#selected_terms <- subset(selected_terms,Term %in% key_terms[[1]])
#selected_terms <- subset(combined_output, Term %in% selected_terms$Term & P.value < 0.05)


selected_terms$group <- factor(
  as.character(selected_terms$module),
  levels = mods
)


# set max pval
quantile(-log(selected_terms$P.value), 0.95)
max_p <- 10

selected_terms$logp <- -log(selected_terms$P.value)
selected_terms$logp <- ifelse(selected_terms$logp > max_p, max_p, selected_terms$logp)

# remove GO Term ID
library(stringr)
selected_terms$Term <- str_replace(selected_terms$Term, " \\s*\\([^\\)]+\\)", "")

selected_terms <- selected_terms %>%
  arrange(group)


selected_terms$wrap <- wrapText(selected_terms$Term, 35)

selected_terms$Term <- factor(
  as.character(selected_terms$Term),
  levels = rev(unique(as.character(selected_terms$Term)))
)

library(viridis)

# GO Term dot plot
p <- selected_terms %>%
  ggplot(aes(x = group, y = Term, color =logp, size=log(Combined.Score))) +
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
#color_df$group<-paste0('M',match(color_df$group,mapping))
lvls <- stringr::str_sort(unique(color_df$group), numeric = TRUE)
color_df$group <- factor(color_df$group, levels = lvls)
# make the colorbar as its own heatmap
color_df$var <- 1
colorbar <- color_df %>%
  ggplot(aes(x=group, y=var, fill=group)) +
  geom_tile() +
  scale_fill_manual(values=color_df$colour) +
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


pdf(paste0('~/Downloads/hdWGCNA_TOM/Consensus_WGNCA_HLA_selected_GO_terms.pdf'), width=13, height=8)
p / colorbar #+ plot_layout(heights=c(20,1))
dev.off()


#######################################
#########  MYELOID PAB PILOT  #########
#######################################


#1343 outlier
#1632,1467 RVF
#1681 1691 1697 NF

M1 <- readRDS('~/Downloads/hdWGCNA_TOM/PAB_data_clean.rds')

M1 <- SetIdent(M1, value = "Names")
M1$group <- M1$orig.ident


M2 <- subset(M1, Names %in% c("Myeloid"))

M2 <- RunPCA(M2)
M2 <- FindNeighbors(M2, dims = 1:50)
M2 <- FindClusters(M2, resolution = 0.5)
M2 <- RunUMAP(M2, dims = 1:10)

#0 is CCR2- rMac
#1 looks like cardiac contam
#2 is
#3 is DC and CCR, HLA high
#4 B cells
#5 is something non immune
#6 EC

M2$Subnames <- M2@active.ident
M2 <- subset(M2, Subnames %in% c('0','1','2','3') )
M2 <- RunUMAP(M2, dims = 1:50)

#0 rMac
#1
#2 
#3 HLA
M2 <- subset(M2, Subnames %in% c('0','3') )
M2 <- RunPCA(M2)
M2 <- RunUMAP(M2, dims = 1:10)

labels <- c('rMac','HLA')

names(labels) <- levels(M2)
M2 <- RenameIdents(M2, labels)
M2$Subnames <- M2@active.ident




#### LOAD BULK DATA AND EMBED INTO ALL PAB
dir <- '/Users/ikuz/Documents/Mouse_PAB_myeloid/output'

samples <- read.csv(file.path(dir, "meta.csv"), header = TRUE)
samples$Batch <- factor(samples$Batch)
samples$Pressure.Loading <- factor(samples$Pressure.Loading)


files <- file.path(dir,'nascent', samples$ID, "abundance.h5")
names(files) <- samples$ID
tx2gene <- read.table('/Users/ikuz/Documents/Mouse_PAB_myeloid/index/t2g.txt',fill=T)
tx2gene = data.frame(TXNAME=tx2gene$V1,GENEID=tx2gene$V3)

txi.kallisto <- tximport(files, type = "kallisto", tx2gene = tx2gene)


rownames(samples) <- samples$ID
bulk <- CreateSeuratObject(counts = txi.kallisto$counts, meta.data = data.frame(samples))

bulk <- SCTransform(bulk, vst.flavor = "v2")
bulk <- RunPCA(bulk, npcs = 10, verbose = FALSE)


M1 <- SplitObject(M1, split.by = "patient")
M1<-PrepSCTIntegration(M1)
features<-SelectIntegrationFeatures(M1)
M1.anchors<-FindIntegrationAnchors(M1,normalization.method = 'SCT',anchor.features = features, reference = c(1,2,3), reduction = "rpca")
M1 <- IntegrateData(anchorset = M1.anchors,normalization.method='SCT')

DefaultAssay(M1) <- "integrated"

M1 <- RunPCA(M1, npcs = 50, verbose = FALSE)
M1 <- RunUMAP(M1, reduction = "pca", dims = 1:30)


anchors <- FindTransferAnchors(
  reference = M1,
  query = bulk,
  normalization.method = "SCT",
  recompute.residuals=FALSE,
  reference.reduction = "pca",
  dims = 1:50,
  k.score =15
)


predictions <- TransferData(anchorset = anchors, refdata = M1$Names, dims = 1:50,k.weight=15)

M1 <- RunUMAP(M1, dims = 1:50, return.model = TRUE)

bulk <- MapQuery(anchorset = anchors, reference = M1, query = bulk,
	refdata = list(celltype = "Names"), reference.reduction = "pca", 
	reduction.model = "umap",transferdata.args = list(k.weight=15))

#score <- MappingScore(anchors)

#bulk$map_score <- score

p1 <- DimPlot(M1, reduction = "umap", group.by = "Names", label = TRUE, label.size = 3, repel = TRUE,raster=TRUE,pt.size=1.5) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(bulk, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE, label.size = 3, pt.size=1.5,repel = TRUE,raster=TRUE) + NoLegend() + ggtitle("Query transferred labels")
pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'RV_PAB_bulkMyeloid_ref_mapped.pdf'), width=10, height=5)
p1 + p2
dev.off()



#### LOAD BULK DATA AND EMBED INTO MYELOID PAB

bulk <- CreateSeuratObject(counts = txi.kallisto$counts, meta.data = data.frame(samples))

bulk <- SCTransform(bulk, vst.flavor = "v2")
bulk <- RunPCA(bulk, npcs = 10, verbose = FALSE)


#M2 <- SplitObject(M2, split.by = "patient")
#M2<-PrepSCTIntegration(M2)
#features<-SelectIntegrationFeatures(M2)
#M2.anchors<-FindIntegrationAnchors(M2,normalization.method = 'SCT',anchor.features = features, reduction = "rpca")
#M2 <- IntegrateData(anchorset = M2.anchors,normalization.method='SCT')

#DefaultAssay(M2) <- "integrated"

#M1 <- RunPCA(M2, npcs = 50, verbose = FALSE)
#M1 <- RunUMAP(M2, reduction = "pca", dims = 1:30)
M3 <- subset(M1, Names %in% c("Myeloid"))

true_myeloid <- colnames(M2)
M3 <- subset(M3,cells=true_myeloid)
M3$Subnames = M2$Subnames

anchors <- FindTransferAnchors(
  reference = M3,
  query = bulk,
  normalization.method = "SCT",
  recompute.residuals=FALSE,
  reference.reduction = "pca",
  dims = 1:50,
  k.score =15
)


predictions <- TransferData(anchorset = anchors, refdata = M3$Subnames, dims = 1:50,k.weight=10)

M3 <- RunUMAP(M3, dims = 1:50, return.model = TRUE)

bulk <- MapQuery(anchorset = anchors, reference = M3, query = bulk,
	refdata = list(celltype = "Subnames"), reference.reduction = "pca", 
	reduction.model = "umap",transferdata.args = list(k.weight=10))

#score <- MappingScore(anchors)

#bulk$map_score <- score


bulk$group <- paste0(bulk$Origin,'_',bulk$Type)
p1 <- DimPlot(M3, reduction = "umap", group.by = "Subnames", label = TRUE, label.size = 3, repel = TRUE,raster=TRUE,pt.size=1.5) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(bulk, reduction = "ref.umap", group.by = "group", label = TRUE, label.size = 3, pt.size=1.5,repel = TRUE,raster=TRUE) + ggtitle("Query transferred labels")
pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'RV_PAB_Myeloid_Only_bulkMyeloid_ref_mapped.pdf'), width=10, height=5)
p1 + p2
dev.off()

bulk_RV = subset(bulk, Origin=='RV')

bulk_RV<-SetIdent(bulk_RV,value='Type')

M3<-SetIdent(M3,value='group')


DefaultAssay(M3) <- 'SCT'

a <- FindMarkers(bulk_RV,ident.1='PAB',ident.2='Sham',min.pct=0.1,logfc.threshold=0)

M3 <- PrepSCTFindMarkers(M3)
b <-  FindMarkers(M3,ident.1='Sev',ident.2='Nor',min.pct=0.25,recorrect_umi=F,logfc.threshold=0)
shared <- intersect(rownames(a),rownames(b))
dataset <- data.frame(bulk=a[shared,]$avg_log2FC,sn=b[shared,]$avg_log2FC)
rownames(dataset) <- shared
labs <- rownames(dataset)

cor(dataset[,1],dataset[,2])

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'RV_PAB_Myeloid_Only_bulkMyeloid_scatter.pdf'), width=10, height=10)
ggplot(dataset, aes(x = sn, y=bulk)) + geom_point() + 
  geom_text_repel(label=labs,max.overlaps=25) + theme_classic()
dev.off()

cat(rownames(subset(dataset, bulk > 0.1 & sn > 0.1)),sep='\n')
cat(rownames(subset(dataset, bulk < 0.1 & sn < 0.1)),sep='\n')



###Compare to human RV


bulk_RV <- AddModuleScore(bulk_RV,list(c('Ciita','Cd74','H2-Ab1','H2-Aa',
	'H2-Eb1','H2-Eb2','H2-Ob','H2-DMb1','H2-DMb2','H2-DMa','H2-Oa')),name='MHCII')

VlnPlot(bulk_RV,'MHCII1')

human2mouse <- read.csv('~/Downloads/hdWGCNA_TOM/human2mouse.csv',header=F)
idx <- match(unique(human2mouse[,2]),human2mouse[,2])
human2mouse<-human2mouse[idx,]
colnames(human2mouse) <-c('human_name', 'mouse_name')
gluc_response <- "VIT;VKORC1L1;ERRFI1;AHCYL1;STEAP4;TRAF3IP2;GALNT15;SERPINE1;JADE1;SLA;CBLB;MT1X;EPS8;CCND3;BMPER;RASSF4;RPS6KA2;ANPEP;C1RL;MAP3K6;IL6R;PDGFRA;MLIP;SCARA5;IL1R1;EBF1;TTC7A;CRISPLD2;SPARCL1;FKBP5;NNMT;LPAR1;SLC1A3;PLA2G5;NID1;ACACB;ZFP36L2;PIK3R5;C3;SCFD2;LPXN;HACL1;SRGAP2;SLC38A2;SLC19A2;S100A10;KLHL29;GADD45B;ZBTB16;ELL2;CORO2B;IGF2R;NFATC4;DERA;SULT1B1;MAFB;BCL6;TMEM236;TBXAS1;NDUFAF2;RGL3;SERPINA3;MCFD2;PTPRS;ELN;PTEN;FMN1;HIF3A;TFCP2L1;PTH1R;SYNE3;CTSS;PTPRG;RNF157;ADAMTS2;C1QTNF1;IMPA2;SH3PXD2B;FLVCR2;EFHD1;AOX1;CERS6;ZHX3;KLF13;ANXA2;IFNGR1;GPX3;NCOA3;SLC39A11;NGF;OSMR;SLC39A14;TGFBR2;TGFBR3;PSMA6;ARHGAP10;MMP14;TBC1D2;SLC7A7;SLC7A8;GFOD1;DPYD;PICK1;FAM20C;COL6A3;PLIN2;ITGA5;MOCS1;ERGIC1;TMEM45A;KANK1;C1S;ADCY3;TFPI;FSTL1;TMEM165;HDAC7;KIAA0513;MTHFD1L;CLMN;PTK2B;PTPN18;GALNT6;GSN;NEGR1;TPK1;CCDC57;TXNRD1;GSR;SUSD1;LHFPL2;MERTK;KLF9;IL18R1"
gluc_response <- str_split(gluc_response[1],';')[[1]]
gluc_response_mouse <- human2mouse$mouse_name[match(gluc_response,human2mouse$human_name)]

bulk_RV <- AddModuleScore(bulk_RV,list(gluc_response_mouse),name='nr3c1')

VlnPlot(bulk_RV,'nr3c11')














human2mouse <- read.csv('~/Downloads/hdWGCNA_TOM/human2mouse.csv',header=F)
idx <- match(unique(human2mouse[,2]),human2mouse[,2])
human2mouse<-human2mouse[idx,]
colnames(human2mouse) <-c('human_name', 'mouse_name')

consensus_modules <- read.csv("~/Downloads/hdWGCNA_TOM/bulk_heart_modules.csv")
consensus_modules <- consensus_modules[,1:3]

idx_match<- match(consensus_modules$gene_name,human2mouse$human_name)

consensus_modules$gene_name <- human2mouse$mouse_name[idx_match]
consensus_modules <- consensus_modules[!is.na(consensus_modules$gene_name),]
rownames(consensus_modules) <- consensus_modules$gene_name

consensus_modules <- subset(consensus_modules, gene_name %in% rownames(M2))
# remove duplicate gene names
consensus_modules <- consensus_modules[match(unique(consensus_modules$gene_name), consensus_modules$gene_name),]
library(dplyr)
score_calc <- consensus_modules %>% group_by(module) %>% group_split()
module_colors <- unique(unlist(lapply(score_calc,'[[','module')))

mapping <- labels2colors(1:100)

module_colors <- paste0('M',match(module_colors,mapping))
M2 <- AddModuleScore(M2,lapply(score_calc,'[[','gene_name'),name="module_score")
cols_current <- colnames(M2@meta.data)
cols_current[startsWith(colnames(M2@meta.data),'module_score')] <- paste0('module_',module_colors)
colnames(M2@meta.data) <- cols_current

M2 <- SetIdent(M2,value='Subnames')

modules_int <- c('M1','M3','M4',"M8")

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'PAB_myeloid_dot_subclust.pdf'), width=5, height=2)

p <- DotPlot(M2,paste0('module_',modules_int),group.by='Subnames',dot.min=0,col.min=-1,col.max=1,scale.min=50,scale.max=100) +
  RotatedAxis() + ylab('')+ xlab('')+
  scale_color_gradient2(high='red', mid='grey95', low='blue') +
  theme(
    panel.border = element_rect(size=1,fill=NA, color='black'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()
) 
p
dev.off()

M2$group <- factor(M2$group,levels = c('Nor','Mod','Sev'))



pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'PAB_myeloid_dot_disease.pdf'), width=5, height=2.5)

p <- DotPlot(M2,paste0('module_',modules_int),group.by='group',dot.min=0,col.min=-1,col.max=1,scale.min=50,scale.max=100) +
  RotatedAxis() + ylab('')+ xlab('')+
  scale_color_gradient2(high='red', mid='grey95', low='blue') +
  theme(
    panel.border = element_rect(size=1,fill=NA, color='black'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()
) 
p
dev.off()


#MHCII

M2 <- AddModuleScore(M2,list(c('Ciita','Cd74','H2-Ab1','H2-Aa',
	'H2-Eb1','H2-Eb2','H2-Ob','H2-DMb1','H2-DMb2','H2-DMa','H2-Oa')),name='MHCII')

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'PAB_myeloid_MHC.pdf'), width=3, height=3)

VlnPlot(subset(M2,Subnames=='HLA'),'MHCII1',group.by='group',pt.size=0)
dev.off()



gluc_response <- "VIT;VKORC1L1;ERRFI1;AHCYL1;STEAP4;TRAF3IP2;GALNT15;SERPINE1;JADE1;SLA;CBLB;MT1X;EPS8;CCND3;BMPER;RASSF4;RPS6KA2;ANPEP;C1RL;MAP3K6;IL6R;PDGFRA;MLIP;SCARA5;IL1R1;EBF1;TTC7A;CRISPLD2;SPARCL1;FKBP5;NNMT;LPAR1;SLC1A3;PLA2G5;NID1;ACACB;ZFP36L2;PIK3R5;C3;SCFD2;LPXN;HACL1;SRGAP2;SLC38A2;SLC19A2;S100A10;KLHL29;GADD45B;ZBTB16;ELL2;CORO2B;IGF2R;NFATC4;DERA;SULT1B1;MAFB;BCL6;TMEM236;TBXAS1;NDUFAF2;RGL3;SERPINA3;MCFD2;PTPRS;ELN;PTEN;FMN1;HIF3A;TFCP2L1;PTH1R;SYNE3;CTSS;PTPRG;RNF157;ADAMTS2;C1QTNF1;IMPA2;SH3PXD2B;FLVCR2;EFHD1;AOX1;CERS6;ZHX3;KLF13;ANXA2;IFNGR1;GPX3;NCOA3;SLC39A11;NGF;OSMR;SLC39A14;TGFBR2;TGFBR3;PSMA6;ARHGAP10;MMP14;TBC1D2;SLC7A7;SLC7A8;GFOD1;DPYD;PICK1;FAM20C;COL6A3;PLIN2;ITGA5;MOCS1;ERGIC1;TMEM45A;KANK1;C1S;ADCY3;TFPI;FSTL1;TMEM165;HDAC7;KIAA0513;MTHFD1L;CLMN;PTK2B;PTPN18;GALNT6;GSN;NEGR1;TPK1;CCDC57;TXNRD1;GSR;SUSD1;LHFPL2;MERTK;KLF9;IL18R1"
gluc_response <- str_split(gluc_response[1],';')[[1]]


gluc_response_mouse <- human2mouse$mouse_name[match(gluc_response,human2mouse$human_name)]


M2 <- AddModuleScore(M2,list(gluc_response_mouse),name='nr3c1')


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'PAB_myeloid_nr3c1.pdf'), width=3, height=3)

VlnPlot(M2,'nr3c11',group.by='group',pt.size=0)
dev.off()

M1$group <- M1$orig.ident
pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'Wnt1.pdf'), width=3, height=3)
DimPlot(M1)
dev.off()


#######################################
#############  FIGURE S6H  ############
#######################################
library(enrichR)


bulk_modules <- consensus_modules
bulk_modules$module <- match(bulk_modules$module,mapping)
dbs <- c('GO_Biological_Process_2023','GO_Cellular_Component_2023','GO_Molecular_Function_2023','ChEA_2022','Reactome_2022')

#Run enrichment by cell type
M2$Names_group <- paste0(M2$Subnames,'_',M2$group)
Idents(M2) <- "Names_group"
combined_set <- data.frame()
combined_output <- data.frame()

mods_idx <- c(1,3,4,8)
cell_types <- unique(M2$Subnames)
comparison <- list(c("Sev","Mod"),c("Mod","Nor"),c("Sev","Nor"))
for (i in mods_idx){
	for (j in cell_types){
		for (k in comparison){
			key_genes <- subset(bulk_modules,module %in% c(i))$gene_name
			key_genes <- key_genes[key_genes %in% rownames(M2)]

			gene_set <- FindMarkers(M2, ident.1 = paste0(j,"_",k[1]), ident.2 = paste0(j,"_",k[2]),features=key_genes,recorrect_umi=F)
			
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


#selected_terms <- subset(combined_output,db=="ChEA_2022")
selected_terms <- subset(combined_output,direction=="down")
selected_terms <- subset(selected_terms,comparison=="Sev_Nor")
selected_terms <- subset(selected_terms,Adjusted.P.value<0.05)
#selected_terms <- subset(selected_terms,module %in% c('M1'))
selected_terms



Idents(M2) <- "group"

library(EnhancedVolcano)

gene_set <- FindMarkers(M2, ident.1 = 'Sev', ident.2 = 'Nor',recorrect_umi=F,features=subset(bulk_modules,module==1)$gene_name)

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'PAB_M1_Myeloid_module_volcano_all_RVF_vs_NF.pdf'), width=8, height=6)

EnhancedVolcano(gene_set,lab=rownames(gene_set),
	x='avg_log2FC',y='p_val_adj',
	FCcutoff = 0.1,pCutoff=0.05) + coord_flip()
dev.off()

gene_set <- FindMarkers(M2, ident.1 = 'Sev', ident.2 = 'Nor',recorrect_umi=F,features=subset(bulk_modules,module==8)$gene_name)

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'PAB_M8_Myeloid_module_volcano_all_RVF_vs_NF.pdf'), width=8, height=6)

EnhancedVolcano(gene_set,lab=rownames(gene_set),
	x='avg_log2FC',y='p_val_adj',
	FCcutoff = 0.1,pCutoff=0.05) + coord_flip()
dev.off()

#######################################
#############  FIGURE S6I  ############
#######################################

M2 <- subset(M1, Names %in% c("CM"))

M2 <- RunPCA(M2)
M2 <- RunHarmony(M2,'patient')
M2 <- FindNeighbors(M2, dims = 1:50,reduction = "harmony")
M2 <- FindClusters(M2, resolution = 0.5,reduction = "harmony")
M2 <- RunUMAP(M2, dims = 1:50,reduction = "harmony")

M2$Names <- M2@active.ident
markers<-FindAllMarkers(M2,recorrect_umi=F)

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'PAB_CM_snUMAP.pdf'), width=5, height=5)
PlotEmbedding(M2,group.by='Names',point_size=1,plot_under=TRUE,plot_theme=umap_theme()+NoLegend(),raster_dpi=400,raster_scale=0.5)
dev.off()



M3 <- readRDS(file = "/Volumes/Extreme\ SSD/Final_Analysis/CellTypes/cm_subclust.rds")

new.cluster.ids <- c("Cm1","Cm2","Cm3","Cm4","Cm5","Cm6","Cm7","Cm8","Cm9","Cm10")
names(new.cluster.ids) <- levels(M3)
M3 <- RenameIdents(M3, new.cluster.ids)

M3$Subnames <- M3@active.ident
M3$SubNames_Groups <- paste(M1$Subnames,M3$group,sep='_')


M3 <- AddModuleScore(M3, features=list(c('MALAT1')),assay="SCT",name="Clust0Score")
M3 <- AddModuleScore(M3, features=list(c('FGF12','SH3RF2','KCNMB2','PRELID2')),assay="SCT",name="Clust1Score")
M3 <- AddModuleScore(M3, features=list(c('TNNT2','TTN','MYBPC3','MYH7')),assay="SCT",name="Clust2Score")
M3 <- AddModuleScore(M3, features=list(c('PALLD','MYO18B','MYPN','ANKRD1')),assay="SCT",name="Clust5Score")
M3 <- AddModuleScore(M3, features=list(c('PDE3A','CDH2','PDLIM5')),assay="SCT",name="Clust4Score")
M3 <- AddModuleScore(M3, features=list(c('AKAP13','OBSCN','LARGE1','THSD4')),assay="SCT",name="Clust3Score")
M3 <- AddModuleScore(M3, features=list(c('PALLD','SORBS2','CAMK2D','CCSER1','PDLIM5')),assay="SCT",name="Clust7Score")
M3 <- AddModuleScore(M3, features=list(c('AC020637.1','LINC02388')),assay="SCT",name="Clust6Score")
M3 <- AddModuleScore(M3, features=list(c('MIR646HG')),assay="SCT",name="Clust8Score")
M3 <- AddModuleScore(M3, features=list(c('GPC5','HS6ST3')),assay="SCT",name="Clust9Score")

#DefaultAssay(M3) <- "RNA"

#M3[["RNA"]] <- split(M3[["RNA"]], f = M3$patient)
#M3[['SCT']] <- NULL
#M3[['decontXcounts']] <- NULL
#M3 <- SCTransform(M3, vst.flavor = "v2")
#M3 <- RunPCA(M3, npcs = 50, verbose = FALSE)
M3 <- SplitObject(M3, split.by = "patient")
M3<-PrepSCTIntegration(M3)
features<-SelectIntegrationFeatures(M3)
M3.anchors<-FindIntegrationAnchors(M3,normalization.method = 'SCT',anchor.features = features, reduction = "rpca")
M3 <- IntegrateData(anchorset = M3.anchors,normalization.method='SCT')

DefaultAssay(M3) <- "integrated"

M3 <- RunPCA(M3, npcs = 50, verbose = FALSE)
M3 <- RunUMAP(M3, reduction = "pca", dims = 1:30)

# RenameGenesSeurat <- function(obj = ls.Seurat[[i]],
#                               newnames = HGNC.updated[[i]]$Suggested.Symbol,
#                               assay = "RNA",
#                               slots = c("data", "counts", "meta.features")) {
#   #
#   message("RenameGenesSeurat, assay: ", assay)
#   warning("Run this before integration and downstream processing. It only attempts to change
#           @counts, @data, and @meta.features in obj@assays$YOUR_ASSAY.", immediate. = TRUE)

#   stopifnot(
#     "Unequal gene name sets: nrow(assayobj) != nrow(newnames):" =
#       length(Features(obj, assay = assay)) == length(newnames)
#   )

#   if (obj@version < 5) warning("obj@version < 5. Old versions are not supported. Update the obj!", immediate. = TRUE)

#   if ("scale.data" %in% slots) {
#     n_genes_sc_dta <- nrow(obj@assays[[assay]]$"scale.data")
#     stopifnot(
#       "scale.data does has different number of genes than newnames!" =
#         n_genes_sc_dta == length(newnames)
#     )
#   }

#   LayersFound <- SeuratObject::Layers(obj@assays[[assay]])
#   #print("Present: ", sort(LayersFound))

#   slots <- sort(intersect(slots, LayersFound))
#   #print("Replaced: ", slots)

#     for (slotX in slots) {
#     print(slotX)
#     if (slotX == "scale.data") browser()
#     nrO <- nrow(SeuratObject::GetAssayData(object = obj, assay = assay, layer = slotX))
#     obj <- .check_and_rename(obj, assay, newnames = newnames, layer.name = slotX)
#     nrN <- nrow(SeuratObject::GetAssayData(object = obj, assay = assay, layer = slotX))
#     stopifnot(nrN == nrO)
#   }
#   return(obj)
# }


# .check_and_rename <- function(obj, assay, newnames, layer.name) {
#   cat(layer.name, fill = TRUE)

#   length_newnames <- length(newnames)
#   length_orig_names <- length(Features(obj, assay = assay))

#   stopifnot(
#     is(obj, "Seurat"),
#     is.character(assay),
#     is.character(layer.name),
#     is.character(newnames),
#     length_orig_names == length_newnames
#   )

#   assayobj <- obj@assays[[assay]]
#   feature.list <- rownames(assayobj@features@.Data)

#   if (length(feature.list) == length(newnames)) {
#     rownames(assayobj@features@.Data) <- newnames
#     nrX <- length(rownames(assayobj@features@.Data))
#   } else {
#     paste0("length feature.list", length(feature.list), "length newnames", length(newnames))
#     stop()
#   }

#   if (layer.name %in% SeuratObject::Layers(assayobj)) {
#     matrix_n <- SeuratObject::LayerData(assayobj, layer = layer.name)
#     nr1 <- nrow(matrix_n)

#     if (all(dim(matrix_n)) > 0) {
#       # browser()
#       stopifnot(nrow(matrix_n) == length(newnames))

#       if ("dgCMatrix" %in% class(matrix_n)) {
#         message(assay, "@", layer.name, " is of type dgeCMatrix!")
#         matrix_n@Dimnames[[1]] <- newnames
#       } else if ("matrix" %in% class(matrix_n)) {
#         message(assay, "@", layer.name, " is of type Matrix!")
#         rownames(matrix_n) <- newnames
#       } else if ("data.frame" %in% class(matrix_n)) {
#         message(assay, "@", layer.name, " is of type data.frame!")
#         rownames(matrix_n) <- newnames
#       } else {
#         warning(">>> No renaming: ", assay, "@", layer.name,
#           " not of type dgeCMatrix / Matrix / data.frame.",
#           immediate. = TRUE
#         )
#       }
#       stopifnot(nr1 == nrow(matrix_n))

#       SeuratObject::LayerData(assayobj, layer = layer.name) <- matrix_n
#       nr3 <- nrow(SeuratObject::LayerData(assayobj, layer = layer.name))
#       stopifnot(nr3 == nrX)
#     }
#   } else {
#     warning(paste(">>>", assay, "@", layer.name, "does not exist!"), immediate. = TRUE)
#   }
#   # obj <- SetAssayData(obj, layer = layer.name, new.data = matrix_n)
#   obj@assays[[assay]] <- assayobj
#   return(obj)
# }

# DefaultAssay(M2) <- "RNA"


# newnames <- human2mouse$human_name[match(rownames(M2),human2mouse$mouse_name)]
# newnames[is.na(newnames)] <- rownames(M2)[is.na(newnames)]

# RenameGenesSeurat(M2,newnames = newnames,assay='RNA')



RNA <- M2@assays$RNA['counts']
newnames <- human2mouse$human_name[match(rownames(RNA),human2mouse$mouse_name)]
newnames[is.na(newnames)] <- rownames(RNA)[is.na(newnames)]
rownames(RNA) <- newnames


M4 <- CreateAssayObject(RNA)
M2[['humanized']] <- M4
DefaultAssay(M2) <- "humanized"
M2[['SCT']] <- NULL

M2 <- SCTransform(M2, vst.flavor = "v2",assay='humanized')
M2 <- RunPCA(M2, npcs = 50, verbose = FALSE)

anchors <- FindTransferAnchors(
  reference = M3,
  query = M2,
  normalization.method = "SCT",
  recompute.residuals=FALSE,
  reference.reduction = "pca",
  dims = 1:50
)

predictions <- TransferData(anchorset = anchors, refdata = M3$Subnames, dims = 1:50)


M3 <- RunUMAP(M3, dims = 1:50, return.model = TRUE)
M2 <- MapQuery(anchorset = anchors, reference = M3, query = M2,
	refdata = list(celltype = "Subnames"), reference.reduction = "pca", reduction.model = "umap")

score <- MappingScore(anchors)

M2$map_score <- score

p1 <- DimPlot(M3, reduction = "umap", group.by = "Subnames", label = TRUE, label.size = 3, repel = TRUE,raster=TRUE,pt.size=1.5) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(M2, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE, label.size = 3, pt.size=1.5,repel = TRUE,raster=TRUE) + NoLegend() + ggtitle("Query transferred labels")
pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'RV_PAB_CM_ref_mapped.pdf'), width=10, height=5)
p1 + p2
dev.off()


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'PAB_CM_ref_mapped.pdf'), width=5, height=5)
PlotEmbedding(M2,group.by='predicted.celltype',reduction = "ref.umap",point_size=0.2,plot_under=TRUE,plot_theme=umap_theme()+NoLegend(),raster_dpi=400,raster_scale=0.5)
dev.off()


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'RV_CM_ref_mapped.pdf'), width=5, height=5)
PlotEmbedding(M3,group.by='Subnames',point_size=0.2,plot_under=TRUE,plot_theme=umap_theme()+NoLegend(),raster_dpi=400,raster_scale=0.5)
dev.off()


FeaturePlot(M2,'map_score',reduction = "ref.umap")

#######################################
#############  FIGURE S6J  ############
#######################################


RV_marks <- FindAllMarkers(M)

RV_marks_sig <- subset(RV_marks, p_val_adj<0.05 & avg_log2FC>0)

idx<-match(RV_marks_sig$cluster,unique(RV_marks_sig$cluster))

marks <- split(RV_marks_sig$gene,idx)

M2 <- AddModuleScore(M2,marks,name='RV_marks')
DefaultAssay(M3) <- "SCT"
M3 <- AddModuleScore(M3,marks,name='RV_marks')

M2$predicted.celltype <- factor(M2$predicted.celltype,levels=levels(M3))

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'PAB_CM_marker_scores.pdf'), width=6, height=4)

DotPlot(M2,c('RV_marks1','RV_marks2','RV_marks3','RV_marks4','RV_marks5'
	,'RV_marks6','RV_marks7','RV_marks8','RV_marks9','RV_marks10'),
	col.min = 0,group.by='predicted.celltype') +
scale_x_discrete(labels=c('Cm1','Cm2','Cm3','Cm4','Cm5','Cm6','Cm7','Cm8','Cm9','Cm10'))
dev.off()


DotPlot(M3,c('RV_marks1','RV_marks2','RV_marks3','RV_marks4','RV_marks5'
	,'RV_marks6','RV_marks7','RV_marks8','RV_marks9','RV_marks10'),col.min = 0) +
scale_x_discrete(labels=c('Cm1','Cm2','Cm3','Cm4','Cm5','Cm6','Cm7','Cm8','Cm9','Cm10'))

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'PAB_CM_scores.pdf'), width=3.75, height=4.25)
VlnPlot(M2,'map_score',group.by='predicted.celltype',pt.size=0)
dev.off()

#######################################
#############  FIGURE S6L  ############
#######################################

M3 <- SetIdent(M3,value = 'group')
M2 <- SetIdent(M2,value = 'group')


a <- FindMarkers(M3,ident.1='RVF',ident.2='NF')
b <-  FindMarkers(M2,ident.1='Sev',ident.2='Nor')
shared <- intersect(rownames(a),rownames(b))
dataset <- data.frame(RV=a[shared,]$avg_log2FC,PAB=b[shared,]$avg_log2FC)
rownames(dataset) <- shared
labs <- rownames(dataset)
#labs[abs(dataset$PAB - dataset$RV)<1] <- NA


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'PAB_vs_RV_CM__dot.pdf'), width=6, height=8)
ggplot(dataset, aes(x = RV, y=PAB)) + geom_point() + 
  geom_text_repel(label=labs,max.overlaps=50) + theme_classic()
dev.off()


a <- FindMarkers(M3,ident.1='RVF',ident.2='NF')
b <-  FindMarkers(M2,ident.1='Mod',ident.2='Nor')
shared <- intersect(rownames(a),rownames(b))
dataset <- data.frame(RV=a[shared,]$avg_log2FC,PAB=b[shared,]$avg_log2FC)
rownames(dataset) <- shared
labs <- rownames(dataset)
#labs[abs(dataset$PAB - dataset$RV)<1] <- NA


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'PAB_vs_RV_CM_Mod_Nor_dot.pdf'), width=6, height=8)
ggplot(dataset, aes(x = RV, y=PAB)) + geom_point() + 
  geom_text_repel(label=labs,max.overlaps=50) + theme_classic()
dev.off()


a <- FindMarkers(M2,ident.1='Sev',ident.2='Nor')
b <-  FindMarkers(M2,ident.1='Mod',ident.2='Nor')
shared <- intersect(rownames(a),rownames(b))
dataset <- data.frame(RV=a[shared,]$avg_log2FC,PAB=b[shared,]$avg_log2FC)
rownames(dataset) <- shared
labs <- rownames(dataset)
#labs[abs(dataset$PAB - dataset$RV)<1] <- NA


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'PAB_CM_Sev_Nor_Mod_Nor_dot.pdf'), width=6, height=8)
ggplot(dataset, aes(x = RV, y=PAB)) + geom_point() + 
  geom_text_repel(label=labs,max.overlaps=50) + theme_classic()
dev.off()


a <- FindMarkers(M2,ident.1='Sev',ident.2='Mod')
b <-  FindMarkers(M2,ident.1='Mod',ident.2='Nor')
shared <- intersect(rownames(a),rownames(b))
dataset <- data.frame(RV=a[shared,]$avg_log2FC,PAB=b[shared,]$avg_log2FC)
rownames(dataset) <- shared
labs <- rownames(dataset)
labs[abs(dataset$PAB - dataset$RV)<1] <- NA


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'PAB_CM_Sev_Mod_Mod_Nor_dot.pdf'), width=6, height=8)
ggplot(dataset, aes(x = RV, y=PAB)) + geom_point() + 
  geom_text_repel(label=labs,max.overlaps=10) + theme_classic()
dev.off()
#######################################
#############  FIGURE S6M  ############
#######################################

consensus_modules <- read.csv("~/Downloads/hdWGCNA_TOM/bulk_heart_modules.csv")
consensus_modules <- consensus_modules[,1:3]
consensus_modules <- subset(consensus_modules, gene_name %in% rownames(M2))
consensus_modules <- consensus_modules[match(unique(consensus_modules$gene_name), consensus_modules$gene_name),]
library(dplyr)
mapping <- labels2colors(1:100)

score_calc <- consensus_modules %>% group_by(module) %>% group_split()
module_colors <- unique(unlist(lapply(score_calc,'[[','module')))
module_colors <- paste0('M',match(module_colors,mapping))


M2 <- AddModuleScore(M2,lapply(score_calc,'[[','gene_name'),name="module_score",ctrl = 50)

cols_current <- colnames(M2@meta.data)
cols_current[startsWith(colnames(M2@meta.data),'module_score')] <- paste0('module_',module_colors)
colnames(M2@meta.data) <- cols_current


M2 <- SetIdent(M2, value = "group")
M2$group <- factor(M2$group,levels=c('Nor','Mod','Sev'))


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'PAB_seurat_dot_CM.pdf'), width=5, height=2.5)

p <- DotPlot(M2,paste0('module_',
  c('M2','M10','M12','M25','M26','M28')),dot.min=0,col.min=0,col.max=2,group.by='group') +
  RotatedAxis() + ylab('')+ xlab('')+
  scale_color_gradient2(high='red', mid='grey95', low='blue') +
  theme(
    panel.border = element_rect(size=1,fill=NA, color='black'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()
) 
p

dev.off()

#######################################
#############  FIGURE S6N  ############
#######################################
library("readxl")

Human.Mito <- read_excel("~/Downloads/hdWGCNA_TOM/Human.MitoCarta3.0.xls", sheet = "A Human MitoCarta3.0")
Mouse.Mito <- read_excel("~/Downloads/hdWGCNA_TOM/Mouse.MitoCarta3.0.xls", sheet = "A Mouse MitoCarta3.0")


M3 <- AddModuleScore(M3,list(Human.Mito$Symbol),name='mito')
M2 <- AddModuleScore(M2,list(union(Human.Mito$Symbol,Mouse.Mito$Symbol)),name='mito',ctrl = 50)

p1<-VlnPlot(M2,'mito1',group.by='group',pt.size=0)
p2<-VlnPlot(M3,'mito1',group.by='group',pt.size=0)

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'PAB_RV_CM_mitocarto.pdf'), width=3, height=4)

p1 / p2
dev.off()

anno <- trimws(unlist(lapply(lapply(lapply(Human.Mito$MitoCarta3.0_MitoPathways,str_split,'>'),'[[',1),'[[',1)))

anno[anno == "Small molecule transport | Signaling"] = "Small molecule transport" 


a<-FindMarkers(M3,ident.1='RVF',ident.2='NF',features = intersect(Human.Mito$Symbol,rownames(M3)))
a$p_val_adj[a$p_val_adj < 1e-50] = 1e-50

keyvals <- anno[match(rownames(a),Human.Mito$Symbol)]
library(colormap)

colors <- colormap(colormap=colormaps$rainbow_soft, nshades=length(unique(keyvals)))

colors <- colors[match(keyvals,unique(keyvals))]
names(colors) <- keyvals


library(EnhancedVolcano)
pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'CM_RV_Mito.pdf'), width=6, height=11)

EnhancedVolcano(a,lab=rownames(a),
  x='avg_log2FC',y='p_val_adj',
  FCcutoff = 0.1,pCutoff=0.05,colCustom = colors,xlim=c(-6.25,3),ylim=c(0,52))
dev.off()
